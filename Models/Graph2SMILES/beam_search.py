import argparse
import datetime
import glob
import numpy as np
import os
import sys
import time
import torch
import torch.distributed as dist
import tree_builder
from tree_builder import NodeWrapper
import preprocess
import tqdm
import pickle
import warnings
import gc
import heapq
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
from models.graph2smiles import Graph2SMILES
from torch.utils.data import DataLoader, SequentialSampler
from torch.utils.data.distributed import DistributedSampler
from train import get_model
from utils import parsing
from utils.data_utils import canonicalize_smiles, load_vocab, G2SDataset, G2SDataset_Beam
from utils.train_utils import log_tensor, log_rank_0, param_count, set_seed, setup_logger


warnings.filterwarnings("ignore", category=UserWarning)
depth_visited = set()


def top_K_search(root, K):
    # Priority queue to keep track of top-K nodes
    top_K = []

    def dfs(node, cumulative_score):
        if node is None:
            return

        # Update the cumulative score
        cumulative_score += float(node.score)

        if node.parent:
            # Wrap the node and push it into the heap
            heapq.heappush(top_K, NodeWrapper(node, cumulative_score))
            # If the heap exceeds size K, remove the smallest element
            if len(top_K) > K:
                heapq.heappop(top_K)

        # Traverse the children
        for child in node.childrens:
            dfs(child, cumulative_score)

    # Start the depth-first search from the root with an initial cumulative score of 0
    dfs(root, 0)

    # Extract the top-K nodes from the heap
    top_K_nodes = [heapq.heappop(top_K).node for _ in range(len(top_K))]
    top_K_nodes.reverse()  # To get the nodes in descending order of cumulative score

    for topK_node in top_K_nodes:
        if topK_node.match:
            return True
        
    return False


def chunk_Nodes(length, chunk_size):
    
    # looping till length l
    for i in range(0, len(length), chunk_size): 
        yield length[i:i + chunk_size]


def get_beam_search_parser():
    parser = argparse.ArgumentParser("Beam Search")
    parsing.add_common_args(parser)
    parsing.add_preprocess_args(parser)
    parsing.add_train_args(parser)
    parsing.add_beam_search_args(parser)
    return parser


def get_predictions(args, model, vocab_tokens, test_loader, device, start):
    all_predictions = []
    with torch.no_grad():
        for test_batch in test_loader:

            test_batch.to(device)
            results = model.predict_step(
                reaction_batch=test_batch,
                batch_size=test_batch.size,
                beam_size=args.beam_size,
                n_best=args.n_best,
                temperature=args.temperature,
                min_length=args.predict_min_len,
                max_length=args.predict_max_len
            )

            for predictions, scores in zip(results["predictions"], results["scores"]):
                smis_with_scores = []
                for prediction, score in zip(predictions, scores):
                    predicted_idx = prediction.detach().cpu().numpy()
                    score = score.detach().cpu().numpy()
                    predicted_tokens = [vocab_tokens[idx] for idx in predicted_idx[:-1]]
                    smi = "".join(predicted_tokens)
                    smis_with_scores.append(f"{smi}_{score}")
                smis_with_scores = ",".join(smis_with_scores)
                all_predictions.append(f"{smis_with_scores}\n")

    return all_predictions


def inference(Node, model, vocab_tokens, device, start, step=0):
    data_dict = preprocess.binarize_g2s_beam_search(args, Node.get_rct_smi(), Node.get_pdt_smi())
    dataset_class = G2SDataset
    args.compute_graph_distance = True
    test_dataset = dataset_class(args, file=data_dict)
    test_dataset.batch(
        batch_type=args.batch_type,
        batch_size=args.predict_batch_size
    )

    if args.local_rank != -1:
        test_sampler = DistributedSampler(test_dataset, shuffle=False)
    else:
        test_sampler = SequentialSampler(test_dataset)

    test_loader = DataLoader(
                    dataset=test_dataset,
                    batch_size=1,
                    sampler=test_sampler,
                    num_workers=args.num_cores,
                    collate_fn=lambda _batch: _batch[0],
                    pin_memory=True
                )

    all_predictions = get_predictions(
        args, model, vocab_tokens, test_loader, device, start)

    for result in all_predictions:
        smi_and_score_list = result.replace("\n", "").split(",")
        # TODO: Assign results to leaf nodes and how to get to the bottom of the tree.
        ct = 0
        for rank, smi_and_score in enumerate(smi_and_score_list):
            smi, score = smi_and_score.split("_")[0], smi_and_score.split("_")[1]
            if not Chem.MolFromSmiles(smi.replace("END", ''), sanitize=False) or ct >= args.topK:
                continue
            if 'END' not in smi:
                rct_smi = Node.get_rct_smi().replace('\n', '')
                smi_list = rct_smi.split('.')
                smi_to_be_replaced = None
                similarity_high = 0
                for original_smi in smi_list:
                    fp1 = FingerprintMols.FingerprintMol(Chem.MolFromSmarts(original_smi))
                    fp2 = FingerprintMols.FingerprintMol(Chem.MolFromSmarts(smi))
                    if DataStructs.TanimotoSimilarity(fp1, fp2) > similarity_high:
                        similarity_high = DataStructs.TanimotoSimilarity(fp1, fp2)
                        smi_to_be_replaced = original_smi

                if smi_to_be_replaced is not None:
                    index_to_be_replaced = smi_list.index(smi_to_be_replaced)
                    smi_list[index_to_be_replaced] = smi
                else:
                    smi_list.append(smi)
                smi = '.'.join(smi_list)

            new_Node = tree_builder.Node((smi, rank, score), Node)
            Node.childrens.append(new_Node)
            ct += 1

    if Node.depth not in depth_visited:
        depth_visited.add(Node.depth + 1)
        log_rank_0(f"Depth: {Node.depth + 1}")

    del data_dict, test_dataset, test_loader, all_predictions
    torch.cuda.empty_cache()
    gc.collect()
        
    if Node.childrens:
        if Node.childrens[0].depth != args.end_step:
            for new_Node in Node.childrens[:args.topK]:
                if not new_Node.terminated():
                    inference(new_Node, model, vocab_tokens, device, start, step+1)
    else:
        Node.match = False
        Node.termination = True
    
    return Node.match


def batch_inference(Node_batch, model, vocab_tokens, device, start, step=0):
    src_list = [Node.get_rct_smi().strip() for Node in Node_batch]
    tgt_list = [Node.get_pdt_smi().strip() for Node in Node_batch]
    data_dict = preprocess.binarize_g2s_beam_search(args, src_list, tgt_list)

    dataset_class = G2SDataset
    args.compute_graph_distance = True
    test_dataset = dataset_class(args, file=data_dict)
    test_dataset.batch(
        batch_type=args.batch_type,
        batch_size=args.predict_batch_size
    )

    if args.local_rank != -1:
        test_sampler = DistributedSampler(test_dataset, shuffle=False)
    else:
        test_sampler = SequentialSampler(test_dataset)

    test_loader = DataLoader(
                    dataset=test_dataset,
                    batch_size=1,
                    sampler=test_sampler,
                    num_workers=args.num_cores,
                    collate_fn=lambda _batch: _batch[0],
                    pin_memory=True
                )

    all_predictions = get_predictions(
        args, model, vocab_tokens, test_loader, device, start)
    
    all_predictions = zip(all_predictions, Node_batch)

    child_Nodes = []

    for result, Node in all_predictions:
        smi_and_score_list = result.replace("\n", "").split(",")
        # TODO: Assign results to leaf nodes and how to get to the bottom of the tree.
        ct = 0
        for rank, smi_and_score in enumerate(smi_and_score_list):
            smi, score = smi_and_score.split("_")[0], smi_and_score.split("_")[1]
            if not Chem.MolFromSmiles(smi.replace("END", ''), sanitize=False) or ct >= args.topK:
                continue
            if 'END' not in smi:
                rct_smi = Node.get_rct_smi().replace('\n', '')
                smi_list = rct_smi.split('.')
                smi_to_be_replaced = None
                similarity_high = 0
                for original_smi in smi_list:
                    fp1 = FingerprintMols.FingerprintMol(Chem.MolFromSmarts(original_smi))
                    fp2 = FingerprintMols.FingerprintMol(Chem.MolFromSmarts(smi))
                    if DataStructs.TanimotoSimilarity(fp1, fp2) > similarity_high:
                        similarity_high = DataStructs.TanimotoSimilarity(fp1, fp2)
                        smi_to_be_replaced = original_smi

                if smi_to_be_replaced is not None:
                    index_to_be_replaced = smi_list.index(smi_to_be_replaced)
                    smi_list[index_to_be_replaced] = smi
                else:
                    smi_list.append(smi)
                smi = '.'.join(smi_list)
                
            new_Node = tree_builder.Node((smi, rank, score), Node)
            if not new_Node.terminated():
                child_Nodes.append(new_Node)
            Node.childrens.append(new_Node)
            ct += 1
        if not Node.childrens:
            Node.match = False
            Node.termination = True

    del data_dict, test_dataset, test_loader, all_predictions
    torch.cuda.empty_cache()
    gc.collect()

    if Node_batch[0].depth not in depth_visited:
        depth_visited.add(Node.depth)
        log_rank_0(f"Depth: {Node.depth + 1}")
        
    if child_Nodes:
        if child_Nodes[0].depth != args.end_step:
            chunked_Nodes = list(chunk_Nodes(child_Nodes, args.batch_search_size))
            for chunked_Node in chunked_Nodes:
                batch_inference(chunked_Node, model, vocab_tokens, device, start, step+1)


def main(args):
    args.device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
    device = args.device
    if args.local_rank != -1:
        dist.init_process_group(backend=args.backend, init_method='env://', timeout=datetime.timedelta(0, 7200))
        torch.cuda.set_device(args.local_rank)
        torch.backends.cudnn.benchmark = True

    if torch.distributed.is_initialized():
        log_rank_0(f"Device rank: {torch.distributed.get_rank()}")
    os.makedirs(os.path.join("./results", args.data_name), exist_ok=True)

    parsing.log_args(args, phase="prediction")

    # initialization ---------------------------------- ckpt parsing
    if args.do_validate:
        checkpoints = glob.glob(os.path.join(args.load_from, "*.pt"))
        checkpoints = sorted(
            checkpoints,
            key=lambda ckpt: int(ckpt.split(".")[-2].split("_")[-1]),
            reverse=True
        )
        checkpoints = [ckpt for ckpt in checkpoints
                       if (args.checkpoint_step_start <= int(ckpt.split(".")[-2].split("_")[0]))
                       and (args.checkpoint_step_end >= int(ckpt.split(".")[-2].split("_")[0]))]
        file_bin = os.path.join(args.processed_data_path, "val.npz")
        file_tgt = os.path.join(args.processed_data_path, "tgt-val.txt")
    elif args.do_beam_search:
        checkpoints = [os.path.join(args.model_path, args.load_from)]
        file_bin = os.path.join(args.processed_data_path, "unseen.npz")
        file_tgt = os.path.join(args.processed_data_path, "tgt-unseen.txt")
        file_src = os.path.join(args.processed_data_path, "src-unseen.txt")
        file_rgt = os.path.join(args.processed_data_path, "rgt-unseen.txt")
    else:
        raise ValueError("Either --do_validate or --do_beam_search need to be specified!")

    tree = tree_builder.initialize_tree(file_src, file_tgt, file_rgt)
    num_rxn_correct = 0
    total_num_rxn = len(tree)
    model = None
    vocab_tokens = None
    start = time.time()
    for ckpt_i, checkpoint in enumerate(checkpoints):
        tree_file = os.path.join(args.test_output_path, f"search_tree_{ckpt_i}.pkl")

        if os.path.exists(tree_file):
            log_rank_0(f"Result file found at {tree_file}, skipping prediction.")
        else:
            log_rank_0(f"Loading from {checkpoint}")
            try:
                state = torch.load(checkpoint)
            except RuntimeError:
                log_rank_0(f"Error loading {checkpoint}, skipping")
                continue            # some weird corrupted files

            pretrain_args = state["args"]
            pretrain_args.load_from = None # Hard code from tracing back to older models

            pretrain_state_dict = state["state_dict"]
            pretrain_args.local_rank = args.local_rank
            if not hasattr(pretrain_args, "n_latent"):
                pretrain_args.n_latent = 1
            args.n_latent = pretrain_args.n_latent
            if not hasattr(pretrain_args, "shared_attention_layer"):
                pretrain_args.shared_attention_layer = 0

            if model is None:
                # initialization ---------------------------------- model
                log_rank_0("Model is None, building model")
                log_rank_0("First logging args for training")
                parsing.log_args(pretrain_args, phase="training")

                # backward
                assert args.model == pretrain_args.model or \
                    pretrain_args.model == "g2s_series_rel", \
                    f"Pretrained model is {pretrain_args.model}!"
                model_class = Graph2SMILES
                args.compute_graph_distance = True

                # initialization ---------------------------------- vocab
                vocab = load_vocab(args)
                vocab_tokens = [k for k, v in sorted(vocab.items(), key=lambda tup: tup[1])]

                model, state = get_model(pretrain_args, model_class, vocab, device)
                if hasattr(model, "module"):
                    model = model.module        # unwrap DDP model to enable accessing model func directly

                log_rank_0(model)
                log_rank_0(f"Number of parameters = {param_count(model)}")

            pretrain_state_dict = {k.replace("module.", ""): v for k, v in pretrain_state_dict.items()}
            model.load_state_dict(pretrain_state_dict)
            log_rank_0(f"Loaded pretrained state_dict from {checkpoint}")
            model.eval()
            
            global depth_visited
            if args.do_batch_search:
                chunked_Nodes = list(chunk_Nodes(tree, args.batch_search_size))
                num_rxn_correct = 0
                for chunked_Node in tqdm.tqdm(chunked_Nodes):
                    batch_inference(chunked_Node, model, vocab_tokens, device, start)
                    num_rxn_correct_batch = 0
                    for node in chunked_Node:
                        if top_K_search(node, args.topK):
                            num_rxn_correct_batch += 1
                    log_rank_0(f"Batch size: {len(chunked_Node)}, number of product found: {num_rxn_correct_batch}")
                    log_rank_0(f"Batch Accuracy: {num_rxn_correct_batch / len(chunked_Node)}")
                    num_rxn_correct += num_rxn_correct_batch
                    depth_visited = set()
            else:
                rxn_searched = 0
                for leaf in tqdm.tqdm(tree):
                    match_result = inference(leaf, model, vocab_tokens, device, start)
                    rxn_searched += 1
                    if match_result:
                        num_rxn_correct += 1
                    depth_visited = set()
                    log_rank_0(f"Accuracy: {num_rxn_correct / rxn_searched}")

            with open(tree_file, 'wb') as f:
                pickle.dump(tree, f)

            if args.do_score:
                accuracy = num_rxn_correct / total_num_rxn
                log_rank_0(f"Accuracy: {accuracy}")
                log_rank_0(f"Number of correct reactions: {num_rxn_correct}")
                log_rank_0(f"Total number of reactions: {total_num_rxn}")
                log_rank_0(f"Elapsed time: {time.time() - start: .2f} s")


if __name__ == "__main__":
    # initialization ---------------------------------- args, logs and devices
    predict_parser = get_beam_search_parser()
    args = predict_parser.parse_args()

    setup_logger(args, warning_off=True)
    np.set_printoptions(threshold=sys.maxsize)
    torch.set_printoptions(profile="full")

    set_seed(args.seed)
    main(args)
