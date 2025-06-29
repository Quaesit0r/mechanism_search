#!/bin/bash

export CUDA_VISIBLE_DEVICES=0,1
export LOAD_FROM=""
export MODEL=graph2smiles
export MAX_REL_POS=4
export ACCUM_COUNT=4
export ENC_PE=none
export ENC_H=256
export ENC_EMB_SCALE=sqrt
export MAX_STEP=850000
# export BATCH_SIZE=8192
# export BATCH_TYPE=tokens_sum
# export ACCUMULATION_COUNT=4
export BATCH_SIZE=4096
export BATCH_TYPE=tokens
export ACCUMULATION_COUNT=8
export ENC_LAYER=4
export REL_BUCKETS=15
export REL_POS=emb_only
export ATTN_LAYER=6
export N_LATENT=1
export MPN_TYPE=dgcn
export DGAT_H=8
export LR=2
export DROPOUT=0.1
export NUM_NODES=1
export NUM_GPU=2


torchrun \
  --standalone \
  --node_rank=0 \
  --nnodes="$NUM_NODES"\
  --nproc_per_node="$NUM_GPU" \
    train.py \
    --backend=nccl \
    --model="$MODEL" \
    --data_name="$DATA_NAME" \
    --log_file="graph2smiles_train_$DATA_NAME" \
    --processed_data_path="$PROCESSED_DATA_PATH" \
    --model_path="$MODEL_PATH" \
    --embed_size=256 \
    --mpn_type="$MPN_TYPE" \
    --dgat_attn_heads="$DGAT_H" \
    --encoder_num_layers="$ENC_LAYER" \
    --encoder_hidden_size="$ENC_H" \
    --encoder_positional_encoding="$ENC_PE" \
    --encoder_emb_scale="$ENC_EMB_SCALE" \
    --attn_enc_num_layers="$ATTN_LAYER" \
    --attn_enc_hidden_size=256 \
    --attn_enc_heads=8 \
    --attn_enc_filter_size=2048 \
    --rel_pos="$REL_POS" \
    --rel_pos_buckets="$REL_BUCKETS" \
    --decoder_num_layers=6 \
    --decoder_hidden_size=256 \
    --decoder_attn_heads=8 \
    --decoder_filter_size=2048 \
    --dropout="$DROPOUT" \
    --attn_dropout="$DROPOUT" \
    --max_relative_positions="$MAX_REL_POS" \
    --seed=42 \
    --epoch=2000 \
    --max_steps="$MAX_STEP" \
    --warmup_steps=8000 \
    --lr="$LR" \
    --weight_decay=0.0 \
    --clip_norm=20.0 \
    --batch_type="$BATCH_TYPE" \
    --train_batch_size=$((BATCH_SIZE / NUM_NODES / NUM_GPU)) \
    --val_batch_size=$((BATCH_SIZE / NUM_NODES / NUM_GPU)) \
    --predict_batch_size=$((BATCH_SIZE / NUM_NODES / NUM_GPU)) \
    --accumulation_count="$ACCUMULATION_COUNT" \
  --num_cores="$NUM_CORES" \
    --beam_size=5 \
    --predict_min_len=1 \
    --predict_max_len=512 \
    --log_iter=100 \
    --eval_iter=5000 \
    --save_iter=10000 \
    --mask_rel_chirality=0 \
    --shared_attention_layer=0 \
    --compute_graph_distance \
    --resume \
    --load_from "$PWD/checkpoints/$DATA_NAME/model.60000_5.pt"

