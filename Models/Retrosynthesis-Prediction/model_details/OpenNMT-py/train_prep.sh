#!/bin/bash

python preprocess.py -train_src data/Mechanism_Prediction/train/train_sources.txt \
-train_tgt data/Mechanism_Prediction/train/train_targets.txt -valid_src data/Mechanism_Prediction/validation/valid_sources.txt \
-valid_tgt data/Mechanism_Prediction/validation/valid_targets.txt -save_data data/Mechanism_Prediction/dataset \
-src_seq_length 1000 -tgt_seq_length 1000 -src_vocab_size 1000 -tgt_vocab_size 1000 -share_vocab