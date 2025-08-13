#!/bin/bash

# Input your checkpoint path
export CHECKPOINT=""
export MODEL=graph2smiles
export BATCH_TYPE=tokens
export BATCH_SIZE=4096

python beam_search.py \
  --do_beam_search \
  --do_score \
  --start_step=1 \
  --end_step=3 \
  --data_name="$DATA_NAME" \
  --model="$MODEL" \
  --load_from="$CHECKPOINT" \
  --log_file="graph2smiles_beam_$DATA_NAME" \
  --processed_data_path="$PROCESSED_DATA_PATH" \
  --model_path="$MODEL_PATH" \
  --test_output_path="$TEST_OUTPUT_PATH" \
  --batch_type="$BATCH_TYPE" \
  --predict_batch_size="$BATCH_SIZE" \
  --accumulation_count=4 \
  --num_cores="$NUM_CORES" \
  --beam_size=30 \
  --n_best=3 \
  --batch_search_size=3000 \
  --do_batch_search=True \
  --topK=3 \
  --predict_min_len=1 \
  --predict_max_len=512 \
  --log_iter=100