#!/bin/bash

python translate.py -model Mechanism_Prediction_Model/Mechanism_Prediction_Model_step_100000.pt \
-src data/Mechanism_Prediction/test/test_sources.txt -output Model_output/output.txt \
-batch_size 128 -replace_unk -max_length 200 -beam_size 1 -n_best 1 -gpu 0