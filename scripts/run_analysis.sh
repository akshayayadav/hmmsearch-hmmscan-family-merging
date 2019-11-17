#!/bin/bash

mkdir /data/family_msa
mkdir /data/family_hmm
hmmsearch_hmmscan_predict_missing_sequences.py

awk '{print FILENAME"/"$0}' /data/results/*/*.predicted_missing_sequences|awk 'BEGIN{FS="/"}{print $4,$6}' >/data/combined-family.predicted_missing_sequences

get_family_mergings.py >family-mergings.txt
