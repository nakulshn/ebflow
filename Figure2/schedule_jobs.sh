#!/usr/bin/env bash

cd /home/ubuntu/ebflow/Figure2

parallel -j 16 :::: joblist_figure2_no_predict_bimodal.txt
