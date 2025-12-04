#!/usr/bin/env bash

cd /home/ubuntu/ebflow/Figure2

parallel -j 32 --joblog joblog.txt --results output :::: joblist_figure2_nonPT_complete.txt
