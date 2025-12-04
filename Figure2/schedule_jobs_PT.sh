#!/usr/bin/env bash

cd /home/ubuntu/ebflow/Figure2

parallel -j 32 --joblog joblog_PT.txt --results output_PT :::: jobs_PT.txt
