#!/usr/bin/env bash

cd /home/ubuntu/ebflow/Figure1

parallel -j 32 --joblog joblog_pt.txt --results output_pt :::: jobs_pt.txt
