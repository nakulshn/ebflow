#!/usr/bin/env bash

cd /home/ubuntu/ebflow/Figure1

parallel -j 32 --joblog joblog.txt --results output :::: jobs_nonpt.txt
