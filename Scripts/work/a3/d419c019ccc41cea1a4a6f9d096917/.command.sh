#!/bin/bash -ue
fastp       -i sample_R1_1.fq.gz       -I sample_R1_2.fq.gz       -o sample_R1_R1.trim.fq.gz       -O sample_R1_R2.trim.fq.gz       -j sample_R1.fastp.json       -h sample_R1.fastp.html       --detect_adapter_for_pe       -w 4
