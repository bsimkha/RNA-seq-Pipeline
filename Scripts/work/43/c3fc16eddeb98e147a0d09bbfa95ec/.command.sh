#!/bin/bash -ue
fastp       -i sample_R2_1.fq.gz       -I sample_R2_2.fq.gz       -o sample_R2_R1.trim.fq.gz       -O sample_R2_R2.trim.fq.gz       -j sample_R2.fastp.json       -h sample_R2.fastp.html       --detect_adapter_for_pe       -w 4
