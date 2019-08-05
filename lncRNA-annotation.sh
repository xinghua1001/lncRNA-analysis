#!/bin/bash
# find the closest protein coding genes for DE lncRNAs using bedtools
bedtools closest -a DE-lnc.bed -b protein_coding.bed -D ref > DE-lnc-closest.txt
