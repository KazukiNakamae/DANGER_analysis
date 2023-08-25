#!/bin/sh
# GetRefs: download reference protein datasets
# specify directory for downloaded datasets
dir=db 
#
cd $dir
# zebrafish
curl -O https://ftp.ensembl.org/pub/current_fasta/danio_rerio/pep/Danio_rerio.GRCz11.pep.all.fa.gz