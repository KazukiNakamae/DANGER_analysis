#!/bin/sh
# batch shell script to run perl script for parsing ggsearch36 results
# Hidemasa Bono

perl=perl # perl command
# zcat=gzcat # command to gzcat
target=Trinity.fasta # name of target organism
pscript=15parseggsearch_cdna.pl # path to the parser script

for i in $@; do
	echo "parsing $i"
	zcat ggsearch_${target}-$i.txt.gz \
	| $perl $pscript > ${target}-$i.txt
done
