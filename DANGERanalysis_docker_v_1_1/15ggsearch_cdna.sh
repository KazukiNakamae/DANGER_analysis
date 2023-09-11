#!/bin/sh
# exec ggsearch for global alignment
# name of query
query=$1
# name of db
db=$2
# threads to use
thre=8
# cutoff threshold (E-value)
evalue=0.1
# number of alignments
d=1
# name of program
ggsearch=ggsearch36
# outfile
out=ggsearch_${query}-${db}.txt
gzip=pigz
# run ggsearch
for file in `\find ./query_files -maxdepth 1 -type f`; do
  echo $file
  ggsearch36 -Q -T $thre -d $d -m10 -E $evalue $file $db >> $out
done
pigz $out;
