#!/bin/bash

# Usage: ./merge_tsv.sh in1.tsv in2.tsv [in3.tsv ...]

# Merge tsv files
# Assumes headers start with a hash
# Retains only header from first file

# check header length
i=0
while [ true ]; do
  if [ $(head -n $i $1 | grep -c -e "^#") -lt $i ]; then
    i=$(($i-1))
    break
  fi
  i=$(($i+1))
done

cat <(head -n1 $1) <(for i in $@; do tail -n +2 $i; done)
