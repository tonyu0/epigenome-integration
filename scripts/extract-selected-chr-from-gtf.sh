#!/bin/bash

if [ $# -lt 2 ]; then
  echo "Usage: $0 <arg1> <arg2>\n<arg1> is a gtf file, <arg2> is a chromosome index you want)\n"
  exit 1
fi


awk -F'\t' -v chrIdx=$2 '$1 == chrIdx'  ${1} > outChr${2}.gtf
