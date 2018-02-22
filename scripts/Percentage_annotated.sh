#!/bin/bash

for f in ../../refseq_entries/*.gb
do
    perl Percentage_Annotated.pl $f
done
