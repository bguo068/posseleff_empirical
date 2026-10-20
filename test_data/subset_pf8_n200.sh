#! /usr/bin/env bash
# this script extract genotype of 200 randomly selected samples
# from MalariaGEN Pf8 

set -ueExo pipefail


for i in {1..14}; do
  chrno=`printf "%02d" $i`
  samples=200_rand_good_samples.txt
  infile=../Pf3D7_${chrno}_v3.passonly.bcf
  outfile=Pf3D7_${chrno}_v3.passonly.subset_n200.bcf
  bcftools view -S $samples -m2 -M2 -v snps --threads 5 $infile -Ou \
    | bcftools view -q 0.01:minor -Ou \
    | bcftools annotate -x ^FMT/GT,FMT/AD --threads 5  -Ob -o $outfile
done
