#! /usr/bin/env bash
# this script shrinks the 200-sample test set to 6 samples: the first 4
# clonal (Fws > 0.95) and the first 2 polyclonal (Fws < 0.8) samples from the
# Fws table of a full run on the 200-sample set (rand200.fws.txt)

set -ueExo pipefail

fws=${1:-rand200.fws.txt}

{
  awk '$2 > 0.95' $fws | head -4
  awk '$2 < 0.8' $fws | head -2
} | cut -f1 > 6_selected_samples.txt

bcftools view -S 6_selected_samples.txt all.passonly.subset_n200.vcf.gz -Ou \
  | bcftools view -i 'MAC>0' -Oz -o all.passonly.subset_n6.vcf.gz
