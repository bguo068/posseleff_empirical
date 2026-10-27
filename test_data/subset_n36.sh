#! /usr/bin/env bash
# this script shrinks the 200-sample test set to 36 samples, picked from the
# Fws table of a full run on the 200-sample set (rand200.fws.txt):
# - the first 30 clonal samples (Fws > 0.95)
# - 4 polyclonal samples with Fws 0.93-0.94, likely to carry a dominant clone
# - the first 2 polyclonal samples with Fws < 0.8 (no dominant clone)

set -ueExo pipefail

fws=${1:-rand200.fws.txt}

{
  awk '$2 > 0.95' $fws | head -30
  for s in RCN13096 QP0060-C SPT26932 QC0558-CW; do
    awk -v s=$s '$1 == s' $fws
  done
  awk '$2 < 0.8' $fws | head -2
} | cut -f1 > 36_selected_samples.txt

bcftools view -S 36_selected_samples.txt all.passonly.subset_n200.vcf.gz -Ou \
  | bcftools view -i 'MAC>0' -Oz -o all.passonly.subset_n36.vcf.gz
