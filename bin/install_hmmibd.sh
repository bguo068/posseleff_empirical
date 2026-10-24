#! /usr/bin/env bash
#
# Build hmmIBD (with this pipeline's recombination rate) into bin/hmmIBD.
# Run with: pixi run install-hmmibd
set -euo pipefail
if [ -x bin/hmmIBD ]; then echo "hmmIBD already available in bin/"; exit 0; fi
rm -rf hmmIBD
git clone https://github.com/glipsnort/hmmIBD.git
cd hmmIBD
# a specific version for reproducibility
git checkout a2f796ef8122d7f6b983ae9ac4c6fba35afcd3aa
sed -i.bak -e 's/const double rec_rate = 7.4e-7/const double rec_rate = 6.67e-7/' hmmIBD.c
${CC:-cc} -o hmmIBD -O3 -Wall hmmIBD.c -lm
cd ..
cp hmmIBD/hmmIBD bin/hmmIBD
rm -rf hmmIBD
