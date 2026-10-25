#! /usr/bin/env bash
#
# Build bcftk into bin/bcftk. Run with: pixi run install-bcftk (in the
# `compile` environment, which provides cargo and libclang).
set -euo pipefail
if [ -x bin/bcftk ]; then echo "bcftk already available in bin/"; exit 0; fi
rm -rf bcftk
git clone git@github.com:bguo068/bcftk.git
cd bcftk
# bindgen (for htslib) needs libclang; take the environment's, not the system's.
LIBCLANG_PATH="$CONDA_PREFIX/lib" cargo build --release
cd ..
cp bcftk/target/release/bcftk bin/
rm -rf bcftk
