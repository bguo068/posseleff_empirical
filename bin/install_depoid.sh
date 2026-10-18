#! /usr/bin/env bash
# 
# run this script with
# pixi run -e deploid bash install_depoid.sh
rm -rf DEploid
git clone --recurse-submodules https://github.com/DEploid-dev/DEploid.git
cd DEploid
git checkout v0.6-beta
git submodule update --init --recursive

bash bootstrap

# "version" file is picked up as a c++ header file.
mv version{,.txt}
mv src/gzstream/version{,.txt}


if [ "$(uname)" = "Darwin" ]; then
    echo "Running on macOS"
    # add the flag to fix INIFINITY not defined error
    make CXXFLAGS='-DINFINITY=__builtin_inff\(\)' -j 10
    # make -j 10
else
    make  -j 10
fi

cp dEploid ../bin/
cd ..
rm -rf dEploid
