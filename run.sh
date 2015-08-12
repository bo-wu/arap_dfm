#! /usr/bin/env bash
cd build
cmake ..
make
cd ..
#bin/morph data/389.off data/399.off
bin/morph data/382.off data/396.off
#bin/morph data/382.off data/382.off
#bin/morph data/399.off data/389.off
