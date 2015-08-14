#! /usr/bin/env bash
cd build_arap
cmake ..
make
cd ..
#bin/morph model_data/389.off model_data/399.off ./model_data/389_399.acr
bin/morph model_data/382.off model_data/399.off ./model_data/382_399.acr
#bin/morph data/382.off data/382.off
#bin/morph data/399.off data/389.off
