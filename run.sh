#! /usr/bin/env bash
cd build_arap
cmake ..
make
cd ..

#bin/morph model_data/382.off model_data/399.off ./model_data/382_399.acr
#bin/morph model_data/389.off model_data/399.off ./model_data/389_399.acr
#bin/morph data/382.off data/382.off
#bin/morph data/399.off data/389.off

#bin/morph ./model_data/aligned_354_plant/354.off ./model_data/aligned_354_plant/plant.OBJ
#bin/morph ./model_data/aligned_354_plant/354.off ./model_data/aligned_354_plant/plant4.obj

#bin/morph ./model_data/aligned_aladdin_deng/alading2b.OBJ ./model_data/aligned_aladdin_deng/dengshen5.OBJ
bin/morph ./model_data/201_208/201.obj ./model_data/201_208/208.obj
#bin/morph ./model_data/dino_397/Dino.off ./model_data/dino_397/397.obj
#bin/morph ./model_data/gecko_seal/gecko.obj ./model_data/gecko_seal/seal.obj
#bin/morph ./model_data/gorilla_horse/gorilla2.obj ./model_data/gorilla_horse/horse9.obj



















#################### just show resolution #################
##data_path='model_data/zhou_yang/'
##for off in `ls $data_path |grep -i obj`
##do
##    echo ${data_path}${off}
##    bin/morph ${data_path}${off} ${data_path}${off} ./model_data/382_399.acr
##done
