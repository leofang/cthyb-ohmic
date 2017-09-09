#!/bin/bash

# prepare Singularity definition
cp cthyb_ohmic.Singularity cthyb_ohmic_temp_Singularity
echo "Build_date  $(date +"%b. %d, %T %Z, %Y")" >> cthyb_ohmic_temp_Singularity

# create the image
singularity create -s 2048 cthyb_ohmic.img
sudo singularity bootstrap cthyb_ohmic.img cthyb_ohmic_temp_Singularity

# clean up
rm cthyb_ohmic_temp_Singularity
