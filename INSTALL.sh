#!/bin/bash

cd $PWD

wget https://users.ics.aalto.fi/~pemartti/fastGEAR/MCRInstallerLinux64bit.zip -P resources --no-check-certificate

mkdir $PWD/resources/mcr
unzip resources/MCRInstallerLinux64bit.zip -d $PWD/resources/mcr
bash $PWD/resources/mcr/install -destinationFolder $PWD/resources/mcr -mode silent -agreeToLicense yes

tar -zvxf $PWD/resources/fastGEARpackageLinux64bit.tar.gz -C resources

conda env create -f workflow/env/install.yaml -n BactPrep
