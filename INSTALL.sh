#!/bin/bash

cd $PWD

mkdir -p resources

# 1) Install MATLAB runtime 
# 1.1) dowload Matlab runtime (mcr)
wget https://users.ics.aalto.fi/~pemartti/fastGEAR/MCRInstallerLinux64bit.zip -P resources --no-check-certificate

# 1.2) unzip mcr and install
mkdir -p $PWD/resources/mcr
unzip resources/MCRInstallerLinux64bit.zip -d $PWD/resources/mcr
bash $PWD/resources/mcr/install -destinationFolder $PWD/resources/mcr -mode silent -agreeToLicense yes

# 2) Install fastGear excutable
# 2.1) download fastGear executable
wget https://users.ics.aalto.fi/~pemartti/fastGEAR/fastGEARpackageLinux64bit.tar.gz -P resources --no-check-certificate

# 2.2) unzip fastGear excutable 
tar -zvxf $PWD/resources/fastGEARpackageLinux64bit.tar.gz -C resources

# 3) create conda environment with all python dependencies installed
conda env create -f workflow/env/install.yaml -n BactPrep
