#!/bin/bash

cd $PWD


# 1) make dir for software if not exist
mkdir -p resources

# 2) Install MATLAB runtime 
# 2.1) dowload Matlab runtime (mcr)
wget https://users.ics.aalto.fi/~pemartti/fastGEAR/MCRInstallerLinux64bit.zip -P resources --no-check-certificate

# 2.2) unzip mcr and install
mkdir -p $PWD/resources/mcr
unzip $PWD/resources/MCRInstallerLinux64bit.zip -d $PWD/resources/mcr
bash $PWD/resources/mcr/install -destinationFolder $PWD/resources/mcr -mode silent -agreeToLicense yes

# 3) Install fastGear excutable
# 3.1) download fastGear executable
wget https://users.ics.aalto.fi/~pemartti/fastGEAR/fastGEARpackageLinux64bit.tar.gz -P resources --no-check-certificate

# 3.2) unzip fastGear excutable 
tar -zvxf $PWD/resources/fastGEARpackageLinux64bit.tar.gz -C resources

# remove zip files from both programs
rm $PWD/resources/MCRInstallerLinux64bit.zip $PWD/resources/fastGEARpackageLinux64bit.tar.gz


# 4) add python file to path
chmod 755 $PWD/start_analysis.py

export PATH="$PWD:$PATH"
