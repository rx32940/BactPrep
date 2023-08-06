#!/bin/bash

cd $PWD

# 1) add python file to path
chmod 755 $PWD/start_analysis.py

ln -s $PWD/start_analysis.py "$(echo $PATH | tr ':' '\n' | grep -m 1 conda)/start_analysis.py"

# 2) make dir for software if not exist
mkdir -p resources

# 3) Install MATLAB runtime 
# 3.1) dowload Matlab runtime (mcr)
wget https://users.ics.aalto.fi/~pemartti/fastGEAR/MCRInstallerLinux64bit.zip -P resources --no-check-certificate

# 3.2) unzip mcr and install
mkdir -p $PWD/resources/mcr
unzip $PWD/resources/MCRInstallerLinux64bit.zip -d $PWD/resources/mcr
bash $PWD/resources/mcr/install -destinationFolder $PWD/resources/mcr -mode silent -agreeToLicense yes

# 4) Install fastGear excutable
# 4.1) download fastGear executable
wget https://users.ics.aalto.fi/~pemartti/fastGEAR/fastGEARpackageLinux64bit.tar.gz -P resources --no-check-certificate

# 4.2) unzip fastGear excutable 
tar -zvxf $PWD/resources/fastGEARpackageLinux64bit.tar.gz -C resources

# remove zip files from both programs
rm $PWD/resources/MCRInstallerLinux64bit.zip $PWD/resources/fastGEARpackageLinux64bit.tar.gz



