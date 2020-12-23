#!/bin/bash

# input login information
login_info=$1 #specify full path
outdir=$2

# #test:
# login_info=/data/Athaliana_RNAseq/login.info
# outdir=/data/Athaliana_RNAseq/

echo "[Credentials]" > ~/.ossutilconfig
echo "language=EN" >> ~/.ossutilconfig
echo "accessKeyID=$(cut -d, -f2 $login_info | head -n1 | tail -n1)" >> ~/.ossutilconfig
echo "accessKeySecret=$(cut -d, -f2 $login_info | head -n2 | tail -n1)" >> ~/.ossutilconfig
OSS_LOC=$(cut -d, -f2 $login_info | head -n3 | tail -n1)
echo "endpoint=$(cut -d, -f2 $login_info | head -n4 | tail -n1).aliyuncs.com" >> ~/.ossutilconfig

/data/Lolium/Softwares/ossutil64 ls ${OSS_LOC}

# or or just download the whole folder
cd $outdir
time /data/Lolium/Softwares/ossutil64 cp -r ${OSS_LOC} .
