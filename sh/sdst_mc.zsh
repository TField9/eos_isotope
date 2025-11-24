#!/bin/zsh

# Usage:
# sdst.zsh [filelist] [ibeg] [iend] [ofn]
# example:
# ./sdst.zsh /cvmfs/ams.cern.ch/Offline/AMSDataDir/DataManagement/DataSetsDesc/Be.B1403/be10.pl1.l1.4400.6_05.txt 1 10 sdst_1_10.root

# 设置 Kerberos 认证
# 设置 Kerberos 环境

export KRB5_CONFIG=/etc/krb5.conf
export KRB5CCNAME=/tmp/krb5cc_condor_${PPID}
KEYTAB="/eos/ams/user/s/selu/mdst/tianye/seludev.keytab"

echo "Trying kinit with keytab: $KEYTAB"
/usr/bin/kinit -k -t $KEYTAB seludev@CERN.CH

# 检查 kinit 是否成功
if [ $? -ne 0 ]; then
    echo "kinit failed! Detailed error:"
    # 尝试获取更详细的错误信息
    /usr/bin/kinit -k -t $KEYTAB seludev@CERN.CH -V
    exit 1
fi

echo "kinit successful. Current ticket:"
klist

# # 测试 EOS 访问
# echo "Testing EOS access..."
# eos root://eosams.cern.ch ls /eos/ams
# if [ $? -ne 0 ]; then
#     echo "EOS access test failed!"
#     klist
#     exit 1
# fi

flist=$1
ibeg=$2
iend=$3
ofn=$4

export Offline=/cvmfs/ams.cern.ch/Offline
export AMSWD=$Offline/vdev #Official 
#export AMSWD=/work/alpha/AMS


export NOCASTOR=1
export PGTRACK=1
export _PGTRACK_=1
export NORFIO=1
export GLIBCXX_USE_CXX11_ABI=0
export NO_RFIOD=1
export NO_NAG=1
export NORFIOD=1
export EOS_MGM_URL=root://eosams.cern.ch
export AMSOPT=/cvmfs/ams.cern.ch/opt
export AMSSRC=$AMSWD
export AMSDataDir=$Offline/AMSDataDir
export UCC=icc64
export UICC=2024
export UROOT=5
export USLC=9

export XRDLIB="/cvmfs/ams.cern.ch/opt/xrootd-4.8.2.el9"
#export ROOTSYS="/cvmfs/ams.cern.ch/Offline/root/Linux/root-v5-34-9-icc64.24-el9"
export ROOTSYS="/cvmfs/ams.cern.ch/Offline/root/Linux/root6-14-04-icc64.24-el9"
export LD_LIBRARY_PATH=.:${ROOTSYS}/lib:${XRDLIB}/lib64
export PATH=${XRDLIB}/bin:$ROOTSYS/bin:.:$PATH:/usr/sbin

pushd ${ROOTSYS}
source $ROOTSYS/bin/thisroot.sh
popd
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/cvmfs/ams.cern.ch/opt/lib64_EL9

export INTEL_LICENSE_FILE=$Offline/intel/licenses
source /cvmfs/projects.cern.ch/intelsw/oneAPI/linux/x86_64/2024/compiler/latest/env/vars.sh

export VERBOSE=1
export GENFIT=1
export DPMJET3=6
export G4MULTITHTREADED=1
export AMSICC=1
export AMSP=1
export G4AMS=1
export CXX=icc
export CC=icc
export FC=ifort
export AMSGeoDir=$Offline/vdev/display/
export amsedcPG=$Offline/vdev/exe/linxx8664icc5.34/amsedcPG

#source $AMSWD/install/g4i10.1.sh
export CVSROOT=/afs/cern.ch/exp/ams/Offline/CVS
export INTELSW=/cvmfs/projects.cern.ch
export INTELDIR=$INTELSW/intelsw/oneAPI/linux/x86_64/2024
export INTELVER=compiler/latest




# Execute the sdst binary with output filename
eos cp /afs/cern.ch/work/s/seludev/private/tianye/sh/param_ML_test.txt ./
eos cp /afs/cern.ch/work/s/seludev/private/tianye/bin/sdst ./
chmod +x sdst
./sdst $flist $ibeg $iend $ofn

sleep 1
sync 
sleep 1
odir="/eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_train"
echo to cp ${ofn} ${odir}

xrdcp --retry 3 ${ofn} root://eosams.cern.ch/${odir}/${ofn}