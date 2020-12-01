set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

########################## 

kernel=$(uname --kernel-name)
machine=$(uname --machine)

if [[ ${machine} != 'x86_64' || ${kernel} != 'Linux'* ]]; then
    bash utilities/error.sh "not Linux x86_64"
    exit 1
fi 

########################## 

if ! which conda; then 
  bash utilities/error.sh "please install conda"
  exit 1
fi 

conda create --name trfermikit python=3.8 

set +o nounset
source activate trfermikit 
set -o nounset

pip install --requirement requirements.txt 

########################## 

mkdir --parents bin
mkdir --parents library 

root=$PWD

bash install/bedtools.sh ${root}

    # download samtools, minimap2, jq, bcftools, etc and install into $PWD/dependencies/bin
    # 	— (reference download_install_*.sh files

