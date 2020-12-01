set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

########################## 

if ! which conda; then 
  echo "please install conda"
fi 

tool="trfermikit"
if conda activate $tool; then 
  echo "conda environment activated"
else 
  conda create --name $tool python=3.8
  conda activate trfermikit
  echo "conda environment created and activated"
fi

pip install --requirement requirements.txt 

########################## 

kernel=$(uname --kernel-name)
machine=$(uname --machine)

if [[ ${machine} == 'x86_64' ]]; then
  if [[ ${kernel} == 'Linux'* ]]; then
    # download bedtools, samtools, minimap2, jq, bcftools, etc and install into $PWD/dependencies/bin
    # 	— (reference download_install_*.sh files
    echo "download bedtools, samtools, minimap2, jq, bcftools, etc and install into $PWD/dependencies/bin"
  else
    echo 'not Linux'
    exit 1
  fi
else
  echo 'not x86_64'
  exit 1
fi

