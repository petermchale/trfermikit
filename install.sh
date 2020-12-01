set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

########################## 

if ! which conda; then 
  bash utilities/error.sh "please install conda"
  exit 1
fi 

conda create --force --name trfermikit python=3.8 colored-traceback colorama cyvcf2

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

