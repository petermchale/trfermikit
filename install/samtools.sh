set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

root=$1

cd library 

version="1.10"
tool="samtools"

if [[ ! -d "${tool}-${version}" ]]; then
  exit 1
  wget https://github.com/samtools/samtools/releases/download/${version}/${tool}-${version}.tar.bz2
  bzip2 -d ${tool}-${version}.tar.bz2
  tar -xvf ${tool}-${version}.tar
  rm ${tool}-${version}.tar
  cd ${tool}-${version}/
  ./configure --prefix=${root}
  make
  make install
else 
  bash ${root}/utilities/info.sh "skipping installation of ${tool}"
fi 





