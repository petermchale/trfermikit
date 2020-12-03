set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

root=$1

cd library 

version="1.9"

if [[ ! -d "htslib-${version}" ]]; then
  wget https://github.com/samtools/htslib/releases/download/${version}/htslib-${version}.tar.bz2
  bzip2 -d htslib-${version}.tar.bz2
  tar -xvf htslib-${version}.tar
  rm htslib-${version}.tar
  cd htslib-${version}/
  ./configure --prefix=${root}
  make
  make install
else 
  bash ${root}/utilities/info.sh "skipping installation of htslib"
fi 
