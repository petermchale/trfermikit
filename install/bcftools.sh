set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

root=$1

cd library 

version="1.10.2"

if [[ ! -d "bcftools-${version}" ]]; then
  # http://www.htslib.org/download/
  wget https://github.com/samtools/bcftools/releases/download/${version}/bcftools-${version}.tar.bz2
  bzip2 -d bcftools-${version}.tar.bz2
  tar -xvf bcftools-${version}.tar
  rm bcftools-${version}.tar
  cd bcftools-${version}
  ./configure --prefix=${root}
  make
  make install
else 
  bash ${root}/utilities/info.sh "skipping installation of bcftools"
fi 




