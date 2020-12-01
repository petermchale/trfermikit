set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

bin=$1 

cd library 

tool="bedtools"

if [[ ! -d ${tool} ]]; then 
  git clone https://github.com/arq5x/${tool}.git
  cd ${tool}
  git checkout 4a3615f5c534a41c22d2f047f6e9d32d162d45a4
  make
  cp bin/bedtools $bin
fi 
