set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

root=$1 

cd library 

if [[ ! -d "minimap2" ]]; then 
  git clone https://github.com/lh3/minimap2
  cd minimap2 
  git checkout cdb7857841fa472f26c7575ba3c51682af5c885a
  make
  cp minimap2 ${root}/bin
else 
  bash ${root}/utilities/info.sh "skipping installation of minimap2"
fi 






