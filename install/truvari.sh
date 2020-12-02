set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

root=$1

cd library 

version="2.0.2"

if [[ ! -d "truvari" ]]; then
  git clone https://github.com/spiralgenetics/truvari
  cd truvari
  python -m pip install --upgrade pip setuptools wheel
  python setup.py sdist bdist_wheel
  pip install dist/Truvari-${version}.tar.gz
else 
  bash ${root}/utilities/info.sh "skipping installation of truvari"
fi 





