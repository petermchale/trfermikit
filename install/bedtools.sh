set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

bin=$1 

git clone https://github.com/arq5x/bedtools2.git
cd bedtools2
git checkout 4a3615f5c534a41c22d2f047f6e9d32d162d45a4
make
cp bin/bedtools $bin
