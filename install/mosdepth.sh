set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

root=$1
mosdepth="${root}/bin/jq"

if [[ ! -f ${mosdepth} ]]; then
  wget -O ${mosdepth} https://github.com/brentp/mosdepth/releases/download/v0.2.6/mosdepth 
  chmod +x ${mosdepth}
else 
  bash ${root}/utilities/info.sh "skipping installation of mosdepth"
fi 

bash ${root}/utilities/info.sh "mosdepth version is: $(${mosdepth} --version)"
