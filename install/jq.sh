set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

root=$1
jq="${root}/bin/jq"

if [[ ! -f ${jq}  ]]; then
  wget -O ${jq} https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64
  chmod +x ${jq}
else
  bash ${root}/utilities/info.sh "skipping installation of jq"
fi 

bash ${root}/utilities/info.sh "jq version is: $(${jq} --version)"


