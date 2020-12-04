set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

root=$1 

cd library 

if [[ ! -f "mysql-8.0.22-linux-glibc2.12-x86_64.tar" ]]; then 
  wget https://dev.mysql.com/get/Downloads/MySQL-8.0/mysql-8.0.22-linux-glibc2.12-x86_64.tar
  tar -xvf mysql-8.0.22-linux-glibc2.12-x86_64.tar
  tar -xvf mysql-8.0.22-linux-glibc2.12-x86_64.tar.xz
  cp mysql-8.0.22-linux-glibc2.12-x86_64/bin/mysql ${root}/bin
else 
  bash ${root}/utilities/info.sh "skipping installation of mysql"
fi 


