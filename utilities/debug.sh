# shellcheck shell=bash

# https://www.shell-tips.com/bash/debug-script/

echo "root is: $root"
exit 1 

log_error () { 
  bash ${root}/utilities/error.sh "${BASH_COMMAND} failed with error code $?"
}
trap log_error ERR

