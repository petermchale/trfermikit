# TODO: move all xtrace commands in this repo into this file
# https://www.shell-tips.com/bash/debug-script/ 

log_error () { 
  bash ${root}/utilities/error.sh "${BASH_COMMAND} failed with error code $?"
}
trap log_error ERR

