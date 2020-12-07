# https://www.shell-tips.com/bash/debug-script/ 

# # TODO: eliminate all xtrace commands in other scripts in favor of the following 
# # (suitably debugged), 
# # which separates xtrace output from non-xtrace output that appears on STDERR
# if [[ -v TRACE ]]; then
#   bash ${root}/utilities/info.sh "Run TRACE mode"
#   exec 4>${output}/xtrace.out
#   BASH_XTRACEFD=4
#   set -o xtrace # same as set -x
# fi

log_error () { 
  bash ${root}/utilities/error.sh "${BASH_COMMAND} failed with error code $?"
}
trap log_error ERR

