log_error () { 
  bash ${root}/utilities/error.sh "${BASH_COMMAND} failed with error code $?"
}
trap log_error ERR

