if [[ -v NOOP ]]; then
  echo "Run NOOP mode"
  set -o noexec # same as set -n
fi

log_error () { 
  bash ${root}/utilities/error.sh "${BASH_COMMAND} failed with error code $?"
}
trap log_error ERR

