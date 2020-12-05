# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --genome-build ) shift; [[ ! $1 =~ ^- ]] && genome_build=$1;;
    --root ) shift; [[ ! $1 =~ ^- ]] && root=$1;;
    --repeats ) shift; [[ ! $1 =~ ^- ]] && repeats=$1;;
    *) bash ${root}/utilities/error.sh "$0: $1 is an invalid flag"; exit 1;;
  esac 
  shift
done

set -o errexit
set -o pipefail
set -o nounset
# set -o noclobber

set -o xtrace
# Must use single quote to prevent variable expansion.
# For example, if double quotes were used, ${LINENO} would take on the value of the current line,
# instead of its value when PS4 is used later in the script
# https://stackoverflow.com/a/6697845/6674256
# ${FOO:+val}    val if $FOO is set
# ${FOO[0]}   element #0 of the FOO array
# https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

table="simpleRepeat"
database="${genome_build}"

chromosome="chrom"
start_coordinate="chromStart"
end_coordinate="chromEnd"
period="period"

# https://genome.ucsc.edu/goldenPath/help/mysql.html
# https://genome.ucsc.edu/cgi-bin/hgTables
${root}/bin/mysql \
    --user=genome \
    --host=genome-mysql.soe.ucsc.edu \
    --port=3306 \
    --column-names \
    --batch \
    --no-auto-rehash \
    --execute="SELECT ${chromosome}, ${start_coordinate}, ${end_coordinate}, ${period} from ${table};" \
    ${database} \
  | ${root}/bin/bgzip --stdout \
  > ${repeats}.tab.gz

