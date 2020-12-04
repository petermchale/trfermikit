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

reference=$1 

[[ -f ${reference}.fa.fai ]] 
echo $? 
exit 1  
# || ${root}/bin/samtools faidx ${reference}.fa

# http://www.htslib.org/doc/faidx.html : 
[[ ! -f ${reference}.genome ]] && cut -f1,2 ${reference}.fa.fai > ${reference}.genome

# -x STR    Preset  []. This option applies multiple options at the same time. It should be applied before other options because options applied later will overwrite the values set by -x.  Available STR are:
# asm5    Long assembly to reference mapping (-k19 -w19 -A1 -B19 -O39,81 -E3,1 -s200 -z200 -N50 --min-occ-floor=100).  Typically, the  alignment  will  not extend to regions with 5% or higher sequence divergence. Only use this preset if the average divergence is far below 5%.
# asm10   Long assembly to reference mapping (-k19 -w19 -A1 -B9 -O16,41 -E2,1 -s200 -z200 -N50 --min-occ-floor=100).  Up to 10% sequence divergence.
# asm20   Long assembly to reference mapping (-k19 -w10 -A1 -B4 -O6,26 -E2,1 -s200 -z200 -N50 --min-occ-floor=100).  Up to 20% sequence divergence.
[[ ! -f ${reference}.mmi ]] && minimap2 -x asm10 -k23 -w11 -d ${reference}.mmi ${reference}.fa