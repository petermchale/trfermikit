# find . -type f -regex ".*\.fq.gz" -ls | awk '{ print $NF }' | xargs rm
# find . -type f -regex ".*/fermikit.unsrt.sam.gz" -ls | awk '{ print $NF }' | xargs rm 
# find . -type f -regex ".*/fermikit.mag.gz" -ls | awk '{ print $NF }' | xargs rm
# find . -type f -regex ".*/fermikit.pre.gz" -ls | awk '{ print $NF }' | xargs rm
find . -type f -regex ".*/fermikit.flt.fmd" -ls | awk '{ print $NF }' | xargs rm
