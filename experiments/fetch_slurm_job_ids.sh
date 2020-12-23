experiment=$1 

rm --force ${experiment}/data/job_ids.txt 

for log in ${experiment}/data/*/slurm*log; do 
  log=$(basename $log)
  log=${log#slurm.} 
  job_id=${log%.log}
  echo $job_id >> ${experiment}/data/job_ids.txt
done
