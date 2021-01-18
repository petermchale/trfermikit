experiment=$1 

for log in ${experiment}/data/*/slurm*log; do 
  log=$(basename $log)
  log=${log#slurm.} 
  job_id=${log%.log}
  echo $job_id 
done
