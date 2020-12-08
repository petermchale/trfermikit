root=$1 
key1=$2 
key2=$3 

fetch_value () { 
  echo $(${root}/bin/jq --raw-output --arg key1 $key1 --arg key2 $key2 '.[$key1][$key2]' ${root}/config.json)
}

