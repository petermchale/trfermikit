#!/usr/bin/env bash

set -o xtrace 

root=$1 
output=$2
key1=$3
key2=$4 

echo $(${root}/bin/jq --raw-output --arg key1 $key1 --arg key2 $key2 '.[$key1][$key2]' ${output}/config.json)

