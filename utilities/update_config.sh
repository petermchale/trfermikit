#!/usr/bin/env bash

set -o xtrace 
set -o errexit 

root=$1 
output=$2
key1=$3
key2=$4 
value=$5

${root}/bin/jq \
    --arg key1 $key1 \
    --arg key2 $key2 \
    --arg value $value \
    '.[$key1][$key2] = $value' \
    $output/config.json 
#   > $output/config.tmp.json
# mv $output/config.tmp.json $output/config.json
