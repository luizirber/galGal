#!/bin/bash

for f in $@
do
  while read line
  do
    if [[ ${line:0:1} == '>' ]]
    then
      echo ">"$(basename $f| cut -d "_" -f2)"|"${line:1:200}
    else
      echo ${line^^}
    fi
  done < $f
done
