#!/bin/bash
shopt -s extglob
while read -r line
do
  set -- $line
  for((i=1;i<=${#};i++))
  do
    s=$(eval echo \${${i}})
    case "$s" in
     +([0-9]).+([0-9]) ) s=$(printf "%.4f " $s);;
     +([-]*)+([0-9]).+([0-9]) ) s=$(printf "%.4f " $s);;
    esac
    printf "%s " $s
  done
  echo
done < "$1"
