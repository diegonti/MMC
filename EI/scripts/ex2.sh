#!/bin/bash

for file in ../test/*.txt
    do
        #head $file$ 
        FIRST=$(head -1 $file) 
        LAST=$(tail -1 $file)
        echo $FIRST, $LAST
        #echo -e "$LAST" "\n$(cat $file | tail -2 | head -n -1)" $"\n$first" > $file
        #sed "1s/*/${LAST}/" $file | sed "$s/*/${FIRST}/"
done
