#!/bin/bash

for file in *.py
do
	forCount=$(grep "for" $file | wc -l)
	echo "The file $file has $forCount for lines"

	if [ $forCount -ge 10 ]; then
		grep "for" $file > $file.forY.txt
	elif [ $forCount -lt 10 ]; then
		grep "for" $file > $file.forN.txt
	fi
done
