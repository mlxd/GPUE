#!/bin/bash
for i in $(ls);
do
	if [ -d $i ];
	then
		mv ./$i $(echo $i|awk '{sub("\\.*0+$","");print}')
	fi
done

