#!/bin/bash
while read line; do
A=$(date '+%y/%m/%d/%H_%M_%S')
echo $A
mkdir -p $A
cp ./gpue ./$A
cp -r ./src ./$A
cp -r ./include ./$A
cp ./Makefile ./$A
cp -r ./py ./$A
cp -r ./bin ./$A
cd ./$A
./gpue $line > result.log
mkdir -p ./images
python ./py/vis.py
cp *.png ./images
cd ./images
ls -tr | grep wfc | grep png > list1.txt
ls -tr | grep cntr | grep png > list2.txt
ls -tr | grep phi | grep png > list3.txt
mencoder mf://@list1.txt -mf w=1280:h=1024:fps=10:type=png -oac copy -ovc copy -o wfc.avi
mencoder mf://@list2.txt -mf w=1280:h=1024:fps=10:type=png -oac copy -ovc copy -o cntr.avi
mencoder mf://@list3.txt -mf w=1280:h=1024:fps=10:type=png -oac copy -ovc copy -o phi.avi
cd ../../../../..
done < ./bin/run_params.conf
