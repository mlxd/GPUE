#GPUE: Split Operator based GPU solver for Nonlinear 
#Schrodinger Equation, Copyright (C) 2011-2015, Lee J. O'Riordan 
#<loriordan@gmail.com>, Tadhg Morgan, Neil Crowley. All rights reserved.
#Redistribution and use in source and binary forms, with or without 
#modification, are permitted provided that the following conditions are 
#met:
#
#1. Redistributions of source code must retain the above copyright 
#notice, this list of conditions and the following disclaimer.
#
#2. Redistributions in binary form must reproduce the above copyright 
#notice, this list of conditions and the following disclaimer in the 
#documentation and/or other materials provided with the distribution.
#
#3. Neither the name of the copyright holder nor the names of its 
#contributors may be used to endorse or promote products derived from 
#this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
#PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
#HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
#SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
#TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
#PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
#LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
#NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#!/bin/bash
i=0
EMAIL=mymail@addr.com
count=0
NAME=$1
PARAMS=$2
declare -a JOBS=(-1 -1 -1 -1 -1 -1 -1 -1)
function run_gpue_test {
	echo $1
}

function run_gpue {
	if [[ $(echo $1 | head -c 1) == "#" ]];then
		return;
	elif [[ $(echo $1 | head -c 1) == "" ]];then 
		return
	fi
	
	if [ -n  "$NAME" ];then
		NAME=$(echo $NAME)_
	fi
	sleep 1
	A=$(date '+%y/%m/%d/%H_%M_%S')
	if [ -d ./$A ]; then
		echo "Exists"
		A=$A-$i
		i=$((i+1))
	fi
	echo "$NAME$A"
	mkdir -p $NAME$A
	cp ./gpue ./$NAME$A; cp -r ./src ./$NAME$A; cp -r ./include ./$NAME$A; cp ./Makefile ./$NAME$A; cp -r ./py ./$NAME$A; cp -r ./bin ./$NAME$A; cp ./wfc_load ./$NAME$A; cp ./wfci_load ./$NAME$A;
	cd ./$NAME$A
	pwd >> result.log
	echo $1 >>result.log
	mail -s "#Started GPU Job# $A" lee.oriordan@oist.jp < result.log
	./gpue $1 2>&1> result.log
	mkdir -p ./images
	#python ./py/vis.py >> result.log
	cp *.png ./images
	cd ./images
	ls | grep wfc_evr | grep _abs | grep png | sort -k3 -t _ -n > list1.txt;mencoder mf://@list1.txt -mf w=1280:h=1024:fps=24:type=png -oac copy -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:mv0:trell:v4mv:cbp:last_pred=3:predia=2:dia=2:vmax_b_frames=2:vb_strategy=1:precmp=2:cmp=2:subcmp=2:preme=2:qns=2:vbitrate=10000000 -o wfc_${PWD##*/}.avi
	rm -rf ./*.png
	#python ./py/hist3d.py
	rm wfc*
	mail -s "#Completed GPU Job# $A" $EMAIL < $(echo $(cat result.log; cat ./Params.dat))
	cd ../../../../..
	sleep 1
}

while read line ; do
	run_gpue "$line" &
	#echo "Running $line"
	JOBS[$count]=$!
	let count+=1
	
	if [ $count -gt 7 ]; then
		wait
		count=0
	fi
done < $PARAMS
