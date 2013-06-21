#
# int_cint.sh - GPUE: Split Operator based GPU solver for Nonlinear 
# Schrodinger Equation, Copyright (C) 2012, Lee J. O'Riordan, Tadhg 
# Morgan, Neil Crowley. 

# This library is free software; you can redistribute it and/or modify 
# it under the terms of the GNU Lesser General Public License as 
# published by the Free Software Foundation; either version 2.1 of the 
# License, or (at your option) any later version. This library is 
# distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public 
# License for more details. You should have received a copy of the GNU 
# Lesser General Public License along with this library; if not, write 
# to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
# Boston, MA 02111-1307 USA 
#

#!/bin/bash
mkdir plots
for i in {2500..6500..25};
do
	j="$(echo $i|awk '{sub("\\.*0+$","");print}')"
	if [ -d $(pwd)/0.0$j ]
	then
		echo "Plotting I_m=0.0$j" >> Plot.log
		cd ./0.0$j
		echo $(pwd)
		ls ./p0
		echo "0.0$j,$(tail -n 1 p0)" >> ./../master_p0.dat
		echo "0.0$j,$(tail -n 1 p1)" >> ./../master_p1.dat
		echo "0.0$j,$(tail -n 1 p2)" >> ./../master_p2.dat
		cp ../population.gp ./
		gnuplot ../population.gp
		mv ./population.png ../plots/0.0$j.png
		cd ./..
	fi
done
