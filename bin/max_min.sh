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

for i in {3500..5500..25};
do
        if [ -d $(echo $(pwd)/0.0${i} | awk ' { sub("\\.*0+$","");print} ') ];
                then
			j=$(echo $(pwd)/0.0${i} | awk ' { sub("\\.*0+$","");print} ')
			cd $j
			cp ../max_min.py ./
			echo $j >> ../MAX_MIN.dat
			python ./max_min.py >> ../MAX_MIN.dat
			echo "" >> ../MAX_MIN.dat
			cd ..
        fi
done

