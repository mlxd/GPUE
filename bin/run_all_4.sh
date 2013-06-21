#
# run_all_4.sh - GPUE: Split Operator based GPU solver for Nonlinear 
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
echo "Run started at $(date +'%Y-%m-%d_%H:%M:%S')" | tee Run.log

for i in {5500..5600..25};
do
        if [ ! -d $(echo $(pwd)/0.0${i} | awk ' { sub("\\.*0+$","");print} ') ];
                then
                        echo "Folder does not exist. Creating $(pwd)/0.0$i" | tee -a Run.log;
                        mkdir $(pwd)/0.0$i;
        		cp $(pwd)/mpi_solver $(pwd)/0.0$i;
        		echo "Running I_m=0.0$i @ $(date +'%Y-%m-%d_%H:%M:%S')" | tee -a Run.log;
		        cd $(pwd)/0.0$i;
		        touch "START_$(date +'%Y-%m-%d_%H:%M:%S')";
		        ./mpi_solver -m 0.0$i -d 0 -f 4.0 -o 1 >> result.log;
		        touch "END_$(date +'%Y-%m-%d_%H:%M:%S')";
		        cd ..
		        echo "Finished I_m=0.0$i @ $(date +'%Y-%m-%d_%H:%M:%S')" | tee -a Run.log;
		else
			echo "Folder exists, skipping 0.0$i" | tee -a Run.log;
	fi
done
echo "Run finished at $(date +'%Y-%m-%d_%H:%M:%S')" | tee -a Run.log
