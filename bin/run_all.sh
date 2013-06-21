#
# run_all.sh - GPUE: Split Operator based GPU solver for Nonlinear 
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
for i in {40..55..1};
do
        if [ ! -d $(pwd)/0.0$i ];
                then
                        echo "Folder does not exist. Creating $(pwd)/0.0$i" | tee -a Run.log;
                        mkdir $(pwd)/0.0$i;
        fi;
        cp $(pwd)/mpi_solver $(pwd)/0.0$i;
        echo "Running I_m=0.0$i @ $(date +'%Y-%m-%d_%H:%M:%S')" | tee -a Run.log;
        cd $(pwd)/0.0$i;
        touch "START_$(date +'%Y-%m-%d_%H:%M:%S')";
        ./mpi_solver -m 0.0$i -l 0.1 -r 0.1 >> result.log;
        touch "END_$(date +'%Y-%m-%d_%H:%M:%S')";
        cd ..
        echo "Finished I_m=0.0$i @ $(date +'%Y-%m-%d_%H:%M:%S')" | tee -a Run.log;
done
echo "Run finished at $(date +'%Y-%m-%d_%H:%M:%S')" | tee -a Run.log
