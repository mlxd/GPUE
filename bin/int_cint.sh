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
	mkdir int counter_int
        cp $(pwd)/mpi_solver ./int;
        cp $(pwd)/mpi_solver ./counter_int;
        echo "Running I_m=0.0492 Counter-Intuitive @ $(date +'%Y-%m-%d_%H:%M:%S')" | tee -a Run.log;
        cd $(pwd)/counter_int;
        touch "START_$(date +'%Y-%m-%d_%H:%M:%S')";
        ./mpi_solver -m 0.0492 -l 0.1 -r 0.1 -o 1 >> result.log;
        touch "END_$(date +'%Y-%m-%d_%H:%M:%S')";
        cd ..
        echo "Finished I_m=0.0492 Counter-Intuitive @ $(date +'%Y-%m-%d_%H:%M:%S')" | tee -a Run.log;

	echo "Running I_m=0.0492 Intuitive @ $(date +'%Y-%m-%d_%H:%M:%S')" | tee -a Run.log;
	cd $(pwd)/int;
        touch "START_$(date +'%Y-%m-%d_%H:%M:%S')";
        ./mpi_solver -m 0.0492 -l 0.1 -r 0.1 -o -1 >> result.log;
        touch "END_$(date +'%Y-%m-%d_%H:%M:%S')";
        cd ..
        echo "Finished I_m=0.0492 Intuitive @ $(date +'%Y-%m-%d_%H:%M:%S')" | tee -a Run.log;
