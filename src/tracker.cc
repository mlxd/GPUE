/*
* tracker.cc - GPUE: Split Operator based GPU solver for Nonlinear
* Schrodinger Equation, Copyright (C) 2012, Lee J. O'Riordan, Tadhg
* Morgan, Neil Crowley.

* This library is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2.1 of the
* License, or (at your option) any later version. This library is
* distributed in the hope that it will be useful, but WITHOUT ANY
* WARRANTY; without even the implied warranty of MERCHANTABILITY or
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
* License for more details. You should have received a copy of the GNU
* Lesser General Public License along with this library; if not, write
* to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
* Boston, MA 02111-1307 USA
*/

#include "../include/tracker.h"
#include "../include/fileIO.h"

void densityMin(double* density, long timestep, int xDim, int yDim, double* X, double* Y){
	int *index;
	int vortexMax = 100;
	int vortexIndex, numVortex;
	int x_vs[256], y_vs[256];
	double r, F_01, F_11, F_12, F_10, F_21;
	index = (int*) malloc(sizeof(int)*2*vortexMax);
	if(timestep==0){
		vortexIndex = 0;
		for(int i=0; i< xDim-1; ++i){
			for(int j=0; i< yDim-1; ++j){
				r = sqrt(X(i)*X(i) + Y(j)*Y(j));
				F_01 = density[(i-1)*yDim + j];
				F_11 = density[i*yDim + j];
				F_21 = density[(i+1)*yDim + j];
				F_10 = density[i*yDim + (j-1)];
				F_12 = density[i*yDim + (j+1)];
				if( (r<1e-4) && (F_11<F_10) && (F_11<F_12) && (F_11<F_01) && (F_11<F_21) ){
					vortexIndex += 1;
					x_vs[vortexIndex] = i;
					y_vs[vortexIndex] = j;
				}
			}
		}
		numVortex = vortexIndex;
	writeOutDouble(buffer, "v_x", x_vs, 256, timestep);		
	writeOutDouble(buffer, "v_y", y_vs, 256, timestep);		
	}
	return index;
}
