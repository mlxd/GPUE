///@cond LICENSE
/*** manip.cu - GPUE: Split Operator based GPU solver for Nonlinear
Schrodinger Equation, Copyright (C) 2011-2015, Lee J. O'Riordan
<loriordan@gmail.com>, Tadhg Morgan, Neil Crowley.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
///@endcond
//##############################################################################
/**
 *  @file    manip.cu
 *  @author  Lee J. O'Riordan (mlxd)
 *  @date    11/08/2015
 *  @version 0.1
 */
 //##############################################################################
#include "../include/manip.h"

void WFC::phaseWinding(double *phi, int winding, double *x, double *y, double dx, double dy, double posx, double posy, int dim){
	for(int ii=0; ii < dim; ++ii){
		for(int jj=0; jj < dim; ++jj){
			phi[ii*dim +jj] = fmod(winding*atan2(y[jj]-(posy-dim/2+1)*dx, x[ii]-(posx-dim/2+1)*dy),2*PI);
		}
	}
}

void WFC::phaseWinding(double *phi, int winding, double *x, double *y, double dx, double dy, double* posx, double* posy, int sites, int dim){
	memset(phi,0,dim*dim);
	for(int a = 0; a< sites; ++a){
		for(int ii=0; ii < dim; ++ii){
			for(int jj=0; jj < dim; ++jj){
				phi[ii*dim +jj] = fmod(phi[ii*dim +jj] + winding*(atan2(y[jj]-(posy[a]-dim/2+1)*dx, x[ii]-(posx[a]-dim/2+1)*dy) ),2*PI);
			}
		}
	}
}

void WFC::applyPhase(double *phi, double2 *wfc, int dim){
	for(int ii=dim-1; ii >= 0; --ii){
		for(int jj=dim-1; jj >= 0; --jj){
			wfc[ii*dim + jj].x *= cos(phi[ii*dim + jj]);
			wfc[ii*dim + jj].y *= -sin(phi[ii*dim + jj]);
		}
	}
}
