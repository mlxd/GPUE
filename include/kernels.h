/*
* kernels.h - GPUE: Split Operator based GPU solver for Nonlinear 
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

#ifndef KERNELS_H
#define KERNELS_H

/* CUDA function declarations */
__global__ void cMult(cufftDoubleComplex*, cufftDoubleComplex*, cufftDoubleComplex*);
__global__ void pinVortex(cufftDoubleComplex*, cufftDoubleComplex*, cufftDoubleComplex*);
__global__ void cMultDensity(double2*, double2*, double2*, double, double,double, int, int);
__global__ void scalarDiv(double2*, double, double2*);
__global__ void scalarDiv1D(double2*, double2*);
__global__ void scalarDiv2D(double2*, double2*);
__global__ void scalarDiv_wfcNorm(double2*, double, double2*, double2*);
__global__ void reduce(double2*, double*);
__device__ double complexMagnitude(double2);
__device__ double complexMagnitudeSquared(double2);
__global__ void multipass(cufftDoubleComplex*, cufftDoubleComplex*, int);

#endif
