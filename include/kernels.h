/*** kernels.h - GPUE: Split Operator based GPU solver for Nonlinear 
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

#ifndef KERNELS_H
#define KERNELS_H
#include<stdio.h>
/* CUDA function declarations */

unsigned int getGid3d3d();;

/**
 * Indexing of threads, blocks
 */
__device__ unsigned int getBid3d3d();
__device__ unsigned int getTid3d3d();

/**
 * Helper functions for complex numbers
 */
__device__ double complexMagnitude(double2);
__device__ double complexMagnitudeSquared(double2);
__device__ double2 conjugate(double2 in);
__device__ double2 realCompMult(double scalar, double2 comp);

/**
 * Multiplication for linear, non-linear and phase-imprinting of the condensate.
 */
__global__ void cMult(cufftDoubleComplex*, cufftDoubleComplex*, cufftDoubleComplex*);
__global__ void cMultPhi(cufftDoubleComplex*, double*, cufftDoubleComplex*);
__global__ void cMultDensity(double2*, double2*, double2*, double, double,double, int, int);

/*
 * Hold vortex at specified position. Not fully implemented. cMultPhi should implement required functionality.
 */
__global__ void pinVortex(cufftDoubleComplex*, cufftDoubleComplex*, cufftDoubleComplex*);

/**
 * FFTW scaling and normalisation routines
 */
__global__ void scalarDiv(double2*, double, double2*);
__global__ void scalarDiv1D(double2*, double2*);
__global__ void scalarDiv2D(double2*, double2*);
__global__ void scalarDiv_wfcNorm(double2*, double, double2*, double2*);


/**
 * Parallel summation.
 */
__global__ void reduce(double2*, double*);
__global__ void multipass(cufftDoubleComplex*, cufftDoubleComplex*, int);

/**
 * Calculate angular momentum. Not fully implemented. Handled in post-processing instead.
 */
__global__ void angularOp(double, double, double2*, double*, double2*);


//####################################################################
/**
 * Non-implemented functions.
 */
__global__ void energyCalc(double2 *wfc, double2 *op, double dt, double2 *energy, int gnd_state, int op_space, double sqrt_omegaz_mass);
inline __device__ double2 braKetMult(double2 in1, double2 in2);
//template<typename T> __global__ void pSumT(T* in1, T* output, int pass);
__global__ void pSum(double* in1, double* output, int pass);
//template<double> __global__ void pSumT(double* in1, double* output, int pass);

#endif
