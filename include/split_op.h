/*** split_op.h - GPUE: Split Operator based GPU solver for Nonlinear 
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

#ifndef SPLIT_OP_H
#define SPLIT_OP_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <ctime>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <ctype.h>
#include <getopt.h>
#include "tracker.h"
#ifdef __linux
	#include<omp.h>
#elif __APPLE__
	//printf("OpenMP support disabled due to Clang/LLVM being behind the trend.",);
#endif

/* Keep track of all params for reading/writing to file*/
extern struct Params *paramS;

/* Error variable & return variables */
cudaError_t err;
cufftResult result;

/* Define operating modes */
int ang_mom = 0;
int gpe = 0;

/* Allocating global variables */
double mass, a_s, omegaX, omegaY, omegaZ;
double xi; //Healing length minimum value defined at central density.

/* Evolution timestep */
double dt, gdt;

/* Grid dimensions vector. xyz are dim length, w is total grid size (x*y*z) */
int xDim, yDim, read_wfc, print, write_it;
long  gsteps, esteps, atoms;
double *x,*y,*xp,*yp,*px,*py,dx,dy,xMax,yMax;

/* CuFFT plans for forward and inverse. May only need to use 1 for both */
cufftHandle plan_2d, plan_1d;

/* Arrays for storing wavefunction, momentum and position op, etc */
cufftDoubleComplex *wfc, *wfc0, *wfc_backup, *GK, *GV_half, *GV, *EK, *EV, *EV_opt, *GxPy, *GyPx, *ExPy, *EyPx, *EappliedField;
double *Energy, *Energy_gpu, *r, *Phi, *V, *V_opt, *K, *xPy, *yPx, *xPy_gpu, *yPx_gpu;

/* CUDA data buffers for FFT */
cufftDoubleComplex *wfc_gpu, *K_gpu, *V_gpu, *par_sum;
double *Phi_gpu;

/* CUDA streams */
cudaStream_t streamA, streamB, streamC, streamD;

/* Scaling the interaction */
double interaction;
double laser_power;

/* Define global dim3 and threads for grid and thread dim calculation */
dim3 grid;
int threads;

/* */
double l;
/* Function declarations */
/*
 * arg1 = Function result code from CUDA CUFFT calls.
 * arg2 = String data for name of function called. Prints value to stdout.
 */
int isError(int, char*); //Checks to see if an error has occurred. 

void parSum(double2* , double2* , int , int , int );
void optLatSetup(struct Vtx::Vortex centre, double* V, struct Vtx::Vortex *vArray, int num_vortices, double theta_opt, double intensity, double* v_opt, double *x, double *y);

double energy_angmom(double* Energy, double* Energy_gpu, double2 *V_op, double2 *K_op, double dx, double dy, double2 *gpuWfc, int gState);

#endif
