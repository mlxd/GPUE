/*
* split_op.h - GPUE: Split Operator based GPU solver for Nonlinear 
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

#ifndef SPLIT_OP_H
#define SPLIT_OP_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <ctype.h>
#include <getopt.h>

/* Error variable & return variables */
cudaError_t err;
cufftResult result;

/* Define operating modes */
int ang_mom = 0;
int gpe = 0;
double mass, omegaX, omegaY, omegaZ;

/* Evolution timestep */
double dt;

/* Grid dimensions vector. xyz are dim length, w is total grid size (x*y*z) */
int xDim, yDim, read_wfc, print;
long  gsteps, esteps, atoms;
double dx,dy;

/* CuFFT plans for forward and inverse. May only need to use 1 for both */
cufftHandle plan_2d, plan_1d;

/* Arrays for storing wavefunction, momentum and position op, etc */
cufftDoubleComplex *wfc, *wfc_backup, *GK, *GV, *EK, *EV, *xPy, *yPx, *GxPy, *GyPx, *ExPy, *EyPx, *EappliedField;
double *Phi, *V, *K;
double2 *r, *r_gpu; 

/* CUDA data buffers for FFT */
cufftDoubleComplex *wfc_gpu, *K_gpu, *V_gpu, *xPy_gpu, *yPx_gpu, *par_sum;


double interaction;

/* Define global dim3 and threads for grid and thread dim calculation */
dim3 grid;
int threads;

/* Function declarations */
/*
 * arg1 = Function result code from CUDA CUFFT calls.
 * arg2 = String data for name of function called. Prints value to stdout.
 */
int isError(int, char*); //Checks to see if an error has occurred. 
void writeOut(char*, char*, double2*, int, int); //Writes out to file

void parSum(double2* , double2* , double2* , int , int , int );

#endif
