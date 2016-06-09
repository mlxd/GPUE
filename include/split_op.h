///@cond LICENSE
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
///@endcond
//##############################################################################
/**
 *  @file    split_op.h
 *  @author  Lee J. O'Riordan (mlxd)
 *  @date    12/11/2015
 *  @version 0.1
 *
 *  @brief Host and device declarations for simulations.
 *
 *  @section DESCRIPTION
 *  These functions and variables are necessary for carrying out the GPUE
 *	simulations. This file will be re-written in an improved form in some
 *	future release.
 */
//##############################################################################

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
#include "ds.h"
#ifdef __linux
	#include<omp.h>
#elif __APPLE__
	//printf("OpenMP support disabled due to Clang/LLVM being behind the trend.",);
#endif

/* Error variable & return variables */
//extern cufftResult result;

/* CuFFT plans for forward and inverse. May only need to use 1 for both */
//extern cufftHandle plan_2d, plan_1d;

/* Arrays for storing wavefunction, momentum and position op, etc */
extern cufftDoubleComplex *wfc, *wfc0, *wfc_backup, *GK, *GV_half, *GV, *EK, *EV, *EV_opt, *GxPy, *GyPx, *ExPy, *EyPx, *EappliedField;
extern double *Energy, *Energy_gpu, *r, *Phi, *V, *V_opt, *K, *xPy, *yPx, *xPy_gpu, *yPx_gpu;

/* CUDA data buffers for FFT */
extern cufftDoubleComplex *wfc_gpu, *K_gpu, *V_gpu, *par_sum;
extern double *Phi_gpu;

/* CUDA streams */
//extern cudaStream_t streamA, streamB, streamC, streamD;

/* Define global dim3 and threads for grid and thread dim calculation */
//extern dim3 grid;
//extern int threads;

/* Function declarations */
/*
 * arg1 = Function result code from CUDA CUFFT calls.
 * arg2 = String data for name of function called. Prints value to stdout.
 */

 /**
 * @brief	Checks if CUDA operation has succeeded. Prints to stdout.
 * @ingroup	data
 * @param	result Function result code of CUDA operation
 * @param	c Descriptor of CUDA operation
 * @return	0 for success. See CUDA failure codes in cuda.h for other values.
 */
int isError(int result, char* c); //Checks to see if an error has occurred.

/**
* @brief	Performs parallel summation and renormalises the wavefunction
* @ingroup	data
* @param	gpuWfc GPU memory location for wavefunction
* @param	gpuParSum GPU memory location for parallel summation memory space
* @param	xDim Length of X dimension
* @param	yDim Length of Y dimension
* @param	threads Number of CUDA threads for operation
* @return	0 for success. See CUDA failure codes in cuda.h for other values.
*/
void parSum(double2* gpuWfc, double2* gpuParSum, int threads, Grid &par,
            Cuda &cupar);

/**
* @brief	Creates the optical lattice to match the vortex lattice constant
* @ingroup	data
* @param	centre Central vortex in condensate
* @param	V Trapping potential for condensate
* @param	vArray Vortex location array
* @param	num_vortices Number of tracked vortices
* @param	theta_opt Offset angle for optical lattice relative to vortex lattice
* @param	intensity Optical lattice amplitude
* @param	v_opt Optical lattice memory address location
* @param	x X grid array
* @param	y Y grid array
*/
void optLatSetup(struct Vtx::Vortex centre, double* V, 
                 struct Vtx::Vortex *vArray, int num_vortices, double theta_opt,
                 double intensity, double* v_opt, double *x, double *y, 
                 Grid &par);

/**
* @brief	Calculates the energy of the condensate. Not implemented.
* @ingroup	data
* @param	Energy Host array of energy calculated values
* @param	Energy_gpu Device array of energy calculated values
* @param	V_op Potential (position space) operator
* @param	K_op Kinetic (momentum space) operator
* @param	dx Increment along x
* @param	dy Increment along y
* @param	gpuWfc Device wavefunction array
* @param	gState Indicate if imaginary or real time evolution
* @return	$\langle \Psi | H | \Psi \rangle$
*/
double energy_angmom(double* Energy, double* Energy_gpu, double2 *V_op, 
                     double2 *K_op, double2 *gpuWfc, 
                     int gState, Grid &par);

#endif
