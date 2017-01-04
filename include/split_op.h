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
#include "constants.h"

#ifdef __linux
	#include<omp.h>
#elif __APPLE__
	//printf("OpenMP support disabled due to Clang/LLVM being behind the trend.",);
#endif

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
void parSum(double2* gpuWfc, double2* gpuParSum, Grid &par,
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
                 Grid &par, Op &opr);

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
