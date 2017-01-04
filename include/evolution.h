///@endcond
//##############################################################################
/**
 *  @file    evolution.h
 *  @author  James Ryan Schloss (leios)
 *  @date    5/31/2016
 *  @version 0.1
 *
 *  @brief function for evolution.
 *
 *  @section DESCRIPTION
 *  These functions and variables are necessary for carrying out the GPUE
 *	simulations. This file will be re-written in an improved form in some
 *	future release.
 */
//##############################################################################

#ifndef EVOLUTION_H
#define EVOLUTION_H

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
#include "split_op.h"
#include "kernels.h"
#include "constants.h"
#include "fileIO.h"
#include "lattice.h"
#include "manip.h"
#include "unit_test.h"


/* Function declarations */
/*
 * arg1 = Function result code from CUDA CUFFT calls.
 * arg2 = String data for name of function called. Prints value to stdout.
 */

// UPDATE LIST LATER
 /**
 * @brief	performs real or imaginary time evolution
 * @ingroup	data
 * @param	result Function result code of CUDA operation
 * @param	c Descriptor of CUDA operation
 * @return	0 for success. See CUDA failure codes in cuda.h for other values.
 */
void evolve_2d(Wave &wave, Op &opr,
            cufftDoubleComplex *gpuParSum, int numSteps, Cuda &cupar,
            unsigned int gstate, Grid &par, 
            std::string buffer);

// UPDATE LIST LATER
 /**
 * @brief       performs real or imaginary time evolution
 * @ingroup     data
 * @param       result Function result code of CUDA operation
 * @param       c Descriptor of CUDA operation
 * @return      0 for success. See CUDA failure codes in cuda.h for other values
 */
void evolve_3d(Wave &wave, Op &opr,
            cufftDoubleComplex *gpuParSum, int numSteps, Cuda &cupar,
            unsigned int gstate, Grid &par, 
            std::string buffer);

#endif
