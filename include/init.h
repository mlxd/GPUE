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

#ifndef INIT_H
#define INIT_H

#include "../include/split_op.h"
#include "../include/kernels.h"
#include "../include/constants.h"
#include "../include/fileIO.h"
#include "../include/tracker.h"
#include "../include/minions.h"
#include "../include/parser.h"
#include "../include/ds.h"
#include "../include/unit_test.h"
#include "../include/operators.h"

#include "../include/lattice.h"
#include "../include/node.h"
#include "../include/edge.h"
#include "../include/manip.h"
#include "../include/vort.h"
#include "../include/evolution.h"
#include <string>
#include <iostream>
/* Function declarations */
/*
 * arg1 = Function result code from CUDA CUFFT calls.
 * arg2 = String data for name of function called. Prints value to stdout.
 */

// UPDATE LIST LATER
 /**
 * @brief	Initializes data structures
 * @ingroup	data
 * @param	Operator class
 * @param	Cuda class
 * @param	Grid class
 * @param	Wave class
 */
int init_2d(Op &opr, Cuda &cupar, Grid &par, Wave &wave);

 /**
 * @brief       Initializes data structures
 * @ingroup     data
 * @param       Operator class
 * @param       Cuda class
 * @param       Grid class
 * @param       Wave class
 */
int init_3d(Op &opr, Cuda &cupar, Grid &par, Wave &wave);

#endif
