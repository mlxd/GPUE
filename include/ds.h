///@cond LICENSE
/*** ds.h - GPUE: Split Operator based GPU solver for Nonlinear
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
 *  @file    ds.h
 *  @author  Lee J. O'Riordan (mlxd)
 *  @date    12/11/2015
 *  @version 0.1
 *
 *  @brief Dastructure for simulation runtime parameters
 *
 *  @section DESCRIPTION
 *  This file keeps track of and generates a parameter file (INI format) of the
 *	simulation parameters used. The resulting file is read in by the Python
 *	post-proc/analysis functions.
 */
 //##############################################################################

#ifndef DS_H
#define DS_H
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <unordered_map>

/**
 * @brief       Class to hold the variable map and grid information
 * @ingroup     data
 */
// NOTE: This is necessary if we ever want to do dynamic grid manipulation.
class Grid{
    // Here we keep our variable map (unordered for performance)
    // and also grid information
    private:
        std::unordered_map<std::string, double> param_double;
        std::unordered_map<std::string, int> param_int;

    // Here we keep the functions to store variables and access grid data
    public:
        void store(std::string, int iparam);
        void store(std::string, double dparam);
        double val_int(std::string id);
        int val_double(std::string id);
};
typedef class Grid Grid;

// NOTE: We may want to combine the Op and Wave class. It depends on the types
//       of functions we need to have in Wave

/**
 * @brief       class to hold all necessary information about the operators
 * @ingroup     data
 */
class Op{
};
typedef class Op Op;

/**
 * @brief       class to hold all necessary information about the wavefunction
 * @ingroup     data
 */
class Wave{
};
typedef class Wave Wave;


// DEPRECATION WARNING

// Gathers all data from command-line parsing
struct Param{ 
    char title[32];
    double data;
};
typedef struct Param Param;

struct Array{
    Param *array;
    size_t length;
    size_t used;
};
typedef struct Array Array;

/**
* @brief	Intialises Array to specified length
* @ingroup	data
* @param	*arr Pointer to parameter storage Array
* @param	initLength Length to initialise Array
*/
void initArr(Array *arr, size_t initLen);
/**
* @brief	Adds data to the parameter storage Array
* @ingroup	data
* @param	*arr Pointer to parameter storage Array
* @param	*t Parameter name to be saved
* @param	d Double value of parameter
*/
void appendData(Array *arr, std::string t, double d);
/**
* @brief	Free all allocated memory from Array
* @ingroup	data
* @param	*arr Pointer to parameter storage Array
*/
void freeArray(Array *arr);
/**
* @brief	Allocate new Param. Unused.
* @ingroup	data
* @param	*arr Pointer to parameter storage Array
* @param	*t Parameter name to be saved
* @param	d Double value of parameter
*/
Param newParam(std::string t,double d);
#endif

// END DEPRECATION WARNING
