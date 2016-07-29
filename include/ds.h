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
 //#############################################################################

#ifndef DS_H
#define DS_H
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <typeinfo>
#include <iostream>

/*----------------------------------------------------------------------------//
* CLASSES
*-----------------------------------------------------------------------------*/

/**
 * @brief       Class to hold the variable map and grid information
 * @ingroup     data
 */
// NOTE: This is necessary if we ever want to do dynamic grid manipulation.
// NOTE: I do not like this integer double split for retrieval. Find solution.
// NOTE: Consider changing unordered maps to switches for performance
// NOTE: Add FileIO to public functions
class Grid{
    // Here we keep our variable map (unordered for performance)
    // and also grid information. Note that dx dy, and dz are in param_double
    private:
        std::unordered_map<std::string, int> param_int;
        std::unordered_map<std::string, double> param_double;
        std::unordered_map<std::string, double*> param_dstar;
        std::unordered_map<std::string, bool> param_bool;
        std::string data_dir;

        // List of all strings for parsing into the appropriate param map
        // 1 -> int, 2 -> double, 3 -> double*
        std::unordered_map<std::string, int> id_list;

    // Here we keep the functions to store variables and access grid data
    public:

        // placing grid parameters in public for now
        double *x, *y, *z, *xp, *yp, *zp;

        // Function to store integer into param_int
        void store(std::string id, int iparam);

        // Function to store double into param_double
        void store(std::string id, double dparam);

        // Function to store double* into param_dstar
        void store(std::string id, double *dsparam);

        // Function to store bool into param_bool
        void store(std::string id, bool bparam);

        // Function to store string into data_dir
        void store(std::string id, std::string sparam);

        // Function to retrieve integer value from param_int
        int ival(std::string id);

        // Function to retrieve double value from param_double
        double dval(std::string id);

        // Function to retrieve double star values from param_dstar
        double *dsval(std::string id);

        // Function to retrieve bool from param_bool
        bool bval(std::string id);

        // Fucntion to retrieve string from data_dir
        std::string sval(std::string id);

        // Function for file writing
        void write(std::string filename);
};
typedef class Grid Grid;

/**
 * @brief        Class to hold CUDA-specific variables and features
 * @ingroup      data
 */
// I will not be using the unordered map for this one because the number of 
// variables stored is small
class Cuda{
    private:
        cudaError_t err;
        cufftResult result;
        cufftHandle plan_1d, plan_2d;
        cudaStream_t streamA, streamB, streamC, streamD;
        dim3 grid;
        int threads;
    public:

        // Functions to store data
        void store(std::string id, cudaError_t errin);
        void store(std::string id, cufftResult resultin);
        void store(std::string id, cufftHandle planin);
        void store(std::string id, cudaStream_t streamin);
        void store(std::string id, dim3 gridin);
        //void store(std::string id, int threadsin);

        // Functions to retrieve data
        cudaError_t cudaError_tval(std::string id);
        cufftResult cufftResultval(std::string id);
        cufftHandle cufftHandleval(std::string id);
        cudaStream_t cudaStream_tval(std::string id);
        dim3 dim3val(std::string id);
        //int ival(std::string id);

};
//typedef class Cuda Cuda;

// NOTE: We may need to store the operators as a 3 1d vectors 
//       instead of 1 3d one.

/**
 * @brief       class to hold all necessary information about the operators
 * @ingroup     data
 */
class Op{
    private:
        std::unordered_map<std::string, double*> Op_dstar;
        std::unordered_map<std::string, cufftDoubleComplex*> Op_cdc;
        // double *V, *V_opt, *K, *xPy, *yPx, *xPy_gpu, *yPx_gpu;
        //cufftDoubleComplex *GK,*GV_half,*GV,*EK,*EV,*EV_opt,*GxPy,*GyPx,
        //                   *ExPy,*EyPx,*EappliedField,*K_gpu,*V_gpu;
    public:

        // Functions to store data
        void store(std::string id, double *data);
        void store(std::string id, cufftDoubleComplex *data);

        // Functions to retrieve data
        double *dsval(std::string id);
        cufftDoubleComplex *cufftDoubleComplexval(std::string id);
};
typedef class Op Op;

/**
 * @brief       class to hold all necessary information about the wavefunction
 * @ingroup     data
 */
class Wave{
    private:
        std::unordered_map<std::string, double*> Wave_dstar;
        std::unordered_map<std::string, cufftDoubleComplex*> Wave_cdc;
        //double *Energy, *energy_gpu, *r, *Phi, *Phi_gpu;
        //cufftDoubleComplex *wfc, *wfc0, *wfc_backup, *wfc_gpu, *par_sum;
    public:

        // functions to store data
        void store(std::string id, double *data);
        void store(std::string id, cufftDoubleComplex *data);


        // Functions to retrieve data
        double *dsval(std::string id);
        cufftDoubleComplex *cufftDoubleComplexval(std::string id);
};
typedef class Wave Wave;

#endif
