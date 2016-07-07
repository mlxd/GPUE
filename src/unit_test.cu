/*
* unit_test.cc - GPUE: Split Operator based GPU solver for Nonlinear 
Schrodinger Equation, Copyright (C) 2011-2015, Lee J. O'Riordan 
<loriordan@gmail.com>, Tadhg Morgan, Neil Crowley. All rights reserved.

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

#include "../include/ds.h"
#include "../include/unit_test.h"
#include "../include/parser.h"
#include <string>
#include <assert.h>
#include <cufft.h>

// Test for the Grid structure with paramters in it 
void parameter_test();

// Test for the parsing function
void parser_test();

/*----------------------------------------------------------------------------//
* MAIN
*-----------------------------------------------------------------------------*/

void test_all(){
    std::cout << "Starting unit tests..." << '\n';
    parameter_test();
}

// Test for the Grid structure with paramters in it
// Initialize all necessary variables and read them back out
void parameter_test(){
    // For this test, we simply need to read in and out stuff from each 
    // class and structure in ds.h / ds.cc
    
    // Certain variables will be used multiple times. 
    double *dstar_var;
    dstar_var = (double *)malloc(sizeof(double) * 5);
    cufftDoubleComplex *cdc_var;
    cdc_var = (cufftDoubleComplex *)malloc(sizeof(cufftDoubleComplex) * 5);
    for (int i = 0; i < 5; ++i){
        dstar_var[i] = (double)i * 0.5;
        cdc_var[i].x = (double)i * 0.5;
        cdc_var[i].y = (double)i * 0.5;
    }

    double dvar = 1.05;
    int ivar = 5;

    // Now testing the Grid class
    Grid grid_test;
    grid_test.store("dstar_var",dstar_var);
    grid_test.store("dvar", dvar);
    grid_test.store("ivar", ivar);

    assert(dstar_var == grid_test.dsval("dstar_var"));
    assert(dvar == grid_test.dval("dvar"));
    assert(ivar == grid_test.ival("ivar"));

    std::cout << "Grid class checked, now checking the Cuda class..." << '\n';

    // Now checking the Cuda class
    // This one will require creating a list of variables...
    cudaError_t err = cudaSuccess;
    cufftHandle plan_1d = 4, plan_2d = 6;
    cudaStream_t streamA = 0, streamB = 0, streamC = 0, streamD = 0;
    cufftResult result = CUFFT_SUCCESS;
    dim3 grid;

    grid.x = 1; grid.y = 2; grid.z = 3;

    // Creating Cuda class to test with
    Cuda cuda_test;

    // Testing the store and value functions
    cuda_test.store("err", err);
    cuda_test.store("result", result);
    cuda_test.store("plan_1d", plan_1d);
    cuda_test.store("plan_2d", plan_2d);
    cuda_test.store("streamA", streamA);
    cuda_test.store("streamB", streamB);
    cuda_test.store("streamC", streamC);
    cuda_test.store("streamD", streamD);
    cuda_test.store("grid", grid);

    assert(err == cuda_test.cudaError_tval("err"));
    assert(result == cuda_test.cufftResultval("result"));
    assert(plan_1d == cuda_test.cufftHandleval("plan_1d"));
    assert(plan_2d == cuda_test.cufftHandleval("plan_2d"));
    assert(streamA == cuda_test.cudaStream_tval("streamA"));
    assert(streamB == cuda_test.cudaStream_tval("streamB"));
    assert(streamC == cuda_test.cudaStream_tval("streamC"));
    assert(streamD == cuda_test.cudaStream_tval("streamD"));
    assert(grid.x == cuda_test.dim3val("grid").x);
    assert(grid.y== cuda_test.dim3val("grid").y);
    assert(grid.y == cuda_test.dim3val("grid").y);

    std::cout << "Cuda class checked, now checking Op class..." << '\n';

    // Now checking the Op class
    // Creating Op class to test with
    Op op_test;
    op_test.store("dstar_var",dstar_var);
    op_test.store("cdc_var",cdc_var);

    assert(dstar_var == op_test.dsval("dstar_var"));
    assert(cdc_var == op_test.cufftDoubleComplexval("cdc_var"));

    std::cout << "Op class checked, now checking Wave class..." << '\n';

    // Now checking the Op class
    // Creating Op class to test with
    Wave wave_test; 
    wave_test.store("dstar_var",dstar_var);
    wave_test.store("cdc_var",cdc_var);

    assert(dstar_var == wave_test.dsval("dstar_var"));
    assert(cdc_var == wave_test.cufftDoubleComplexval("cdc_var"));

    std::cout << "All data structures checked" << '\n';

}
