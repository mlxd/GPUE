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
#include "../include/evolution.h"
#include "../include/init.h"
#include <string.h>
#include <assert.h>
#include <cufft.h>
#include <vector>

// Test for the Grid structure with paramters in it 
void parameter_test();

// Test for the parsing function
void parser_test();

// Testing the evolve function in evolution.cu
void evolve_test();

// Kernel testing will be added later

/*----------------------------------------------------------------------------//
* MAIN
*-----------------------------------------------------------------------------*/

void test_all(){
    std::cout << "Starting unit tests..." << '\n';
    parameter_test();
    parser_test();
    evolve_test();

    std::cout << "All tests completed. GPUE passed." << '\n';
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
    bool bvar = true;

    // Now testing the Grid class
    Grid grid_test;
    grid_test.store("dstar_var",dstar_var);
    grid_test.store("dvar", dvar);
    grid_test.store("ivar", ivar);
    grid_test.store("bvar", bvar);

    assert(dstar_var == grid_test.dsval("dstar_var"));
    assert(dvar == grid_test.dval("dvar"));
    assert(ivar == grid_test.ival("ivar"));
    assert(bvar == grid_test.bval("bvar"));

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

// Test for the parsing function
void parser_test(){

    // Testing the command-line parser with defaults and with modifications
    std::cout << "Testing command-line parser with no arguments..." << '\n';

    // First testing default values in and out of the parser function
    char *fake_noargv[] = {NULL};
    Grid noarg_grid;
    noarg_grid = parseArgs(0,fake_noargv);

    // Checking contents of noarg_grid:
    assert(noarg_grid.ival("xDim") == 256);
    assert(noarg_grid.ival("yDim") == 256);
    assert(noarg_grid.ival("zDim") == 256);
    assert(noarg_grid.dval("omega") == 0);
    assert(noarg_grid.dval("gammaY") == 1.0);
    assert(noarg_grid.dval("gsteps") == 1e4);
    assert(noarg_grid.dval("esteps") == 1000);
    assert(noarg_grid.dval("gdt") == 1e-4);
    assert(noarg_grid.dval("dt") == 1e-4);
    assert(noarg_grid.ival("device") == 0);
    assert(noarg_grid.ival("atoms") == 1);
    assert(noarg_grid.bval("read_wfc") == false);
    assert(noarg_grid.ival("printSteps") == 100);
    assert(noarg_grid.dval("winding") == 0);
    assert(noarg_grid.bval("corotating") == false);
    assert(noarg_grid.bval("gpe") == false);
    assert(noarg_grid.dval("omegaZ") == 0);
    assert(noarg_grid.dval("int_scaling") == 0);
    assert(noarg_grid.dval("laser_power") == 0);
    assert(noarg_grid.dval("angle_sweep") == 0);
    assert(noarg_grid.ival("kick_it") == 0);
    assert(noarg_grid.bval("write_it") == false);
    assert(noarg_grid.dval("x0_shift") == 0);
    assert(noarg_grid.dval("y0_shift") == 0);
    assert(noarg_grid.dval("sepMinEpsilon") == 0);
    assert(noarg_grid.bval("graph") == false);
    assert(noarg_grid.bval("unit_test") == false);

    // Now testing all values specified by command-line arguments
    std::cout << "Testing command-line parser with all arguments..." << '\n';
    std::vector<std::string> argarray(10);

    char *fake_fullargv[] = {strdup("./gpue"), strdup("-d"), strdup("0"), strdup("-e"), strdup("1000"), strdup("-G"), strdup("1"), strdup("-g"), strdup("1e4"), strdup("-i"), strdup("0"), strdup("-k"), strdup("0"), strdup("-L"), strdup("0"), strdup("-n"), strdup("1"), strdup("-O"), strdup("0"), strdup("-o"), strdup("0"), strdup("-P"), strdup("0"), strdup("-p"), strdup("100"), strdup("-S"), strdup("0"), strdup("-T"), strdup("1e-4"), strdup("-t"), strdup("1e-4"), strdup("-U"), strdup("0"), strdup("-V"), strdup("0"), strdup("-W"), strdup("-w"), strdup("0"), strdup("-X"), strdup("1.0"), strdup("-x"), strdup("256"), strdup("-Y"), strdup("1.0"), strdup("-y"), strdup("256"), strdup("-r"), strdup("-l"), strdup("-a"), strdup("-s"), NULL};
    int fake_argc = sizeof(fake_fullargv) / sizeof(char *) - 1;

    // Now to read into gpue and see what happens
    Grid fullarg_grid;
    fullarg_grid = parseArgs(fake_argc, fake_fullargv);

    // Checking contents of fullarg_grid:
    assert(fullarg_grid.ival("xDim") == 256);
    assert(fullarg_grid.ival("yDim") == 256);
    assert(fullarg_grid.ival("zDim") == 256);
    assert(fullarg_grid.dval("omega") == 0);
    assert(fullarg_grid.dval("gammaY") == 1.0);
    assert(fullarg_grid.dval("gsteps") == 1e4);
    assert(fullarg_grid.dval("esteps") == 1000);
    assert(fullarg_grid.dval("gdt") == 1e-4);
    assert(fullarg_grid.dval("dt") == 1e-4);
    assert(fullarg_grid.ival("device") == 0);
    assert(fullarg_grid.ival("atoms") == 1);
    assert(fullarg_grid.bval("read_wfc") == true);
    assert(fullarg_grid.ival("printSteps") == 100);
    assert(fullarg_grid.dval("winding") == 0);
    assert(fullarg_grid.bval("corotating") == true);
    assert(fullarg_grid.bval("gpe") == true);
    assert(fullarg_grid.dval("omegaZ") == 0);
    assert(fullarg_grid.dval("int_scaling") == 0);
    assert(fullarg_grid.dval("laser_power") == 0);
    assert(fullarg_grid.dval("angle_sweep") == 0);
    assert(fullarg_grid.ival("kick_it") == 0);
    assert(fullarg_grid.bval("write_it") == true);
    assert(fullarg_grid.dval("x0_shift") == 0);
    assert(fullarg_grid.dval("y0_shift") == 0);
    assert(fullarg_grid.dval("sepMinEpsilon") == 0);
    assert(fullarg_grid.bval("graph") == true);
    assert(fullarg_grid.dval("omegaY") == 1.0);
    assert(fullarg_grid.dval("omegaX") == 1.0);
    assert(fullarg_grid.bval("unit_test") == false);

}

// Testing the evolve function in evolution.cu
void evolve_test(){
    // First, we need to create all the necessary data structures for the
    // The evolve function, FOLLOWING INIT.CU

    std::cout << "Testing the evolve function" << '\n';

    char * fake_argv[] = {strdup("./gpue"), strdup("-d"), strdup("0"), strdup("-e"), strdup("2.01e4"), strdup("-G"), strdup("1.0"), strdup("-g"), strdup("0"), strdup("-i"), strdup("1.0"), strdup("-k"), strdup("0"), strdup("-L"), strdup("0"), strdup("-n"), strdup("1e6"), strdup("-O"), strdup("0.0"), strdup("-o"), strdup("0.0"), strdup("-P"), strdup("0.0"), strdup("-p"), strdup("1000"), strdup("-S"), strdup("0.0"), strdup("-T"), strdup("1e-4"), strdup("-t"), strdup("1e-4"), strdup("-U"), strdup("0"), strdup("-V"), strdup("0"), strdup("-w"), strdup("0.0"), strdup("-X"), strdup("1.0"), strdup("-x"), strdup("256"), strdup("-Y"), strdup("1.0"), strdup("-y"), strdup("256"), strdup("-W"), strdup("-D"), strdup("data"), NULL};
    int fake_argc = sizeof(fake_argv) / sizeof(char *) - 1;

    // Now to read into gpue and see what happens
    Grid par;
    par = parseArgs(fake_argc, fake_argv);

    Wave wave;
    Op opr;
    Cuda cupar;

    std::cout << "omegaX is: " << par.dval("omegaX") << '\n';
    std::cout << "x / yDim are: " << par.ival("xDim") << '\t' 
              << par.ival("yDim") << '\n';
    int device = par.ival("device");
    cudaSetDevice(device);

    std::string buffer;

    //************************************************************//
    /*
    * Initialise the Params data structure to track params and variables
    */
    //************************************************************//

    initialise(opr, cupar, par, wave);

    // Re-establishing variables from parsed Grid class
    double dx = par.dval("dx");
    double dy = par.dval("dy");
    double *x = par.dsval("x");
    double *y = par.dsval("y");
    double *V_opt = opr.dsval("V_opt");
    double *xPy = opr.dsval("xPy");
    double *yPx = opr.dsval("yPx");
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    bool read_wfc = par.bval("read_wfc");
    int gsteps = par.ival("gsteps");
    int esteps = par.ival("esteps");
    cufftDoubleComplex *wfc = wave.cufftDoubleComplexval("wfc");
    cufftDoubleComplex *V_gpu = opr.cufftDoubleComplexval("V_gpu");
    cufftDoubleComplex *GK = opr.cufftDoubleComplexval("GK");
    cufftDoubleComplex *GV = opr.cufftDoubleComplexval("GV");
    cufftDoubleComplex *EV = opr.cufftDoubleComplexval("EV");
    cufftDoubleComplex *EK = opr.cufftDoubleComplexval("EK");
    cufftDoubleComplex *ExPy = opr.cufftDoubleComplexval("ExPy");
    cufftDoubleComplex *EyPx = opr.cufftDoubleComplexval("EyPx");
    cufftDoubleComplex *wfc_gpu = wave.cufftDoubleComplexval("wfc_gpu");
    cufftDoubleComplex *K_gpu = opr.cufftDoubleComplexval("K_gpu");
    cufftDoubleComplex *xPy_gpu = opr.cufftDoubleComplexval("xPy_gpu");
    cufftDoubleComplex *yPx_gpu = opr.cufftDoubleComplexval("yPx_gpu");
    cufftDoubleComplex *par_sum = wave.cufftDoubleComplexval("par_sum");
    cudaError_t err = cupar.cudaError_tval("err");

    std::cout << "variables re-established" << '\n';
    std::cout << read_wfc << '\n';

    std::cout << "omegaY is: " << par.ival("omegaY") << '\t'
              << "omegaX is: " << par.dval("omegaX") << '\n';
   
    if(gsteps > 0){
        err=cudaMemcpy(K_gpu, GK, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess)
            exit(1);
        err=cudaMemcpy(V_gpu, GV, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess)
            exit(1);
        err=cudaMemcpy(xPy_gpu, xPy, sizeof(double)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess)
            exit(1);
        err=cudaMemcpy(yPx_gpu, yPx, sizeof(double)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess)
            exit(1);
        err=cudaMemcpy(wfc_gpu, wfc, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess)
            exit(1);
    
        evolve(wave, opr, par_sum,
               gsteps, cupar, 0, 0, par, buffer);
        wfc = wave.cufftDoubleComplexval("wfc");
        wfc_gpu = wave.cufftDoubleComplexval("wfc_gpu");
        cudaMemcpy(wfc, wfc_gpu, sizeof(cufftDoubleComplex)*xDim*yDim,
                   cudaMemcpyDeviceToHost);
    }

    std::cout << GV[0].x << '\t' << GK[0].x << '\t'
              << xPy[0] << '\t' << yPx[0] << '\n';

    //free(GV); free(GK); free(xPy); free(yPx);

    // Re-initializing wfc after evolution
    wfc = wave.cufftDoubleComplexval("wfc");
    wfc_gpu = wave.cufftDoubleComplexval("wfc_gpu");

    std::cout << "evolution started..." << '\n';
    std::cout << "esteps: " << esteps << '\n';

    //************************************************************//
    /*
    * Evolution
    */
    //************************************************************//
    if(esteps > 0){
        err=cudaMemcpy(xPy_gpu, ExPy, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << 1 << '\n';
            exit(1);
        }
        err=cudaMemcpy(yPx_gpu, EyPx, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << 2 << '\n';
            exit(1);
        }
        err=cudaMemcpy(xPy_gpu, ExPy, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << 3 << '\n';
            exit(1);
        }
        err=cudaMemcpy(yPx_gpu, EyPx, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << 4 << '\n';
            exit(1);
        }
        err=cudaMemcpy(K_gpu, EK, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << 5 << '\n';
            exit(1);
        }
        err=cudaMemcpy(V_gpu, EV, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << 6 << '\n';
            exit(1);
        }
        err=cudaMemcpy(wfc_gpu, wfc, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << 7 << '\n';
            exit(1);
        }

        evolve(wave, opr, par_sum,
               esteps, cupar, 1, 0, par, buffer);

        wfc = wave.cufftDoubleComplexval("wfc");
        wfc_gpu = wave.cufftDoubleComplexval("wfc_gpu");
    }

    std::cout << "done evolving" << '\n';
    free(EV); free(EK); free(ExPy); free(EyPx);
    free(x);free(y);
    cudaFree(wfc_gpu); cudaFree(K_gpu); cudaFree(V_gpu); cudaFree(yPx_gpu);
    cudaFree(xPy_gpu); cudaFree(par_sum);

    // At this point, we have a wavefunction that is testable, which we will be
    // doing in much the same way as in the linear/perf branch of GPUE.
    // For this, we must recreate the en.py file in a testable format in cpp
    // Note that we could be using the GPUs for this, but because it is a unit
    // test and we do not care that much about perfomance, we will be using the 
    // CPU instead. We may later add in the appropriate GPU kernels.

    // We first need to grab the wavefunctions from the evolve function
    // After evolution

    wfc = wave.cufftDoubleComplexval("wfc");
    wfc_gpu = wave.cufftDoubleComplexval("wfc_gpu");

}

// Performs simple trapezoidal integral -- following python notation
// Note: Because of the shape of the wfc array, we may need to create a temp 
//       array to do the actual integration here. Look into it!
double trapz(double *array, int dimension, double dx){
    double integral = 0;
    for (int i = 1; i < dimension; ++i){
        integral += (array[i-1] + array[i]) * 0.5 * dx;
    }
    return integral;
}
