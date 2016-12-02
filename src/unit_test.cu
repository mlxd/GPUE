/*
* unit_test.cu - GPUE: Split Operator based GPU solver for Nonlinear 
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

// Testing the evolve_2d function in evolution.cu
void evolve_2d_test();

// Testing the parSum function
void parSum_test();

// Simple test of grid / cuda stuff
void grid_test();

// Kernel testing will be added later

/*----------------------------------------------------------------------------//
* MAIN
*-----------------------------------------------------------------------------*/

void test_all(){
    std::cout << "Starting unit tests..." << '\n';
    //parameter_test();
    //parser_test();
    //evolve_2d_test();
    //grid_test();
    parSum_test();

    std::cout << "All tests completed. GPUE passed." << '\n';
}

// Simple test of CUDA grid stuff
void grid_test(){

    std::cout << "testing grid / threads and stuff" << '\n';

    int max_threads = 128;

    int xDim = 1024;
    int yDim = 1024;
    int zDim = 1;

    int xD = 1, yD = 1, zD = 1;

    int gsize = xDim * yDim;

    // Now to set up the CUDA grid / threads
    dim3 block;
    dim3 grid;

    if (xDim <= max_threads){
        block.x = xDim;
        block.y = 1;
        block.z = 1;

        xD = 1;
        yD = yDim;
        zD = 1;
    } 
    else{
        int count = 0;
        int dim_tmp = xDim;
        while (dim_tmp > max_threads){
            count++;
            dim_tmp /= 2;
        }

        std::cout << "count is: " << count << '\n';

        block.x = dim_tmp;
        block.y = 1;
        block.z = 1;
        xD = pow(2,count);
        yD = yDim;
        zD = 1;
    }

    std::cout << "threads in x are: " << block.x << '\n';
    std::cout << "dimensions are: " << xD << '\t' << yD << '\t' << zD << '\n';

    grid.x=xD; 
    grid.y=yD; 
    grid.z=zD; 

    int total_threads = block.x * block.y * block.z;

    // Now we need to initialize our double * and send it to the gpu
    double *host_array, *device_array;
    host_array = (double *) malloc(sizeof(double)*gsize);
    cudaMalloc((void**) &device_array, sizeof(double)*gsize);

    // initializing 2d array
    for (int i = 0; i < gsize; i++){
        host_array[i] = -1;
    }

    // Now to copy to device
    cudaMemcpy(device_array, host_array,
               sizeof(double)*gsize,
               cudaMemcpyHostToDevice);

    // Test
    thread_test<<<grid,block>>>(device_array,device_array);

    // Now to copy back and print
    cudaMemcpy(host_array, device_array,
               sizeof(double)*gsize,
               cudaMemcpyDeviceToHost);
    
    
/*
    for (int i = 0; i < gsize; i++){
        std::cout << i << '\t' <<  host_array[i] << '\n';
    }
*/
    std::cout << "1024 x 1024 is: " << host_array[gsize-1] << '\n';
    assert(host_array[gsize-1] == 1024*1024-1);

    std::cout << "2d grid tests completed. now for 3d cases" << '\n';

    // Insert 3d tests

    std::cout << "Grid test concluded" << '\n';
}

// Test of the parSum function in 3d
void parSum_test(){

    // Setting error
    cudaError_t err;

    // first, we need to initialize the Grid and Cuda classes
    Grid par;
    Cuda cupar;

    // 2D test first

    // For now, we will assume an 8x8 array for summing
    dim3 threads(32, 1, 1);
    int total_threads = threads.x*threads.y*threads.z;
    int xDim = 64;
    int yDim = 64;
    int zDim = 1;

    par.store("dimnum", 2);
    par.store("xDim", xDim);
    par.store("yDim", yDim);
    par.store("zDim", zDim);
    par.store("dx",1.0);
    par.store("dy",1.0);
    par.store("dz",1.0);
    cupar.store("threads",threads);

    // Now we need to initialize the grid for the getGid3d3d kernel
    int gsize = xDim*yDim;
    dim3 grid;
    grid.x = 2;
    grid.y = yDim;

    cupar.store("grid", grid);

    // now we need to initialize the wfc to all 1's;
    double2 *wfc, *host_sum;
    wfc = (cufftDoubleComplex *) malloc(sizeof(cufftDoubleComplex) * gsize);
    host_sum = (cufftDoubleComplex *) 
               malloc(sizeof(cufftDoubleComplex) * gsize / total_threads);

    // init wfc
    for (int i = 0; i < gsize; i++){
        wfc[i].x = 1;
        wfc[i].y = 0;
    }

    double2 *gpu_wfc;
    cudaMalloc((void**) &gpu_wfc, sizeof(cufftDoubleComplex)*gsize);

    // copying wfc to device
    err = cudaMemcpy(gpu_wfc, wfc, sizeof(cufftDoubleComplex)*gsize,
                     cudaMemcpyHostToDevice);

    if (err!=cudaSuccess){
        std::cout << "ERROR: Could not copy wfc to device!" << '\n';
    }

    // Creating parsum on device
    double2 *par_sum;
    cudaMalloc((void**) &par_sum, sizeof(cufftDoubleComplex)*gsize/total_threads);

    parSum(gpu_wfc, par_sum, par, cupar);

    // copying parsum back
    err = cudaMemcpy(host_sum, par_sum, 
                     sizeof(cufftDoubleComplex)*gsize / total_threads, 
                     cudaMemcpyDeviceToHost);
    if (err!=cudaSuccess){
        std::cout << err << '\n';
        std::cout << "ERROR: Could not copy par_sum to the host!" << '\n';
        exit(1);
    }

    // The output value should be 4096
    std::cout << "2d parSum is:" << '\n';
    std::cout << host_sum[0].x << " + " << host_sum[0].y << " i" << '\n';

    if (host_sum[0].x != 4096){
        std::cout << "parSum 2d test has failed! Sum is: "
                  << host_sum[0].x << '\n';
        assert((int)host_sum[0].x == 4096);
    }

    // Now for the 3d case
    // For now, we will assume a 4x4x4 array for summing
    par.store("dimnum", 3);
    par.store("xDim", 16);
    par.store("yDim", 16);
    par.store("zDim", 16);
    par.store("dx",1.0);
    par.store("dy",1.0);
    par.store("dz",1.0);
    cupar.store("threads",threads);

    // Now we need to initialize the grid for the getGid3d3d kernel
    grid.x = gsize;
    grid.y = 1;
    grid.z = 1;

    cupar.store("grid", grid);

    // copying host wfc back to device
    err = cudaMemcpy(gpu_wfc, wfc, sizeof(cufftDoubleComplex)*gsize,
                     cudaMemcpyHostToDevice);

    parSum(gpu_wfc, par_sum, par, cupar);

    // copying parsum back
    err = cudaMemcpy(host_sum, par_sum, 
                     sizeof(cufftDoubleComplex)*gsize / total_threads, 
                     cudaMemcpyDeviceToHost);
    if (err!=cudaSuccess){
        std::cout << "ERROR: Could not copy par_sum to the host!" << '\n';
        exit(1);
    }

    std::cout << "3d parSum is:" << '\n';
    std::cout << host_sum[0].x << " + " << host_sum[0].y << " i" << '\n';

    if (host_sum[0].x != 4096){
        std::cout << "parSum 3d test has failed!" << '\n';
        assert((int)host_sum[0].x == 4096);
    }

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

    // I apologize for the mess... If you have a better way of creating the 
    // char ** for this without running into memory issues, let me know!
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

// Testing the evolve_2d function in evolution.cu
void evolve_2d_test(){
    // First, we need to create all the necessary data structures for the
    // The evolve_2d function, FOLLOWING INIT.CU

    std::cout << "Testing the evolve_2d function" << '\n';

    // Note: the omega_z value (-o flag) is arbitrary
    char * fake_argv[] = {strdup("./gpue"), strdup("-d"), strdup("0"), strdup("-e"), strdup("2.01e4"), strdup("-G"), strdup("1.0"), strdup("-g"), strdup("0"), strdup("-i"), strdup("1.0"), strdup("-k"), strdup("0"), strdup("-L"), strdup("0"), strdup("-n"), strdup("1e6"), strdup("-O"), strdup("0.0"), strdup("-o"), strdup("10.0"), strdup("-P"), strdup("0.0"), strdup("-p"), strdup("1000"), strdup("-S"), strdup("0.0"), strdup("-T"), strdup("1e-4"), strdup("-t"), strdup("1e-4"), strdup("-U"), strdup("0"), strdup("-V"), strdup("0"), strdup("-w"), strdup("0.0"), strdup("-X"), strdup("1.0"), strdup("-x"), strdup("256"), strdup("-Y"), strdup("1.0"), strdup("-y"), strdup("256"), strdup("-W"), strdup("-D"), strdup("data"), NULL};
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

    init_2d(opr, cupar, par, wave);

    // Re-establishing variables from parsed Grid class
    double dx = par.dval("dx");
    double dy = par.dval("dy");
    double *x = par.dsval("x");
    double *y = par.dsval("y");
    double *V_opt = opr.dsval("V_opt");
    double *pAy = opr.dsval("pAy");
    double *pAx = opr.dsval("pAx");
    double *pAy_gpu = opr.dsval("pAy_gpu");
    double *pAx_gpu = opr.dsval("pAx_gpu");
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
    cufftDoubleComplex *EpAy = opr.cufftDoubleComplexval("EpAy");
    cufftDoubleComplex *EpAx = opr.cufftDoubleComplexval("EpAx");
    cufftDoubleComplex *GpAx = opr.cufftDoubleComplexval("GpAx");
    cufftDoubleComplex *GpAy = opr.cufftDoubleComplexval("GpAy");
    cufftDoubleComplex *wfc_gpu = wave.cufftDoubleComplexval("wfc_gpu");
    cufftDoubleComplex *K_gpu = opr.cufftDoubleComplexval("K_gpu");
    cufftDoubleComplex *par_sum = wave.cufftDoubleComplexval("par_sum");
    cudaError_t err = cupar.cudaError_tval("err");

    std::cout << "variables re-established" << '\n';
    std::cout << read_wfc << '\n';

    std::cout << "omegaY is: " << par.ival("omegaY") << '\t'
              << "omegaX is: " << par.dval("omegaX") << '\n';

/*
    for (int i = 0; i < xDim * yDim; ++i){
        std::cout << i << '\t' << wfc[i].x << '\t' << wfc[i].y << '\n';
    }
*/

    std::cout << "gsteps: " << gsteps << '\n';
   
    if(gsteps > 0){
        err=cudaMemcpy(K_gpu, GK, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy K_gpu to device" << '\n';
            exit(1);
        }
        err=cudaMemcpy(V_gpu, GV, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy V_gpu to device" << '\n';
            exit(1);
        }
        err=cudaMemcpy(pAy_gpu, GpAy, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy pAy_gpu to device" << '\n';
            exit(1);
        }
        err=cudaMemcpy(pAx_gpu, GpAx, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy pAx_gpu to device" << '\n';
            exit(1);
        }
        err=cudaMemcpy(wfc_gpu, wfc, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy wfc_gpu to device" << '\n';
            exit(1);
        }
    
        evolve_2d(wave, opr, par_sum,
               gsteps, cupar, 0, par, buffer);
        wfc = wave.cufftDoubleComplexval("wfc");
        wfc_gpu = wave.cufftDoubleComplexval("wfc_gpu");
        cudaMemcpy(wfc, wfc_gpu, sizeof(cufftDoubleComplex)*xDim*yDim,
                   cudaMemcpyDeviceToHost);
    }

    std::cout << GV[0].x << '\t' << GK[0].x << '\t'
              << pAy[0] << '\t' << pAx[0] << '\n';

    //free(GV); free(GK); free(pAy); free(pAx);

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
        err=cudaMemcpy(pAy_gpu, EpAy, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy pAy_gpu to device" << '\n';
            exit(1);
        }
        err=cudaMemcpy(pAx_gpu, EpAx, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy pAx_gpu to device" << '\n';
            exit(1);
        }
        err=cudaMemcpy(K_gpu, EK, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy K_gpu to device" << '\n';
            exit(1);
        }
        err=cudaMemcpy(V_gpu, EV, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy V_gpu to device" << '\n';
            exit(1);
        }
        err=cudaMemcpy(wfc_gpu, wfc, sizeof(cufftDoubleComplex)*xDim*yDim,
                       cudaMemcpyHostToDevice);
        if(err!=cudaSuccess){
            std::cout << "ERROR: Could not copy wfc_gpu to device" << '\n';
            exit(1);
        }

        evolve_2d(wave, opr, par_sum,
               esteps, cupar, 1, par, buffer);

    }

    std::cout << "done evolving, checking result" << '\n';

    // At this point, we have a wavefunction that is testable, which we will be
    // doing in much the same way as in the linear/perf branch of GPUE.
    // For this, we must recreate the en.py file in a testable format in cpp
    // Note that we could be using the GPUs for this, but because it is a unit
    // test and we do not care that much about perfomance, we will be using the 
    // CPU instead. We may later add in the appropriate GPU kernels.

    // We first need to grab the wavefunctions from the evolve_2d function
    // After evolution
    wfc = wave.cufftDoubleComplexval("wfc");
    wfc_gpu = wave.cufftDoubleComplexval("wfc_gpu");
    unsigned int gSize = xDim * yDim;

    // Now to grab K and V, note that these are different than the values used 
    // for K / V_gpu or for E / G K / V in the evolve_2d function
    // The additional 0 in the gpu variable name indicate this (sorry)
    double *K_0_gpu = opr.dsval("K");
    double *K = opr.dsval("K");
    double *V_0_gpu = opr.dsval("V");
    double *V = opr.dsval("V");

    // Now we need som CUDA specific variables for the kernels later on...
    int threads = par.ival("threads");
    dim3 grid = cupar.dim3val("grid");

    // Momentum-space (p) wavefunction
    double2 *wfc_p = wfc;
    double2 *wfc_p_gpu = wfc_gpu;

    // Conjugate (c) wavefunction
    double2 *wfc_c = wfc;
    double2 *wfc_c_gpu = wfc_gpu;

    // Energies
    double2 *Energy_1, *Energy_2, *Energy_k, *Energy_v;
    Energy_1 = wfc_gpu;
    Energy_2 = wfc_gpu;

    // Plan for 2d FFT
    cufftHandle plan_2d = cupar.cufftHandleval("plan_2d");

    std::cout << "allocating space on device..." << '\n';

    // Allocating space on GPU
    cudaMalloc((void **) &wfc_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void **) &K_0_gpu, sizeof(double) * gSize);
    cudaMalloc((void **) &V_0_gpu, sizeof(double) * gSize);
    cudaMalloc((void **) &wfc_p_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void **) &wfc_c_gpu, sizeof(cufftDoubleComplex) * gSize);
    cudaMalloc((void **) &par_sum, sizeof(cufftDoubleComplex)*(gSize/threads));

    std::cout << "copying contents... " << '\n';

    // Copy variables over to device
    cudaMemcpy(wfc_gpu, wfc, sizeof(cufftDoubleComplex) * gSize,
               cudaMemcpyHostToDevice);
    std::cout << "wfc copied..." << '\n';
    cudaMemcpy(K_0_gpu, K, sizeof(cufftDoubleComplex) * gSize,
               cudaMemcpyHostToDevice);
    std::cout << "K copied..." << '\n';
    cudaMemcpy(V_0_gpu, GV, sizeof(cufftDoubleComplex) * gSize,
               cudaMemcpyHostToDevice);
    std::cout << "V copied..." << '\n';
    cudaMemcpy(wfc_p_gpu, wfc_p, sizeof(cufftDoubleComplex) * gSize,
               cudaMemcpyHostToDevice);
    std::cout << "wfc_p copied..." << '\n';
    cudaMemcpy(wfc_c_gpu, wfc_c, sizeof(cufftDoubleComplex) * gSize,
               cudaMemcpyHostToDevice);
    std::cout << "wfc_c copied..." << '\n';

    std::cout << "performing energy calculations..." << '\n';


    // In the example python code, it was necessary to reshape everything, 
    // But let's see what happens if I don't do that here...

    // FFT for the wfc in momentum-space
    cufftExecZ2Z(plan_2d, wfc_gpu, wfc_p, CUFFT_FORWARD);

    // Conjugate for the wfc
    vecConjugate<<<grid,threads>>>(wfc_gpu, wfc_c);

    // K * wfc
    vecMult<<<grid,threads>>>(wfc_gpu,K_0_gpu,wfc_p);
    cufftExecZ2Z(plan_2d, wfc_p, Energy_1, CUFFT_INVERSE); 

    vecMult<<<grid,threads>>>(wfc_gpu, V_0_gpu, Energy_2);

/*
    for (int i = 0; i < xDim * yDim; ++i){
        std::cout << Energy_1[i].y << '\t' << Energy_2[i].x << '\n';
    }
*/

    //std::cout << wfc_gpu[0].x << '\t' << wfc_gpu[0].y << '\n';

    free(EV); free(EK); free(EpAy); free(EpAx);
    free(x);free(y);
    cudaFree(wfc_gpu); cudaFree(K_gpu); cudaFree(V_gpu); cudaFree(pAx_gpu);
    cudaFree(pAy_gpu); cudaFree(par_sum);

    std::cout << "Evolution test complete." << '\n';
    std::cout << "EVOLUTION TEST UNFINISHED!" << '\n';
    
}

/*
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
*/
