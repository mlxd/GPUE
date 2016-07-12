/*** ds.cc - GPUE: Split Operator based GPU solver for Nonlinear 
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
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSEDvd AND ON ANY THEORY OF 
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "../include/ds.h"

// Function to store integer into Grid->param_int
void Grid::store(std::string id, int iparam){
    param_int[id] = iparam;
}

// Function to store double into Grid->param_double
void Grid::store(std::string id, double dparam){
    param_double[id] = dparam;
}

// Function to store double* into param_dstar
void Grid::store(std::string id, double *dsparam){
    param_dstar[id] = dsparam;
}

// Function to store bool into param_bool
void Grid::store(std::string id, bool bparam){
    param_bool[id] = bparam;
}
// Function to retrieve integer from Grid->param_int
int Grid::ival(std::string id){
    return param_int[id];
}

// Function to retrieve double from Grid->param_double
double Grid::dval(std::string id){
    return param_double[id];
}

// Function to retrieve double star values from param_dstar
double *Grid::dsval(std::string id){
    return param_dstar[id];
}

// Function to retrieve bool values from param_bool
bool Grid::bval(std::string id){
    return param_bool[id];
}

// Function for file writing (to replace writeOutParam)
void Grid::write(std::string filename){
    std::ofstream output;
    output.open(filename);
    // We simply iterate through the int and double param maps
    for (auto item : param_double){
        output << item.first << "=" << item.second << '\n';
    }

    for (auto item : param_int){
        output << item.first << "=" << item.second << '\n';
    }

    output.close();
}

// Functions to store data in Cuda class
void Cuda::store(std::string id, cudaError_t errin){
    err = errin;
}

void Cuda::store(std::string id, cufftResult resultin){
    result = resultin;
}

void Cuda::store(std::string id, cufftHandle planin){
    if (id == "plan_1d"){
        plan_1d = planin;
    }
    else if (id == "plan_2d"){
        plan_2d = planin;
    }
    else{
        std::cout << "Error: plan not found!" << '\n';
    }
}

void Cuda::store(std::string id, cudaStream_t streamin){
    if (id == "streamA"){
        streamA = streamin;
    }
    else if (id == "streamB"){
        streamB = streamin;
    }
    else if (id == "streamC"){
        streamC = streamin;
    }
    else if (id == "streamD"){
        streamD = streamin;
    }
    else{
        std::cout << "Error: stream not found!" << '\n';
    }
}

void Cuda::store(std::string id, dim3 gridin){
    grid = gridin;
}

// Functions to retrieve data from Cuda class
// Note: There are definitely more elegant ways to do this.
cudaError_t Cuda::cudaError_tval(std::string id){
    return err;
}

cufftResult Cuda::cufftResultval(std::string id){
    return result;
}

// Returns nothing if called incorrectly.
cufftHandle Cuda::cufftHandleval(std::string id){
    if (id == "plan_1d"){
        return plan_1d;
    }
    else if (id == "plan_2d"){
        return plan_2d;
    }
    else{
        std::cout << "Error: plan not found!" << '\n';
        exit(1);
    }
}

// Returns nothing if called incorrectly
cudaStream_t Cuda::cudaStream_tval(std::string id){
    if (id == "streamA"){
        return streamA;
    }
    else if (id == "streamB"){
        return streamB;
    }
    else if (id == "streamC"){
        return streamC;
    }
    else if (id == "streamD"){
        return streamD;
    }
    else{
        std::cout << "Error: stream not found!" << '\n';
        exit(1);
    }
}

dim3 Cuda::dim3val(std::string id){
    return grid;
}

// Functions to store data in the Op class
void Op::store(std::string id, double *data){
    Op_dstar[id] = data;
}

void Op::store(std::string id, cufftDoubleComplex *data){
    Op_cdc[id] = data;
}

// Functions to retrieve data from the Op class
double *Op::dsval(std::string id){
    return Op_dstar[id];
}

cufftDoubleComplex *Op::cufftDoubleComplexval(std::string id){
    return Op_cdc[id];
}

// Functions to store data in the Wave class
void Wave::store(std::string id, double *data){
    Wave_dstar[id] = data;
}

void Wave::store(std::string id, cufftDoubleComplex *data){
    Wave_cdc[id] = data;
}

// Functions to retrieve data from the Wave class
double *Wave::dsval(std::string id){
    return Wave_dstar[id];
}

cufftDoubleComplex *Wave::cufftDoubleComplexval(std::string id){
    return Wave_cdc[id];
}
