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
#include "../include/operators.h"

/*----------------------------------------------------------------------------//
* GRID
*-----------------------------------------------------------------------------*/

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

// Function to store string into data_dir
void Grid::store(std::string id, std::string sparam){
    param_string[id] = sparam;
}

// Two boolean functions to check whether a string exists in 
// param_double or param_dstar
bool Grid::is_double(std::string id){
    auto it = param_double.find(id);
    if (it != param_double.end()){
        return true;
    }
    else {
        return false;
    }
}

bool Grid::is_dstar(std::string id){
    auto it = param_dstar.find(id);
    if (it != param_dstar.end()){
        return true;
    }
    else {
        return false;
    }
}

// Function to retrieve integer from Grid->param_int
int Grid::ival(std::string id){
    return param_int[id];
}

// Function to retrieve double from Grid->param_double
double Grid::dval(std::string id){
    auto it = param_double.find(id);
    if (it == param_double.end()){
        std::cout << "ERROR: could not find string " << id 
                  << " in Grid::param_double." << '\n';
        assert(it != param_double.end());
    }
    return it->second;
}

// Function to retrieve double star values from param_dstar
double *Grid::dsval(std::string id){
    auto it = param_dstar.find(id);
    if (it == param_dstar.end()){
        std::cout << "ERROR: could not find string " << id
                  << " in Grid::param_dstar." << '\n';
        assert(it != param_dstar.end());
    }
    return it->second;
}

// Function to retrieve bool values from param_bool
bool Grid::bval(std::string id){
    auto it = param_bool.find(id);
    if (it == param_bool.end()){
        std::cout << "ERROR: could not find string " << id 
                  << " in Grid::param_bool." << '\n';
        assert(it != param_bool.end());
    }
    return it->second;
}

// Function to retrieve string from data_dir
// Note: There is only one string value in the Grid struct... 
//       We must add an unordered map for strings if further strings are desired
std::string Grid::sval(std::string id){
    auto it = param_string.find(id);
    if (it == param_string.end()){
        std::cout << "ERROR: could not find string " << id 
                  << " in Grid::param_string." << '\n';
        assert(it != param_string.end());
    }
    return it->second;
}

// Function for file writing (to replace fileIO::writeOutParam)
void Grid::write(std::string filename){
    std::ofstream output;
    output.open(filename);
    // We simply iterate through the int and double param maps
    for (auto item : param_double){
        output << item.first << "=" << item.second << '\n';
        std::cout << item.first << "=" << item.second << '\n';
    }

    for (auto item : param_int){
        output << item.first << "=" << item.second << '\n';
        std::cout << item.first << "=" << item.second << '\n';
    }

    output.close();
}

// Function to print all available variables
void Grid::print_map(){
    for (auto item : param_double){
        std::cout << item.first << '\n';
    }
    for (auto item : param_dstar){
       std::cout << item.first << '\n';
    }
}

/*----------------------------------------------------------------------------//
* CUDA
*-----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------//
* OP
*-----------------------------------------------------------------------------*/

// Functions to store data in the Op class
void Op::store(std::string id, double *data){
    Op_dstar[id] = data;
}

void Op::store(std::string id, cufftDoubleComplex *data){
    Op_cdc[id] = data;
}

// Functions to retrieve data from the Op class
double *Op::dsval(std::string id){
    auto it = Op_dstar.find(id);
    if (it == Op_dstar.end()){
        std::cout << "ERROR: could not find string " << id 
                  << " in Op::Op_dstar." << '\n';
        assert(it != Op_dstar.end());
    }
    return it->second;
}

cufftDoubleComplex *Op::cufftDoubleComplexval(std::string id){
    auto it = Op_cdc.find(id);
    if (it == Op_cdc.end()){
        std::cout << "ERROR: could not find string " << id 
                  << " in Op::Op_cdc." << '\n';
        assert(it != Op_cdc.end());
    }
    return it->second;
}

// Map for function pointers and keys K and V
Op::functionPtr Op::K_fn(std::string id){
    auto it = Op_K_fns.find(id);
    if (it == Op_K_fns.end()){
        std::cout << "ERROR: could not find string " << id 
                  << " in Op::Op_K_fns. Did you mean: " << '\n';
        for (auto item : Op_K_fns){
            std::cout << item.first << '\n';
        }
    }
    return it->second;
}

Op::functionPtr Op::V_fn(std::string id){
    auto it = Op_V_fns.find(id);
    if (it == Op_V_fns.end()){
        std::cout << "ERROR: could not find string " << id 
                  << " in Op::Op_V_fns. Did you mean: " << '\n';
        for (auto item : Op_V_fns){
            std::cout << item.first << '\n';
        }
    }
    return it->second;
}

// Map for function pointers for pAx, pAy, and pAz
Op::functionPtr Op::pAx_fn(std::string id){
    auto it = Op_pAx_fns.find(id);
    if (it == Op_pAx_fns.end()){
        std::cout << "ERROR: could not find string " << id 
                  << " in Op::Op_pAx_fns. Did you mean: " << '\n';
        for (auto item : Op_pAx_fns){
            std::cout << item.first << '\n';
        }
    }
    return it->second;
}

Op::functionPtr Op::pAy_fn(std::string id){
    auto it = Op_pAy_fns.find(id);
    if (it == Op_pAy_fns.end()){
        std::cout << "ERROR: could not find string " << id 
                  << " in Op::Op_pAy_fns. Did you mean: " << '\n';
        for (auto item : Op_pAy_fns){
            std::cout << item.first << '\n';
        }
    }
    return it->second;
}

Op::functionPtr Op::pAz_fn(std::string id){
    auto it = Op_pAz_fns.find(id);
    if (it == Op_pAz_fns.end()){
        std::cout << "ERROR: could not find string " << id 
                  << " in Op::Op_pAz_fns. Did you mean: " << '\n';
        for (auto item : Op_pAz_fns){
            std::cout << item.first << '\n';
        }
    }
    return it->second;
}

Op::functionPtr Op::Ax_fn(std::string id){
    auto it = Op_Ax_fns.find(id);
    if (it == Op_Ax_fns.end()){
        std::cout << "ERROR: could not find string " << id 
                  << " in Op::Op_Ax_fns. Did you mean: " << '\n';
        for (auto item : Op_Ax_fns){
            std::cout << item.first << '\n';
        }
    }
    return it->second;
}

Op::functionPtr Op::Ay_fn(std::string id){
    auto it = Op_Ay_fns.find(id);
    if (it == Op_Ay_fns.end()){
        std::cout << "ERROR: could not find string " << id 
                  << " in Op::Op_Ay_fns. Did you mean: " << '\n';
        for (auto item : Op_Ay_fns){
            std::cout << item.first << '\n';
        }
    }
    return it->second;
}

Op::functionPtr Op::Az_fn(std::string id){
    auto it = Op_Az_fns.find(id);
    if (it == Op_Az_fns.end()){
        std::cout << "ERROR: could not find string " << id 
                  << " in Op::Op_Az_fns. Did you mean: " << '\n';
        for (auto item : Op_Az_fns){
            std::cout << item.first << '\n';
        }
    }
    return it->second;
}
// Function to set functionPtrs for K and V
// Unfortunately, these must be set one at a time.
void Op::set_fns(){

    Op_K_fns["rotation_K"] = rotation_K;
    Op_K_fns["rotation_gauge_K"] = rotation_gauge_K;
    Op_V_fns["harmonic_V"] = harmonic_V;
    Op_V_fns["harmonic_gauge_V"] = harmonic_gauge_V;
    Op_pAy_fns["rotation"] = rotation_pAy;
    Op_pAx_fns["rotation"] = rotation_pAx;
    Op_Ax_fns["rotation"] = rotation_Ax;
    Op_Ay_fns["rotation"] = rotation_Ay;
    Op_Ax_fns["dynamic"] = dynamic_Ax;
    Op_Ay_fns["dynamic"] = dynamic_Ay;
    Op_pAx_fns["rotation_squared"] = rotation_Ax;
    Op_pAy_fns["rotation_squared"] = rotation_Ay;
    Op_Ay_fns["rotation_squared"] = rotation_squared_Ay;
    Op_Ax_fns["rotation_squared"] = rotation_squared_Ax;
    

    //Op_K_fns.emplace("rotation_K", rotation_K);
    //Op_V_fns.emplace("harmonic_V", harmonic_V);
    
}

/*----------------------------------------------------------------------------//
* WAVE
*-----------------------------------------------------------------------------*/

// Functions to store data in the Wave class
void Wave::store(std::string id, double *data){
    Wave_dstar[id] = data;
}

void Wave::store(std::string id, cufftDoubleComplex *data){
    Wave_cdc[id] = data;
}

// Functions to retrieve data from the Wave class
double *Wave::dsval(std::string id){
    auto it = Wave_dstar.find(id);
    if (it == Wave_dstar.end()){
        std::cout << "ERROR: could not find string " << id 
                  << " in Wave::Wave_dstar." << '\n';
        assert(it != Wave_dstar.end());
    }
    return it->second;
}

cufftDoubleComplex *Wave::cufftDoubleComplexval(std::string id){
    auto it = Wave_cdc.find(id);
    if (it == Wave_cdc.end()){
        std::cout << "ERROR: could not find string " << id 
                  << " in Wave::Wave_cdc." << '\n';
        assert(it != Wave_cdc.end());
    }
    return it->second;
}

/*----------------------------------------------------------------------------//
* MISC
*-----------------------------------------------------------------------------*/

/*
// Template function to print all values in map
template <typename T> void print_map(std::unordered_map<std::string, T> map){
    std::cout << "Contents of map are: " << '\n';
    std::cout << "key: " << '\t' << "element: " << '\n';
    for (auto element : map){
        std::cout << element.first << '\t' << element.second << '\n';
    }
}
*/
