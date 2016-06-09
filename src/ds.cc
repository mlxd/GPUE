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
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "../include/ds.h"

// Unfortunately, I only know how to do this as a macro...
// Turns variable name into string
#define STR(a) #a

// Function to store integer into Grid->param_int
void Grid::store(std::string id, int iparam){
    param_int[id] = iparam;
    id_list[id] = 1;
}

// Function to store double into Grid->param_double
void Grid::store(std::string id, double dparam){
    param_double[id] = dparam;
    id_list[id] = 2;
}

// Function to store double* into param_dstar
void Grid::store(std::string id, double *dsparam){
    param_dstar[id] = dsparam;
    id_list[id] = 3;
}

/*
// Defining template and function to return values from parameter maps.
template <typename arbitrary> 
arbitrary Grid::val(std::string id){
    switch (id_list[id]){
        case 1:
            return param_int[id];
            break;
        case 2:
            return param_double[id];
            break;
        case 3:
            return param_dstar[id];
            break;
    }
}
*/

// Function to retrieve integer from Grid->param_int
int Grid::ival(std::string id){
    return param_int[id];
}

// Function to retrieve double from Grid->param_double
double Grid::dval(std::string id){
    return param_double[id];
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

// Failed template attempt for the Cuda class. I spent too much time on it, may
// return to it later
/*
// Not necessary because everything is public for now...
// Function to evaluate arbitrary value from Cuda class
template <typename arbitrary = int> 
arbitrary Cuda::val(std::string id){
    switch (id){
        case "err":
            return err;
            break;
        case "result":
            return result;
            break;
        case "plan_1d":
            return plan_1d;
            break;
        case "plan_2d":
            return plan_2d;
            break;
        case "streamA":
            return streamA;
            break;
        case "streamB":
            return streamB;
            break;
        case "streamC":
            return streamC;
            break;
        case "streamD":
            return streamD;
            break;
        case "grid":
            return grid;
            break;
        case "threads":
            return threads;
            break;
    }
}

template <class arbitrary> arbitrary Cuda::store(arbitrary data){
    this->data = data;
}
*/

// Functions to store data in Cuda class
void Cuda::store(std::string id, cudaError_t errin){
    err = errin;
}
void Cuda::store(std::string id, cufftHandle planin){
    if (id == "plan_1d"){
        plan_1d = planin;
    }
    else if (id == "plan_2d"){
        plan_2d = planin;
    }
    else{
        std::cout << "ERROR: plan not found!" << '\n';
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
        std::cout << "ERROR: stream not found!" << '\n';
    }
}
void Cuda::store(std::string id, dim3 gridin){
    grid = gridin;
}
/*
void Cuda::store(std::string id, int threadsin){
    threads = threadsin;
}
*/

// Functions to retrieve data
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
        std::cout << "ERROR: plan not found!" << '\n';
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
        std::cout << "ERROR: stream not found!" << '\n';
    }
}
dim3 Cuda::dim3val(std::string id){
    return grid;
}
int Cuda::ival(std::string id){
    return threads;
}
/*----------------------------------------------------------------------------//
* DEPRECATION WARNING
*-----------------------------------------------------------------------------*/
// 
// void initArr(Array *arr, size_t initLen){
// 	arr->array = (Param*) malloc(initLen*sizeof(Param));
// 	arr->used = 0;
// 	arr->length = initLen;
// }
// 
// void appendData(Array *arr, std::string t, double d){
// 	Param p = newParam(t,d);
// 	if(arr->used == arr->length){
// 		arr->length *= 2;
// 		arr->array = (Param*)realloc(arr->array, arr->length*sizeof(Param));
// 	}
// 	arr->array[arr->used] = p;
// 	arr->used = arr->used + 1;
// }
// 
// void freeArray(Array *arr){
// 	free(arr->array);
// 	arr->array = NULL;
// 	arr->used = 0;
// 	arr->length = 0;
// }
// 
// Param newParam(std::string t,double d){
// 	Param p;
// 	strcpy(p.title,t.c_str());
// 	p.data = d;
// 	return p;
// }
