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

// Function to store integer into Grid->param_int
void Grid::store(std::string id, int iparam){
    param_int[id] = iparam;
}

// Function to store double into Grid->param_double
void Grid::store(std::string id, double dparam){
    param_double[id] = dparam;
}

// Function to retrieve integer from Grid->param_int
int Grid::ival(std::string id){
    return param_int[id];
}

//Function to retrieve double from Grid->param_double
double Grid::dval(std::string id){
    return param_double[id];
}

/*----------------------------------------------------------------------------//
* DEPRECATION WARNING
*-----------------------------------------------------------------------------*/

void initArr(Array *arr, size_t initLen){
	arr->array = (Param*) malloc(initLen*sizeof(Param));
	arr->used = 0;
	arr->length = initLen;
}

void appendData(Array *arr, std::string t, double d){
	Param p = newParam(t,d);
	if(arr->used == arr->length){
		arr->length *= 2;
		arr->array = (Param*)realloc(arr->array, arr->length*sizeof(Param));
	}
	arr->array[arr->used] = p;
	arr->used = arr->used + 1;
}

void freeArray(Array *arr){
	free(arr->array);
	arr->array = NULL;
	arr->used = 0;
	arr->length = 0;
}

Param newParam(std::string t,double d){
	Param p;
	strcpy(p.title,t.c_str());
	p.data = d;
	return p;
}
