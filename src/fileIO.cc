/*
* fileIO.c - GPUE: Split Operator based GPU solver for Nonlinear 
* Schrodinger Equation, Copyright (C) 2012, Lee J. O'Riordan, Tadhg 
* Morgan, Neil Crowley. 

* This library is free software; you can redistribute it and/or modify 
* it under the terms of the GNU Lesser General Public License as 
* published by the Free Software Foundation; either version 2.1 of the 
* License, or (at your option) any later version. This library is 
* distributed in the hope that it will be useful, but WITHOUT ANY 
* WARRANTY; without even the implied warranty of MERCHANTABILITY or 
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public 
* License for more details. You should have received a copy of the GNU 
* Lesser General Public License along with this library; if not, write 
* to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
* Boston, MA 02111-1307 USA 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime.h>


int readIn(char* fileR, char* fileI, int xDim, int yDim, double2 *arr){
	FILE *f;
	f = fopen(fileR,"r");
	int i = 0;
	double line;
	while(i < xDim*yDim){
		fscanf(f,"%lf",&line);
		arr[i].x = line;
		++i;
	}
	f = fopen(fileI,"r");
	i = 0;
	while(i < xDim*yDim){
		fscanf(f,"%lf",&line);
		arr[i].y = line;
		++i;
	}
	fclose(f);
	return 0;
}


void writeOut(char* buffer, char *file, double2 *data, int length, int step){
	FILE *f;
	sprintf (buffer, "%s_%d", file, step);
	f = fopen (buffer,"w");
	int i;
	for (i = 0; i < length; i++) 
		fprintf (f, "%.16e\n",data[i].x);
	fclose (f);
	
	sprintf (buffer, "%si_%d", file, step);
	f = fopen (buffer,"w");
	for (i = 0; i < length; i++) 
		fprintf (f, "%.16e\n",data[i].y);
	fclose (f);
}

void writeOutDouble(char* buffer, char *file, double *data, int length, int step){
	FILE *f;
	sprintf (buffer, "%s_%d", file, step);
	f = fopen (buffer,"w");
	int i;
	for (i = 0; i < length; i++) 
		fprintf (f, "%.16e\n",data[i]);
	fclose (f);
}


int readState(char* name){
	FILE *f;
	f = fopen(name,"r");
	fclose(f); 
	return 0;
}
