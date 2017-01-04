///@endcond
//##############################################################################
/**
 *  @file    split_op.h
 *  @author  James Ryan Schloss (leios)
 *  @date    5/25/2016
 *  @version 0.1
 *
 *  @brief command line parser file.
 *
 */
//##############################################################################

#ifndef PARSER_H
#define PARSER_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fstream>
#include <sys/stat.h>
#include <algorithm>
#include "../include/ds.h"
#include "../include/unit_test.h"

/**
* @brief	Parses command-line input, creates initial grid
*/
Grid parseArgs(int argc, char** argv);

#endif
