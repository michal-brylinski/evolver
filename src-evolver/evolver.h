/*
===============================================================================
           ___ ___         __                   
   .-----.|   |   |.-----.|  |.--.--.-----.----.
   |  -__||   |   ||  _  ||  ||  |  |  -__|   _|
   |_____| \_____/ |_____||__| \___/|_____|__|  
                                                  
   eVolver - protein sequence generator

   Computational Systems Biology Group
   Department of Biological Sciences
   Center for Computation & Technology
   Louisiana State University
   407 Choppin Hall, Baton Rouge, LA 70803, USA

   http://www.brylinski.org

   Report bugs to michal@brylinski.org

   Copyright 2013 Michal Brylinski

   This file is part of eVolver.

   eVolver is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   eVolver is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with eVolver. If not, see <http://www.gnu.org/licenses/>.

===============================================================================
*/


#ifndef __EVOLVER_H_
#define __EVOLVER_H_

#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <ctime>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_siman.h>
#include <gsl/gsl_ieee_utils.h>

#include "size.h"
#include "data.h"
#include "dfire.h"
#include "params.h"
#include "parser.h"
#include "energy.h"
#include "siman.h"
#include "shuffle.h"
#include "output.h"

#define N_TRIES 200
#define ITERS_FIXED_T 2000
#define STEP_SIZE 1.0
#define K 1.0
#define T_INITIAL 5000.0
#define MU_T 1.002
#define T_MIN 5.0e-3

gsl_siman_params_t params = { N_TRIES, ITERS_FIXED_T, STEP_SIZE, K, T_INITIAL, MU_T, T_MIN };

double aa_comp[MAXRES];
double aa_freq[MAXRES];
double p_pot[MAXRES][2];
double p_bur[MAXBUR][MAXRES];
double p_str[MAXSTR][MAXRES];
double p_frq[MAXSTR];

bool profseq_opt;

int verb_opt;

int dfire1[MAXDF1];
double dfire2[MAXDF4][MAXRES][MAXRES];
int dfire3[MAXDF2];
double dfire4[MAXDF3][MAXRES][MAXRES];

int res;

int bin_ca[MAXPRO][MAXPRO];
int bin_sc[MAXPRO][MAXPRO];

int sequence[MAXPRO];
int stride[MAXPRO];
int burial[MAXPRO];

double profstr[MAXPRO][MAXSTR];
double profseq[MAXPRO][MAXRES];

double weights[MAXCMP][3];

vector<std::string> pdbcontent;

double dif1;

#endif
