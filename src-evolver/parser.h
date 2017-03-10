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


#ifndef __PARSER_H_
#define __PARSER_H_

#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <algorithm>

#include "size.h"
#include "data.h"

using namespace std;

extern double aa_freq[MAXRES];

extern int dfire1[MAXDF1];
extern int dfire3[MAXDF2];

extern int res;

extern int bin_ca[MAXPRO][MAXPRO];
extern int bin_sc[MAXPRO][MAXPRO];

extern int sequence[MAXPRO];
extern int stride[MAXPRO];
extern int burial[MAXPRO];

extern double profstr[MAXPRO][MAXSTR];
extern double profseq[MAXPRO][MAXRES];

extern vector<std::string> pdbcontent;

void read_pdb( std::string );

void read_stride( std::string );

void read_profstr( std::string );

void read_profseq( std::string );

#endif
