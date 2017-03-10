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


#ifndef __DFIRE_H_
#define __DFIRE_H_

#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "size.h"
#include "data.h"

using namespace std;

extern int dfire1[MAXDF1];
extern double dfire2[MAXDF4][MAXRES][MAXRES];
extern int dfire3[MAXDF2];
extern double dfire4[MAXDF3][MAXRES][MAXRES];

void read_potential_sc( std::string );

void read_potential_ca( std::string );

#endif
