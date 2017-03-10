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


#include "params.h"

using namespace std;

void read_params( std::string params_name )
{
 for ( int i = 0; i < MAXRES; i++ )
 {
  aa_comp[i] = 0.0;
  
  for ( int j = 0; j < 2; j++ )
   p_pot[i][j] = 0.0;
  
  for ( int j = 0; j < MAXBUR; j++ )
   p_bur[j][i] = 0.0;
  
  for ( int j = 0; j < MAXSTR; j++ )
   p_str[j][i] = 0.0;
 }
 
 for ( int i = 0; i < MAXCMP; i++ )
  for ( int j = 0; j < 3; j++ )
   weights[i][j] = 0.0;
 
 string line1;
 
 ifstream params_file( params_name.c_str() );
 
 if ( !params_file.is_open() ) { cout << "Cannot open " << params_name << endl; exit(EXIT_FAILURE); }
 
 while (getline(params_file,line1))
  if ( line1.length() > 5 )
  {
   if ( line1.substr(0,4) == "WGHT")
   {
    weights[atoi(line1.substr(5,1).c_str())][0] = atof(line1.substr(6,5).c_str());
    weights[atoi(line1.substr(5,1).c_str())][1] = atof(line1.substr(11,5).c_str());
    weights[atoi(line1.substr(5,1).c_str())][2] = atof(line1.substr(16,9).c_str());
   }
   
   else if ( line1.substr(0,4) == "COMP")
    aa_comp[one2number(line1.substr(5,1).c_str())] = atof(line1.substr(6,9).c_str());
   
   else if ( line1.substr(0,4) == "AAFQ")
    aa_freq[one2number(line1.substr(5,1).c_str())] = atof(line1.substr(6,7).c_str());
   
   else if ( line1.substr(0,4) == "POTS")
   {
    p_pot[one2number(line1.substr(5,1).c_str())][0] = atof(line1.substr(6,10).c_str());
    p_pot[one2number(line1.substr(5,1).c_str())][1] = atof(line1.substr(16,9).c_str());
   }
   
   else if ( line1.substr(0,4) == "PBUR")
    p_bur[burial2number(line1.substr(5,1).c_str())][one2number(line1.substr(7,1).c_str())] = atof(line1.substr(8,8).c_str());
   
   else if ( line1.substr(0,4) == "SSTR")
    p_str[stride2number(line1.substr(5,1).c_str())][one2number(line1.substr(7,1).c_str())] = atof(line1.substr(8,8).c_str());
   
   else if ( line1.substr(0,4) == "FREQ")
    p_frq[atoi(line1.substr(5,1).c_str())] = atof(line1.substr(6,9).c_str());
  }
 
 params_file.close();
}
