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


#include "shuffle.h"

using namespace std;

void seq_shuffle( int t_seq[] )
{
 const gsl_rng * r = gsl_rng_alloc (gsl_rng_env_setup());
 
 gsl_ran_shuffle (r, t_seq, res, sizeof (int));
}

void seq_random( int t_seq[] )
{
 int r_seq[MAXPRO];
 
 int n2 = 0;
 
 for ( int i = 0; i < MAXRES; i++ )
 {
  int n1 = int ( aa_freq[i] * double ( res ) );
  
  for ( int j = 0; j <= n1; j++ )
   r_seq[n2++] = i;
 }
 
 const gsl_rng * r = gsl_rng_alloc (gsl_rng_env_setup());
 
 gsl_ran_shuffle (r, r_seq, n2, sizeof (int));
 
 for ( int i = 0; i < res; i++ )
  t_seq[i] = r_seq[i];
}
