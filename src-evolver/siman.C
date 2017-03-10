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


#include "siman.h"

using namespace std;

double Etsp( void *xp )
{
 int *seq_t = (int *) xp;
 
 double erg3[MAXCMP];
 
 erg3[0] = get_energy_burial(seq_t);
 erg3[1] = get_energy_stride(seq_t);
 erg3[2] = get_energy_profile(seq_t);
 erg3[3] = get_energy_dfire_ca(seq_t);
 erg3[4] = get_energy_dfire_sc(seq_t);
 erg3[5] = get_energy_pot(seq_t);
 
 for ( int i = 0; i < MAXCMP; i++ )
  Enormalize( erg3[i], weights[i][0], weights[i][1] );
 
 erg3[5] = 1.0 - erg3[5];
 
 double erg4 = 0.0;
 
 for ( int i = 0; i < MAXCMP; i++ )
 {
  erg3[i] *= weights[i][2];
  
  erg4 += erg3[i];
 }
 
 dif1++;
 
 return 1000.0 / erg4;
}

double Mtsp( void *xp, void *yp )
{
 int *seq_t1 = (int *) xp;
 int *seq_t2 = (int *) yp;
 
 double dist = 0;
 
 for ( int i = 0; i < res; i++ )
  dist += ( ( seq_t1[i] == seq_t2[i] ) ? 0 : 1 );
 
 return dist;
}

void Stsp( const gsl_rng * r, void *xp, double step_size )
{
 int x1, x2, dummy;
 int *seq_t = (int *) xp;
 
 step_size = 0;
 
 x1 = (gsl_rng_get (r) % res);
 do {
   x2 = (gsl_rng_get (r) % res);
 } while (x2 == x1);
 
 dummy = seq_t[x1];
 seq_t[x1] = seq_t[x2];
 seq_t[x2] = dummy;
}

void Ptsp( void *xp )
{
 if ( verb_opt == 1 )
  cout << " STEP";
 
 else if ( verb_opt == 2 )
 {
  int *seq_t = (int *) xp;
  
  string seq2 = "";
  
  for ( int i = 0; i < res; i++ )
   seq2 += number2one(seq_t[i]);
  
  double erg5[MAXCMP];
  
  erg5[0] = get_energy_burial(seq_t);
  erg5[1] = get_energy_stride(seq_t);
  erg5[2] = get_energy_profile(seq_t);
  erg5[3] = get_energy_dfire_ca(seq_t);
  erg5[4] = get_energy_dfire_sc(seq_t);
  erg5[5] = get_energy_pot(seq_t);
  
  for ( int i = 0; i < MAXCMP; i++ )
   Enormalize( erg5[i], weights[i][0], weights[i][1] );
  
  erg5[5] = 1.0 - erg5[5];
  
  double erg6 = 0.0;
  
  for ( int i = 0; i < MAXCMP; i++ )
  {
   erg5[i] *= weights[i][2];
   
   erg6 += erg5[i];
  }
  
  if ( profseq_opt )
   cout << " STEP " << get_profile_score(seq_t) << " " << get_sequence_score(seq_t) << " " << erg6 << " " << seq2;
  
  else
   cout << " STEP " << get_sequence_score(seq_t) << " " << erg6 << " " << seq2;
 }
}
