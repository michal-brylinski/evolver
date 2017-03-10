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


#include "energy.h"

double get_energy_burial( int t_seq[] )
{
 double ene1 = 0.0;
 
 for ( int i = 0; i < res; i++ )
  ene1 += p_bur[burial[i]][t_seq[i]];
 
 ene1 /= double ( res );
 
 return ene1;
}

double get_energy_stride( int t_seq[] )
{
 double ene1 = 0.0;
 
 for ( int i = 0; i < res; i++ )
  ene1 += p_str[stride[i]][t_seq[i]];
 
 ene1 /= double ( res );
 
 return ene1;
}

double get_energy_profile( int t_seq[] )
{
 double ene1 = 0.0;
 double ene2 = 0.0;
 
 for ( int i = 0; i < res; i++ )
 {
  ene1 += profstr[i][number2pinak(t_seq[i])];
  
  double ene3 = 0.0;
  
  for ( int j = 0; j < MAXSTR; j++ )
   ene3 += profstr[i][j];
  
  if ( ene3 > 0.0 )
   ene2 += 1.0;
 }
 
 ene1 /= ene2;
 
 return ene1;
}

double get_energy_dfire_ca( int t_seq[] )
{
 double ene1 = 0.0;
 
 for ( int i = 0; i < res - 1; i++ )
  for ( int j = i + 1; j < res; j++ )
   ene1 -= dfire4[ bin_ca[i][j] ][ t_seq[i] ][ t_seq[j] ];
 
 ene1 = ene1 / ( 1.8031 * res - 18.669 );
 
 return ene1;
}

double get_energy_dfire_sc( int t_seq[] )
{
 double ene1 = 0.0;
 
 for ( int i = 0; i < res - 1; i++ )
  for ( int j = i + 1; j < res; j++ )
   ene1 -= dfire2[ bin_sc[i][j] ][ t_seq[i] ][ t_seq[j] ];
 
 ene1 = ene1 / ( 1.5766 * res - 44.863 );
 
 return ene1;
}

double get_energy_pot( int t_seq[] )
{
 double ene1 = 0.0;
 double ene2 = 0.0;
 
 for ( int i = 0; i < MAXRES; i++ )
 {
  int n1[MAXPRO];
  
  int n2 = 0;
  
  for ( int j = 0; j < res; j++ )
   if ( t_seq[j] == i )
    n1[n2++] = j;
  
  if ( n2 > 2 )
  {
   double m1 = res;
   
   double r1 = exp( ( -1.0 * n2 ) / m1 );
   
   double pot1 = 0.0;
   
   for ( int k = 0; k < n2 - 1; k++ )
    for ( int l = k + 1; l < n2; l++ )
     pot1 += exp(  -1.0 * abs( n1[k] - n1[l] ) * ( n2 / m1 ) );
   
   double pot0 = n2 * ( n2 - 1.0 ) * r1 * ( m1 - 1.0 / ( 1.0 - r1 ) ) / ( m1 * ( m1 - 1.0 ) * ( 1.0 - r1 ) );
   
   double pot2 = pot1 - pot0;
   
   double dev1 = sqrt( ( r1 * r1 * n2 * n2 * ( 1.0 - n2 / m1 ) * ( 1.0 - n2 / m1 ) ) / ( m1 * ( 1.0 - r1 * r1 ) ) );
   
   double z1 = 0.0;
   
   if ( dev1 != 0 )
    z1 = pot2 / dev1;
   
   ene1 += ( 0.5 * gsl_sf_pow_int( ( z1 - p_pot[i][0] ) / p_pot[i][1], 2 ) - log( 1.0 / ( p_pot[i][1] * sqrt( 2.0 * PI ) ) ) );
   
   ene2 += 1.0;
  }
 }
 
 ene1 /= ene2;
 
 return ene1;
}

double get_profile_score( int t_seq[] )
{
 double sco1 = 0.0;
 
 for ( int i = 0; i < res; i++ )
  sco1 += profseq[i][t_seq[i]];
 
 sco1 /= double ( res );
 
 return sco1;
}

double get_sequence_score( int t_seq[] )
{
 double sco1 = 0.0;
 
 for ( int i = 0; i < res; i++ )
  if ( t_seq[i] == sequence[i] )
   sco1++;
 
 sco1 /= double ( res );
 
 return sco1;
}

void Enormalize( double & v, double lb, double ub )
{
 v = ( ( v - lb ) / ( ub - lb ) );
 
 if ( v < 0.0 )
  v = 0.0;
 if ( v > 1.0 )
  v = 1.0;
}
