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


#include "parser.h"

using namespace std;

void read_pdb( std::string pdb_name )
{
 string line1;
 
 res = 0;
 
 double q_ca[MAXPRO][3];
 double q_cb[MAXPRO][3];
 double q_sc[MAXPRO][3];
 
 double q_ms1[MAXPRO][4];
 
 for ( int i = 0; i < MAXPRO; i++ )
  for ( int j = 0; j < 4; j++ )
   q_ms1[i][j] = 0.0;
 
 ifstream pdb_file( pdb_name.c_str() );
 
 if ( !pdb_file.is_open() ) { cout << "Cannot open " << pdb_name << endl; exit(EXIT_FAILURE); }
 
 while (getline(pdb_file,line1))
  if ( line1.length() > 53 )
   if ( line1.substr(0,6) == "ATOM  " )
   {
    int residue1 = three2number(line1.substr(17,3));
    int residue2 = atoi(line1.substr(22,4).c_str());
    
    string atom1 = line1.substr(12,4);
    string atom2 = line1.substr(13,1);
    
    double x1 = atof(line1.substr(30,8).c_str());
    double y1 = atof(line1.substr(38,8).c_str());
    double z1 = atof(line1.substr(46,8).c_str());
    
    if ( atom1 == " CA " )
    {
     q_ca[residue2-1][0] = x1;
     q_ca[residue2-1][1] = y1;
     q_ca[residue2-1][2] = z1;
     
     sequence[residue2-1] = residue1;
     
     res++;
     
     if ( residue1 == 5 )
     {
      q_sc[residue2-1][0] = x1;
      q_sc[residue2-1][1] = y1;
      q_sc[residue2-1][2] = z1;
      
      q_cb[residue2-1][0] = x1;
      q_cb[residue2-1][1] = y1;
      q_cb[residue2-1][2] = z1;
     }
    }
    else if ( atom1 != " N  " && atom1 != " CA " && atom1 != " C  " && atom1 != " O  " )
    {
     double q_ms2 = 0.0;
     
          if ( atom2 == "C" ) { q_ms2 = 12.010700; }
     else if ( atom2 == "N" ) { q_ms2 = 14.006700; }
     else if ( atom2 == "O" ) { q_ms2 = 15.999400; }
     else if ( atom2 == "P" ) { q_ms2 = 30.973761; }
     else if ( atom2 == "S" ) { q_ms2 = 32.065000; }
     
     q_ms1[residue2-1][0] += q_ms2 * x1;
     q_ms1[residue2-1][1] += q_ms2 * y1;
     q_ms1[residue2-1][2] += q_ms2 * z1;
     
     q_ms1[residue2-1][3] += q_ms2;
     
     if ( atom1 == " CB " )
     {
      q_cb[residue2-1][0] = x1;
      q_cb[residue2-1][1] = y1;
      q_cb[residue2-1][2] = z1;
     }
    }
    
    if ( line1.substr(12,4) == " N  " || line1.substr(12,4) == " CA " || line1.substr(12,4) == " C  " || line1.substr(12,4) == " O  " )
     pdbcontent.push_back( line1.substr(0,54) );
   }
 
 pdb_file.close();
 
 for ( int i = 0; i < res; i++ )
  if ( sequence[i] != 5 )
   for ( int j = 0; j < 3; j++ )
    q_sc[i][j] = q_ms1[i][j] / q_ms1[i][3];
 
 for ( int i = 0; i < res; i++ )
 {
  int k = 0;
  
  for ( int j = 0; j < res; j++ )
  {
   double r1 = sqrt( pow( q_ca[i][0] - q_ca[j][0], 2.0 ) + pow( q_ca[i][1] - q_ca[j][1], 2.0 ) + pow( q_ca[i][2] - q_ca[j][2], 2.0 ) );
   double r2 = sqrt( pow( q_sc[i][0] - q_sc[j][0], 2.0 ) + pow( q_sc[i][1] - q_sc[j][1], 2.0 ) + pow( q_sc[i][2] - q_sc[j][2], 2.0 ) );
   
   int r3 = int ( 2.0 * r1 ) - 1;
   int r4 = int ( 2.0 * r2 ) - 1;
   
        if ( r3 < 0 )          { r3 = 0; }
   else if ( r3 > MAXDF2 - 1 ) { r3 = MAXDF2 - 1; }
   
        if ( r4 < 0 )          { r4 = 0; }
   else if ( r4 > MAXDF1 - 1 ) { r4 = MAXDF1 - 1; }
   
   bin_ca[i][j] = dfire3[r3];
   bin_sc[i][j] = dfire1[r4];
   
   if ( bin_ca[i][j] > MAXDF3 - 1 )
    bin_ca[i][j] = MAXDF3 - 1;
   
   if ( bin_sc[i][j] > MAXDF4 - 1 )
    bin_sc[i][j] = MAXDF4 - 1;
   
   if ( i != j )
   {
    double r5 = sqrt( pow( q_cb[i][0] - q_cb[j][0], 2.0 ) + pow( q_cb[i][1] - q_cb[j][1], 2.0 ) + pow( q_cb[i][2] - q_cb[j][2], 2.0 ) );
    
    if ( r5 <= 14.0 )
     k++;
   }
  }
  
  if ( k < 27 )
   burial[i] = 0;
  
  else if ( k < 34 )
   burial[i] = 1;
  
  else if ( k < 40 )
   burial[i] = 2;
  
  else if ( k < 47 )
   burial[i] = 3;
  
  else if ( k < 55 )
   burial[i] = 4;
  
  else if ( k < 66 )
   burial[i] = 5;
  
  else
   burial[i] = 6;
 }
}

void read_stride( std::string sec_name )
{
 string line1;
 
 ifstream sec_file( sec_name.c_str() );
 
 if ( !sec_file.is_open() ) { cout << "Cannot open " << sec_name << endl; exit(EXIT_FAILURE); }
 
 while (getline(sec_file,line1))
  if ( line1.length() > 39 )
   if ( line1.substr(0,4) == "ASG ")
   {
    string sec1 = line1.substr(24,1);
    
    std::transform(sec1.begin(), sec1.end(), sec1.begin(), ::toupper);
    
    int sec2 = atoi(line1.substr(10,5).c_str()) - 1;
    
    stride[sec2] = stride2number(sec1);
   }
 
 sec_file.close();
}

void read_profstr( std::string prf_name )
{
 string line1;
 
 ifstream prf_file( prf_name.c_str() );
 
 if ( !prf_file.is_open() ) { cout << "Cannot open " << prf_name << endl; exit(EXIT_FAILURE); }
 
 while (getline(prf_file,line1))
  if ( line1.length() > 74 )
  {
   int prf1 = atoi(line1.substr(0,4).c_str()) - 1;
   
   profstr[prf1][0] = atof(line1.substr(12,9).c_str());
   profstr[prf1][1] = atof(line1.substr(21,9).c_str());
   profstr[prf1][2] = atof(line1.substr(30,9).c_str());
   profstr[prf1][3] = atof(line1.substr(39,9).c_str());
   profstr[prf1][4] = atof(line1.substr(48,9).c_str());
   profstr[prf1][5] = atof(line1.substr(57,9).c_str());
   profstr[prf1][6] = atof(line1.substr(66,9).c_str());
  }
 
 prf_file.close();
}

void read_profseq( std::string prf_name )
{
 string line1;
 
 ifstream prf_file( prf_name.c_str() );
 
 if ( !prf_file.is_open() ) { cout << "Cannot open " << prf_name << endl; exit(EXIT_FAILURE); }
 
 int l1 = 0;
 int l2 = 0;
 int l3 = 0;
 int l4 = 0;
 
 while (getline(prf_file,line1))
  if ( line1.length() > 0 )
   if ( l1++ > 0 )
   {
    istringstream iss(line1);
    
    string line2;
    
    while ( getline( iss, line2, ' ' ) )
     if ( line2.length() > 0 )
     {
      if ( l2 == 2 )
       l4 = atoi(line2.c_str());
      
      else if ( l2 > 2 )
      {
       int l5 = 0;
       
            if ( l2 == 11 ) l5 = 0;
       else if ( l2 == 3  ) l5 = 1;
       else if ( l2 == 18 ) l5 = 2;
       else if ( l2 == 17 ) l5 = 3;
       else if ( l2 == 5  ) l5 = 4;
       else if ( l2 == 12 ) l5 = 5;
       else if ( l2 == 19 ) l5 = 6;
       else if ( l2 == 6  ) l5 = 7;
       else if ( l2 == 21 ) l5 = 8;
       else if ( l2 == 7  ) l5 = 9;
       else if ( l2 == 4  ) l5 = 10;
       else if ( l2 == 16 ) l5 = 11;
       else if ( l2 == 22 ) l5 = 12;
       else if ( l2 == 15 ) l5 = 13;
       else if ( l2 == 20 ) l5 = 14;
       else if ( l2 == 14 ) l5 = 15;
       else if ( l2 == 13 ) l5 = 16;
       else if ( l2 == 8  ) l5 = 17;
       else if ( l2 == 9  ) l5 = 18;
       else if ( l2 == 10 ) l5 = 19;
       
       double l6 = atof(line2.c_str()) * double ( l4 );
       
       profseq[l3][l5] = ( l6 + aa_freq[l5] * sqrt( double ( l4 ) ) ) / ( double ( l4 ) + sqrt( double ( l4 ) ) );
      }
      
      if ( ++l2 >= MAXRES + 3 )
      {
       l3++;
       
       l2 = 0;
      }
     }
   }
 
 prf_file.close();
}
