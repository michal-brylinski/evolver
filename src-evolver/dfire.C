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


#include "dfire.h"

using namespace std;

void read_potential_sc( std::string dfire1_name )
{
 for ( int i = 0; i < MAXDF1; i++ )
  dfire1[i] = 0;
 
 for ( int i = 0; i < MAXDF4; i++ )
  for ( int j = 0; j < MAXRES; j++ )
   for ( int k = 0; k < MAXRES; k++ )
    dfire2[i][j][k] = 0.0;
 
 string line1;
 
 ifstream dfire1_file( dfire1_name.c_str() );
 
 if ( !dfire1_file.is_open() ) { cout << "Cannot open " << dfire1_name << endl; exit(EXIT_FAILURE); }
 
 int l1 = 0;
 
 while (getline(dfire1_file,line1))
  if ( line1.length() > 5 )
  {
   if ( l1 > 1 && l1 < MAXDF1 + 2 )
   {
    istringstream iss(line1);
    
    string line2;
    
    int l2 = 0;
    
    while ( getline( iss, line2, ' ' ) )
     if ( l2++ < 1 )
      dfire1[l1-2] = atoi(line2.c_str()) - 1;
   }
   else if ( l1 > MAXDF1 + 1 )
   {
    istringstream iss(line1);
    
    string line2;
    
    int l2 = 0;
    
    int b1 = 0;
    int b2 = 0;
    int b3 = 0;
    
    double b4 = 0.0;
    
    while ( getline( iss, line2, ' ' ) )
    {
     if ( l2 == 3 )
     {
      b1 = atoi(line2.c_str()) - 1;
     }
     else if ( l2 == 0 )
     {
      b2 = three2number(line2.c_str());
     }
     else if ( l2 == 1 )
     {
      b3 = three2number(line2.c_str());
     }
     else if ( l2 == 2 )
     {
      b4 = atof(line2.c_str());
     }
     
     l2++;
    }
    
    dfire2[b1][b2][b3] = b4;
    dfire2[b1][b3][b2] = b4;
   }
   
   l1++;
  }
 
 dfire1_file.close();
}

void read_potential_ca( std::string dfire2_name )
{
 for ( int i = 0; i < MAXDF1; i++ )
  dfire3[i] = 0;
 
 for ( int i = 0; i < MAXDF3; i++ )
  for ( int j = 0; j < MAXRES; j++ )
   for ( int k = 0; k < MAXRES; k++ )
    dfire4[i][j][k] = 0.0;
 
 string line1;
 
 ifstream dfire2_file( dfire2_name.c_str() );
 
 if ( !dfire2_file.is_open() ) { cout << "Cannot open " << dfire2_name << endl; exit(EXIT_FAILURE); }
 
 int l1 = 0;
 
 while (getline(dfire2_file,line1))
  if ( line1.length() > 5 )
  {
   if ( l1 > 1 && l1 < MAXDF2 + 2 )
   {
    istringstream iss(line1);
    
    string line2;
    
    int l2 = 0;
    
    while ( getline( iss, line2, ' ' ) )
     if ( l2++ < 1 )
      dfire3[l1-2] = atoi(line2.c_str()) - 1;
   }
   else if ( l1 > MAXDF2 + 1 )
   {
    istringstream iss(line1);
    
    string line2;
    
    int l2 = 0;
    
    int b1 = 0;
    int b2 = 0;
    int b3 = 0;
    
    double b4 = 0.0;
    
    while ( getline( iss, line2, ' ' ) )
    {
     if ( l2 == 2 )
     {
      b1 = atoi(line2.c_str()) - 1;
     }
     else if ( l2 == 0 )
     {
      b2 = three2number(line2.c_str());
     }
     else if ( l2 == 1 )
     {
      b3 = three2number(line2.c_str());
     }
     else if ( l2 == 3 )
     {
      b4 = atof(line2.c_str());
     }
     
     l2++;
    }
    
    dfire4[b1][b2][b3] = b4;
    dfire4[b1][b3][b2] = b4;
   }
   
   l1++;
  }
 
 dfire2_file.close();
}
