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


#include "data.h"

using namespace std;

std::string one2three( std::string resnam1 )
{
      if ( resnam1 == "A" ) { return "ALA"; }
 else if ( resnam1 == "C" ) { return "CYS"; }
 else if ( resnam1 == "D" ) { return "ASP"; }
 else if ( resnam1 == "E" ) { return "GLU"; }
 else if ( resnam1 == "F" ) { return "PHE"; }
 else if ( resnam1 == "G" ) { return "GLY"; }
 else if ( resnam1 == "H" ) { return "HIS"; }
 else if ( resnam1 == "I" ) { return "ILE"; }
 else if ( resnam1 == "K" ) { return "LYS"; }
 else if ( resnam1 == "L" ) { return "LEU"; }
 else if ( resnam1 == "M" ) { return "MET"; }
 else if ( resnam1 == "N" ) { return "ASN"; }
 else if ( resnam1 == "P" ) { return "PRO"; }
 else if ( resnam1 == "Q" ) { return "GLN"; }
 else if ( resnam1 == "R" ) { return "ARG"; }
 else if ( resnam1 == "S" ) { return "SER"; }
 else if ( resnam1 == "T" ) { return "THR"; }
 else if ( resnam1 == "V" ) { return "VAL"; }
 else if ( resnam1 == "W" ) { return "TRP"; }
 else if ( resnam1 == "Y" ) { return "TYR"; }
 
 else
 {
  cout << "Unknown residue passed to one2three: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

std::string three2one( std::string resnam1 )
{
      if ( resnam1 == "ALA" ) { return "A"; }
 else if ( resnam1 == "CYS" ) { return "C"; }
 else if ( resnam1 == "ASP" ) { return "D"; }
 else if ( resnam1 == "GLU" ) { return "E"; }
 else if ( resnam1 == "PHE" ) { return "F"; }
 else if ( resnam1 == "GLY" ) { return "G"; }
 else if ( resnam1 == "HIS" ) { return "H"; }
 else if ( resnam1 == "ILE" ) { return "I"; }
 else if ( resnam1 == "LYS" ) { return "K"; }
 else if ( resnam1 == "LEU" ) { return "L"; }
 else if ( resnam1 == "MET" ) { return "M"; }
 else if ( resnam1 == "ASN" ) { return "N"; }
 else if ( resnam1 == "PRO" ) { return "P"; }
 else if ( resnam1 == "GLN" ) { return "Q"; }
 else if ( resnam1 == "ARG" ) { return "R"; }
 else if ( resnam1 == "SER" ) { return "S"; }
 else if ( resnam1 == "THR" ) { return "T"; }
 else if ( resnam1 == "VAL" ) { return "V"; }
 else if ( resnam1 == "TRP" ) { return "W"; }
 else if ( resnam1 == "TYR" ) { return "Y"; }
 
 else
 {
  cout << "Unknown residue passed to three2one: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

int three2number( std::string resnam1 )
{
      if ( resnam1 == "ALA" ) { return  0; }
 else if ( resnam1 == "CYS" ) { return  1; }
 else if ( resnam1 == "ASP" ) { return  2; }
 else if ( resnam1 == "GLU" ) { return  3; }
 else if ( resnam1 == "PHE" ) { return  4; }
 else if ( resnam1 == "GLY" ) { return  5; }
 else if ( resnam1 == "HIS" ) { return  6; }
 else if ( resnam1 == "ILE" ) { return  7; }
 else if ( resnam1 == "LYS" ) { return  8; }
 else if ( resnam1 == "LEU" ) { return  9; }
 else if ( resnam1 == "MET" ) { return 10; }
 else if ( resnam1 == "ASN" ) { return 11; }
 else if ( resnam1 == "PRO" ) { return 12; }
 else if ( resnam1 == "GLN" ) { return 13; }
 else if ( resnam1 == "ARG" ) { return 14; }
 else if ( resnam1 == "SER" ) { return 15; }
 else if ( resnam1 == "THR" ) { return 16; }
 else if ( resnam1 == "VAL" ) { return 17; }
 else if ( resnam1 == "TRP" ) { return 18; }
 else if ( resnam1 == "TYR" ) { return 19; }
 
 else
 {
  cout << "Unknown residue passed to three2number: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

int one2number( std::string resnam1 )
{
      if ( resnam1 == "A" ) { return  0; }
 else if ( resnam1 == "C" ) { return  1; }
 else if ( resnam1 == "D" ) { return  2; }
 else if ( resnam1 == "E" ) { return  3; }
 else if ( resnam1 == "F" ) { return  4; }
 else if ( resnam1 == "G" ) { return  5; }
 else if ( resnam1 == "H" ) { return  6; }
 else if ( resnam1 == "I" ) { return  7; }
 else if ( resnam1 == "K" ) { return  8; }
 else if ( resnam1 == "L" ) { return  9; }
 else if ( resnam1 == "M" ) { return 10; }
 else if ( resnam1 == "N" ) { return 11; }
 else if ( resnam1 == "P" ) { return 12; }
 else if ( resnam1 == "Q" ) { return 13; }
 else if ( resnam1 == "R" ) { return 14; }
 else if ( resnam1 == "S" ) { return 15; }
 else if ( resnam1 == "T" ) { return 16; }
 else if ( resnam1 == "V" ) { return 17; }
 else if ( resnam1 == "W" ) { return 18; }
 else if ( resnam1 == "Y" ) { return 19; }
 
 else
 {
  cout << "Unknown residue passed to one2number: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

std::string number2three( int resnam1 )
{
      if ( resnam1 ==  0 ) { return "ALA"; }
 else if ( resnam1 ==  1 ) { return "CYS"; }
 else if ( resnam1 ==  2 ) { return "ASP"; }
 else if ( resnam1 ==  3 ) { return "GLU"; }
 else if ( resnam1 ==  4 ) { return "PHE"; }
 else if ( resnam1 ==  5 ) { return "GLY"; }
 else if ( resnam1 ==  6 ) { return "HIS"; }
 else if ( resnam1 ==  7 ) { return "ILE"; }
 else if ( resnam1 ==  8 ) { return "LYS"; }
 else if ( resnam1 ==  9 ) { return "LEU"; }
 else if ( resnam1 == 10 ) { return "MET"; }
 else if ( resnam1 == 11 ) { return "ASN"; }
 else if ( resnam1 == 12 ) { return "PRO"; }
 else if ( resnam1 == 13 ) { return "GLN"; }
 else if ( resnam1 == 14 ) { return "ARG"; }
 else if ( resnam1 == 15 ) { return "SER"; }
 else if ( resnam1 == 16 ) { return "THR"; }
 else if ( resnam1 == 17 ) { return "VAL"; }
 else if ( resnam1 == 18 ) { return "TRP"; }
 else if ( resnam1 == 19 ) { return "TYR"; }
 
 else
 {
  cout << "Unknown residue passed to number2three: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

std::string number2one( int resnam1 )
{
      if ( resnam1 ==  0 ) { return "A"; }
 else if ( resnam1 ==  1 ) { return "C"; }
 else if ( resnam1 ==  2 ) { return "D"; }
 else if ( resnam1 ==  3 ) { return "E"; }
 else if ( resnam1 ==  4 ) { return "F"; }
 else if ( resnam1 ==  5 ) { return "G"; }
 else if ( resnam1 ==  6 ) { return "H"; }
 else if ( resnam1 ==  7 ) { return "I"; }
 else if ( resnam1 ==  8 ) { return "K"; }
 else if ( resnam1 ==  9 ) { return "L"; }
 else if ( resnam1 == 10 ) { return "M"; }
 else if ( resnam1 == 11 ) { return "N"; }
 else if ( resnam1 == 12 ) { return "P"; }
 else if ( resnam1 == 13 ) { return "Q"; }
 else if ( resnam1 == 14 ) { return "R"; }
 else if ( resnam1 == 15 ) { return "S"; }
 else if ( resnam1 == 16 ) { return "T"; }
 else if ( resnam1 == 17 ) { return "V"; }
 else if ( resnam1 == 18 ) { return "W"; }
 else if ( resnam1 == 19 ) { return "Y"; }
 
 else
 {
  cout << "Unknown residue passed to number2one: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

int burial2number( std::string resnam1 )
{
      if ( resnam1 == "A" ) { return  0; }
 else if ( resnam1 == "B" ) { return  1; }
 else if ( resnam1 == "C" ) { return  2; }
 else if ( resnam1 == "D" ) { return  3; }
 else if ( resnam1 == "E" ) { return  4; }
 else if ( resnam1 == "F" ) { return  5; }
 else if ( resnam1 == "G" ) { return  6; }
 
 else
 {
  cout << "Unknown state passed to burial2number: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

int stride2number( std::string resnam1 )
{
      if ( resnam1 == "B" ) { return  0; }
 else if ( resnam1 == "C" ) { return  1; }
 else if ( resnam1 == "E" ) { return  2; }
 else if ( resnam1 == "G" ) { return  3; }
 else if ( resnam1 == "H" ) { return  4; }
 else if ( resnam1 == "I" ) { return  5; }
 else if ( resnam1 == "T" ) { return  6; }
 
 else
 {
  cout << "Unknown state passed to stride2number: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}

int number2pinak( int resnam1 )
{
      if ( resnam1 ==  0 ) { return 0; }
 else if ( resnam1 ==  1 ) { return 0; }
 else if ( resnam1 ==  2 ) { return 2; }
 else if ( resnam1 ==  3 ) { return 2; }
 else if ( resnam1 ==  4 ) { return 5; }
 else if ( resnam1 ==  5 ) { return 1; }
 else if ( resnam1 ==  6 ) { return 6; }
 else if ( resnam1 ==  7 ) { return 0; }
 else if ( resnam1 ==  8 ) { return 4; }
 else if ( resnam1 ==  9 ) { return 0; }
 else if ( resnam1 == 10 ) { return 0; }
 else if ( resnam1 == 11 ) { return 3; }
 else if ( resnam1 == 12 ) { return 5; }
 else if ( resnam1 == 13 ) { return 3; }
 else if ( resnam1 == 14 ) { return 4; }
 else if ( resnam1 == 15 ) { return 1; }
 else if ( resnam1 == 16 ) { return 1; }
 else if ( resnam1 == 17 ) { return 0; }
 else if ( resnam1 == 18 ) { return 5; }
 else if ( resnam1 == 19 ) { return 5; }
 
 else
 {
  cout << "Unknown residue passed to number2pinak: " << resnam1 << endl << endl;
  exit(EXIT_FAILURE);
 }
 
 exit(EXIT_FAILURE);
}
