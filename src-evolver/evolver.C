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


#include "evolver.h"

using namespace std;

int main(int argc, char *argv[])
{
 srand(time(0));
 
 time_t t_start, t_end;
 
 time(&t_start);
 
 dif1 = 0.0;
 
 cout << "------------------------------------------------------------" << endl
      << "                          evolver" << endl
      << "                        version 1.0" << endl
      << "                 protein sequence generator" << endl
      << "            report bugs to michal@brylinski.org" << endl
      << "------------------------------------------------------------" << endl << endl;
 
 if ( argc < 2 )
 {
  cout << " Usage and options:" << endl << endl
       << "  evolver -p <target protein structure, PDB>" << endl
       << "          -s <STRIDE assignment>" << endl
       << "          -r <structure profile>" << endl
       << "          -o <output filename>" << endl
       << "          -f <output format, default fasta>" << endl
       << "          -q <initial sequence, default 0>" << endl
       << "          -v <verbosity level, default 0>" << endl << endl
       << " Options for the output format: " << endl << endl
       << "  fasta" << endl
       << "  pdb" << endl << endl
       << " Options for the initial sequence:" << endl << endl
       << "  0 - from input PDB" << endl
       << "  1 - from input PDB but shuffled" << endl
       << "  2 - random sequence" << endl << endl
       << " Options for the verbosity level: " << endl << endl
       << "  0 - no progress from SA" << endl
       << "  1 - reduced progress" << endl
       << "  2 - full progress" << endl << endl;
  
  exit(EXIT_SUCCESS);
 }
 
 string protein_name;
 bool protein_opt = false;
 
 string profseq_name;
 profseq_opt = false;
 
 string profstr_name;
 bool profstr_opt = false;
 
 string stride_name;
 bool stride_opt = false;
 
 string output_name;
 bool output_opt = false;
 
 int seq_ini = 0;
 
 string out_format = "fasta";
 
 verb_opt = 0;
 
 for ( int i = 0; i < argc; i++ )
 {
  if ( !strcmp(argv[i],"-p")  && i < argc ) { protein_name = string(argv[i+1]); protein_opt = true; }
  if ( !strcmp(argv[i],"-s")  && i < argc ) { stride_name  = string(argv[i+1]); stride_opt  = true; }
  if ( !strcmp(argv[i],"-r")  && i < argc ) { profstr_name = string(argv[i+1]); profstr_opt = true; }
  if ( !strcmp(argv[i],"-a")  && i < argc ) { profseq_name = string(argv[i+1]); profseq_opt = true; }
  if ( !strcmp(argv[i],"-o")  && i < argc ) { output_name  = string(argv[i+1]); output_opt  = true; }
  if ( !strcmp(argv[i],"-q")  && i < argc ) { seq_ini      = atoi(argv[i+1]);                       }
  if ( !strcmp(argv[i],"-f")  && i < argc ) { out_format   = string(argv[i+1]);                     }
  if ( !strcmp(argv[i],"-v")  && i < argc ) { verb_opt     = atoi(argv[i+1]);                       }
 }
 
 if ( !protein_opt )
 {
  cout << "Provide target protein structure" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !profstr_opt )
 {
  cout << "Provide structure profile" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !stride_opt )
 {
  cout << "Provide STRIDE assignment" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !output_opt )
 {
  cout << "Provide output filename" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( seq_ini != 0 && seq_ini != 1 && seq_ini != 2 )
 {
  cout << "Available options for the initial sequence are: 0, 1 and 2" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( out_format.compare("fasta") != 0 && out_format.compare("pdb") != 0 )
 {
  cout << "Available options for the output format are: fasta and pdb" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( verb_opt != 0 && verb_opt != 1 && verb_opt != 2 )
 {
  cout << "Available verbosity levels are: 0, 1 and 2" << endl;
  exit(EXIT_FAILURE);
 }
 
 char * path1;
 
 path1 = getenv("EVOLDAT"); if ( path1==NULL ) { cout << "EVOLDAT is not set" << endl; exit(EXIT_FAILURE); }
 
 string data_path;
 data_path = getenv("EVOLDAT");
 
 read_params( data_path+"/params.dat" );
 
 read_potential_sc( data_path+"/fort.21_scm" );
 
 read_potential_ca( data_path+"/fort.21_1.61_CA" );
 
 read_pdb( protein_name );
 
 read_stride( stride_name );
 
 read_profstr( profstr_name );
 
 if ( profseq_opt )
  read_profseq( profseq_name );
 
 int seq1[MAXPRO];
 
 for ( int i = 0; i < res; i++ )
  seq1[i] = sequence[i];
 
 if ( seq_ini == 1 )
  seq_shuffle(seq1);
 
 else if ( seq_ini == 2 )
  seq_random(seq1);
 
 double erg1[MAXCMP];
 
 erg1[0] = get_energy_burial(seq1);
 erg1[1] = get_energy_stride(seq1);
 erg1[2] = get_energy_profile(seq1);
 erg1[3] = get_energy_dfire_ca(seq1);
 erg1[4] = get_energy_dfire_sc(seq1);
 erg1[5] = get_energy_pot(seq1);
 
 for ( int i = 0; i < MAXCMP; i++ )
  Enormalize( erg1[i], weights[i][0], weights[i][1] );
 
 erg1[5] = 1.0 - erg1[5];
 
 double erg2 = 0.0;
 
 for ( int i = 0; i < MAXCMP; i++ )
 {
  erg1[i] *= weights[i][2];
  
  erg2 += erg1[i];
 }
 
 cout << "Initial score: " << fixed << setw(9) << setprecision(4) << erg2  << endl << endl;
 
 const gsl_rng * r = gsl_rng_alloc (gsl_rng_env_setup());
 
 gsl_ieee_env_setup ();
 
 if ( verb_opt == 0 )
 {
  cout << "Running simulated annealing ... " << flush;
  
  gsl_siman_solve(r, seq1, Etsp, Stsp, Mtsp, NULL, NULL, NULL, NULL, MAXPRO*sizeof(int), params);
  
  cout << "done" << endl;
 }
 
 else
 {
  cout << "Running simulated annealing: " << endl << endl;
  
  gsl_siman_solve(r, seq1, Etsp, Stsp, Mtsp, Ptsp, NULL, NULL, NULL, MAXPRO*sizeof(int), params);
 }
 
 erg1[0] = get_energy_burial(seq1);
 erg1[1] = get_energy_stride(seq1);
 erg1[2] = get_energy_profile(seq1);
 erg1[3] = get_energy_dfire_ca(seq1);
 erg1[4] = get_energy_dfire_sc(seq1);
 erg1[5] = get_energy_pot(seq1);
 
 for ( int i = 0; i < MAXCMP; i++ )
  Enormalize( erg1[i], weights[i][0], weights[i][1] );
 
 erg1[5] = 1.0 - erg1[5];
 
 erg2 = 0.0;
 
 for ( int i = 0; i < MAXCMP; i++ )
 {
  erg1[i] *= weights[i][2];
  
  erg2 += erg1[i];
 }
 
 if ( out_format.compare("fasta") == 0 )
  out_fasta(output_name, seq1, erg2);
 
 else if ( out_format.compare("pdb") == 0 )
  out_pdb(output_name, seq1, erg2);
 
 pdbcontent.clear();
 
 cout << endl << "Final score:   " << fixed << setw(9) << setprecision(4) << erg2  << endl << endl;
 
 if ( seq_ini == 0 )
  cout << "SID to the initial sequence: " << fixed << setprecision(1) << ( get_sequence_score(seq1) * 100.0 ) << "%" << endl << endl;
 
 cout << "Results written to " << output_name  << endl << endl;
 
 time(&t_end);
 
 double dif2 = difftime(t_end, t_start);
 
 cout << "Simulation time in sec:         " << fixed << setprecision(0) << setw(10) << dif2 << endl;
 cout << "Function evaluations (total):   " << fixed << setprecision(0) << setw(10) << dif1 << endl;
 cout << "Function evaluations (per sec): " << fixed << setprecision(0) << setw(10) << dif1 / dif2 << endl;
 
 return 0;
}
