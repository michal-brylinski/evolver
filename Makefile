#===============================================================================
#           ___ ___         __                   
#   .-----.|   |   |.-----.|  |.--.--.-----.----.
#   |  -__||   |   ||  _  ||  ||  |  |  -__|   _|
#   |_____| \_____/ |_____||__| \___/|_____|__|  
#                                                  
#   eVolver - protein sequence generator
#
#   Computational Systems Biology Group
#   Department of Biological Sciences
#   Center for Computation & Technology
#   Louisiana State University
#   407 Choppin Hall, Baton Rouge, LA 70803, USA
#
#   http://www.brylinski.org
#
#   Report bugs to michal@brylinski.org
#
#   Copyright 2013 Michal Brylinski
#
#   This file is part of eVolver.
#
#   eVolver is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   eVolver is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with eVolver. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================


EXE = evolver eprofile

OBJ = data.o dfire.o energy.o evolver.o output.o params.o parser.o shuffle.o siman.o

CXX = g++

SH = sh

PERL = perl

CPPFLAGS = -O2 -Wall -I../lib -I.

LDFLAGS = -lm -lgsl -lgslcblas -L../lib -L.

default: $(EXE)

all: $(EXE)

evolver: $(OBJ)
	$(CXX) -o $@ $(OBJ) $(LDFLAGS)
	@mv evolver ../bin/

data.o: data.C
	$(CXX) $(CPPFLAGS) -c -o data.o data.C

dfire.o: dfire.C
	$(CXX) $(CPPFLAGS) -c -o dfire.o dfire.C

energy.o: energy.C
	$(CXX) $(CPPFLAGS) -c -o energy.o energy.C

evolver.o: evolver.C
	$(CXX) $(CPPFLAGS) -c -o evolver.o evolver.C

output.o: output.C
	$(CXX) $(CPPFLAGS) -c -o output.o output.C

params.o: params.C
	$(CXX) $(CPPFLAGS) -c -o params.o params.C

parser.o: parser.C
	$(CXX) $(CPPFLAGS) -c -o parser.o parser.C

shuffle.o: shuffle.C
	$(CXX) $(CPPFLAGS) -c -o shuffle.o shuffle.C

siman.o: siman.C
	$(CXX) $(CPPFLAGS) -c -o siman.o siman.C

eprofile:
	$(SH) eprofile.shar
	@chmod +x eprofile
	@mv eprofile ../bin/

clean:
	@(rm -f ${EXE} ../bin/evolver ../bin/eprofile ${OBJ})

