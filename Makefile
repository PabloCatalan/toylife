# /******************************************************************************
#  *   Copyright (C) 2018  Pablo Catalan                                        *
#  *                                                                            *
#  *  This program is free software: you can redistribute it and/or modify      *
#  *  it under the terms of the GNU General Public License as published by      *
#  *  the Free Software Foundation, either version 3 of the License, or         *
#  *  (at your option) any later version.                                       *
#  *                                                                            *
#  *  This program is distributed in the hope that it will be useful,           *
#  *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
#  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
#  *  GNU General Public License for more details.                              *
#  *                                                                            *
#  *  You should have received a copy of the GNU General Public License         *
#  *  along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
#  ******************************************************************************/

############## Operating-System and Processor Architecture Detection #########
# (many thanks to Juan Antonio García-Martín for this part of the code)
# This section attempts to auto-detect the OS. This is used only for convenience
# and allows the same Makefile to be used on multiple operating systems
# without modification.
# You can circumvent auto-detection by setting these environment variables:
#   RNA_MAKE_OS   -- the operating system name (e.g. Linux, Windows, Mac)
#   RNA_MAKE_ARCH -- the target architecture (e.g. x86 or x86_64)
# (These can be set as environment variables or on the MAKE command-line. 
#    e.g. `make all RNA_MAKE_OS=Linux`)
    ifneq (${RNA_MAKE_OS},) 
      #if RNA_MAKE_OS is NOT blank, use it as the OS
      OPSYSTEM=$(RNA_MAKE_OS)
    else ifeq (${OPSYSTEM},)
      # IF both RNA_MAKE_OS and OPSYSTEM are blank, use the `uname` command 
      #   (if available) to determine the OS.
      #   Replace 'UNKNOWN' with default OS if desired. 
      OPSYSTEM_RAW:=$(shell uname -s 2>/dev/null || echo UNKNOWN) 
      # Perform some replacements to normalize the output of uname on various systems.
      # OS_REPLACEMENTS= CYGWIN%=Windows MSYS%=Windows Darwin=Mac GNU%=Linux
      OPSYSTEM := $(OPSYSTEM_RAW:CYGWIN%=Windows)
      OPSYSTEM := $(OPSYSTEM:MSYS%=Windows)
      OPSYSTEM := $(OPSYSTEM:Darwin=Mac)
      OPSYSTEM := $(OPSYSTEM:GNU%=Linux)
      $(if $(DEBUG), $(info Make: Operating System: $(OPSYSTEM)))
      export OPSYSTEM #make it available for recursive calls to make, so auto-detection is performed only once.
    endif


CXX     	= g++
CXXFLAGS        = -fPIC -std=c++0x -O4 -DNDEBUG 
    ifeq (${OPSYSTEM},Linux)
      ############ LINUX ##########################################################
      LDFLAGS         = -lstdc++ -lgcc -lm
    else ifeq (${OPSYSTEM},Mac)
      ############ MAC ############################################################
      LDFLAGS         = -lstdc++ -lgcc -lm
    else ifeq (,${OPSYSTEM})
      ############ NO OS ########################################################
	  $(error No Operating system defined!!!)
    else ifeq (${OPSYSTEM},Windows)
      ############ WINDOWS ########################################################
      LDFLAGS         = -lstdc++ -lgcc -lm
   else
	  ############ UNKNOWN OS ###################################################
      $(error Unknown Operating system defined: $(OPSYSTEM))
    endif

all: example

example: example.o helper_functions.o toy_plugin.o
	$(CXX) $(CXXFLAGS) example.o helper_functions.o toy_plugin.o $(LDFLAGS) -o example

example.o: example.cpp helper_functions.o toy_plugin.o
	$(CXX) $(CXXFLAGS) -c example.cpp -o example.o

helper_functions.o: helper_functions.cpp
	$(CXX) $(CXXFLAGS) -c helper_functions.cpp -o helper_functions.o

toy_plugin.o: toy_plugin.cpp
	$(CXX) $(CXXFLAGS) -c toy_plugin.cpp -o toy_plugin.o

entrapmentplugin.o: entrapmentplugin.cpp
	$(CXX) $(CXXFLAGS) -c entrapmentplugin.cpp -o entrapmentplugin.o

clean:
	rm -rf example *.o 

