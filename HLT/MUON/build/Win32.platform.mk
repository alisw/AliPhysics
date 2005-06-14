############################################################################
#
# Author: Artur Szostak
# Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
#
############################################################################

CXX = g++
INCLUDE_PREFIX   = -I
MACRO_PREFIX     = -D
LIBPATH_PREFIX   = -L
LIBRARY_PREFIX   = -l
CXX_DEBUG        = -g3 -ggdb3
OPTIMIZATION     = -O3

BINARY_CXXFLAGS  = -fexceptions -Wall
LIBRARY_CXXFLAGS  = -fexceptions -Wall -fPIC -shared

BIN_EXT = .exe
LIB_EXT = .dll
LIB_PREFIX =

MACROS += WIN32

