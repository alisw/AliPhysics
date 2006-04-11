////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONBASE_H
#define ALIHLTMUONBASE_H

#include "TError.h"
#include "TMath.h"
#include "Utils.hpp"


// I prefer to have the option of compiling in the ROOT Assert statements.
#ifndef COMPILE_ROOT_ASSERT
#	ifndef DEBUG
		// If we are not building the DEBUG version then remove the
		// compilation of Assert.
#		undef Assert
#		define Assert(code) 
#	endif // DEBUG
#endif // COMPILE_ROOT_ASSERT


#endif // ALIHLTMUONBASE_H
