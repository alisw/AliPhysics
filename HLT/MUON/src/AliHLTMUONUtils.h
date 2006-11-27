#ifndef ALIHLTHLTMUONUTILITIES_H
#define ALIHLTHLTMUONUTILITIES_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
// We define some useful macros used in the code.
//
////////////////////////////////////////////////////////////////////////////////


/* Since c++ is missing a finally "keyword" we define one. Its usage is identical
   to a try..finally statement in Java etc.. however, since it is officialy a macro
   one must use the ( ) brackets instead of { }
 */

// if the compiler supports __finally use it otherwise make our own
#if defined(__BORLANDC__)

#	define finally(str) __finally{str}

#else

#	define finally(code) \
		catch(...) \
		{ \
			code \
			throw; \
		}; \
		code

#endif // __BORLANDC__


/* Define logical operators that are easier to read.
   and = &&
   or = ||
   not = !
 */
#if !defined(__GNUC__) && !defined(__CINT__)

#	define and &&
#	define or ||
#	define not !

#endif // GCC


#ifdef DEBUG

#	ifdef __ROOT__
#		include <TError.h>
#	endif // __ROOT__
#       ifndef Assert
#		include <cassert>
		// Define assert with a capital first letter to maintain coding style. 
#		define Assert(statement) assert(statement);
#       endif

	// Any code that should be compiled in only when the DEBUG macro is specified
	// can be enclosed with this DebugCode macro.
#	define DebugCode(code) code

#else // DEBUG

#	ifdef __ROOT__
#		include <TError.h>
#	endif // __ROOT__
#       ifndef Assert
#		define Assert(statement)
#       endif

#	define DebugCode(code)

#endif // DEBUG


#ifdef DEBUG
extern int gAliHLTMUONDebugLevel;
#endif // DEBUG


/* Here we define the DebugMsg(level, message) macro for easy embedding of debug information
   into a program. The output is only generated in programs compiled with the DEBUG macro
   defined. Here is a usage example:
   
       // statements...
       DebugMsg(2, "some debug information.");
       // statements...
       
   One can also use C++ ostream operators << like so:
   
       // statements...
       int x, y;
       DebugMsg(2, "x = " << x << " and y = " << y );
       // statements...
 */
#ifdef USE_ALILOG

#	include "AliLog.h"

	// We are defining this DebugMsg_PREECODE macro to use in DebugMsg in such a way
	// so that the code is removed when the LOG_NO_DEBUG macro is specified but 
	// compiled otherwise.
#ifdef __sun
#       include <sstream>
#else
#	include <Rstrstream.h>
#endif
#	ifndef LOG_NO_DEBUG
#		define __DebugMsg_PREECODE__(message) std::ostringstream os; os << message;
#	else // LOG_NO_DEBUG
#		define __DebugMsg_PREECODE__(message) 
#	endif // LOG_NO_DEBUG

#	define DebugMsg(level, message) \
		{ \
			__DebugMsg_PREECODE__(message) \
			AliDebugGeneral("dHLT", level, os.str().c_str()); \
		}

#else // USE_ALILOG

#	include <iostream>
#	define DebugMsg(level, message) \
		DebugCode( if (dHLT::gAliHLTMUONDebugLevel > level) std::cout << message << std::endl; )

#endif // USE_ALILOG


#endif // ALIHLTHLTMUONUTILITIES_H
