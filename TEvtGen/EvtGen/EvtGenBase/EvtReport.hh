//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtGen/EvtReport.hh
//
// Description:
//
// Modification history:
//
// Author:      Simon Patton
// Created:     Mon Jun  3 12:45:03 EDT 1996
//
//------------------------------------------------------------------------

#if !defined(TOOLBOX_FUNCTIONS_HH)
#define TOOLBOX_FUNCTIONS_HH

#if !defined(FILENAME_ONLY) /* relative path includes */

// system include files
#include <iostream>

// user include files

#else /* filename-only includes */
#include <iostream>
#include <types.h>
#endif /* filename-only includes */
// system include files

// user include files

// forward declarations

//
// constants, enums and typedefs
//
namespace Severity {
  enum Enum {
    Emergency,           // fatal
    Alert,               // requires immediate action
    Critical,            // serious
    Error,
    Warning,
    Notice,              // "normal but significant"
    Info,                // informational
    Debug                // debug
  };
};

// function declaration
std::ostream& report( Severity::Enum severity ,
                 const char* facility = 0 ) ;

// inline function definitions

#endif /* TOOLBOX_FUNCTIONS_HH */

