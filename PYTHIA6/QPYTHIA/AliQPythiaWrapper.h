#ifndef ALIQPYTHIAWRAPPER_H
#define ALIQPYTHIAWRAPPER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <Rtypes.h>
#include <TError.h>

class AliFastGlauber;

class AliQPythiaWrapper {
 public:
    AliQPythiaWrapper() {
	// Default constructor. The static data member is initialized 
	// in the implementation file
    }
    AliQPythiaWrapper(const AliQPythiaWrapper & /*rn*/) {
	// Copy constructor: no copy allowed for the object
	::Fatal("Copy constructor","Not allowed\n");
    }
    virtual ~AliQPythiaWrapper() {
	// Destructor
    }
    AliQPythiaWrapper & operator=(const AliQPythiaWrapper& /*rn*/) {
	// Assignment operator: no assignment allowed
	::Fatal("Assignment operator","Not allowed\n");
	return (*this);
    }
  
private:

  ClassDef(AliQPythiaWrapper, 0)  // Wrappers for C++ functionalities needed by the QPythia Fortran Code
};

#endif 

