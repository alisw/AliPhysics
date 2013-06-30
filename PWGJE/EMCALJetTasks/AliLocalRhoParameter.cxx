// $Id$
//
// Container and interface class to store and return local energy density
// parameter rho. The global rho value can be obtained via the members of AliRhoParameter,
// in addition this class provides rho as a function (TF1) of one parameter (denoted as phi).
// A local value of rho is evaluated in area phi-r, phi+r where r is a user defined
// parameter (e.g. a jet radius)
// Functions are implemented inline for optimization.
//
// Author: Redmer Alexander Bertens, Utrecht University, Utrecht, Netherlands
//         (rbertens@cern.ch, rbertens@nikhef.nl, r.a.bertens@uu.nl)

#include "AliRhoParameter.h"
#include "AliLocalRhoParameter.h"

ClassImp(AliLocalRhoParameter)

//________________________________________________________________________
AliLocalRhoParameter::AliLocalRhoParameter() : 
  AliRhoParameter(),
  fLocalRho(0) 
{ 
  // Constructor for root IO. 
}

//________________________________________________________________________
AliLocalRhoParameter::AliLocalRhoParameter(const char* name, Double_t val) : 
  AliRhoParameter(name, val), 
  fLocalRho(0x0)
{ 
  // Constructor
}
