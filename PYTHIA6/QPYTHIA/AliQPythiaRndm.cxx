/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliPythiaRndm.cxx 31620 2009-03-20 15:31:17Z morsch $ */

//-----------------------------------------------------------------------------
//   Class: AliPythiaRndm
//   Responsibilities: Interface to Root random number generator 
//                     from Fortran (re-implements FINCTION PYR from PYTHIA)
//                     Very similar to AliHijingRndm
//   Collaborators: AliPythia and AliGenPythia classes
//   Example:
//
//   root> AliPythia::Instance();
//   root> AliPythiaRndm::SetPythiaRandom(new TRandom3()); 
//   root> AliPythiaRndm::GetPythiaRandom()->SetSeed(0);
//   root> cout<<"Seed "<< AliPythiaRndm::GetPythiaRandom()->GetSeed() <<endl;
//
//-----------------------------------------------------------------------------

#include <TMath.h>
#include <TRandom.h>

#include "AliPythiaRndm.h"

TRandom * AliPythiaRndm::fgPythiaRandom=0;

ClassImp(AliPythiaRndm)


//_______________________________________________________________________
void AliPythiaRndm::SetPythiaRandom(TRandom *ran) {
  //
  // Sets the pointer to an existing random numbers generator
  //
  if(ran) fgPythiaRandom=ran;
  else fgPythiaRandom=gRandom;
}

//_______________________________________________________________________
TRandom * AliPythiaRndm::GetPythiaRandom() {
  //
  // Retrieves the pointer to the random numbers generator
  //
  return fgPythiaRandom;
}

//_______________________________________________________________________
#define pyr        pyr_
#define pygauss    pygauss_
#define pyrset     pyrset_
#define pyrget     pyrget_

extern "C" {
    Double_t pyr(Int_t*) 
    {
	// Wrapper to FUNCTION PYR from PYTHIA
	// Uses static method to retrieve the pointer to the (C++) generator
	Double_t r;
	do r=AliPythiaRndm::GetPythiaRandom()->Rndm();
	while(0 >= r || r >= 1);
	return r;
    }
    
    Double_t pygauss(Double_t x0, Double_t sig)
    {
	Double_t s = 2.;
	Double_t v1 = 0.;
	Double_t v2 = 0.;
	
	while (s > 1.) {
	    v1 = 2. * pyr(0) - 1.;
	    v2 = 2. * pyr(0) - 1.;
	    s = v1 * v1 + v2 * v2;
	}
	return v1 * TMath::Sqrt(-2. * TMath::Log(s) / s) * sig + x0;
    }

    void pyrset(Int_t*,Int_t*) {}
    void pyrget(Int_t*,Int_t*) {}
}
