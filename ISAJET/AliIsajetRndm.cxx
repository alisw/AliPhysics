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

/* $Id$ */

//-----------------------------------------------------------------------------
//   Class: AliHerwigRndm
//   Responsibilities: Interface to Root random number generator 
//                     from Fortran (re-implements FINCTION RLU_HERWIG 
//                     from HERWIG)
//   Note: Since AliGenHerwig belongs to another module (THerwig) one cannot
//         pass the ponter to the generator via static variable
//   Collaborators: AliGenHerwig class
//   Example:
//
//   root> AliGenHerwig *gener = new AliGenHerwig(-1);
//   root> AliHerwigRndm::SetHerwigRandom(new TRandom3());
//   root> AliHerwigRndm::GetHerwigRandom()->SetSeed(0);
//   root> cout<<"Seed "<< AliHerwigRndm::GetHerwigRandom()->GetSeed() <<endl;
//-----------------------------------------------------------------------------

#include <TRandom.h>

#include "AliIsajetRndm.h"

TRandom * AliIsajetRndm::fgIsajetRandom=0;

ClassImp(AliIsajetRndm)

//_______________________________________________________________________
void AliIsajetRndm::SetIsajetRandom(TRandom *ran) {
  //
  // Sets the pointer to an existing random numbers generator
  //
  if(ran) fgIsajetRandom=ran;
  else fgIsajetRandom=gRandom;
}

//_______________________________________________________________________
TRandom * AliIsajetRndm::GetIsajetRandom() {
  //
  // Retrieves the pointer to the random numbers generator
  //
  return fgIsajetRandom;
}

//_______________________________________________________________________
# define ranf    ranf_

extern "C" {
  Double_t ranf() 
  {
    // Wrapper to FUNCTION HWR from HERWIG
    // Uses static method to retrieve the pointer to the (C++) generator
      Double_t r;
      do r=AliIsajetRndm::GetIsajetRandom()->Rndm(); 
      while(0 >= r || r >= 1);
      return r;
  }
}
