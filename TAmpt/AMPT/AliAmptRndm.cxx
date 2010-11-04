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

#include <TRandom.h>

#include "AliAmptRndm.h"

TRandom * AliAmptRndm::fgAmptRandom=0;

ClassImp(AliAmptRndm)

//_______________________________________________________________________
void AliAmptRndm::SetAmptRandom(TRandom *ran) {
  //
  // Sets the pointer to an existing random numbers generator
  //
  if(ran) fgAmptRandom=ran;
  else fgAmptRandom=gRandom;
}

//_______________________________________________________________________
TRandom * AliAmptRndm::GetAmptRandom() {
  //
  // Retrieves the pointer to the random numbers generator
  //
  return fgAmptRandom;
}

//_______________________________________________________________________
# define rluget_ampt rluget_ampt_
# define rluset_ampt rluset_ampt_
# define rlu_ampt    rlu_ampt_
# define ranart      ranart_
# define ran1        ran1_
# define rlu         rlu_

extern "C" {
  void rluget_ampt(Int_t & /*lfn*/, Int_t & /*move*/)
  {printf("Dummy version of rluget_ampt reached\n");}

  void rluset_ampt(Int_t & /*lfn*/, Int_t & /*move*/)
  {printf("Dummy version of rluset_ampt reached\n");}

  Float_t rlu_ampt(Int_t & /*idum*/) 
  {
    // Wrapper to FINCTION RLU_AMPT from AMPT
    // Uses static method to retrieve the pointer to the (C++) generator
      Double_t r;
      do r=AliAmptRndm::GetAmptRandom()->Rndm(); 
      while(0 >= r || r >= 1);
      return (Float_t)r;
  }

  Float_t ranart(Int_t &idum) 
  {
    return rlu_ampt(idum);
  }

  Float_t ran1(Int_t &idum) 
  {
    return rlu_ampt(idum);
  }

  Float_t rlu(Int_t &idum) 
  {
    return rlu_ampt(idum);
  }
}
