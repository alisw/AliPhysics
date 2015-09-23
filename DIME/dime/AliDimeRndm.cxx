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
//   Class: AliDimeRndm
//   Responsibilities: Interface to ROOT random number generator 
//                     from Fortran (re-implements function rn from DIME)
//
//   Note: Since AliGenDime belongs to another module (TDime) one cannot
//         pass the ponter to the generator via static variable
//   Collaborators: AliGenDime class
//   Example:
//
//   root> AliGenDime *gener = new AliGenDime(-1);
//   root> AliDimeRndm::SetDimeRandom(new TRandom3());
//   root> AliDimeRndm::GetDimeRandom()->SetSeed(0);
//   root> cout<<"Seed "<< AliDimeRndm::GetDimeRandom()->GetSeed() <<endl;
//-----------------------------------------------------------------------------
// NOTE! TRandom3 algorithm is being used (Mersenne-Twister),
// not TRandom which is a bad quality Linear Congruential generator in ROOT
//
// Author: Mikael.Mieskolainen@cern.ch


#include <TRandom3.h>

#include "AliDimeRndm.h"

TRandom3* AliDimeRndm::fgDimeRandom=0;

ClassImp(AliDimeRndm)

//_______________________________________________________________________
void AliDimeRndm::SetDimeRandom(TRandom3* ran) {
  //
  // Sets the pointer to an existing random numbers generator
  //
  if (ran) {
    fgDimeRandom = ran;
    printf("AliDimeRndm::SetDimeRandom(): Seed: %d\n", fgDimeRandom->GetSeed());
  } else {
    printf("AliDimeRndm::SetDimeRandom(): Problem with TRandom3* object\n");
    // Do nothing
  }
}

//_______________________________________________________________________
TRandom3* AliDimeRndm::GetDimeRandom() {
  //
  // Retrieves the pointer to the random numbers generator
  //
  return fgDimeRandom;
}

//_______________________________________________________________________
# define rn rn_

extern "C" {

  Double_t rn(Int_t& /*idum*/)
   {
  // Wrapper to function rn() from DIME
  // Uses static method to retrieve the pointer to the (C++) generator
      Double_t r = 0;
      do {
        r = AliDimeRndm::GetDimeRandom()->Rndm();
        //printf("Random number %0.5f \n", r);
      } while (0 >= r || r >= 1);
      
      return r;
  }

}
