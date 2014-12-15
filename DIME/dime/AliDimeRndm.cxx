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
//   Responsibilities: Interface to Root random number generator 
//                     from Fortran (re-implements FINCTION RLU_HIJING 
//                     from HIJING)
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

#include <TRandom.h>

#include "AliDimeRndm.h"

TRandom * AliDimeRndm::fgDimeRandom=0;

ClassImp(AliDimeRndm)

//_______________________________________________________________________
void AliDimeRndm::SetDimeRandom(TRandom *ran) {
  //
  // Sets the pointer to an existing random numbers generator
  //
  if(ran) fgDimeRandom=ran;
  else fgDimeRandom=gRandom;
}

//_______________________________________________________________________
TRandom * AliDimeRndm::GetDimeRandom() {
  //
  // Retrieves the pointer to the random numbers generator
  //
  return fgDimeRandom;
}

//_______________________________________________________________________
//# define rluget_hijing rluget_hijing_
//# define rluset_hijing rluset_hijing_
//# define rlu_hijing    rlu_hijing_

extern "C" {
  
  //  void rluget_hijing(Int_t & /*lfn*/, Int_t & /*move*/)
  //{printf("Dummy version of rluget_hijing reached\n");}

  //void rluset_hijing(Int_t & /*lfn*/, Int_t & /*move*/)
  //{printf("Dummy version of rluset_hijing reached\n");}

  //Float_t rlu_hijing(Int_t & /*idum*/) 
  // {
  // Wrapper to FINCTION RLU_HIJING from HIJING
    // Uses static method to retrieve the pointer to the (C++) generator
  //    Double_t r;
  //    do r=AliDimeRndm::GetDimeRandom()->Rndm(); 
  //    while(0 >= r || r >= 1);
  //    return (Float_t)r;
  // }
}
