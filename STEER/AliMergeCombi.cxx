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

/*
$Log$
Revision 1.2  2001/10/04 15:50:23  jchudoba
Implement non default combination

Revision 1.1  2001/09/19 06:22:13  jchudoba
Class to generate combinations for merging

*/

////////////////////////////////////////////////////////////////////////
//
// AliMergeCombi.cxx
//
// returns combinations of input event numbers
//
////////////////////////////////////////////////////////////////////////

#include "AliMergeCombi.h"

ClassImp(AliMergeCombi)

//_______________________________________________________________________
AliMergeCombi::AliMergeCombi():
  fDim(1),
  fSperb(1),
  fCounter(0)
{
  //
  // default ctor
  //
}

//_______________________________________________________________________
AliMergeCombi::AliMergeCombi(Int_t dim, Int_t sperb):
  fDim(dim),
  fSperb(sperb),
  fCounter(0)
{
  //
  // Standard ctor
  //
}

//_______________________________________________________________________
AliMergeCombi::~AliMergeCombi()
{
  // default dtor
}

//_______________________________________________________________________
Bool_t AliMergeCombi::Combination(Int_t /* evNumber */ [], Int_t delta[])
{
  delta[0] = 1;
  for (Int_t i=1; i<fDim; i++) {
    if (fCounter%fSperb == 0) {
      delta[i] = 1;
    } else {
      delta[i] = 0;
    }
    fCounter++;
  }      
  return kTRUE;
}
