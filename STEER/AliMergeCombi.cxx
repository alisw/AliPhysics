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

AliMergeCombi::AliMergeCombi()
{
// default ctor

}

////////////////////////////////////////////////////////////////////////
AliMergeCombi::AliMergeCombi(Int_t dim, Int_t sperb)
{
// default ctor
  fDim = dim;
  fSperb = sperb;
  fLastDelta = new Int_t[fDim];
  for (Int_t i=0; i<fDim; i++) fLastDelta[i] = -1;
}

////////////////////////////////////////////////////////////////////////
AliMergeCombi::~AliMergeCombi()
{
// default dtor
  if (fLastDelta) {
    delete fLastDelta;
    fLastDelta = 0;
  }
}

////////////////////////////////////////////////////////////////////////
Bool_t AliMergeCombi::Combination(Int_t evNumber[], Int_t delta[])
{
   for (Int_t i=0; i<fDim; i++) fLastDelta[i] = delta[i] = 1;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////
