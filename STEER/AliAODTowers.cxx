/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

//-------------------------------------------------------------------------
//     AOD class to store tower data
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include "AliAODTowers.h"

ClassImp(AliAODTowers)

  AliAODTowers::AliAODTowers() : TNamed(), fNTowers(0), fIndex(0), fAmplitude(0), fIsSorted(kTRUE), fType(kUndef)
{
  // default constructor
}

AliAODTowers::AliAODTowers(const char* name, const char* title, AODTwrs_t ttype) : TNamed(name, title), fNTowers(0), fIndex(0), fAmplitude(0), fIsSorted(kTRUE), fType(ttype)
{
  // TNamed constructor
}

void AliAODTowers::CreateContainer(Short_t nTowers)
{
  // function that creates container to store tower data

  DeleteContainer();
  
  if (nTowers <= 0) {
    fNTowers = 0;
    return;
  }

  fNTowers = nTowers;

  fIndex = new Short_t[fNTowers];
  fAmplitude = new Double32_t[fNTowers];
}

AliAODTowers::~AliAODTowers()
{
  // destructor

  DeleteContainer();
}

void AliAODTowers::DeleteContainer()
{
  // deletes allocated memory

  if (fIndex)
  {
    delete[] fIndex;
    fIndex = 0;
  }

  if (fAmplitude)
  {
    delete[] fAmplitude;
    fAmplitude = 0;
  }

  fNTowers = 0;
  fIsSorted = kFALSE;
}

Bool_t AliAODTowers::SetTower(Short_t pos, Short_t index, Double32_t amplitude)
{
  // Sets a tower at the given position

  if (pos>=0 && pos < fNTowers) {
    fIndex[pos] = index;
    fAmplitude[pos] = amplitude;
    fIsSorted = kFALSE;
    return kTRUE;
  } else {
    return kFALSE;
  }
}
