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

#include "AliEMCALFixedWindowClusterInfo.h"

#include <TNamed.h>
#include <TArrayI.h>

ClassImp(AliEMCALFixedWindowClusterInfo)

AliEMCALFixedWindowClusterInfo::AliEMCALFixedWindowClusterInfo()
:TNamed(),
fSize(0),
lastElement(0),
fIds(new TArrayI()),
fIndexes(new TArrayI()),
fPhi(new TArrayI()),
fEta(new TArrayI())
{
}

AliEMCALFixedWindowClusterInfo::AliEMCALFixedWindowClusterInfo(const char* name, Int_t size)
:TNamed(name, name),
fSize(size),
lastElement(0),
fIds(new TArrayI(size)),
fIndexes(new TArrayI(size)),
fPhi(new TArrayI(size)),
fEta(new TArrayI(size))
{
  fIds->Reset(-1);
  fIndexes->Reset(-1);
  fPhi->Reset(-1);
  fEta->Reset(-1);
}

AliEMCALFixedWindowClusterInfo::AliEMCALFixedWindowClusterInfo(const TString& name, Int_t size)
:TNamed(name, name),
fSize(size),
lastElement(0),
fIds(new TArrayI(size)),
fIndexes(new TArrayI(size)),
fPhi(new TArrayI(size)),
fEta(new TArrayI(size))
{
  fIds->Reset(-1);
  fIndexes->Reset(-1);
  fPhi->Reset(-1);
  fEta->Reset(-1);
}

AliEMCALFixedWindowClusterInfo::~AliEMCALFixedWindowClusterInfo()
{
  delete fIds;
  delete fIndexes;
  delete fPhi;
  delete fEta;
}

Bool_t AliEMCALFixedWindowClusterInfo::GetInfoFromId(Int_t idclus, Int_t &index, Int_t &eta, Int_t &phi)
{
  Int_t i = GetPositionFromId(idclus);
  if (i > -1)
  {
    index = fIndexes->At(i);
    eta = fEta->At(i);
    phi = fPhi->At(i);
    return kTRUE;
  }
  return kFALSE;
}

Bool_t AliEMCALFixedWindowClusterInfo::GetInfoFromIndex(Int_t index, Int_t &idclus, Int_t &eta, Int_t &phi)
{
  Int_t i = GetPositionFromIndex(index);
  if (i > -1)
  {
    idclus = fIds->At(i);
    eta = fEta->At(i);
    phi = fPhi->At(i);
    return kTRUE;
  }
  return kFALSE;
}

Bool_t AliEMCALFixedWindowClusterInfo::SetIndexFromId(Int_t idclus, Int_t index)
{
  Int_t i = GetPositionFromId(idclus);
  if (i > -1)
  {
    fIndexes->AddAt(index, i);
    return kTRUE;
  }
  return kFALSE;
}

void AliEMCALFixedWindowClusterInfo::Add(Int_t idclus, Int_t index, Int_t eta, Int_t phi)
{
  if (lastElement >= fSize)
    Expand(fSize * 2 + 1);
  
  fIds->AddAt(idclus, lastElement);
  fIndexes->AddAt(index, lastElement);
  fPhi->AddAt(phi, lastElement);
  fEta->AddAt(eta, lastElement);
  lastElement++;
}

Bool_t AliEMCALFixedWindowClusterInfo::RemoveId(Int_t idclus)
{
  Int_t i = GetPositionFromId(idclus);
  if (i > -1)
  {
    fIds->AddAt(-1, i);
    fIndexes->AddAt(-1, i);
    fPhi->AddAt(-1, i);
    fEta->AddAt(-1, i);
    return kTRUE;
  }
  return kFALSE;
}

Bool_t AliEMCALFixedWindowClusterInfo::RemoveIndex(Int_t index)
{
  Int_t i = GetPositionFromIndex(index);
  if (i > -1)
  {
    fIds->AddAt(-1, i);
    fIndexes->AddAt(-1, i);
    fPhi->AddAt(-1, i);
    fEta->AddAt(-1, i);
    return kTRUE;
  }
  return kFALSE;
}

void AliEMCALFixedWindowClusterInfo::Expand(Int_t size)
{
  fIds->Set(size);
  fIndexes->Set(size);
  fPhi->Set(size);
  fEta->Set(size);
  
  for (Int_t i = fSize; i < size; i++)
  {
    fIds->AddAt(-1, i);
    fIndexes->AddAt(-1, i);
    fPhi->AddAt(-1, i);
    fEta->AddAt(-1, i);
  }
  
  fSize = size;
}

void AliEMCALFixedWindowClusterInfo::Clear(Option_t* option)
{
  fIds->Reset(-1);
  fIndexes->Reset(-1);
  fPhi->Reset(-1);
  fEta->Reset(-1);
  lastElement = 0;
}

Int_t AliEMCALFixedWindowClusterInfo::GetSize()
{
  return fSize;
}

Bool_t AliEMCALFixedWindowClusterInfo::ContainsId(Int_t idclus)
{
  if(GetPositionFromId(idclus) > -1)
    return kTRUE;
  else
    return kFALSE;
}

Bool_t AliEMCALFixedWindowClusterInfo::ContainsIndex(Int_t index)
{
  if(GetPositionFromIndex(index) > -1)
    return kTRUE;
  else
    return kFALSE;
}

Int_t AliEMCALFixedWindowClusterInfo::GetPositionFromId(Int_t idclus)
{
  for (Int_t i = 0; i < lastElement; i++)
  {
    if (fIds->At(i) == idclus)
      return i;
  }
  return -1;
}

Int_t AliEMCALFixedWindowClusterInfo::GetPositionFromIndex(Int_t index)
{
  for (Int_t i = 0; i < lastElement; i++)
  {
    if (fIndexes->At(i) == index)
      return i;
  }
  return -1;
}

Int_t AliEMCALFixedWindowClusterInfo::GetLastElementId()
{
  return lastElement;
}
