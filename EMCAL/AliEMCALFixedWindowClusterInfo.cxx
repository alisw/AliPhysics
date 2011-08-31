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

// --- Root ---
#include <TNamed.h>
#include <TArrayI.h>

#include "AliEMCALFixedWindowClusterInfo.h"

ClassImp(AliEMCALFixedWindowClusterInfo)

//___________________________________________________________________________________________________________
AliEMCALFixedWindowClusterInfo::AliEMCALFixedWindowClusterInfo()
  : TNamed(),
    fLastPos(0),
    fIds(new TArrayI()),
    fIndexes(new TArrayI()),
    fPhi(new TArrayI()),
    fEta(new TArrayI())
{
  // Constructor
}

//___________________________________________________________________________________________________________
AliEMCALFixedWindowClusterInfo::AliEMCALFixedWindowClusterInfo(const char* name, Int_t size)
  : TNamed(name, name),
    fLastPos(0),
    fIds(new TArrayI(size)),
    fIndexes(new TArrayI(size)),
    fPhi(new TArrayI(size)),
    fEta(new TArrayI(size))
{
  // Constructor
  
  fIds->Reset(-1);
  fIndexes->Reset(-1);
  fPhi->Reset(-1);
  fEta->Reset(-1);
}

//___________________________________________________________________________________________________________
AliEMCALFixedWindowClusterInfo::AliEMCALFixedWindowClusterInfo(const TString& name, Int_t size)
  : TNamed(name, name),
    fLastPos(0),
    fIds(new TArrayI(size)),
    fIndexes(new TArrayI(size)),
    fPhi(new TArrayI(size)),
    fEta(new TArrayI(size))
{
  // Constructor
  
  fIds->Reset(-1);
  fIndexes->Reset(-1);
  fPhi->Reset(-1);
  fEta->Reset(-1);
}

//___________________________________________________________________________________________________________
AliEMCALFixedWindowClusterInfo::~AliEMCALFixedWindowClusterInfo()
{
  // Destructor

  delete fIds;
  delete fIndexes;
  delete fPhi;
  delete fEta;
}

//___________________________________________________________________________________________________________
Bool_t AliEMCALFixedWindowClusterInfo::GetInfoFromId(Int_t idclus, Int_t &index, Int_t &eta, Int_t &phi) const
{
  // Search for the unique ID and return index, eta and phi
  
  Int_t i = GetPositionFromId(idclus);
  if (i > -1){
    index = fIndexes->At(i);
    eta = fEta->At(i);
    phi = fPhi->At(i);
    return kTRUE;
  }
  return kFALSE;
}

//___________________________________________________________________________________________________________
Bool_t AliEMCALFixedWindowClusterInfo::GetInfoFromIndex(Int_t index, Int_t &idclus, Int_t &eta, Int_t &phi) const
{
  // Search for the index and return unique ID, eta and phi
  
  Int_t i = GetPositionFromIndex(index);
  if (i > -1){
    idclus = fIds->At(i);
    eta = fEta->At(i);
    phi = fPhi->At(i);
    return kTRUE;
  }
  return kFALSE;
}

//___________________________________________________________________________________________________________
Bool_t AliEMCALFixedWindowClusterInfo::SetIndexFromId(Int_t idclus, Int_t index)
{
  // Search for the unique ID and set index
  
  Int_t i = GetPositionFromId(idclus);
  if (i > -1){
    fIndexes->AddAt(index, i);
    return kTRUE;
  }
  return kFALSE;
}

//___________________________________________________________________________________________________________
void AliEMCALFixedWindowClusterInfo::Add(Int_t idclus, Int_t index, Int_t eta, Int_t phi)
{
  // Add a new element to the arrays; if necessary expand the arrays
  
  if (fLastPos >= fIds->GetSize())
    Expand(fIds->GetSize() * 2 + 1);
  
  fIds->AddAt(idclus, fLastPos);
  fIndexes->AddAt(index, fLastPos);
  fPhi->AddAt(phi, fLastPos);
  fEta->AddAt(eta, fLastPos);
  fLastPos++;
}

//___________________________________________________________________________________________________________
Bool_t AliEMCALFixedWindowClusterInfo::RemoveId(Int_t idclus)
{
  // Search for the unique ID and remove the corresponding element from the arrays
  
  Int_t i = GetPositionFromId(idclus);
  if (i > -1){
    fIds->AddAt(-1, i);
    fIndexes->AddAt(-1, i);
    fPhi->AddAt(-1, i);
    fEta->AddAt(-1, i);
    return kTRUE;
  }
  return kFALSE;
}

//___________________________________________________________________________________________________________
Bool_t AliEMCALFixedWindowClusterInfo::RemoveIndex(Int_t index)
{
  // Search for and index and remove the corresponding element from the arrays
  
  Int_t i = GetPositionFromIndex(index);
  if (i > -1){
    fIds->AddAt(-1, i);
    fIndexes->AddAt(-1, i);
    fPhi->AddAt(-1, i);
    fEta->AddAt(-1, i);
    return kTRUE;
  }
  return kFALSE;
}

//___________________________________________________________________________________________________________
void AliEMCALFixedWindowClusterInfo::Expand(Int_t size)
{
  // Expand the arrays
  
  Int_t oldsize = fIds->GetSize();
  
  fIds->Set(size);
  fIndexes->Set(size);
  fPhi->Set(size);
  fEta->Set(size);
  
  for (Int_t i = oldsize; i < size; i++){
    fIds->AddAt(-1, i);
    fIndexes->AddAt(-1, i);
    fPhi->AddAt(-1, i);
    fEta->AddAt(-1, i);
  }
}

//___________________________________________________________________________________________________________
void AliEMCALFixedWindowClusterInfo::Clear(Option_t* /*option*/)
{
  // Clear the arrays
  
  fIds->Reset(-1);
  fIndexes->Reset(-1);
  fPhi->Reset(-1);
  fEta->Reset(-1);
  fLastPos = 0;
}

//___________________________________________________________________________________________________________
Int_t AliEMCALFixedWindowClusterInfo::GetSize() const
{
  // Return the size of the arrays
  
  return fIds->GetSize();
}

//___________________________________________________________________________________________________________
Bool_t AliEMCALFixedWindowClusterInfo::ContainsId(Int_t idclus) const
{
  // Search for the unique ID and return whether or not it exists
  
  if(GetPositionFromId(idclus) > -1)
    return kTRUE;
  else
    return kFALSE;
}

//___________________________________________________________________________________________________________
Bool_t AliEMCALFixedWindowClusterInfo::ContainsIndex(Int_t index) const
{
  // Search for the index and return whether or not it exists
  
  if(GetPositionFromIndex(index) > -1)
    return kTRUE;
  else
    return kFALSE;
}

//___________________________________________________________________________________________________________
Int_t AliEMCALFixedWindowClusterInfo::GetPositionFromId(Int_t idclus) const
{
  // Return the the position corresponding to a given unique ID
  
  for (Int_t i = 0; i < fLastPos; i++){
    if (fIds->At(i) == idclus)
      return i;
  }
  return -1;
}

//___________________________________________________________________________________________________________
Int_t AliEMCALFixedWindowClusterInfo::GetPositionFromIndex(Int_t index) const
{
  // Return the the position corresponding to a given index
  
  for (Int_t i = 0; i < fLastPos; i++){
    if (fIndexes->At(i) == index)
      return i;
  }
  return -1;
}

//___________________________________________________________________________________________________________
Int_t AliEMCALFixedWindowClusterInfo::GetLastElementPosition() const
{
  // Return the position of the last non-empty element in the arrays
  
  return fLastPos;
}
