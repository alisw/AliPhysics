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

#include <TObjArray.h>
#include <TVector.h>

#include "AliMUONTransientDigit.h"
#include "AliLog.h"

ClassImp(AliMUONTransientDigit)

//____________________________________________________________________________
AliMUONTransientDigit::AliMUONTransientDigit()
  : AliMUONDigit()
{
// Constructor

  fTrackList=0;
}
 
//____________________________________________________________________________

  AliMUONTransientDigit::AliMUONTransientDigit(const AliMUONTransientDigit& digit)
    : AliMUONDigit(digit)
{
// Protected copy constructor

  AliFatal( "Not implemented.");
}


AliMUONTransientDigit::AliMUONTransientDigit(Int_t ich, Int_t *digits): 
    AliMUONDigit(digits)
{
    //
    // Creates a MUON digit list object
    //
    fChamber     = ich;
    fTrackList   = new TObjArray(5);   
    // 5 is arbitrary number, just to decrease default 16
}

////////////////////////////////////////////////////////////////////////
AliMUONTransientDigit::~AliMUONTransientDigit() 
{
  fTrackList->Delete();
  delete fTrackList;
}

////////////////////////////////////////////////////////////////////////
AliMUONTransientDigit& 
AliMUONTransientDigit::operator =(const AliMUONTransientDigit& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}

////////////////////////////////////////////////////////////////////////
void AliMUONTransientDigit::AddToTrackList(Int_t track, Int_t charge)
{
  TVector *pTrInfo = new TVector(3);
  TVector &trInfo = *pTrInfo;
  trInfo(0) = track;
  trInfo(1) = charge;
  fTrackList->Add(pTrInfo);
}

////////////////////////////////////////////////////////////////////////
void AliMUONTransientDigit::UpdateTrackList(Int_t track, Int_t charge)
{
  Int_t lastEntry = fTrackList->GetLast();
  TVector *pVect = static_cast<TVector*>(fTrackList->At(lastEntry));
  if ( static_cast<Int_t>((*pVect)(0)) == track) {
    (*pVect)(1) += charge;  // update charge
  } else {
    AddToTrackList(track,charge);
  }
}

////////////////////////////////////////////////////////////////////////
Int_t AliMUONTransientDigit::GetTrack(Int_t i) const
{
  if (i > fTrackList->GetEntriesFast()) return 0;
  TVector *pVect = static_cast<TVector*>(fTrackList->At(i));
  return static_cast<Int_t>((*pVect)(0));
}


////////////////////////////////////////////////////////////////////////
Int_t AliMUONTransientDigit::GetCharge(Int_t i) const
{
  if (i > fTrackList->GetEntriesFast()) return 0;
  TVector *pVect = static_cast<TVector*>(fTrackList->At(i));
  return static_cast<Int_t>((*pVect)(1));
}

