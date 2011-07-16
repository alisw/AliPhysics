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

//-----------------------------------------------------------------
//   Implementation of the base class for track residuals
//
//
//-----------------------------------------------------------------

#include "AliTrackResiduals.h"

#include "AliAlignObj.h"
#include "AliAlignObjParams.h"
#include "AliTrackPointArray.h"

ClassImp(AliTrackResiduals)

//_____________________________________________________________________________
AliTrackResiduals::AliTrackResiduals():
  fN(0),
  fLast(0),
  fAlignObj(0),
  fVolArray(0),
  fTrackArray(0),
  fChi2(0),
  fNdf(0),
  fMinNPoints(0),
  fIsOwner(kTRUE)
{
  // Default constructor
  for (Int_t ipar=0; ipar<6; ipar++){
    fBFixed[ipar] = kFALSE;
    fFixed[ipar]  = 0.;
  }  
}

//_____________________________________________________________________________
AliTrackResiduals::AliTrackResiduals(Int_t ntracks):
  fN(ntracks),
  fLast(0),
  fAlignObj(0),
  fVolArray(0),
  fTrackArray(0),
  fChi2(0),
  fNdf(0),
  fMinNPoints(0),
  fIsOwner(kTRUE)
{
  // Constructor
  if (ntracks > 0) {
    fVolArray = new AliTrackPointArray*[ntracks];
    fTrackArray = new AliTrackPointArray*[ntracks];
    for (Int_t itrack = 0; itrack < ntracks; itrack++)
      fVolArray[itrack] = fTrackArray[itrack] = 0x0;
  }

  for (Int_t ipar=0; ipar<6; ipar++){
    fBFixed[ipar] = kFALSE;
    fFixed[ipar]  = 0.;
  }  
}

//_____________________________________________________________________________
AliTrackResiduals::AliTrackResiduals(const AliTrackResiduals &res):
  TObject(res),
  fN(res.fN),
  fLast(res.fLast),
  fAlignObj(0),
  fVolArray(0),
  fTrackArray(0),
  fChi2(res.fChi2),
  fNdf(res.fNdf),
  fMinNPoints(res.fMinNPoints),
  fIsOwner(kTRUE)
{
  // Copy constructor
  // By default the created copy owns the track point arrays
  if (res.fAlignObj)
    fAlignObj = (AliAlignObj *)res.fAlignObj->Clone();

  if (fN > 0) {
    fVolArray = new AliTrackPointArray*[fN];
    fTrackArray = new AliTrackPointArray*[fN];
    for (Int_t itrack = 0; itrack < fN; itrack++)
      {
	if (res.fVolArray[itrack])
	  fVolArray[itrack] = new AliTrackPointArray(*res.fVolArray[itrack]);
	else
	  fVolArray = 0x0;
	if (res.fTrackArray[itrack])
	  fTrackArray[itrack] = new AliTrackPointArray(*res.fTrackArray[itrack]);
	else
	  fTrackArray = 0x0;
      }
  }
  for(Int_t i=0;i<6;i++) { 
    fBFixed[i]=res.fBFixed[i];
    fFixed[i]=res.fFixed[i];
  }
}

//_____________________________________________________________________________
AliTrackResiduals &AliTrackResiduals::operator =(const AliTrackResiduals& res)
{
  // assignment operator
  // Does not copy the track point arrays
  if(this==&res) return *this;
  ((TObject *)this)->operator=(res);

  fN = res.fN;
  fLast = res.fLast;
  fChi2 = res.fChi2;
  fNdf  = res.fNdf;
  fMinNPoints = res.fMinNPoints;
  fIsOwner = kFALSE;
  fAlignObj = res.fAlignObj;

  fVolArray = res.fVolArray;
  fTrackArray = res.fTrackArray;

  for(Int_t i=0;i<6;i++) { 
    fBFixed[i]=res.fBFixed[i];
    fFixed[i]=res.fFixed[i];
  }

  return *this;
}

//_____________________________________________________________________________
AliTrackResiduals::~AliTrackResiduals()
{
  // Destructor
  if (fAlignObj) delete fAlignObj;
  DeleteTrackPointArrays();
}

//_____________________________________________________________________________
void AliTrackResiduals::SetNTracks(Int_t ntracks)
{
  // Set new size for the track point arrays.
  // Delete the old arrays and allocate the
  // new ones.
  DeleteTrackPointArrays();

  fN = ntracks;
  fLast = 0;
  fChi2 = 0;
  fNdf  = 0;
  fIsOwner = kTRUE;

  if (ntracks > 0) {
    fVolArray = new AliTrackPointArray*[ntracks];
    fTrackArray = new AliTrackPointArray*[ntracks];
    for (Int_t itrack = 0; itrack < ntracks; itrack++)
      fVolArray[itrack] = fTrackArray[itrack] = 0x0;
  }
  else {
    fVolArray = fTrackArray = 0x0;
  }
}

//_____________________________________________________________________________
Bool_t AliTrackResiduals::AddTrackPointArrays(AliTrackPointArray *volarray, AliTrackPointArray *trackarray)
{
  // Adds pair of track space point and
  // track extrapolation point arrays
  if (!fVolArray || !fTrackArray) return kFALSE;

  if (!volarray || !trackarray) return kFALSE;

  if (volarray->GetNPoints() < fMinNPoints) return kFALSE;

  if (fLast >= fN) return kFALSE;

  fVolArray[fLast] = volarray;
  fTrackArray[fLast] = trackarray;
  fLast++;

  return kTRUE;
}

//_____________________________________________________________________________
void AliTrackResiduals::InitAlignObj()
{
  // Create the alignment object 
  // to be updated
  if (fAlignObj) delete fAlignObj;
  fAlignObj = new AliAlignObjParams;
}


//_____________________________________________________________________________
Bool_t AliTrackResiduals::GetTrackPointArrays(Int_t i, AliTrackPointArray* &volarray, AliTrackPointArray* &trackarray) const
{
  // Provide an access to a pair of track point arrays
  // with given index
  if (i >= fLast) {
    volarray = trackarray = 0x0;
    return kFALSE;
  }
  else {
    volarray = fVolArray[i];
    trackarray = fTrackArray[i];
    return kTRUE;
  }
}

//_____________________________________________________________________________
void AliTrackResiduals::DeleteTrackPointArrays()
{
  // Deletes the track point arrays only in case
  // the object is their owner.
  // Called by the destructor and SetNTracks methods.
  if (fIsOwner) {
    if (fVolArray) {
      for (Int_t itrack = 0; itrack < fN; itrack++) {
	if (fVolArray[itrack]) delete fVolArray[itrack];
      }
      delete [] fVolArray;
    }
    if (fTrackArray) {
      for (Int_t itrack = 0; itrack < fN; itrack++) {
	if (fTrackArray[itrack]) delete fTrackArray[itrack];
      }
      delete [] fTrackArray;
    }
  }
}

//_____________________________________________________
Int_t AliTrackResiduals::GetNFreeParam(){ 
  Int_t unfixedparam=6;
  for(Int_t j=0;j<6;j++){
    if(fBFixed[j]==kTRUE)unfixedparam--;
  }
  return unfixedparam;
}
