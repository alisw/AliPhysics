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
//   Implementation of the base class for fast track fitters
//
//
//-----------------------------------------------------------------

#include <TMatrixDSym.h>
#include <TArrayI.h>

#include "AliTrackFitter.h"
#include "AliTrackPointArray.h"
#include "AliLog.h"

ClassImp(AliTrackFitter)

//_____________________________________________________________________________
AliTrackFitter::AliTrackFitter() :
  TObject(),
  fCov(0),
  fPoints(0),
  fPVolId(0),
  fPTrack(0),
  fChi2(0),
  fNdf(0),
  fMinNPoints(0),
  fIsOwner(kFALSE)
{
  // default constructor
  //
  for (Int_t i=0;i<6;i++) fParams[i] = 0;
}

//_____________________________________________________________________________
AliTrackFitter::AliTrackFitter(AliTrackPointArray *array, Bool_t owner) :
  TObject(),
  fCov(new TMatrixDSym(6)),
  fPoints(0),
  fPVolId(0),
  fPTrack(0),
  fChi2(0),
  fNdf(0),
  fMinNPoints(0),
  fIsOwner(kFALSE)
  
{
  // constructor from space points array
  //
  for (Int_t i=0;i<6;i++) fParams[i] = 0;
  SetTrackPointArray(array,owner);
}

//_____________________________________________________________________________
AliTrackFitter::AliTrackFitter(const AliTrackFitter &fitter):
  TObject(fitter),
  fCov(new TMatrixDSym(*fitter.fCov)),
  fPoints(0),
  fPVolId(0),
  fPTrack(0),
  fChi2(fitter.fChi2),
  fNdf(fitter.fNdf),
  fMinNPoints(fitter.fMinNPoints),
  fIsOwner(kFALSE)
{
  // Copy constructor
  //
  SetTrackPointArray(fitter.fPoints,fitter.fIsOwner);
  for (Int_t i=0;i<6;i++) fParams[i] = fitter.fParams[i];
}

//_____________________________________________________________________________
AliTrackFitter &AliTrackFitter::operator =(const AliTrackFitter& fitter)
{
  // assignment operator
  //
  if(this==&fitter) return *this;

  SetTrackPointArray(fitter.fPoints);
  for (Int_t i=0;i<6;i++) fParams[i] = fitter.fParams[i];
  fCov = new TMatrixDSym(*fitter.fCov);
  fChi2 = fitter.fChi2;
  fNdf = fitter.fNdf;
  fMinNPoints = fitter.fMinNPoints;
  fIsOwner = kFALSE;
  
  return *this;
}

//_____________________________________________________________________________
AliTrackFitter::~AliTrackFitter()
{
  if (fIsOwner)
    delete fPoints;
  delete fCov;
}

//_____________________________________________________________________________
void AliTrackFitter::Reset()
{
  for (Int_t i=0;i<6;i++) fParams[i] = 0;
  delete fCov;
  fCov = new TMatrixDSym(6);
  fPVolId = fPTrack = 0;
  fChi2 = 0;
  fNdf = 0;
}

void AliTrackFitter::SetTrackPointArray(AliTrackPointArray *array, Bool_t owner)
{
  // Load space points from array
  // By default we don't copy them but
  // just put the pointers to them
  if (!array) {
    AliWarning("Invalid pointer to the space-points array !");
    if (fIsOwner) delete fPoints;
    fPoints = NULL;
    return;
  }

  Reset();

  if (fIsOwner) delete fPoints;

  if (owner) {
    fPoints = new AliTrackPointArray(*array);
    fIsOwner = kTRUE;
  }
  else {
    fPoints = array;
    fIsOwner = kFALSE;
  }
}

Bool_t AliTrackFitter::FindVolId(const TArrayI *array, UShort_t volid) const
{
  // The method is used to check whenever
  // the volume id (volid) is contained in
  // a array of integers
  Int_t nVolIds = array->GetSize();
  if (nVolIds == 0) return kFALSE;

  Bool_t found = kFALSE;
  for (Int_t iVolId = 0; iVolId < nVolIds; iVolId++) {
    if ((*array)[iVolId] == volid) {
      found = kTRUE;
      break;
    }
  }

  return found;
}
