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

#include "AliTrackFitter.h"
#include "AliTrackPointArray.h"

ClassImp(AliTrackFitter)

//_____________________________________________________________________________
AliTrackFitter::AliTrackFitter()
{
  // default constructor
  //
  for (Int_t i=0;i<6;i++) fParams[i] = 0;
  fCov = 0;
  fPoints = 0;
  fIsOwner = kFALSE;
}

//_____________________________________________________________________________
AliTrackFitter::AliTrackFitter(AliTrackPointArray *array, Bool_t owner)
{
  // constructor from space points array
  //
  for (Int_t i=0;i<6;i++) fParams[i] = 0;
  fCov = new TMatrixDSym(6);
  fIsOwner = kFALSE;
  SetTrackPointArray(array,owner);
}

//_____________________________________________________________________________
AliTrackFitter::AliTrackFitter(const AliTrackFitter &fitter):
  TObject(fitter)
{
  // Copy constructor
  //
  for (Int_t i=0;i<6;i++) fParams[i] = fitter.fParams[i];
  fCov = new TMatrixDSym(*fitter.fCov);
  fIsOwner = kFALSE;
  SetTrackPointArray(fitter.fPoints,fitter.fIsOwner);
}

//_____________________________________________________________________________
AliTrackFitter &AliTrackFitter::operator =(const AliTrackFitter& fitter)
{
  // assignment operator
  //
  if(this==&fitter) return *this;

  for (Int_t i=0;i<6;i++) fParams[i] = fitter.fParams[i];
  fCov = new TMatrixDSym(*fitter.fCov);
  fIsOwner = kFALSE;
  SetTrackPointArray(fitter.fPoints);
  
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
}

void AliTrackFitter::SetTrackPointArray(AliTrackPointArray *array, Bool_t owner)
{
  // Load space points from array
  // By default we don't copy them but
  // just put the pointers to them
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
