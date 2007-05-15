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

Bool_t AliTrackFitter::Fit(const TArrayI *volIds,const TArrayI *volIdsFit,
AliGeomManager::ELayerID layerRangeMin,
AliGeomManager::ELayerID layerRangeMax)
{
  //-------------------------------------------------------------------
  //
  //                      Fit the track points. 
  //
  // volIds    - the array of IDs of volumes where the residuals 
  //             will be calculated.
  // volIdsFit - the array of IDs of volumes having the points
  //              that will be fitted
  // 
  // If volIdsFit==0, the IDs of volumes having the points to fit
  // are taken in the range defined by the two last parameters
  // 
  //
  // The function fills two track-point arrays: fPVolId and fPTrack.
  // The first one contains the track points from the volumes with IDs  
  // taken from the "volIds". The second array is filled with 
  // the intersection points between the fitted track and the volumes
  // the points from the first arry belong to.
  //
  // The two arrays are used for calculation of the residuals
  // and for the construction of a chi2 function to be minimized 
  // in the alignment procedures. 
  //
  //--------------------------------------------------------------------

  Int_t npoints=fPoints->GetNPoints();
  if (npoints<fMinNPoints) return kFALSE;

  // Fast counting the points
  Int_t countFit=0;
  Int_t countPnt=0;

  Int_t fst=-1;
  Int_t lst=-1;
  if (volIdsFit != 0x0) {
     for (Int_t i=0; i<npoints; i++) {
         if (FindVolId(volIds,   fPoints->GetVolumeID()[i])) countPnt++;
         if (FindVolId(volIdsFit,fPoints->GetVolumeID()[i])) {
            countFit++;
            if (fst<0) fst=i;
            lst=i;
         }
     }
  } else {
     for (Int_t i=0; i<npoints; i++) {
         UShort_t id=fPoints->GetVolumeID()[i]; 
         if (FindVolId(volIds,id)) countPnt++;
         if (id < AliGeomManager::LayerToVolUID(layerRangeMin,0)) continue;
	 if (id > AliGeomManager::LayerToVolUID(layerRangeMax,
		  AliGeomManager::LayerSize(layerRangeMax))) continue;
         countFit++;
         if (fst<0) fst=i;
         lst=i;
     }
  }
  if (countPnt==0) return kFALSE;
  if (countFit<fMinNPoints) return kFALSE;



  //************* Fit the selected track points

  if (!Begin(fst,lst)) return kFALSE;

  AliTrackPoint p;
  if (volIdsFit != 0x0) {
     for (Int_t i=0; i<npoints; i++) {
         if (!FindVolId(volIdsFit,fPoints->GetVolumeID()[i])) continue;
         fPoints->GetPoint(p,i);
         if (!AddPoint(&p)) return kFALSE;
     }
  } else {
     for (Int_t i=0; i<npoints; i++) {
         UShort_t id=fPoints->GetVolumeID()[i]; 
         if (id < AliGeomManager::LayerToVolUID(layerRangeMin,0)) continue;
	 if (id > AliGeomManager::LayerToVolUID(layerRangeMax,
		  AliGeomManager::LayerSize(layerRangeMax))) continue;
         fPoints->GetPoint(p,i);
         if (!AddPoint(&p)) continue;
     }
  }

  if (!Update()) return kFALSE;




  //************* Calculate the intersection points

  fPVolId = new AliTrackPointArray(countPnt);
  fPTrack = new AliTrackPointArray(countPnt);

  Int_t n=0;
  AliTrackPoint p2;
  for (Int_t i=0; i<npoints; i++) {
      if (!FindVolId(volIds,fPoints->GetVolumeID()[i])) continue;
      fPoints->GetPoint(p,i);
      if (GetPCA(p,p2)) {
	fPVolId->AddPoint(n,&p);
	fPTrack->AddPoint(n,&p2);
        n++;
      } else {
	delete fPVolId;
	fPVolId=0;
	delete fPTrack;
	fPTrack=0;
	return kFALSE;
      }
  }
  
  return kTRUE;
}
