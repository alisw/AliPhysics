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

// ---------------------
// Class AliMUONTrackK
// ---------------------
// Reconstructed track in the muons system based on the extended 
// Kalman filter approach
// Author: Alexander Zinchenko, JINR Dubna

#include <Riostream.h>
#include <TClonesArray.h>
#include <TArrayD.h>
#include <TMatrixD.h>

#include "AliMUONTrackK.h"
#include "AliMUON.h"
#include "AliMUONConstants.h"

#include "AliMUONTrackReconstructorK.h"
#include "AliMagF.h"
#include "AliMUONSegment.h"
#include "AliMUONHitForRec.h"
#include "AliMUONRawCluster.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONEventRecoCombi.h"
#include "AliMUONDetElement.h"
#include "AliRun.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMUONTrackK) // Class implementation in ROOT context
/// \endcond

const Int_t AliMUONTrackK::fgkSize = 5;
const Int_t AliMUONTrackK::fgkNSigma = 12; 
const Double_t AliMUONTrackK::fgkChi2max = 100; 
const Int_t AliMUONTrackK::fgkTriesMax = 10000; 
const Double_t AliMUONTrackK::fgkEpsilon = 0.002; 

void mnvertLocal(Double_t* a, Int_t l, Int_t m, Int_t n, Int_t& ifail); // from AliMUONTrack

Int_t AliMUONTrackK::fgDebug = -1; //-1;
Int_t AliMUONTrackK::fgNOfPoints = 0; 
AliMUONTrackReconstructorK* AliMUONTrackK::fgTrackReconstructor = NULL; 
TClonesArray* AliMUONTrackK::fgHitForRec = NULL; 
AliMUONEventRecoCombi *AliMUONTrackK::fgCombi = NULL;

  //__________________________________________________________________________
AliMUONTrackK::AliMUONTrackK()
  : AliMUONTrack(),
    fStartSegment(0x0),
    fPosition(0.),
    fPositionNew(0.),
    fChi2(0.),
    fTrackHits(0x0),
    fNmbTrackHits(0),
    fTrackDir(1),
    fBPFlag(kFALSE),
    fRecover(0),
    fSkipHit(0x0),
    fTrackPar(0x0),
    fTrackParNew(0x0),
    fCovariance(0x0),
    fWeight(0x0),
    fParExtrap(0x0),
    fParFilter(0x0), 
    fParSmooth(0x0),
    fCovExtrap(0x0),
    fCovFilter(0x0),
    fJacob(0x0),
    fNSteps(0),
    fSteps(0x0),
    fChi2Array(0x0),
    fChi2Smooth(0x0)
{
/// Default constructor

  fgTrackReconstructor = NULL; // pointer to event reconstructor
  fgHitForRec = NULL; // pointer to points
  fgNOfPoints = 0; // number of points

  return;
}

  //__________________________________________________________________________
AliMUONTrackK::AliMUONTrackK(AliMUONTrackReconstructorK *TrackReconstructor, TClonesArray *hitForRec)
  : AliMUONTrack(),
    fStartSegment(0x0),
    fPosition(0.),
    fPositionNew(0.),
    fChi2(0.),
    fTrackHits(0x0),
    fNmbTrackHits(0),
    fTrackDir(1),
    fBPFlag(kFALSE),
    fRecover(0),
    fSkipHit(0x0),
    fTrackPar(0x0),
    fTrackParNew(0x0),
    fCovariance(0x0),
    fWeight(0x0),
    fParExtrap(0x0),
    fParFilter(0x0), 
    fParSmooth(0x0),
    fCovExtrap(0x0),
    fCovFilter(0x0),
    fJacob(0x0),
    fNSteps(0),
    fSteps(0x0),
    fChi2Array(0x0),
    fChi2Smooth(0x0)
{
/// Constructor

  if (!TrackReconstructor) return;
  fgTrackReconstructor = TrackReconstructor; // pointer to event reconstructor
  fgHitForRec = hitForRec; // pointer to points
  fgNOfPoints = fgHitForRec->GetEntriesFast(); // number of points
  fgCombi = NULL;
  if (fgTrackReconstructor->GetTrackMethod() == 3) fgCombi = AliMUONEventRecoCombi::Instance();
}

  //__________________________________________________________________________
AliMUONTrackK::AliMUONTrackK(AliMUONSegment *segment)
  //: AliMUONTrack(segment, segment, fgTrackReconstructor)
  : AliMUONTrack(NULL, segment),
    fStartSegment(segment),
    fPosition(0.),
    fPositionNew(0.),
    fChi2(0.),
    fTrackHits(new TObjArray(13)),
    fNmbTrackHits(2),
    fTrackDir(1),
    fBPFlag(kFALSE),
    fRecover(0),
    fSkipHit(0x0),
    fTrackPar(new TMatrixD(fgkSize,1)),
    fTrackParNew(new TMatrixD(fgkSize,1)),
    fCovariance(new TMatrixD(fgkSize,fgkSize)),
    fWeight(new TMatrixD(fgkSize,fgkSize)),
    fParExtrap(new TObjArray(15)),
    fParFilter(new TObjArray(15)), 
    fParSmooth(0x0),
    fCovExtrap(new TObjArray(15)),
    fCovFilter(new TObjArray(15)),
    fJacob(new TObjArray(15)),
    fNSteps(0),
    fSteps(new TArrayD(15)),
    fChi2Array(new TArrayD(13)),
    fChi2Smooth(0x0)
{
/// Constructor from a segment
  Double_t dX, dY, dZ;
  AliMUONHitForRec *hit1, *hit2;
  AliMUONRawCluster *clus;
  TClonesArray *rawclusters;

  // Pointers to hits from the segment
  hit1 = segment->GetHitForRec1();
  hit2 = segment->GetHitForRec2();
  hit1->SetNTrackHits(hit1->GetNTrackHits()+1); // mark hit as being on track
  hit2->SetNTrackHits(hit2->GetNTrackHits()+1); // mark hit as being on track
  // check sorting in Z
  if (TMath::Abs(hit1->GetZ()) > TMath::Abs(hit2->GetZ())) {
    hit1 = hit2;
    hit2 = segment->GetHitForRec1();
  }

  // Fill array of track parameters
  if (hit1->GetChamberNumber() > 7) {
    // last tracking station
    (*fTrackPar)(0,0) = hit1->GetBendingCoor(); // y
    (*fTrackPar)(1,0) = hit1->GetNonBendingCoor(); // x
    fPosition = hit1->GetZ(); // z
    fTrackHits->Add(hit2); // add hit 2
    fTrackHits->Add(hit1); // add hit 1
    //AZ(Z->-Z) fTrackDir = -1;
    fTrackDir = 1;
  } else {
    // last but one tracking station
    (*fTrackPar)(0,0) = hit2->GetBendingCoor(); // y
    (*fTrackPar)(1,0) = hit2->GetNonBendingCoor(); // x
    fPosition = hit2->GetZ(); // z
    fTrackHits->Add(hit1); // add hit 1
    fTrackHits->Add(hit2); // add hit 2
    //AZ(Z->-Z) fTrackDir = 1;
    fTrackDir = -1;
  }
  dZ = hit2->GetZ() - hit1->GetZ();
  dY = hit2->GetBendingCoor() - hit1->GetBendingCoor();
  dX = hit2->GetNonBendingCoor() - hit1->GetNonBendingCoor();
  (*fTrackPar)(2,0) = TMath::ATan2(dY,dZ); // alpha
  if ((*fTrackPar)(2,0) < 0.) (*fTrackPar)(2,0) += 2*TMath::Pi(); // from 0 to 2*pi
  (*fTrackPar)(3,0) = TMath::ATan2(-dX,dZ/TMath::Cos((*fTrackPar)(2,0))); // beta
  (*fTrackPar)(2,0) -= TMath::Pi();
  (*fTrackPar)(4,0) = 1/fgTrackReconstructor->GetBendingMomentumFromImpactParam(segment->GetBendingImpact()); // 1/Pt
  (*fTrackPar)(4,0) *= TMath::Cos((*fTrackPar)(3,0)); // 1/p
  // Evaluate covariance (and weight) matrix
  EvalCovariance(dZ);

  if (fgDebug < 0 ) return;
  cout << fgTrackReconstructor->GetNRecTracks()-1 << " " << fgTrackReconstructor->GetBendingMomentumFromImpactParam(segment->GetBendingImpact()) << " " << 1/(*fTrackPar)(4,0) << " ";
    // from raw clusters
  for (Int_t i=0; i<2; i++) {
    hit1 = (AliMUONHitForRec*) ((*fTrackHits)[i]);
    rawclusters = fgTrackReconstructor->GetMUONData()->RawClusters(hit1->GetChamberNumber());
    clus = (AliMUONRawCluster*) rawclusters->UncheckedAt(hit1->GetHitNumber());
    cout << clus->GetTrack(1);
    if (clus->GetTrack(2) != -1) cout << " " << clus->GetTrack(2);
    if (i == 0) cout << " <--> ";
  }
  cout << " @ " << fStartSegment->GetHitForRec1()->GetChamberNumber() << endl;
  
}

  //__________________________________________________________________________
AliMUONTrackK::~AliMUONTrackK()
{
/// Destructor

  if (fTrackHits) {
    //cout << fNmbTrackHits << endl;
    for (Int_t i = 0; i < fNmbTrackHits; i++) {
      AliMUONHitForRec *hit = (AliMUONHitForRec*) fTrackHits->UncheckedAt(i);
      hit->SetNTrackHits(hit->GetNTrackHits()-1);
    }
    delete fTrackHits; // delete the TObjArray of pointers to TrackHit's
    fTrackHits = NULL;
  }
  if (fTrackPar) {
    delete fTrackPar; delete fTrackParNew; delete fCovariance;
    delete fWeight;
  } 

  if (fSteps) delete fSteps; 
  if (fChi2Array) delete fChi2Array;
  if (fChi2Smooth) delete fChi2Smooth;
  if (fParSmooth) {fParSmooth->Delete(); delete fParSmooth; }
  // Delete only matrices not shared by several tracks
  if (!fParExtrap) return;

  Int_t id = 0;
  for (Int_t i=0; i<fNSteps; i++) {
    //if (fParExtrap->UncheckedAt(i)->GetUniqueID() > 1) 
    //  fParExtrap->UncheckedAt(i)->SetUniqueID(fParExtrap->RemoveAt(i)->GetUniqueID()-1);
    id = fParExtrap->UncheckedAt(i)->GetUniqueID();
    if (id > 1) { 
      fParExtrap->UncheckedAt(i)->SetUniqueID(id-1);
      fParExtrap->RemoveAt(i);
    }
    //if (fParFilter->UncheckedAt(i)->GetUniqueID() > 1) 
    //  fParFilter->UncheckedAt(i)->SetUniqueID(fParFilter->RemoveAt(i)->GetUniqueID()-1);
    id = fParFilter->UncheckedAt(i)->GetUniqueID();
    if (id > 1) { 
      fParFilter->UncheckedAt(i)->SetUniqueID(id-1);
      fParFilter->RemoveAt(i);
    }
    //if (fCovExtrap->UncheckedAt(i)->GetUniqueID() > 1) 
    //  fCovExtrap->UncheckedAt(i)->SetUniqueID(fCovExtrap->RemoveAt(i)->GetUniqueID()-1);
    id = fCovExtrap->UncheckedAt(i)->GetUniqueID();
    if (id > 1) { 
      fCovExtrap->UncheckedAt(i)->SetUniqueID(id-1);
      fCovExtrap->RemoveAt(i);
    }

    //if (fCovFilter->UncheckedAt(i)->GetUniqueID() > 1) 
    //  fCovFilter->UncheckedAt(i)->SetUniqueID(fCovFilter->RemoveAt(i)->GetUniqueID()-1);
    id = fCovFilter->UncheckedAt(i)->GetUniqueID();
    if (id > 1) { 
      fCovFilter->UncheckedAt(i)->SetUniqueID(id-1);
      fCovFilter->RemoveAt(i);
    }
    //if (fJacob->UncheckedAt(i)->GetUniqueID() > 1) 
    //  fJacob->UncheckedAt(i)->SetUniqueID(fJacob->RemoveAt(i)->GetUniqueID()-1);
    id = fJacob->UncheckedAt(i)->GetUniqueID();
    if (id > 1) { 
      fJacob->UncheckedAt(i)->SetUniqueID(id-1);
      fJacob->RemoveAt(i);
    }
  }
  /*
  for (Int_t i=0; i<fNSteps; i++) {
    if (fParExtrap->UncheckedAt(i)) ((TMatrixD*)fParExtrap->UncheckedAt(i))->Delete();
    if (fParFilter->UncheckedAt(i)) ((TMatrixD*)fParFilter->UncheckedAt(i))->Delete();
    if (fCovExtrap->UncheckedAt(i)) ((TMatrixD*)fCovExtrap->UncheckedAt(i))->Delete();
    cout << fCovFilter->UncheckedAt(i) << " " << (*fSteps)[i] << endl;
    if (fCovFilter->UncheckedAt(i)) ((TMatrixD*)fCovFilter->UncheckedAt(i))->Delete();
    if (fJacob->UncheckedAt(i)) ((TMatrixD*)fJacob->UncheckedAt(i))->Delete();
  }
  */
  fParExtrap->Delete(); fParFilter->Delete();
  fCovExtrap->Delete(); fCovFilter->Delete();
  fJacob->Delete();
  delete fParExtrap; delete fParFilter;
  delete fCovExtrap; delete fCovFilter;
  delete fJacob;
}

  //__________________________________________________________________________
AliMUONTrackK & AliMUONTrackK::operator=(const AliMUONTrackK& source)
{
/// Assignment operator

  // Members
  if(&source == this) return *this;

  // base class assignement
  //AZ TObject::operator=(source);
  AliMUONTrack::operator=(source);

  fStartSegment = source.fStartSegment;
  fNmbTrackHits = source.fNmbTrackHits;
  fChi2 = source.fChi2;
  fPosition = source.fPosition;
  fPositionNew = source.fPositionNew;
  fTrackDir = source.fTrackDir;
  fBPFlag = source.fBPFlag;
  fRecover = source.fRecover;
  fSkipHit = source.fSkipHit;

  // Pointers
  fTrackHits = new TObjArray(*source.fTrackHits);
  //source.fTrackHits->Dump();
  //fTrackHits->Dump();
  for (Int_t i = 0; i < fNmbTrackHits; i++) {
    AliMUONHitForRec *hit = (AliMUONHitForRec*) fTrackHits->UncheckedAt(i);
    hit->SetNTrackHits(hit->GetNTrackHits()+1);
  }
  
  fTrackPar = new TMatrixD(*source.fTrackPar); // track parameters
  fTrackParNew = new TMatrixD(*source.fTrackParNew); // track parameters
  fCovariance = new TMatrixD(*source.fCovariance); // covariance matrix
  fWeight = new TMatrixD(*source.fWeight); // weight matrix (inverse of covariance)

  // For smoother
  fParExtrap = new TObjArray(*source.fParExtrap);
  fParFilter = new TObjArray(*source.fParFilter);
  fCovExtrap = new TObjArray(*source.fCovExtrap);
  fCovFilter = new TObjArray(*source.fCovFilter);
  fJacob = new TObjArray(*source.fJacob);
  fSteps = new TArrayD(*source.fSteps);
  fNSteps = source.fNSteps;
  fChi2Array = NULL;
  if (source.fChi2Array) fChi2Array = new TArrayD(*source.fChi2Array);
  fChi2Smooth = NULL;
  fParSmooth = NULL;

  for (Int_t i=0; i<fParExtrap->GetEntriesFast(); i++) {
    fParExtrap->UncheckedAt(i)->SetUniqueID(fParExtrap->UncheckedAt(i)->GetUniqueID()+1);
    fParFilter->UncheckedAt(i)->SetUniqueID(fParFilter->UncheckedAt(i)->GetUniqueID()+1);
    fCovExtrap->UncheckedAt(i)->SetUniqueID(fCovExtrap->UncheckedAt(i)->GetUniqueID()+1);
    fCovFilter->UncheckedAt(i)->SetUniqueID(fCovFilter->UncheckedAt(i)->GetUniqueID()+1);
    fJacob->UncheckedAt(i)->SetUniqueID(fJacob->UncheckedAt(i)->GetUniqueID()+1);
  }
  return *this;
}

  //__________________________________________________________________________
void AliMUONTrackK::EvalCovariance(Double_t dZ)
{
/// Evaluate covariance (and weight) matrix for track candidate
  Double_t sigmaB, sigmaNonB, tanA, tanB, dAdY, rad, dBdX, dBdY;

  dZ = -dZ;
  sigmaB = fgTrackReconstructor->GetBendingResolution(); // bending resolution
  sigmaNonB = fgTrackReconstructor->GetNonBendingResolution(); // non-bending resolution

  (*fWeight)(0,0) = sigmaB*sigmaB; // <yy>

  (*fWeight)(1,1) = sigmaNonB*sigmaNonB; // <xx>

  tanA = TMath::Tan((*fTrackPar)(2,0));
  dAdY = 1/(1+tanA*tanA)/dZ;
  (*fWeight)(2,2) = dAdY*dAdY*(*fWeight)(0,0)*2; // <aa>
  (*fWeight)(0,2) = dAdY*(*fWeight)(0,0); // <ya>
  (*fWeight)(2,0) = (*fWeight)(0,2);

  rad = dZ/TMath::Cos((*fTrackPar)(2,0));
  tanB = TMath::Tan((*fTrackPar)(3,0));
  dBdX = 1/(1+tanB*tanB)/rad;
  dBdY = 0; // neglect
  (*fWeight)(3,3) = dBdX*dBdX*(*fWeight)(1,1)*2; // <bb>
  (*fWeight)(1,3) = dBdX*(*fWeight)(1,1); // <xb>
  (*fWeight)(3,1) = (*fWeight)(1,3);

  (*fWeight)(4,4) = ((*fTrackPar)(4,0)*0.5)*((*fTrackPar)(4,0)*0.5); // error 50%

  // check whether the Invert method returns flag if matrix cannot be inverted,
  // and do not calculate the Determinant in that case !!!!
  if (fWeight->Determinant() != 0) {
    // fWeight->Invert();
    Int_t ifail;
    mnvertLocal(&((*fWeight)(0,0)), fgkSize,fgkSize,fgkSize,ifail);
  } else {
    AliWarning(" Determinant fWeight=0:");
  }
  return;
}

  //__________________________________________________________________________
Bool_t AliMUONTrackK::KalmanFilter(Int_t ichamBeg, Int_t ichamEnd, Bool_t Back, Double_t zDipole1, Double_t zDipole2)
{
/// Follows track through detector stations 
  Bool_t miss, success;
  Int_t ichamb, iFB, iMin, iMax, dChamb, ichambOK;
  Int_t ihit, firstIndx = -1, lastIndx = -1, currIndx = -1, iz0 = -1, iz = -1;
  Double_t zEnd, dChi2;
  AliMUONHitForRec *hitAdd, *firstHit = NULL, *lastHit = NULL, *hit;
  AliMUONRawCluster *clus;
  TClonesArray *rawclusters;
  hit = 0; clus = 0; rawclusters = 0;

  miss = success = kTRUE;
  Int_t endOfProp = 0;
  //iFB = (ichamEnd == ichamBeg) ? fTrackDir : TMath::Sign(1,ichamEnd-ichamBeg);
  iFB = (ichamEnd == ichamBeg) ? -fTrackDir : TMath::Sign(1,ichamEnd-ichamBeg);
  iMin = TMath::Min(ichamEnd,ichamBeg);
  iMax = TMath::Max(ichamEnd,ichamBeg);
  ichamb = ichamBeg;
  ichambOK = ichamb;

  if (Back) {
    // backpropagation
    currIndx = 1; 
    if (((AliMUONHitForRec*)fTrackHits->First())->GetChamberNumber() != ichamb) currIndx = 0;
  } else if (fRecover) {
    hit = GetHitLastOk();
    currIndx = fTrackHits->IndexOf(hit);
    if (currIndx < 0) hit = fStartSegment->GetHitForRec1(); // for station 3
    Back = kTRUE;
    ichamb = hit->GetChamberNumber();
    if (hit == fSkipHit || fRecover == 2 && currIndx >= 0) {
      // start from the last point or outlier
      // remove the skipped hit
      fTrackHits->Remove(fSkipHit); // remove hit
      fNmbTrackHits --;
      fSkipHit->SetNTrackHits(fSkipHit->GetNTrackHits()-1); // unmark hit 
      if (fRecover == 1) {
	// recovery
	Back = kFALSE;
	fRecover = 0; 
	ichambOK = ((AliMUONHitForRec*)((*fTrackHits)[fNmbTrackHits-1]))->GetChamberNumber();
	//if (ichambOK >= 8 && ((AliMUONHitForRec*)((*fTrackHits)[0]))->GetChamberNumber() == 6) ichambOK = 6;
	if (fgTrackReconstructor->GetTrackMethod() == 3 && 
	    fSkipHit->GetHitNumber() < 0) {
	  iz0 = fgCombi->IZfromHit(fSkipHit);
	  currIndx = -1;
	}
	else currIndx = fgHitForRec->IndexOf(fSkipHit);
      } else {
	// outlier
	fTrackHits->Compress();
      }
    } // if (hit == fSkipHit)
    else if (currIndx < 0) currIndx = fTrackHits->IndexOf(hit);
  } // else if (fRecover) 
  else {
    // Get indices of the 1'st and last hits on the track candidate
    firstHit = (AliMUONHitForRec*) fTrackHits->First();
    lastHit = (AliMUONHitForRec*) fTrackHits->Last();
    if (fgTrackReconstructor->GetTrackMethod() == 3 && 
	lastHit->GetHitNumber() < 0) iz0 = fgCombi->IZfromHit(lastHit);
    else {
      firstIndx = fgHitForRec->IndexOf(firstHit);
      lastIndx = fgHitForRec->IndexOf(lastHit);
      currIndx = TMath::Abs (TMath::Max(firstIndx*iFB,lastIndx*iFB));
    }
  }

  if (iz0 < 0) iz0 = iFB;
  while (ichamb >= iMin && ichamb <= iMax) {
  // Find the closest hit in Z, not belonging to the current plane
    if (Back) {
      // backpropagation
      if (currIndx < fNmbTrackHits) {
	hitAdd = (AliMUONHitForRec*) fTrackHits->UncheckedAt(currIndx);
	zEnd = hitAdd->GetZ();
	//AZ(z->-z) } else zEnd = -9999;
      } else zEnd = 9999;
    } else {
      //AZ(Z->-Z) zEnd = -9999;
      zEnd = 9999;
      for (ihit=currIndx+iFB; ihit>=0 && ihit<fgNOfPoints; ihit+=iFB) {
	hitAdd = (AliMUONHitForRec*) ((*fgHitForRec)[ihit]);
	//if (TMath::Abs(hitAdd->GetZ()-fPosition) > 0.1) {
	if (TMath::Abs(hitAdd->GetZ()-fPosition) > 0.5) {
	  zEnd = hitAdd->GetZ();
	  currIndx = ihit;
	  break;
	}
      }

      // Combined cluster / track finder
      if (zEnd > 999 && iFB < 0 && fgTrackReconstructor->GetTrackMethod() == 3) {
	currIndx = -2;
	AliMUONHitForRec hitTmp;
	for (iz = iz0 - iFB; iz < fgCombi->Nz(); iz++) {
	  if (TMath::Abs(fgCombi->Z(iz)-fPosition) < 0.5) continue;
	  Int_t *pDEatZ = fgCombi->DEatZ(iz);
	  //cout << iz << " " << fgCombi->Z(iz) << endl;
	  zEnd = fgCombi->Z(iz);
	  iz0 = iz;
	  AliMUONDetElement *detElem = fgCombi->DetElem(pDEatZ[0]);
	  hitAdd = &hitTmp;
	  hitAdd->SetChamberNumber(detElem->Chamber());
	  //hitAdd = (AliMUONHitForRec*) detElem->HitsForRec()->First();
	  if (hitAdd) break;
	}
      }
    }
    if (zEnd>999 && ichamb==ichamEnd) endOfProp = 1; // end-of-propagation
    else {
      // Check if there is a chamber without hits
      if (zEnd>999 || TMath::Abs(hitAdd->GetChamberNumber()-ichamb) > 1) {
	if (!Back && zEnd<999) currIndx -= iFB;
	ichamb += iFB;
	zEnd = AliMUONConstants::DefaultChamberZ(ichamb);
	miss = kTRUE;
      } else {
	ichamb = hitAdd->GetChamberNumber();
	miss = kFALSE;
      }
    }
    if (ichamb<iMin || ichamb>iMax) break;
    // Check for missing station 
    if (!Back) {
      dChamb = TMath::Abs(ichamb-ichambOK); 
      //cout << dChamb << " " << ichambOK << " " << fgNOfPoints << endl;
      Int_t dStatMiss = TMath::Abs (ichamb/2 - ichambOK/2);
      if (zEnd > 999) dStatMiss++;
      if (dStatMiss > 1) {
      //if (dStatMiss == 2 && ichambOK/2 != 3 || dStatMiss > 2) { // AZ - missing st. 3
	// missing station - abandon track
	//cout << dChamb << " " << ichambOK << " " << fgNOfPoints << " " << 1/(*fTrackPar)(4,0) << endl;
	if (fgDebug >= 10) {
	  for (Int_t i1=0; i1<fgNOfPoints; i1++) {
	    cout << " Hit #" << ((AliMUONHitForRec*)((*fgHitForRec)[i1]))->GetChamberNumber() << " ";
	    cout << ((AliMUONHitForRec*)((*fgHitForRec)[i1]))->GetBendingCoor() << " ";
	    cout << ((AliMUONHitForRec*)((*fgHitForRec)[i1]))->GetNonBendingCoor() << " ";
	    cout << ((AliMUONHitForRec*)((*fgHitForRec)[i1]))->GetZ() << " " << " ";
	    cout << ((AliMUONHitForRec*)((*fgHitForRec)[i1]))->GetTTRTrack() << endl;
	  }
	  cout << endl;
	  cout << fNmbTrackHits << endl;
	  for (Int_t i1=0; i1<fNmbTrackHits; i1++) {
	    hit = (AliMUONHitForRec*) ((*fTrackHits)[i1]);
	    printf(" * %d %10.4f %10.4f %10.4f", 
	           hit->GetChamberNumber(), hit->GetBendingCoor(), 
		   hit->GetNonBendingCoor(), hit->GetZ());
	    // from raw clusters
	    rawclusters = fgTrackReconstructor->GetMUONData()->RawClusters(hit->GetChamberNumber());
	    clus = (AliMUONRawCluster*) rawclusters->UncheckedAt(hit->GetHitNumber());
	    printf("%3d", clus->GetTrack(1)-1); 
	    if (clus->GetTrack(2) != 0)
	      printf("%3d \n", clus->GetTrack(2)-1);
	    else 
	      printf("\n");
	    
	  }
	} // if (fgDebug >= 10) 
	if (fNmbTrackHits>2 && fRecover==0) Recover(); // try to recover track later
	return kFALSE;
      } // if (dStatMiss > 1)
    } // if (!Back)
    if (endOfProp != 0) break;

    // propagate to the found Z

    // Check if track steps into dipole
    //AZ(z->-z) if (fPosition>zDipole2 && zEnd<zDipole2) {
    if (fPosition<zDipole2 && zEnd>zDipole2) {
      //LinearPropagation(zDipole2-zBeg); 
      ParPropagation(zDipole2); 
      MSThin(1); // multiple scattering in the chamber
      WeightPropagation(zDipole2, kTRUE); // propagate weight matrix
      fPosition = fPositionNew;
      *fTrackPar = *fTrackParNew; 
      //MagnetPropagation(zEnd); 
      ParPropagation(zEnd); 
      WeightPropagation(zEnd, kTRUE);
      fPosition = fPositionNew;
    } 
    // Check if track steps out of dipole
    //AZ(z->-z) else if (fPosition>zDipole1 && zEnd<zDipole1) {
    else if (fPosition<zDipole1 && zEnd>zDipole1) {
      //MagnetPropagation(zDipole1-zBeg); 
      ParPropagation(zDipole1); 
      MSThin(1); // multiple scattering in the chamber
      WeightPropagation(zDipole1, kTRUE);
      fPosition = fPositionNew;
      *fTrackPar = *fTrackParNew; 
      //LinearPropagation(zEnd-zDipole1); 
      ParPropagation(zEnd); 
      WeightPropagation(zEnd, kTRUE);
      fPosition = fPositionNew;
    } else {
      ParPropagation(zEnd);
      //MSThin(1); // multiple scattering in the chamber
      if (TMath::Abs(zEnd-fPosition) > 5) MSThin(1); // multiple scattering in the chamber
      WeightPropagation(zEnd, kTRUE);
      fPosition = fPositionNew;
    }

    // Add measurement
    if (fRecover != 0 && hitAdd == fSkipHit && !miss) {
      // recovered track - remove the hit
      miss = kTRUE;
      ichamb = fSkipHit->GetChamberNumber();
      // remove the skipped hit
      fTrackHits->Remove(fSkipHit); 
      fNmbTrackHits --;
      //AZ fSkipHit->SetNTrackHits(fSkipHit->GetNTrackHits()-1); // unmark hit 
      Back = kFALSE;
      fRecover = 0;
      ichambOK = ((AliMUONHitForRec*)((*fTrackHits)[fNmbTrackHits-1]))->GetChamberNumber();
      //? if (fgTrackReconstructor->GetTrackMethod() != 3) currIndx = fgHitForRec->IndexOf(fSkipHit);
      currIndx = fgHitForRec->IndexOf(fSkipHit);
    }

    if (Back && !miss) {
      // backward propagator
      if (currIndx) {
	TMatrixD pointWeight(fgkSize,fgkSize);
	TMatrixD point(fgkSize,1);
	TMatrixD trackParTmp = point;
	point(0,0) = hitAdd->GetBendingCoor();
	point(1,0) = hitAdd->GetNonBendingCoor();
	pointWeight(0,0) = 1/hitAdd->GetBendingReso2();
	pointWeight(1,1) = 1/hitAdd->GetNonBendingReso2();
	TryPoint(point,pointWeight,trackParTmp,dChi2);
	*fTrackPar = trackParTmp;
	*fWeight += pointWeight; 
	if (fTrackDir > 0) AddMatrices (this, dChi2, hitAdd);
	fChi2 += dChi2; // Chi2
      } else *fTrackPar = *fTrackParNew; // adjust end point to the last hit
      if (ichamb==ichamEnd) break; 
      currIndx ++;
    } else {
      // forward propagator
      if (miss || !FindPoint(ichamb, zEnd, currIndx, iFB, hitAdd, iz)) {
	// missing point
	*fTrackPar = *fTrackParNew; 
      } else {
	//add point
	fTrackHits->Add(hitAdd); // add hit
	fNmbTrackHits ++;
	hitAdd->SetNTrackHits(hitAdd->GetNTrackHits()+1); // mark hit as being on track
	ichambOK = ichamb;
	currIndx = fgHitForRec->IndexOf(hitAdd); // Check
      }
    }
  } // while
  if (fgDebug > 0) cout << fNmbTrackHits << " " << fChi2 << " " << 1/(*fTrackPar)(4,0) << " " << fPosition << endl;
  if (1/TMath::Abs((*fTrackPar)(4,0)) < fgTrackReconstructor->GetMinBendingMomentum()) success = kFALSE; // p < p_min
  if (GetRecover() < 0) success = kFALSE;
  return success;
}

  //__________________________________________________________________________
void AliMUONTrackK::ParPropagation(Double_t zEnd)
{
/// Propagation of track parameters to zEnd
  Int_t iFB, nTries;
  Double_t dZ, step, distance, charge;
  Double_t vGeant3[7], vGeant3New[7];
  AliMUONTrackParam trackParam;

  nTries = 0;
  // First step using linear extrapolation
  dZ = zEnd - fPosition;
  fPositionNew = fPosition;
  *fTrackParNew = *fTrackPar;
  if (dZ == 0) return; //AZ ???
  iFB = (Int_t)TMath::Sign(Double_t(1.0),dZ);
  step = dZ/TMath::Cos((*fTrackPar)(2,0))/TMath::Cos((*fTrackPar)(3,0)); // linear estimate
  //AZ(z->-z) charge = iFB*TMath::Sign(Double_t(1.0),(*fTrackPar)(4,0));
  charge = -iFB*TMath::Sign(Double_t(1.0),(*fTrackPar)(4,0));
  SetGeantParam(vGeant3,iFB);
  //fTrackParNew->Print();

  // Check if overstep
  do {
    step = TMath::Abs(step);
    // Propagate parameters
    AliMUONTrackExtrap::ExtrapOneStepRungekutta(charge,step,vGeant3,vGeant3New);
    //extrap_onestep_rungekutta(charge,step,vGeant3,vGeant3New);
    distance = zEnd - vGeant3New[2];
    step *= dZ/(vGeant3New[2]-fPositionNew);
    nTries ++;
  } while (distance*iFB < 0 && TMath::Abs(distance) > fgkEpsilon);

  GetFromGeantParam(vGeant3New,iFB);
  //fTrackParNew->Print();

  // Position adjustment (until within tolerance)
  while (TMath::Abs(distance) > fgkEpsilon) {
    dZ = zEnd - fPositionNew;
    iFB = (Int_t)TMath::Sign(Double_t(1.0),dZ);
    step = dZ/TMath::Cos((*fTrackParNew)(2,0))/TMath::Cos((*fTrackParNew)(3,0));
    step = TMath::Abs(step);
    SetGeantParam(vGeant3,iFB);
    do {
      // binary search
      // Propagate parameters
      AliMUONTrackExtrap::ExtrapOneStepRungekutta(charge,step,vGeant3,vGeant3New);
      //extrap_onestep_rungekutta(charge,step,vGeant3,vGeant3New);
      distance = zEnd - vGeant3New[2];
      step /= 2;
      nTries ++;
      if (nTries > fgkTriesMax) AliError(Form(" Too many tries: %d", nTries));
    } while (distance*iFB < 0);

    GetFromGeantParam(vGeant3New,iFB);
  }
  //cout << nTries << endl;
  //fTrackParNew->Print();
  return;
}

  //__________________________________________________________________________
void AliMUONTrackK::WeightPropagation(Double_t zEnd, Bool_t smooth)
{
/// Propagation of the weight matrix
/// W = DtWD, where D is Jacobian 
  Int_t i, j;
  Double_t dPar;

  TMatrixD jacob(fgkSize,fgkSize);
  jacob = 0;

  // Save initial and propagated parameters
  TMatrixD trackPar0 = *fTrackPar;
  TMatrixD trackParNew0 = *fTrackParNew;

  // Get covariance matrix
  *fCovariance = *fWeight;
  // check whether the Invert method returns flag if matrix cannot be inverted,
  // and do not calculate the Determinant in that case !!!!
  if (fCovariance->Determinant() != 0) {
    Int_t ifail;
    mnvertLocal(&((*fCovariance)(0,0)), fgkSize,fgkSize,fgkSize,ifail);
    //fCovariance->Print();
  } else {
    AliWarning(" Determinant fCovariance=0:");
  }

  // Loop over parameters to find change of the propagated vs initial ones
  for (i=0; i<fgkSize; i++) {
    dPar = TMath::Sqrt((*fCovariance)(i,i));
    *fTrackPar = trackPar0;
    if (i == 4) dPar *= TMath::Sign(1.,-trackPar0(4,0)); // 1/p
    (*fTrackPar)(i,0) += dPar;
    ParPropagation(zEnd);
    for (j=0; j<fgkSize; j++) {
      jacob(j,i) = ((*fTrackParNew)(j,0)-trackParNew0(j,0))/dPar;
    }
  }

  //trackParNew0.Print();
  //TMatrixD par1(jacob,TMatrixD::kMult,trackPar0); //
  //par1.Print();
  TMatrixD jacob0 = jacob;
  if (jacob.Determinant() != 0) {
    jacob.Invert();
  } else {
    AliWarning(" Determinant jacob=0:");
  }
  TMatrixD weight1(*fWeight,TMatrixD::kMult,jacob); // WD
  *fWeight = TMatrixD(jacob,TMatrixD::kTransposeMult,weight1); // DtWD

  // Restore initial and propagated parameters
  *fTrackPar = trackPar0;
  *fTrackParNew = trackParNew0;

  // Save for smoother
  if (!smooth) return; // do not use smoother
  if (fTrackDir < 0) return; // only for propagation towards int. point
  TMatrixD *tmp = new TMatrixD(*fTrackParNew); // extrapolated parameters
  fParExtrap->Add(tmp);
  tmp->SetUniqueID(1);
  tmp = new TMatrixD(*fTrackParNew); // filtered parameters (if no measurement will be found)
  fParFilter->Add(tmp);
  tmp->SetUniqueID(1);
  *fCovariance = *fWeight;
  if (fCovariance->Determinant() != 0) {
    Int_t ifail;
    mnvertLocal(&((*fCovariance)(0,0)), fgkSize,fgkSize,fgkSize,ifail);
  } else {
    AliWarning(" Determinant fCovariance=0:");
  }
  tmp = new TMatrixD(*fCovariance); // extrapolated covariance
  fCovExtrap->Add(tmp);
  tmp->SetUniqueID(1);
  tmp = new TMatrixD(*fCovariance); // filtered covariance (if no measurement will be found)
  fCovFilter->Add(tmp);
  tmp->SetUniqueID(1);
  tmp = new TMatrixD(jacob0); // Jacobian
  fJacob->Add(tmp);
  tmp->SetUniqueID(1);
  if (fSteps->fN <= fNSteps) fSteps->Set(fSteps->fN+10);
  fSteps->AddAt(fPositionNew,fNSteps++); 
  if (fgDebug > 0) cout << " WeightPropagation " << fNSteps << " " << fPositionNew << endl;
  return;
}

  //__________________________________________________________________________
Bool_t AliMUONTrackK::FindPoint(Int_t ichamb, Double_t zEnd, Int_t currIndx, Int_t iFB, AliMUONHitForRec *&hitAdd, Int_t iz)
{
/// Picks up point within a window for the chamber No ichamb 
/// Split the track if there are more than 1 hit
  Int_t ihit, nRecTracks;
  Double_t windowB, windowNonB, dChi2Tmp=0, dChi2, y, x, savePosition=0;
  TClonesArray *trackPtr;
  AliMUONHitForRec *hit, *hitLoop;
  AliMUONTrackK *trackK;
  AliMUONDetElement *detElem = NULL;

  Bool_t ok = kFALSE;

  if (fgTrackReconstructor->GetTrackMethod() == 3 && iz >= 0) {
    // Combined cluster / track finder
    // Check if extrapolated track passes thru det. elems. at iz
    Int_t *pDEatZ = fgCombi->DEatZ(iz);
    Int_t nDetElem = pDEatZ[-1];
    //cout << fgCombi->Z(iz) << " " << nDetElem << endl;
    for (Int_t i = 0; i < nDetElem; i++) {
      detElem = fgCombi->DetElem(pDEatZ[i]);
      if (detElem->Inside((*fTrackParNew)(1,0), (*fTrackParNew)(0,0), fPosition)) {
	detElem->ClusterReco((*fTrackParNew)(1,0), (*fTrackParNew)(0,0));
	hitAdd = (AliMUONHitForRec*) detElem->HitsForRec()->First();
	ok = kTRUE;
	break;
      }
    }
    if (!ok) return ok; // outside det. elem. 
    ok = kFALSE;
  }

  //sigmaB = fgTrackReconstructor->GetBendingResolution(); // bending resolution
  //sigmaNonB = fgTrackReconstructor->GetNonBendingResolution(); // non-bending resolution
  *fCovariance = *fWeight;
  // check whether the Invert method returns flag if matrix cannot be inverted,
  // and do not calculate the Determinant in that case !!!!
  if (fCovariance->Determinant() != 0) {
      Int_t ifail;
      mnvertLocal(&((*fCovariance)(0,0)), fgkSize,fgkSize,fgkSize,ifail);
  } else {
    AliWarning(" Determinant fCovariance=0:");
  }
  //windowB = fgkNSigma*TMath::Sqrt((*fCovariance)(0,0)+sigmaB*sigmaB);
  //windowNonB = fgkNSigma*TMath::Sqrt((*fCovariance)(1,1)+sigmaNonB*sigmaNonB);
  // Loop over all hits and take hits from the chamber
  TMatrixD pointWeight(fgkSize,fgkSize);
  TMatrixD saveWeight = pointWeight;
  TMatrixD pointWeightTmp = pointWeight;
  TMatrixD point(fgkSize,1);
  TMatrixD trackPar = point;
  TMatrixD trackParTmp = point;
  Int_t nHitsOK = 0, ihitB = currIndx, ihitE = fgNOfPoints, iDhit = iFB;
  Double_t zLast;
  zLast = ((AliMUONHitForRec*)fTrackHits->Last())->GetZ();
  if (fgTrackReconstructor->GetTrackMethod() == 3 && detElem) {
    ihitB = 0;
    ihitE = detElem->NHitsForRec();
    iDhit = 1;
  }

  TArrayD branchChi2(20);
  for (ihit = ihitB; ihit >= 0 && ihit < ihitE; ihit+=iDhit) {
    if (fgTrackReconstructor->GetTrackMethod() != 3 || !detElem) 
      hit = (AliMUONHitForRec*) ((*fgHitForRec)[ihit]);
    else hit = (AliMUONHitForRec*) detElem->HitsForRec()->UncheckedAt(ihit);
    if (hit->GetChamberNumber() == ichamb) {
      //if (TMath::Abs(hit->GetZ()-zEnd) < 0.1) {
      if (TMath::Abs(hit->GetZ()-zEnd) < 0.5) {
	y = hit->GetBendingCoor();
	x = hit->GetNonBendingCoor();
	if (hit->GetBendingReso2() < 0) {
	  // Combined cluster / track finder
	  hit->SetBendingReso2(fgTrackReconstructor->GetBendingResolution()*
			       fgTrackReconstructor->GetBendingResolution());
	  hit->SetNonBendingReso2(fgTrackReconstructor->GetNonBendingResolution()*
				  fgTrackReconstructor->GetNonBendingResolution());
	}
	windowB = fgkNSigma*TMath::Sqrt((*fCovariance)(0,0)+hit->GetBendingReso2());
	windowNonB = fgkNSigma*TMath::Sqrt((*fCovariance)(1,1)+hit->GetNonBendingReso2());

	// windowB = TMath::Min (windowB,5.);
	if (fgkNSigma > 6) windowB = TMath::Min (windowB,5.);
	/*
	if (TMath::Abs(hit->GetZ()-zLast) < 50) {
	  windowB = TMath::Min (windowB,0.5);
	  windowNonB = TMath::Min (windowNonB,3.);
	} else if (TMath::Abs(hit->GetZ()-zLast) < 200) {
	  windowB = TMath::Min (windowB,1.5);
	  windowNonB = TMath::Min (windowNonB,3.);
	} else if (TMath::Abs(hit->GetZ()-zLast) < 350) {
	  windowB = TMath::Min (windowB,4.);
	  windowNonB = TMath::Min (windowNonB,6.);
	}
	*/

	// Test
	/*
	if (TMath::Abs(hit->GetZ()-zLast) < 50) {
	  windowB = 5*TMath::Sqrt((*fCovariance)(0,0)+hit->GetBendingReso2());
	  windowNonB = 5*TMath::Sqrt((*fCovariance)(1,1)+hit->GetNonBendingReso2());
	} else {
	  windowB = fgkNSigma*TMath::Sqrt((*fCovariance)(0,0)+hit->GetBendingReso2());
	  windowNonB = fgkNSigma*TMath::Sqrt((*fCovariance)(1,1)+hit->GetNonBendingReso2());
	}
	*/

	//cout << TMath::Abs((*fTrackParNew)(0,0)-y) << " " << windowB << " " << TMath::Abs((*fTrackParNew)(1,0)-x) << " " << windowNonB << " " << zEnd << endl;
	if (TMath::Abs((*fTrackParNew)(0,0)-y) <= windowB &&
	    TMath::Abs((*fTrackParNew)(1,0)-x) <= windowNonB) {
	//if (TMath::Abs((*fTrackParNew)(0,0)-y) <= windowB &&
	  //    TMath::Abs((*fTrackParNew)(1,0)-x) <= windowNonB &&
	  //  hit->GetTrackRefSignal() == 1) { // just for test
	  // Vector of measurements and covariance matrix
	  //fprintf(lun1,"%3d %3d %10.4f %10.4f \n", gAlice->GetEvNumber(), ichamb, x, y);
	  if (TMath::Abs(hit->GetZ()-zEnd) > 0.05) {
	    // Adjust position: for multiple hits in the chamber or misalignment (Z as a function of X or Y)
	    //AliWarning(Form(" *** adjust %f %f ", zEnd, hit->GetZ()));
	    zEnd = hit->GetZ();
	    *fTrackPar = *fTrackParNew;
	    ParPropagation(zEnd);
	    WeightPropagation(zEnd, kTRUE);
	    fPosition = fPositionNew;
	    *fTrackPar = *fTrackParNew;
	    // Get covariance
	    *fCovariance = *fWeight;
	    if (fCovariance->Determinant() != 0) {
	      Int_t ifail;
	      mnvertLocal(&((*fCovariance)(0,0)), fgkSize,fgkSize,fgkSize,ifail);
	    } else {
	      AliWarning(" Determinant fCovariance=0:" );
	    }
	  }
	  point.Zero();
	  point(0,0) = y;
	  point(1,0) = x;
	  pointWeight(0,0) = 1/hit->GetBendingReso2();
	  pointWeight(1,1) = 1/hit->GetNonBendingReso2();
	  TryPoint(point,pointWeight,trackPar,dChi2);
	  if (TMath::Abs(1./(trackPar)(4,0)) < fgTrackReconstructor->GetMinBendingMomentum()) continue; // p < p_min - next hit
	  // if (TMath::Sign(1.,(trackPar)(4,0)*(*fTrackPar)(4,0)) < 0) continue; // change of sign
	  ok = kTRUE;
	  nHitsOK++;
	  //if (nHitsOK > -1) {
	  if (nHitsOK == 1) {
	    // Save current members
	    saveWeight = *fWeight;
	    savePosition = fPosition;
	    // temporary storage for the current track
	    dChi2Tmp = dChi2;
	    trackParTmp = trackPar;
	    pointWeightTmp = pointWeight;
	    hitAdd = hit;
	    if (fgDebug > 0) cout << " Added point: " << x << " " << y << " " << dChi2 << endl;
	    branchChi2[0] = dChi2;
	  } else {
	    // branching: create a new track
	    trackPtr = fgTrackReconstructor->GetRecTracksPtr();
	    nRecTracks = fgTrackReconstructor->GetNRecTracks();
	    trackK = new ((*trackPtr)[nRecTracks]) AliMUONTrackK(NULL, NULL); 
	    *trackK = *this;
	    fgTrackReconstructor->SetNRecTracks(nRecTracks+1);
	    if (fgDebug > 0) cout << " ******** New track: " << ichamb << " " << hit->GetTTRTrack() << " " << 1/(trackPar)(4,0) << " " << hit->GetBendingCoor() << " " << hit->GetNonBendingCoor() << " " << fNmbTrackHits << " " << nRecTracks << endl;
	    trackK->fRecover = 0;
	    *(trackK->fTrackPar) = trackPar;
	    *(trackK->fWeight) += pointWeight; 
	    trackK->fChi2 += dChi2;
	    // Check
	    /*
	    *(trackK->fCovariance) = *(trackK->fWeight);
	    if (trackK->fCovariance->Determinant() != 0) {
	      Int_t ifail;
	      mnvertLocal(&((*(trackK->fCovariance))(0,0)), fgkSize,fgkSize,fgkSize,ifail);
	    }
	    cout << (*(trackK->fCovariance))(0,0) << " " << (*(trackK->fCovariance))(1,1) << " " << (*fCovariance)(0,0) << " " << (*fCovariance)(1,1) << endl;
	    */
	    // Add filtered matrices for the smoother
	    if (fTrackDir > 0) {
	      if (nHitsOK > 2) { // check for position adjustment
		for (Int_t i=trackK->fNSteps-1; i>=0; i--) {
		  if (TMath::Abs(hit->GetZ()-(*trackK->fSteps)[i]) > 0.1) {
		  //if (TMath::Abs(hit->GetZ()-(trackK->fSteps)[i]) > 0.1) {
		    RemoveMatrices(trackK);
		    AliError(Form(" *** Position adjustment 1: %f %f", hit->GetZ(), (*trackK->fSteps)[i]));
		  }
		  else break;
		}
	      }
	      AddMatrices (trackK, dChi2, hit);
	    }
	    // Mark hits as being on 2 tracks
	    for (Int_t i=0; i<fNmbTrackHits; i++) {
	      hitLoop = (AliMUONHitForRec*) ((*fTrackHits)[i]);
	      //AZ hitLoop->SetNTrackHits(hitLoop->GetNTrackHits()+1); 
	      if (fgDebug >=10) {
		cout << " ** ";
		cout << hitLoop->GetChamberNumber() << " ";
		cout << hitLoop->GetBendingCoor() << " ";
		cout << hitLoop->GetNonBendingCoor() << " ";
		cout << hitLoop->GetZ() << " " << " ";
		cout << hitLoop->GetTTRTrack() << endl;
		printf(" ** %d %10.4f %10.4f %10.4f\n", 
		       hitLoop->GetChamberNumber(), hitLoop->GetBendingCoor(), 
		       hitLoop->GetNonBendingCoor(), hitLoop->GetZ());
	      }
	    }
	    //add point
	    trackK->fTrackHits->Add(hit); // add hit
	    trackK->fNmbTrackHits ++;
	    hit->SetNTrackHits(hit->GetNTrackHits()+1); // mark hit as being on track
	    if (ichamb == 9) {
	      // the last chamber
	      trackK->fTrackDir = 1;
	      trackK->fBPFlag = kTRUE; 
	    }
	    if (nHitsOK > branchChi2.GetSize()) branchChi2.Set(branchChi2.GetSize()+10);
	    branchChi2[nHitsOK-1] = dChi2;
	  }
	}
      }
    } else break; // different chamber
  } // for (ihit=currIndx;
  if (ok) {
    // Restore members
    *fTrackPar = trackParTmp;
    *fWeight = saveWeight;
    *fWeight += pointWeightTmp; 
    fChi2 += dChi2Tmp; // Chi2
    fPosition = savePosition;
    // Add filtered matrices for the smoother
    if (fTrackDir > 0) {
      for (Int_t i=fNSteps-1; i>=0; i--) {
	if (TMath::Abs(fPosition-(*fSteps)[i]) > 0.1) {
	//if (TMath::Abs(fPosition-fSteps[i]) > 0.1) {
	  RemoveMatrices(this);
	  if (fgDebug > 0) cout << " *** Position adjustment 2 " << fPosition << " " << (*fSteps)[i] << endl;
	}
	else break;
      } // for (Int_t i=fNSteps-1;
      AddMatrices (this, dChi2Tmp, hitAdd);
      /*
      if (nHitsOK > 1) {
	for (Int_t i=0; i<fNSteps; i++) cout << (*fSteps)[i] << " "; cout << endl;
	for (Int_t i=0; i<trackK->fNSteps; i++) cout << (*trackK->fSteps)[i] << " "; cout << endl;
      }
      */
    } // if (fTrackDir > 0)
    // Check for maximum number of branches - exclude excessive
    if (nHitsOK > 1) CheckBranches(branchChi2, nHitsOK); 
  }
  return ok;
}

  //__________________________________________________________________________
void AliMUONTrackK::CheckBranches(TArrayD &branchChi2, Int_t nBranch)
{
/// Check for maximum number of branches - exclude excessive

  Int_t nBranchMax = 5;
  if (nBranch <= nBranchMax) return;

  Double_t *chi2 = branchChi2.GetArray();
  Int_t *indx = new Int_t [nBranch];
  TMath::Sort (nBranch, chi2, indx, kFALSE);
  TClonesArray *trackPtr = fgTrackReconstructor->GetRecTracksPtr();
  Int_t nRecTracks = fgTrackReconstructor->GetNRecTracks();
  Int_t ibeg = nRecTracks - nBranch;

  // Discard excessive branches with higher Chi2 contribution
  for (Int_t i = nBranchMax; i < nBranch; ++i) {
    if (indx[i] == 0) {
      // Discard current track
      SetRecover(-1);
      continue;
    }
    Int_t j = ibeg + indx[i];
    AliMUONTrackK *trackK = (AliMUONTrackK*) trackPtr->UncheckedAt(j);
    trackK->SetRecover(-1);
  }
  delete [] indx;
}

  //__________________________________________________________________________
void AliMUONTrackK::TryPoint(TMatrixD &point, const TMatrixD &pointWeight, TMatrixD &trackParTmp, Double_t &dChi2)
{
/// Adds a measurement point (modifies track parameters and computes
/// change of Chi2)

  // Solving linear system (W+U)p' = U(m-p) + (W+U)p
  TMatrixD wu = *fWeight;
  wu += pointWeight; // W+U
  trackParTmp = point;
  trackParTmp -= *fTrackParNew; // m-p
  TMatrixD right(pointWeight,TMatrixD::kMult,trackParTmp); // U(m-p)
  TMatrixD right1(wu,TMatrixD::kMult,*fTrackParNew); // (W+U)p
  right += right1; // U(m-p) + (W+U)p

  // check whether the Invert method returns flag if matrix cannot be inverted,
  // and do not calculate the Determinant in that case !!!!
  if (wu.Determinant() != 0) {
    Int_t ifail;
    mnvertLocal(&((wu)(0,0)), fgkSize,fgkSize,fgkSize,ifail);
  } else {
    AliWarning(" Determinant wu=0:");
  }
  trackParTmp = TMatrixD(wu,TMatrixD::kMult,right); 

  right1 = trackParTmp;
  right1 -= point; // p'-m
  point = trackParTmp;
  point -= *fTrackParNew; // p'-p
  right = TMatrixD(*fWeight,TMatrixD::kMult,point); // W(p'-p)
  TMatrixD value(point,TMatrixD::kTransposeMult,right); // (p'-p)'W(p'-p)
  dChi2 = value(0,0);
  right = TMatrixD(pointWeight,TMatrixD::kMult,right1); // U(p'-m)
  value = TMatrixD(right1,TMatrixD::kTransposeMult,right); // (p'-m)'U(p'-m)
  dChi2 += value(0,0);
  return;
}

  //__________________________________________________________________________
void AliMUONTrackK::MSThin(Int_t sign)
{
/// Adds multiple scattering in a thin layer (only angles are affected)
  Double_t cosAlph, cosBeta, momentum, velo, path, theta0;

  // check whether the Invert method returns flag if matrix cannot be inverted,
  // and do not calculate the Determinant in that case !!!!
  if (fWeight->Determinant() != 0) {
    Int_t ifail;
    mnvertLocal(&((*fWeight)(0,0)), fgkSize,fgkSize,fgkSize,ifail);
  } else {
    AliWarning(" Determinant fWeight=0:");
  }

  cosAlph = TMath::Cos((*fTrackParNew)(2,0));
  cosBeta = TMath::Cos((*fTrackParNew)(3,0));
  momentum = 1/(*fTrackParNew)(4,0); // particle momentum
  //velo = momentum/TMath::Sqrt(momentum*momentum+muonMass*muonMass); // velocity/c for muon hypothesis
  velo = 1; // relativistic
  path = TMath::Abs(fgTrackReconstructor->GetChamberThicknessInX0()/cosAlph/cosBeta); // path length
  theta0 = 0.0136/velo/momentum*TMath::Sqrt(path)*(1+0.038*TMath::Log(path)); // projected scattering angle

  (*fWeight)(2,2) += sign*theta0/cosBeta*theta0/cosBeta; // alpha
  (*fWeight)(3,3) += sign*theta0*theta0; // beta
  Int_t ifail;
  mnvertLocal(&((*fWeight)(0,0)), fgkSize,fgkSize,fgkSize,ifail);
  return;
}

  //__________________________________________________________________________
void AliMUONTrackK::StartBack(void)
{
/// Starts backpropagator
  
  fBPFlag = kTRUE;
  fChi2 = 0;
  for (Int_t i=0; i<fgkSize; i++) {
    for (Int_t j=0; j<fgkSize; j++) {
      if (j==i) (*fWeight)(i,i) /= 100;
      //if (j==i) (*fWeight)(i,i) /= fNmbTrackHits*fNmbTrackHits;
      else (*fWeight)(j,i) = 0;
    }
  }
  // Sort hits on track in descending order in abs(z)
  SortHits(0, fTrackHits);
}

  //__________________________________________________________________________
void AliMUONTrackK::SortHits(Int_t iflag, TObjArray *array)
{
/// Sort hits in Z if the seed segment in the last but one station
/// (if iflag==0 in descending order in abs(z), if !=0 - unsort)
  
  if (iflag && ((AliMUONHitForRec*)(array->UncheckedAt(0)))->GetChamberNumber() == 6) return;
  Double_t z = 0, zmax = TMath::Abs(((AliMUONHitForRec*)(array->UncheckedAt(0)))->GetZ());
  Int_t i = 1, entries = array->GetEntriesFast(); 
  for ( ; i<entries; i++) {
    if (iflag) {
      if (((AliMUONHitForRec*)(array->UncheckedAt(i)))->GetChamberNumber() == 6) break;
    } else {
      z = TMath::Abs(((AliMUONHitForRec*)(array->UncheckedAt(i)))->GetZ());
      if (z < zmax) break;
      zmax = z;
      if (fgDebug >= 10) cout << " " << zmax << " " << z << " " << i << endl;
    }
  }
  if (!iflag) i--;
  for (Int_t j=0; j<=(i-1)/2; j++) {
    TObject *hit = array->UncheckedAt(j);
    array->AddAt(array->UncheckedAt(i-j),j);
    array->AddAt(hit,i-j);
  }
  if (fgDebug >= 10) {
    for (i=0; i<entries; i++) 
      cout << ((AliMUONHitForRec*)(array->UncheckedAt(i)))->GetZ() << " ";
    cout << " - Sort" << endl;
  }
}

  //__________________________________________________________________________
void AliMUONTrackK::SetGeantParam(Double_t *VGeant3, Int_t iFB)
{
/// Set vector of Geant3 parameters pointed to by "VGeant3"
/// from track parameters 

  VGeant3[0] = (*fTrackParNew)(1,0); // X
  VGeant3[1] = (*fTrackParNew)(0,0); // Y
  VGeant3[2] = fPositionNew; // Z
  VGeant3[3] = iFB*TMath::Sin((*fTrackParNew)(3,0)); // Px/Ptot
  VGeant3[4] = iFB*TMath::Cos((*fTrackParNew)(3,0))*TMath::Sin((*fTrackParNew)(2,0)); // Py/Ptot
  VGeant3[5] = iFB*TMath::Sqrt(1.0-VGeant3[3]*VGeant3[3]-VGeant3[4]*VGeant3[4]); // Pz/Ptot
  VGeant3[6] = 1/TMath::Abs((*fTrackParNew)(4,0)); // Ptot
}

  //__________________________________________________________________________
void AliMUONTrackK::GetFromGeantParam(Double_t *VGeant3, Int_t iFB)
{
/// Get track parameters from vector of Geant3 parameters pointed 
/// to by "VGeant3"

  fPositionNew = VGeant3[2]; // Z
  (*fTrackParNew)(0,0) = VGeant3[1]; // Y 
  (*fTrackParNew)(1,0) = VGeant3[0]; // X
  (*fTrackParNew)(3,0) = TMath::ASin(iFB*VGeant3[3]); // beta
  (*fTrackParNew)(2,0) = TMath::ASin(iFB*VGeant3[4]/TMath::Cos((*fTrackParNew)(3,0))); // alpha
  (*fTrackParNew)(4,0) = 1/VGeant3[6]*TMath::Sign(Double_t(1.0),(*fTrackPar)(4,0)); // 1/Ptot
}

  //__________________________________________________________________________
void AliMUONTrackK::SetTrackQuality(Int_t iChi2)
{
/// Computes "track quality" from Chi2 (if iChi2==0) or vice versa

  if (fChi2 > 500) {
    AliWarning(Form(" *** Too high Chi2: %f ", fChi2));
    fChi2 = 500;
  }
  if (iChi2 == 0) fChi2 = fNmbTrackHits + (500.-fChi2)/501;
  else fChi2 = 500 - (fChi2-fNmbTrackHits)*501;
}

  //__________________________________________________________________________
Int_t AliMUONTrackK::Compare(const TObject* trackK) const
{
/// "Compare" function to sort with decreasing "track quality".
/// Returns +1 (0, -1) if quality of current track
/// is smaller than (equal to, larger than) quality of trackK

  if (fChi2 < ((AliMUONTrackK*)trackK)->fChi2) return(+1);
  else if (fChi2 == ((AliMUONTrackK*)trackK)->fChi2) return(0);
  else return(-1);
}

  //__________________________________________________________________________
Bool_t AliMUONTrackK::KeepTrack(AliMUONTrackK* track0) const
{
/// Check whether or not to keep current track 
/// (keep, if it has less than half of common hits with track0)
  Int_t hitsInCommon, nHits0, i, j, nTrackHits2;
  AliMUONHitForRec *hit0, *hit1;

  hitsInCommon = 0;
  nHits0 = track0->fNmbTrackHits;
  nTrackHits2 = fNmbTrackHits/2;

  for (i=0; i<nHits0; i++) {
    // Check if hit belongs to several tracks
    hit0 = (AliMUONHitForRec*) (*track0->fTrackHits)[i]; 
    if (hit0->GetNTrackHits() == 1) continue; 
    for (j=0; j<fNmbTrackHits; j++) {
      hit1 = (AliMUONHitForRec*) (*fTrackHits)[j]; 
      if (hit1->GetNTrackHits() == 1) continue; 
      if (hit0 == hit1) {
	hitsInCommon++;
	if (hitsInCommon >= nTrackHits2) return kFALSE;
	break;
      }
    } // for (j=0; 
  } // for (i=0; 
  return kTRUE;
}

  //__________________________________________________________________________
void AliMUONTrackK::Kill(void)
{
/// Kill track candidate
  fgTrackReconstructor->GetRecTracksPtr()->Remove(this);
}

  //__________________________________________________________________________
void AliMUONTrackK::FillMUONTrack(void)
{
/// Compute track parameters at hit positions (as for AliMUONTrack)

  // Set Chi2
  SetFitFMin(fChi2);

  // Set track parameters at vertex
  AliMUONTrackParam trackParam;
  SetTrackParam(&trackParam, fTrackPar, fPosition);
  SetTrackParamAtVertex(&trackParam);

  // Set track parameters at hits
  for (Int_t i = fNmbTrackHits-1; i>=0; i--) {
    if ((*fChi2Smooth)[i] < 0) {
      // Propagate through last chambers
      AliMUONTrackExtrap::ExtrapToZ(&trackParam, ((AliMUONHitForRec*)((*fTrackHits)[i]))->GetZ());
    } else {
      // Take saved info
      SetTrackParam(&trackParam, (TMatrixD*)fParSmooth->UncheckedAt(i), ((AliMUONHitForRec*)((*fTrackHits)[i]))->GetZ());
    }
    AddTrackParamAtHit(&trackParam,(AliMUONHitForRec*)fTrackHits->UncheckedAt(i));
    // Fill array of HitForRec's
    AddHitForRecAtHit((AliMUONHitForRec*)fTrackHits->UncheckedAt(i)); 
  }
}

  //__________________________________________________________________________
void AliMUONTrackK::SetTrackParam(AliMUONTrackParam *trackParam, TMatrixD *par, Double_t z)
{
/// Fill AliMUONTrackParam object

  trackParam->SetBendingCoor((*par)(0,0));
  trackParam->SetNonBendingCoor((*par)(1,0));
  trackParam->SetBendingSlope(TMath::Tan((*par)(2,0)));
  trackParam->SetNonBendingSlope(TMath::Tan((*par)(3,0))/TMath::Cos((*par)(2,0)));
  trackParam->SetInverseBendingMomentum((*par)(4,0)/TMath::Cos((*par)(3,0)));
  trackParam->SetZ(z);
}

  //__________________________________________________________________________
void AliMUONTrackK::Branson(void)
{
/// Propagates track to the vertex thru absorber using Branson correction
/// (makes use of the AliMUONTrackParam class)
 
  //AliMUONTrackParam *trackParam = new AliMUONTrackParam();
  AliMUONTrackParam trackParam = AliMUONTrackParam();
  /*
  trackParam->SetBendingCoor((*fTrackPar)(0,0));
  trackParam->SetNonBendingCoor((*fTrackPar)(1,0));
  trackParam->SetBendingSlope(TMath::Tan((*fTrackPar)(2,0)));
  trackParam->SetNonBendingSlope(TMath::Tan((*fTrackPar)(3,0))/TMath::Cos((*fTrackPar)(2,0)));
  trackParam->SetInverseBendingMomentum((*fTrackPar)(4,0)/TMath::Cos((*fTrackPar)(3,0)));
  trackParam->SetZ(fPosition);
  */
  SetTrackParam(&trackParam, fTrackPar, fPosition);

  AliMUONTrackExtrap::ExtrapToVertex(&trackParam, Double_t(0.), Double_t(0.), Double_t(0.));

  (*fTrackPar)(0,0) = trackParam.GetBendingCoor();
  (*fTrackPar)(1,0) = trackParam.GetNonBendingCoor();
  (*fTrackPar)(2,0) = TMath::ATan(trackParam.GetBendingSlope());
  (*fTrackPar)(3,0) = TMath::ATan(TMath::Cos((*fTrackPar)(2,0))*trackParam.GetNonBendingSlope());
  (*fTrackPar)(4,0) = TMath::Cos((*fTrackPar)(3,0))*trackParam.GetInverseBendingMomentum();
  fPosition = trackParam.GetZ();
  //delete trackParam;
  if (fgDebug > 0) cout << 1/(*fTrackPar)(4,0) << " " << fPosition << " " << (*fTrackPar)(0,0) << endl;

  // Get covariance matrix
  *fCovariance = *fWeight;
  if (fCovariance->Determinant() != 0) {
    Int_t ifail;
    mnvertLocal(&((*fCovariance)(0,0)), fgkSize,fgkSize,fgkSize,ifail);
  } else {
    AliWarning(" Determinant fCovariance=0:");
  }
}

  //__________________________________________________________________________
void AliMUONTrackK::GoToZ(Double_t zEnd)
{
/// Propagates track to given Z

  ParPropagation(zEnd);
  MSThin(1); // multiple scattering in the chamber
  WeightPropagation(zEnd, kFALSE);
  fPosition = fPositionNew;
  *fTrackPar = *fTrackParNew; 
}

  //__________________________________________________________________________
void AliMUONTrackK::GoToVertex(Int_t iflag)
{
/// Version 3.08
/// Propagates track to the vertex
/// All material constants are taken from AliRoot

    static Double_t x01[5] = { 24.282,  // C
  			       24.282,  // C
  			       11.274,  // Concrete
  			        1.758,  // Fe 
  			        1.758}; // Fe (cm)
  // inner part theta < 3 degrees
    static Double_t x02[5] = { 30413,  // Air
			       24.282, // C
			       11.274, // Concrete
			       1.758,  // Fe
			       0.369}; // W (cm)
  // z positions of the materials inside the absober outer part theta > 3 degres
  static Double_t zPos[10] = {-90, -105, -315, -443, -468};

  Double_t dZ, r0Norm, x0, deltaP, dChi2, pTotal, pOld;
  AliMUONHitForRec *hit;
  AliMUONRawCluster *clus;
  TClonesArray *rawclusters;

  // First step to the rear end of the absorber
  Double_t zRear = -503;
  GoToZ(zRear);
  Double_t tan3 = TMath::Tan(3./180*TMath::Pi());

  // Go through absorber
  pOld = 1/(*fTrackPar)(4,0);
  Double_t r0Rear = (*fTrackPar)(0,0)*(*fTrackPar)(0,0) + 
                    (*fTrackPar)(1,0)*(*fTrackPar)(1,0);
  r0Rear = TMath::Sqrt(r0Rear)/TMath::Abs(fPosition)/tan3;
  r0Norm = r0Rear;
  Double_t p0, cos25, cos60;
  if (!iflag) goto vertex;

  for (Int_t i=4; i>=0; i--) {
    ParPropagation(zPos[i]);
    WeightPropagation(zPos[i], kFALSE);
    dZ = TMath::Abs (fPositionNew-fPosition);
    if (r0Norm > 1) x0 = x01[i];
    else x0 = x02[i];
    MSLine(dZ,x0); // multiple scattering in the medium (linear approximation)
    fPosition = fPositionNew;
    *fTrackPar = *fTrackParNew; 
    r0Norm = (*fTrackPar)(0,0)*(*fTrackPar)(0,0) + 
             (*fTrackPar)(1,0)*(*fTrackPar)(1,0);
    r0Norm = TMath::Sqrt(r0Norm)/TMath::Abs(fPosition)/tan3;
  }
  // Correct momentum for energy losses
  pTotal = 1/TMath::Abs((*fTrackPar)(4,0));
  p0 = pTotal;
  cos25 = TMath::Cos(2.5/180*TMath::Pi());
  cos60 = TMath::Cos(6.0/180*TMath::Pi());
  for (Int_t j=0; j<1; j++) {
    /*
    if (r0Rear > 1) {
      if (p0 < 20) {
	deltaP = 2.164 + 0.145e-1*p0 - 0.417e-3*p0*p0;
      } else {
	deltaP = 2.275 + 0.102e-2*p0 - 0.674e-6*p0*p0;
      }
    } else {
      if (p0 < 20) {
	deltaP = 2.581 + 0.188e-1*p0 - 0.398e-3*p0*p0;
      } else {
	deltaP = 2.727 + 0.356e-2*p0 + 0.242e-5*p0*p0;
      }
    }
    */
    /*
    if (r0Rear < 1) {
      //W
      if (p0<15) {
	deltaP = 2.737 + 0.0494*p0 - 0.001123*p0*p0;
      } else {
	deltaP = 3.0643 + 0.01346*p0;
      }
      deltaP *= 0.95;
    } else {
      //Pb
      if (p0<15) {
	deltaP  = 2.1380 + 0.0351*p0 - 0.000853*p0*p0;
      } else {
	deltaP = 2.407 + 0.00702*p0;
      }
      deltaP *= 0.95;
    }
    */
    /*
    if (r0Rear < 1) {
      //W
      if (p0<18) {
	deltaP = 2.439 + 0.806e-1*p0 - 0.500e-2*p0*p0 + 0.106e-3*p0*p0*p0;
      } else {
	deltaP = 2.767 + 0.742e-2*p0 - 0.196e-4*p0*p0 + 0.403e-7*p0*p0*p0;
      }
      //deltaP += 0.2;
      deltaP *= cos25;
    } else {
      //Pb
      if (p0<18) {
	deltaP  = 2.209 + 0.800e-2*p0;
      } else {
	deltaP = 2.285 + 0.141e-2*p0 - 0.446e-6*p0*p0;
      }
      //deltaP += 0.2;
      deltaP *= cos60;
    }
    deltaP *= 1.1;
    */
    //*
    if (r0Rear  < 1) {
      if (p0 < 20) {
	deltaP = 2.5938 + 0.0570 * p0 - 0.001151 * p0 * p0;
      } else {
	deltaP = 3.0714 + 0.011767 * p0;
      }
      deltaP *= 0.75; 
    } else {
      if (p0 < 20) {
	deltaP  = 2.1207 + 0.05478 * p0 - 0.00145079 * p0 * p0;
      } else { 
	deltaP = 2.6069 + 0.0051705 * p0;
      }
      deltaP *= 0.9; 
    }
    //*/

    p0 = pTotal + deltaP/TMath::Abs(TMath::Cos((*fTrackPar)(2,0))/TMath::Cos((*fTrackPar)(3,0)));
  }
  (*fTrackPar)(4,0) = 1/p0*TMath::Sign((Double_t)1.,(*fTrackPar)(4,0));

  // Go to the vertex
vertex:
  ParPropagation((Double_t)0.);
  WeightPropagation((Double_t)0., kFALSE);
  fPosition = fPositionNew;
  //*fTrackPar = *fTrackParNew; 
  // Add vertex as a hit
  TMatrixD pointWeight(fgkSize,fgkSize);
  TMatrixD point(fgkSize,1);
  TMatrixD trackParTmp = point;
  point(0,0) = 0; // vertex coordinate - should be taken somewhere
  point(1,0) = 0; // vertex coordinate - should be taken somewhere
  pointWeight(0,0) = 1/1.e-3/1.e-3; // 10 um error
  pointWeight(1,1) = 1/1.e-3/1.e-3; // 10 um error
  TryPoint(point,pointWeight,trackParTmp,dChi2);
  *fTrackPar = trackParTmp;
  *fWeight += pointWeight; 
  fChi2 += dChi2; // Chi2
  if (fgDebug < 0) return; // no output

  cout << pOld << " " << 1/(*fTrackPar)(4,0) << " " << dChi2 << " " << fChi2 << " " << fNmbTrackHits << endl;
  for (Int_t i1=0; i1<fNmbTrackHits; i1++) {
    hit =  (AliMUONHitForRec*) ((*fTrackHits)[i1]);
    printf ("%5d", hit->GetChamberNumber()); 
  }
  cout << endl;
  if (fgDebug > 0) {
    for (Int_t i1=0; i1<fNmbTrackHits; i1++) {
      hit =  (AliMUONHitForRec*) ((*fTrackHits)[i1]);
      //cout << ((AliMUONHitForRec*)((*fTrackHits)[i1]))->GetHitNumber() << " ";
      //cout << ((AliMUONHitForRec*)((*fTrackHits)[i1]))->GetZ() << " ";
      printf ("%5d", fgHitForRec->IndexOf(hit)); 
    }
    cout << endl;
  }

  // from raw clusters
  for (Int_t i1=0; i1<fNmbTrackHits; i1++) {
    hit =  (AliMUONHitForRec*) ((*fTrackHits)[i1]);
    if (hit->GetHitNumber() < 0) { // combined cluster / track finder
      Int_t index = -hit->GetHitNumber() / 100000;
      Int_t iPos = -hit->GetHitNumber() - index * 100000;
      clus = (AliMUONRawCluster*) fgCombi->DetElem(index-1)->RawClusters()->UncheckedAt(iPos);
    } else {
      rawclusters = fgTrackReconstructor->GetMUONData()->RawClusters(hit->GetChamberNumber());
      clus = (AliMUONRawCluster*) rawclusters->UncheckedAt(hit->GetHitNumber());
    }
    printf ("%5d", clus->GetTrack(1)%10000000); 
    
    cout << endl;
    for (Int_t i1=0; i1<fNmbTrackHits; i1++) {
      hit =  (AliMUONHitForRec*) ((*fTrackHits)[i1]);
      if (hit->GetHitNumber() < 0) { // combined cluster / track finder
	Int_t index = -hit->GetHitNumber() / 100000;
	Int_t iPos = -hit->GetHitNumber() - index * 100000;
	clus = (AliMUONRawCluster*) fgCombi->DetElem(index-1)->RawClusters()->UncheckedAt(iPos);
      } else {
	rawclusters = fgTrackReconstructor->GetMUONData()->RawClusters(hit->GetChamberNumber());
	clus = (AliMUONRawCluster*) rawclusters->UncheckedAt(hit->GetHitNumber());
      }
      if (clus->GetTrack(2) != -1) printf ("%5d", clus->GetTrack(2)%10000000);
      else printf ("%5s", "    ");
    }
  }
  cout << endl;
  for (Int_t i1=0; i1<fNmbTrackHits; i1++) {
    //cout << ((AliMUONHitForRec*)((*fTrackHits)[i1]))->GetHitNumber() << " ";
    cout << ((AliMUONHitForRec*)((*fTrackHits)[i1]))->GetZ() << " ";
    //cout << fgHitForRec->IndexOf(((AliMUONHitForRec*)((*fTrackHits)[i1]))) << " ";
  }
  cout << endl;
  for (Int_t i1=0; i1<fNmbTrackHits; i1++) printf("%8.4f", (*fChi2Smooth)[i1]);
  cout << endl;
  cout << "---------------------------------------------------" << endl;

  // Get covariance matrix
  /* Not needed - covariance matrix is not interesting to anybody
  *fCovariance = *fWeight;
  if (fCovariance->Determinant() != 0) {
    Int_t ifail;
    mnvertLocal(&((*fCovariance)(0,0)), fgkSize,fgkSize,fgkSize,ifail);
  } else {
    AliWarning(" Determinant fCovariance=0:" );
  }
  */
}

  //__________________________________________________________________________
void AliMUONTrackK::MSLine(Double_t dZ, Double_t x0)
{
/// Adds multiple scattering in a thick layer for linear propagation

  Double_t cosAlph = TMath::Cos((*fTrackPar)(2,0));
  Double_t tanAlph = TMath::Tan((*fTrackPar)(2,0));
  Double_t cosBeta = TMath::Cos((*fTrackPar)(3,0));
  Double_t sinBeta;
  sinBeta = TMath::Sin((*fTrackPar)(3,0));
  Double_t tanBeta = TMath::Tan((*fTrackPar)(3,0));
  Double_t momentum = 1/(*fTrackPar)(4,0);
  Double_t velo = 1; // relativistic velocity
  Double_t step = TMath::Abs(dZ/cosAlph/cosBeta); // step length

  // Projected scattering angle
  Double_t theta0 = 0.0136/velo/momentum/TMath::Sqrt(x0)*(1+0.038*TMath::Log(step/x0)); 
  Double_t theta02 = theta0*theta0;
  Double_t dl2 = step*step/2*theta02;
  Double_t dl3 = dl2*step*2/3;

  //Derivatives
  Double_t dYdT = 1/cosAlph;
  Double_t dYdB = 0; //(*fTrackPar)(2,0)*sinBeta/cosAlph;
  Double_t dXdT = tanAlph*tanBeta;
  //Double_t dXdB = (1+(*fTrackPar)(2,0)*tanAlph*sinBeta*sinBeta)/cosBeta;
  Double_t dXdB = 1/cosBeta;
  Double_t dAdT = 1/cosBeta;
  Double_t dAdB = 0; //(*fTrackPar)(2,0)*tanBeta;

  // Get covariance matrix
  *fCovariance = *fWeight;
  if (fCovariance->Determinant() != 0) {
    //   fCovariance->Invert();
    Int_t ifail;
    mnvertLocal(&((*fCovariance)(0,0)), fgkSize,fgkSize,fgkSize,ifail);
  } else {
    AliWarning(" Determinant fCovariance=0:" );
  }

  (*fCovariance)(0,0) += dl3*(dYdT*dYdT+dYdB*dYdB); // <yy>
  (*fCovariance)(1,1) += dl3*(dXdT*dXdT+dXdB*dXdB); // <xx>
  (*fCovariance)(2,2) += theta02*step*(dAdT*dAdT+dAdB*dAdB); // <aa>
  (*fCovariance)(3,3) += theta02*step; // <bb>

  (*fCovariance)(0,1) += dl3*(dYdT*dXdT+dYdB*dXdB); // <yx>
  (*fCovariance)(1,0) = (*fCovariance)(0,1);

  (*fCovariance)(0,2) += dl2*(dYdT*dAdT+dYdB*dAdB); // <ya>
  (*fCovariance)(2,0) = (*fCovariance)(0,2);

  (*fCovariance)(0,3) += dl2*dYdB; // <yb>
  (*fCovariance)(3,0) = (*fCovariance)(0,3);

  (*fCovariance)(1,2) += dl2*(-dXdT*dAdT+dXdB*dAdB); // <xa>
  (*fCovariance)(2,1) = (*fCovariance)(1,2);

  (*fCovariance)(1,3) += dl2*dXdB; // <xb>
  (*fCovariance)(3,1) = (*fCovariance)(1,3);

  (*fCovariance)(2,3) += theta02*step*dAdB; // <ab>
  (*fCovariance)(3,2) = (*fCovariance)(2,3);

  // Get weight matrix
  *fWeight = *fCovariance;
  if (fWeight->Determinant() != 0) {
    Int_t ifail;
    mnvertLocal(&((*fWeight)(0,0)), fgkSize,fgkSize,fgkSize,ifail);
  } else {
    AliWarning(" Determinant fWeight=0:");
  }
}
 
  //__________________________________________________________________________
Bool_t AliMUONTrackK::Recover(void)
{
/// Adds new failed track(s) which can be tried to be recovered
  Int_t nRecTracks;
  TClonesArray *trackPtr;
  AliMUONTrackK *trackK;

  if (fgDebug > 0) cout << " ******** Enter Recover " << endl;
  trackPtr = fgTrackReconstructor->GetRecTracksPtr();

  // Remove hit with the highest chi2
  Double_t chi2 = 0;
  if (fgDebug > 0) {
    for (Int_t i=0; i<fNmbTrackHits; i++) {
      chi2 = fChi2Smooth ? (*fChi2Smooth)[i] : (*fChi2Array)[i];
      printf("%10.4f", chi2);
    }
    printf("\n");
    for (Int_t i=0; i<fNmbTrackHits; i++) {
      printf("%10d", ((AliMUONHitForRec*)fTrackHits->UncheckedAt(i))->GetChamberNumber());
    }
    printf("\n");
  }
  Double_t chi2max = 0;
  Int_t imax = 0;
  for (Int_t i=0; i<fNmbTrackHits; i++) {
    chi2 = fChi2Smooth ? (*fChi2Smooth)[i] : (*fChi2Array)[i];
    if (chi2 < chi2max) continue;
    chi2max = chi2;
    imax = i;
  }
  //if (chi2max < 10) return kFALSE; // !!!
  //if (chi2max < 25) imax = fNmbTrackHits - 1;
  if (chi2max < 15) imax = fNmbTrackHits - 1; // discard the last point
  // Check if the outlier is not from the seed segment
  AliMUONHitForRec *skipHit = (AliMUONHitForRec*) fTrackHits->UncheckedAt(imax);
  if (skipHit == fStartSegment->GetHitForRec1() || skipHit == fStartSegment->GetHitForRec2()) {
    //DropBranches(fStartSegment); // drop all tracks with the same seed segment
    return kFALSE; // to be changed probably
  }
  
  // Make a copy of track hit collection
  TObjArray *hits = new TObjArray(*fTrackHits);
  Int_t imax0;
  imax0 = imax;

  // Hits after the found one will be removed
  if (GetStation0() == 3 && skipHit->GetChamberNumber() >= 7) {
    SortHits(1, fTrackHits); // unsort hits
    imax = fTrackHits->IndexOf(skipHit);
  }
  Int_t nTrackHits = fNmbTrackHits;
  for (Int_t i=imax+1; i<nTrackHits; i++) {
    AliMUONHitForRec *hit = (AliMUONHitForRec*) fTrackHits->UncheckedAt(i);
    fTrackHits->Remove(hit);
    hit->SetNTrackHits(hit->GetNTrackHits()-1); // unmark hit 
    fNmbTrackHits--;
  }

  // Check if the track candidate doesn't exist yet
  if (ExistDouble()) { delete hits; return kFALSE; }

  //DropBranches(imax0, hits); // drop branches downstream the discarded hit
  delete hits;

  nRecTracks = fgTrackReconstructor->GetNRecTracks();
  skipHit = (AliMUONHitForRec*) ((*fTrackHits)[fNmbTrackHits-1]);
  // Remove all saved steps and smoother matrices after the skipped hit 
  RemoveMatrices(skipHit->GetZ());

  //AZ(z->-z) if (skipHit->GetZ() > fStartSegment->GetHitForRec2()->GetZ() || !fNSteps) {
  if (TMath::Abs(skipHit->GetZ()) > TMath::Abs(fStartSegment->GetHitForRec2()->GetZ()) || !fNSteps) {
    // Propagation toward high Z or skipped hit next to segment - 
    // start track from segment 
    trackK = new ((*trackPtr)[nRecTracks]) AliMUONTrackK(fStartSegment); 
    fgTrackReconstructor->SetNRecTracks(nRecTracks+1);
    trackK->fRecover = 1;
    trackK->fSkipHit = skipHit;
    trackK->fNmbTrackHits = fNmbTrackHits;
    delete trackK->fTrackHits; // not efficient ?
    trackK->fTrackHits = new TObjArray(*fTrackHits);
    if (fgDebug > 0) cout << nRecTracks << " " << trackK->fRecover << endl;
    return kTRUE;
  } 

  trackK = new ((*trackPtr)[nRecTracks]) AliMUONTrackK(NULL, NULL);
  *trackK = *this;
  fgTrackReconstructor->SetNRecTracks(nRecTracks+1);
  //AZ(z->-z) trackK->fTrackDir = -1;
  trackK->fTrackDir = 1;
  trackK->fRecover = 1;
  trackK->fSkipHit = skipHit;
  Int_t iD = trackK->fParFilter->Last()->GetUniqueID();
  if (iD > 1) { 
    trackK->fParFilter->Last()->SetUniqueID(iD-1);
    CreateMatrix(trackK->fParFilter); 
  }
  *((TMatrixD*)trackK->fParFilter->Last()) = *((TMatrixD*)fParExtrap->Last());
  trackK->fParFilter->Last()->SetUniqueID(1);
  *(trackK->fTrackPar) = *((TMatrixD*)trackK->fParFilter->Last());
  iD = trackK->fCovFilter->Last()->GetUniqueID();
  if (iD > 1) { 
    trackK->fCovFilter->Last()->SetUniqueID(iD-1);
    CreateMatrix(trackK->fCovFilter); 
  }
  *((TMatrixD*)trackK->fCovFilter->Last()) = *((TMatrixD*)fCovExtrap->Last());
  trackK->fCovFilter->Last()->SetUniqueID(1);
  *(trackK->fWeight) = *((TMatrixD*)trackK->fCovFilter->Last());
  if (trackK->fWeight->Determinant() != 0) {
    Int_t ifail;
    mnvertLocal(&((*(trackK->fWeight))(0,0)), fgkSize,fgkSize,fgkSize,ifail);
  } else {
    AliWarning(" Determinant fWeight=0:");
  }
  trackK->fPosition = trackK->fPositionNew = (*fSteps)[fNSteps-1];
  trackK->fChi2 = 0;
  for (Int_t i=0; i<fNmbTrackHits-1; i++) trackK->fChi2 += (*fChi2Array)[i];
  if (fgDebug > 0) cout << nRecTracks << " " << trackK->fRecover << endl;
  return kTRUE;
}

  //__________________________________________________________________________
void AliMUONTrackK::AddMatrices(AliMUONTrackK *trackK, Double_t dChi2, AliMUONHitForRec *hitAdd)
{
/// Adds matrices for the smoother and keep Chi2 for the point
/// Track parameters
  //trackK->fParFilter->Last()->Print();
  Int_t iD = trackK->fParFilter->Last()->GetUniqueID();
  if (iD > 1) { 
    trackK->fParFilter->Last()->SetUniqueID(iD-1);
    CreateMatrix(trackK->fParFilter); 
    iD = 1; 
  }
  *((TMatrixD*)(trackK->fParFilter->Last())) = *(trackK->fTrackPar);
  trackK->fParFilter->Last()->SetUniqueID(iD);
  if (fgDebug > 1) {
    cout << " Add matrices" << " " << fPosition << " " << dChi2 << " " << fgHitForRec->IndexOf(hitAdd) << endl;
    //trackK->fTrackPar->Print();
    //trackK->fTrackParNew->Print();
    trackK->fParFilter->Last()->Print();
    cout << " Add matrices" << endl;
  }
  // Covariance
  *(trackK->fCovariance) = *(trackK->fWeight);
  if (trackK->fCovariance->Determinant() != 0) {
    Int_t ifail;
    mnvertLocal(&((*(trackK->fCovariance))(0,0)), fgkSize,fgkSize,fgkSize,ifail);
  } else {
    AliWarning(" Determinant fCovariance=0:");
  }
  iD = trackK->fCovFilter->Last()->GetUniqueID();
  if (iD > 1) { 
    trackK->fCovFilter->Last()->SetUniqueID(iD-1);
    CreateMatrix(trackK->fCovFilter); 
    iD = 1; 
  }
  *((TMatrixD*)(trackK->fCovFilter->Last())) = *(trackK->fCovariance);
  trackK->fCovFilter->Last()->SetUniqueID(iD);

  // Save Chi2-increment for point
  Int_t indx = trackK->fTrackHits->IndexOf(hitAdd);
  if (indx < 0) indx = fNmbTrackHits;
  if (trackK->fChi2Array->fN <= indx) trackK->fChi2Array->Set(indx+10); 
  trackK->fChi2Array->AddAt(dChi2,indx);
}

  //__________________________________________________________________________
void AliMUONTrackK::CreateMatrix(TObjArray *objArray) const
{
/// Create new matrix and add it to TObjArray 

  TMatrixD *matrix = (TMatrixD*) objArray->First();
  TMatrixD *tmp = new TMatrixD(*matrix);
  objArray->AddAtAndExpand(tmp,objArray->GetLast());
  //cout << " !!! " << tmp << " " << tmp->GetUniqueID() << " " << tmp->GetNoElements() << endl;
}

  //__________________________________________________________________________
void AliMUONTrackK::RemoveMatrices(Double_t zEnd)
{
/// Remove matrices (and saved steps) in the smoother part with abs(z) < abs(zEnd)

  for (Int_t i=fNSteps-1; i>=0; i--) {
    if (fgDebug > 1) printf("%10.4f %10.4f \n", (*fSteps)[i], zEnd);
    if (TMath::Abs((*fSteps)[i]) > TMath::Abs(zEnd) - 0.01) break;
    RemoveMatrices(this);
  } // for (Int_t i=fNSteps-1;
}

  //__________________________________________________________________________
void AliMUONTrackK::RemoveMatrices(AliMUONTrackK* trackK)
{
/// Remove last saved matrices and steps in the smoother part 

  trackK->fNSteps--;
  Int_t i = trackK->fNSteps;

  Int_t id = 0;
  // Delete only matrices not shared by several tracks
  id = trackK->fParExtrap->Last()->GetUniqueID();
  if (id > 1) {
    trackK->fParExtrap->Last()->SetUniqueID(id-1);
    trackK->fParExtrap->RemoveAt(i);
  }
  else ((TMatrixD*)(trackK->fParExtrap->RemoveAt(i)))->Delete();
  id = fParFilter->Last()->GetUniqueID();
  if (id > 1) { 
    trackK->fParFilter->Last()->SetUniqueID(id-1);
    trackK->fParFilter->RemoveAt(i);
  }
  else ((TMatrixD*)(trackK->fParFilter->RemoveAt(i)))->Delete();
  id = trackK->fCovExtrap->Last()->GetUniqueID();
  if (id > 1) { 
    trackK->fCovExtrap->Last()->SetUniqueID(id-1);
    trackK->fCovExtrap->RemoveAt(i);
  }
  else ((TMatrixD*)(trackK->fCovExtrap->RemoveAt(i)))->Delete();
  id = trackK->fCovFilter->Last()->GetUniqueID();
  if (id > 1) { 
    trackK->fCovFilter->Last()->SetUniqueID(id-1);
    trackK->fCovFilter->RemoveAt(i);
  }
  else ((TMatrixD*)(trackK->fCovFilter->RemoveAt(i)))->Delete();
  id = trackK->fJacob->Last()->GetUniqueID();
  if (id > 1) { 
    trackK->fJacob->Last()->SetUniqueID(id-1);
    trackK->fJacob->RemoveAt(i);
  }
  else ((TMatrixD*)(trackK->fJacob->RemoveAt(i)))->Delete();
}

  //__________________________________________________________________________
Bool_t AliMUONTrackK::Smooth(void)
{
/// Apply smoother
  Int_t ihit = fNmbTrackHits - 1;
  AliMUONHitForRec *hit = (AliMUONHitForRec*) (*fTrackHits)[ihit];
  fChi2Smooth = new TArrayD(fNmbTrackHits);
  fChi2Smooth->Reset(-1);
  fChi2 = 0;
  fParSmooth = new TObjArray(15);
  fParSmooth->Clear();

  if (fgDebug > 0) {
    cout << " ******** Enter Smooth " << endl;
    cout << (*fSteps)[fNSteps-1] << " " << fPosition << " " << fPositionNew << endl; 
    /*
    for (Int_t i=fNSteps-1; i>=0; i--) cout << (*fSteps)[i] << " ";
    cout << endl;
    for (Int_t i=fNSteps-1; i>=0; i--) {cout << i << " " << (*fSteps)[i]; ((TMatrixD*)fParFilter->UncheckedAt(i))->Print(); ((TMatrixD*)fParExtrap->UncheckedAt(i))->Print(); ((TMatrixD*)fJacob->UncheckedAt(i))->Print();}
    */
    for (Int_t i=ihit; i>=0; i--) cout << ((AliMUONHitForRec*)(*fTrackHits)[i])->GetZ() << " ";
    cout << endl;
  }

  // Find last point corresponding to the last hit
  Int_t iLast = fNSteps - 1;
  for ( ; iLast>=0; iLast--) {
    //AZ(z->-z) if ((*fSteps)[iLast] + 0.01 > ((AliMUONHitForRec*)(*fTrackHits)[fNmbTrackHits-1])->GetZ()) break;
    if (TMath::Abs((*fSteps)[iLast]) + 0.01 > TMath::Abs(((AliMUONHitForRec*)(*fTrackHits)[fNmbTrackHits-1])->GetZ())) break;
  }

  TMatrixD parSmooth = *((TMatrixD*)(fParFilter->UncheckedAt(iLast))); // last filtered = first smoothed
  //parSmooth.Dump();
  TMatrixD covSmooth = *((TMatrixD*)(fCovFilter->UncheckedAt(iLast))); // last filtered = first smoothed
  TMatrixD tmp=covSmooth, weight=covSmooth, gain = covSmooth;
  TMatrixD tmpPar = *fTrackPar;
  //parSmooth.Print(); ((TMatrixD*)(fParExtrap->Last()))->Print(); 

  Bool_t found;
  Double_t chi2max = 0;
  for (Int_t i=iLast+1; i>0; i--) {
    if (i == iLast + 1) goto L33; // temporary fix

    // Smoother gain matrix
    weight = *((TMatrixD*)(fCovExtrap->UncheckedAt(i)));
    if (weight.Determinant() != 0) {
      Int_t ifail;
      mnvertLocal(&(weight(0,0)), fgkSize,fgkSize,fgkSize,ifail);
    } else {
      AliWarning(" Determinant weight=0:");
    }
    // Fk'Wkk+1
    tmp = TMatrixD(*((TMatrixD*)(fJacob->UncheckedAt(i))),TMatrixD::kTransposeMult,weight);
    gain = TMatrixD(*((TMatrixD*)(fCovFilter->UncheckedAt(i-1))),TMatrixD::kMult,tmp);
    //cout << (*fSteps)[i-1] << " " << (*fSteps)[i] << endl; gain.Print();

    // Smoothed parameter vector
    //((TMatrixD*)(fParExtrap->UncheckedAt(i)))->Dump();
    parSmooth -= *((TMatrixD*)(fParExtrap->UncheckedAt(i)));
    tmpPar = *((TMatrixD*)(fParFilter->UncheckedAt(i-1)));
    tmpPar += TMatrixD(gain,TMatrixD::kMult,parSmooth);
    parSmooth = tmpPar;

    // Smoothed covariance
    covSmooth -= *((TMatrixD*)(fCovExtrap->UncheckedAt(i)));
    weight = TMatrixD(TMatrixD::kTransposed,gain);
    tmp = TMatrixD(covSmooth,TMatrixD::kMult,weight);
    covSmooth = *((TMatrixD*)(fCovFilter->UncheckedAt(i-1)));
    covSmooth += TMatrixD(gain,TMatrixD::kMult,tmp);

    // Check if there was a measurement at given z
    found = kFALSE;
    for ( ; ihit>=0; ihit--) { 
      hit = (AliMUONHitForRec*) (*fTrackHits)[ihit];
      if (TMath::Abs(hit->GetZ()-(*fSteps)[i-1]) < 0.1) { found = kTRUE; break; }
      //AZ(z->-z) else if (ihit < fNmbTrackHits-1 && hit->GetZ() > (*fSteps)[i-1]) { ihit++; break; }
      else if (ihit < fNmbTrackHits-1 && TMath::Abs(hit->GetZ()) > TMath::Abs((*fSteps)[i-1])) { ihit++; break; }
    }
    if (!found) continue; // no measurement - skip the rest
    else if (fgDebug > 1) cout << i << " " << ihit << " " << hit->GetZ() << endl;
    if (ihit == 0) continue; // the first hit - skip the rest

L33:
    // Smoothed residuals
    tmpPar = 0;
    tmpPar(0,0) = hit->GetBendingCoor() - parSmooth(0,0); // bending coordinate
    tmpPar(1,0) = hit->GetNonBendingCoor() - parSmooth(1,0); // non-bending coordinate
    if (fgDebug > 1) {
      cout << hit->GetBendingCoor() << " " << parSmooth(0,0) << " " << tmpPar(0,0) << endl;
      cout << hit->GetNonBendingCoor() << " " << parSmooth(1,0) << " " << tmpPar(1,0) << endl;
    }
    // Cov. matrix of smoothed residuals
    tmp = 0;
    tmp(0,0) = hit->GetBendingReso2() - covSmooth(0,0);
    tmp(1,1) = hit->GetNonBendingReso2() - covSmooth(1,1);
    tmp(0,1) = tmp(1,0) = -covSmooth(0,1);
    tmp(2,2) = tmp(3,3) = tmp(4,4) = 1;

    // Calculate Chi2 of smoothed residuals
    if (tmp.Determinant() != 0) {
      Int_t ifail;
      mnvertLocal(&(tmp(0,0)), fgkSize,fgkSize,fgkSize,ifail);
    } else {
      AliWarning(" Determinant tmp=0:");
    }
    TMatrixD vector(tmp,TMatrixD::kMult,tmpPar);
    TMatrixD chi2(tmpPar,TMatrixD::kTransposeMult,vector);
    if (fgDebug > 1) chi2.Print();
    (*fChi2Smooth)[ihit] = chi2(0,0);
    if (chi2(0,0) > chi2max) chi2max = chi2(0,0); 
    fChi2 += chi2(0,0);
    if (chi2(0,0) < 0) { 
      //chi2.Print(); 
      AliError(Form(" *** chi2 < 0: %d %d ", i, iLast));
    }
    // Save smoothed parameters
    TMatrixD *par = new TMatrixD(parSmooth);
    fParSmooth->AddAtAndExpand(par, ihit);

  } // for (Int_t i=iLast+1;

  //if (chi2max > 16) { 
  //if (chi2max > 25) { 
  //if (chi2max > 50) { 
  //if (chi2max > 100) { 
  if (chi2max > fgkChi2max) { 
    //if (Recover()) DropBranches(); 
    //Recover();
    Outlier();
    return kFALSE; 
  }
  return kTRUE;
}
 
  //__________________________________________________________________________
void AliMUONTrackK::Outlier()
{
/// Adds new track with removed hit having the highest chi2

  if (fgDebug > 0) {
    cout << " ******** Enter Outlier " << endl;
    for (Int_t i=0; i<fNmbTrackHits; i++) printf("%10.4f", (*fChi2Smooth)[i]);
    printf("\n");
    for (Int_t i=0; i<fNmbTrackHits; i++) {
      printf("%10d", ((AliMUONHitForRec*)fTrackHits->UncheckedAt(i))->GetChamberNumber());
    }
    printf("\n");
  }

  Double_t chi2max = 0;
  Int_t imax = 0;
  for (Int_t i=0; i<fNmbTrackHits; i++) {
    if ((*fChi2Smooth)[i] < chi2max) continue;
    chi2max = (*fChi2Smooth)[i];
    imax = i;
  }
  // Check if the outlier is not from the seed segment
  AliMUONHitForRec *hit = (AliMUONHitForRec*) fTrackHits->UncheckedAt(imax);
  if (hit == fStartSegment->GetHitForRec1() || hit == fStartSegment->GetHitForRec2()) return; // to be changed probably

  // Check for missing station
  Int_t ok = 1;
  if (imax == 0) {
    if (((AliMUONHitForRec*)fTrackHits->UncheckedAt(1))->GetChamberNumber() < 8) ok--; 
  } else if (imax == fNmbTrackHits-1) {
    if (((AliMUONHitForRec*)fTrackHits->UncheckedAt(fNmbTrackHits-2))->GetChamberNumber() > 1) ok--; 
  } 
  else if (((AliMUONHitForRec*)fTrackHits->UncheckedAt(imax-1))->GetChamberNumber()/2 - ((AliMUONHitForRec*)fTrackHits->UncheckedAt(imax+1))->GetChamberNumber()/2 > 1) ok--;
  if (!ok) { Recover(); return; } // try to recover track
  //AZ if (!ok) { if (fgDebug >= 0) cout << imax << endl; DropBranches(imax, 0); return; } 

  // Remove saved steps and smoother matrices after the outlier
  RemoveMatrices(hit->GetZ()); 
  
  // Check for possible double track candidates
  //if (ExistDouble(hit)) return;

  TClonesArray *trackPtr = fgTrackReconstructor->GetRecTracksPtr();
  Int_t nRecTracks = fgTrackReconstructor->GetNRecTracks();

  AliMUONTrackK *trackK = 0;
  if (!fNSteps || GetStation0() == 3 && hit->GetChamberNumber() > 7) {
    // start track from segment 
    trackK = new ((*trackPtr)[nRecTracks]) AliMUONTrackK(fStartSegment); 
    fgTrackReconstructor->SetNRecTracks(nRecTracks+1);
    trackK->fRecover = 2;
    trackK->fSkipHit = hit;
    trackK->fNmbTrackHits = fNmbTrackHits;

    hit = (AliMUONHitForRec*) trackK->fTrackHits->UncheckedAt(0);
    hit->SetNTrackHits(hit->GetNTrackHits()-1);
    hit = (AliMUONHitForRec*) trackK->fTrackHits->UncheckedAt(1);
    hit->SetNTrackHits(hit->GetNTrackHits()-1);
    delete trackK->fTrackHits; // not efficient ?
    trackK->fTrackHits = new TObjArray(*fTrackHits);
    for (Int_t i = 0; i < fNmbTrackHits; i++) {
      hit = (AliMUONHitForRec*) fTrackHits->UncheckedAt(i);
      hit->SetNTrackHits(hit->GetNTrackHits()+1);
    }

    if (GetStation0() == 3) trackK->SortHits(1, trackK->fTrackHits);
    if (fgDebug > 0) cout << nRecTracks << " Outlier" << endl;
    return;
  } 
  trackK = new ((*trackPtr)[nRecTracks]) AliMUONTrackK(NULL, NULL);
  *trackK = *this;
  fgTrackReconstructor->SetNRecTracks(nRecTracks+1);
  trackK->fTrackDir = 1;
  trackK->fRecover = 2;
  trackK->fSkipHit = hit;
  Int_t iD = trackK->fParFilter->Last()->GetUniqueID();
  if (iD > 1) { 
    trackK->fParFilter->Last()->SetUniqueID(iD-1);
    CreateMatrix(trackK->fParFilter); 
  }
  *((TMatrixD*)trackK->fParFilter->Last()) = *((TMatrixD*)fParExtrap->Last());
  trackK->fParFilter->Last()->SetUniqueID(1);
  *(trackK->fTrackPar) = *((TMatrixD*)trackK->fParFilter->Last());
  iD = trackK->fCovFilter->Last()->GetUniqueID();
  if (iD > 1) { 
    trackK->fCovFilter->Last()->SetUniqueID(iD-1);
    CreateMatrix(trackK->fCovFilter); 
  }
  *((TMatrixD*)trackK->fCovFilter->Last()) = *((TMatrixD*)fCovExtrap->Last());
  trackK->fCovFilter->Last()->SetUniqueID(1);
  *(trackK->fWeight) = *((TMatrixD*)trackK->fCovFilter->Last());
  if (trackK->fWeight->Determinant() != 0) {
    Int_t ifail;
    mnvertLocal(&((*(trackK->fWeight))(0,0)), fgkSize,fgkSize,fgkSize,ifail);
  } else {
    AliWarning(" Determinant fWeight=0:");
  }
  trackK->fPosition = trackK->fPositionNew = (*fSteps)[fNSteps-1];
  trackK->fChi2 = 0;
  for (Int_t i=0; i<fNmbTrackHits-1; i++) trackK->fChi2 += (*fChi2Array)[i];
  if (fgDebug > 0) cout << nRecTracks << " Outlier " << trackK->fChi2 << endl;
}

  //__________________________________________________________________________
Double_t AliMUONTrackK::GetChi2PerPoint(Int_t iPoint) const
{
/// Return Chi2 at point
  return fChi2Smooth ? (*fChi2Smooth)[iPoint] : (*fChi2Array)[iPoint];
  //return 0.;
}

  //__________________________________________________________________________
void AliMUONTrackK::Print(FILE *lun) const
{
/// Print out track information

  Int_t flag = 1;
  AliMUONHitForRec *hit = 0; 
    // from raw clusters
  AliMUONRawCluster *clus = 0;
  TClonesArray *rawclusters = 0;
  for (Int_t i1=0; i1<fNmbTrackHits; i1++) {
    hit =  (AliMUONHitForRec*) ((*fTrackHits)[i1]);
    rawclusters = fgTrackReconstructor->GetMUONData()->RawClusters(hit->GetChamberNumber());
    clus = (AliMUONRawCluster*) rawclusters->UncheckedAt(hit->GetHitNumber());
    if (TMath::Abs(clus->GetTrack(1)-1) < 2) {
      if (clus->GetTrack(2)) flag = 2;
      continue;
    }
    if (clus->GetTrack(2) && TMath::Abs(clus->GetTrack(2)-1) < 2) {
      flag = 3;
      continue;
    }
    flag = 0;
    break;
    
    Int_t sig[2]={1,1}, tid[2]={0};
    for (Int_t i1=0; i1<fNmbTrackHits; i1++) {
      if (GetChi2PerPoint(i1) < -0.1) continue;
      hit =  (AliMUONHitForRec*) ((*fTrackHits)[i1]);
      rawclusters = fgTrackReconstructor->GetMUONData()->RawClusters(hit->GetChamberNumber());
      clus = (AliMUONRawCluster*) rawclusters->UncheckedAt(hit->GetHitNumber());
      for (Int_t j=0; j<2; j++) {
	tid[j] = clus->GetTrack(j+1) - 1;
	if (clus->GetTrack(j+1) < 0) { sig[j] = 0; tid[j] = 999; }
      }
      fprintf(lun,"%3d %3d %10.4f", gAlice->GetEvNumber(), hit->GetChamberNumber(), GetChi2PerPoint(i1));
      if (!(clus->GetTrack(2))) fprintf(lun, "%3d %3d", sig[0], tid[0]); // simple cluster
      else { // track overlap
	fprintf(lun, "%3d %3d", TMath::Max(sig[0],sig[1]), TMath::Min(tid[0],tid[1]));
	//if (tid[0] < 2) flag *= 2;
	//else if (tid[1] < 2) flag *= 3;
      }
      fprintf (lun, "%3d \n", flag); 
    }
  }
}

  //__________________________________________________________________________
void AliMUONTrackK::DropBranches(Int_t imax, TObjArray *hits)
{
/// Drop branches downstream of the skipped hit 
  Int_t nRecTracks;
  TClonesArray *trackPtr;
  AliMUONTrackK *trackK;

  trackPtr = fgTrackReconstructor->GetRecTracksPtr();
  nRecTracks = fgTrackReconstructor->GetNRecTracks();
  Int_t icand = trackPtr->IndexOf(this);
  if (!hits) hits = fTrackHits; 

  // Check if the track candidate doesn't exist yet
  for (Int_t i=icand+1; i<nRecTracks; i++) {
    trackK = (AliMUONTrackK*) ((*trackPtr)[i]);
    if (trackK->fNmbTrackHits == 2 && trackK->GetRecover() == 0) continue;
    if (trackK->GetRecover() < 0) continue;

    if (trackK->fNmbTrackHits >= imax + 1) {
      for (Int_t j=0; j<=imax; j++) {
	//if (j != fNmbTrackHits-1 && (*trackK->fTrackHits)[j] != (*fTrackHits)[j]) break;
	if ((*trackK->fTrackHits)[j] != (*hits)[j]) break;
	if (j == imax) {
	  if (hits != fTrackHits) {
	    // Drop all branches downstream the hit (for Recover)
	    trackK->SetRecover(-1);
	    if (fgDebug >= 0 )cout << " Recover " << i << endl;
	    continue;
	  }
	  // Check if the removal of the hit breaks the track
	  Int_t ok = 1;
	  if (imax == 0) {
	    if (((AliMUONHitForRec*)trackK->fTrackHits->UncheckedAt(1))->GetChamberNumber() < 8) ok--; }
	  else if (imax == trackK->fNmbTrackHits-1) continue;
	    // else if (imax == trackK->fNmbTrackHits-1) {
	    //if (((AliMUONHitForRec*)trackK->fTrackHits->UncheckedAt(trackK->fNmbTrackHits-2))->GetChamberNumber() > 1) ok--; 
	    //} 
	  else if (((AliMUONHitForRec*)trackK->fTrackHits->UncheckedAt(imax-1))->GetChamberNumber()/2 - ((AliMUONHitForRec*)trackK->fTrackHits->UncheckedAt(imax+1))->GetChamberNumber()/2 > 1) ok--;
	  if (!ok) trackK->SetRecover(-1);
	}
      } // for (Int_t j=0;
    }
  } // for (Int_t i=0;
}

  //__________________________________________________________________________
void AliMUONTrackK::DropBranches(AliMUONSegment *segment)
{
/// Drop all candidates with the same seed segment
  Int_t nRecTracks;
  TClonesArray *trackPtr;
  AliMUONTrackK *trackK;

  trackPtr = fgTrackReconstructor->GetRecTracksPtr();
  nRecTracks = fgTrackReconstructor->GetNRecTracks();
  Int_t icand = trackPtr->IndexOf(this);

  for (Int_t i=icand+1; i<nRecTracks; i++) {
    trackK = (AliMUONTrackK*) ((*trackPtr)[i]);
    if (trackK->fNmbTrackHits == 2 && trackK->GetRecover() == 0) continue;
    if (trackK->GetRecover() < 0) continue;
    if (trackK->fStartSegment == segment) trackK->SetRecover(-1);
  }
  if (fgDebug >= 0) cout << " Drop segment " << endl; 
}

  //__________________________________________________________________________
AliMUONHitForRec* AliMUONTrackK::GetHitLastOk(void)
{
/// Return the hit where track stopped

  if (!fNSteps) return (AliMUONHitForRec*)((*fTrackHits)[1]);
  return fSkipHit;
}

  //__________________________________________________________________________
Int_t AliMUONTrackK::GetStation0(void)
{
/// Return seed station number
  return fStartSegment->GetHitForRec1()->GetChamberNumber() / 2;
}

  //__________________________________________________________________________
Bool_t AliMUONTrackK::ExistDouble(AliMUONHitForRec *hit)
{
/// Check if the track will make a double after outlier removal

  TClonesArray *trackPtr = fgTrackReconstructor->GetRecTracksPtr();
  Int_t nRecTracks = fgTrackReconstructor->GetNRecTracks();
  TObjArray *hitArray = new TObjArray(*fTrackHits);
  TObjArray *hitArray1 = new TObjArray(*hitArray);
  hitArray1->Remove(hit);
  hitArray1->Compress();

  Bool_t same = kFALSE;
  for (Int_t i=0; i<nRecTracks; i++) {
    AliMUONTrackK *trackK = (AliMUONTrackK*) ((*trackPtr)[i]);
    if (trackK->fNmbTrackHits == 2 && trackK->GetRecover() == 0) continue;
    if (trackK == this) continue;
    if (trackK->fNmbTrackHits == fNmbTrackHits || trackK->fNmbTrackHits == fNmbTrackHits-1) {
      TObjArray *hits = new TObjArray(*trackK->fTrackHits);
      same = kTRUE;
      if (trackK->fNmbTrackHits == fNmbTrackHits) {
	for (Int_t j=0; j<fNmbTrackHits; j++) {
	  if (hits->UncheckedAt(j) != hitArray->UncheckedAt(j)) { same = kFALSE; break; }
	}
	if (same) { delete hits; break; }
	if (trackK->fSkipHit) {
	  TObjArray *hits1 = new TObjArray(*hits);
	  if (hits1->Remove(trackK->fSkipHit) > 0) {
	    hits1->Compress();
	    same = kTRUE;
	    for (Int_t j=0; j<fNmbTrackHits-1; j++) {
	      if (hits1->UncheckedAt(j) != hitArray1->UncheckedAt(j)) { same = kFALSE; break; }
	    }
	    if (same) { delete hits1; break; }
	  }
	  delete hits1;
	}
      } else {
	// Check with removed outlier
	same = kTRUE;
	for (Int_t j=0; j<fNmbTrackHits-1; j++) {
	  if (hits->UncheckedAt(j) != hitArray1->UncheckedAt(j)) { same = kFALSE; break; }
	}
	if (same) { delete hits; break; }
      } 
      delete hits;
    }
  } // for (Int_t i=0; i<nRecTracks;
  delete hitArray; delete hitArray1;
  if (same && fgDebug >= 0) cout << " Same" << endl;
  return same;
}

  //__________________________________________________________________________
Bool_t AliMUONTrackK::ExistDouble(void)
{
/// Check if the track will make a double after recovery

  TClonesArray *trackPtr = fgTrackReconstructor->GetRecTracksPtr();
  Int_t nRecTracks = fgTrackReconstructor->GetNRecTracks();

  TObjArray *hitArray = new TObjArray(*fTrackHits);
  if (GetStation0() == 3) SortHits(0, hitArray); // sort
  //if (GetStation0() == 3) SortHits(1, hitArray); // unsort

  Bool_t same = kFALSE;
  for (Int_t i=0; i<nRecTracks; i++) {
    AliMUONTrackK *trackK = (AliMUONTrackK*) ((*trackPtr)[i]);
    if (trackK->fNmbTrackHits == 2 && trackK->GetRecover() == 0) continue;
    if (trackK == this) continue;
    //AZ if (trackK->GetRecover() < 0) continue; //
    if (trackK->fNmbTrackHits >= fNmbTrackHits) {
      TObjArray *hits = new TObjArray(*trackK->fTrackHits);
      if (trackK->GetStation0() == 3) SortHits(0, hits); // sort
      //if (trackK->GetStation0() == 3) SortHits(1, hits); // unsort
      for (Int_t j=0; j<fNmbTrackHits; j++) {
	//cout << fNmbTrackHits << " " << i << " " << j << " " << (*hitArray)[j] << " " << (*hits)[j] << " " << trackK->fSkipHit << endl;
	if (j != fNmbTrackHits-1 && (*hitArray)[j] != (*hits)[j]) break; 
	if (j == fNmbTrackHits-1) {
	  if (trackK->fSkipHit && TMath::Abs(((AliMUONHitForRec*)((*hitArray)[j]))->GetZ()-trackK->fSkipHit->GetZ()) < 0.5) same = kTRUE; 
	  //if (trackK->fNmbTrackHits > fNmbTrackHits && 
	  //if (trackK->fSkipHit) cout << ((AliMUONHitForRec*)((*hitArray)[j]))->GetZ() << " " << trackK->fSkipHit->GetZ() << endl;
	}
      } // for (Int_t j=0;
      delete hits;
      if (same) break; 
    }
  } // for (Int_t i=0;
  delete hitArray;
  return same;
}
