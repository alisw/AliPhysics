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

#include <stdlib.h> // for exit()

#include <Riostream.h>
#include <TClonesArray.h>
#include <TMatrixD.h>

#include "AliMUONTrackK.h"
#include "AliCallf77.h"
#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONEventReconstructor.h"
#include "AliMUONSegment.h"
#include "AliMUONHitForRec.h"
#include "AliMUONRawCluster.h"
#include "AliMUONTrackParam.h"
#include "AliRun.h"
#include "AliLog.h"
//#include "AliMagF.h"

const Int_t AliMUONTrackK::fgkSize = 5;
const Int_t AliMUONTrackK::fgkNSigma = 4; 
const Int_t AliMUONTrackK::fgkTriesMax = 10000; 
const Double_t AliMUONTrackK::fgkEpsilon = 0.002; 

void mnvertLocalK(Double_t* a, Int_t l, Int_t m, Int_t n, Int_t& ifail);

ClassImp(AliMUONTrackK) // Class implementation in ROOT context

  // A few calls in Fortran or from Fortran (extrap.F).
#ifndef WIN32 
# define extrap_onestep_helix extrap_onestep_helix_
# define extrap_onestep_helix3 extrap_onestep_helix3_
# define extrap_onestep_rungekutta extrap_onestep_rungekutta_
# define gufld_double gufld_double_
#else 
# define extrap_onestep_helix EXTRAP_ONESTEP_HELIX
# define extrap_onestep_helix3 EXTRAP_ONESTEP_HELIX3
# define extrap_onestep_rungekutta EXTRAP_ONESTEP_RUNGEKUTTA
# define gufld_double GUFLD_DOUBLE
#endif 

extern "C" {
  void type_of_call extrap_onestep_helix
  (Double_t &Charge, Double_t &StepLength, Double_t *VGeant3, Double_t *VGeant3New);

  void type_of_call extrap_onestep_helix3
  (Double_t &Field, Double_t &StepLength, Double_t *VGeant3, Double_t *VGeant3New);

  void type_of_call extrap_onestep_rungekutta
  (Double_t &Charge, Double_t &StepLength, Double_t *VGeant3, Double_t *VGeant3New);

  void type_of_call gufld_double(Double_t *Position, Double_t *Field);
    /*  void type_of_call gufld_double(Double_t *Position, Double_t *Field) {
    // interface to "gAlice->Field()->Field" for arguments in double precision
    Float_t x[3], b[3];
    x[0] = Position[0]; x[1] = Position[1]; x[2] = Position[2];
    gAlice->Field()->Field(x, b);
    Field[0] = b[0]; Field[1] = b[1]; Field[2] = b[2];
  }
    */
}

Int_t AliMUONTrackK::fgNOfPoints = 0; 
AliMUON* AliMUONTrackK::fgMUON = NULL;
AliMUONEventReconstructor* AliMUONTrackK::fgEventReconstructor = NULL; 
TClonesArray* AliMUONTrackK::fgHitForRec = NULL; 

  //__________________________________________________________________________
AliMUONTrackK::AliMUONTrackK()
  : TObject()
{
  // Default constructor

  fgEventReconstructor = NULL; // pointer to event reconstructor
  fgMUON = NULL; // pointer to Muon module
  fgHitForRec = NULL; // pointer to points
  fgNOfPoints = 0; // number of points

  fStartSegment = NULL;
  fTrackHitsPtr = NULL;
  fNTrackHits = 0;
  fTrackPar = NULL;
  fTrackParNew = NULL;
  fCovariance = NULL;
  fWeight = NULL;
  fSkipHit = NULL;

  return;
}

  //__________________________________________________________________________
AliMUONTrackK::AliMUONTrackK(AliMUONEventReconstructor *EventReconstructor, TClonesArray *hitForRec)
  : TObject()
{
  // Constructor

  fgEventReconstructor = EventReconstructor; // pointer to event reconstructor
  fgMUON = (AliMUON*) gAlice->GetModule("MUON"); // pointer to Muon module
  fgHitForRec = hitForRec; // pointer to points
  fgNOfPoints = fgHitForRec->GetEntriesFast(); // number of points

  fStartSegment = NULL;
  fTrackHitsPtr = NULL;
  fNTrackHits = 0;
  fChi2 = 0;
  fTrackPar = NULL;
  fTrackParNew = NULL;
  fCovariance = NULL;
  fWeight = NULL;
  fSkipHit = NULL;

  return;
}

  //__________________________________________________________________________
AliMUONTrackK::AliMUONTrackK(AliMUONSegment *segment)
  : TObject()
{
  // Constructor from a segment
  Double_t dX, dY, dZ;
  AliMUONHitForRec *hit1, *hit2;
  AliMUONRawCluster *clus;
  TClonesArray *rawclusters;

  fStartSegment = segment;
  fRecover = 0;
  // Pointers to hits from the segment
  hit1 = segment->GetHitForRec1();
  hit2 = segment->GetHitForRec2();
  hit1->SetNTrackHits(hit1->GetNTrackHits()+1); // mark hit as being on track
  hit2->SetNTrackHits(hit2->GetNTrackHits()+1); // mark hit as being on track
  // check sorting in Z
  if (hit1->GetZ() > hit2->GetZ()) {
    hit1 = hit2;
    hit2 = segment->GetHitForRec1();
  }
  // memory allocation for the TObjArray of pointers to reconstructed TrackHit's
  fTrackHitsPtr = new TObjArray(10);
  fNTrackHits = 2;
  fChi2 = 0;
  fBPFlag = kFALSE;
  fTrackPar = new TMatrixD(fgkSize,1); // track parameters
  fTrackParNew = new TMatrixD(fgkSize,1); // track parameters
  fCovariance = new TMatrixD(fgkSize,fgkSize); // covariance matrix
  fWeight = new TMatrixD(fgkSize,fgkSize); // weight matrix (inverse of covariance)

  // Fill array of track parameters
  if (hit1->GetChamberNumber() > 7) {
    // last tracking station
    (*fTrackPar)(0,0) = hit1->GetBendingCoor(); // y
    (*fTrackPar)(1,0) = hit1->GetNonBendingCoor(); // x
    fPosition = hit1->GetZ(); // z
    fTrackHitsPtr->Add((TObjArray*)hit2); // add hit 2
    fTrackHitsPtr->Add((TObjArray*)hit1); // add hit 1
    fTrackDir = -1;
  } else {
    // last but one tracking station
    (*fTrackPar)(0,0) = hit2->GetBendingCoor(); // y
    (*fTrackPar)(1,0) = hit2->GetNonBendingCoor(); // x
    fPosition = hit2->GetZ(); // z
    fTrackHitsPtr->Add((TObjArray*)hit1); // add hit 1
    fTrackHitsPtr->Add((TObjArray*)hit2); // add hit 2
    fTrackDir = 1;
  }
  dZ = hit2->GetZ() - hit1->GetZ();
  dY = hit2->GetBendingCoor() - hit1->GetBendingCoor();
  dX = hit2->GetNonBendingCoor() - hit1->GetNonBendingCoor();
  (*fTrackPar)(2,0) = TMath::ATan2(dY,dZ); // alpha
  (*fTrackPar)(3,0) = TMath::ATan2(dX,dZ/TMath::Cos((*fTrackPar)(2,0))); // beta
  (*fTrackPar)(4,0) = 1/fgEventReconstructor->GetBendingMomentumFromImpactParam(segment->GetBendingImpact()); // 1/Pt
  (*fTrackPar)(4,0) *= TMath::Cos((*fTrackPar)(3,0)); // 1/p
  cout << fgEventReconstructor->GetBendingMomentumFromImpactParam(segment->GetBendingImpact()) << " " << 1/(*fTrackPar)(4,0) << " ";
  if (fgEventReconstructor->GetRecGeantHits()) { 
    // from GEANT hits
    cout << ((AliMUONHitForRec*)((*fTrackHitsPtr)[0]))->GetTHTrack() << "<-->" << ((AliMUONHitForRec*)((*fTrackHitsPtr)[1]))->GetTHTrack() << endl;
  } else {
    // from raw clusters
    for (Int_t i=0; i<2; i++) {
      hit1 = (AliMUONHitForRec*) ((*fTrackHitsPtr)[i]);
      rawclusters = fgMUON->GetMUONData()->RawClusters(hit1->GetChamberNumber());
      clus = (AliMUONRawCluster*) rawclusters->UncheckedAt(hit1->GetHitNumber());
      cout << clus->GetTrack(1)-1;
      if (clus->GetTrack(2) != 0) cout << " " << clus->GetTrack(2)-1;
      if (i == 0) cout << " <--> ";
    }
    cout << endl;
  }
  // Evaluate covariance (and weight) matrix
  EvalCovariance(dZ);

  return;
}

  //__________________________________________________________________________
AliMUONTrackK::~AliMUONTrackK()
{
  // Destructor

  if (fTrackHitsPtr) {
    delete fTrackHitsPtr; // delete the TObjArray of pointers to TrackHit's
    fTrackHitsPtr = NULL;
  }
  delete fTrackPar; delete fTrackParNew; delete fCovariance;
  delete fWeight; 
}

  //__________________________________________________________________________
AliMUONTrackK::AliMUONTrackK (const AliMUONTrackK& source)
  : TObject(source)
{
// Protected copy constructor

  AliFatal("Not implemented.");
}

  //__________________________________________________________________________
AliMUONTrackK & AliMUONTrackK::operator=(const AliMUONTrackK& source)
{
  // Assignment operator
  // Members
  if(&source == this) return *this;

  // base class assignement
  TObject::operator=(source);

  fStartSegment = source.fStartSegment;
  fNTrackHits = source.fNTrackHits;
  fChi2 = source.fChi2;
  fPosition = source.fPosition;
  fPositionNew = source.fPositionNew;
  fTrackDir = source.fTrackDir;
  fBPFlag = source.fBPFlag;
  fRecover = source.fRecover;
  fSkipHit = source.fSkipHit;

  // Pointers
  fTrackHitsPtr = new TObjArray(*source.fTrackHitsPtr);
  //source.fTrackHitsPtr->Dump();
  //fTrackHitsPtr->Dump();
  
  fTrackPar = new TMatrixD(*source.fTrackPar); // track parameters
  fTrackParNew = new TMatrixD(*source.fTrackParNew); // track parameters
  fCovariance = new TMatrixD(*source.fCovariance); // covariance matrix
  fWeight = new TMatrixD(*source.fWeight); // weight matrix (inverse of covariance)

  return *this;
}

  //__________________________________________________________________________
void AliMUONTrackK::EvalCovariance(Double_t dZ)
{
  // Evaluate covariance (and weight) matrix for track candidate
  Double_t sigmaB, sigmaNonB, tanA, tanB, dAdY, rad, dBdX, dBdY;

  sigmaB = fgEventReconstructor->GetBendingResolution(); // bending resolution
  sigmaNonB = fgEventReconstructor->GetNonBendingResolution(); // non-bending resolution

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

  //(*fWeight)(4,4) = ((*fTrackPar)(4,0)*0.2)*((*fTrackPar)(4,0)*0.2); // error 20%
  (*fWeight)(4,4) = ((*fTrackPar)(4,0)*0.5)*((*fTrackPar)(4,0)*0.5); // error 50%

  // check whether the Invert method returns flag if matrix cannot be inverted,
  // and do not calculate the Determinant in that case !!!!
  if (fWeight->Determinant() != 0) {

    // fWeight->Invert();

    Int_t ifailWeight;
    mnvertLocalK(&((*fWeight)(0,0)), fgkSize,fgkSize,fgkSize,ifailWeight);
  } else {
    AliWarning(" Determinant fWeight=0:");
  }
  return;
}

  //__________________________________________________________________________
Bool_t AliMUONTrackK::KalmanFilter(Int_t ichamBeg, Int_t ichamEnd, Bool_t Back, Double_t zDipole1, Double_t zDipole2)
{
  // Follows track through detector stations 
  Bool_t miss, success;
  Int_t ichamb, iFB, iMin, iMax, dChamb, ichambOK, i;
  Int_t ihit, firstIndx, lastIndx, currIndx, dChambMiss, iDindx=0;
  Double_t zEnd, dChi2;
  AliMUONHitForRec *hitAdd, *firstHit, *lastHit, *hit;
  AliMUONRawCluster *clus;
  TClonesArray *rawclusters;
  hit = 0; clus = 0; rawclusters = 0;

  miss = kTRUE;
  success = kTRUE;
  Int_t endOfProp = 0;
  iFB = TMath::Sign(1,ichamEnd-ichamBeg);
  iMin = TMath::Min(ichamEnd,ichamBeg);
  iMax = TMath::Max(ichamEnd,ichamBeg);
  ichamb = ichamBeg;
  ichambOK = ichamb;

  // Get indices of the 1'st and last hits on the track candidate
  firstHit = (AliMUONHitForRec*) fTrackHitsPtr->First();
  lastHit = (AliMUONHitForRec*) fTrackHitsPtr->Last();
  firstIndx = fgHitForRec->IndexOf(firstHit);
  lastIndx = fgHitForRec->IndexOf(lastHit);
  currIndx = TMath::Abs (TMath::Max(firstIndx*iFB,lastIndx*iFB));
  if (Back) {
    // backpropagation
    currIndx = 2; 
    iDindx = 1;
    if (fRecover != 0) {
      // find hit with the highest Z
      Double_t zbeg = 0;
      for (i=0; i<fNTrackHits; i++) {
        hitAdd = (AliMUONHitForRec*) ((*fTrackHitsPtr)[i]);
        zEnd = hitAdd->GetZ();
	if (zEnd > zbeg) zbeg = zEnd;
	else {
	  currIndx = fNTrackHits - i + 2; //???
	  break;
	}
      } //for (Int_t i=0;
    }
  } else if (fRecover != 0) {
    Back = kTRUE; // dirty trick
    iDindx = -1;
    if (ichamBeg == 7 || ichamBeg == 8) currIndx = fNTrackHits - 2;
    else {
      Double_t zbeg = ((AliMUONHitForRec*)((*fTrackHitsPtr)[0]))->GetZ();
      for (i=1; i<fNTrackHits; i++) {
        hitAdd = (AliMUONHitForRec*) ((*fTrackHitsPtr)[i]);
        zEnd = hitAdd->GetZ();
	if (zEnd < zbeg) break;
      } //for (Int_t i=1;
      currIndx = fNTrackHits - i; //???
    }
  }

  while (ichamb>=iMin && ichamb<=iMax) {
  // Find the closest hit in Z, not belonging to the current plane
    if (Back) {
      // backpropagation
      hitAdd = (AliMUONHitForRec*) ((*fTrackHitsPtr)[fNTrackHits-currIndx]);
      zEnd = hitAdd->GetZ();
    } else {
      zEnd = -9999;
      for (ihit=currIndx+iFB; ihit>=0 && ihit<fgNOfPoints; ihit+=iFB) {
	hitAdd = (AliMUONHitForRec*) ((*fgHitForRec)[ihit]);
	//if (TMath::Abs(hitAdd->GetZ()-fPosition) > 0.1) {
	if (TMath::Abs(hitAdd->GetZ()-fPosition) > 0.5) {
	  zEnd = hitAdd->GetZ();
	  currIndx = ihit;
	  break;
	}
      }
    }
    if (zEnd<-999 && ichamb==ichamEnd) endOfProp = 1; // end-of-propagation
    else {
      // Check if there is a missing chamber
      if (zEnd<-999 || TMath::Abs(hitAdd->GetChamberNumber()-ichamb) > 1) {
	if (!Back && zEnd>-999) currIndx -= iFB;
	ichamb += iFB;
	zEnd = (&(fgMUON->Chamber(ichamb)))->Z();
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
      if (dChamb > 1) {
	dChambMiss = endOfProp;
        //Check if (iFB > 0) dChambMiss++;
        if (iFB > 0) {
	  if (TMath::Odd(ichambOK)) dChambMiss++;
	  else dChambMiss--;
	}
	//cout << dChamb << " " << ichambOK << " " << fgNOfPoints << endl;
	if (TMath::Odd(ichambOK) && dChamb > 3-dChambMiss) {
	  // missing station - abandon track
	  //cout << dChamb << " " << ichambOK << " " << fgNOfPoints << " " << 1/(*fTrackPar)(4,0) << endl;
	  /*
	  for (Int_t i1=0; i1<fgNOfPoints; i1++) {
	    cout << " Hit #" << ((AliMUONHitForRec*)((*fgHitForRec)[i1]))->GetChamberNumber() << " ";
	    cout << ((AliMUONHitForRec*)((*fgHitForRec)[i1]))->GetBendingCoor() << " ";
	    cout << ((AliMUONHitForRec*)((*fgHitForRec)[i1]))->GetNonBendingCoor() << " ";
	    cout << ((AliMUONHitForRec*)((*fgHitForRec)[i1]))->GetZ() << " " << " ";
	    cout << ((AliMUONHitForRec*)((*fgHitForRec)[i1]))->GetTHTrack() << endl;
	  }
	  //cout << endl;
	  */
	  /*
	  cout << fNTrackHits << endl;
	  for (Int_t i1=0; i1<fNTrackHits; i1++) {
	    hit = (AliMUONHitForRec*) ((*fTrackHitsPtr)[i1]);
	    printf(" * %d %10.4f %10.4f %10.4f", 
		   hit->GetChamberNumber(), hit->GetBendingCoor(), 
		   hit->GetNonBendingCoor(), hit->GetZ());
	    if (fgEventReconstructor->GetRecGeantHits()) { 
	      // from GEANT hits
	      printf(" %3d %3d \n", hit->GetGeantSignal(), hit->GetTHTrack());
	    } else {
	      // from raw clusters
	      rawclusters = fgMUON->RawClustAddress(hit->GetChamberNumber());
	      clus = (AliMUONRawCluster*) rawclusters->UncheckedAt(hit->GetHitNumber());
	      printf("%3d", clus->fTracks[1]-1); 
	      if (clus->fTracks[2] != 0) printf("%3d \n", clus->fTracks[2]-1);
	      else printf("\n");
	    }
	  }
	  */
	  if (fNTrackHits>2 && fRecover==0 && !(ichambOK==((AliMUONHitForRec*)((*fTrackHitsPtr)[0]))->GetChamberNumber())) {
	    // try to recover track later
	    Recover();
	  } 
	  return kFALSE;
	}
	//Check else if (TMath::Even(ichambOK) && dChamb > 2-endOfProp) {
	else if (TMath::Even(ichambOK) && dChamb > 2-dChambMiss) {
 	  // missing station - abandon track
	  //cout << dChamb << " " << ichambOK << " " << fgNOfPoints << " " << 1/(*fTrackPar)(4,0) << endl;
	  /*
	  for (Int_t i1=0; i1<fgNOfPoints; i1++) {
	    cout << " Hit #" << ((AliMUONHitForRec*)((*fgHitForRec)[i1]))->GetChamberNumber() << " ";
	    cout << ((AliMUONHitForRec*)((*fgHitForRec)[i1]))->GetBendingCoor() << " ";
	    cout << ((AliMUONHitForRec*)((*fgHitForRec)[i1]))->GetNonBendingCoor() << " ";
	    cout << ((AliMUONHitForRec*)((*fgHitForRec)[i1]))->GetZ() << " " << " ";
	    cout << ((AliMUONHitForRec*)((*fgHitForRec)[i1]))->GetTHTrack() << endl;
	  }
	  //cout << endl;
	  */
	  /*
	  cout << fNTrackHits << endl;
	  for (Int_t i1=0; i1<fNTrackHits; i1++) {
	    hit = (AliMUONHitForRec*) ((*fTrackHitsPtr)[i1]);
	    printf(" * %d %10.4f %10.4f %10.4f", 
		   hit->GetChamberNumber(), hit->GetBendingCoor(), 
		   hit->GetNonBendingCoor(), hit->GetZ());
	    if (fgEventReconstructor->GetRecGeantHits()) { 
	      // from GEANT hits
	      printf(" %3d %3d \n", hit->GetGeantSignal(), hit->GetTHTrack());
	    } else {
	      // from raw clusters
	      rawclusters = fgMUON->RawClustAddress(hit->GetChamberNumber());
	      clus = (AliMUONRawCluster*) rawclusters->UncheckedAt(hit->GetHitNumber());
	      printf("%3d", clus->fTracks[1]-1); 
	      if (clus->fTracks[2] != 0) printf("%3d \n", clus->fTracks[2]-1);
	      else printf("\n");
	    }
	  }
	  */
	  if (fNTrackHits>2 && fRecover==0 && !(ichambOK==((AliMUONHitForRec*)((*fTrackHitsPtr)[0]))->GetChamberNumber())) {
	    // try to recover track later
	    Recover();
	  } 
	  return kFALSE;
	}
      }
    }
    if (endOfProp != 0) break;

    // propagate to the found Z

    // Check if track steps into dipole
    if (fPosition>zDipole2 && zEnd<zDipole2) {
      //LinearPropagation(zDipole2-zBeg); 
      ParPropagation(zDipole2); 
      MSThin(1); // multiple scattering in the chamber
      WeightPropagation(zDipole2); // propagate weight matrix
      fPosition = fPositionNew;
      *fTrackPar = *fTrackParNew; 
      //MagnetPropagation(zEnd); 
      ParPropagation(zEnd); 
      WeightPropagation(zEnd);
      fPosition = fPositionNew;
    } 
    // Check if track steps out of dipole
    else if (fPosition>zDipole1 && zEnd<zDipole1) {
      //MagnetPropagation(zDipole1-zBeg); 
      ParPropagation(zDipole1); 
      MSThin(1); // multiple scattering in the chamber
      WeightPropagation(zDipole1);
      fPosition = fPositionNew;
      *fTrackPar = *fTrackParNew; 
      //LinearPropagation(zEnd-zDipole1); 
      ParPropagation(zEnd); 
      WeightPropagation(zEnd);
      fPosition = fPositionNew;
    } else {
      ParPropagation(zEnd);
      //MSThin(1); // multiple scattering in the chamber
      if (TMath::Abs(zEnd-fPosition) > 5) MSThin(1); // multiple scattering in the chamber
      WeightPropagation(zEnd);
      fPosition = fPositionNew;
    }

    // Add measurement
    if (fRecover != 0 && hitAdd == fSkipHit && !miss) {
      // recovered track - remove the hit
      miss = kTRUE;
      ichamb = hitAdd->GetChamberNumber();
      if (fRecover == 1) {
	// remove the last hit
	fTrackHitsPtr->Remove((TObjArray*)hitAdd); // remove hit
	fNTrackHits --;
	hitAdd->SetNTrackHits(hitAdd->GetNTrackHits()-1); // unmark hit 
      } else {
	// remove the hits
	for (i=fNTrackHits-1; i>1; i--) {
	  hitAdd = (AliMUONHitForRec*)((*fTrackHitsPtr)[i]);
	  fTrackHitsPtr->Remove((TObjArray*)hitAdd); // remove hit
	  hitAdd->SetNTrackHits(hitAdd->GetNTrackHits()-1); // unmark hit 
	  fNTrackHits --;
	  if (hitAdd == fSkipHit) break;
	} // for (i=fNTrackHits-1;
      }
      Back = kFALSE;
      fRecover =0; // ????????? Dec-17-2001
      ichambOK = ((AliMUONHitForRec*)((*fTrackHitsPtr)[fNTrackHits-1]))->GetChamberNumber();
      currIndx = fgHitForRec->IndexOf(fSkipHit);
    }

    if (Back && !miss) {
      // backward propagator
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
      fChi2 += dChi2; // Chi2
      if (ichamb==ichamEnd) break; 
      currIndx += iDindx;
    } else {
      // forward propagator
      if (miss || !FindPoint(ichamb,zEnd,currIndx,iFB,hitAdd)) {
	// missing point
	*fTrackPar = *fTrackParNew; 
      } else {
	//add point
	fTrackHitsPtr->Add((TObjArray*)hitAdd); // add hit
	fNTrackHits ++;
	hitAdd->SetNTrackHits(hitAdd->GetNTrackHits()+1); // mark hit as being on track
	ichambOK = ichamb;
	currIndx = fgHitForRec->IndexOf(hitAdd); // Check
      }
    }
  } // while
  cout << fNTrackHits << " " << fChi2 << " " << 1/(*fTrackPar)(4,0) << " " << fPosition << endl;
  return success;
}

  //__________________________________________________________________________
void AliMUONTrackK::ParPropagation(Double_t zEnd)
{
  // Propagation of track parameters to zEnd
  Int_t iFB, nTries;
  Double_t dZ, step, distance, charge;
  Double_t vGeant3[7], vGeant3New[7];

  nTries = 0;
  // First step using linear extrapolation
  dZ = zEnd - fPosition;
  iFB = (Int_t)TMath::Sign(Double_t(1.0),dZ);
  step = dZ/TMath::Cos((*fTrackPar)(2,0))/TMath::Cos((*fTrackPar)(3,0)); // linear estimate
  charge = iFB*TMath::Sign(Double_t(1.0),(*fTrackPar)(4,0));
  fPositionNew = fPosition;
  *fTrackParNew = *fTrackPar;
  SetGeantParam(vGeant3,iFB);

  // Check if overstep
  do {
    step = TMath::Abs(step);
    // Propagate parameters
    extrap_onestep_rungekutta(charge,step,vGeant3,vGeant3New);
    distance = zEnd - vGeant3New[2];
    step *= dZ/(vGeant3New[2]-fPositionNew);
    nTries ++;
  } while (distance*iFB < 0 && TMath::Abs(distance) > fgkEpsilon);

  GetFromGeantParam(vGeant3New,iFB);

  // Position ajustment (until within tolerance)
  while (TMath::Abs(distance) > fgkEpsilon) {
    dZ = zEnd - fPositionNew;
    iFB = (Int_t)TMath::Sign(Double_t(1.0),dZ);
    step = dZ/TMath::Cos((*fTrackParNew)(2,0))/TMath::Cos((*fTrackParNew)(3,0));
    step = TMath::Abs(step);
    SetGeantParam(vGeant3,iFB);
    do {
      // binary search
      // Propagate parameters
      extrap_onestep_rungekutta(charge,step,vGeant3,vGeant3New);
      distance = zEnd - vGeant3New[2];
      step /= 2;
      nTries ++;
      if (nTries > fgkTriesMax) {
	cout << " ***** ParPropagation: too many tries " << nTries << endl;
	exit(0);
      }
    } while (distance*iFB < 0);

    GetFromGeantParam(vGeant3New,iFB);
  }
  //cout << nTries << endl;
  return;
}
/*
  //__________________________________________________________________________
void AliMUONTrackK::WeightPropagation(void)
{
  // Propagation of the weight matrix
  // W = DtWD, where D is Jacobian 

  // !!! not implemented TMatrixD weight1(*fJacob,TMatrixD::kAtBA,*fWeight); // DtWD
  TMatrixD weight1(*fWeight,TMatrixD::kMult,*fJacob); // WD
  *fWeight = TMatrixD(*fJacob,TMatrixD::kTransposeMult,weight1); // DtWD
  return;
}
*/
  //__________________________________________________________________________
void AliMUONTrackK::WeightPropagation(Double_t zEnd)
{
  // Propagation of the weight matrix
  // W = DtWD, where D is Jacobian 
  Int_t i, j;
  Double_t dPar;

  TMatrixD jacob(fgkSize,fgkSize);
  jacob = 0;

  // Save initial and propagated parameters
  TMatrixD trackPar0 = *fTrackPar;
  TMatrixD trackParNew0 = *fTrackParNew;
  Double_t savePosition = fPositionNew;

  // Get covariance matrix
  *fCovariance = *fWeight;
  // check whether the Invert method returns flag if matrix cannot be inverted,
  // and do not calculate the Determinant in that case !!!!
  if (fCovariance->Determinant() != 0) {
    //   fCovariance->Invert();
    Int_t ifailCov;
    mnvertLocalK(&((*fCovariance)(0,0)), fgkSize,fgkSize,fgkSize,ifailCov);
  } else {
    AliWarning(" Determinant fCovariance=0:");
  }

  // Loop over parameters to find change of the initial vs propagated ones
  zEnd = fPosition;
  fPosition = fPositionNew;
  for (i=0; i<fgkSize; i++) {
    dPar = TMath::Sqrt((*fCovariance)(i,i));
    *fTrackPar = trackParNew0;
    (*fTrackPar)(i,0) += dPar;
    ParPropagation(zEnd);
    for (j=0; j<fgkSize; j++) {
      jacob(j,i) = ((*fTrackParNew)(j,0)-trackPar0(j,0))/dPar;
    }
  }

  //jacob->Print();
  //trackParNew0.Print();
  //TMatrixD par1(jacob,TMatrixD::kMult,trackPar0); //
  //par1.Print();
  /*
  if (jacob.Determinant() != 0) {
    //  jacob.Invert();
  } else {
    cout << " ***** Warning in WeightPropagation: Determinant jacob=0:" << endl;
  }
  */
  TMatrixD weight1(*fWeight,TMatrixD::kMult,jacob); // WD
  *fWeight = TMatrixD(jacob,TMatrixD::kTransposeMult,weight1); // DtWD
  //fWeight->Print();

  // Restore initial and propagated parameters
  *fTrackPar = trackPar0;
  *fTrackParNew = trackParNew0;
  fPosition = zEnd;
  fPositionNew = savePosition;
  return;
}

  //__________________________________________________________________________
Bool_t AliMUONTrackK::FindPoint(Int_t ichamb, Double_t zEnd, Int_t currIndx, Int_t iFB, AliMUONHitForRec *&hitAdd)
{
  // Picks up point within a window for the chamber No ichamb 
  // Split the track if there are more than 1 hit
  Int_t ihit, nRecTracks;
  Double_t windowB, windowNonB, dChi2Tmp=0, dChi2, y, x, savePosition=0;
  TClonesArray *trackPtr;
  AliMUONHitForRec *hit, *hitLoop;
  AliMUONTrackK *trackK;

  Bool_t ok = kFALSE;
  //sigmaB = fgEventReconstructor->GetBendingResolution(); // bending resolution
  //sigmaNonB = fgEventReconstructor->GetNonBendingResolution(); // non-bending resolution
  *fCovariance = *fWeight;
  // check whether the Invert method returns flag if matrix cannot be inverted,
  // and do not calculate the Determinant in that case !!!!
  if (fCovariance->Determinant() != 0) {
    //  fCovariance->Invert();

      Int_t ifailCov;
      mnvertLocalK(&((*fCovariance)(0,0)), fgkSize,fgkSize,fgkSize,ifailCov);
  } else {
    AliWarning("Determinant fCovariance=0:");
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
  Int_t nHitsOK = 0;

  for (ihit=currIndx; ihit>=0 && ihit<fgNOfPoints; ihit+=iFB) {
    hit = (AliMUONHitForRec*) ((*fgHitForRec)[ihit]);
    if (hit->GetChamberNumber() == ichamb) {
      //if (TMath::Abs(hit->GetZ()-zEnd) < 0.1) {
      if (TMath::Abs(hit->GetZ()-zEnd) < 0.5) {
        if (TMath::Abs(hit->GetZ()-zEnd) > 0.1) {
	  // adjust position: for multiple hits in the chamber
	  // (mostly (only?) for GEANT hits)
	  zEnd = hit->GetZ();
	  *fTrackPar = *fTrackParNew;
	  ParPropagation(zEnd);
	  WeightPropagation(zEnd);
	  fPosition = fPositionNew;
	  *fTrackPar = *fTrackParNew;
	  // Get covariance
	  *fCovariance = *fWeight;
	  if (fCovariance->Determinant() != 0) {
	    //fCovariance->Invert();
	    Int_t ifailCov;
	    mnvertLocalK(&((*fCovariance)(0,0)), fgkSize,fgkSize,fgkSize,ifailCov);
	  } else {
	    AliWarning("Determinant fCovariance=0:");
	  }
	}
	y = hit->GetBendingCoor();
	x = hit->GetNonBendingCoor();
	windowB = fgkNSigma*TMath::Sqrt((*fCovariance)(0,0)+hit->GetBendingReso2());
	windowNonB = fgkNSigma*TMath::Sqrt((*fCovariance)(1,1)+hit->GetNonBendingReso2());
	if (TMath::Abs((*fTrackParNew)(0,0)-y) <= windowB &&
	    TMath::Abs((*fTrackParNew)(1,0)-x) <= windowNonB) {
	  // Vector of measurements and covariance matrix
	  point.Zero();
	  point(0,0) = y;
	  point(1,0) = x;
	  pointWeight(0,0) = 1/hit->GetBendingReso2();
	  pointWeight(1,1) = 1/hit->GetNonBendingReso2();
	  TryPoint(point,pointWeight,trackPar,dChi2);
	  if (TMath::Abs(1./(trackPar)(4,0)) < fgEventReconstructor->GetMinBendingMomentum()) continue; // p < p_min - next hit
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
	  } else {
	    // branching: create a new track
	    trackPtr = fgEventReconstructor->GetRecTracksPtr();
	    nRecTracks = fgEventReconstructor->GetNRecTracks();
	    trackK = new ((*trackPtr)[nRecTracks])
	    	     AliMUONTrackK(*this); // dummy copy constructor
	    *trackK = *this;
	    fgEventReconstructor->SetNRecTracks(nRecTracks+1);
	    //cout << " ******** New track: " << ichamb << " " << hit->GetTHTrack() << " " << 1/(trackPar)(4,0) << " " << hit->GetBendingCoor() << " " << fNTrackHits << " " << nRecTracks << endl;
	    trackK->fRecover = 0;
	    *(trackK->fTrackPar) = trackPar;
	    *(trackK->fWeight) += pointWeight; 
	    trackK->fChi2 += dChi2;
	    // Mark hits as being on 2 tracks
	    for (Int_t i=0; i<fNTrackHits; i++) {
	      hitLoop = (AliMUONHitForRec*) ((*fTrackHitsPtr)[i]);
	      hitLoop->SetNTrackHits(hitLoop->GetNTrackHits()+1); 
	      /*
	      cout << " ** ";
	      cout << hitLoop->GetChamberNumber() << " ";
	      cout << hitLoop->GetBendingCoor() << " ";
	      cout << hitLoop->GetNonBendingCoor() << " ";
	      cout << hitLoop->GetZ() << " " << " ";
	      cout << hitLoop->GetGeantSignal() << " " << " ";
	      cout << hitLoop->GetTHTrack() << endl;
	      printf(" ** %d %10.4f %10.4f %10.4f %d %d \n", 
		     hitLoop->GetChamberNumber(), hitLoop->GetBendingCoor(), 
		     hitLoop->GetNonBendingCoor(), hitLoop->GetZ(), 
		     hitLoop->GetGeantSignal(), hitLoop->GetTHTrack());
	      */
	    }
	    //add point
	    trackK->fTrackHitsPtr->Add((TObjArray*)hit); // add hit
	    trackK->fNTrackHits ++;
	    hit->SetNTrackHits(hit->GetNTrackHits()+1); // mark hit as being on track
	    if (ichamb == 9) {
	      // the last chamber
	      trackK->fTrackDir = -1;
	      trackK->fBPFlag = kTRUE; 
	    }
	  }
	}
      }
    } else break; // different chamber
  } // for (ihit=currIndx;
  if (ok) {
    *fTrackPar = trackParTmp;
    *fWeight = saveWeight;
    *fWeight += pointWeightTmp; 
    fChi2 += dChi2Tmp; // Chi2
    // Restore members
    fPosition = savePosition;
  }
  return ok;
}

  //__________________________________________________________________________
void AliMUONTrackK::TryPoint(TMatrixD &point, const TMatrixD &pointWeight, TMatrixD &trackParTmp, Double_t &dChi2)
{
  // Adds a measurement point (modifies track parameters and computes
  // change of Chi2)

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

    //  wu.Invert();
     Int_t ifailWU;
      mnvertLocalK(&((wu)(0,0)), fgkSize,fgkSize,fgkSize,ifailWU);
  } else {
    AliWarning("Determinant wu=0:");
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
  // Adds multiple scattering in a thin layer (only angles are affected)
  Double_t cosAlph, cosBeta, momentum, velo, path, theta0;

  // check whether the Invert method returns flag if matrix cannot be inverted,
  // and do not calculate the Determinant in that case !!!!
  if (fWeight->Determinant() != 0) {
    //fWeight->Invert(); // covariance

    Int_t ifailWeight;
    mnvertLocalK(&((*fWeight)(0,0)), fgkSize,fgkSize,fgkSize,ifailWeight);
  } else {
    AliWarning("Determinant fWeight=0:");
  }

  cosAlph = TMath::Cos((*fTrackParNew)(2,0));
  cosBeta = TMath::Cos((*fTrackParNew)(3,0));
  momentum = 1/(*fTrackParNew)(4,0); // particle momentum
  //velo = momentum/TMath::Sqrt(momentum*momentum+muonMass*muonMass); // velocity/c for muon hypothesis
  velo = 1; // relativistic
  path = fgEventReconstructor->GetChamberThicknessInX0()/cosAlph/cosBeta; // path length
  theta0 = 0.0136/velo/momentum*TMath::Sqrt(path)*(1+0.038*TMath::Log(path)); // projected scattering angle

  (*fWeight)(2,2) += sign*theta0/cosBeta*theta0/cosBeta; // alpha
  (*fWeight)(3,3) += sign*theta0*theta0; // beta
  //fWeight->Invert(); // weight

  Int_t ifailWeight;
  mnvertLocalK(&((*fWeight)(0,0)), fgkSize,fgkSize,fgkSize,ifailWeight);
  return;
}
  //__________________________________________________________________________
void AliMUONTrackK::StartBack(void)
{
  // Starts backpropagator
  
  fBPFlag = kTRUE;
  fChi2 = 0;
  for (Int_t i=0; i<fgkSize; i++) {
    for (Int_t j=0; j<fgkSize; j++) {
      if (j==i) (*fWeight)(i,i) /= 100;
      //if (j==i) (*fWeight)(i,i) /= fNTrackHits*fNTrackHits;
      else (*fWeight)(j,i) = 0;
    }
  }
}

  //__________________________________________________________________________
void AliMUONTrackK::SetGeantParam(Double_t *VGeant3, Int_t iFB)
{
  // Set vector of Geant3 parameters pointed to by "VGeant3"
  // from track parameters 

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
  // Get track parameters from vector of Geant3 parameters pointed 
  // to by "VGeant3"

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
  // Computes "track quality" from Chi2 (if iChi2==0) or vice versa

  if (fChi2 > 250) {
    cout << " ***** Too high Chi2: " << fChi2 << endl;
    fChi2 = 250;
    //   exit(0);
  }
  if (iChi2 == 0) fChi2 = fNTrackHits + (250.-fChi2)/251;
  else fChi2 = 250 - (fChi2-fNTrackHits)*251;
}

  //__________________________________________________________________________
Int_t AliMUONTrackK::Compare(const TObject* trackK) const
{
  // "Compare" function to sort with decreasing "track quality".
  // Returns +1 (0, -1) if quality of current track
  // is smaller than (equal to, larger than) quality of trackK

  if (fChi2 < ((AliMUONTrackK*)trackK)->fChi2) return(+1);
  else if (fChi2 == ((AliMUONTrackK*)trackK)->fChi2) return(0);
  else return(-1);
}

  //__________________________________________________________________________
Bool_t AliMUONTrackK::KeepTrack(AliMUONTrackK* track0) const
{
  // Check whether or not to keep current track 
  // (keep, if it has less than half of common hits with track0)
  Int_t hitsInCommon, nHits0, i, j, nTrackHits2;
  AliMUONHitForRec *hit0, *hit1;

  hitsInCommon = 0;
  nHits0 = track0->fNTrackHits;
  nTrackHits2 = fNTrackHits/2;

  for (i=0; i<nHits0; i++) {
    // Check if hit belongs to several tracks
    hit0 = (AliMUONHitForRec*) (*track0->fTrackHitsPtr)[i]; 
    if (hit0->GetNTrackHits() == 1) continue; 
    for (j=0; j<fNTrackHits; j++) {
      hit1 = (AliMUONHitForRec*) (*fTrackHitsPtr)[j]; 
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
  // Kill track candidate
  Int_t i;
  AliMUONHitForRec *hit;

  if (fTrackHitsPtr) {
    // Remove track mark from hits
    for (i=0; i<fNTrackHits; i++) {
      hit = (AliMUONHitForRec*) (*fTrackHitsPtr)[i]; 
      hit->SetNTrackHits(hit->GetNTrackHits()-1); 
    }
  }
  fgEventReconstructor->GetRecTracksPtr()->Remove(this);
}

  //__________________________________________________________________________
void AliMUONTrackK::Branson(void)
{
  // Propagates track to the vertex thru absorber using Branson correction
  // (makes use of the AliMUONTrackParam class)
 
  AliMUONTrackParam *trackParam = new AliMUONTrackParam();
  trackParam->SetBendingCoor((*fTrackPar)(0,0));
  trackParam->SetNonBendingCoor((*fTrackPar)(1,0));
  trackParam->SetBendingSlope(TMath::Tan((*fTrackPar)(2,0)));
  trackParam->SetNonBendingSlope(TMath::Tan((*fTrackPar)(3,0))/TMath::Cos((*fTrackPar)(2,0)));
  trackParam->SetInverseBendingMomentum((*fTrackPar)(4,0)/TMath::Cos((*fTrackPar)(3,0)));
  trackParam->SetZ(fPosition);

  trackParam->ExtrapToVertex(0.,0.,0.);

  (*fTrackPar)(0,0) = trackParam->GetBendingCoor();
  (*fTrackPar)(1,0) = trackParam->GetNonBendingCoor();
  (*fTrackPar)(2,0) = TMath::ATan(trackParam->GetBendingSlope());
  (*fTrackPar)(3,0) = TMath::ATan(TMath::Cos((*fTrackPar)(2,0))*trackParam->GetNonBendingSlope());
  (*fTrackPar)(4,0) = TMath::Cos((*fTrackPar)(3,0))*trackParam->GetInverseBendingMomentum();
  fPosition = trackParam->GetZ();
  delete trackParam;
  cout << 1/(*fTrackPar)(4,0) << " " << fPosition << " " << (*fTrackPar)(0,0) << endl;

  // Get covariance matrix
  *fCovariance = *fWeight;
  if (fCovariance->Determinant() != 0) {
    //    fCovariance->Invert();

      Int_t ifailCov;
      mnvertLocalK(&((*fCovariance)(0,0)), fgkSize,fgkSize,fgkSize,ifailCov);
  } else {
    AliWarning("Determinant fCovariance=0:");
  }
}

  //__________________________________________________________________________
void AliMUONTrackK::GoToZ(Double_t zEnd)
{
  // Propagates track to given Z

  ParPropagation(zEnd);
  MSThin(1); // multiple scattering in the chamber
  WeightPropagation(zEnd);
  fPosition = fPositionNew;
  *fTrackPar = *fTrackParNew; 
}

  //__________________________________________________________________________
void AliMUONTrackK::GoToVertex(void)
{
  // Version 3.08
  // Propagates track to the vertex
  // All material constants are taken from AliRoot

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
  static Double_t zPos[10] = {90, 105, 315, 443, 468};
  // R > 1
  // R < 1

  Double_t dZ, r0Norm, x0, deltaP, dChi2, pTotal, pOld;
  AliMUONHitForRec *hit;
  AliMUONRawCluster *clus;
  TClonesArray *rawclusters;

  // First step to the rear end of the absorber
  Double_t zRear = 503;
  GoToZ(zRear);
  Double_t tan3 = TMath::Tan(3./180*TMath::Pi());

  // Go through absorber
  pOld = 1/(*fTrackPar)(4,0);
  Double_t r0Rear = (*fTrackPar)(0,0)*(*fTrackPar)(0,0) + 
                    (*fTrackPar)(1,0)*(*fTrackPar)(1,0);
  r0Rear = TMath::Sqrt(r0Rear)/fPosition/tan3;
  r0Norm = r0Rear;
  for (Int_t i=4; i>=0; i--) {
    ParPropagation(zPos[i]);
    WeightPropagation(zPos[i]);
    dZ = TMath::Abs (fPositionNew-fPosition);
    if (r0Norm > 1) x0 = x01[i];
    else x0 = x02[i];
    MSLine(dZ,x0); // multiple scattering in the medium (linear approximation)
    fPosition = fPositionNew;
    *fTrackPar = *fTrackParNew; 
    r0Norm = (*fTrackPar)(0,0)*(*fTrackPar)(0,0) + 
             (*fTrackPar)(1,0)*(*fTrackPar)(1,0);
    r0Norm = TMath::Sqrt(r0Norm)/fPosition/tan3;
  }
  // Correct momentum for energy losses
  pTotal = 1/TMath::Abs((*fTrackPar)(4,0));
  Double_t p0 = pTotal;
  for (Int_t j=0; j<2; j++) {
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
    if (r0Rear < 1) {
      //W
      if (p0<15) {
	deltaP = 2.737 + 0.0494*p0 - 0.001123*p0*p0;
      } else {
	deltaP = 3.0643 + 0.01346*p0;
      }
    } else {
      //Pb
      if (p0<15) {
	deltaP  = 2.1380 + 0.0351*p0 - 0.000853*p0*p0;
      } else {
	deltaP = 2.407 + 0.00702*p0;
      }
    }

    p0 = pTotal + deltaP/TMath::Cos((*fTrackPar)(2,0))/TMath::Cos((*fTrackPar)(3,0));
  }
  (*fTrackPar)(4,0) = 1/p0*TMath::Sign((Double_t)1.,(*fTrackPar)(4,0));

  // Go to the vertex
  ParPropagation((Double_t)0.);
  WeightPropagation((Double_t)0.);
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
  cout << pOld << " " << 1/(*fTrackPar)(4,0) << " " << dChi2 << " " << fChi2 << " " << fNTrackHits << endl;
  for (Int_t i1=0; i1<fNTrackHits; i1++) {
    hit =  (AliMUONHitForRec*) ((*fTrackHitsPtr)[i1]);
    printf ("%4d", hit->GetChamberNumber()); 
    //cout << ((AliMUONHitForRec*)((*fTrackHitsPtr)[i1]))->GetChamberNumber() << " ";
  }
  cout << endl;
  for (Int_t i1=0; i1<fNTrackHits; i1++) {
    hit =  (AliMUONHitForRec*) ((*fTrackHitsPtr)[i1]);
    //cout << ((AliMUONHitForRec*)((*fTrackHitsPtr)[i1]))->GetHitNumber() << " ";
    //cout << ((AliMUONHitForRec*)((*fTrackHitsPtr)[i1]))->GetZ() << " ";
    printf ("%4d", fgHitForRec->IndexOf(hit)); 
    //cout << fgHitForRec->IndexOf(((AliMUONHitForRec*)((*fTrackHitsPtr)[i1]))) << " ";
  }
  cout << endl;
  if (fgEventReconstructor->GetRecGeantHits()) { 
      // from GEANT hits
    for (Int_t i1=0; i1<fNTrackHits; i1++) {
      hit =  (AliMUONHitForRec*) ((*fTrackHitsPtr)[i1]);
      cout << hit->GetTHTrack() + hit->GetGeantSignal()*10000 << " ";
    }
  } else {
    // from raw clusters
    for (Int_t i1=0; i1<fNTrackHits; i1++) {
      hit =  (AliMUONHitForRec*) ((*fTrackHitsPtr)[i1]);
      rawclusters = fgMUON->GetMUONData()->RawClusters(hit->GetChamberNumber());
      clus = (AliMUONRawCluster*) rawclusters->UncheckedAt(hit->GetHitNumber());
      printf ("%4d", clus->GetTrack(1) - 1); 
      //cout << clus->fTracks[1] - 1 << " ";
    }
    cout << endl;
    for (Int_t i1=0; i1<fNTrackHits; i1++) {
      hit =  (AliMUONHitForRec*) ((*fTrackHitsPtr)[i1]);
      rawclusters = fgMUON->GetMUONData()->RawClusters(hit->GetChamberNumber());
      clus = (AliMUONRawCluster*) rawclusters->UncheckedAt(hit->GetHitNumber());
      if (clus->GetTrack(2) != 0) printf ("%4d", clus->GetTrack(2) - 1);
      else printf ("%4s", "   ");
      //if (clus->fTracks[2] != 0) cout << clus->fTracks[2] - 1 << " ";
    }
  }
  cout << endl;
  for (Int_t i1=0; i1<fNTrackHits; i1++) {
    //cout << ((AliMUONHitForRec*)((*fTrackHitsPtr)[i1]))->GetHitNumber() << " ";
    cout << ((AliMUONHitForRec*)((*fTrackHitsPtr)[i1]))->GetZ() << " ";
    //cout << fgHitForRec->IndexOf(((AliMUONHitForRec*)((*fTrackHitsPtr)[i1]))) << " ";
  }
  cout << endl;
  cout << "---------------------------------------------------" << endl;

  // Get covariance matrix
  *fCovariance = *fWeight;
  if (fCovariance->Determinant() != 0) {
    //   fCovariance->Invert();

      Int_t ifailCov;
      mnvertLocalK(&((*fCovariance)(0,0)), fgkSize,fgkSize,fgkSize,ifailCov);
  } else {
   AliWarning("Determinant fCovariance=0:" );
  }
}

  //__________________________________________________________________________
void AliMUONTrackK::MSLine(Double_t dZ, Double_t x0)
{
  // Adds multiple scattering in a thick layer for linear propagation

  Double_t cosAlph = TMath::Cos((*fTrackPar)(2,0));
  Double_t tanAlph = TMath::Tan((*fTrackPar)(2,0));
  Double_t cosBeta = TMath::Cos((*fTrackPar)(3,0));
  Double_t sinBeta;
  sinBeta = TMath::Sin((*fTrackPar)(3,0));
  Double_t tanBeta = TMath::Tan((*fTrackPar)(3,0));
  Double_t momentum = 1/(*fTrackPar)(4,0);
  Double_t velo = 1; // relativistic velocity
  Double_t step = TMath::Abs(dZ)/cosAlph/cosBeta; // step length

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

       Int_t ifailCov;
       mnvertLocalK(&((*fCovariance)(0,0)), fgkSize,fgkSize,fgkSize,ifailCov);
  } else {
    AliWarning("Determinant fCovariance=0:" );
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

  (*fCovariance)(1,2) += dl2*(dXdT*dAdT+dXdB*dAdB); // <xa>
  (*fCovariance)(2,1) = (*fCovariance)(1,2);

  (*fCovariance)(1,3) += dl2*dXdB; // <xb>
  (*fCovariance)(3,1) = (*fCovariance)(1,3);

  (*fCovariance)(2,3) += theta02*step*dAdB; // <ab>
  (*fCovariance)(3,2) = (*fCovariance)(2,3);

  // Get weight matrix
  *fWeight = *fCovariance;
  if (fWeight->Determinant() != 0) {
    //  fWeight->Invert();

       Int_t ifailWeight;
       mnvertLocalK(&((*fWeight)(0,0)), fgkSize,fgkSize,fgkSize,ifailWeight);
  } else {
    AliWarning("Determinant fWeight=0:");
  }
}
 
  //__________________________________________________________________________
void AliMUONTrackK::Recover(void)
{
  // Adds new failed track(s) which can be tried to be recovered
  Int_t nRecTracks, ichamb;
  TClonesArray *trackPtr;
  AliMUONTrackK *trackK;

  //cout << " ******** Enter Recover " << endl;
  //return;
  trackPtr = fgEventReconstructor->GetRecTracksPtr();

  // The last hit will be removed
  nRecTracks = fgEventReconstructor->GetNRecTracks();

  // Check if the track candidate doesn't exist yet
  for (Int_t i=0; i<nRecTracks; i++) {
    trackK = (AliMUONTrackK*) ((*trackPtr)[i]);
    if (trackK->fNTrackHits == 2 && trackK->GetRecover() == 0) continue;
    if (trackK == this) continue;
    //if (trackK->GetRecover() != 1) continue;
    if (trackK->fNTrackHits >= fNTrackHits-1) {
      /*
      for (Int_t j=0; j<fNTrackHits-1; j++) {
	if ((*trackK->fTrackHitsPtr)[j] != ((*fTrackHitsPtr)[j])) break;
	return;
      } // for (Int_t j=0;
      */
      if ((*trackK->fTrackHitsPtr)[0] == ((*fTrackHitsPtr)[0])) return;
    }
  } // for (Int_t i=0;

  cout << " ******** Enter Recover " << endl;
  trackK = new ((*trackPtr)[nRecTracks]) AliMUONTrackK(fStartSegment); 
  fgEventReconstructor->SetNRecTracks(nRecTracks+1);
  trackK->fRecover = 1;
  trackK->fSkipHit = (AliMUONHitForRec*) ((*fTrackHitsPtr)[fNTrackHits-1]);
  trackK->fNTrackHits = fNTrackHits;
  delete trackK->fTrackHitsPtr; // not efficient ?
  trackK->fTrackHitsPtr = new TObjArray(*fTrackHitsPtr);
  cout << nRecTracks << " " << trackK->fRecover << endl;

  // The hit before missing chamber will be removed
  Int_t ichamBeg = ((AliMUONHitForRec*)((*fTrackHitsPtr)[0]))->GetChamberNumber();
  Int_t indxSkip = -1;
  if (ichamBeg == 9) {
    // segment in the last station
    // look for the missing chamber
    for (Int_t i=1; i<fNTrackHits; i++) {
      ichamb = ((AliMUONHitForRec*)((*fTrackHitsPtr)[i]))->GetChamberNumber();
      if (TMath::Abs(ichamBeg-ichamb)>1 && i>2) {
	indxSkip = i;
	break;
      }
      ichamBeg = ichamb;
    } // for (Int_t i=1;
  } else {
    // in the last but one station
    for (Int_t i=1; i<fNTrackHits; i++) {
      ichamb = ((AliMUONHitForRec*)((*fTrackHitsPtr)[i]))->GetChamberNumber();
      if (TMath::Abs(ichamBeg-ichamb)>1 && ichamb<4) {
	indxSkip = i;
	break;
      }
      ichamBeg = ichamb;
    } // for (Int_t i=1;
  }
  if (indxSkip < 0) return;
  
  // Check if the track candidate doesn't exist yet
  for (Int_t i=0; i<nRecTracks; i++) {
    trackK = (AliMUONTrackK*) ((*trackPtr)[i]);
    if (trackK->fNTrackHits == 2 && trackK->GetRecover() == 0) continue;
    if (trackK == this) continue;
    //if (trackK->GetRecover() != 1) continue;
    if (trackK->fNTrackHits >= indxSkip-1) {
      /*
      for (Int_t j=0; j<indxSkip-1; j++) {
	if ((*trackK->fTrackHitsPtr)[j] != ((*fTrackHitsPtr)[j])) break;
	return;
      } // for (Int_t j=0;
      */
      if ((*trackK->fTrackHitsPtr)[0] == ((*fTrackHitsPtr)[0])) return;
    }
  } // for (Int_t i=0;

  nRecTracks = fgEventReconstructor->GetNRecTracks();
  trackK = new ((*trackPtr)[nRecTracks]) AliMUONTrackK(fStartSegment); 
  fgEventReconstructor->SetNRecTracks(nRecTracks+1);
  trackK->fRecover = 2;
  trackK->fSkipHit = (AliMUONHitForRec*) ((*fTrackHitsPtr)[indxSkip-1]);
  trackK->fNTrackHits = fNTrackHits;
  delete trackK->fTrackHitsPtr; // not efficient ?
  trackK->fTrackHitsPtr = new TObjArray(*fTrackHitsPtr);
  cout << nRecTracks << " " << trackK->fRecover << endl;
}

//______________________________________________________________________________
 void mnvertLocalK(Double_t *a, Int_t l, Int_t, Int_t n, Int_t &ifail)
{
//*-*-*-*-*-*-*-*-*-*-*-*Inverts a symmetric matrix*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                    ==========================
//*-*        inverts a symmetric matrix.   matrix is first scaled to
//*-*        have all ones on the diagonal (equivalent to change of units)
//*-*        but no pivoting is done since matrix is positive-definite.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  // taken from TMinuit package of Root (l>=n)
  // fVERTs, fVERTq and fVERTpp changed to localVERTs, localVERTq and localVERTpp
  //  Double_t localVERTs[n], localVERTq[n], localVERTpp[n];
  Double_t * localVERTs = new Double_t[n];
  Double_t * localVERTq = new Double_t[n];
  Double_t * localVERTpp = new Double_t[n];
  // fMaxint changed to localMaxint
  Int_t localMaxint = n;

    /* System generated locals */
    Int_t aOffset;

    /* Local variables */
    Double_t si;
    Int_t i, j, k, kp1, km1;

    /* Parameter adjustments */
    aOffset = l + 1;
    a -= aOffset;

    /* Function Body */
    ifail = 0;
    if (n < 1) goto L100;
    if (n > localMaxint) goto L100;
//*-*-                  scale matrix by sqrt of diag elements
    for (i = 1; i <= n; ++i) {
        si = a[i + i*l];
        if (si <= 0) goto L100;
        localVERTs[i-1] = 1 / TMath::Sqrt(si);
    }
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) {
            a[i + j*l] = a[i + j*l]*localVERTs[i-1]*localVERTs[j-1];
        }
    }
//*-*-                                       . . . start main loop . . . .
    for (i = 1; i <= n; ++i) {
        k = i;
//*-*-                  preparation for elimination step1
        if (a[k + k*l] != 0) localVERTq[k-1] = 1 / a[k + k*l];
        else goto L100;
        localVERTpp[k-1] = 1;
        a[k + k*l] = 0;
        kp1 = k + 1;
        km1 = k - 1;
        if (km1 < 0) goto L100;
        else if (km1 == 0) goto L50;
        else               goto L40;
L40:
        for (j = 1; j <= km1; ++j) {
            localVERTpp[j-1] = a[j + k*l];
            localVERTq[j-1]  = a[j + k*l]*localVERTq[k-1];
            a[j + k*l]   = 0;
        }
L50:
        if (k - n < 0) goto L51;
        else if (k - n == 0) goto L60;
        else                goto L100;
L51:
        for (j = kp1; j <= n; ++j) {
            localVERTpp[j-1] = a[k + j*l];
            localVERTq[j-1]  = -a[k + j*l]*localVERTq[k-1];
            a[k + j*l]   = 0;
        }
//*-*-                  elimination proper
L60:
        for (j = 1; j <= n; ++j) {
            for (k = j; k <= n; ++k) { a[j + k*l] += localVERTpp[j-1]*localVERTq[k-1]; }
        }
    }
//*-*-                  elements of left diagonal and unscaling
    for (j = 1; j <= n; ++j) {
        for (k = 1; k <= j; ++k) {
            a[k + j*l] = a[k + j*l]*localVERTs[k-1]*localVERTs[j-1];
            a[j + k*l] = a[k + j*l];
        }
    }
    delete localVERTs;
    delete localVERTq;
    delete localVERTpp;
    return;
//*-*-                  failure return
L100:
    delete localVERTs;
    delete localVERTq;
    delete localVERTpp;
    ifail = 1;
} /* mnvertLocal */
