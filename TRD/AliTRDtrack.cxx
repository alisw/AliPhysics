/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use	, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

#include <Riostream.h>

#include <TMath.h>
#include <TVector2.h>

#include "AliTracker.h"
#include "AliESDtrack.h"

#include "AliTRDgeometry.h" 
#include "AliTRDcluster.h" 
#include "AliTRDtrack.h"
#include "AliTRDtracklet.h"

// A. Bercuci - used for PID calculations
#include "AliTRDcalibDB.h"
#include "Cal/AliTRDCalPID.h"

ClassImp(AliTRDtrack)

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Represents a reconstructed TRD track                                     //
//  Local TRD Kalman track                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________
AliTRDtrack::AliTRDtrack()
  :AliKalmanTrack()
  ,fSeedLab(-1)
  ,fdEdx(0)
  ,fDE(0)
  ,fClusterOwner(kFALSE)
  ,fStopped(kFALSE)
  ,fLhElectron(0)
  ,fNWrong(0)
  ,fNRotate(0)
  ,fNCross(0)
  ,fNExpected(0)
  ,fNLast(0)
  ,fNExpectedLast(0)
  ,fNdedx(0)
  ,fChi2Last(1e10)
  ,fBackupTrack(0x0)
{
  //
  // AliTRDtrack default constructor
  //

  for (Int_t i = 0; i < kNplane; i++) {
    for (Int_t j = 0; j < kNslice; j++) {
      fdEdxPlane[i][j] = 0.0;
    }
    fTimBinPlane[i] = -1;
    fMom[i]         = -1.;
    fSnp[i]         = 0.;
    fTgl[i]         = 0.;
  }

  for (UInt_t i = 0; i < kMAXCLUSTERSPERTRACK; i++) {
    fIndex[i]       = 0;
    fIndexBackup[i] = 0;
    fdQdl[i]        = 0;
    fClusters[i]    = 0x0;
  }

  for (Int_t i = 0; i < 3; i++) {
    fBudget[i] = 0;
  }

}

//_____________________________________________________________________________
AliTRDtrack::AliTRDtrack(AliTRDcluster *c, Int_t index
                       , const Double_t p[5], const Double_t cov[15] 
                       , Double_t x, Double_t alpha)
  :AliKalmanTrack()
  ,fSeedLab(-1)
  ,fdEdx(0)
  ,fDE(0)
  ,fClusterOwner(kFALSE)
  ,fStopped(kFALSE)
  ,fLhElectron(0)
  ,fNWrong(0)
  ,fNRotate(0)
  ,fNCross(0)
  ,fNExpected(0)
  ,fNLast(0)
  ,fNExpectedLast(0)
  ,fNdedx(0)
  ,fChi2Last(1e10)
  ,fBackupTrack(0x0)
{
  //
  // The main AliTRDtrack constructor.
  //

  Double_t cnv   = 1.0 / (GetBz() * kB2C);

  Double_t pp[5] = { p[0]    
                   , p[1]
                   , x*p[4] - p[2]
                   , p[3]
                   , p[4]*cnv      };

  Double_t c22 = x*x*cov[14] - 2*x*cov[12] + cov[ 5];
  Double_t c32 =   x*cov[13] -     cov[ 8];
  Double_t c20 =   x*cov[10] -     cov[ 3];
  Double_t c21 =   x*cov[11] -     cov[ 4];
  Double_t c42 =   x*cov[14] -     cov[12];

  Double_t cc[15] = { cov[ 0]
                    , cov[ 1],     cov[ 2]
                    , c20,         c21,         c22
                    , cov[ 6],     cov[ 7],     c32,     cov[ 9]
                    , cov[10]*cnv, cov[11]*cnv, c42*cnv, cov[13]*cnv, cov[14]*cnv*cnv };

  Set(x,alpha,pp,cc);
  SetNumberOfClusters(1);
  fIndex[0]    = index;
  fClusters[0] = c;

  for (Int_t i = 0; i < kNplane; i++) {
    for (Int_t j = 0; j < kNslice; j++) {
      fdEdxPlane[i][j] = 0.0;
    }
    fTimBinPlane[i] = -1;
    fMom[i]         = -1.;
    fSnp[i]         =  0.;
    fTgl[i]         =  0.;
  }

  Double_t q = TMath::Abs(c->GetQ());
  Double_t s = GetSnp();
  Double_t t = GetTgl();
  if (s*s < 1) {
    q *= TMath::Sqrt((1-s*s)/(1+t*t));
  }

  fdQdl[0] = q;  	
  for (UInt_t i = 1; i < kMAXCLUSTERSPERTRACK; i++) {
    fdQdl[i]        = 0;
    fIndex[i]       = 0;
    fIndexBackup[i] = 0;
    fClusters[i]    = 0x0;
  }

  for (Int_t i = 0; i < 3;i++) {
    fBudget[i]      = 0;
  }

}                              
           
//_____________________________________________________________________________
AliTRDtrack::AliTRDtrack(const AliTRDtrack &t/*, const Bool_t owner*/)
  :AliKalmanTrack(t) 
  ,fSeedLab(t.GetSeedLabel())
  ,fdEdx(t.fdEdx)
  ,fDE(t.fDE)
  ,fClusterOwner(kTRUE)
  ,fStopped(t.fStopped)
  ,fLhElectron(0)
  ,fNWrong(t.fNWrong)
  ,fNRotate(t.fNRotate)
  ,fNCross(t.fNCross)
  ,fNExpected(t.fNExpected)
  ,fNLast(t.fNLast)
  ,fNExpectedLast(t.fNExpectedLast)
  ,fNdedx(t.fNdedx)
  ,fChi2Last(t.fChi2Last)
  ,fBackupTrack(0x0)
{
  //
  // Copy constructor.
  //

  for (Int_t i = 0; i < kNplane; i++) {
    for (Int_t j = 0; j < kNslice; j++) {
      fdEdxPlane[i][j] = t.fdEdxPlane[i][j];
    }
    fTimBinPlane[i] = t.fTimBinPlane[i];
    fTracklets[i]   = t.fTracklets[i];
    fMom[i]         = t.fMom[i];
    fSnp[i]         = t.fSnp[i];
    fTgl[i]         = t.fTgl[i];
  }

  Int_t n = t.GetNumberOfClusters(); 
  SetNumberOfClusters(n);

  for (Int_t i = 0; i < n; i++) {
    fIndex[i]       = t.fIndex[i];
    fIndexBackup[i] = t.fIndex[i];
    fdQdl[i]        = t.fdQdl[i];
    if (fClusterOwner && t.fClusters[i]) {
      fClusters[i] = new AliTRDcluster(*(t.fClusters[i]));
    }
    else {
      fClusters[i] = t.fClusters[i];
    }
  }

  for (UInt_t i = n; i < kMAXCLUSTERSPERTRACK; i++) {
    fdQdl[i]        = 0;
    fIndex[i]       = 0;
    fIndexBackup[i] = 0;
    fClusters[i]    = 0x0;
  }

  for (Int_t i = 0; i < 3;i++) {
    fBudget[i]      = t.fBudget[i];
  }

}                                

//_____________________________________________________________________________
AliTRDtrack::AliTRDtrack(const AliKalmanTrack &t, Double_t /*alpha*/) 
  :AliKalmanTrack(t) 
  ,fSeedLab(-1)
  ,fdEdx(t.GetPIDsignal())
  ,fDE(0)
  ,fClusterOwner(kFALSE)
  ,fStopped(kFALSE)
  ,fLhElectron(0.0)
  ,fNWrong(0)
  ,fNRotate(0)
  ,fNCross(0)
  ,fNExpected(0)
  ,fNLast(0)
  ,fNExpectedLast(0)
  ,fNdedx(0)
  ,fChi2Last(0.0)
  ,fBackupTrack(0x0)
{
  //
  // Constructor from AliTPCtrack or AliITStrack
  //

  SetLabel(t.GetLabel());
  SetChi2(0.0);
  SetMass(t.GetMass());
  SetNumberOfClusters(0);

  for (Int_t i = 0; i < kNplane; i++) {
    for (Int_t j = 0; j < kNslice; j++) {
      fdEdxPlane[i][j] = 0.0;
    }
    fTimBinPlane[i] = -1;
    fMom[i]         = -1.;
    fSnp[i]         =  0.;
    fTgl[i]         =  0.;
  }

  for (UInt_t i = 0; i < kMAXCLUSTERSPERTRACK; i++) {
    fdQdl[i]        = 0;
    fIndex[i]       = 0;
    fIndexBackup[i] = 0;
    fClusters[i]    = 0x0;
  }
  
  for (Int_t i = 0; i < 3; i++) { 
    fBudget[i]      = 0;
  }

}              

//_____________________________________________________________________________
AliTRDtrack::AliTRDtrack(const AliESDtrack &t)
  :AliKalmanTrack()
  ,fSeedLab(-1)
  ,fdEdx(t.GetTRDsignal())
  ,fDE(0)
  ,fClusterOwner(kFALSE)
  ,fStopped(kFALSE)
  ,fLhElectron(0)
  ,fNWrong(0)
  ,fNRotate(0)
  ,fNCross(0)
  ,fNExpected(0)
  ,fNLast(0)
  ,fNExpectedLast(0)
  ,fNdedx(0)
  ,fChi2Last(1e10)
  ,fBackupTrack(0x0)
{
  //
  // Constructor from AliESDtrack
  //

  SetLabel(t.GetLabel());
  SetChi2(0.0);
  SetMass(t.GetMass());
  SetNumberOfClusters(t.GetTRDclusters(fIndex)); 

  Int_t ncl = t.GetTRDclusters(fIndexBackup);
  for (UInt_t i = ncl; i < kMAXCLUSTERSPERTRACK; i++) {
    fIndexBackup[i] = 0;
    fIndex[i]       = 0;
  }

  for (Int_t i = 0; i < kNplane; i++) {
    for (Int_t j = 0; j < kNslice; j++) {
      fdEdxPlane[i][j] = t.GetTRDsignals(i,j);
    }
    fTimBinPlane[i] = t.GetTRDTimBin(i);
    fMom[i]         = -1.;
    fSnp[i]         =  0.;
    fTgl[i]         =  0.;
  }

  const AliExternalTrackParam *par = &t;
  if (t.GetStatus() & AliESDtrack::kTRDbackup) { 
    par = t.GetOuterParam();
    if (!par) {
      AliError("No backup info!"); 
      par = &t;
    }
  }
  Set(par->GetX(),par->GetAlpha(),par->GetParameter(),par->GetCovariance());

  for (UInt_t i = 0; i < kMAXCLUSTERSPERTRACK; i++) {
    fdQdl[i]     = 0;
    fClusters[i] = 0x0;
  }

  for (Int_t i = 0; i < 3; i++) {
    fBudget[i] = 0;
  }

  if ((t.GetStatus() & AliESDtrack::kTIME) == 0) {
    return;
  }

  StartTimeIntegral();
  Double_t times[10]; 
  t.GetIntegratedTimes(times); 
  SetIntegratedTimes(times);
  SetIntegratedLength(t.GetIntegratedLength());

}  

//____________________________________________________________________________
AliTRDtrack::~AliTRDtrack()
{
  //
  // Destructor
  //

  if (fBackupTrack) {
    delete fBackupTrack;
  }
  fBackupTrack = 0x0;
  if (fClusterOwner) {
    UInt_t icluster = 0;
    while ((icluster < kMAXCLUSTERSPERTRACK) && fClusters[icluster]) {
      delete fClusters[icluster];
      fClusters[icluster] = 0x0;
      icluster++;
    }
  }

}

//____________________________________________________________________________
Float_t AliTRDtrack::StatusForTOF()
{
  //
  // Defines the status of the TOF extrapolation
  //

  // Definition of res ????
  Float_t res = (0.2 + 0.8 * (fN / (fNExpected + 5.0)))
                     * (0.4 + 0.6 * fTracklets[5].GetN() / 20.0);
  res *= (0.25 + 0.8 * 40.0 / (40.0 + fBudget[2]));
  return res;

  // This part of the function is never reached ????
  // What defines these parameters ????
  //Int_t status = 0;
  //if (GetNumberOfClusters()                <  20)  return 0;
  //if ((fN                                  > 110) && 
  //    (fChi2/(Float_t(fN))                 <   3)) return 3; // Gold
  //if ((fNLast                              >  30) &&
  //    (fChi2Last/(Float_t(fNLast))         <   3)) return 3; // Gold
  //if ((fNLast                              >  20) &&
  //    (fChi2Last/(Float_t(fNLast))         <   2)) return 3; // Gold
  //if ((fNLast/(fNExpectedLast+3.0)         > 0.8) && 
  //    (fChi2Last/Float_t(fNLast)           <   5) && 
  //    (fNLast                              >  20)) return 2; // Silver
  //if ((fNLast                              >   5) &&
  //    (((fNLast+1.0)/(fNExpectedLast+1.0)) > 0.8) && 
  //    (fChi2Last/(fNLast-5.0)              <   6)) return 1; 
  //
  //return status;

}

//_____________________________________________________________________________
Int_t AliTRDtrack::Compare(const TObject *o) const 
{
  //
  // Compares tracks according to their Y2 or curvature
  //

  AliTRDtrack *t = (AliTRDtrack *) o;

  Double_t co = TMath::Abs(t->GetC());
  Double_t c  = TMath::Abs(GetC());  

  if      (c > co) {
    return 1;
  }
  else if (c < co) {
    return -1;
  }
  return 0;

}                

//_____________________________________________________________________________
void AliTRDtrack::CookdEdx(Double_t low, Double_t up) 
{
  //
  // Calculates the truncated dE/dx within the "low" and "up" cuts.
  //

  // Array to sort the dEdx values according to amplitude
  Float_t sorted[kMAXCLUSTERSPERTRACK];
  fdEdx = 0.;
	
  // Require at least 10 clusters for a dedx measurement
  if (fNdedx < 10) return;

  // Can fdQdl be negative ????
  for (Int_t i = 0; i < fNdedx; i++) {
    sorted[i] = TMath::Abs(fdQdl[i]);
  }
  // Sort the dedx values by amplitude
  Int_t *index = new Int_t[fNdedx];
  TMath::Sort(fNdedx, sorted, index, kFALSE);

  // Sum up the truncated charge between lower and upper bounds 
  Int_t nl = Int_t(low * fNdedx);
  Int_t nu = Int_t( up * fNdedx);
  for (Int_t i = nl; i <= nu; i++) {
    fdEdx += sorted[index[i]];
  }
  fdEdx /= (nu - nl + 1.0);

  delete[] index;

}                     

//_____________________________________________________________________________
void AliTRDtrack::CookdEdxTimBin(const Int_t/* tid*/)
{
  //
  // Set fdEdxPlane and fTimBinPlane and also get the 
  // Time bin for Max. Cluster
  //
  // Authors:
  // Prashant Shukla (shukla@physi.uni-heidelberg.de)
  // Alexandru Bercuci (A.Bercuci@gsi.de)
  //

  Double_t  maxcharge[AliESDtrack::kNPlane]; // max charge in chamber
  // Number of clusters attached to track per chamber and slice
  Int_t     nCluster[AliESDtrack::kNPlane][AliESDtrack::kNSlice];
  // Number of time bins in chamber
  Int_t ntb = AliTRDcalibDB::Instance()->GetNumberOfTimeBins();
  Int_t plane;                  // Plane of current cluster
  Int_t tb;                     // Time bin of current cluster
  Int_t slice;                  // Current slice
  AliTRDcluster *cluster = 0x0; // Pointer to current cluster

  // Reset class and local counters/variables
  for (Int_t iPlane = 0; iPlane < AliESDtrack::kNPlane; iPlane++) {
    fTimBinPlane[iPlane] = -1;
    maxcharge[iPlane]    = 0.;
    for (Int_t iSlice = 0; iSlice < AliESDtrack::kNSlice; iSlice++) {
      fdEdxPlane[iPlane][iSlice] = 0.;
      nCluster[iPlane][iSlice]   = 0;
    }
  }

  // Start looping over clusters attached to this track
  for (Int_t iClus = 0; iClus < GetNumberOfClusters(); iClus++) {

    cluster = fClusters[iClus]; //(AliTRDcluster*)tracker->GetCluster(fIndex[iClus]);
    if(!cluster) continue;

    // Read info from current cluster
    plane  = AliTRDgeometry::GetPlane(cluster->GetDetector());
    if (plane < 0 || plane >= AliESDtrack::kNPlane) {
      AliError(Form("Wrong plane %d", plane));
      continue;
    }

    tb     = cluster->GetLocalTimeBin();
    if ((tb == 0) || (tb >= ntb)) {
      AliWarning(Form("time bin 0 or > %d in cluster %d", ntb, iClus));
      AliInfo(Form("dQ/dl %f", fdQdl[iClus]));
      continue;
    }
	
    slice = tb * AliESDtrack::kNSlice / ntb;

    fdEdxPlane[plane][slice] += fdQdl[iClus];
    if (fdQdl[iClus] > maxcharge[plane]) {
      maxcharge[plane]    = fdQdl[iClus];
      fTimBinPlane[plane] = tb;
    }

    nCluster[plane][slice]++;

  } // End of loop over cluster
	
  // Normalize fdEdxPlane to number of clusters and set track segments
  for (Int_t iPlane = 0; iPlane < AliESDtrack::kNPlane; iPlane++) {
    for (Int_t iSlice = 0; iSlice < AliESDtrack::kNSlice; iSlice++) {
      if (nCluster[iPlane][iSlice]) {
        fdEdxPlane[iPlane][iSlice] /= nCluster[iPlane][iSlice];
      }
    }
  }

}

//_____________________________________________________________________________
void	AliTRDtrack::CookdEdxNN(Float_t *dedx){

  // This function calcuates the input for the neural networks 
  // which are used for the PID. 

	Int_t ntb = AliTRDcalibDB::Instance()->GetNumberOfTimeBins();//number of time bins in chamber
	Int_t plane;   // plane of current cluster
	Int_t tb;      // time bin of current cluster
	Int_t slice;   // curent slice
	AliTRDcluster *cluster = 0x0; // pointer to current cluster
	const Int_t lMLPscale = 16000; // scaling of the MLP input to be smaller than 1

	// Reset class and local contors/variables
	for (Int_t iPlane = 0; iPlane < AliESDtrack::kNPlane; iPlane++){
	  for (Int_t iSlice = 0; iSlice < kNMLPslice; iSlice++) {
	    *(dedx + (kNMLPslice * iPlane) + iSlice) = 0.;
	  }
	}

	// start looping over clusters attached to this track
	for (Int_t iClus = 0; iClus < GetNumberOfClusters(); iClus++) {
	  cluster = fClusters[iClus]; //(AliTRDcluster*)tracker->GetCluster(fIndex[iClus]);
	  if(!cluster) continue;
	  
	  // Read info from current cluster
	  plane  = AliTRDgeometry::GetPlane(cluster->GetDetector());
	  if (plane < 0 || plane >= AliESDtrack::kNPlane) {
	    AliError(Form("Wrong plane %d", plane));
	    continue;
	  }

	  tb     = cluster->GetLocalTimeBin();
	  if(tb == 0 || tb >= ntb){
	    AliWarning(Form("time bin 0 or > %d in cluster %d", ntb, iClus));
	    AliInfo(Form("dQ/dl %f", fdQdl[iClus]));
	    continue;
	  }

	  slice   = tb * kNMLPslice / ntb;
	  
	  *(dedx+(kNMLPslice * plane) + slice) += fdQdl[iClus]/lMLPscale;
	} // End of loop over cluster
	 
}


//_____________________________________________________________________________
void	AliTRDtrack::SetTrackSegmentDirMom(const Int_t plane)
{
  //
  // Set the momenta for a track segment in a given plane
  //

  if ((plane <       0) || 
      (plane>= kNplane)) {
    AliError(Form("Trying to access out of range plane (%d)", plane));
    return;
  }
	
  fSnp[plane] = GetSnp();
  fTgl[plane] = GetTgl();
  Double_t p[3]; 
  GetPxPyPz(p);
  fMom[plane] = TMath::Sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);

}

//_____________________________________________________________________________
Float_t	AliTRDtrack::GetTrackLengthPlane(Int_t plane) const
{
  //
  // Calculate the track length of a track segment in a given plane
  //

  if ((plane < 0) || (plane >= kNplane)) return 0.;

  return (AliTRDgeometry::AmThick() + AliTRDgeometry::DrThick())
        / TMath::Sqrt((1.0 - fSnp[plane]*fSnp[plane]) 
                    / (1.0 + fTgl[plane]*fTgl[plane]));

}

//_____________________________________________________________________________
Bool_t AliTRDtrack::CookPID(Int_t &pidQuality)
{
  //
  // This function calculates the PID probabilities based on TRD signals
  //
  // The method produces probabilities based on the charge
  // and the position of the maximum time bin in each layer.
  // The dE/dx information can be used as global charge or 2 to 3
  // slices. Check AliTRDCalPID and AliTRDCalPIDRefMaker for the actual
  // implementation.
  //
  // Author
  // Alex Bercuci (A.Bercuci@gsi.de) 2nd May 2007
  //

  AliTRDcalibDB *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliError("No access to calibration data");
    return kFALSE;
  }
	
  // Retrieve the CDB container class with the probability distributions
  const AliTRDCalPID *pd = calibration->GetPIDObject(fPIDmethod == kNN ? 0 : 1);
  if (!pd) {
    AliError("No access to AliTRDCalPID");
    return kFALSE;
  }

  // Calculate the input for the NN if fPIDmethod is kNN
  Float_t ldEdxNN[AliTRDCalPID::kNPlane * kNMLPslice], *dedx = 0x0;
  if(fPIDmethod == kNN) {
    CookdEdxNN(&ldEdxNN[0]);
  }

  // Sets the a priori probabilities
  for(int ispec=0; ispec<AliPID::kSPECIES; ispec++) {
    fPID[ispec] = 1./AliPID::kSPECIES;	
  }


  if(AliPID::kSPECIES != 5){
    AliError("Probabilities array defined only for 5 values. Please modify !!");
    return kFALSE;
  }


  pidQuality = 0;
  Float_t length;  // track segment length in chamber

  // Skip tracks which have no TRD signal at all
  if (fdEdx == 0.) return kFALSE;
	
  for (Int_t iPlane = 0; iPlane < AliTRDgeometry::kNplan; iPlane++) {

    length = (AliTRDgeometry::AmThick() + AliTRDgeometry::DrThick())/TMath::Sqrt((1. - fSnp[iPlane]*fSnp[iPlane]) / (1. + fTgl[iPlane]*fTgl[iPlane]));

    // check data
    if((fdEdxPlane[iPlane][0] + fdEdxPlane[iPlane][1] + fdEdxPlane[iPlane][2]) <=  0.
       || fTimBinPlane[iPlane] == -1.) continue;

    // this track segment has fulfilled all requierments
    pidQuality++;

    // Get the probabilities for the different particle species
    for (Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++) {
      switch(fPIDmethod){
      case kLQ:
	dedx = fdEdxPlane[iPlane];
	break;
      case kNN:
	dedx = &ldEdxNN[iPlane*kNMLPslice];
	break;
      }
      fPID[iSpecies] *= pd->GetProbability(iSpecies, fMom[iPlane], dedx, length, iPlane);
    }
  }
  if(pidQuality == 0) return kTRUE;

  // normalize probabilities
  Double_t probTotal = 0.;
  for (Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++) probTotal += fPID[iSpecies];

  if(probTotal <= 0.) {
    AliWarning("The total probability over all species <= 0. This may be caused by some error in the reference histograms.");
    return kFALSE;
  }
  for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++) fPID[iSpecies] /= probTotal;
	
  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliTRDtrack::PropagateTo(Double_t xk, Double_t xx0, Double_t xrho)
{
  //
  // Propagates this track to a reference plane defined by "xk" [cm] 
  // correcting for the mean crossed material.
  //
  // "xx0"  - thickness/rad.length [units of the radiation length] 
  // "xrho" - thickness*density    [g/cm^2] 
  // 

  if (xk == GetX()) {
    return kTRUE;
  }

  Double_t oldX = GetX();
  Double_t oldY = GetY();
  Double_t oldZ = GetZ();

  Double_t bz   = GetBz();

  if (!AliExternalTrackParam::PropagateTo(xk,bz)) {
    return kFALSE;
  }

  Double_t x = GetX();
  Double_t y = GetY();
  Double_t z = GetZ();

  if (oldX < xk) {
    xrho = -xrho;
    if (IsStartedTimeIntegral()) {
      Double_t l2  = TMath::Sqrt((x-oldX)*(x-oldX) + (y-oldY)*(y-oldY) + (z-oldZ)*(z-oldZ));
      Double_t crv = GetC();
      if (TMath::Abs(l2*crv) > 0.0001) {
        // Make correction for curvature if neccesary
        l2 = 0.5 * TMath::Sqrt((x-oldX)*(x-oldX) + (y-oldY)*(y-oldY));
        l2 = 2.0 * TMath::ASin(l2 * crv) / crv;
        l2 = TMath::Sqrt(l2*l2 + (z-oldZ)*(z-oldZ));
      }
      AddTimeStep(l2);
    }
  }

  if (!AliExternalTrackParam::CorrectForMeanMaterial(xx0,xrho,GetMass())) { 
    return kFALSE;
  }

  {

    // Energy losses************************
    Double_t p2    = (1.0 + GetTgl()*GetTgl()) / (GetSigned1Pt()*GetSigned1Pt());
    Double_t beta2 = p2 / (p2 + GetMass()*GetMass());
    if (beta2<1.0e-10 || (5940.0 * beta2/(1.0 - beta2 + 1.0e-10) - beta2) < 0.0) {
      return kFALSE;
    }

    Double_t dE    = 0.153e-3 / beta2 
                   * (log(5940.0 * beta2/(1.0 - beta2 + 1.0e-10)) - beta2)
                   * xrho;
    fBudget[0] += xrho;

    /*
    // Suspicious part - think about it ?
    Double_t kinE =  TMath::Sqrt(p2);
    if (dE > 0.8*kinE) dE = 0.8 * kinE;  //      
    if (dE < 0)        dE = 0.0;         // Not valid region for Bethe bloch 
    */
 
    fDE += dE;

    /*
    // Suspicious ! I.B.
    Double_t sigmade = 0.07 * TMath::Sqrt(TMath::Abs(dE));   // Energy loss fluctuation 
    Double_t sigmac2 = sigmade*sigmade*fC*fC*(p2+GetMass()*GetMass())/(p2*p2);
    fCcc += sigmac2;
    fCee += fX*fX * sigmac2;  
    */

  }

  return kTRUE;            

}     

//_____________________________________________________________________________
Bool_t AliTRDtrack::Update(const AliTRDcluster *c, Double_t chisq, Int_t index
                         , Double_t h01)
{
  //
  // Assignes found cluster to the track and updates track information
  //

  Bool_t fNoTilt = kTRUE;
  if (TMath::Abs(h01) > 0.003) {
    fNoTilt = kFALSE;
  }

  // Add angular effect to the error contribution -  MI
  Float_t tangent2 = GetSnp()*GetSnp();
  if (tangent2 < 0.90000) {
    tangent2 = tangent2 / (1.0 - tangent2);
  }
  //Float_t errang = tangent2 * 0.04;

  Double_t p[2]   = {c->GetY(), c->GetZ() };
  //Double_t cov[3] = {c->GetSigmaY2()+errang, 0.0, c->GetSigmaZ2()*100.0 };
  Double_t sy2    = c->GetSigmaY2() * 4.0;
  Double_t sz2    = c->GetSigmaZ2() * 4.0;
  Double_t cov[3] = { sy2 + h01*h01*sz2, h01*(sy2-sz2), sz2 + h01*h01*sy2 };

  if (!AliExternalTrackParam::Update(p,cov)) {
    return kFALSE;
  }

  Int_t n   = GetNumberOfClusters();
  fIndex[n] = index;
  SetNumberOfClusters(n+1);

  SetChi2(GetChi2()+chisq);

  return kTRUE;

}        

//_____________________________________________________________________________
Int_t AliTRDtrack::UpdateMI(AliTRDcluster *c, Double_t chisq, Int_t index
                          , Double_t h01, Int_t /*plane*/, Int_t tid)
{
  //
  // Assignes found cluster to the track and updates track information
  //

  Bool_t fNoTilt = kTRUE;
  if (TMath::Abs(h01) > 0.003) {
    fNoTilt = kFALSE;
  }

  // Add angular effect to the error contribution and make correction  -  MI
  Double_t tangent2 = GetSnp()*GetSnp();
  if (tangent2 < 0.90000){
    tangent2 = tangent2 / (1.0-tangent2);
  }
  Double_t tangent = TMath::Sqrt(tangent2);
  if (GetSnp() < 0) {
    tangent *= -1;
  }

  //
  // Is the following still needed ????
  //

  //  Double_t correction = 0*plane;
  /*
  Double_t errang = tangent2*0.04;  //
  Double_t errsys =0.025*0.025*20;  //systematic error part 

  Float_t extend =1;
  if (c->GetNPads()==4) extend=2;
  */
  //if (c->GetNPads()==5)  extend=3;
  //if (c->GetNPads()==6)  extend=3;
  //if (c->GetQ()<15) return 1;

  /*
  if (corrector!=0){
  //if (0){
    correction = corrector->GetCorrection(plane,c->GetLocalTimeBin(),tangent);
    if (TMath::Abs(correction)>0){
      //if we have info 
      errang     = corrector->GetSigma(plane,c->GetLocalTimeBin(),tangent);
      errang    *= errang;      
      errang    += tangent2*0.04;
    }
  }
  */
  //
  //Double_t padlength = TMath::Sqrt(c->GetSigmaZ2()*12.);
  /*
    {
      Double_t dy=c->GetY() - GetY(), dz=c->GetZ() - GetZ();     
      printf("%e %e %e %e\n",dy,dz,padlength/2,h01);
    }
  */

  Double_t p[2]   = { c->GetY(), c->GetZ() };
  //Double_t cov[3]={ (c->GetSigmaY2()+errang+errsys)*extend, 0.0, c->GetSigmaZ2()*10000.0 };
  Double_t sy2    = c->GetSigmaY2() * 4.0;
  Double_t sz2    = c->GetSigmaZ2() * 4.0;
  Double_t cov[3] = { sy2 + h01*h01*sz2, h01*(sy2-sz2), sz2 + h01*h01*sy2 };

  if (!AliExternalTrackParam::Update(p,cov)) {
    return kFALSE;
  }

  // Register cluster to track
  Int_t n      = GetNumberOfClusters();
  fIndex[n]    = index;
  fClusters[n] = c;
  SetNumberOfClusters(n+1);
  SetChi2(GetChi2() + chisq);

  return kTRUE;

}                     

// //_____________________________________________________________________________
// Int_t AliTRDtrack::UpdateMI(const AliTRDtracklet &tracklet)
// {
//   //
//   // Assignes found tracklet to the track and updates track information
//   //
//   // Can this be removed ????
//   //
//
//   Double_t r00=(tracklet.GetTrackletSigma2()), r01=0., r11= 10000.;
//   r00+=fCyy; r01+=fCzy; r11+=fCzz;
//   //
//   Double_t det=r00*r11 - r01*r01;
//   Double_t tmp=r00; r00=r11/det; r11=tmp/det; r01=-r01/det;
//   //

//   Double_t dy=tracklet.GetY() - fY, dz=tracklet.GetZ() - fZ;

  
//   Double_t s00 = tracklet.GetTrackletSigma2();  // error pad
//   Double_t s11 = 100000;   // error pad-row
//   Float_t  h01 = tracklet.GetTilt();
//   //
//   //  r00 = fCyy + 2*fCzy*h01 + fCzz*h01*h01+s00;
//   r00 = fCyy + fCzz*h01*h01+s00;
//   //  r01 = fCzy + fCzz*h01;
//   r01 = fCzy ;
//   r11 = fCzz + s11;
//   det = r00*r11 - r01*r01;
//   // inverse matrix
//   tmp=r00; r00=r11/det; r11=tmp/det; r01=-r01/det;

//   Double_t k00=fCyy*r00+fCzy*r01, k01=fCyy*r01+fCzy*r11;
//   Double_t k10=fCzy*r00+fCzz*r01, k11=fCzy*r01+fCzz*r11;
//   Double_t k20=fCey*r00+fCez*r01, k21=fCey*r01+fCez*r11;
//   Double_t k30=fCty*r00+fCtz*r01, k31=fCty*r01+fCtz*r11;
//   Double_t k40=fCcy*r00+fCcz*r01, k41=fCcy*r01+fCcz*r11;
  
//   // K matrix
// //   k00=fCyy*r00+fCzy*(r01+h01*r00),k01=fCyy*r01+fCzy*(r11+h01*r01);
// //   k10=fCzy*r00+fCzz*(r01+h01*r00),k11=fCzy*r01+fCzz*(r11+h01*r01);
// //   k20=fCey*r00+fCez*(r01+h01*r00),k21=fCey*r01+fCez*(r11+h01*r01);
// //   k30=fCty*r00+fCtz*(r01+h01*r00),k31=fCty*r01+fCtz*(r11+h01*r01);
// //   k40=fCcy*r00+fCcz*(r01+h01*r00),k41=fCcy*r01+fCcz*(r11+h01*r01);  
//   //
//   //Update measurement
//   Double_t cur=fC + k40*dy + k41*dz, eta=fE + k20*dy + k21*dz;  
//   //  cur=fC + k40*dy + k41*dz; eta=fE + k20*dy + k21*dz;
//   if (TMath::Abs(cur*fX-eta) >= 0.90000) {
//     //Int_t n=GetNumberOfClusters();
//     //      if (n>4) cerr<<n<<" AliTRDtrack warning: Filtering failed !\n";
//     return 0;
//   }                           
// //   k01+=h01*k00;
// //   k11+=h01*k10;
// //   k21+=h01*k20;
// //   k31+=h01*k30;
// //   k41+=h01*k40;  


//   fY += k00*dy + k01*dz;
//   fZ += k10*dy + k11*dz;
//   fE  = eta;
//   fT += k30*dy + k31*dz;
//   fC  = cur;
    
  
//   //Update covariance
//   //
//   //
//   Double_t oldyy = fCyy, oldzz = fCzz; //, oldee=fCee, oldcc =fCcc;
//   Double_t oldzy = fCzy, oldey = fCey, oldty=fCty, oldcy =fCcy;
//   Double_t oldez = fCez, oldtz = fCtz, oldcz=fCcz;
//   //Double_t oldte = fCte, oldce = fCce;
//   //Double_t oldct = fCct;

//   fCyy-=k00*oldyy+k01*oldzy;   
//   fCzy-=k10*oldyy+k11*oldzy;
//   fCey-=k20*oldyy+k21*oldzy;   
//   fCty-=k30*oldyy+k31*oldzy;
//   fCcy-=k40*oldyy+k41*oldzy;  
//   //
//   fCzz-=k10*oldzy+k11*oldzz;
//   fCez-=k20*oldzy+k21*oldzz;   
//   fCtz-=k30*oldzy+k31*oldzz;
//   fCcz-=k40*oldzy+k41*oldzz;
//   //
//   fCee-=k20*oldey+k21*oldez;   
//   fCte-=k30*oldey+k31*oldez;
//   fCce-=k40*oldey+k41*oldez;
//   //
//   fCtt-=k30*oldty+k31*oldtz;
//   fCct-=k40*oldty+k41*oldtz;
//   //

//   fCcc-=k40*oldcy+k41*oldcz;                 
//   //
  
//   //Int_t n=GetNumberOfClusters();
//   //fIndex[n]=index;
//   //SetNumberOfClusters(n+1);

//   //SetChi2(GetChi2()+chisq);
//   //  cerr<<"in update: fIndex["<<fN<<"] = "<<index<<endl;

//   return 1;      

// }                     

//_____________________________________________________________________________
Bool_t AliTRDtrack::Rotate(Double_t alpha, Bool_t absolute)
{
  //
  // Rotates track parameters in R*phi plane
  // if absolute rotation alpha is in global system
  // otherwise alpha rotation is relative to the current rotation angle
  //  

  if (absolute) {
    alpha -= GetAlpha();
  }
  else{
    fNRotate++;
  }

  return AliExternalTrackParam::Rotate(GetAlpha()+alpha);

}                         

//_____________________________________________________________________________
Double_t AliTRDtrack::GetPredictedChi2(const AliTRDcluster *c, Double_t h01) const
{
  //
  // Returns the track chi2
  //  

  Double_t p[2]   = { c->GetY(), c->GetZ() };
  Double_t sy2    = c->GetSigmaY2() * 4.0;
  Double_t sz2    = c->GetSigmaZ2() * 4.0;
  Double_t cov[3] = { sy2 + h01*h01*sz2, h01*(sy2-sz2), sz2 + h01*h01*sy2 };

  return AliExternalTrackParam::GetPredictedChi2(p,cov);

  //
  // Can the following be removed ????
  //
  /*
  Bool_t fNoTilt = kTRUE;
  if(TMath::Abs(h01) > 0.003) fNoTilt = kFALSE;

  return (c->GetY() - GetY())*(c->GetY() - GetY())/c->GetSigmaY2();
  */

  /*
  Double_t chi2, dy, r00, r01, r11;

  if(fNoTilt) {
    dy=c->GetY() - fY;
    r00=c->GetSigmaY2();    
    chi2 = (dy*dy)/r00;    
  }
  else {
    Double_t padlength = TMath::Sqrt(c->GetSigmaZ2()*12);
    //
    r00=c->GetSigmaY2(); r01=0.; r11=c->GetSigmaZ2();
    r00+=fCyy; r01+=fCzy; r11+=fCzz;

    Double_t det=r00*r11 - r01*r01;
    if (TMath::Abs(det) < 1.e-10) {
      Int_t n=GetNumberOfClusters(); 
      if (n>4) cerr<<n<<" AliTRDtrack warning: Singular matrix !\n";
      return 1e10;
    }
    Double_t tmp=r00; r00=r11; r11=tmp; r01=-r01;
    Double_t dy=c->GetY() - fY, dz=c->GetZ() - fZ;
    Double_t tiltdz = dz;
    if (TMath::Abs(tiltdz)>padlength/2.) {
      tiltdz = TMath::Sign(padlength/2,dz);
    }
    //    dy=dy+h01*dz;
    dy=dy+h01*tiltdz;

    chi2 = (dy*r00*dy + 2*r01*dy*dz + dz*r11*dz)/det; 
  }

  return chi2;
  */

}

//_____________________________________________________________________________
void AliTRDtrack::MakeBackupTrack()
{
  //
  // Creates a backup track
  //

  if (fBackupTrack) {
    delete fBackupTrack;
  }
  fBackupTrack = new AliTRDtrack(*this);
  
}

//_____________________________________________________________________________
Int_t AliTRDtrack::GetProlongation(Double_t xk, Double_t &y, Double_t &z)
{
  //
  // Find a prolongation at given x
  // Return 0 if it does not exist
  //  

  Double_t bz = GetBz();

  if (!AliExternalTrackParam::GetYAt(xk,bz,y)) {
    return 0;
  }
  if (!AliExternalTrackParam::GetZAt(xk,bz,z)) {
    return 0;
  }

  return 1;  

}

//_____________________________________________________________________________
Int_t AliTRDtrack::PropagateToX(Double_t xr, Double_t step)
{
  //
  // Propagate track to given x  position 
  // Works inside of the 20 degree segmentation
  // (local cooordinate frame for TRD , TPC, TOF)
  // 
  // Material budget from geo manager
  // 

  Double_t xyz0[3];
  Double_t xyz1[3];
  Double_t y;
  Double_t z;

  const Double_t kAlphac  = TMath::Pi()/9.0;   
  const Double_t kTalphac = TMath::Tan(kAlphac*0.5);

  // Critical alpha  - cross sector indication
  Double_t dir = (GetX() > xr) ? -1.0 : 1.0;

  // Direction +-
  for (Double_t x = GetX()+dir*step; dir*x < dir*xr; x += dir*step) {

    GetXYZ(xyz0);	
    GetProlongation(x,y,z);
    xyz1[0] = x * TMath::Cos(GetAlpha()) + y * TMath::Sin(GetAlpha()); 
    xyz1[1] = x * TMath::Sin(GetAlpha()) - y * TMath::Cos(GetAlpha());
    xyz1[2] = z;
    Double_t param[7];
    AliTracker::MeanMaterialBudget(xyz0,xyz1,param);

    if ((param[0] > 0) &&
        (param[1] > 0)) {
      PropagateTo(x,param[1],param[0]*param[4]);
    }

    if (GetY() >  GetX()*kTalphac) {
      Rotate(-kAlphac);
    }
    if (GetY() < -GetX()*kTalphac) {
      Rotate( kAlphac);
    }

  }

  PropagateTo(xr);

  return 0;

}

//_____________________________________________________________________________
Int_t   AliTRDtrack::PropagateToR(Double_t r,Double_t step)
{
  //
  // Propagate track to the radial position
  // Rotation always connected to the last track position
  //

  Double_t xyz0[3];
  Double_t xyz1[3];
  Double_t y;
  Double_t z; 

  Double_t radius = TMath::Sqrt(GetX()*GetX() + GetY()*GetY());
  // Direction +-
  Double_t dir    = (radius > r) ? -1.0 : 1.0;   

  for (Double_t x = radius+dir*step; dir*x < dir*r; x += dir*step) {

    GetXYZ(xyz0);	
    Double_t alpha = TMath::ATan2(xyz0[1],xyz0[0]);
    Rotate(alpha,kTRUE);
    GetXYZ(xyz0);	
    GetProlongation(x,y,z);
    xyz1[0] = x * TMath::Cos(alpha) + y * TMath::Sin(alpha); 
    xyz1[1] = x * TMath::Sin(alpha) - y * TMath::Cos(alpha);
    xyz1[2] = z;
    Double_t param[7];
    AliTracker::MeanMaterialBudget(xyz0,xyz1,param);
    if (param[1] <= 0) {
      param[1] = 100000000;
    }
    PropagateTo(x,param[1],param[0]*param[4]);

  } 

  GetXYZ(xyz0);	
  Double_t alpha = TMath::ATan2(xyz0[1],xyz0[0]);
  Rotate(alpha,kTRUE);
  GetXYZ(xyz0);	
  GetProlongation(r,y,z);
  xyz1[0] = r * TMath::Cos(alpha) + y * TMath::Sin(alpha); 
  xyz1[1] = r * TMath::Sin(alpha) - y * TMath::Cos(alpha);
  xyz1[2] = z;
  Double_t param[7];
  AliTracker::MeanMaterialBudget(xyz0,xyz1,param);

  if (param[1] <= 0) {
    param[1] = 100000000;
  }
  PropagateTo(r,param[1],param[0]*param[4]);

  return 0;

}

//_____________________________________________________________________________
Int_t AliTRDtrack::GetSector() const
{
  //
  // Return the current sector
  //

  return Int_t(TVector2::Phi_0_2pi(GetAlpha()) / AliTRDgeometry::GetAlpha())
             % AliTRDgeometry::kNsect;

}

//_____________________________________________________________________________
void AliTRDtrack::SetSampledEdx(Float_t q, Int_t i)    
{
  //
  // The sampled energy loss
  //

  Double_t s = GetSnp();
  Double_t t = GetTgl();
  q *= TMath::Sqrt((1.0 - s*s) / (1.0 + t*t));
  fdQdl[i] = q;

}     

//_____________________________________________________________________________
void AliTRDtrack::SetSampledEdx(Float_t q) 
{
  //
  // The sampled energy loss
  //

  Double_t s = GetSnp();
  Double_t t = GetTgl();
  q *= TMath::Sqrt((1.0 - s*s) / (1.0 + t*t));
  fdQdl[fNdedx] = q;
  fNdedx++;

}     

//_____________________________________________________________________________
Double_t AliTRDtrack::GetBz() const 
{
  //
  // Returns Bz component of the magnetic field (kG)
  //

  if (AliTracker::UniformField()) {
    return AliTracker::GetBz();
  }
  Double_t r[3]; 
  GetXYZ(r);

  return AliTracker::GetBz(r);

}
