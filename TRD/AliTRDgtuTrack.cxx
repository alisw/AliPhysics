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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
//  TRD module global track (GTU)                                            //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TMatrixD.h>
#include <TObjArray.h>

#include "AliLog.h"

#include "AliTRDgeometry.h"
#include "AliTRDcalibDB.h"
#include "AliTRDltuTracklet.h"
#include "AliTRDgtuTrack.h"
#include "Cal/AliTRDCalPIDLQ.h"

ClassImp(AliTRDgtuTrack)

//_____________________________________________________________________________
AliTRDgtuTrack::AliTRDgtuTrack()
  :TObject()
  ,fTracklets(new TObjArray(400))
  ,fYproj(0)
  ,fZproj(0)
  ,fSlope(0)
  ,fDetector(-1)
  ,fNtracklets(0)
  ,fNplanes(0)
  ,fNclusters(0)
  ,fPt(0)
  ,fPhi(0)
  ,fEta(0)
  ,fLabel(-1)
  ,fPID(0)
  ,fIsElectron(kFALSE)
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliTRDgtuTrack::AliTRDgtuTrack(const AliTRDgtuTrack &t)
  :TObject(t)
  ,fTracklets(NULL)
  ,fYproj(t.fYproj)
  ,fZproj(t.fZproj)
  ,fSlope(t.fSlope)
  ,fDetector(t.fDetector)
  ,fNtracklets(t.fNtracklets)
  ,fNplanes(t.fNplanes)
  ,fNclusters(t.fNclusters)
  ,fPt(t.fPt)
  ,fPhi(t.fPhi)
  ,fEta(t.fEta)
  ,fLabel(t.fLabel)
  ,fPID(t.fPID)
  ,fIsElectron(t.fIsElectron)
{
  //
  // Copy contructor
  //

}

//_____________________________________________________________________________
AliTRDgtuTrack &AliTRDgtuTrack::operator=(const AliTRDgtuTrack &t)
{
  //
  // Assignment operator
  //

  if (this != &t) ((AliTRDgtuTrack &) t).Copy(*this); 

  return *this;

} 

//_____________________________________________________________________________
void AliTRDgtuTrack::Copy(TObject &t) const
{
  //
  // Copy function
  //

  ((AliTRDgtuTrack &) t).fTracklets  = NULL;
  ((AliTRDgtuTrack &) t).fYproj      = fYproj;
  ((AliTRDgtuTrack &) t).fZproj      = fZproj;
  ((AliTRDgtuTrack &) t).fSlope      = fSlope;
  ((AliTRDgtuTrack &) t).fDetector   = fDetector;
  ((AliTRDgtuTrack &) t).fNtracklets = fNtracklets;
  ((AliTRDgtuTrack &) t).fNplanes    = fNplanes;
  ((AliTRDgtuTrack &) t).fNclusters  = fNclusters;
  ((AliTRDgtuTrack &) t).fPt         = fPt;
  ((AliTRDgtuTrack &) t).fPhi        = fPhi;
  ((AliTRDgtuTrack &) t).fEta        = fEta;
  ((AliTRDgtuTrack &) t).fLabel      = fLabel;
  ((AliTRDgtuTrack &) t).fPID        = fPID;
  ((AliTRDgtuTrack &) t).fIsElectron = fIsElectron;

}

//_____________________________________________________________________________
AliTRDgtuTrack::~AliTRDgtuTrack()
{
  //
  // Destructor
  //

  if (fTracklets) {
    //fTracklets->Delete();
    fTracklets->Clear();
    delete fTracklets;
    fTracklets = 0;
  }

}

//_____________________________________________________________________________
void AliTRDgtuTrack::AddTracklet(AliTRDltuTracklet *trk)
{
  //
  // Add a LTU tracklet to this track
  //

  Tracklets()->Add(trk);

}

//_____________________________________________________________________________
AliTRDltuTracklet *AliTRDgtuTrack::GetTracklet(Int_t pos) const
{
  //
  // Return LTU tracklet at position "pos"
  //

  if (fTracklets == 0) {
    return 0;
  }
  void *trk = fTracklets->UncheckedAt(pos);
  if (trk == 0) {
    return 0;
  }

  return (AliTRDltuTracklet *) trk;

}

//_____________________________________________________________________________
Int_t AliTRDgtuTrack::Compare(const TObject *o) const
{
  //
  // Compare function for sorting the tracks
  //

  AliTRDgtuTrack *gtutrack = (AliTRDgtuTrack *) o;

  if (fYproj <  gtutrack->GetYproj()) {
    return -1;
  }
  if (fYproj == gtutrack->GetYproj()) {
    return  0;
  }

  return +1;

}

//_____________________________________________________________________________
void AliTRDgtuTrack::Reset()
{
  //
  // Reset the track information
  //

  fYproj      = 0.0;
  fZproj      = 0.0;
  fSlope      = 0.0;
  fDetector   = -1;
  fNtracklets = 0;
  fNplanes    = 0;
  fNclusters  = 0;
  fPt         = 0.0;
  fPhi        = 0.0;
  fEta        = 0.0;
  fLabel      = -1;
  fPID        = 0.0;
  fIsElectron = kFALSE;
  
}

//_____________________________________________________________________________
void AliTRDgtuTrack::Track(Float_t xpl, Float_t field)
{
  //
  // Calculate the kinematics of the found track
  //

  Int_t  i    = 0;
  Int_t  iDet = 0;
  Int_t  nDet = 0;
  Bool_t newDetector;

  AliTRDltuTracklet *trk;
  Int_t nTracklets = GetNtracklets();
  Float_t fC[kNmaxTrk][3];            // X, Y, Z  coordinates of segments

  fYproj      = 0.0;
  fZproj      = 0.0;
  fSlope      = 0.0;
  fNclusters  = 0;
  fNplanes    = 0;
  fNtracklets = GetNtracklets();
  Int_t inDetector[kNplan];
  for (i = 0; i < kNplan; i++) {
    inDetector[i] = -1;
  }

  for (i = 0; i < nTracklets; i++) {

    trk         = GetTracklet(i);
    fYproj     += trk->GetYproj(xpl);
    fZproj     += trk->GetZproj(xpl);
    fSlope     += trk->GetSlope();
    fNclusters += trk->GetNclusters();
    iDet        = trk->GetDetector();

    newDetector = kTRUE;
    for (Int_t id = 0; id < nDet; id++) {
      if (iDet == inDetector[id]) {
	newDetector = kFALSE;
	break;
      }
    }
    if (newDetector) {
      inDetector[nDet++] = iDet;
      fNplanes++;
    }

    fC[i][0] = trk->GetTime0();
    fC[i][1] = trk->GetOffset();
    fC[i][2] = trk->GetRowz();

  }
  fYproj /= (Float_t) nTracklets;
  fZproj /= (Float_t) nTracklets;
  fSlope /= (Float_t) nTracklets;

  Float_t x[kNmaxTrk+1];
  Float_t y[kNmaxTrk+1];
  Float_t z[kNmaxTrk+1];
  Bool_t  count[kNmaxTrk];
  for (i = 0; i < kNmaxTrk; i++) {
    count[i] = kFALSE;
  }

  Int_t iXmin = -1;
  Int_t j     =  0;
  x[0] = y[0] = z[0] = 0.0;
  while (j < nTracklets) {
    iXmin = -1;
    for (i = 0; i < nTracklets; i++) {
      if (count[i]) continue;
      if (iXmin == -1) {
	iXmin = i;
	continue;
      }
      if (fC[i][0] < fC[iXmin][0]) {
        iXmin = i;
      }
    }
    x[j+1] = fC[iXmin][0];
    y[j+1] = fC[iXmin][1];
    z[j+1] = fC[iXmin][2];
    j++;
    count[iXmin] = kTRUE;
  }

  TMatrixD smatrix(2,2);
  TMatrixD sums(2,1);
  TMatrixD res(2,1);
  Double_t xv, yv;
  Float_t a, b;

  smatrix.Zero();
  sums.Zero();
  for (i = 0; i < nTracklets; i++) {
    xv = (Double_t) x[i+1];
    yv = (Double_t) y[i+1];
    smatrix(0,0) += 1.0;
    smatrix(1,1) += xv*xv;
    smatrix(0,1) += xv;
    smatrix(1,0) += xv;
    sums(0,0)    += yv;
    sums(1,0)    += xv*yv;
  }
  res = smatrix.Invert() * sums;
  a   = res(0,0);
  b   = res(1,0);

  Float_t dist  = AliTRDgeometry::GetTime0(1) - AliTRDgeometry::GetTime0(0);
  Float_t fx1   = x[1]          + dist * (Float_t) (nTracklets-1) / 6.0;
  Float_t fy1   = a + b * fx1;
  Float_t fx2   = x[nTracklets] - dist * (Float_t) (nTracklets-1) / 6.0;
  Float_t fy2   = a + b * fx2;
  Float_t d12   = TMath::Sqrt((fx2-fx1)*(fx2-fx1)+(fy2-fy1)*(fy2-fy1));
  Float_t alpha = TMath::ATan(fy2/fx2) - TMath::ATan(fy1/fx1);
  Float_t r     = (d12 / 2.0) / TMath::Sin(alpha);

  fPt = 0.3 * field * 0.01 * r;

  Float_t d1 = fx1*fx1 + fy1*fy1;
  Float_t d2 = fx2*fx2 + fy2*fy2;
  Float_t d  = fx1*fy2 - fx2*fy1;
  
  Float_t xc = (d1*fy2 - d2*fy1) / (2.0*d);
  Float_t yc = (d2*fx1 - d1*fx2) / (2.0*d);

  if (yc != 0.0) {
    fPhi = TMath::ATan(xc/yc);
  } 
  else {
    fPhi = TMath::PiOver2();
  }

  fPhi *= 180.0/TMath::Pi();

  smatrix.Zero();
  sums.Zero();
  for (i = 0; i < nTracklets+1; i++) {
    xv = (Double_t) z[i];
    yv = (Double_t) x[i];
    smatrix(0,0) += 1.0;
    smatrix(1,1) += xv*xv;
    smatrix(0,1) += xv;
    smatrix(1,0) += xv;
    sums(0,0)    += yv;
    sums(1,0)    += xv*yv;
  }
  res = smatrix.Invert() * sums;
  a   = res(0,0);
  b   = res(1,0);
  Float_t theta = TMath::ATan(b);
  
  if (theta < 0.0) {
    theta = TMath::Pi() + theta;
  }
  if (theta == 0.0) {
    fEta = 0.0;
  } 
  else {
    fEta = -TMath::Log(TMath::Tan(theta/2.0));
  }
  
}

//_____________________________________________________________________________
void AliTRDgtuTrack::MakePID()
{
  //
  // Electron likelihood signal
  //

  Int_t i = 0;

  AliTRDcalibDB *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliError("No instance of AliTRDcalibDB.");
    return;  
  }
  const AliTRDCalPIDLQ *pd = calibration->GetPIDLQObject();
  
  AliTRDltuTracklet *trk;
  Int_t   nTracklets = GetNtracklets();
  Int_t   det;
  Int_t   pla;
  Float_t sl;
  Float_t th;
  Float_t q;
  Float_t probPio = 1.0;
  Float_t probEle = 1.0;

  for (i = 0; i < nTracklets; i++) {

    trk = GetTracklet(i);

    sl  = TMath::Abs(trk->GetSlope());     // Tracklet inclination in X-y plane
    th  = trk->GetRowz()/trk->GetTime0();  // Tracklet inclination in X-z plane
    th  = TMath::ATan(TMath::Abs(th));

    q   = trk->GetQ() 
        * TMath::Cos(sl/180.0*TMath::Pi()) 
        * TMath::Cos(th/180.0*TMath::Pi());

    det = trk->GetDetector();
    pla = trk->GetPlane(det);

    // Unclear yet factor to match the PS distributions = 5.8
    // not explained only by the tail filter ...
    // AliRoot v4-03-07 , v4-03-Release
    //q = q * 5.8;

    // HEAD28Mar06
    // Temporary (B. Vulpescu):
    // The charge distributions do not match the new changes in simulation (A. Bercuci),
    // which are nevertheless now in agreement with the beam tests.
    // Some tricks will be used to still have reasonable results
    // To match the existing charge distributions, the charge per layer has to be modified
    // as folows:
    /*
      if (k == 0) {
      // electrons
      q = 4.3 * q + 95.0;
      } else {
      // others
      q = 4.2 * q + 70.0;
      }
    */
    // Since at tracking time we have no information on the particle type, we will modify
    // instead the charge distributions accordingly. This will slow down the sampling.
    // The modified distributions are in TRDdEdxHistogramsV1_BV.root and the CDB has
    // been regenerated with AliTRDCreateDummyCDB.C
    // The new PIDLQ data base has the version :
    // I-AliCDBLocal::Get: CDB object retrieved: path "TRD/Calib/PIDLQ"; run range [0,0];
    // version v0_s1

    // HEAD28Mar06 new distributions: factor = 6.0

    probEle *= pd->GetProbability(0,TMath::Abs(fPt),q*6.0);
    probPio *= pd->GetProbability(2,TMath::Abs(fPt),q*6.0);

  }

  if ((probEle + probPio) > 0.0) {
    fPID = probEle/(probEle + probPio);
  } 
  else {
    fPID = 0.0;
  }

  // Thresholds for LQ cut at 90% electron efficiency (from AliRoot v4-03-07)
  // P [GeV/c]  fPIDcut (between 0 and 1)
  // 2          0.925
  // 3          0.915
  // 4          0.875
  // 5          0.855
  // 6          0.845
  // 8          0.785
  // 10         0.735
  //
  // PIDcut = 0.978 - 0.024 * P[GeV/c]
  //Float_t PIDcut = 0.978 - 0.024 * TMath::Abs(fPt);

  // HEAD28Mar06 with modified distributions (A. Bercuci changes, P. Shukla distributions)
  //Float_t fPIDcut = 0.829 - 0.032 * TMath::Abs(fPt);

  // HEAD28Mar06 with new generated distributions (pol2 fit)
  Float_t absPt   = TMath::Abs(fPt);
  //Float_t fPIDcut = 0.9575 - 0.0832 * absPt + 0.0037 * absPt*absPt;  // 800 bins in dEdx
  Float_t fPIDcut = 0.9482 - 0.0795 * absPt + 0.0035 * absPt*absPt;    // 200 bins in dEdx

  fIsElectron = kFALSE;
  if (fPID >= fPIDcut) {
    fIsElectron = kTRUE;
  }

}

//_____________________________________________________________________________
void AliTRDgtuTrack::CookLabel()
{
  //
  // Cook the track label from tracklets labels
  //

  AliTRDltuTracklet *trk;

  const Int_t kMaxTracks = 10;
  Int_t trackLabel[kMaxTracks];
  Int_t trackCount[kMaxTracks];
  for (Int_t it = 0; it < kMaxTracks; it++) {
    trackLabel[it] = -1;
    trackCount[it] =  0;
  }

  Bool_t counted;
  Int_t  label;
  Int_t  nTracks = 0;
  for (Int_t itrk = 0; itrk < fNtracklets; itrk++) {

    trk = GetTracklet(itrk);

    if (trk->GetLabel() == -1) {
      continue;
    }

    label = trk->GetLabel();

    counted = kFALSE;
    for (Int_t it = 0; it < nTracks; it++) {
      if (label == trackLabel[it]) {
	trackCount[it]++;
	counted = kTRUE;
	break;
      }
    }
    if (!counted) {
      trackLabel[nTracks] = label;
      trackCount[nTracks]++;
      nTracks++;
      if (nTracks == kMaxTracks) {
	Warning("CookLabel","Too many tracks for this tracklet.");
	nTracks--;
	break;
      }
    }
    
  }

  Float_t frac = 4.0 / 5.0;
  for (Int_t it = 0; it < kMaxTracks; it++) {
    if (trackCount[it] >= (Int_t) (frac*fNtracklets)) {
      fLabel = trackLabel[it];
      break;
    }
  }

}

//_____________________________________________________________________________
void AliTRDgtuTrack::ResetTracklets() 
{ 
  //
  // Resets the list of tracklets
  //

  if (fTracklets) {
    //fTracklets->Delete(); 
    fTracklets->Clear();
  }

}

//_____________________________________________________________________________
TObjArray *AliTRDgtuTrack::Tracklets() 
{ 
  //
  // Returns the list of tracklets
  //

  if (!fTracklets) {
    fTracklets = new TObjArray(400); 
  }

  return fTracklets; 

}

//_____________________________________________________________________________
Int_t AliTRDgtuTrack::GetNtracklets() const 
{
  //
  // Returns the number of tracklets
  //

  if (fTracklets) {
    return fTracklets->GetEntriesFast();
  }

  return 0;

}
