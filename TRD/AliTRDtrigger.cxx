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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD trigger class                                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TTree.h>
#include <TBranch.h>
#include <TMatrixD.h>

#include "AliLog.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliHeader.h"
#include "AliRawReader.h"

#include "AliTRDdigitsManager.h"
#include "AliTRDgeometry.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCommonParam.h"
#include "Cal/AliTRDCalPIDLQ.h"
#include "AliTRDrawData.h"

#include "AliTRDtrigger.h"
#include "AliTRDmcmTracklet.h"
#include "AliTRDtrigParam.h"
#include "AliTRDmcm.h"

ClassImp(AliTRDtrigger)
ClassImp(AliTRDltuTracklet)
ClassImp(AliTRDgtuTrack)
ClassImp(AliTRDmodule)

//_____________________________________________________________________________
AliTRDltuTracklet::AliTRDltuTracklet(Int_t det, 
				     Int_t row, 
				     Float_t rowz,
				     Float_t slope, 
				     Float_t offset, 
				     Float_t time, 
				     Int_t ncl,
				     Int_t label,
				     Float_t q) 
{
  //
  // AliTRDltuTracklet constructor
  //

  fDetector  = det;
  fRow       = row;
  fRowz      = rowz;
  fSlope     = slope;
  fX         = time;
  fY         = offset;
  fNclusters = ncl;
  fLabel     = label;
  fQ         = q;

}

//_____________________________________________________________________________
Float_t AliTRDltuTracklet::GetPt(Float_t field)
{
  // transverse momentum calculation
  // curvature R = (fX*fX + fY*fY) / (2 * sin(alpha))
  // alpha = angle deviation from "infinite momentum"
  //
  // consistent with AliTRDmcmTracklet::GetPt(...)

  Float_t InfSlope = TMath::ATan(fY/fX)/TMath::Pi()*180.0;    
  Float_t alpha = fSlope - InfSlope;
  Float_t R = TMath::Sqrt(fX*fX + fY*fY)/(2.0*TMath::Sin(alpha/180.0*TMath::Pi()));
  
  Float_t Pt = 0.3 * field * 0.01 * R;
 
  return Pt;
 
}

//_____________________________________________________________________________
Int_t AliTRDltuTracklet::Compare(const TObject * o) const
{
  //
  // compare two LTU tracklets according to the intercept point Y1
  //

  AliTRDltuTracklet *ltutrk = (AliTRDltuTracklet*)o;

  if (fRow      != ltutrk->fRow)      return +1;
  if (fDetector != ltutrk->fDetector) return +1;

  if (fY <  ltutrk->fY) return -1;
  if (fY == ltutrk->fY) return  0;

  return 1;

}

//_____________________________________________________________________________
Float_t AliTRDltuTracklet::GetYproj(Float_t xpl)
{

  Float_t Yproj;

  Yproj = fY + TMath::Tan(fSlope/180.0*TMath::Pi()) * (xpl - fX);

  return Yproj;

}

//_____________________________________________________________________________
Float_t AliTRDltuTracklet::GetZproj(Float_t xpl)
{

  Float_t Zproj;

  Zproj = fRowz * xpl / fX;

  return Zproj;

}

//_____________________________________________________________________________
AliTRDgtuTrack::AliTRDgtuTrack()
{

  fYproj = 0.0;
  fZproj = 0.0;
  fSlope = 0.0;
  fDetector = -1;
  fNtracklets = 0;
  fNplanes = 0;
  fNclusters = 0;
  fPt = 0.0;
  fPhi = 0.0;
  fEta = 0.0;
  fLabel = -1;
  fPID = 0.0;
  fIsElectron = kFALSE;

  fTracklets = new TObjArray(400);

}

//_____________________________________________________________________________
AliTRDgtuTrack::AliTRDgtuTrack(const AliTRDgtuTrack& track):
  TObject(track),
  fYproj(track.fYproj),
  fZproj(track.fZproj),
  fSlope(track.fSlope),
  fDetector(track.fDetector),
  fNtracklets(track.fNtracklets),
  fNplanes(track.fNplanes),
  fNclusters(track.fNclusters),
  fPt(track.fPt),
  fPhi(track.fPhi),
  fEta(track.fEta),
  fLabel(track.fLabel),
  fPID(track.fPID),
  fIsElectron(track.fIsElectron)
{
  //
  // copy contructor
  //

  fTracklets = NULL;

}

//_____________________________________________________________________________
void AliTRDgtuTrack::AddTracklet(AliTRDltuTracklet *trk)
{

  Tracklets()->Add(trk);

}

//_____________________________________________________________________________
AliTRDltuTracklet* AliTRDgtuTrack::GetTracklet(Int_t pos)
{

  if (fTracklets == 0) return 0;
  void * trk = fTracklets->UncheckedAt(pos);
  if (trk == 0) return 0;

  return (AliTRDltuTracklet*)trk;

}

//_____________________________________________________________________________
Int_t AliTRDgtuTrack::Compare(const TObject * o) const
{

  AliTRDgtuTrack *gtutrack = (AliTRDgtuTrack*)o;

  if (fYproj <  gtutrack->GetYproj()) return -1;
  if (fYproj == gtutrack->GetYproj()) return  0;

  return +1;

}

//_____________________________________________________________________________
void AliTRDgtuTrack::Reset()
{

  fYproj = 0.0;
  fZproj = 0.0;
  fSlope = 0.0;
  fDetector = -1;
  fNtracklets = 0;
  fNplanes = 0;
  fNclusters = 0;
  fPt = 0.0;
  fPhi = 0.0;
  fEta = 0.0;
  fLabel = -1;
  fPID = 0.0;
  fIsElectron = kFALSE;
  
}

//_____________________________________________________________________________
void AliTRDgtuTrack::Track(Float_t xpl, Float_t field)
{

  AliTRDltuTracklet *trk;
  Int_t nTracklets = GetNtracklets();
  Float_t fC[kNmaxTrk][3];            // X, Y, Z  coordinates of segments

  fYproj = 0.0;
  fZproj = 0.0;
  fSlope = 0.0;
  fNclusters = 0;
  fNplanes = 0;
  fNtracklets = GetNtracklets();
  Int_t InDetector[kNplan];
  for (Int_t i = 0; i < kNplan; i++) InDetector[i] = -1;
  Int_t iDet, nDet = 0;
  Bool_t NewDetector;
  for (Int_t i = 0; i < nTracklets; i++) {

    trk = GetTracklet(i);
    fYproj += trk->GetYproj(xpl);
    fZproj += trk->GetZproj(xpl);
    fSlope += trk->GetSlope();
    fNclusters += trk->GetNclusters();
    iDet = trk->GetDetector();

    NewDetector = kTRUE;
    for (Int_t id = 0; id < nDet; id++) {
      if (iDet == InDetector[id]) {
	NewDetector = kFALSE;
	break;
      }
    }
    if (NewDetector) {
      InDetector[nDet++] = iDet;
      fNplanes++;
    }

    fC[i][0] = trk->GetTime0();
    fC[i][1] = trk->GetOffset();
    fC[i][2] = trk->GetRowz();

  }
  fYproj /= (Float_t)nTracklets;
  fZproj /= (Float_t)nTracklets;
  fSlope /= (Float_t)nTracklets;

  Float_t X[kNmaxTrk+1], Y[kNmaxTrk+1], Z[kNmaxTrk+1];
  Bool_t count[kNmaxTrk];
  for (Int_t i = 0; i < kNmaxTrk; i++) count[i] = kFALSE;

  Int_t iXmin = -1;
  Int_t j = 0;
  X[0] = Y[0] = Z[0] = 0.0;
  while (j < nTracklets) {
    iXmin = -1;
    for (Int_t i = 0; i < nTracklets; i++) {
      if (count[i]) continue;
      if (iXmin == -1) {
	iXmin = i;
	continue;
      }
      if (fC[i][0] < fC[iXmin][0]) iXmin = i;
    }
    X[j+1] = fC[iXmin][0];
    Y[j+1] = fC[iXmin][1];
    Z[j+1] = fC[iXmin][2];
    j++;
    count[iXmin] = kTRUE;
  }

  TMatrixD smatrix(2,2);
  TMatrixD sums(2,1);
  TMatrixD res(2,1);
  Double_t x, y;
  Float_t A, B;

  smatrix.Zero();
  sums.Zero();
  for (Int_t i = 0; i < nTracklets; i++) {
    x = (Double_t)X[i+1];
    y = (Double_t)Y[i+1];
    smatrix(0,0) += 1.0;
    smatrix(1,1) += x*x;
    smatrix(0,1) += x;
    smatrix(1,0) += x;
    sums(0,0)    += y;
    sums(1,0)    += x*y;
  }
  res = smatrix.Invert() * sums;
  A = res(0,0);
  B = res(1,0);

  Float_t dist = AliTRDgeometry::GetTime0(1) - AliTRDgeometry::GetTime0(0);

  Float_t fX1 = X[1]          + dist * (Float_t)(nTracklets-1)/6.0;
  Float_t fY1 = A + B * fX1;
  Float_t fX2 = X[nTracklets] - dist * (Float_t)(nTracklets-1)/6.0;
  Float_t fY2 = A + B * fX2;
  Float_t D12 = TMath::Sqrt((fX2-fX1)*(fX2-fX1)+(fY2-fY1)*(fY2-fY1));
  Float_t Alpha = TMath::ATan(fY2/fX2) - TMath::ATan(fY1/fX1);
  Float_t R = (D12/2.0)/TMath::Sin(Alpha);

  fPt = 0.3 * field * 0.01 * R;

  Float_t D1 = fX1*fX1+fY1*fY1;
  Float_t D2 = fX2*fX2+fY2*fY2;
  Float_t D  = fX1*fY2-fX2*fY1;
  
  Float_t Xc = (D1*fY2-D2*fY1)/(2*D);
  Float_t Yc = (D2*fX1-D1*fX2)/(2*D);

  if (Yc != 0.0) {
    fPhi = TMath::ATan(Xc/Yc);
  } else {
    fPhi = TMath::PiOver2();
  }

  fPhi *= 180.0/TMath::Pi();

  smatrix.Zero();
  sums.Zero();
  for (Int_t i = 0; i < nTracklets+1; i++) {
    x = (Double_t)Z[i];
    y = (Double_t)X[i];
    smatrix(0,0) += 1.0;
    smatrix(1,1) += x*x;
    smatrix(0,1) += x;
    smatrix(1,0) += x;
    sums(0,0)    += y;
    sums(1,0)    += x*y;
  }
  res = smatrix.Invert() * sums;
  A = res(0,0);
  B = res(1,0);
  Float_t theta = TMath::ATan(B);
  
  if (theta < 0.0) theta = TMath::Pi() + theta;
  
  if (theta == 0.0) {
    fEta = 0.0;
  } else {
    fEta = -TMath::Log(TMath::Tan(theta/2.0));
  }
  
}

//_____________________________________________________________________________
void AliTRDgtuTrack::MakePID()
{

  //
  // Electron likelihood signal
  //

  AliTRDcalibDB* calibration = AliTRDcalibDB::Instance();
  if (!calibration)
  {
    Error("MakePID","No instance of AliTRDcalibDB.");
    return;  
  }
  const AliTRDCalPIDLQ *pd = calibration->GetPIDLQObject();
  
  AliTRDltuTracklet *trk;
  Int_t nTracklets = GetNtracklets();
  Int_t det, pla;
  Float_t sl, th, q, probPio = 1.0, probEle = 1.0;
  for (Int_t i = 0; i < nTracklets; i++) {

    trk = GetTracklet(i);

    sl = TMath::Abs(trk->GetSlope());     // tracklet inclination in X-y plane
    th = trk->GetRowz()/trk->GetTime0();  // tracklet inclination in X-z plane
    th = TMath::ATan(TMath::Abs(th));

    q = trk->GetQ() 
      * TMath::Cos(sl/180.0*TMath::Pi()) 
      * TMath::Cos(th/180.0*TMath::Pi());

    det = trk->GetDetector();
    pla = trk->GetPlane(det);

    // unclear yet factor to match the PS distributions = 5.8
    // not explained only by the tail filter ...

    // AliRoot v4-03-07 , v4-03-Release
    //q = q * 5.8;

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

    probEle *= pd->GetProbability(0,TMath::Abs(fPt),q);
    probPio *= pd->GetProbability(2,TMath::Abs(fPt),q);

  }

  if ((probEle+probPio) > 0.0) {
    fPID = probEle/(probEle+probPio);
  } else {
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
  Float_t PIDcut = 0.829 - 0.032 * TMath::Abs(fPt);

  fIsElectron = kFALSE;
  if (fPID >= PIDcut) {
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
    trackCount[it] = 0;
  }

  Bool_t counted;
  Int_t label, nTracks = 0;
  for (Int_t itrk = 0; itrk < fNtracklets; itrk++) {

    trk = GetTracklet(itrk);

    if (trk->GetLabel() == -1) continue;

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

  Float_t frac = 4.0/5.0;
  for (Int_t it = 0; it < kMaxTracks; it++) {
    if (trackCount[it] >= (Int_t)(frac*fNtracklets)) {
      fLabel = trackLabel[it];
      break;
    }
  }

}

//_____________________________________________________________________________
AliTRDmodule::AliTRDmodule(AliTRDtrigParam *trigp) 
{

  //
  // AliTRDmodule default constructor
  //

  fDeltaY     = trigp->GetDeltaY();
  fDeltaS     = trigp->GetDeltaS();
  fXprojPlane = trigp->GetXprojPlane();
  fField      = trigp->GetField();
  fLTUtrk     = 0;
  fGTUtrk     = 0;
  fTracklets  = new TObjArray(400);
  fTracks     = new TObjArray(400);

}

//_____________________________________________________________________________
void AliTRDmodule::Reset() 
{

  ResetTracklets();
  ResetTracks();

  fLTUtrk     = 0;
  fGTUtrk     = 0;
  fTracklets  = new TObjArray(400);
  fTracks     = new TObjArray(400);

}

//_____________________________________________________________________________
void AliTRDmodule::ResetTracks() 
{

  if (fTracks) {

    AliTRDgtuTrack *trk;
    for (Int_t i = 0; i < GetNtracks(); i++) {
      
      trk = GetTrack(i);
      trk->Reset();
      
    }

    fTracks->Delete();

  }

}

//_____________________________________________________________________________
AliTRDgtuTrack* AliTRDmodule::GetTrack(Int_t pos)
{

  if (fTracks == 0) return 0;
  void * trk = fTracks->UncheckedAt(pos);
  if (trk == 0) return 0;

  return (AliTRDgtuTrack*)trk;

}

//_____________________________________________________________________________
void AliTRDmodule::RemoveTrack(Int_t pos)
{

  if (fTracks == 0) return;
  fTracks->RemoveAt(pos);
  fTracks->Compress();

}

//_____________________________________________________________________________
void AliTRDmodule::AddTracklet(Int_t det, 
			       Int_t row, 
			       Float_t rowz,
			       Float_t slope, 
			       Float_t offset, 
			       Float_t time, 
			       Int_t ncl,
			       Int_t label,
			       Float_t q) 
{
  
  fLTUtrk = new AliTRDltuTracklet(det,row,rowz,slope,offset,time,ncl,label,q);
  
  Tracklets()->Add(fLTUtrk);

}

//_____________________________________________________________________________
AliTRDltuTracklet* AliTRDmodule::GetTracklet(Int_t pos)
{

  if (fTracklets == 0) return 0;
  void * trk = fTracklets->UncheckedAt(pos);
  if (trk == 0) return 0;

  return (AliTRDltuTracklet*)trk;

}

//_____________________________________________________________________________
void AliTRDmodule::RemoveTracklet(Int_t pos)
{

  if (fTracklets == 0) return;
  fTracklets->RemoveAt(pos);
  fTracklets->Compress();

}

//_____________________________________________________________________________
void AliTRDmodule::RemoveMultipleTracklets()
{

  Float_t OffDiffMin = 0.5;  // [cm]

  AliTRDltuTracklet *trk;
  Int_t Det1, Det2, Row1, Row2, Ncl1, Ncl2, Label1, Label2;
  Float_t Off1, Off2;
  Int_t itrk = 0;
  while (itrk < (GetNtracklets()-1)) {

    trk = GetTracklet(itrk  );

    Det1 = trk->GetDetector();
    Row1 = trk->GetRow();
    Off1 = trk->GetOffset();
    Ncl1 = trk->GetNclusters();
    Label1 = trk->GetLabel();

    trk = GetTracklet(itrk+1);

    Det2 = trk->GetDetector();
    Row2 = trk->GetRow();
    Off2 = trk->GetOffset();
    Ncl2 = trk->GetNclusters();
    Label2 = trk->GetLabel();

    if (Det1 == Det2 && Row1 == Row2) {
      if ((Off2-Off1) < OffDiffMin) {
	if (Ncl1 < Ncl2) {
	  RemoveTracklet(itrk  );
	} else {
	  RemoveTracklet(itrk+1);
	}
      }
    }

    itrk++;

  }

}

//_____________________________________________________________________________
void AliTRDmodule::SortZ(Int_t cha)
{

  InitZLUT();

  AliTRDltuTracklet *trk;
  Int_t row, pla, det;
  for (Int_t iTrk = 0; iTrk < GetNtracklets(); iTrk++) {

    trk = GetTracklet(iTrk);

    row = trk->GetRow();
    det = trk->GetDetector();
    pla = trk->GetPlane(det);

    for (Int_t iZchan = 0; iZchan < kNsubZchan; iZchan++) {
      if (fZChannelMap[cha][iZchan][pla][row] == 1) {
	fZtrkid[pla][fZnchan[pla][iZchan]][iZchan] = iTrk;
	fZnchan[pla][iZchan]++;
      }
    }

  }

}

//_____________________________________________________________________________
void AliTRDmodule::InitZLUT()
{

  for (Int_t iPlan = 0; iPlan < AliTRDgeometry::Nplan(); iPlan++) {
    for (Int_t i = 0; i < kNsubZchan; i++) {
      fZnchan[iPlan][i] = 0;
      for (Int_t j = 0; j < kNmaxZchan; j++) {
        fZtrkid[iPlan][j][i] = -1;
      }
    }
  }

}

//_____________________________________________________________________________
void AliTRDmodule::FindTracks() 
{

  for (Int_t iZchan = 0; iZchan < kNsubZchan; iZchan++) {
    FindTracksCombi(iZchan);
  }

}

//_____________________________________________________________________________
void AliTRDmodule::FindTracksCombi(Int_t zchan) 
{

  //
  // find tracks by pure combinatorics...
  //
  
  static Int_t TrkTrack[12];
  
  Int_t nTracklets, nPlanes;
  Int_t Ntrk1, TrkId1, Ntrk2, TrkId2;
  Float_t Y1, Y1min, Y1max, S1, Z1, S1min, S1max, Y2, S2, Z2;
  AliTRDltuTracklet *Trk1;
  AliTRDltuTracklet *Trk2;
  AliTRDltuTracklet *Trk ;

  Bool_t IsPlane[kNplan];

  for (Int_t iPlan1 = 0; iPlan1 < kNplan; iPlan1++) {

    Ntrk1 = fZnchan[iPlan1][zchan];

    for (Int_t iTrk1 = 0; iTrk1 < Ntrk1; iTrk1++) {

      for (Int_t iPlan = 0; iPlan < kNplan; iPlan++) IsPlane[iPlan] = kFALSE;

      TrkId1 = fZtrkid[iPlan1][iTrk1][zchan];

      nTracklets = 0;
      for (Int_t iList = 0; iList < kNmaxTrk; iList++) {
	TrkTrack[iList] = -1;
      }

      TrkTrack[nTracklets++] = TrkId1;

      IsPlane[iPlan1] = kTRUE;

      Trk1 = GetTracklet(TrkId1);

      Y1    = Trk1->GetYproj(fXprojPlane);
      Y1min = Y1 - fDeltaY;
      Y1max = Y1 + fDeltaY;
      S1    = Trk1->GetSlope();
      S1min = S1 - fDeltaS;
      S1max = S1 + fDeltaS;
      Z1    = Trk1->GetZproj(fXprojPlane);      

      for (Int_t iPlan2 = 0; iPlan2 < kNplan; iPlan2++) {

	if (iPlan2 == iPlan1) continue;

	Ntrk2 = fZnchan[iPlan2][zchan];

	for (Int_t iTrk2 = 0; iTrk2 < Ntrk2; iTrk2++) {

	  TrkId2 = fZtrkid[iPlan2][iTrk2][zchan];

	  if (TrkId2 == TrkId1) continue;

	  Trk2 = GetTracklet(TrkId2);

	  Y2 = Trk2->GetYproj(fXprojPlane);
	  S2 = Trk2->GetSlope();
	  Z2 = Trk2->GetZproj(fXprojPlane);

	  if ((Y1min < Y2 && Y2 < Y1max) && 
	      (S1min < S2 && S2 < S1max)) {

	    if (nTracklets >= kNmaxTrk) {
	      Warning("FindTracksCombi","Too many tracklets for this track.");
	    } else {
	      TrkTrack[nTracklets++] = TrkId2;
	      IsPlane[iPlan2] = kTRUE;
	    }

	  }

	}  // end trk 2

      }  // end plan 2

      nPlanes = 0;
      for (Int_t iPlan = 0; iPlan < kNplan; iPlan++) {
	nPlanes += (Int_t)IsPlane[iPlan];
      }
      
      if (nPlanes >= 4) {

	Int_t Cha1, Cha2, Npoints1, Npoints2;
	for (Int_t iList = 0; iList < (nTracklets-1); iList++) {

	  if (TrkTrack[iList] == -1 || TrkTrack[iList+1] == -1) continue;
	  Trk1 = GetTracklet(TrkTrack[iList  ]);
	  Trk2 = GetTracklet(TrkTrack[iList+1]);

	  Cha1 = Trk1->GetDetector();
	  Cha2 = Trk2->GetDetector();
	  if (Cha1 != Cha2) continue;

	  Npoints1 = Trk1->GetNclusters();
	  Npoints2 = Trk2->GetNclusters();

	  if (Npoints1 == Npoints2) {
	    TrkTrack[iList] = -1;
	  } else {
	    if (Npoints1 > Npoints2) TrkTrack[iList+1] = -1;
	    if (Npoints1 < Npoints2) TrkTrack[iList  ] = -1;
	  }

	}

	fGTUtrk = new AliTRDgtuTrack();
	for (Int_t iList = 0; iList < nTracklets; iList++) {
	  if (TrkTrack[iList] == -1) continue;
	  Trk = GetTracklet(TrkTrack[iList]);
	  fGTUtrk->AddTracklet(Trk);
	}
	fGTUtrk->Track(fXprojPlane,fField);
	AddTrack();
      }
           
    }  // end trk 1

  }  // end plan 1

}

//_____________________________________________________________________________
void AliTRDmodule::AddTrack() 
{
  
  Tracks()->Add(fGTUtrk);

}

//_____________________________________________________________________________
void AliTRDmodule::RemoveMultipleTracks()
{

  AliTRDgtuTrack *trk1;
  AliTRDgtuTrack *trk2;

  Float_t Yproj1, Yproj2, Alpha1, Alpha2;
  Int_t ntrk1, ntrk2;
  Int_t iTrack = 0;

  while (iTrack < (GetNtracks()-1)) {

    trk1 = GetTrack(iTrack  );
    trk2 = GetTrack(iTrack+1);

    ntrk1  = trk1->GetNtracklets();
    Yproj1 = trk1->GetYproj();
    Alpha1 = trk1->GetSlope();
    ntrk2  = trk2->GetNtracklets();
    Yproj2 = trk2->GetYproj();
    Alpha2 = trk2->GetSlope();

    if (TMath::Abs(Yproj1-Yproj2) < fDeltaY && TMath::Abs(Alpha1-Alpha2) < fDeltaS) {
      if (ntrk1 < ntrk2) {
	RemoveTrack(iTrack  );
      } else {
	RemoveTrack(iTrack+1);
      }
    } else {
      iTrack++;
    }
    
  }

}

//_____________________________________________________________________________
AliTRDtrigger::AliTRDtrigger():
  TNamed(),
  fTracks("AliTRDgtuTrack",0)
{
  //
  // AliTRDtrigger default constructor
  //

  fDigitsManager = NULL;
  fTrackletTree = NULL;
  fTracklets     = NULL;

  fNROB = 0;
  fTrigParam = NULL;
  fMCM = NULL;
  fTrk = NULL;
  fGTUtrk = NULL;

  fNtracklets = 0;

  fDigits = NULL;
  fTrack0 = NULL;
  fTrack1 = NULL;
  fTrack2 = NULL;

  fModule = NULL;

  fNPrimary = 0;

}

//_____________________________________________________________________________
AliTRDtrigger::AliTRDtrigger(const Text_t *name, const Text_t *title):
  TNamed(name,title),
  fTracks("AliTRDgtuTrack",1000)
{
  //
  // AliTRDtrigger constructor
  //

  fDigitsManager = new AliTRDdigitsManager();
  fTrackletTree = NULL;
  fTracklets = new TObjArray(400);

  fNROB = 0;
  fTrigParam = NULL;
  fMCM = NULL;
  fTrk = NULL;
  fGTUtrk = NULL;

  fNtracklets = 0;

  fDigits = NULL;
  fTrack0 = NULL;
  fTrack1 = NULL;
  fTrack2 = NULL;

  fModule = NULL;

  fNPrimary = 0;

}

//_____________________________________________________________________________
AliTRDtrigger::AliTRDtrigger(const AliTRDtrigger &p):TNamed(p)
{
  //
  // AliTRDtrigger copy constructor
  //

  ((AliTRDtrigger &) p).Copy(*this);

}

///_____________________________________________________________________________
AliTRDtrigger::~AliTRDtrigger()
{
  //
  // AliTRDtrigger destructor
  //

  if (fTracklets) {
    fTracklets->Delete();
    delete fTracklets;
  }

  fTracks.Delete();

}

//_____________________________________________________________________________
AliTRDtrigger &AliTRDtrigger::operator=(const AliTRDtrigger &p)
{
  //
  // Assignment operator
  //

  if (this != &p) ((AliTRDtrigger &) p).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDtrigger::Copy(TObject &) const
{
  //
  // Copy function
  //

  AliFatal("Not implemented");

}

//_____________________________________________________________________________
void AliTRDtrigger::Init()
{

  fModule = new AliTRDmodule(fTrigParam); 
  /*
  AliHeader *header = fRunLoader->GetHeader();
  fNPrimary = header->GetNprimary();
  */
  fTracks.Clear();

}

//_____________________________________________________________________________
Bool_t AliTRDtrigger::Open(const Char_t *name, Int_t nEvent)
{
  //
  // Opens the AliROOT file.
  //

  TString evfoldname = AliConfig::GetDefaultEventFolderName();
  fRunLoader = AliRunLoader::GetRunLoader(evfoldname);

  if (!fRunLoader)
    fRunLoader = AliRunLoader::Open(name);

  if (!fRunLoader) {
    Error("Open","Can not open session for file %s.",name);
    return kFALSE;
  }

  // Open input

  if (fRunLoader->GetAliRun() == 0x0) fRunLoader->LoadgAlice();
  gAlice = fRunLoader->GetAliRun();

  if (!(gAlice)) {
    fRunLoader->LoadgAlice();
    gAlice = fRunLoader->GetAliRun();
    if (!(gAlice)) {
      Error("Open","Could not find AliRun object.");
      return kFALSE;
    }
  }

  // Import the Trees for the event nEvent in the file
  fRunLoader->GetEvent(nEvent);

  // Open output

  TObjArray *ioArray = 0;

  AliLoader* loader = fRunLoader->GetLoader("TRDLoader");
  loader->MakeTree("T");
  fTrackletTree = loader->TreeT();
  fTrackletTree->Branch("TRDmcmTracklet","TObjArray",&ioArray,32000,0);

  fRunLoader->LoadHeader();

  Init();

  return kTRUE;

}


//_____________________________________________________________________________
Bool_t AliTRDtrigger::ReadDigits() 
{
  //
  // Reads the digits arrays from the input aliroot file
  //

  if (!fRunLoader) {
    Error("ReadDigits","Can not find the Run Loader");
    return kFALSE;
  }

  AliLoader* loader = fRunLoader->GetLoader("TRDLoader");
  if (!loader->TreeD()) loader->LoadDigits();

  return (fDigitsManager->ReadDigits(loader->TreeD()));

}

//_____________________________________________________________________________
Bool_t AliTRDtrigger::ReadDigits(AliRawReader* rawReader)
{
  //
  // Reads the digits arrays from the ddl file
  //

  AliTRDrawData *raw = new AliTRDrawData();
  raw->SetDebug(1);

  fDigitsManager = raw->Raw2Digits(rawReader);

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDtrigger::ReadTracklets(AliRunLoader *rl) 
{
  //
  // Reads the tracklets find the tracks
  //

  Int_t idet;

  AliLoader *loader = rl->GetLoader("TRDLoader");
  loader->LoadTracks();
  fTrackletTree = loader->TreeT();

  TBranch *branch = fTrackletTree->GetBranch("TRDmcmTracklet");
  if (!branch) {
    Error("ReadTracklets","Can't get the branch !");
    return kFALSE;
  }
  TObjArray *tracklets = new TObjArray(400);
  branch->SetAddress(&tracklets);

  Int_t nEntries = (Int_t) fTrackletTree->GetEntries();
  Int_t iEntry, itrk;
  Int_t iStack, iStackPrev = -1;
  
  for (iEntry = 0; iEntry < nEntries; iEntry++) {    
    fTrackletTree->GetEvent(iEntry);
    
    for (itrk = 0; itrk < tracklets->GetEntriesFast(); itrk++){

      fTrk = (AliTRDmcmTracklet*)tracklets->UncheckedAt(itrk);
      
      idet = fTrk->GetDetector();

      iStack = idet / (AliTRDgeometry::Nplan());
      if (iStackPrev != iStack) {
	if (iStackPrev == -1) {
	  iStackPrev = iStack;
	} else {
	  MakeTracks(idet-AliTRDgeometry::Nplan());
	  ResetTracklets();
	  iStackPrev = iStack;
	}
      }
      
      Tracklets()->Add(fTrk);

      if (iEntry == (nEntries-1) && itrk == (tracklets->GetEntriesFast()-1)) {
	idet++;
	MakeTracks(idet-AliTRDgeometry::Nplan());
	ResetTracklets();
      }

    }

  }

  loader->UnloadTracks();

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDtrigger::MakeTracklets()
{

  AliTRDcalibDB* calibration = AliTRDcalibDB::Instance();
  if (!calibration)
  {
    Error("MakeTracklets","No instance of AliTRDcalibDB.");
    return kFALSE;  
  }
  
  AliTRDCommonParam* commonParam = AliTRDCommonParam::Instance();
  if (!commonParam)
  {
    Error("MakeTracklets","No common params.");
    return kFALSE;
  }
    
  AliTRDgeometry *geo = AliTRDgeometry::GetGeometry(fRunLoader);

  Int_t    chamBeg = 0;
  Int_t    chamEnd = AliTRDgeometry::Ncham();
  Int_t    planBeg = 0;
  Int_t    planEnd = AliTRDgeometry::Nplan();
  Int_t    sectBeg = 0;
  Int_t    sectEnd = AliTRDgeometry::Nsect();

  fMCM = new AliTRDmcm(fTrigParam,0);

  Int_t time, col, row, col1, col2;
  Float_t amp;
  Int_t idet, iStack, iStackPrev;
  iStack     = -1;
  iStackPrev = -1;
  for (Int_t isect = sectBeg; isect < sectEnd; isect++) {

    for (Int_t icham = chamBeg; icham < chamEnd; icham++) {

      // number of ROBs in the chamber
      if( icham == 2 ) {
	fNROB = 6;
      } else {
	fNROB = 8;
      }

      for (Int_t iplan = planBeg; iplan < planEnd; iplan++) {

        idet = geo->GetDetector(iplan,icham,isect);
	ResetTracklets();
	/*
	iStack = idet / (AliTRDgeometry::Nplan());
	if (iStackPrev != iStack) {
	  if (iStackPrev == -1) {
	    iStackPrev = iStack;
	  } else {
	    MakeTracks(idet-AliTRDgeometry::Nplan());
	    ResetTracklets();
	    iStackPrev = iStack;
	  }
	}
	*/
        Int_t    nRowMax     = commonParam->GetRowMax(iplan,icham,isect);
	Int_t    nColMax     = commonParam->GetColMax(iplan);
        Int_t    nTimeTotal  = calibration->GetNumberOfTimeBins();

        // Get the digits
        fDigits = fDigitsManager->GetDigits(idet);
        fDigits->Expand();
        fTrack0 = fDigitsManager->GetDictionary(idet,0);
        fTrack0->Expand();
        fTrack1 = fDigitsManager->GetDictionary(idet,1);
        fTrack1->Expand();
        fTrack2 = fDigitsManager->GetDictionary(idet,2); 
        fTrack2->Expand();


	for (Int_t iRob = 0; iRob < fNROB; iRob++) {

	  for (Int_t iMcm = 0; iMcm < kNMCM; iMcm++) {

	    fMCM->Reset();

	    fMCM->SetRobId(iRob);
	    fMCM->SetChaId(idet);

	    SetMCMcoordinates(iMcm);

	    row = fMCM->GetRow();

	    if (row < 0 || row > nRowMax) {
	      Error("MakeTracklets","MCM row number out of range.");
	    }

	    fMCM->GetColRange(col1,col2);
	    
            for (time = 0; time < nTimeTotal; time++) {
	      for (col = col1; col < col2; col++) {
		if (col >= 0 && col < nColMax) {
		  amp = TMath::Abs(fDigits->GetDataUnchecked(row,col,time));
		} else {
		  amp = 0.0;
		}
		fMCM->SetADC(col-col1,time,amp);

	      }
	    }

	    if (fTrigParam->GetTailCancelation()) {
	      fMCM->Filter(fTrigParam->GetNexponential(),fTrigParam->GetFilterType());
	    }
	    
	    if (fMCM->Run()) {

	      for (Int_t iSeed = 0; iSeed < kMaxTrackletsPerMCM; iSeed++) {
		
		if (fMCM->GetSeedCol()[iSeed] < 0) continue;

		if ( fTrigParam->GetDebugLevel() > 1 ) 
		  printf("Add tracklet %d in col %02d \n",fNtracklets,fMCM->GetSeedCol()[iSeed]);

		if ( fTrigParam->GetDebugLevel() == -1 ) {
		  printf("Add tracklet %d in col %02d \n",fNtracklets,fMCM->GetSeedCol()[iSeed]);
		  for (time = 0; time < nTimeTotal; time++) {
		    for (col = 0; col < kMcmCol; col++) {		    
		      printf("%03.0f  ",fMCM->GetADC(col,time));
		    }
		    printf("\n");
		  }
		}

		AddTracklet(idet,row,iSeed,fNtracklets++);

	      }

	    }

	  }

      
	}

	// Compress the arrays
        fDigits->Compress(1,0);
        fTrack0->Compress(1,0);
        fTrack1->Compress(1,0);
        fTrack2->Compress(1,0);

	WriteTracklets(idet);

     }
    }
  }
  /*
  idet++;
  MakeTracks(idet-AliTRDgeometry::Nplan());
  ResetTracklets();
  */
  return kTRUE;

}

//_____________________________________________________________________________
void AliTRDtrigger::SetMCMcoordinates(Int_t imcm)
{

  Int_t robid = fMCM->GetRobId();

  // setting the Row and Col range

  const Int_t kNcolRob = 2;  // number of ROBs per chamber in column direction
  const Int_t kNmcmRob = 4;  // number of MCMs per ROB in column/row direction

  Int_t mcmid = imcm%(kNmcmRob*kNmcmRob);

  if (robid%kNcolRob == 0) {

    if ( mcmid%kNmcmRob == 0 ) {
      fMCM->SetColRange(18*0-1,18*1-1+2+1);
    }
    if ( mcmid%kNmcmRob == 1 ) {
      fMCM->SetColRange(18*1-1,18*2-1+2+1);
    }
    if ( mcmid%kNmcmRob == 2 ) {
      fMCM->SetColRange(18*2-1,18*3-1+2+1);
    }
    if ( mcmid%kNmcmRob == 3 ) {
      fMCM->SetColRange(18*3-1,18*4-1+2+1);
    }

  } else {

    if ( mcmid%kNmcmRob == 0 ) {
      fMCM->SetColRange(18*4-1,18*5-1+2+1);
    }
    if ( mcmid%kNmcmRob == 1 ) {
      fMCM->SetColRange(18*5-1,18*6-1+2+1);
    }
    if ( mcmid%kNmcmRob == 2 ) {
      fMCM->SetColRange(18*6-1,18*7-1+2+1);
    }
    if ( mcmid%kNmcmRob == 3 ) {
      fMCM->SetColRange(18*7-1,18*8-1+2+1);
    }

  } 

  fMCM->SetRow(kNmcmRob*(robid/kNcolRob)+mcmid/kNmcmRob);

}

//_____________________________________________________________________________
void AliTRDtrigger::AddTracklet(Int_t det, Int_t row, Int_t seed, Int_t n)
{

  Float_t field = fTrigParam->GetField();
  AliTRDgeometry *geo = (AliTRDgeometry*)AliTRDgeometry::GetGeometry(fRunLoader);

  AliTRDcalibDB* calibration = AliTRDcalibDB::Instance();
  if (!calibration)
  {
    Error("AddTracklets","No instance of AliTRDcalibDB.");
    return;  
  }
  
  Int_t nTimeTotal  = calibration->GetNumberOfTimeBins();

  fTrk = new AliTRDmcmTracklet(det,row,n);

  Int_t iCol, iCol1, iCol2, track[3];
  iCol = fMCM->GetSeedCol()[seed];  // 0....20 (MCM)
  fMCM->GetColRange(iCol1,iCol2);   // range in the pad plane
	    
  Float_t Amp[3];
  for (Int_t iTime = 0; iTime < nTimeTotal; iTime++) {

    Amp[0] = fMCM->GetADC(iCol-1,iTime);
    Amp[1] = fMCM->GetADC(iCol  ,iTime);
    Amp[2] = fMCM->GetADC(iCol+1,iTime);

    // extract track contribution only from the central pad
    track[0] = fTrack0->GetDataUnchecked(row,iCol+iCol1,iTime);
    track[1] = fTrack1->GetDataUnchecked(row,iCol+iCol1,iTime);
    track[2] = fTrack2->GetDataUnchecked(row,iCol+iCol1,iTime);

    if (fMCM->IsCluster(iCol,iTime)) {

      fTrk->AddCluster(iCol+iCol1,iTime,Amp,track);

    } else if ((iCol+1+1) < kMcmCol) {

      Amp[0] = fMCM->GetADC(iCol-1+1,iTime);
      Amp[1] = fMCM->GetADC(iCol  +1,iTime);
      Amp[2] = fMCM->GetADC(iCol+1+1,iTime);

      if (fMCM->IsCluster(iCol+1,iTime)) {

	// extract track contribution only from the central pad
	track[0] = fTrack0->GetDataUnchecked(row,iCol+1+iCol1,iTime);
	track[1] = fTrack1->GetDataUnchecked(row,iCol+1+iCol1,iTime);
	track[2] = fTrack2->GetDataUnchecked(row,iCol+1+iCol1,iTime);

	fTrk->AddCluster(iCol+1+iCol1,iTime,Amp,track);

      }

    } else {
    }

  }

  fTrk->CookLabel(0.8);  
  /*
  if (fTrk->GetLabel() >= fNPrimary) {
    Info("AddTracklet","Only primaries are stored!");
    return;
  }
  */
  // LTU Pt cut
  fTrk->MakeTrackletGraph(geo,field);
  fTrk->MakeClusAmpGraph();
  if (TMath::Abs(fTrk->GetPt()) < fTrigParam->GetLtuPtCut()) return;
      
  Tracklets()->Add(fTrk);

}

//_____________________________________________________________________________
Bool_t AliTRDtrigger::WriteTracklets(Int_t det) 
{
  //
  // Fills TRDmcmTracklet branch in the tree with the Tracklets 
  // found in detector = det. For det=-1 writes the tree. 
  //

  if ((det < -1) || (det >= AliTRDgeometry::Ndet())) {
    Error("WriteTracklets","Unexpected detector index %d.",det);
    return kFALSE;
  }

  TBranch *branch = fTrackletTree->GetBranch("TRDmcmTracklet");
  if (!branch) {
    TObjArray *ioArray = 0;
    branch = fTrackletTree->Branch("TRDmcmTracklet","TObjArray",&ioArray,32000,0);
  }

  if ((det >= 0) && (det < AliTRDgeometry::Ndet())) {

    Int_t nTracklets = Tracklets()->GetEntriesFast();
    TObjArray *detTracklets = new TObjArray(400);

    for (Int_t i = 0; i < nTracklets; i++) {
      AliTRDmcmTracklet *trk = (AliTRDmcmTracklet *) Tracklets()->UncheckedAt(i);
      
      if (det == trk->GetDetector()) {
        detTracklets->AddLast(trk);
      }
      else {
      }
    }

    branch->SetAddress(&detTracklets);
    fTrackletTree->Fill();

    delete detTracklets;

    return kTRUE;

  }

  if (det == -1) {

    Info("WriteTracklets","Writing the Tracklet tree %s for event %d."
	 ,fTrackletTree->GetName(),fRunLoader->GetEventNumber());

    AliLoader* loader = fRunLoader->GetLoader("TRDLoader");
    loader->WriteTracks("OVERWRITE");
    
    return kTRUE;  

  }

  return kFALSE;

}

//_____________________________________________________________________________
void AliTRDtrigger::MakeTracks(Int_t det)
{
  //
  // Create GTU tracks per module (stack of 6 chambers)
  //
  
  fModule->Reset();

  AliTRDCommonParam* commonParam = AliTRDCommonParam::Instance();
  if (!commonParam)
  {
    Error("MakeTracks","No common params.");
    return;
  }
    
  Int_t nRowMax, iplan, icham, isect, row;

  AliTRDgeometry *geo = (AliTRDgeometry*)AliTRDgeometry::GetGeometry(fRunLoader);

  if ((det < 0) || (det >= AliTRDgeometry::Ndet())) {
    Error("MakeTracks","Unexpected detector index %d.",det);
    return;
  }
  
  Int_t nTracklets = Tracklets()->GetEntriesFast();
  
  AliTRDmcmTracklet *trk;
  for (Int_t i = 0; i < nTracklets; i++) {
    
    trk = (AliTRDmcmTracklet *) Tracklets()->UncheckedAt(i);
    
    iplan = geo->GetPlane(trk->GetDetector());
    icham = geo->GetChamber(trk->GetDetector());
    isect = geo->GetSector(trk->GetDetector());

    nRowMax = commonParam->GetRowMax(iplan,icham,isect);
    row = trk->GetRow();

    fModule->AddTracklet(trk->GetDetector(),
			 row,
			 trk->GetRowz(),
			 trk->GetSlope(),
			 trk->GetOffset(),
			 trk->GetTime0(),
			 trk->GetNclusters(),
			 trk->GetLabel(),
			 trk->GetdQdl());
    
  }

  fModule->SortTracklets();
  fModule->RemoveMultipleTracklets();
  fModule->SortZ((Int_t)geo->GetChamber(det));
  fModule->FindTracks();
  fModule->SortTracks();
  fModule->RemoveMultipleTracks();

  Int_t nModTracks = fModule->GetNtracks();
  AliTRDgtuTrack *gtutrk;
  for (Int_t i = 0; i < nModTracks; i++) {
    gtutrk = (AliTRDgtuTrack*)fModule->GetTrack(i);
    if (TMath::Abs(gtutrk->GetPt()) < fTrigParam->GetGtuPtCut()) continue;
    gtutrk->CookLabel();
    gtutrk->MakePID();
    AddTrack(gtutrk,det);
  }
  
}


