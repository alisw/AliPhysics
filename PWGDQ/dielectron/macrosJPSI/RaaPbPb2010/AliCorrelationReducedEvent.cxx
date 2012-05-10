#ifndef ALICORRELATIONREDUCEDEVENT_H
#include "AliCorrelationReducedEvent.h"
#endif

#include <TMath.h>
#include <TClonesArray.h>

ClassImp(AliCorrelationReducedTrack)
ClassImp(AliCorrelationReducedPair)
ClassImp(AliCorrelationReducedEvent)
ClassImp(AliCorrelationReducedEventFriend)
ClassImp(AliCorrelationReducedCaloCluster)

TClonesArray* AliCorrelationReducedEvent::fgTracks = 0;
TClonesArray* AliCorrelationReducedEvent::fgCandidates = 0;
TClonesArray* AliCorrelationReducedEvent::fgCaloClusters = 0;

//_______________________________________________________________________________
AliCorrelationReducedTrack::AliCorrelationReducedTrack() :
  fTrackId(-1),
  fStatus(0),
  fGlobalPhi(0.0),
  fGlobalPt(0.0),
  fGlobalEta(0.0),
  fTPCPhi(0.0),
  fTPCPt(0.0),
  fTPCEta(0.0),
  fMomentumInner(0.0),
  fDCA(),
  fITSclusterMap(0),
  fITSsignal(0.0),
  fTPCNcls(0),
  fTPCCrossedRows(0),
  fTPCNclsF(0),
  fTPCNclsIter1(0),
  fTPCClusterMap(0),
  fTPCsignal(0),
  fTPCnSig(),
  fTOFbeta(0.0),
  fTOFnSig(),
  fTRDntracklets(),
  fTRDpid(),
  fCaloClusterId(-999),
  fBayesPID(),
  fFlags(0)
{
  //
  // Constructor
  //
  fDCA[0] = 0.0; fDCA[1]=0.0;
  for(Int_t i=0; i<4; ++i) {fTPCnSig[i]=-999.; fTOFnSig[i]=-999.;} 
  for(Int_t i=0; i<3; ++i) {fBayesPID[i]=-999.;}
  fTRDpid[0]=-999.; fTRDpid[1]=-999.;
}


//_______________________________________________________________________________
AliCorrelationReducedTrack::~AliCorrelationReducedTrack()
{
  //
  // De-Constructor
  //
}


//_______________________________________________________________________________
AliCorrelationReducedPair::AliCorrelationReducedPair() :
  fCandidateId(-1),
  fPairType(-1), 
  fLegIds(),
  fMass(),
  fPhi(0.0),
  fPt(0.0),
  fEta(0.0),
  fLxy(0.0),
  fOpeningAngle(-1.0),
  //fOnTheFly(kFALSE),
  fMCid(0)
{
  //
  // Constructor
  //
  fLegIds[0] = -1; fLegIds[1] = -1;
  fMass[0]=-999.; fMass[1]=-999.; fMass[2]=-999.;
}


//_______________________________________________________________________________
AliCorrelationReducedPair::AliCorrelationReducedPair(const AliCorrelationReducedPair &c) :
  TObject(c),
  fCandidateId(c.CandidateId()),
  fPairType(c.PairType()),
  fLegIds(),
  fMass(),
  fPhi(c.Phi()),
  fPt(c.Pt()),
  fEta(c.Eta()),
  fLxy(c.Lxy()),
  fOpeningAngle(c.OpeningAngle()),
  //fOnTheFly(c.IsOnTheFly()),
  fMCid(c.MCid())
{
  //
  // copy constructor
  //
  fLegIds[0] = c.LegId(0);
  fLegIds[1] = c.LegId(1);
  fMass[0] = c.Mass(0); fMass[1] = c.Mass(1); fMass[2] = c.Mass(2);
}


//_______________________________________________________________________________
AliCorrelationReducedPair::~AliCorrelationReducedPair() {
  //
  // destructor
  //
}

//____________________________________________________________________________
AliCorrelationReducedEvent::AliCorrelationReducedEvent() :
  fRunNo(0),
  fBC(0),
  fTriggerMask(0),
  fIsPhysicsSelection(kTRUE),
  fVtx(),
  fNVtxContributors(0),
  fVtxTPC(),
  fNVtxTPCContributors(0),
  fCentrality(),
  fCentQuality(0),
  fNV0candidates(),
  fNDielectronCandidates(0),
  fNtracks(),
  fSPDntracklets(0),
  fVZEROMult(),
  fZDCnEnergy(),
  fNCaloClusters(0)
//fFMDMult()
{
  //
  // Constructor
  //
  for(Int_t i=0; i<3; ++i) {fVtx[i]=-999.; fVtxTPC[i]=-999.;}
  for(Int_t i=0; i<4; ++i) fCentrality[i]=-1.;
  fNV0candidates[0]=0; fNV0candidates[1]=0;
  fNtracks[0]=0; fNtracks[1]=0;
  for(Int_t i=0; i<64; ++i) fVZEROMult[i] = 0.0;
  for(Int_t i=0; i<8; ++i) fZDCnEnergy[i]=0.0;
  
  if(!fgTracks) fgTracks = new TClonesArray("AliCorrelationReducedTrack", 100000);
  fTracks = fgTracks;
  if(!fgCandidates) fgCandidates = new TClonesArray("AliCorrelationReducedPair", 100000);
  fCandidates = fgCandidates;
  if(!fgCaloClusters) fgCaloClusters = new TClonesArray("AliCorrelationReducedCaloCluster", 50000);
  fCaloClusters = fgCaloClusters;
}


//____________________________________________________________________________
AliCorrelationReducedEvent::~AliCorrelationReducedEvent()
{
  //
  // De-Constructor
  //
  ClearEvent();
}


//____________________________________________________________________________
Float_t AliCorrelationReducedEvent::MultVZEROA() const
{
  //
  // Total VZERO multiplicity in A side
  //
  Float_t mult=0.0;
  for(Int_t i=32;i<64;++i)
    mult += fVZEROMult[i];
  return mult;
}


//____________________________________________________________________________
Float_t AliCorrelationReducedEvent::MultVZEROC() const
{
  //
  // Total VZERO multiplicity in C side
  //
  Float_t mult=0.0;
  for(Int_t i=0;i<32;++i)
    mult += fVZEROMult[i];
  return mult;
}


//____________________________________________________________________________
Float_t AliCorrelationReducedEvent::MultVZERO() const
{
  //
  // Total VZERO multiplicity
  //
  return MultVZEROA()+MultVZEROC();
}


//____________________________________________________________________________
Float_t AliCorrelationReducedEvent::MultRingVZEROA(Int_t ring) const 
{
  //
  //  VZERO multiplicity in a given ring on A side
  //
  if(ring<0 || ring>3) return -1.0;

  Float_t mult=0.0;
  for(Int_t i=32+8*ring; i<32+8*(ring+1); ++i)
    mult += fVZEROMult[i];
  return mult;
}


//____________________________________________________________________________
Float_t AliCorrelationReducedEvent::MultRingVZEROC(Int_t ring) const 
{
  //
  //  VZERO multiplicity in a given ring on C side
  //
  if(ring<0 || ring>3) return -1.0;

  Float_t mult=0.0;
  for(Int_t i=8*ring; i<8*(ring+1); ++i)
    mult += fVZEROMult[i];
  return mult;
}

//_____________________________________________________________________________
void AliCorrelationReducedEvent::ClearEvent() {
  //
  // clear the event
  //
  fTracks->Clear("C");
  fCandidates->Clear("C");
  fCaloClusters->Clear("C");
  fRunNo = 0;
  fBC = 0;
  fTriggerMask = 0;
  fIsPhysicsSelection = kTRUE;
  fCentQuality = 0;
  fNV0candidates[0] = 0; fNV0candidates[1] = 0;
  fNDielectronCandidates = 0;
  fNtracks[0] = 0; fNtracks[1] = 0;
  for(Int_t i=0; i<3; ++i) {fVtx[i]=-999.; fVtxTPC[i]=-999.; fCentrality[i]=-1.;}
  for(Int_t i=0; i<64; ++i) fVZEROMult[i] = 0.0;
  for(Int_t i=0; i<8; ++i) fZDCnEnergy[i]=0.0;
}


//_______________________________________________________________________________
AliCorrelationReducedEventFriend::AliCorrelationReducedEventFriend() :
 fQvector(),
 fEventPlaneStatus()
{
  //
  // default constructor
  //
  for(Int_t idet=0; idet<kNdetectors; ++idet) {
    for(Int_t ih=0; ih<fgkNMaxHarmonics; ++ih) {
      fEventPlaneStatus[idet][ih] = 0;
      for(Int_t ic=0; ic<2; ++ic)
	fQvector[idet][ih][ic] = 0.0;
    }
  }
}


//____________________________________________________________________________
AliCorrelationReducedEventFriend::~AliCorrelationReducedEventFriend()
{
  //
  // De-Constructor
  //
  ClearEvent();
}


//_____________________________________________________________________________
void AliCorrelationReducedEventFriend::ClearEvent() {
  //
  // clear the event
  //
  for(Int_t idet=0; idet<kNdetectors; ++idet) {
    for(Int_t ih=0; ih<fgkNMaxHarmonics; ++ih) {
      fEventPlaneStatus[idet][ih] = 0;
      for(Int_t ic=0; ic<2; ++ic)
	fQvector[idet][ih][ic] = 0.0;
    }
  }
}


//____________________________________________________________________________
void AliCorrelationReducedEventFriend::CopyEvent(AliCorrelationReducedEventFriend* event) {
  //
  // copy information from another event to this one
  //
  for(Int_t idet=0; idet<kNdetectors; ++idet) {
    for(Int_t ih=0; ih<fgkNMaxHarmonics; ++ih) {
      fQvector[idet][ih][0] = event->Qx(idet, ih+1);
      fQvector[idet][ih][1] = event->Qy(idet, ih+1);
      fEventPlaneStatus[idet][ih] = event->GetEventPlaneStatus(idet, ih+1);
    }
  }
}


//_____________________________________________________________________________
AliCorrelationReducedCaloCluster::AliCorrelationReducedCaloCluster() :
 fType(kUndefined),
 fEnergy(-999.),
 fTrackDx(-999.),
 fTrackDz(-999.)
{
  //
  // default constructor
  //
}


//_____________________________________________________________________________
AliCorrelationReducedCaloCluster::~AliCorrelationReducedCaloCluster()
{
  //
  // destructor
  //
}


//_______________________________________________________________________________
void AliCorrelationReducedEvent::GetQvector(Double_t Qvec[][2], Int_t det) {
  //
  // Get the event plane for a specified detector
  //
  if(det==AliCorrelationReducedEventFriend::kTPC || 
     det==AliCorrelationReducedEventFriend::kTPCptWeights ||
     det==AliCorrelationReducedEventFriend::kTPCpos ||
     det==AliCorrelationReducedEventFriend::kTPCneg) {
    GetTPCQvector(Qvec, det);
    return;
  }
  if(det==AliCorrelationReducedEventFriend::kVZEROA ||
     det==AliCorrelationReducedEventFriend::kVZEROC) {
    GetVZEROQvector(Qvec, det);   
    return;
  }
  if(det==AliCorrelationReducedEventFriend::kZDCA ||
     det==AliCorrelationReducedEventFriend::kZDCC) {
    GetZDCQvector(Qvec, det);
    return;
  }
  if(det==AliCorrelationReducedEventFriend::kFMD) {
    //TODO implementation
    return;
  }
  return;
}


//_________________________________________________________________
void AliCorrelationReducedEvent::GetTPCQvector(Double_t Qvec[][2], Int_t det) {
  //
  // Construct the event plane using tracks in the barrel
  //
  if(!(det==AliCorrelationReducedEventFriend::kTPC ||
       det==AliCorrelationReducedEventFriend::kTPCpos ||
       det==AliCorrelationReducedEventFriend::kTPCneg))
    return;
  AliCorrelationReducedTrack* track=0x0;
  Double_t pt=0.0; Double_t x=0.0; Double_t y=0.0; 
  Short_t charge=0;
  TIter nextTrack(fTracks);
  while((track=static_cast<AliCorrelationReducedTrack*>(nextTrack()))) {
    if(!track->UsedForQvector()) continue;
    charge = track->Charge();
    if(det==AliCorrelationReducedEventFriend::kTPCpos && charge<0) continue;
    if(det==AliCorrelationReducedEventFriend::kTPCneg && charge>0) continue;
    pt = 1.0;
    if(det==AliCorrelationReducedEventFriend::kTPCptWeights) {
      pt = track->Pt();
      if(pt>2.0) pt = 2.0;    // pt is the weight used for the event plane
    }
    x = TMath::Cos(track->Phi());
    y = TMath::Sin(track->Phi());
    //  1st harmonic  
    Qvec[0][0] += pt*x;
    Qvec[0][1] += pt*y;
    //  2nd harmonic
    Qvec[1][0] += pt*(2.0*TMath::Power(x,2.0)-1);
    Qvec[1][1] += pt*(2.0*x*y);
    //  3rd harmonic
    Qvec[2][0] += pt*(4.0*TMath::Power(x,3.0)-3.0*x);
    Qvec[2][1] += pt*(3.0*y-4.0*TMath::Power(y,3.0));
    //  4th harmonic
    Qvec[3][0] += pt*(1.0-8.0*TMath::Power(x*y,2.0));
    Qvec[3][1] += pt*(4.0*x*y-8.0*x*TMath::Power(y,3.0));
    //  5th harmonic
    Qvec[4][0] += pt*(16.0*TMath::Power(x,5.0)-20.0*TMath::Power(x, 3.0)+5.0*x);
    Qvec[4][1] += pt*(16.0*TMath::Power(y,5.0)-20.0*TMath::Power(y, 3.0)+5.0*y);
    //  6th harmonic
    Qvec[5][0] += pt*(32.0*TMath::Power(x,6.0)-48.0*TMath::Power(x, 4.0)+18.0*TMath::Power(x, 2.0)-1.0);
    Qvec[5][1] += pt*(x*y*(32.0*TMath::Power(y,4.0)-32.0*TMath::Power(y, 2.0)+6.0)); 
  }
}


//_________________________________________________________________
void AliCorrelationReducedEvent::GetVZEROQvector(Double_t Qvec[][2], Int_t det) {
  //
  // Get the reaction plane from the VZERO detector for a given harmonic
  //
  GetVZEROQvector(Qvec, det, fVZEROMult);
}


//_________________________________________________________________
void AliCorrelationReducedEvent::GetVZEROQvector(Double_t Qvec[][2], Int_t det, Float_t* vzeroMult) {
  //
  // Get the reaction plane from the VZERO detector for a given harmonic
  //
  //  Q{x,y} = SUM_i mult(i) * {cos(n*phi_i), sin(n*phi_i)} 
  //  phi_i - phi angle of the VZERO sector i
  //          Each sector covers 45 degrees(8 sectors per ring). Middle of sector 0 is at 45/2
  //        channel 0: 22.5
  //                1: 22.5+45
  //                2: 22.5+45*2
  //               ...
  //        at the next ring continues the same
  //        channel 8: 22.5
  //        channel 9: 22.5 + 45
  //       
  if(!(det==AliCorrelationReducedEventFriend::kVZEROA ||
       det==AliCorrelationReducedEventFriend::kVZEROC))
    return; 
  
  const Double_t kX[8] = {0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268, 0.38268, 0.92388};    // cosines of the angles of the VZERO sectors (8 per ring)
  const Double_t kY[8] = {0.38268, 0.92388, 0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268};    // sines     -- " --
  Int_t phi;
  
  for(Int_t iChannel=0; iChannel<64; ++iChannel) {
    if(iChannel<32 && det==AliCorrelationReducedEventFriend::kVZEROA) continue;
    if(iChannel>=32 && det==AliCorrelationReducedEventFriend::kVZEROC) continue;
    phi=iChannel%8;
    //  1st harmonic  
    Qvec[0][0] += vzeroMult[iChannel]*kX[phi];
    Qvec[0][1] += vzeroMult[iChannel]*kY[phi];
    //  2nd harmonic
    Qvec[1][0] += vzeroMult[iChannel]*(2.0*TMath::Power(kX[phi],2.0)-1);
    Qvec[1][1] += vzeroMult[iChannel]*(2.0*kX[phi]*kY[phi]);
    //  3rd harmonic
    Qvec[2][0] += vzeroMult[iChannel]*(4.0*TMath::Power(kX[phi],3.0)-3.0*kX[phi]);
    Qvec[2][1] += vzeroMult[iChannel]*(3.0*kY[phi]-4.0*TMath::Power(kY[phi],3.0));
    //  4th harmonic
    Qvec[3][0] += vzeroMult[iChannel]*(1.0-8.0*TMath::Power(kX[phi]*kY[phi],2.0));
    Qvec[3][1] += vzeroMult[iChannel]*(4.0*kX[phi]*kY[phi]-8.0*kX[phi]*TMath::Power(kY[phi],3.0));
    //  5th harmonic
    Qvec[4][0] += vzeroMult[iChannel]*(16.0*TMath::Power(kX[phi],5.0)-20.0*TMath::Power(kX[phi], 3.0)+5.0*kX[phi]);
    Qvec[4][1] += vzeroMult[iChannel]*(16.0*TMath::Power(kY[phi],5.0)-20.0*TMath::Power(kY[phi], 3.0)+5.0*kY[phi]);
    //  6th harmonic
    Qvec[5][0] += vzeroMult[iChannel]*(32.0*TMath::Power(kX[phi],6.0)-48.0*TMath::Power(kX[phi], 4.0)+18.0*TMath::Power(kX[phi], 2.0)-1.0);
    Qvec[5][1] += vzeroMult[iChannel]*(kX[phi]*kY[phi]*(32.0*TMath::Power(kY[phi],4.0)-32.0*TMath::Power(kY[phi], 2.0)+6.0));
  }    // end loop over channels 
}


//_________________________________________________________________
void AliCorrelationReducedEvent::GetZDCQvector(Double_t Qvec[][2], Int_t det) {
  //
  // Construct the event plane using the ZDC
  // ZDC has 2 side (A and C) with 4 calorimeters on each side  
  // The XY position of each calorimeter is specified by the 
  // zdcNTowerCenters_x and zdcNTowerCenters_y arrays
  if(!(det==AliCorrelationReducedEventFriend::kZDCA || 
       det==AliCorrelationReducedEventFriend::kZDCC)) return;       // bad detector option
  const Float_t zdcTowerCenter = 1.75;
  const Float_t zdcNTowerCenters_x[4] = {-zdcTowerCenter,  zdcTowerCenter, -zdcTowerCenter, zdcTowerCenter};
  const Float_t zdcNTowerCenters_y[4] = {-zdcTowerCenter, -zdcTowerCenter,  zdcTowerCenter, zdcTowerCenter};

  Qvec[0][0] = 0.0; Qvec[0][1] = 0.0;   // first harmonic Q-vector
  Float_t zdcNCentroidSum = 0;
  Float_t zdcN_alpha = 0.5;
    
  for(Int_t i=0; i<4; ++i) {
    if(fZDCnEnergy[i+(det==AliCorrelationReducedEventFriend::kZDCA ? 4 : 0)]>0.0) {
      Float_t zdcNenergy_alpha = TMath::Power(fZDCnEnergy[i+(det==AliCorrelationReducedEventFriend::kZDCA ? 4 : 0)], zdcN_alpha);
      Qvec[0][0] += zdcNTowerCenters_x[i-1]*zdcNenergy_alpha;
      Qvec[0][1] += zdcNTowerCenters_y[i-1]*zdcNenergy_alpha;
      zdcNCentroidSum += zdcNenergy_alpha;
    }
  }   // end loop over channels
  
  if(zdcNCentroidSum>0.0) {
    Qvec[0][0] /= zdcNCentroidSum;
    Qvec[0][1] /= zdcNCentroidSum;
  }
}
