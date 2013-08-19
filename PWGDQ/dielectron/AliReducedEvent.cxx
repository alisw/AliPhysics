/*
***********************************************************
  Implementation of reduced ESD information classes for
  quick analysis.
  Contact: i.c.arsene@gsi.de, i.c.arsene@cern.ch
  2012/06/21
  *********************************************************
*/

#ifndef ALIREDUCEDEVENT_H
#include "AliReducedEvent.h"
#endif

#include <TMath.h>
#include <TClonesArray.h>
#include <TClass.h>

ClassImp(AliReducedTrack)
ClassImp(AliReducedPair)
ClassImp(AliReducedFMD)
ClassImp(AliReducedEvent)
ClassImp(AliReducedEventFriend)
ClassImp(AliReducedCaloCluster)

TClonesArray* AliReducedEvent::fgTracks = 0;
TClonesArray* AliReducedEvent::fgCandidates = 0;
TClonesArray* AliReducedEvent::fgCaloClusters = 0;
TClonesArray* AliReducedEvent::fgFMD1 = 0;
TClonesArray* AliReducedEvent::fgFMD2I = 0;
TClonesArray* AliReducedEvent::fgFMD2O = 0;
TClonesArray* AliReducedEvent::fgFMD3I = 0;
TClonesArray* AliReducedEvent::fgFMD3O = 0;

//_______________________________________________________________________________
AliReducedTrack::AliReducedTrack() :
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
  fTrackLength(0.0),
  fITSclusterMap(0),
  fITSsignal(0.0),
  fITSnSig(),
  fITSchi2(0.0),
  fTPCNcls(0),
  fTPCCrossedRows(0),
  fTPCNclsF(0),
  fTPCNclsIter1(0),
  fTPCClusterMap(0),
  fTPCsignal(0),
  fTPCsignalN(0),
  fTPCnSig(),
  fTPCchi2(0.0),
  fTOFbeta(0.0),
  fTOFnSig(),
  fTOFdeltaBC(0),
  fTRDntracklets(),
  fTRDpid(),
  fTRDpidLQ2D(),
  fCaloClusterId(-999),
  fBayesPID(),
  fFlags(0),
  fMoreFlags(0)
{
  //
  // Constructor
  //
  fDCA[0] = 0.0; fDCA[1]=0.0;
  for(Int_t i=0; i<4; ++i) {fTPCnSig[i]=-999.; fTOFnSig[i]=-999.; fITSnSig[i]=-999.;} 
  for(Int_t i=0; i<3; ++i) {fBayesPID[i]=-999.;}
  fTRDpid[0]=-999.; fTRDpid[1]=-999.;
  fTRDpidLQ2D[0] = -999.; fTRDpidLQ2D[1] = -999.;
}


//_______________________________________________________________________________
AliReducedTrack::~AliReducedTrack()
{
  //
  // De-Constructor
  //
}


//_______________________________________________________________________________
AliReducedPair::AliReducedPair() :
  fCandidateId(-1),
  fPairType(-1), 
  fLegIds(),
  fMass(),
  fPhi(0.0),
  fPt(0.0),
  fEta(0.0),
  fLxy(0.0),
  fLxyErr(0.0),
  fPointingAngle(0.0),
  fChisquare(0.0),
  fMCid(0)
{
  //
  // Constructor
  //
  fLegIds[0] = -1; fLegIds[1] = -1;
  fMass[0]=-999.; fMass[1]=-999.; fMass[2]=-999.; fMass[3]=-999.;
}


//_______________________________________________________________________________
AliReducedPair::AliReducedPair(const AliReducedPair &c) :
  TObject(c),
  fCandidateId(c.CandidateId()),
  fPairType(c.PairType()),
  fLegIds(),
  fMass(),
  fPhi(c.Phi()),
  fPt(c.Pt()),
  fEta(c.Eta()),
  fLxy(c.Lxy()),
  fLxyErr(c.LxyErr()),
  fPointingAngle(c.PointingAngle()),
  fChisquare(c.Chi2()),
  fMCid(c.MCid())
{
  //
  // copy constructor
  //
  fLegIds[0] = c.LegId(0);
  fLegIds[1] = c.LegId(1);
  fMass[0] = c.Mass(0); fMass[1] = c.Mass(1); fMass[2] = c.Mass(2); fMass[3] = c.Mass(3);
}


//_______________________________________________________________________________
AliReducedPair::~AliReducedPair() {
  //
  // destructor
  //
}


//_______________________________________________________________________________
AliReducedFMD::AliReducedFMD() :
  fMultiplicity(0),
  //fEta(0.0),
  fId(0)
{
  AliReducedFMD::Class()->IgnoreTObjectStreamer();
  //
  // Constructor
  //
}


//_______________________________________________________________________________
AliReducedFMD::~AliReducedFMD()
{
  //
  // De-Constructor
  //
}


//____________________________________________________________________________
AliReducedEvent::AliReducedEvent() :
  TObject(),
  fEventTag(0),
  fEventNumberInFile(0),
  fL0TriggerInputs(0),
  fL1TriggerInputs(0),
  fL2TriggerInputs(0),
  fRunNo(0),
  fBC(0),
  fTimeStamp(0),
  fEventType(0),
  fTriggerMask(0),
  fIsPhysicsSelection(kTRUE),
  fIsSPDPileup(kFALSE),
  fIsSPDPileupMultBins(kFALSE),
  fIRIntClosestIntMap(),
  fIsFMDReduced(kTRUE),
  fVtx(),
  fNVtxContributors(0),
  fVtxTPC(),
  fNVtxTPCContributors(0),
  fNpileupSPD(0),
  fNpileupTracks(0),
  fNPMDtracks(0),
  fNTRDtracks(0),
  fNTRDtracklets(0),
  fCentrality(),
  fCentQuality(0),
  fNV0candidates(),
  fNDielectronCandidates(0),
  fNtracks(),
  fSPDntracklets(0),
  fVZEROMult(),
  fZDCnEnergy(),
  fZDCpEnergy(),
  fT0amplitude(),
  fT0TOF(),
  fT0TOFbest(),
  fT0zVertex(0),
  fT0start(0),
  fT0pileup(kFALSE),
  fT0sattelite(kFALSE),
  fTracks(0x0),
  fCandidates(0x0),
  fNFMDchannels(),
  fFMDtotalMult(),
  fFMD1(0x0),
  fFMD2I(0x0),
  fFMD2O(0x0),
  fFMD3I(0x0),
  fFMD3O(0x0),
  fNCaloClusters(0),
  fCaloClusters(0x0)
{
  //
  // Constructor
  //
  for(Int_t i=0; i<2; ++i) fIRIntClosestIntMap[i] = 0;
  for(Int_t i=0; i<3; ++i) {fVtx[i]=-999.; fVtxTPC[i]=-999.;}
  for(Int_t i=0; i<4; ++i) fCentrality[i]=-1.;
  fNV0candidates[0]=0; fNV0candidates[1]=0;
  fNtracks[0]=0; fNtracks[1]=0;
  for(Int_t i=0; i<32; ++i) fSPDntrackletsEta[i]=0;
  for(Int_t i=0; i<32; ++i) fNtracksPerTrackingFlag[i]=0;
  for(Int_t i=0; i<64; ++i) fVZEROMult[i] = 0.0;
  for(Int_t i=0; i<10; ++i) fZDCnEnergy[i]=0.0;
  for(Int_t i=0; i<10; ++i) fZDCpEnergy[i]=0.0;
  for(Int_t i=0; i<26; ++i) fT0amplitude[i]=0.0;
  for(Int_t i=0; i<3; ++i)  fT0TOF[i]=0.0;
  for(Int_t i=0; i<3; ++i)  fT0TOFbest[i]=0.0;
  for(Int_t i=0; i<5; ++i)  fNFMDchannels[i]=0.0;
  for(Int_t i=0; i<5; ++i)  fFMDtotalMult[i]=0.0;
}


//____________________________________________________________________________
AliReducedEvent::AliReducedEvent(const Char_t* /*name*/) :
  TObject(),
  fEventTag(0),
  fEventNumberInFile(0),
  fL0TriggerInputs(0),
  fL1TriggerInputs(0),
  fL2TriggerInputs(0),
  fRunNo(0),
  fBC(0),
  fTimeStamp(0),
  fEventType(0),
  fTriggerMask(0),
  fIsPhysicsSelection(kTRUE),
  fIsSPDPileup(kFALSE),
  fIsSPDPileupMultBins(kFALSE),
  fIRIntClosestIntMap(),
  fIsFMDReduced(kTRUE),
  fVtx(),
  fNVtxContributors(0),
  fVtxTPC(),
  fNVtxTPCContributors(0),
  fNpileupSPD(0),
  fNpileupTracks(0),
  fNPMDtracks(0),
  fNTRDtracks(0),
  fNTRDtracklets(0),
  fCentrality(),
  fCentQuality(0),
  fNV0candidates(),
  fNDielectronCandidates(0),
  fNtracks(),
  fSPDntracklets(0),
  fVZEROMult(),
  fZDCnEnergy(),
  fZDCpEnergy(),
  fT0amplitude(),
  fT0TOF(),
  fT0TOFbest(),
  fT0zVertex(0),
  fT0start(0),
  fT0pileup(kFALSE),
  fT0sattelite(kFALSE),
  fTracks(0x0),
  fCandidates(0x0),
  fNFMDchannels(),
  fFMDtotalMult(),
  fFMD1(0x0),
  fFMD2I(0x0),
  fFMD2O(0x0),
  fFMD3I(0x0),
  fFMD3O(0x0),
  fNCaloClusters(0),
  fCaloClusters(0x0)
//fFMDMult()
{
  //
  // Constructor
  //
  for(Int_t i=0; i<2; ++i) fIRIntClosestIntMap[i] = 0;
  for(Int_t i=0; i<3; ++i) {fVtx[i]=-999.; fVtxTPC[i]=-999.;}
  for(Int_t i=0; i<4; ++i) fCentrality[i]=-1.;
  fNV0candidates[0]=0; fNV0candidates[1]=0;
  fNtracks[0]=0; fNtracks[1]=0;
  for(Int_t i=0; i<32; ++i) fSPDntrackletsEta[i]=0;
  for(Int_t i=0; i<32; ++i) fNtracksPerTrackingFlag[i]=0;
  for(Int_t i=0; i<64; ++i) fVZEROMult[i] = 0.0;
  for(Int_t i=0; i<10; ++i) fZDCnEnergy[i]=0.0;
  for(Int_t i=0; i<10; ++i) fZDCpEnergy[i]=0.0;
  for(Int_t i=0; i<26; ++i) fT0amplitude[i]=0.0;
  for(Int_t i=0; i<3; ++i)  fT0TOF[i]=0.0;
  for(Int_t i=0; i<3; ++i)  fT0TOFbest[i]=0.0;
  for(Int_t i=0; i<5; ++i)  fNFMDchannels[i]=0.0;
  for(Int_t i=0; i<5; ++i)  fFMDtotalMult[i]=0.0;
  
  if(!fgTracks) fgTracks = new TClonesArray("AliReducedTrack", 100000);
  fTracks = fgTracks;
  if(!fgCandidates) fgCandidates = new TClonesArray("AliReducedPair", 100000);
  fCandidates = fgCandidates;
  if(!fgCaloClusters) fgCaloClusters = new TClonesArray("AliReducedCaloCluster", 50000);
  fCaloClusters = fgCaloClusters;
  if(!fgFMD1) fgFMD1 = new TClonesArray("AliReducedFMD", 10240);
  fFMD1 = fgFMD1;
  if(!fgFMD2I) fgFMD2I = new TClonesArray("AliReducedFMD", 10240);
  fFMD2I = fgFMD2I;
  if(!fgFMD2O) fgFMD2O = new TClonesArray("AliReducedFMD", 10240);
  fFMD2O = fgFMD2O;
  if(!fgFMD3I) fgFMD3I = new TClonesArray("AliReducedFMD", 10240);
  fFMD3I = fgFMD3I;
  if(!fgFMD3O) fgFMD3O = new TClonesArray("AliReducedFMD", 10240);
  fFMD3O = fgFMD3O;
}


//____________________________________________________________________________
AliReducedEvent::~AliReducedEvent()
{
  //
  // De-Constructor
  //
  //ClearEvent();
}


//____________________________________________________________________________
Float_t AliReducedEvent::MultVZEROA() const
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
Float_t AliReducedEvent::MultVZEROC() const
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
Float_t AliReducedEvent::MultVZERO() const
{
  //
  // Total VZERO multiplicity
  //
  return MultVZEROA()+MultVZEROC();
}


//____________________________________________________________________________
Float_t AliReducedEvent::MultRingVZEROA(Int_t ring) const 
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
Float_t AliReducedEvent::MultRingVZEROC(Int_t ring) const 
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

//____________________________________________________________________________
Float_t AliReducedEvent::AmplitudeTZEROA() const
{
  //
  // Total TZERO multiplicity in A side
  //
  Float_t mult=0.0;
  for(Int_t i=12;i<24;++i)
    mult += fT0amplitude[i];
  return mult;
}


//____________________________________________________________________________
Float_t AliReducedEvent::AmplitudeTZEROC() const
{
  //
  // Total TZERO multiplicity in C side
  //
  Float_t mult=0.0;
  for(Int_t i=0;i<12;++i)
    mult += fT0amplitude[i];
  return mult;
}


//____________________________________________________________________________
Float_t AliReducedEvent::AmplitudeTZERO() const
{
  //
  // Total TZERO multiplicity
  //
  return AmplitudeTZEROA()+AmplitudeTZEROC();
}



//____________________________________________________________________________
Float_t AliReducedEvent::EnergyZDCA() const
{
  //
  // Total ZDC energy in A side
  //
  Float_t energy=0.0;
  for(Int_t i=6;i<10;++i){
      if(fZDCnEnergy[i]>0.) {
  Float_t zdcNenergyAlpha = TMath::Power(fZDCnEnergy[i], gZdcNalpha);
    energy += zdcNenergyAlpha;}
  }
  return energy;
}


//____________________________________________________________________________
Float_t AliReducedEvent::EnergyZDCC() const
{
  //
  // Total ZDC energy in C side
  //
  Float_t energy=0.0;
  for(Int_t i=1;i<5;++i){
      if(fZDCnEnergy[i]>0.) {
    Float_t zdcNenergyAlpha = TMath::Power(fZDCnEnergy[i], gZdcNalpha);
    energy += zdcNenergyAlpha;}
  }
  return energy;
}


//____________________________________________________________________________
Float_t AliReducedEvent::EnergyZDC() const
{
  //
  // Total ZDC energy   
  //
  return EnergyZDCA()+EnergyZDCC();
}


//____________________________________________________________________________
Float_t AliReducedEvent::EnergyZDCn(Int_t ch) const
{
  //
  // Total ZDC energy in channel
  //
  Float_t energy=0;
  if(ch<0 || ch>9) return -999.;
      if(fZDCnEnergy[ch]>0.) {
  			energy = TMath::Power(fZDCnEnergy[ch], gZdcNalpha);
		}
  return energy;

}


//_____________________________________________________________________________
void AliReducedEvent::ClearEvent() {
  //
  // clear the event
  //
  if(fTracks) fTracks->Clear("C");
  if(fCandidates) fCandidates->Clear("C");
  if(fCaloClusters) fCaloClusters->Clear("C");
  if(fFMD1) fFMD1->Clear("C");
  if(fFMD2I) fFMD2I->Clear("C");
  if(fFMD2O) fFMD2O->Clear("C");
  if(fFMD3I) fFMD3I->Clear("C");
  if(fFMD3O) fFMD3O->Clear("C");
  fEventTag = 0;
  fEventNumberInFile = -999;
  fL0TriggerInputs=0;
  fL1TriggerInputs=0;
  fL2TriggerInputs=0;
  fIRIntClosestIntMap[0] = 0; fIRIntClosestIntMap[1] = 0;
  fRunNo = 0;
  fBC = 0;
  fTimeStamp = 0;
  fEventType = 0;
  fTriggerMask = 0;
  fIsPhysicsSelection = kTRUE;
  fIsSPDPileup = kFALSE;
  fIsSPDPileupMultBins = kFALSE;
  fIsFMDReduced = kTRUE;
  fNVtxContributors = 0;
  fNVtxTPCContributors = 0;
  fCentQuality = 0;
  for(Int_t i=0;i<4;++i) fCentrality[i] = -1.0;
  fNV0candidates[0] = 0; fNV0candidates[1] = 0;
  fNpileupSPD=0;
  fNpileupTracks=0;
  fNPMDtracks=0;
  fNTRDtracks=0;
  fNTRDtracklets=0;
  fNDielectronCandidates = 0;
  fNtracks[0] = 0; fNtracks[1] = 0;
  for(Int_t i=0; i<32; ++i) fSPDntrackletsEta[i] = 0;
  for(Int_t i=0; i<32; ++i) fNtracksPerTrackingFlag[i] = 0;
  for(Int_t i=0; i<3; ++i) {fVtx[i]=-999.; fVtxTPC[i]=-999.;}
  for(Int_t i=0; i<64; ++i) fVZEROMult[i] = 0.0;
  for(Int_t i=0; i<10; ++i) fZDCnEnergy[i]=0.0;
  for(Int_t i=0; i<10; ++i) fZDCpEnergy[i]=0.0;
  for(Int_t i=0; i<26; ++i) fT0amplitude[i]=0.0;
  for(Int_t i=0; i<3; ++i)  fT0TOF[i]=0.0;
  for(Int_t i=0; i<3; ++i)  fT0TOFbest[i]=0.0;
  fT0pileup = kFALSE;
  fT0zVertex = -999.;
  fT0start = -999.;
  fT0sattelite = kFALSE;
  for(Int_t i=0; i<5; ++i)  fNFMDchannels[i]=0.0;
  for(Int_t i=0; i<5; ++i)  fFMDtotalMult[i]=0.0;
}


//_______________________________________________________________________________
AliReducedEventFriend::AliReducedEventFriend() :
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
AliReducedEventFriend::~AliReducedEventFriend()
{
  //
  // De-Constructor
  //
  ClearEvent();
}


//_____________________________________________________________________________
void AliReducedEventFriend::ClearEvent() {
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
void AliReducedEventFriend::CopyEvent(const AliReducedEventFriend* event) {
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
AliReducedCaloCluster::AliReducedCaloCluster() :
 fType(kUndefined),
 fEnergy(-999.),
 fTrackDx(-999.),
 fTrackDz(-999.),
 fM20(-999.),
 fM02(-999.),
 fDispersion(-999.)
{
  //
  // default constructor
  //
}


//_____________________________________________________________________________
AliReducedCaloCluster::~AliReducedCaloCluster()
{
  //
  // destructor
  //
}


//_______________________________________________________________________________
void AliReducedEvent::GetQvector(Double_t Qvec[][2], Int_t det,
                                 Float_t etaMin/*=-0.8*/, Float_t etaMax/*=+0.8*/,
				 Bool_t (*IsTrackSelected)(AliReducedTrack*)/*=NULL*/) {
  //
  // Get the event plane for a specified detector
  //
  if(det==AliReducedEventFriend::kTPC || 
     det==AliReducedEventFriend::kTPCptWeights ||
     det==AliReducedEventFriend::kTPCpos ||
     det==AliReducedEventFriend::kTPCneg) {
    GetTPCQvector(Qvec, det, etaMin, etaMax, IsTrackSelected);
    return;
  }
  if(det==AliReducedEventFriend::kVZEROA ||
     det==AliReducedEventFriend::kVZEROC) {
    GetVZEROQvector(Qvec, det);   
    return;
  }
  if(det==AliReducedEventFriend::kZDCA ||
     det==AliReducedEventFriend::kZDCC) {
    GetZDCQvector(Qvec, det);
    return;
  }
  if(det==AliReducedEventFriend::kFMD) {
    //TODO implementation
    return;
  }
  return;
}


//_________________________________________________________________
Int_t AliReducedEvent::GetTPCQvector(Double_t Qvec[][2], Int_t det, 
                                     Float_t etaMin/*=-0.8*/, Float_t etaMax/*=+0.8*/,
				     Bool_t (*IsTrackSelected)(AliReducedTrack*)/*=NULL*/) {
  //
  // Construct the event plane using tracks in the barrel
  //
  if(!(det==AliReducedEventFriend::kTPC ||
       det==AliReducedEventFriend::kTPCpos ||
       det==AliReducedEventFriend::kTPCneg))
    return 0;
  Int_t nUsedTracks = 0;
  Short_t charge = 0;
  AliReducedTrack* track=0x0;
  Double_t weight=0.0; Double_t absWeight = 0.0; Double_t x=0.0; Double_t y=0.0; 
  TIter nextTrack(fTracks);
  while((track=static_cast<AliReducedTrack*>(nextTrack()))) {
    if(track->Eta()<etaMin) continue;
    if(track->Eta()>etaMax) continue;
    charge = track->Charge();
    if(det==AliReducedEventFriend::kTPCpos && charge<0) continue;
    if(det==AliReducedEventFriend::kTPCneg && charge>0) continue;
    
    if(IsTrackSelected && !IsTrackSelected(track)) continue;
    absWeight = 1.0;
    if(det==AliReducedEventFriend::kTPCptWeights) {
      absWeight = track->Pt();
      if(absWeight>2.0) absWeight = 2.0;    // pt is the weight used for the event plane
    }
    weight = absWeight;
        
    ++nUsedTracks;
    x = TMath::Cos(track->Phi());
    y = TMath::Sin(track->Phi());
    //  1st harmonic  
    Qvec[0][0] += weight*x;
    Qvec[0][1] += weight*y;
    //  2nd harmonic
    Qvec[1][0] += absWeight*(2.0*TMath::Power(x,2.0)-1);
    Qvec[1][1] += absWeight*(2.0*x*y);
    //  3rd harmonic
    Qvec[2][0] += weight*(4.0*TMath::Power(x,3.0)-3.0*x);
    Qvec[2][1] += weight*(3.0*y-4.0*TMath::Power(y,3.0));
    //  4th harmonic
    Qvec[3][0] += absWeight*(1.0-8.0*TMath::Power(x*y,2.0));
    Qvec[3][1] += absWeight*(4.0*x*y-8.0*x*TMath::Power(y,3.0));
    //  5th harmonic
    Qvec[4][0] += weight*(16.0*TMath::Power(x,5.0)-20.0*TMath::Power(x, 3.0)+5.0*x);
    Qvec[4][1] += weight*(16.0*TMath::Power(y,5.0)-20.0*TMath::Power(y, 3.0)+5.0*y);
    //  6th harmonic
    Qvec[5][0] += absWeight*(32.0*TMath::Power(x,6.0)-48.0*TMath::Power(x, 4.0)+18.0*TMath::Power(x, 2.0)-1.0);
    Qvec[5][1] += absWeight*(x*y*(32.0*TMath::Power(y,4.0)-32.0*TMath::Power(y, 2.0)+6.0)); 
  }
  return nUsedTracks;
}


//____________________________________________________________________________
void AliReducedEvent::SubtractParticleFromQvector(
	AliReducedTrack* particle, Double_t Qvec[][2], Int_t det, 
        Float_t etaMin/*=-0.8*/, Float_t etaMax/*=+0.8*/,
	Bool_t (*IsTrackSelected)(AliReducedTrack*)/*=NULL*/) {
  //
  // subtract a particle from the event Q-vector
  //
  Float_t eta = particle->Eta();
  if(eta<etaMin) return;
  if(eta>etaMax) return;
  
  Float_t charge = particle->Charge();
  if(det==AliReducedEventFriend::kTPCpos && charge<0) return;
  if(det==AliReducedEventFriend::kTPCneg && charge>0) return;
  
  if(IsTrackSelected && !IsTrackSelected(particle)) return;
  
  Double_t weight=0.0; Double_t absWeight = 0.0;
  if(det==AliReducedEventFriend::kTPCptWeights) {
    absWeight = particle->Pt();
    if(absWeight>2.0) absWeight = 2.0;
  }
  weight = absWeight;
  //  if(eta<0.0) weight *= -1.0;
    
  Float_t x = TMath::Cos(particle->Phi());
  Float_t y = TMath::Sin(particle->Phi());
  
  //  1st harmonic  
  Qvec[0][0] -= weight*x;
  Qvec[0][1] -= weight*y;
  //  2nd harmonic
  Qvec[1][0] -= absWeight*(2.0*TMath::Power(x,2.0)-1);
  Qvec[1][1] -= absWeight*(2.0*x*y);
  //  3rd harmonic
  Qvec[2][0] -= weight*(4.0*TMath::Power(x,3.0)-3.0*x);
  Qvec[2][1] -= weight*(3.0*y-4.0*TMath::Power(y,3.0));
  //  4th harmonic
  Qvec[3][0] -= absWeight*(1.0-8.0*TMath::Power(x*y,2.0));
  Qvec[3][1] -= absWeight*(4.0*x*y-8.0*x*TMath::Power(y,3.0));
  //  5th harmonic
  Qvec[4][0] -= weight*(16.0*TMath::Power(x,5.0)-20.0*TMath::Power(x, 3.0)+5.0*x);
  Qvec[4][1] -= weight*(16.0*TMath::Power(y,5.0)-20.0*TMath::Power(y, 3.0)+5.0*y);
  //  6th harmonic
  Qvec[5][0] -= absWeight*(32.0*TMath::Power(x,6.0)-48.0*TMath::Power(x, 4.0)+18.0*TMath::Power(x, 2.0)-1.0);
  Qvec[5][1] -= absWeight*(x*y*(32.0*TMath::Power(y,4.0)-32.0*TMath::Power(y, 2.0)+6.0)); 
  
  return;
}


//_________________________________________________________________
void AliReducedEvent::GetVZEROQvector(Double_t Qvec[][2], Int_t det) {
  //
  // Get the reaction plane from the VZERO detector for a given harmonic
  //
  GetVZEROQvector(Qvec, det, fVZEROMult);
}


//_________________________________________________________________
void AliReducedEvent::GetVZEROQvector(Double_t Qvec[][2], Int_t det, Float_t* vzeroMult) {
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
  if(!(det==AliReducedEventFriend::kVZEROA ||
       det==AliReducedEventFriend::kVZEROC))
    return; 
  
  const Double_t kX[8] = {0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268, 0.38268, 0.92388};    // cosines of the angles of the VZERO sectors (8 per ring)
  const Double_t kY[8] = {0.38268, 0.92388, 0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268};    // sines     -- " --
  Int_t phi;
  
  for(Int_t iChannel=0; iChannel<64; ++iChannel) {
    if(iChannel<32 && det==AliReducedEventFriend::kVZEROA) continue;
    if(iChannel>=32 && det==AliReducedEventFriend::kVZEROC) continue;
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
void AliReducedEvent::GetZDCQvector(Double_t Qvec[][2], Int_t det) const {
  //
  // Get the reaction plane from the ZDC detector for a given harmonic
  //
  GetZDCQvector(Qvec, det, fZDCnEnergy);
}


//_________________________________________________________________
void AliReducedEvent::GetZDCQvector(Double_t Qvec[][2], Int_t det, const Float_t* zdcEnergy) const {
  //
  // Construct the event plane using the ZDC
  // ZDC has 2 side (A and C) with 4 calorimeters on each side  
  // The XY position of each calorimeter is specified by the 
  // zdcNTowerCenters_x and zdcNTowerCenters_y arrays
  if(!(det==AliReducedEventFriend::kZDCA ||
       det==AliReducedEventFriend::kZDCC ))
    return; 
 

  //Int_t minTower,maxTower;
  //Float_t totalEnergy = 0.0;

  const Float_t zdcTowerCenter = 1.75;


  Float_t zdcNTowerCentersX[4] = {0.0};
  Float_t zdcNTowerCentersY[4] = {0.0};


  if(det==AliReducedEventFriend::kZDCA){
    zdcNTowerCentersX[0] =  zdcTowerCenter; zdcNTowerCentersX[1] =-zdcTowerCenter;
    zdcNTowerCentersX[2] =  zdcTowerCenter; zdcNTowerCentersX[3] =-zdcTowerCenter;
  }

  if(det==AliReducedEventFriend::kZDCC){
    zdcNTowerCentersX[0] = -zdcTowerCenter; zdcNTowerCentersX[1] = zdcTowerCenter;
    zdcNTowerCentersX[2] = -zdcTowerCenter; zdcNTowerCentersX[3] = zdcTowerCenter;
  }

  zdcNTowerCentersY[0] = -zdcTowerCenter; zdcNTowerCentersY[1] =-zdcTowerCenter;
  zdcNTowerCentersY[2] =  zdcTowerCenter; zdcNTowerCentersY[3] = zdcTowerCenter;

 
  for(Int_t ih=0; ih<6; ih++) {Qvec[ih][0] = 0.0; Qvec[ih][1] = 0.0;}   // harmonic Q-vector
  Float_t zdcNCentroidSum = 0;
    
  for(Int_t iTow=(det==AliReducedEventFriend::kZDCA ? 6 : 1); iTow<(det==AliReducedEventFriend::kZDCA ? 10 : 5); ++iTow) {
    //if(fZDCnEnergy[i+(det==AliReducedEventFriend::kZDCA ? 4 : 0)]>0.0) {
    //  1st harmonic  
    Qvec[0][0] += zdcEnergy[iTow]*zdcNTowerCentersX[(iTow-1)%5];
    Qvec[0][1] += zdcEnergy[iTow]*zdcNTowerCentersY[(iTow-1)%5];
    //  2nd harmonic
    Qvec[1][0] += zdcEnergy[iTow]*(2.0*TMath::Power(zdcNTowerCentersX[(iTow-1)%5],2.0)-1);
    Qvec[1][1] += zdcEnergy[iTow]*(2.0*zdcNTowerCentersX[(iTow-1)%5]*zdcNTowerCentersY[(iTow-1)%5]);
    //  3rd harmonic
    Qvec[2][0] += zdcEnergy[iTow]*(4.0*TMath::Power(zdcNTowerCentersX[(iTow-1)%5],3.0)-3.0*zdcNTowerCentersX[(iTow-1)%5]);
    Qvec[2][1] += zdcEnergy[iTow]*(3.0*zdcNTowerCentersY[(iTow-1)%5]-4.0*TMath::Power(zdcNTowerCentersY[(iTow-1)%5],3.0));
    //  4th harmonic
    Qvec[3][0] += zdcEnergy[iTow]*(1.0-8.0*TMath::Power(zdcNTowerCentersX[(iTow-1)%5]*zdcNTowerCentersY[(iTow-1)%5],2.0));
    Qvec[3][1] += zdcEnergy[iTow]*(4.0*zdcNTowerCentersX[(iTow-1)%5]*zdcNTowerCentersY[(iTow-1)%5]-8.0*zdcNTowerCentersX[(iTow-1)%5]*TMath::Power(zdcNTowerCentersY[(iTow-1)%5],3.0));
    //  5th harmonic
    Qvec[4][0] += zdcEnergy[iTow]*(16.0*TMath::Power(zdcNTowerCentersX[(iTow-1)%5],5.0)-20.0*TMath::Power(zdcNTowerCentersX[(iTow-1)%5], 3.0)+5.0*zdcNTowerCentersX[(iTow-1)%5]);
    Qvec[4][1] += zdcEnergy[iTow]*(16.0*TMath::Power(zdcNTowerCentersY[(iTow-1)%5],5.0)-20.0*TMath::Power(zdcNTowerCentersY[(iTow-1)%5], 3.0)+5.0*zdcNTowerCentersY[(iTow-1)%5]);
    //  6th harmonic
    Qvec[5][0] += zdcEnergy[iTow]*(32.0*TMath::Power(zdcNTowerCentersX[(iTow-1)%5],6.0)-48.0*TMath::Power(zdcNTowerCentersX[(iTow-1)%5], 4.0)+18.0*TMath::Power(zdcNTowerCentersX[(iTow-1)%5], 2.0)-1.0);
    Qvec[5][1] += zdcEnergy[iTow]*(zdcNTowerCentersX[(iTow-1)%5]*zdcNTowerCentersY[(iTow-1)%5]*(32.0*TMath::Power(zdcNTowerCentersY[(iTow-1)%5],4.0)-32.0*TMath::Power(zdcNTowerCentersY[(iTow-1)%5], 2.0)+6.0));

    zdcNCentroidSum += zdcEnergy[iTow];

    if(zdcNCentroidSum>0.0) {
      Qvec[0][0] /= zdcNCentroidSum;
      Qvec[0][1] /= zdcNCentroidSum;
      Qvec[1][0] /= zdcNCentroidSum;
      Qvec[1][1] /= zdcNCentroidSum;
      Qvec[2][0] /= zdcNCentroidSum;
      Qvec[2][1] /= zdcNCentroidSum;
      Qvec[3][0] /= zdcNCentroidSum;
      Qvec[3][1] /= zdcNCentroidSum;
      Qvec[4][0] /= zdcNCentroidSum;
      Qvec[4][1] /= zdcNCentroidSum;
      Qvec[5][0] /= zdcNCentroidSum;
      Qvec[5][1] /= zdcNCentroidSum;
    }
  }   // end loop over channels
  return;
}
