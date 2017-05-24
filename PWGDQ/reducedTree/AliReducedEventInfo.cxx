/*
***********************************************************
  Implementation of AliReducedEventInfo class.
  Contact: iarsene@cern.ch
  2015/04/15
  *********************************************************
*/

#ifndef ALIREDUCEDEVENTINFO_H
#include "AliReducedEventInfo.h"
#endif

#include <iostream>
#include "AliReducedBaseEvent.h"
#include "AliReducedCaloClusterInfo.h"
#include "AliReducedFMDInfo.h"
#include "AliReducedPairInfo.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedEventPlaneInfo.h"

ClassImp(AliReducedEventInfo)

TClonesArray* AliReducedEventInfo::fgCaloClusters = 0;
const Float_t AliReducedEventInfo::fgkZdcNalpha = 1.0;
TClonesArray* AliReducedEventInfo::fgFMD = 0;

//____________________________________________________________________________
AliReducedEventInfo::AliReducedEventInfo() :
  AliReducedBaseEvent(),
  fEventNumberInFile(0),
  fL0TriggerInputs(0),
  fL1TriggerInputs(0),
  fL2TriggerInputs(0),
  fBC(0),
  fTimeStamp(0),
  fEventType(0),
  fTriggerMask(0),
  fMultiplicityEstimators(),
  fMultiplicityEstimatorPercentiles(),
  fIsPhysicsSelection(kTRUE),
  fIsSPDPileup(kFALSE),
  fIsSPDPileupMultBins(kFALSE),
  fIRIntClosestIntMap(),
  fVtxCovMatrix(),
  fVtxTPC(),
  fNVtxTPCContributors(0),
  fVtxSPD(),
  fNVtxSPDContributors(0),
  fNpileupSPD(0),
  fNpileupTracks(0),
  fNTPCclusters(0),
  fNPMDtracks(0),
  fNTRDtracks(0),
  fNTRDtracklets(0),
  fSPDntracklets(0),
  fSPDnSingle(0),
  fVZEROMult(),
  fVZEROTotalMult(),
  fZDCnEnergy(),
  fZDCpEnergy(),
  fZDCnTotalEnergy(),
  fZDCpTotalEnergy(),
  fT0amplitude(),
  fT0TOF(),
  fT0TOFbest(),
  fT0zVertex(0),
  fT0start(0),
  fT0pileup(kFALSE),
  fT0sattelite(kFALSE),
  fNCaloClusters(0),
  fCaloClusters(0x0),
  fFMD(0x0),
  //fEventPlane(0x0)
  fEventPlane()
{
  //
  // Constructor
  //
  for(Int_t i=0; i<2; ++i) fIRIntClosestIntMap[i] = 0;
  for(Int_t i=0; i<6; ++i) fVtxCovMatrix[i]=0.;
  for(Int_t i=0; i<3; ++i) fVtxTPC[i]=-999.;
  for(Int_t i=0; i<3; ++i) fVtxSPD[i]=-999.;
  for(Int_t i=0; i<10; ++i) fMultiplicityEstimators[i]=-999.;
  for(Int_t i=0; i<10; ++i) fMultiplicityEstimatorPercentiles[i]=-999.;
  for(Int_t i=0; i<32; ++i) fSPDntrackletsEta[i]=0;
  for(Int_t i=0; i<2; ++i) fSPDFiredChips[i]=0;
  for(Int_t i=0; i<6; ++i) fITSClusters[i]=0;
  for(Int_t i=0; i<32; ++i) fNtracksPerTrackingFlag[i]=0;
  for(Int_t i=0; i<64; ++i) fVZEROMult[i] = 0.0;
  for(Int_t i=0; i<2; ++i) fVZEROTotalMult[i] = 0.0;
  for(Int_t i=0; i<10; ++i) fZDCnEnergy[i]=0.0;
  for(Int_t i=0; i<10; ++i) fZDCpEnergy[i]=0.0;
  for(Int_t i=0; i<2; ++i) fZDCnTotalEnergy[i]=0.0;
  for(Int_t i=0; i<2; ++i) fZDCpTotalEnergy[i]=0.0;
  for(Int_t i=0; i<26; ++i) fT0amplitude[i]=0.0;
  for(Int_t i=0; i<3; ++i)  fT0TOF[i]=0.0;
  for(Int_t i=0; i<3; ++i)  fT0TOFbest[i]=0.0;
}


//____________________________________________________________________________
AliReducedEventInfo::AliReducedEventInfo(const Char_t* name, Int_t trackOption /*=kNoInit*/) :
  AliReducedBaseEvent(name, trackOption),
  fEventNumberInFile(0),
  fL0TriggerInputs(0),
  fL1TriggerInputs(0),
  fL2TriggerInputs(0),
  fBC(0),
  fTimeStamp(0),
  fEventType(0),
  fTriggerMask(0),
  fMultiplicityEstimators(),
  fMultiplicityEstimatorPercentiles(),
  fIsPhysicsSelection(kTRUE),
  fIsSPDPileup(kFALSE),
  fIsSPDPileupMultBins(kFALSE),
  fIRIntClosestIntMap(),
  fVtxTPC(),
  fNVtxTPCContributors(0),
  fVtxSPD(),
  fNVtxSPDContributors(0),
  fNpileupSPD(0),
  fNpileupTracks(0),
  fNTPCclusters(0),
  fNPMDtracks(0),
  fNTRDtracks(0),
  fNTRDtracklets(0),
  fSPDntracklets(0),
  fSPDnSingle(0),
  fVZEROMult(),
  fVZEROTotalMult(),
  fZDCnEnergy(),
  fZDCpEnergy(),
  fZDCnTotalEnergy(),
  fZDCpTotalEnergy(),
  fT0amplitude(),
  fT0TOF(),
  fT0TOFbest(),
  fT0zVertex(0),
  fT0start(0),
  fT0pileup(kFALSE),
  fT0sattelite(kFALSE),
  fNCaloClusters(0),
  fCaloClusters(0x0),
  fFMD(0x0),
  //fEventPlane(0x0)
  fEventPlane()
{
  //
  // Constructor
  //
  for(Int_t i=0; i<2; ++i) fIRIntClosestIntMap[i] = 0;
  for(Int_t i=0; i<3; ++i) fVtxTPC[i]=-999.;
  for(Int_t i=0; i<3; ++i) fVtxSPD[i]=-999.;
  for(Int_t i=0; i<10; ++i) fMultiplicityEstimators[i]=-999.;
  for(Int_t i=0; i<10; ++i) fMultiplicityEstimatorPercentiles[i]=-999.;
  for(Int_t i=0; i<32; ++i) fSPDntrackletsEta[i]=0;
  for(Int_t i=0; i<2; ++i) fSPDFiredChips[i]=0;
  for(Int_t i=0; i<6; ++i) fITSClusters[i]=0;
  for(Int_t i=0; i<32; ++i) fNtracksPerTrackingFlag[i]=0;
  for(Int_t i=0; i<64; ++i) fVZEROMult[i] = 0.0;
  for(Int_t i=0; i<2; ++i) fVZEROTotalMult[i] = 0.0;
  for(Int_t i=0; i<10; ++i) fZDCnEnergy[i]=0.0;
  for(Int_t i=0; i<10; ++i) fZDCpEnergy[i]=0.0;
  for(Int_t i=0; i<2; ++i) fZDCnTotalEnergy[i]=0.0;
  for(Int_t i=0; i<2; ++i) fZDCpTotalEnergy[i]=0.0;
  for(Int_t i=0; i<26; ++i) fT0amplitude[i]=0.0;
  for(Int_t i=0; i<3; ++i)  fT0TOF[i]=0.0;
  for(Int_t i=0; i<3; ++i)  fT0TOFbest[i]=0.0;
  
  if(!fgCaloClusters) fgCaloClusters = new TClonesArray("AliReducedCaloClusterInfo", 50000);
  fCaloClusters = fgCaloClusters;
  if(!fgFMD) fgFMD = new TClonesArray("AliReducedFMDInfo", 40000);
  fFMD = fgFMD;
}


//____________________________________________________________________________
AliReducedEventInfo::~AliReducedEventInfo()
{
  //
  // De-Constructor
  //
}

//_____________________________________________________________________________
void AliReducedEventInfo::ClearEvent() {
  //
  // clear the event
  //
  AliReducedBaseEvent::ClearEvent();
  
  if(fCaloClusters) fCaloClusters->Clear("C");
  if(fFMD) fFMD->Clear("C");
  fEventNumberInFile = -999;
  fL0TriggerInputs=0;
  fL1TriggerInputs=0;
  fL2TriggerInputs=0;
  fIRIntClosestIntMap[0] = 0; fIRIntClosestIntMap[1] = 0;
  fBC = 0;
  fTimeStamp = 0;
  fEventType = 0;
  fTriggerMask = 0;
  fIsPhysicsSelection = kTRUE;
  fIsSPDPileup = kFALSE;
  fIsSPDPileupMultBins = kFALSE;
  fNVtxTPCContributors = 0;
  fNVtxSPDContributors = 0;
  fNpileupSPD=0;
  fNpileupTracks=0;
  fNPMDtracks=0;
  fNTRDtracks=0;
  fNTRDtracklets=0;
  fNtracks[0] = 0; fNtracks[1] = 0;
  for(Int_t i=0; i<32; ++i) fSPDntrackletsEta[i] = 0;
  fSPDnSingle = 0;
  for(Int_t i=0; i<2; ++i) fSPDFiredChips[i]=0;
  for(Int_t i=0; i<6; ++i) fITSClusters[i]=0;
  for(Int_t i=0; i<32; ++i) fNtracksPerTrackingFlag[i] = 0;
  for(Int_t i=0; i<3; ++i) fVtxTPC[i]=-999.;
  for(Int_t i=0; i<3; ++i) fVtxSPD[i]=-999.;
  for(Int_t i=0; i<10; ++i) fMultiplicityEstimators[i]=-999.;
  for(Int_t i=0; i<10; ++i) fMultiplicityEstimatorPercentiles[i]=-999.;
  for(Int_t i=0; i<64; ++i) fVZEROMult[i] = 0.0;
  for(Int_t i=0; i<2; ++i) fVZEROTotalMult[i] = 0.0;
  for(Int_t i=0; i<10; ++i) fZDCnEnergy[i]=0.0;
  for(Int_t i=0; i<10; ++i) fZDCpEnergy[i]=0.0;
  for(Int_t i=0; i<2; ++i) fZDCnTotalEnergy[i]=0.0;
  for(Int_t i=0; i<2; ++i) fZDCpTotalEnergy[i]=0.0;
  for(Int_t i=0; i<26; ++i) fT0amplitude[i]=0.0;
  for(Int_t i=0; i<3; ++i)  fT0TOF[i]=0.0;
  for(Int_t i=0; i<3; ++i)  fT0TOFbest[i]=0.0;
  fT0pileup = kFALSE;
  fT0zVertex = -999.;
  fT0start = -999.;
  fT0sattelite = kFALSE;
}

//_______________________________________________________________________________
void AliReducedEventInfo::GetQvector(Double_t Qvec[][2], Int_t det,
                                 Float_t etaMin/*=-0.8*/, Float_t etaMax/*=+0.8*/,
				 Bool_t (*IsTrackSelected)(AliReducedTrackInfo*)/*=NULL*/) {
  //
  // Get the event plane for a specified detector
  //
  if(det==AliReducedEventPlaneInfo::kTPC || 
     det==AliReducedEventPlaneInfo::kTPCptWeights ||
     det==AliReducedEventPlaneInfo::kTPCpos ||
     det==AliReducedEventPlaneInfo::kTPCneg) {
    GetTPCQvector(Qvec, det, etaMin, etaMax, IsTrackSelected);
    return;
  }
  if(det==AliReducedEventPlaneInfo::kVZEROA ||
     det==AliReducedEventPlaneInfo::kVZEROC) {
    GetVZEROQvector(Qvec, det);   
    return;
  }
  if(det==AliReducedEventPlaneInfo::kZDCA ||
     det==AliReducedEventPlaneInfo::kZDCC) {
    GetZDCQvector(Qvec, det);
    return;
  }
  if(det==AliReducedEventPlaneInfo::kFMD) {
    //TODO implementation
    return;
  }
  return;
}


//_________________________________________________________________
Int_t AliReducedEventInfo::GetTPCQvector(Double_t Qvec[][2], Int_t det, 
                                     Float_t etaMin/*=-0.8*/, Float_t etaMax/*=+0.8*/,
				     Bool_t (*IsTrackSelected)(AliReducedTrackInfo*)/*=NULL*/) {
  //
  // Construct the event plane using tracks in the barrel
  //
  if(!(det==AliReducedEventPlaneInfo::kTPC ||
       det==AliReducedEventPlaneInfo::kTPCpos ||
       det==AliReducedEventPlaneInfo::kTPCneg))
    return 0;
  Int_t nUsedTracks = 0;
  Short_t charge = 0;
  AliReducedTrackInfo* track=0x0;
  Double_t weight=0.0; Double_t absWeight = 0.0; Double_t x=0.0; Double_t y=0.0; 
  TIter nextTrack(fTracks);
  while((track=static_cast<AliReducedTrackInfo*>(nextTrack()))) {
    if(track->Eta()<etaMin) continue;
    if(track->Eta()>etaMax) continue;
    charge = track->Charge();
    if(det==AliReducedEventPlaneInfo::kTPCpos && charge<0) continue;
    if(det==AliReducedEventPlaneInfo::kTPCneg && charge>0) continue;
    
    if(IsTrackSelected && !IsTrackSelected(track)) continue;
    absWeight = 1.0;
    if(det==AliReducedEventPlaneInfo::kTPCptWeights) {
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
void AliReducedEventInfo::SubtractParticleFromQvector(
	AliReducedTrackInfo* particle, Double_t Qvec[][2], Int_t det, 
        Float_t etaMin/*=-0.8*/, Float_t etaMax/*=+0.8*/,
	Bool_t (*IsTrackSelected)(AliReducedTrackInfo*)/*=NULL*/) {
  //
  // subtract a particle from the event Q-vector
  //
  //if(!particle) {std::cout<<"subtraction particle doesn't exist"<<std::endl; return;}
  Float_t eta = particle->Eta();
  if(eta<etaMin) return;
  if(eta>etaMax) return;
  
  Float_t charge = particle->Charge();
  if(det==AliReducedEventPlaneInfo::kTPCpos && charge<0) return;
  if(det==AliReducedEventPlaneInfo::kTPCneg && charge>0) return;
  
  if(IsTrackSelected && !IsTrackSelected(particle)) return;
  
  Double_t weight=0.0; Double_t absWeight = 0.0;
  if(det==AliReducedEventPlaneInfo::kTPCptWeights) {
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
void AliReducedEventInfo::GetVZEROQvector(Double_t Qvec[][2], Int_t det) {
  //
  // Get the reaction plane from the VZERO detector for a given harmonic
  //
  GetVZEROQvector(Qvec, det, fVZEROMult);
}


//_________________________________________________________________
void AliReducedEventInfo::GetVZEROQvector(Double_t Qvec[][2], Int_t det, Float_t* vzeroMult) {
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
  if(!(det==AliReducedEventPlaneInfo::kVZEROA ||
       det==AliReducedEventPlaneInfo::kVZEROC))
    return; 
  
  const Double_t kX[8] = {0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268, 0.38268, 0.92388};    // cosines of the angles of the VZERO sectors (8 per ring)
  const Double_t kY[8] = {0.38268, 0.92388, 0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268};    // sines     -- " --
  Int_t phi;
  
  for(Int_t iChannel=0; iChannel<64; ++iChannel) {
    if(iChannel<32 && det==AliReducedEventPlaneInfo::kVZEROA) continue;
    if(iChannel>=32 && det==AliReducedEventPlaneInfo::kVZEROC) continue;
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
void AliReducedEventInfo::GetZDCQvector(Double_t Qvec[][2], Int_t det) const {
  //
  // Get the reaction plane from the ZDC detector for a given harmonic
  //
  GetZDCQvector(Qvec, det, fZDCnEnergy);
}


//_________________________________________________________________
void AliReducedEventInfo::GetZDCQvector(Double_t Qvec[][2], Int_t det, const Float_t* zdcEnergy) const {
  //
  // Construct the event plane using the ZDC
  // ZDC has 2 side (A and C) with 4 calorimeters on each side  
  // The XY position of each calorimeter is specified by the 
  // zdcNTowerCenters_x and zdcNTowerCenters_y arrays
  if(!(det==AliReducedEventPlaneInfo::kZDCA ||
       det==AliReducedEventPlaneInfo::kZDCC ))
    return; 
 

  //Int_t minTower,maxTower;
  //Float_t totalEnergy = 0.0;

  const Float_t zdcTowerCenter = 1.75;


  Float_t zdcNTowerCentersX[4] = {0.0};
  Float_t zdcNTowerCentersY[4] = {0.0};


  if(det==AliReducedEventPlaneInfo::kZDCA){
    zdcNTowerCentersX[0] =  zdcTowerCenter; zdcNTowerCentersX[1] =-zdcTowerCenter;
    zdcNTowerCentersX[2] =  zdcTowerCenter; zdcNTowerCentersX[3] =-zdcTowerCenter;
  }

  if(det==AliReducedEventPlaneInfo::kZDCC){
    zdcNTowerCentersX[0] = -zdcTowerCenter; zdcNTowerCentersX[1] = zdcTowerCenter;
    zdcNTowerCentersX[2] = -zdcTowerCenter; zdcNTowerCentersX[3] = zdcTowerCenter;
  }

  zdcNTowerCentersY[0] = -zdcTowerCenter; zdcNTowerCentersY[1] =-zdcTowerCenter;
  zdcNTowerCentersY[2] =  zdcTowerCenter; zdcNTowerCentersY[3] = zdcTowerCenter;

 
  for(Int_t ih=0; ih<6; ih++) {Qvec[ih][0] = 0.0; Qvec[ih][1] = 0.0;}   // harmonic Q-vector
  Float_t zdcNCentroidSum = 0;
    
  for(Int_t iTow=(det==AliReducedEventPlaneInfo::kZDCA ? 6 : 1); iTow<(det==AliReducedEventPlaneInfo::kZDCA ? 10 : 5); ++iTow) {
    //if(fZDCnEnergy[i+(det==AliReducedEventPlaneInfo::kZDCA ? 4 : 0)]>0.0) {
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

//____________________________________________________________________________
Float_t AliReducedEventInfo::EnergyZDCA() const
{
  //
  // Total ZDC energy in A side
  //
  /*Float_t energy=0.0;
  for(Int_t i=6;i<10;++i){
      if(fZDCnEnergy[i]>0.) {
  Float_t zdcNenergyAlpha = TMath::Power(fZDCnEnergy[i], fgkZdcNalpha);
    energy += zdcNenergyAlpha;}
  }
  return energy;*/
  return fZDCnTotalEnergy[0];
}


//____________________________________________________________________________
Float_t AliReducedEventInfo::EnergyZDCC() const
{
  //
  // Total ZDC energy in C side
  //
  /*Float_t energy=0.0;
  for(Int_t i=1;i<5;++i){
      if(fZDCnEnergy[i]>0.) {
    Float_t zdcNenergyAlpha = TMath::Power(fZDCnEnergy[i], fgkZdcNalpha);
    energy += zdcNenergyAlpha;}
  }
  return energy;*/
  return fZDCnTotalEnergy[1];
}


//____________________________________________________________________________
Float_t AliReducedEventInfo::EnergyZDC() const
{
  //
  // Total ZDC energy   
  //
  return EnergyZDCA()+EnergyZDCC();
}


//____________________________________________________________________________
Float_t AliReducedEventInfo::EnergyZDCn(Int_t ch) const
{
  //
  // Total ZDC energy in channel
  //
  Float_t energy=0;
  if(ch<0 || ch>9) return -999.;
      if(fZDCnEnergy[ch]>0.) {
  			energy = TMath::Power(fZDCnEnergy[ch], fgkZdcNalpha);
		}
  return energy;

}

//____________________________________________________________________________
Float_t AliReducedEventInfo::MultVZEROA() const
{
  //
  // Total VZERO multiplicity in A side
  //
  /*Float_t mult=0.0;
  for(Int_t i=32;i<64;++i)
    mult += fVZEROMult[i];
  return mult;*/
  return fVZEROTotalMult[0];
}


//____________________________________________________________________________
Float_t AliReducedEventInfo::MultVZEROC() const
{
  //
  // Total VZERO multiplicity in C side
  //
  /*Float_t mult=0.0;
  for(Int_t i=0;i<32;++i)
    mult += fVZEROMult[i];
  return mult;*/
  return fVZEROTotalMult[1];
}


//____________________________________________________________________________
Float_t AliReducedEventInfo::MultVZERO() const
{
  //
  // Total VZERO multiplicity
  //
  return MultVZEROA()+MultVZEROC();
}


//____________________________________________________________________________
Float_t AliReducedEventInfo::MultRingVZEROA(Int_t ring) const 
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
Float_t AliReducedEventInfo::MultRingVZEROC(Int_t ring) const 
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
Float_t AliReducedEventInfo::AmplitudeTZEROA() const
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
Float_t AliReducedEventInfo::AmplitudeTZEROC() const
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
Float_t AliReducedEventInfo::AmplitudeTZERO() const
{
  //
  // Total TZERO multiplicity
  //
  return AmplitudeTZEROA()+AmplitudeTZEROC();
}

//_________________________________________________________________________________
Double_t AliReducedEventInfo::GetQvectorFMD(Int_t c, Double_t etamin, Double_t etamax){
  Double_t q=0.0;
  AliReducedFMDInfo* fmd = new AliReducedFMDInfo();
  TIter nextFMD(fFMD);
  for(Int_t it=0; it<fFMD->GetEntriesFast(); ++it) {
    fmd = (AliReducedFMDInfo*)nextFMD();
    if(!fmd) continue;
    if(fmd->Eta()>etamin&&fmd->Eta()<etamax){
      if(c==0) q+=fmd->Multiplicity()*TMath::Cos(2.*fmd->Phi());
      if(c==1) q+=fmd->Multiplicity()*TMath::Sin(2.*fmd->Phi());
    }
  }
  return q;
}

