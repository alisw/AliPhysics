/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//
// Dielectron helper functions wrapped in a namespace
// 
//
// Authors: 
//   Jens Wiechula <Jens.Wiechula@cern.ch> 
//   Frederick Kramer <Frederick.Kramer@cern.ch> 
// 




#include <TError.h>
#include <TMath.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TVectorD.h>
#include <TF1.h>
#include <TRandom.h>
#include <TProfile.h>

#include <AliVEvent.h>
#include <AliVParticle.h>
#include <AliKFParticle.h>
#include <AliESDtrackCuts.h>
#include <AliESDEvent.h>
#include <AliMCEvent.h>
#include <AliAODEvent.h>
#include <AliAODTracklets.h>
#include <AliMultiplicity.h>
#include <AliStack.h>

#include "AliDielectronVarCuts.h"
#include "AliDielectronTrackCuts.h"
#include "AliDielectronVarManager.h"
#include "AliDielectronHelper.h"

//_____________________________________________________________________________
TVectorD* AliDielectronHelper::MakeLogBinning(Int_t nbinsX, Double_t xmin, Double_t xmax)
{
  //
  // Make logarithmic binning
  // the user has to delete the array afterwards!!!
  //
  
  //check limits
  if (xmin<1e-20 || xmax<1e-20){
    Error("AliDielectronHelper::MakeLogBinning","For Log binning xmin and xmax must be > 1e-20. Using linear binning instead!");
    return AliDielectronHelper::MakeLinBinning(nbinsX, xmin, xmax);
  }
  if (xmax<xmin){
    Double_t tmp=xmin;
    xmin=xmax;
    xmax=tmp;
  }
  TVectorD *binLim=new TVectorD(nbinsX+1);
  Double_t first=xmin;
  Double_t last=xmax;
  Double_t expMax=TMath::Log(last/first);
  for (Int_t i=0; i<nbinsX+1; ++i){
    (*binLim)[i]=first*TMath::Exp(expMax/nbinsX*(Double_t)i);
  }
  return binLim;
}

//_____________________________________________________________________________
TVectorD* AliDielectronHelper::MakeLinBinning(Int_t nbinsX, Double_t xmin, Double_t xmax)
{
  //
  // Make linear binning
  // the user has to delete the array afterwards!!!
  //
  if (xmax<xmin){
    Double_t tmp=xmin;
    xmin=xmax;
    xmax=tmp;
  }
  TVectorD *binLim=new TVectorD(nbinsX+1);
  Double_t first=xmin;
  Double_t last=xmax;
  Double_t binWidth=(last-first)/nbinsX;
  for (Int_t i=0; i<nbinsX+1; ++i){
    (*binLim)[i]=first+binWidth*(Double_t)i;
  }
  return binLim;
}

//_____________________________________________________________________________
TVectorD* AliDielectronHelper::MakeArbitraryBinning(const char* bins)
{
  //
  // Make arbitrary binning, bins separated by a ','
  //
  TString limits(bins);
  if (limits.IsNull()){
    Error("AliDielectronHelper::MakeArbitraryBinning","Bin Limit string is empty, cannot add the variable");
    return 0x0;
  }
  
  TObjArray *arr=limits.Tokenize(",");
  Int_t nLimits=arr->GetEntries();
  if (nLimits<2){
    Error("AliDielectronHelper::MakeArbitraryBinning","Need at leas 2 bin limits, cannot add the variable");
    delete arr;
    return 0x0;
  }
  
  TVectorD *binLimits=new TVectorD(nLimits);
  for (Int_t iLim=0; iLim<nLimits; ++iLim){
    (*binLimits)[iLim]=(static_cast<TObjString*>(arr->At(iLim)))->GetString().Atof();
  }
  
  delete arr;
  return binLimits;
}


//_____________________________________________________________________________
Int_t AliDielectronHelper::GetNch(const AliMCEvent *ev, Double_t etaRange){
  // determination of Nch
  if (!ev || ev->IsA()!=AliMCEvent::Class()) return -1;

  AliStack *stack = ((AliMCEvent*)ev)->Stack();

  if (!stack) return -1;

  Int_t nParticles = stack->GetNtrack();
  Int_t nCh = 0;

  // count..
  for (Int_t iMc = 0; iMc < nParticles; ++iMc) {
    if (!stack->IsPhysicalPrimary(iMc)) continue;

    TParticle* particle = stack->Particle(iMc);
    if (!particle) continue;
    if (particle->GetPDG()->Charge() == 0) continue;

    Float_t eta = particle->Eta();
    if (TMath::Abs(eta) < TMath::Abs(etaRange)) nCh++;
  }

  return nCh;
}


//_____________________________________________________________________________
Int_t AliDielectronHelper::GetNaccTrcklts(const AliVEvent *ev, Double_t etaRange){
  // Compute the collision multiplicity based on AOD or ESD tracklets
  // Code taken from: AliAnalysisTaskMuonCollisionMultiplicity::ComputeMultiplicity()

  if (!ev) return -1;

  Int_t nTracklets = 0;
  Int_t nAcc = 0;
  
  if (ev->IsA() == AliAODEvent::Class()) {
    AliAODTracklets *tracklets = ((AliAODEvent*)ev)->GetTracklets();
    nTracklets = tracklets->GetNumberOfTracklets();
    for (Int_t nn = 0; nn < nTracklets; nn++) {
      Double_t theta = tracklets->GetTheta(nn);
      Double_t eta = -TMath::Log(TMath::Tan(theta/2.0));
      if (TMath::Abs(eta) < etaRange) nAcc++;
    }
  } else if (ev->IsA() == AliESDEvent::Class()) {
    nTracklets = ((AliESDEvent*)ev)->GetMultiplicity()->GetNumberOfTracklets();
    for (Int_t nn = 0; nn < nTracklets; nn++) {
      Double_t eta = ((AliESDEvent*)ev)->GetMultiplicity()->GetEta(nn);
      if (TMath::Abs(eta) < etaRange) nAcc++;
    }
  } else return -1;

  return nAcc;
}



//________________________________________________________________
Double_t AliDielectronHelper::GetNaccTrckltsCorrected(const AliVEvent *event, Double_t uncorrectedNacc, Double_t vtxZ, Int_t type) {
  //
  // Correct the number of accepted tracklets based on the period and event vertex
  //

  Int_t runNo = event->GetRunNumber();

  Int_t period = -1;   // 0-LHC10b, 1-LHC10c, 2-LHC10d, 3-LHC10e
  Double_t refMult = 0.0;   // reference multiplicity
  
  if(runNo>114930 && runNo<117223) period = 0;
  if(runNo>119158 && runNo<120830) period = 1;
  if(runNo>122373 && runNo<126438) period = 2;
  if(runNo>127711 && runNo<130841) period = 3;
  if(period<0 || period>3) return uncorrectedNacc;

  if(type<0 || type>8) return uncorrectedNacc;
  if(type == 0) refMult = 5.0;         // SPD tracklets in |eta|<0.5 
  if(type == 1) refMult = 9.5;         // SPD tracklets in |eta|<1.0
  if(type == 2) refMult = 13.0;        // SPD tracklets in |eta|<1.6
  if(type == 3) refMult = 6.0;         // ITSTPC+ in |eta|<0.5
  if(type == 4) refMult = 12.0;        // ITSTPC+ in |eta|<1.0
  if(type == 5) refMult = 16.0;        // ITSTPC+ in |eta|<1.6
  if(type == 6) refMult = 6.0;         // ITSSA+ in |eta|<0.5
  if(type == 7) refMult = 12.0;        // ITSSA+ in |eta|<1.0
  if(type == 8) refMult = 15.0;        // ITSSA+ in |eta|<1.6

  if(TMath::Abs(vtxZ)>10.0) return uncorrectedNacc;

  TProfile* estimatorAvg = AliDielectronVarManager::GetEstimatorHistogram(period, type);
  if(!estimatorAvg) return uncorrectedNacc;

  Double_t localAvg = estimatorAvg->GetBinContent(estimatorAvg->FindBin(vtxZ));

  Double_t deltaM = uncorrectedNacc*(refMult/localAvg - 1);

  Double_t correctedNacc = uncorrectedNacc + (deltaM>0 ? 1 : -1) * gRandom->Poisson(TMath::Abs(deltaM));

  return correctedNacc;
}

//_____________________________________________________________________________
Int_t AliDielectronHelper::GetNacc(const AliVEvent *ev){
  // put a robust Nacc definition here

  if (!ev) return -1;
  
  AliDielectronVarCuts *varCuts = new AliDielectronVarCuts("VarCuts","VarCuts");
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.9,   0.9);
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,     70.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,  -0.5,   0.5);   //noKinks
    
  AliDielectronTrackCuts *trkCuts = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  trkCuts->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kAny);
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE);

  Int_t nRecoTracks = ev->GetNumberOfTracks();
  Int_t nAcc = 0;
  
  for (Int_t iTrack = 0; iTrack < nRecoTracks; iTrack++) {
    AliVParticle* candidate = ev->GetTrack(iTrack);
    if (!candidate) continue;
    AliVTrack *track        = static_cast<AliVTrack*>(candidate);
    if (!track) continue;
    if (varCuts->IsSelected(track) && trkCuts->IsSelected(track)) 
      nAcc++;
  }
  
  delete varCuts;
  delete trkCuts;

  return nAcc;
}

//_____________________________________________________________________________
void AliDielectronHelper::RotateKFParticle(AliKFParticle * kfParticle,Double_t angle, const AliVEvent * const ev){
  // Before rotate needs to be moved to position 0,0,0, ; move back after rotation
  if (!kfParticle) return;
  Double_t dx = 0.;
  Double_t dy = 0.;
  Double_t dz = 0.;

  if (ev){
    dx = ev->GetPrimaryVertex()->GetX()-0.;
    dy = ev->GetPrimaryVertex()->GetY()-0.;
    dz = ev->GetPrimaryVertex()->GetZ()-0.;
  }
  
  kfParticle->X() = kfParticle->GetX() - dx;
  kfParticle->Y() = kfParticle->GetY() - dy;
  kfParticle->Z() = kfParticle->GetZ() - dz;
  
  
  // Rotate the kf particle
  Double_t c = cos(angle);
  Double_t s = sin(angle);
  
  Double_t mA[8][ 8];
  for( Int_t i=0; i<8; i++ ){
    for( Int_t j=0; j<8; j++){
      mA[i][j] = 0;
    }
  }
  for( int i=0; i<8; i++ ){
    mA[i][i] = 1;
  }
  mA[0][0] =  c;  mA[0][1] = s;
  mA[1][0] = -s;  mA[1][1] = c;
  mA[3][3] =  c;  mA[3][4] = s;
  mA[4][3] = -s;  mA[4][4] = c;
  
  Double_t mAC[8][8];
  Double_t mAp[8];
  
  for( Int_t i=0; i<8; i++ ){
    mAp[i] = 0;
    for( Int_t k=0; k<8; k++){
      mAp[i]+=mA[i][k] * kfParticle->GetParameter(k);
    }
  }
  
  for( Int_t i=0; i<8; i++){
    kfParticle->Parameter(i) = mAp[i];
  }
  
  for( Int_t i=0; i<8; i++ ){
    for( Int_t j=0; j<8; j++ ){
      mAC[i][j] = 0;
      for( Int_t k=0; k<8; k++ ){
        mAC[i][j]+= mA[i][k] * kfParticle->GetCovariance(k,j);
      }
    }
  }
  
  for( Int_t i=0; i<8; i++ ){
    for( Int_t j=0; j<=i; j++ ){
      Double_t xx = 0;
      for( Int_t k=0; k<8; k++){
        xx+= mAC[i][k]*mA[j][k];
      }
      kfParticle->Covariance(i,j) = xx;
    }
  }
  
  kfParticle->X() = kfParticle->GetX() + dx;
  kfParticle->Y() = kfParticle->GetY() + dy;
  kfParticle->Z() = kfParticle->GetZ() + dz;
  
}

//_____________________________________________________________________________
Int_t AliDielectronHelper::GetNMothers(const AliMCEvent *ev, Double_t etaRange, Int_t pdgMother, Int_t pdgDaughter, Int_t prim){
  // counting number of mother particles generated in given eta range and 2 particle decay
  if (!ev || ev->IsA()!=AliMCEvent::Class()) return -1;
  
  AliStack *stack = ((AliMCEvent*)ev)->Stack();
  
  if (!stack) return -1;
  
  Int_t nParticles = stack->GetNtrack();
  Int_t nMothers   = 0;
  
  // count..
  for (Int_t iMc = 0; iMc < nParticles; ++iMc) {
    
    TParticle* particle = stack->Particle(iMc);
    if (!particle) continue;
    if (particle->GetPdgCode() != pdgMother)               continue;
    if (TMath::Abs(particle->Eta()) > TMath::Abs(etaRange)) continue;

    if (particle->GetNDaughters() != 2)                 continue;
    // 1st daugther
    if (particle->GetFirstDaughter()>=nParticles ||
	particle->GetFirstDaughter()<0             ) continue;
    
    TParticle* dau1 = stack->Particle(particle->GetFirstDaughter());
    if (TMath::Abs(dau1->GetPdgCode()) != pdgDaughter)     continue;
    if (TMath::Abs(dau1->Eta()) > TMath::Abs(etaRange)) continue;
    
    // 2nd daughter
    if (particle->GetLastDaughter()>=nParticles ||
	particle->GetLastDaughter()<0             ) continue;

    TParticle* dau2 = stack->Particle(particle->GetLastDaughter());
    if (TMath::Abs(dau2->GetPdgCode()) != pdgDaughter)     continue;
    if (TMath::Abs(dau2->Eta()) > TMath::Abs(etaRange)) continue;
    
    // primary
    if (prim != -1) {
      if(particle->IsPrimary() != prim) continue;
    }
    nMothers++;
  }
  return nMothers;
}


