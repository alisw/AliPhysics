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
//   Julian Book <Julian.Book@cern.ch>

#include <iostream>
#include <vector>

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
#include <AliESDtrack.h>
#include <AliESDEvent.h>
#include <AliMCEvent.h>
#include <AliAODEvent.h>
#include <AliAODTracklets.h>
#include <AliMultiplicity.h>

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
void AliDielectronHelper::GetMaxPtAndPhi(const AliVEvent *ev, Double_t &ptMax, Double_t &phiOfptMax){
  //
  // find the highest pt and its phi in the event
  //
  for(Int_t itrk=0; itrk<ev->GetNumberOfTracks(); itrk++) {
    AliVParticle *part= ev->GetTrack(itrk);
    if(part && part->Pt() > ptMax) {
      ptMax      = part->Pt();
      phiOfptMax = part->Phi();
    }
  }


}

//_____________________________________________________________________________
Int_t AliDielectronHelper::GetNch(const AliMCEvent *ev, Double_t etaRange, Bool_t excludeJpsiDaughters){
  // determination of Nch
  // Obsolete function see CountMCtracks
  if (!ev || ev->IsA()!=AliMCEvent::Class()) return -1;
  Int_t nParticles = ev->GetNumberOfTracks();
  Int_t nCh = 0;
  if(!ev->Particle(1) && !ev->GetTrack(1))  return -1;

  Bool_t isAOD = kFALSE;
  if(ev->GetTrack(1)->IsA() == AliAODMCParticle::Class()) isAOD = kTRUE;
  if(isAOD){
    for (Int_t iParticle = 1; iParticle < nParticles; iParticle++) {
      AliAODMCParticle *mcPar = (AliAODMCParticle*) ev->GetTrack(iParticle);
      if(!mcPar) continue;
      if(!mcPar->IsPhysicalPrimary()) continue;
      if(mcPar->Charge() == 0) continue;

      if(excludeJpsiDaughters){
        Int_t iMother = mcPar->GetMother();
        if(iMother > 0){
          AliAODMCParticle *mcParMother = (AliAODMCParticle*) ev->GetTrack(iMother);
          if(mcParMother)
          if(TMath::Abs(mcParMother->GetPdgCode()) == 443) continue;
        }
      }
      Float_t eta = mcPar->Eta();
      if (TMath::Abs(eta) < TMath::Abs(etaRange)) nCh++;
    }
  }
  // count.. for ESD input
  else{
    for (Int_t iMc = 0; iMc < nParticles; ++iMc) {
      if (!ev->IsPhysicalPrimary(iMc)) continue;

      TParticle* particle = ev->Particle(iMc);
      if (!particle) continue;
      if (particle->GetPDG()->Charge() == 0) continue;

      if(excludeJpsiDaughters){
        Int_t iMother = particle->GetMother(0);
        if( ! (iMother < 0) ) {
          TParticle* mother = ev->Particle(iMother);
          if( mother && TMath::Abs( mother->GetPdgCode() ) == 443  ) continue;
        }
      }
      Float_t eta = particle->Eta();
      if (TMath::Abs(eta) < TMath::Abs(etaRange)) nCh++;
    }
  }
  return nCh;
}

//_____________________________________________________________________________
Double_t AliDielectronHelper::GetTransverseSpherocity(const AliVEvent *ev){

  //GetNtracks
  const Int_t Ntracks = ev->GetNumberOfTracks();

  //fPx fPy
  AliVParticle *part=0x0;
  Double_t Px=-99;
  Double_t Py=-99;

  if(Ntracks < 10)
    return -1;
  Double_t stepSize=0.1;
  Double_t RetTransverseSpherocity = 1000;
  Double_t sumpt = 0;
  Int_t steplimit = 360/stepSize;
  for(Int_t i = 0; i < steplimit; ++i) {
    //Divide the whole azimuth into segments and do the projection on these segments (below)
    Double_t phiparam = ((TMath::Pi()) * i * stepSize) / 180;
    Double_t nx = TMath::Cos(phiparam); // x component of a unitary vector n
    Double_t ny = TMath::Sin(phiparam); // y component of a unitary vector n

    Double_t num = 0;
    for(Int_t j = 0; j < Ntracks; j++) {
      part= ev->GetTrack(j);
      Px = part->Px(); 
      Py = part->Py();
      num += TMath::Abs(ny*Px - nx*Py);

      if(i==0)
        sumpt += TMath::Sqrt(Px*Px + Py*Py);
    }

    Double_t pFull = TMath::Power((num/sumpt), 2); //Projection of sp. on the segment
    if(pFull < RetTransverseSpherocity)  //Select the lowest projection
      RetTransverseSpherocity = pFull;
  };

  RetTransverseSpherocity *= TMath::Pi()*TMath::Pi()/4.0;

  return RetTransverseSpherocity;
}


//_____________________________________________________________________________
Double_t AliDielectronHelper::GetTransverseSpherocityTracks(const AliVEvent *ev){

  //GetNtracks
  const Int_t Ntracks = ev->GetNumberOfTracks();

  //fPx fPy
  AliVParticle *part=0x0;
  Double_t Px=-99;
  Double_t Py=-99;


  if(Ntracks < 10)
    return -1;

  Double_t RetTransverseSpherocity = 1000;
  Double_t sumpt = 0;
  //const Double_t pt = 1;
  for(Int_t i = 0; i < Ntracks; i++) {

    part = ev->GetTrack(i);
    Px = part->Px(); 
    Py = part->Py();
    Double_t pt = part->Pt();
    Double_t nx = Px / pt; // x component of a unitary vector n
    Double_t ny = Py / pt; // y component of a unitary vector n

    Double_t num = 0;
    for(Int_t j = 0; j < Ntracks; j++) {
      part = ev->GetTrack(j);
      Px = part->Px(); 
      Py = part->Py();
      num += TMath::Abs(ny*Px - nx*Py);

      if(i==0)
        sumpt += part->Pt();
    }

    Double_t pFull = TMath::Power((num/sumpt), 2); //Projection of sp. on the segment
    if(pFull < RetTransverseSpherocity)  { //Select the lowest projection
      RetTransverseSpherocity = pFull;
    };
  };

  RetTransverseSpherocity *= TMath::Pi()*TMath::Pi()/4.0;

  return RetTransverseSpherocity;
 
}

//_____________________________________________________________________________
Bool_t AliDielectronHelper::CheckESDtrack(AliESDtrack *track){

  AliESDtrackCuts* TrackCutsTPCRefit = new AliESDtrackCuts("AliESDtrackCutsTPCRefit","AliESDtrackCutsTPCRefit");
  TrackCutsTPCRefit = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts(); //If not running, set to kFALSE;
  TrackCutsTPCRefit->SetRequireTPCRefit(kTRUE);
  TrackCutsTPCRefit->SetRequireITSRefit(kTRUE);
  TrackCutsTPCRefit->SetEtaRange(-0.8,0.8);

  if(TrackCutsTPCRefit->AcceptTrack(track)) return kTRUE;
  else return kFALSE;

}


//_____________________________________________________________________________
std::vector<AliESDtrack*> AliDielectronHelper::GetESDtracks(const AliESDEvent *ev){

  //GetNtracks
  const Int_t Ntracks = ev->GetNumberOfTracks();

  AliESDtrackCuts* TrackCutsTPCRefit = new AliESDtrackCuts("AliESDtrackCutsTPCRefit","AliESDtrackCutsTPCRefit");
  TrackCutsTPCRefit = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts(); //If not running, set to kFALSE;
  TrackCutsTPCRefit->SetRequireTPCRefit(kTRUE);
  TrackCutsTPCRefit->SetRequireITSRefit(kTRUE);
  TrackCutsTPCRefit->SetEtaRange(-0.8,0.8);

  AliESDtrack *track = 0x0;
  std::vector<AliESDtrack *> tracks;
  for(Int_t i=0; i<Ntracks; i++){

    track = ev->GetTrack(i); 
    if(TrackCutsTPCRefit->AcceptTrack(track)) tracks.push_back(track);
  }

  delete TrackCutsTPCRefit;
  TrackCutsTPCRefit=0x0;
  return tracks;

}

//_____________________________________________________________________________
Double_t AliDielectronHelper::GetTransverseSpherocityESD(const AliESDEvent *ev){ 

  std::vector<AliESDtrack *> tracks;
  AliESDtrack *track = 0x0;
  tracks = GetESDtracks(ev);

  const int Ntracks = tracks.size();
  if(Ntracks<4)
    return -1.;

  //fPx fPy
  Double_t Px=-99;
  Double_t Py=-99;

  Double_t stepSize=0.1;
  Double_t RetTransverseSpherocity = 1000;
  Double_t sumpt = 0;
  Int_t steplimit = 360/stepSize;
  for(Int_t i = 0; i < steplimit; ++i) {
    //Divide the whole azimuth into segments and do the projection on these segments (below)
    Double_t phiparam = ((TMath::Pi()) * i * stepSize) / 180;
    Double_t nx = TMath::Cos(phiparam); // x component of a unitary vector n
    Double_t ny = TMath::Sin(phiparam); // y component of a unitary vector n

    Double_t num = 0;
    for(Int_t j = 0; j < Ntracks; j++) {
      
      track = tracks.at(j);
      Px = track->Px(); 
      Py = track->Py();
      num += TMath::Abs(ny*Px - nx*Py);

      if(i==0)
        sumpt += TMath::Sqrt(Px*Px + Py*Py);
    }

    Double_t pFull = TMath::Power((num/sumpt), 2); //Projection of sp. on the segment
    if(pFull < RetTransverseSpherocity)  //Select the lowest projection
      RetTransverseSpherocity = pFull;
  };

  RetTransverseSpherocity *= TMath::Pi()*TMath::Pi()/4.0;

  return RetTransverseSpherocity;
}


//_____________________________________________________________________________
Double_t AliDielectronHelper::GetTransverseSpherocityESDtracks(const AliESDEvent *ev){

  std::vector<AliESDtrack *> tracks;
  AliESDtrack *track = 0x0;
  tracks = GetESDtracks(ev);

  const int Ntracks = tracks.size();
  if(Ntracks<10)
    return -1.;

  //fPx fPy
  Double_t Px=-99.;
  Double_t Py=-99.;

  Double_t RetTransverseSpherocity = 1000.;
  Double_t sumpt = 0;
  Int_t firstTrack = 0;
  //const Double_t pt = 1;
  for(Int_t i = 0; i < Ntracks; i++) {

    track = tracks.at(i);
    Px = track->Px();
    Py = track->Py();
    Double_t pt = track->Pt();
    Double_t nx = Px / pt; // x component of a unitary vector n
    Double_t ny = Py / pt; // y component of a unitary vector n

    Double_t num = 0;
    for(Int_t j = 0; j < Ntracks; j++) {
      track = tracks.at(j);
      Px = track->Px();
      Py = track->Py();
      //Px = TMath::Cos(part->Phi());
      //Py = TMath::Sin(part->Phi());
      num += TMath::Abs(ny*Px - nx*Py);

      if(i==0)
        sumpt += TMath::Sqrt(Px*Px + Py*Py);
    }

    Double_t pFull = TMath::Power((num/sumpt), 2); //Projection of sp. on the segment
    if(pFull < RetTransverseSpherocity)  { //Select the lowest projection
      RetTransverseSpherocity = pFull;
    };
  };

  RetTransverseSpherocity *= TMath::Pi()*TMath::Pi()/4.0;
  if(RetTransverseSpherocity < 0.) return -1.1;
  if(RetTransverseSpherocity > 1.) return -0.9;
  return RetTransverseSpherocity;

}


//_____________________________________________________________________________
Double_t AliDielectronHelper::GetTransverseSpherocityESDwoPtWeight(const AliESDEvent *ev){

  std::vector<AliESDtrack *> tracks;
  AliESDtrack *track = 0x0;
  tracks = GetESDtracks(ev);

  const int Ntracks = tracks.size();
  if(Ntracks<4)
    return -1.;

  //fPx fPy
  Double_t Px=-99;
  Double_t Py=-99;

  Double_t stepSize=0.1;
  Double_t RetTransverseSpherocity = 1000;
  Double_t sumpt = 0;
  Int_t steplimit = 360/stepSize;
  for(Int_t i = 0; i < steplimit; ++i) {
    //Divide the whole azimuth into segments and do the projection on these segments (below)
    Double_t phiparam = ((TMath::Pi()) * i * stepSize) / 180;
    Double_t nx = TMath::Cos(phiparam); // x component of a unitary vector n
    Double_t ny = TMath::Sin(phiparam); // y component of a unitary vector n

    Double_t num = 0;
    for(Int_t j = 0; j < Ntracks; j++) {
      
      track = tracks.at(j);
      Px = track->Px() / track->Pt(); 
      Py = track->Py() / track->Pt();
      num += TMath::Abs(ny*Px - nx*Py);

      if(i==0)
        sumpt += TMath::Sqrt(Px*Px + Py*Py);
    }

    Double_t pFull = TMath::Power((num/sumpt), 2); //Projection of sp. on the segment
    if(pFull < RetTransverseSpherocity)  //Select the lowest projection
      RetTransverseSpherocity = pFull;
  };

  RetTransverseSpherocity *= TMath::Pi()*TMath::Pi()/4.0;

  return RetTransverseSpherocity;
}



//_____________________________________________________________________________
Double_t AliDielectronHelper::GetTransverseSpherocityESDtracksWoPtWeight(const AliESDEvent *ev){

  std::vector<AliESDtrack *> tracks;
  AliESDtrack *track = 0x0;
  tracks = GetESDtracks(ev);

  const int Ntracks = tracks.size();
  if(Ntracks<10)
    return -1.;

  //fPx fPy
  Double_t Px=-99.;
  Double_t Py=-99.;

  Double_t RetTransverseSpherocity = 1000.;
  Double_t sumpt = 0;
  Int_t firstTrack = 0;
  //const Double_t pt = 1;
  for(Int_t i = 0; i < Ntracks; i++) {

    track = tracks.at(i);
    Px = track->Px();
    Py = track->Py();
    Double_t pt = track->Pt();
    Double_t nx = Px / pt; // x component of a unitary vector n
    Double_t ny = Py / pt; // y component of a unitary vector n

    Double_t num = 0;
    for(Int_t j = 0; j < Ntracks; j++) {
      track = tracks.at(j);
      Px = track->Px() / track->Pt();
      Py = track->Py() / track->Pt();
      //Px = TMath::Cos(part->Phi());
      //Py = TMath::Sin(part->Phi());
      num += TMath::Abs(ny*Px - nx*Py);

      if(i==0)
        sumpt += TMath::Sqrt( Px*Px + Py*Py);
    }

    Double_t pFull = TMath::Power((num/sumpt), 2); //Projection of sp. on the segment
    if(pFull < RetTransverseSpherocity)  { //Select the lowest projection
      RetTransverseSpherocity = pFull;
    };
  };


  RetTransverseSpherocity *= TMath::Pi()*TMath::Pi()/4.0;
  if(RetTransverseSpherocity < 0.) return -1.1;
  if(RetTransverseSpherocity > 1.) return -0.9;

  return RetTransverseSpherocity;

}


//_____________________________________________________________________________
Int_t AliDielectronHelper::GetNaccTrcklts(const AliVEvent *ev, Double_t etaRange)
{
  // Compute the collision multiplicity based on AOD or ESD tracklets
  // Code taken from: AliAnalysisTaskMuonCollisionMultiplicity::ComputeMultiplicity()
  if (!ev) return -1;

  Int_t nTracklets = 0;
  Int_t nAcc = 0;

  if (ev->IsA() == AliAODEvent::Class()) {
    AliAODTracklets *tracklets = ((AliAODEvent*)ev)->GetTracklets();
    if (!tracklets) return -1;
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
Double_t AliDielectronHelper::GetNaccTrckltsCorrected(const AliVEvent *event, Double_t uncorrectedNacc, Double_t vtxZ, Int_t type)
{
  //
  // Correct the number of accepted tracklets based on the period and event vertex
  //

  Int_t runNo = event->GetRunNumber();

  Int_t period = -1;   // 0-LHC10b, 1-LHC10c, 2-LHC10d, 3-LHC10e
  Double_t refMult = 0.0;   // reference multiplicity

  //  if(runNo>114930 && runNo<117223) period = 0; Org Fred Kramer Analysis
  if(runNo>114736 && runNo<117224) period = 0; // LHC10b PASS4 2015
  //  if(runNo>119158 && runNo<120830) period = 1; Org Fred Kramer Analysis
  if(runNo>118358 && runNo<121041) period = 1; // LHC10c PASS4 2015
  //  if(runNo>122373 && runNo<126438) period = 2; Org Fred Kramer Analysis
  if(runNo>121693 && runNo<126439) period = 2; // LHC10d PASS4 2015
  //  if(runNo>127711 && runNo<130841) period = 3; Org Fred Kramer Analysis
  if(runNo>127102 && runNo<130851) period = 3; // LHC10e PASS4 2015
  //pp 13 TeV CJ
  if(runNo>254124 && runNo<264035) period = 4; // LHC16l,16k,16i,16j,16o pass1 2016 CJ analysis

  if(period<0 || period>4) return uncorrectedNacc;

  if(type<0 || type>8) return uncorrectedNacc;
  if(type == 0) refMult = 5.0;         // SPD tracklets in |eta|<0.5
  if(type == 1 && period!=4) refMult = 9.5;         // SPD tracklets in |eta|<1.0
  if(type == 1 && period==4) refMult = 12.0;        // SPD tracklets in |eta|<1.0 refMulti for period 4 (LHC16l) //CJ analysis
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

  AliDielectronVarCuts varCuts;
  varCuts.AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCuts.AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCuts.AddCut(AliDielectronVarManager::kEta,         -0.9,   0.9);
  varCuts.AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  varCuts.AddCut(AliDielectronVarManager::kNclsTPC,     70.0, 160.0);
  varCuts.AddCut(AliDielectronVarManager::kKinkIndex0,  -0.5,   0.5);   //noKinks

  AliDielectronTrackCuts trkCuts;
  trkCuts.SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kAny);
  trkCuts.SetRequireITSRefit(kTRUE);
  trkCuts.SetRequireTPCRefit(kTRUE);

  Int_t nRecoTracks = ev->GetNumberOfTracks();
  Int_t nAcc = 0;

  // @TODO: IsSelected() sets the VarManager fgFillMap to the 'trkCuts' fillmap.
  // Since 'trkCuts' is a stack object and will be deleted, the VarManager fgFillMap goes into an undefined state.
  for (Int_t iTrack = 0; iTrack < nRecoTracks; iTrack++) {
    AliVTrack *track        = static_cast<AliVTrack*>(ev->GetTrack(iTrack));
    if (!track) continue;
    if (!trkCuts.IsSelected(track)) continue;
    if (!varCuts.IsSelected(track)) continue;
    nAcc++;
  }

  return nAcc;
}

//_____________________________________________________________________________
Double_t AliDielectronHelper::GetITSTPCMatchEff(const AliVEvent *ev, Double_t *efficiencies, Bool_t bEventPlane, Bool_t useV0Cep)
{
  // recalulate the its-tpc matching efficiecy
  if (!ev) return -1;


  // AliDielectronVarCuts varCutsTPC;
  // varCutsTPC.AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  // varCutsTPC.AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  // varCutsTPC.AddCut(AliDielectronVarManager::kEta,         -0.9,   0.9);
  // varCutsTPC.AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  // varCutsTPC.AddCut(AliDielectronVarManager::kNclsTPC,     50.0, 160.0);
  AliDielectronTrackCuts trkCutsTPC;
  trkCutsTPC.SetRequireTPCRefit(kTRUE);

  // AliDielectronVarCuts varCutsITS;
  // varCutsITS.AddCut(AliDielectronVarManager::kEta,         -0.9,   0.9);
  AliDielectronTrackCuts trkCutsITS;
  trkCutsITS.SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kAny);
  trkCutsITS.SetRequireITSRefit(kTRUE);

  // Eventplane angle Cuts

  // AliDielectronVarCuts varCutsEPinP;
  // AliDielectronVarCuts varCutsEPoutP;
  // if(bEventPlane){
  //   varCutsEPinP.AddCut(AliDielectronVarManager::kQnDeltaPhiTrackTPCrpH2, -2.35, -0.79, kTRUE);
  //   varCutsEPinP.AddCut(AliDielectronVarManager::kQnDeltaPhiTrackTPCrpH2, 0.79, 2.35, kTRUE);
  //   varCutsEPoutP.AddCut(AliDielectronVarManager::kQnDeltaPhiTrackTPCrpH2, -0.79, 0.79, kTRUE);
  //   varCutsEPoutP.AddCut(AliDielectronVarManager::kQnDeltaPhiTrackTPCrpH2, -3.15, -2.35, kTRUE);
  //   varCutsEPoutP.AddCut(AliDielectronVarManager::kQnDeltaPhiTrackTPCrpH2, 2.35, 3.15, kTRUE);
  // }

  // @TODO: IsSelected() sets the VarManager fgFillMap to the 'trkCutsITS' fillmap.
  // Since 'trkCutsITS' is a stack object and will be deleted, the VarManager fgFillMap goes into an undefined state.

  Int_t nRecoTracks = ev->GetNumberOfTracks();
  Double_t nTPC = 0, nITS = 0;
  Double_t nTPCinP = 0, nITSinP = 0;
  Double_t nTPCoutP = 0, nITSoutP = 0;

  Float_t impactParXY = -99.;
  Float_t impactParZ = -99.;
  Float_t eta = -99.;
  Float_t nClsTPC = -99.;
  Float_t tpcChi2Cl = -99.;
  Float_t eventplane = -99.;

  TString epDetector = "TPC";
  if(useV0Cep)  epDetector = "V0C";



  if(bEventPlane){
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    if(man)
    if( AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask =
        dynamic_cast<AliAnalysisTaskFlowVectorCorrections*> (man->GetTask("FlowQnVectorCorrections")) )
      if(flowQnVectorTask != NULL){
        AliQnCorrectionsManager *flowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
        TList *qnlist = flowQnVectorMgr->GetQnVectorList();
        if(qnlist != NULL){
          const AliQnCorrectionsQnVector *qVecQnFrameworkTPC = AliDielectronQnEPcorrection::GetQnVectorFromList( qnlist, epDetector.Data(), "latest", "latest" );
          if(qVecQnFrameworkTPC != NULL){
            TVector2 qVectorTPC( qVecQnFrameworkTPC->Qx(2), qVecQnFrameworkTPC->Qy(2) );
            eventplane = TVector2::Phi_mpi_pi(qVectorTPC.Phi())/2;
          }
        }
      }
  }

  for (Int_t iTrack = 0; iTrack < nRecoTracks; iTrack++) {
    AliVTrack *track        = static_cast<AliVTrack*>(ev->GetTrack(iTrack));
    if (!track) continue;

    // ITS cuts begin
    if(!trkCutsITS.IsSelected(track)) continue;
    //
    if(track->IsA() == AliESDtrack::Class())track->GetImpactParameters(impactParXY, impactParZ);
    else{
      AliAODTrack *trackAOD = (AliAODTrack*) track;
      Double_t xyz[2] = {-100.,-100.};
      Double_t dcaRes[3] = {-100.,-100.,-100.};
      AliDielectronVarManager::GetDCA(trackAOD, xyz, dcaRes);
      impactParXY = xyz[0];
      impactParZ = xyz[1];
    }
    if(TMath::Abs(impactParXY) > 1.0) continue;
    if(TMath::Abs(impactParZ) > 3.0) continue;
    eta = track->Eta();
    if(TMath::Abs(eta) > 0.9) continue;
    //
    // if(!varCutsITS.IsSelected(track)) continue;
    //
    // ITS cuts end
    nITS+=1.;

    // Eventplane cuts for ITS
    if(bEventPlane){
      if(AliDielectronHelper::IsInPlane(track->Phi(),eventplane))  nITSinP++;
      if(AliDielectronHelper::IsOutOfPlane(track->Phi(),eventplane))  nITSoutP++;
    }


    // TPC cuts begin
    if(!trkCutsTPC.IsSelected(track)) continue;
    track->GetImpactParameters(impactParXY, impactParZ);
    if(TMath::Abs(impactParXY) > 1.0) continue;
    if(TMath::Abs(impactParZ) > 3.0) continue;
    eta = track->Eta();
    if(TMath::Abs(eta) > 0.9) continue;
    nClsTPC = track->GetTPCNcls();
    if(nClsTPC < 50. || nClsTPC > 160.) continue;
    tpcChi2Cl = nClsTPC > 0 ? track->GetTPCchi2() / nClsTPC : -1.;
    if(tpcChi2Cl < 0. || tpcChi2Cl > 4.) continue;
    // TPC cuts end


    // if(!varCutsTPC.IsSelected(track)) continue;
    nTPC+=1.;

    if(bEventPlane){
      if(AliDielectronHelper::IsInPlane(track->Phi(),eventplane))  nTPCinP++;
      if(AliDielectronHelper::IsOutOfPlane(track->Phi(),eventplane))  nTPCoutP++;
    }
  }
  if( bEventPlane && efficiencies){
    efficiencies[0] = nITSinP > 0. ? nTPCinP/nITSinP : -1.;
    efficiencies[1] = nITSoutP > 0. ? nTPCoutP/nITSoutP : -1.;
  }
  //  printf(" tracks TPC %.3e ITS %.3e = %.5f \n",nTPC,nITS,(nITS>0. ? nTPC/nITS : -1));


  return (nITS>0. ? nTPC/nITS : -1);

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
Int_t AliDielectronHelper::GetNMothers(const AliMCEvent *ev, Double_t etaRange, Int_t pdgMother, Int_t pdgDaughter, Int_t prim)
{
  // Obsolete function see CountMCtracks
  // counting number of mother particles generated in given eta range and 2 particle decay
  if (!ev || ev->IsA()!=AliMCEvent::Class()) return -1;

  if(!ev->Particle(1) && !ev->GetTrack(1))  return -1;

  Bool_t isAOD = kFALSE;
  if(ev->GetTrack(1)->IsA() == AliAODMCParticle::Class()) isAOD = kTRUE;

  Int_t nParticles = ev->GetNumberOfTracks();
  Int_t nMothers   = 0;

  // count.. for ESDs
  if(isAOD)
  {
    for (Int_t iMc = 1; iMc < nParticles; ++iMc) {
      AliAODMCParticle *mcPar = (AliAODMCParticle*) ev->GetTrack(iMc);
      if(!mcPar) continue;
      if(mcPar->GetPdgCode() != pdgMother) continue;
      if(TMath::Abs(mcPar->Eta()) > TMath::Abs(etaRange)) continue;

      // Daughters
      if(mcPar->GetNDaughters() != 2) continue;
      if(mcPar->GetDaughterFirst() >= nParticles || mcPar->GetDaughterFirst() < 0) continue;
      AliAODMCParticle *mcParDau1 = (AliAODMCParticle*) ev->GetTrack(mcPar->GetDaughterFirst());
      if(!mcParDau1) continue;
      if(TMath::Abs(mcParDau1->GetPdgCode()) != pdgDaughter) continue;
      if(TMath::Abs(mcParDau1->Eta()) > TMath::Abs(etaRange)) continue;

      if(mcPar->GetDaughterLast() >= nParticles || mcPar->GetDaughterLast() < 0) continue;
      AliAODMCParticle *mcParDau2 = (AliAODMCParticle*) ev->GetTrack(mcPar->GetDaughterLast());
      if(!mcParDau2) continue;
      if(TMath::Abs(mcParDau2->GetPdgCode()) != pdgDaughter) continue;
      if(TMath::Abs(mcParDau2->Eta()) > TMath::Abs(etaRange)) continue;
      if(prim != -1 && mcPar->IsPhysicalPrimary() != prim) continue;
      nMothers++;
    }
  }
  else{
    for (Int_t iMc = 0; iMc < nParticles; ++iMc) {

      TParticle* particle = ev->Particle(iMc);
      if (!particle) continue;
      if (particle->GetPdgCode() != pdgMother)               continue;
      if (TMath::Abs(particle->Eta()) > TMath::Abs(etaRange)) continue;

      if (particle->GetNDaughters() != 2)                 continue;
      // 1st daugther
      if (particle->GetFirstDaughter()>=nParticles ||
  	particle->GetFirstDaughter()<0             ) continue;

      TParticle* dau1 = ev->Particle(particle->GetFirstDaughter());
      if (TMath::Abs(dau1->GetPdgCode()) != pdgDaughter)     continue;
      if (TMath::Abs(dau1->Eta()) > TMath::Abs(etaRange)) continue;

      // 2nd daughter
      if (particle->GetLastDaughter()>=nParticles ||
  	particle->GetLastDaughter()<0             ) continue;

      TParticle* dau2 = ev->Particle(particle->GetLastDaughter());
      if (TMath::Abs(dau2->GetPdgCode()) != pdgDaughter)     continue;
      if (TMath::Abs(dau2->Eta()) > TMath::Abs(etaRange)) continue;

      // primary
      if (prim != -1) {
        if(particle->IsPrimary() != prim) continue;
      }
      nMothers++;
    }
  }
  return nMothers;
}

//______________________________________________
Bool_t AliDielectronHelper::IsInPlane(Float_t trackPhi, Float_t eventPhi){
  Float_t deltaPhi = TVector2::Phi_mpi_pi(trackPhi - eventPhi);
  if(TMath::Abs(deltaPhi) < TMath::Pi()/4 || (TMath::Abs(deltaPhi) < TMath::Pi() && TMath::Abs(deltaPhi) > TMath::Pi()*3/4) )  return kTRUE;
  else return kFALSE;
}

Bool_t AliDielectronHelper::IsOutOfPlane(Float_t trackPhi, Float_t eventPhi){
  Float_t deltaPhi = TVector2::Phi_mpi_pi(trackPhi - eventPhi);
  if(TMath::Abs(deltaPhi) > TMath::Pi()/4 && TMath::Abs(deltaPhi) < TMath::Pi()*3/4)  return kTRUE;
  else return kFALSE;
}

//______________________________________________
void AliDielectronHelper::CountMCtracks(const AliMCEvent *ev, Double_t numbers[11], Int_t pdgMother, Int_t pdgDaughter)
{
  /* Count the number of charged tracks, jpsi mothers and so on in one loop (more efficient then running GetNch and GetNMothers several times from the VarManager) */
  // 0 - N charged tracks, eta 1.6, incl jpsi daughters
  // 1 - N charged tracks, eta 1.6, excl jpsi daughters
  // 2 - N charged tracks, eta 0.5, incl jpsi daughters
  // 3 - N charged tracks, eta 0.5, excl jpsi daughters
  // 4 - N charged tracks, eta 1.0, incl jpsi daughters
  // 5 - N charged tracks, eta 1.0, excl jpsi daughters
  // 6 - N charged tracks, eta 0.9, incl jpsi daughters
  // 7 - N charged tracks, eta 0.9, excl jpsi daughters

  // 8 - N Jpsis, eta 0.9
  // 9 - N Prompt Jpsis, eta 0.9
  // 10 - N Non-Prompt Jpsis, eta 0.9

  for (Int_t i = 0; i < 11; i++) {
    numbers[i] = 0;
  }

  if (!ev || ev->IsA()!=AliMCEvent::Class()) return;

  if(!ev->Particle(1) && !ev->GetTrack(1))  return;

  Int_t nParticles = ev->GetNumberOfTracks();
  for (Int_t iMc = 1; iMc < nParticles; ++iMc) {
    AliVParticle *mcPar = (AliVParticle*) ev->GetTrack(iMc);
    // AliAODMCParticle *mcPar = (AliAODMCParticle*) ev->GetTrack(iMc);
    if(!mcPar) continue;
    Float_t eta = mcPar->Eta();
    if(mcPar->IsPhysicalPrimary() && mcPar->Charge() != 0){
      // count nCharged
      if(TMath::Abs(eta) < TMath::Abs(1.6))  numbers[0]++;
      if(TMath::Abs(eta) < TMath::Abs(0.5))  numbers[2]++;
      if(TMath::Abs(eta) < TMath::Abs(1.0))  numbers[4]++;
      if(TMath::Abs(eta) < TMath::Abs(0.9))  numbers[6]++;
      // count nCharged without jpsis daughters
      Int_t iMother = mcPar->GetMother();
      if(iMother > 0){
        AliVParticle *mcParMother = (AliVParticle*) ev->GetTrack(iMother);
        // AliAODMCParticle *mcParMother = (AliAODMCParticle*) ev->GetTrack(iMother);
        if(mcParMother)
        if(TMath::Abs(mcParMother->PdgCode()) != pdgMother){
        // if(TMath::Abs(mcParMother->GetPdgCode()) != pdgMother){
          if(TMath::Abs(eta) < TMath::Abs(1.6))  numbers[1]++;
          if(TMath::Abs(eta) < TMath::Abs(0.5))  numbers[3]++;
          if(TMath::Abs(eta) < TMath::Abs(1.0))  numbers[5]++;
          if(TMath::Abs(eta) < TMath::Abs(0.9))  numbers[7]++;
        }
      }
    } // End of counting nCharged
    // njpsis counting
    if(mcPar->PdgCode() != pdgMother) continue;
    // if(mcPar->GetPdgCode() != pdgMother) continue;

    // Daughters
    if(mcPar->GetNDaughters() != 2) continue;
    if(mcPar->GetDaughterFirst() >= nParticles || mcPar->GetDaughterFirst() < 0) continue;

    AliVParticle *mcParDau1 = (AliVParticle*) ev->GetTrack(mcPar->GetDaughterFirst());
    // AliAODMCParticle *mcParDau1 = (AliAODMCParticle*) ev->GetTrack(mcPar->GetDaughterFirst());
    if(!mcParDau1) continue;
    if(TMath::Abs(mcParDau1->PdgCode()) != pdgDaughter) continue;
    // if(TMath::Abs(mcParDau1->GetPdgCode()) != pdgDaughter) continue;

    if(mcPar->GetDaughterLast() >= nParticles || mcPar->GetDaughterLast() < 0) continue;
    AliVParticle *mcParDau2 = (AliVParticle*) ev->GetTrack(mcPar->GetDaughterLast());
    // AliAODMCParticle *mcParDau2 = (AliAODMCParticle*) ev->GetTrack(mcPar->GetDaughterLast());
    if(!mcParDau2) continue;
    if(TMath::Abs(mcParDau2->PdgCode()) != pdgDaughter) continue;
    // if(TMath::Abs(mcParDau2->GetPdgCode()) != pdgDaughter) continue;

    if((TMath::Abs(eta) > TMath::Abs(0.9)) && (TMath::Abs(mcParDau1->Eta()) > TMath::Abs(0.9)) && (TMath::Abs(mcParDau2->Eta()) > TMath::Abs(0.9)))
    {
      numbers[8]++;
      //prompt
      if(mcPar->IsPhysicalPrimary()) numbers[9]++;
      //non-prompt
      if(!mcPar->IsPhysicalPrimary()) numbers[10]++;
    }
  }
}
