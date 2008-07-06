/**************************************************************************
 * Author: Panos Christakoglou.                                           *
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

//-----------------------------------------------------------------
//                 AliProtonAnalysis class
//   This is the class to deal with the proton analysis
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------
#include <Riostream.h>
#include <TFile.h>
#include <TSystem.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TParticle.h>

#include "AliProtonAnalysis.h"

#include <AliExternalTrackParam.h>
#include <AliAODEvent.h>
#include <AliESDEvent.h>
#include <AliLog.h>
#include <AliPID.h>

#include <AliStack.h>

ClassImp(AliProtonAnalysis)

//____________________________________________________________________//
AliProtonAnalysis::AliProtonAnalysis() : 
  TObject(), 
  fNBinsY(0), fMinY(0), fMaxY(0),
  fNBinsPt(0), fMinPt(0), fMaxPt(0),
  fMinTPCClusters(0), fMinITSClusters(0),
  fMaxChi2PerTPCCluster(0), fMaxChi2PerITSCluster(0),
  fMaxCov11(0), fMaxCov22(0), fMaxCov33(0), fMaxCov44(0), fMaxCov55(0),
  fMaxSigmaToVertex(0),
  fMinTPCClustersFlag(kFALSE), fMinITSClustersFlag(kFALSE),
  fMaxChi2PerTPCClusterFlag(kFALSE), fMaxChi2PerITSClusterFlag(kFALSE),
  fMaxCov11Flag(kFALSE), fMaxCov22Flag(kFALSE), fMaxCov33Flag(kFALSE), fMaxCov44Flag(kFALSE), fMaxCov55Flag(kFALSE),
  fMaxSigmaToVertexFlag(kFALSE),
  fITSRefitFlag(kFALSE), fTPCRefitFlag(kFALSE),
  fFunctionProbabilityFlag(kFALSE), 
  fElectronFunction(0), fMuonFunction(0),
  fPionFunction(0), fKaonFunction(0), fProtonFunction(0),
  fUseTPCOnly(kFALSE), fHistEvents(0), fHistYPtProtons(0), fHistYPtAntiProtons(0) {
  //Default constructor
  for(Int_t i = 0; i < 5; i++) fPartFrac[i] = 0.0;
}

//____________________________________________________________________//
AliProtonAnalysis::AliProtonAnalysis(Int_t nbinsY, Float_t fLowY, Float_t fHighY,Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt) : 
  TObject(),
  fNBinsY(nbinsY), fMinY(fLowY), fMaxY(fHighY),
  fNBinsPt(nbinsPt), fMinPt(fLowPt), fMaxPt(fHighPt),
  fMinTPCClusters(0), fMinITSClusters(0),
  fMaxChi2PerTPCCluster(0), fMaxChi2PerITSCluster(0),
  fMaxCov11(0), fMaxCov22(0), fMaxCov33(0), fMaxCov44(0), fMaxCov55(0),
  fMaxSigmaToVertex(0),
  fMinTPCClustersFlag(kFALSE), fMinITSClustersFlag(kFALSE),
  fMaxChi2PerTPCClusterFlag(kFALSE), fMaxChi2PerITSClusterFlag(kFALSE),
  fMaxCov11Flag(kFALSE), fMaxCov22Flag(kFALSE), fMaxCov33Flag(kFALSE), fMaxCov44Flag(kFALSE), fMaxCov55Flag(kFALSE),
  fMaxSigmaToVertexFlag(kFALSE),
  fITSRefitFlag(kFALSE), fTPCRefitFlag(kFALSE),
  fFunctionProbabilityFlag(kFALSE), 
  fElectronFunction(0), fMuonFunction(0),
  fPionFunction(0), fKaonFunction(0), fProtonFunction(0),
  fUseTPCOnly(kFALSE), fHistEvents(0), fHistYPtProtons(0), fHistYPtAntiProtons(0) {
  //Default constructor

  fHistEvents = new TH1I("fHistEvents","Analyzed events",1,0,1);

  fHistYPtProtons = new TH2F("fHistYPtProtons","y-Pt Protons",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
  fHistYPtProtons->SetStats(kTRUE);
  fHistYPtProtons->GetYaxis()->SetTitle("P_{T} [GeV]");
  fHistYPtProtons->GetXaxis()->SetTitle("y");
  fHistYPtProtons->GetXaxis()->SetTitleColor(1);

  fHistYPtAntiProtons = new TH2F("fHistYPtAntiProtons","y-Pt Antiprotons",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
  fHistYPtAntiProtons->SetStats(kTRUE);
  fHistYPtAntiProtons->GetYaxis()->SetTitle("P_{T} [GeV]");
  fHistYPtAntiProtons->GetXaxis()->SetTitle("y");
  fHistYPtAntiProtons->GetXaxis()->SetTitleColor(1);
} 

//____________________________________________________________________//
AliProtonAnalysis::~AliProtonAnalysis() {
  //Default destructor
  
}

//____________________________________________________________________//
void AliProtonAnalysis::InitHistograms(Int_t nbinsY, Float_t fLowY, Float_t fHighY, Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt) {
  fNBinsY = nbinsY;
  fMinY = fLowY;
  fMaxY = fHighY;
  fNBinsPt = nbinsPt;
  fMinPt = fLowPt;
  fMaxPt = fHighPt;

  fHistEvents = new TH1I("fHistEvents","Anallyzed events",1,0,1);

  fHistYPtProtons = new TH2F("fHistYPtProtons","y-Pt Protons",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
  fHistYPtProtons->SetStats(kTRUE);
  fHistYPtProtons->GetYaxis()->SetTitle("P_{T} [GeV]");
  fHistYPtProtons->GetXaxis()->SetTitle("y");
  fHistYPtProtons->GetXaxis()->SetTitleColor(1);

  fHistYPtAntiProtons = new TH2F("fHistYPtAntiProtons","y-Pt Antiprotons",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
  fHistYPtAntiProtons->SetStats(kTRUE);
  fHistYPtAntiProtons->GetYaxis()->SetTitle("P_{T} [GeV]");
  fHistYPtAntiProtons->GetXaxis()->SetTitle("y");
  fHistYPtAntiProtons->GetXaxis()->SetTitleColor(1);
}

//____________________________________________________________________//
Bool_t AliProtonAnalysis::ReadFromFile(const char* filename) {
  Bool_t status = kTRUE;

  TFile *file = TFile::Open(filename);
  if(!file) {
    cout<<"Could not find the input file "<<filename<<endl;
    status = kFALSE;
  }

  TList *list = (TList *)file->Get("outputList1");
  if(list) {
    cout<<"Retrieving objects from the list "<<list->GetName()<<"..."<<endl; 
    fHistYPtProtons = (TH2F *)list->At(0);
    fHistYPtAntiProtons = (TH2F *)list->At(1);
    fHistEvents = (TH1I *)list->At(2);
  }
  else if(!list) {
    cout<<"Retrieving objects from the file... "<<endl;
    fHistYPtProtons = (TH2F *)file->Get("fHistYPtProtons");
    fHistYPtAntiProtons = (TH2F *)file->Get("fHistYPtAntiProtons");
    fHistEvents = (TH1I *)file->Get("fHistEvents");
  }
  if((!fHistYPtProtons)||(!fHistYPtAntiProtons)||(!fHistEvents)) {
    cout<<"Input containers were not found!!!"<<endl;
    status = kFALSE;
  }
  else {
    fHistYPtProtons->Sumw2();
    fHistYPtAntiProtons->Sumw2();
  }

  return status;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetProtonYHistogram() {
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();
  TH1D *fYProtons = (TH1D *)fHistYPtProtons->ProjectionX("fYProtons",0,fHistYPtProtons->GetYaxis()->GetNbins(),"e"); 
  fYProtons->SetStats(kFALSE);
  fYProtons->GetYaxis()->SetTitle("(1/N_{events})(dN/dy)");
  fYProtons->SetTitle("dN/dy protons");
  fYProtons->SetMarkerStyle(kFullCircle);
  fYProtons->SetMarkerColor(4);
  if(nAnalyzedEvents > 0)
    fYProtons->Scale(1./nAnalyzedEvents);

  return fYProtons;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetAntiProtonYHistogram() {
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();
  TH1D *fYAntiProtons = (TH1D *)fHistYPtAntiProtons->ProjectionX("fYAntiProtons",0,fHistYPtAntiProtons->GetYaxis()->GetNbins(),"e"); 
  fYAntiProtons->SetStats(kFALSE);
  fYAntiProtons->GetYaxis()->SetTitle("(1/N_{events})(dN/dy)");
  fYAntiProtons->SetTitle("dN/dy antiprotons");
  fYAntiProtons->SetMarkerStyle(kFullCircle);
  fYAntiProtons->SetMarkerColor(4);
  if(nAnalyzedEvents > 0)
    fYAntiProtons->Scale(1./nAnalyzedEvents);

  return fYAntiProtons;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetProtonPtHistogram() {
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();
  TH1D *fPtProtons = (TH1D *)fHistYPtProtons->ProjectionY("fPtProtons",0,fHistYPtProtons->GetXaxis()->GetNbins(),"e"); 
  fPtProtons->SetStats(kFALSE);
  fPtProtons->GetYaxis()->SetTitle("(1/N_{events})(dN/dP_{T})");
  fPtProtons->SetTitle("dN/dPt protons");
  fPtProtons->SetMarkerStyle(kFullCircle);
  fPtProtons->SetMarkerColor(4);
  if(nAnalyzedEvents > 0)
    fPtProtons->Scale(1./nAnalyzedEvents);

  return fPtProtons;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetAntiProtonPtHistogram() {
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();
  TH1D *fPtAntiProtons = (TH1D *)fHistYPtAntiProtons->ProjectionY("fPtAntiProtons",0,fHistYPtProtons->GetXaxis()->GetNbins(),"e"); 
  fPtAntiProtons->SetStats(kFALSE);
  fPtAntiProtons->GetYaxis()->SetTitle("(1/N_{events})(dN/dP_{T})");
  fPtAntiProtons->SetTitle("dN/dPt antiprotons");
  fPtAntiProtons->SetMarkerStyle(kFullCircle);
  fPtAntiProtons->SetMarkerColor(4);
  if(nAnalyzedEvents > 0)
    fPtAntiProtons->Scale(1./nAnalyzedEvents);

  return fPtAntiProtons;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetYRatioHistogram() {
  TH1D *fYProtons = GetProtonYHistogram();
  TH1D *fYAntiProtons = GetAntiProtonYHistogram();
  
  TH1D *hRatioY = new TH1D("hRatioY","",fYProtons->GetNbinsX(),fYProtons->GetXaxis()->GetXmin(),fYProtons->GetXaxis()->GetXmax());
  hRatioY->Divide(fYAntiProtons,fYProtons,1.0,1.0);
  hRatioY->SetMarkerStyle(kFullCircle);
  hRatioY->SetMarkerColor(4);
  hRatioY->GetYaxis()->SetTitle("#bar{p}/p");
  hRatioY->GetYaxis()->SetTitleOffset(1.4);
  hRatioY->GetXaxis()->SetTitle("y");
  hRatioY->GetXaxis()->SetTitleColor(1);
  hRatioY->SetStats(kFALSE);

  return hRatioY;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetPtRatioHistogram() {
  TH1D *fPtProtons = GetProtonPtHistogram();
  TH1D *fPtAntiProtons = GetAntiProtonPtHistogram();
  
  TH1D *hRatioPt = new TH1D("hRatioPt","",fPtProtons->GetNbinsX(),fPtProtons->GetXaxis()->GetXmin(),fPtProtons->GetXaxis()->GetXmax());
  hRatioPt->Divide(fPtAntiProtons,fPtProtons,1.0,1.0);
  hRatioPt->SetMarkerStyle(kFullCircle);
  hRatioPt->SetMarkerColor(4);
  hRatioPt->GetYaxis()->SetTitle("#bar{p}/p");
  hRatioPt->GetYaxis()->SetTitleOffset(1.4);
  hRatioPt->GetXaxis()->SetTitle("P_{T} [GeV/c]");
  hRatioPt->GetXaxis()->SetTitleColor(1);
  hRatioPt->SetStats(kFALSE);

  return hRatioPt;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetYAsymmetryHistogram() {
  TH1D *fYProtons = GetProtonYHistogram();
  TH1D *fYAntiProtons = GetAntiProtonYHistogram();
  
  TH1D *hsum = new TH1D("hsumY","",fYProtons->GetNbinsX(),fYProtons->GetXaxis()->GetXmin(),fYProtons->GetXaxis()->GetXmax());
  hsum->Add(fYProtons,fYAntiProtons,1.0,1.0);

  TH1D *hdiff = new TH1D("hdiffY","",fYProtons->GetNbinsX(),fYProtons->GetXaxis()->GetXmin(),fYProtons->GetXaxis()->GetXmax());
  hdiff->Add(fYProtons,fYAntiProtons,1.0,-1.0);

  TH1D *hAsymmetryY = new TH1D("hAsymmetryY","",fYProtons->GetNbinsX(),fYProtons->GetXaxis()->GetXmin(),fYProtons->GetXaxis()->GetXmax());
  hAsymmetryY->Divide(hdiff,hsum,2.0,1.);
  hAsymmetryY->SetMarkerStyle(kFullCircle);
  hAsymmetryY->SetMarkerColor(4);
  hAsymmetryY->GetYaxis()->SetTitle("A_{p}");
  hAsymmetryY->GetYaxis()->SetTitleOffset(1.4);
  hAsymmetryY->GetXaxis()->SetTitle("y");
  hAsymmetryY->GetXaxis()->SetTitleColor(1);
  hAsymmetryY->SetStats(kFALSE);

  return hAsymmetryY;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetPtAsymmetryHistogram() {
  TH1D *fPtProtons = GetProtonPtHistogram();
  TH1D *fPtAntiProtons = GetAntiProtonPtHistogram();
  
  TH1D *hsum = new TH1D("hsumPt","",fPtProtons->GetNbinsX(),fPtProtons->GetXaxis()->GetXmin(),fPtProtons->GetXaxis()->GetXmax());
  hsum->Add(fPtProtons,fPtAntiProtons,1.0,1.0);

  TH1D *hdiff = new TH1D("hdiffPt","",fPtProtons->GetNbinsX(),fPtProtons->GetXaxis()->GetXmin(),fPtProtons->GetXaxis()->GetXmax());
  hdiff->Add(fPtProtons,fPtAntiProtons,1.0,-1.0);

  TH1D *hAsymmetryPt = new TH1D("hAsymmetryPt","",fPtProtons->GetNbinsX(),fPtProtons->GetXaxis()->GetXmin(),fPtProtons->GetXaxis()->GetXmax());
  hAsymmetryPt->Divide(hdiff,hsum,2.0,1.);
  hAsymmetryPt->SetMarkerStyle(kFullCircle);
  hAsymmetryPt->SetMarkerColor(4);
  hAsymmetryPt->GetYaxis()->SetTitle("A_{p}");
  hAsymmetryPt->GetYaxis()->SetTitleOffset(1.4);
  hAsymmetryPt->GetXaxis()->SetTitle("P_{T} [GeV/c]");
  hAsymmetryPt->GetXaxis()->SetTitleColor(1);
  hAsymmetryPt->SetStats(kFALSE);

  return hAsymmetryPt;
}

//____________________________________________________________________//
Double_t AliProtonAnalysis::GetParticleFraction(Int_t i, Double_t p) {
  Double_t partFrac=0;
  if(fFunctionProbabilityFlag) {
    if(i == 0) partFrac = fElectronFunction->Eval(p);
    if(i == 1) partFrac = fMuonFunction->Eval(p);
    if(i == 2) partFrac = fPionFunction->Eval(p);
    if(i == 3) partFrac = fKaonFunction->Eval(p);
    if(i == 4) partFrac = fProtonFunction->Eval(p);
  }
  else partFrac = fPartFrac[i];

  return partFrac;
}

//____________________________________________________________________//
void AliProtonAnalysis::Analyze(AliESDEvent* fESD) {
  //Main analysis part - ESD
  fHistEvents->Fill(0); //number of analyzed events
  Double_t Pt = 0.0, P = 0.0;
  Int_t nGoodTracks = fESD->GetNumberOfTracks();
  for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    Double_t probability[5];

    if(IsAccepted(track)) {	
      if(fUseTPCOnly) {
	AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
	if(!tpcTrack) continue;
	Pt = tpcTrack->Pt();
	P = tpcTrack->P();
	
	//pid
	track->GetTPCpid(probability);
	Double_t rcc = 0.0;
	for(Int_t i = 0; i < AliPID::kSPECIES; i++) 
	  rcc += probability[i]*GetParticleFraction(i,P);
	if(rcc == 0.0) continue;
	Double_t w[5];
	for(Int_t i = 0; i < AliPID::kSPECIES; i++) 
	  w[i] = probability[i]*GetParticleFraction(i,P)/rcc;
	Long64_t fParticleType = TMath::LocMax(AliPID::kSPECIES,w);
	if(fParticleType == 4) {
	  if(tpcTrack->Charge() > 0) 
	    fHistYPtProtons->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz()),Pt);
	  else if(tpcTrack->Charge() < 0) 
	    fHistYPtAntiProtons->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz()),Pt);
	}//proton check
      }//TPC only tracks
      else if(!fUseTPCOnly) {
	Pt = track->Pt();
	P = track->P();
	
	//pid
	track->GetESDpid(probability);
	Double_t rcc = 0.0;
	for(Int_t i = 0; i < AliPID::kSPECIES; i++) 
	  rcc += probability[i]*GetParticleFraction(i,P);
	if(rcc == 0.0) continue;
	Double_t w[5];
	for(Int_t i = 0; i < AliPID::kSPECIES; i++) 
	  w[i] = probability[i]*GetParticleFraction(i,P)/rcc;
	Long64_t fParticleType = TMath::LocMax(AliPID::kSPECIES,w);
	if(fParticleType == 4) {
	  //cout<<"(Anti)protons found..."<<endl;
	  if(track->Charge() > 0) 
	    fHistYPtProtons->Fill(Rapidity(track->Px(),track->Py(),track->Pz()),Pt);
	  else if(track->Charge() < 0) 
	    fHistYPtAntiProtons->Fill(Rapidity(track->Px(),track->Py(),track->Pz()),Pt);
	}//proton check 
      }//combined tracking
    }//cuts
  }//track loop 
}

//____________________________________________________________________//
void AliProtonAnalysis::Analyze(AliAODEvent* fAOD) {
  //Main analysis part - AOD
  fHistEvents->Fill(0); //number of analyzed events
  Int_t nGoodTracks = fAOD->GetNumberOfTracks();
  for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {
    AliAODTrack* track = fAOD->GetTrack(iTracks);
    Double_t Pt = track->Pt();
    Double_t P = track->P();
    
    //pid
    Double_t probability[10];
    track->GetPID(probability);
    Double_t rcc = 0.0;
    for(Int_t i = 0; i < AliPID::kSPECIESN; i++) rcc += probability[i]*GetParticleFraction(i,P);
    if(rcc == 0.0) continue;
    Double_t w[10];
    for(Int_t i = 0; i < AliPID::kSPECIESN; i++) w[i] = probability[i]*GetParticleFraction(i,P)/rcc;
    Long64_t fParticleType = TMath::LocMax(AliPID::kSPECIESN,w);
    if(fParticleType == 4) {
      if(track->Charge() > 0) fHistYPtProtons->Fill(track->Y(fParticleType),Pt);
      else if(track->Charge() < 0) fHistYPtAntiProtons->Fill(track->Y(fParticleType),Pt);
    }//proton check
  }//track loop 
}

//____________________________________________________________________//
void AliProtonAnalysis::Analyze(AliStack* stack) {
  //Main analysis part - MC
  fHistEvents->Fill(0); //number of analyzed events
  for(Int_t i = 0; i < stack->GetNprimary(); i++) {
    TParticle *particle = stack->Particle(i);
    if(particle->Pt() < 0.1) continue;
    if(TMath::Abs(particle->Eta()) > 1.0) continue;
    Int_t pdgcode = particle->GetPdgCode();
    if(pdgcode == 2212) fHistYPtProtons->Fill(Rapidity(particle->Px(),
						       particle->Py(),
						       particle->Pz()),
					      particle->Pt());
    if(pdgcode == -2212) fHistYPtAntiProtons->Fill(Rapidity(particle->Px(),
							    particle->Py(),
							    particle->Pz()),
						   particle->Pt());
  }//particle loop                                                                  
}

//____________________________________________________________________//
Bool_t AliProtonAnalysis::IsAccepted(AliESDtrack* track) {
  // Checks if the track is excluded from the cuts
  Double_t Pt = 0.0, Px = 0.0, Py = 0.0, Pz = 0.0;
  if(fUseTPCOnly) {
    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(!tpcTrack) {
      Pt = 0.0; Px = 0.0; Py = 0.0; Pz = 0.0;
    }
    else {
      Pt = tpcTrack->Pt();
      Px = tpcTrack->Px();
      Py = tpcTrack->Py();
      Pz = tpcTrack->Pz();
    }
  }
  else{
    Pt = track->Pt();
    Px = track->Px();
    Py = track->Py();
    Pz = track->Pz();
  }
     
  Int_t  fIdxInt[200];
  Int_t nClustersITS = track->GetITSclusters(fIdxInt);
  Int_t nClustersTPC = track->GetTPCclusters(fIdxInt);

  Float_t chi2PerClusterITS = -1;
  Float_t chi2PerClusterTPC = -1;
  if (nClustersTPC!=0)
    chi2PerClusterTPC = track->GetTPCchi2()/Float_t(nClustersTPC);

  Double_t extCov[15];
  track->GetExternalCovariance(extCov);

  if(fMinITSClustersFlag)
    if(nClustersITS < fMinITSClusters) return kFALSE;
  if(fMinTPCClustersFlag)
    if(nClustersTPC < fMinTPCClusters) return kFALSE;
  if(fMaxChi2PerTPCClusterFlag)
    if(chi2PerClusterTPC > fMaxChi2PerTPCCluster) return kFALSE; 
  if(fMaxChi2PerITSClusterFlag)
    if(chi2PerClusterITS > fMaxChi2PerITSCluster) return kFALSE; 
  if(fMaxCov11Flag)
    if(extCov[0] > fMaxCov11) return kFALSE;
  if(fMaxCov22Flag)
    if(extCov[2] > fMaxCov22) return kFALSE;
  if(fMaxCov33Flag)
    if(extCov[5] > fMaxCov33) return kFALSE;
  if(fMaxCov44Flag)
    if(extCov[9] > fMaxCov44) return kFALSE;
  if(fMaxCov55Flag)
    if(extCov[14] > fMaxCov55) return kFALSE;
  if(fMaxSigmaToVertexFlag)
    if(GetSigmaToVertex(track) > fMaxSigmaToVertex) return kFALSE;
  if(fITSRefitFlag)
    if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) return kFALSE;
  if(fTPCRefitFlag)
    if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) return kFALSE;

  if((Pt < fMinPt) || (Pt > fMaxPt)) return kFALSE;
  if((Rapidity(Px,Py,Pz) < fMinY) || (Rapidity(Px,Py,Pz) > fMaxY)) return kFALSE;

  return kTRUE;
}

//____________________________________________________________________//
Float_t AliProtonAnalysis::GetSigmaToVertex(AliESDtrack* esdTrack) {
  // Calculates the number of sigma to the vertex.
  
  Float_t b[2];
  Float_t bRes[2];
  Float_t bCov[3];
  esdTrack->GetImpactParameters(b,bCov);
  if (bCov[0]<=0 || bCov[2]<=0) {
    //AliDebug(1, "Estimated b resolution lower or equal zero!");
    bCov[0]=0; bCov[2]=0;
  }
  bRes[0] = TMath::Sqrt(bCov[0]);
  bRes[1] = TMath::Sqrt(bCov[2]);
  
  if (bRes[0] == 0 || bRes[1] ==0) return -1;
  
  Float_t d = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));
  
  if (TMath::Exp(-d * d / 2) < 1e-10) return 1000;
  
  d = TMath::ErfInverse(1 - TMath::Exp(-d * d / 2)) * TMath::Sqrt(2);
  
  return d;
}

//____________________________________________________________________//
Double_t AliProtonAnalysis::Rapidity(Double_t Px, Double_t Py, Double_t Pz) {
  //returns the rapidity of the proton - to be removed
  Double_t fMass = 9.38270000000000048e-01;
  
  Double_t P = TMath::Sqrt(TMath::Power(Px,2) + 
                           TMath::Power(Py,2) + 
			   TMath::Power(Pz,2));
  Double_t energy = TMath::Sqrt(P*P + fMass*fMass);
  Double_t y = -999;
  if(energy != Pz) 
    y = 0.5*TMath::Log((energy + Pz)/(energy - Pz));

  return y;
}

//____________________________________________________________________//
Bool_t AliProtonAnalysis::PrintMean(TH1 *hist, Double_t edge) {
  //calculates the mean value of the ratio/asymmetry within \pm edge
  Double_t sum = 0.0;
  Int_t nentries = 0;
  //calculate the mean
  for(Int_t i = 0; i < hist->GetXaxis()->GetNbins(); i++) {
    Double_t x = hist->GetBinCenter(i+1);
    Double_t y = hist->GetBinContent(i+1);
    if(TMath::Abs(x) < edge) {
      sum += y;
      nentries += 1;
    }
  }
  Double_t mean = 0.0;
  if(nentries != 0)
    mean = sum/nentries;

  //calculate the error
  for(Int_t i = 0; i < hist->GetXaxis()->GetNbins(); i++) {
    Double_t x = hist->GetBinCenter(i+1);
    Double_t y = hist->GetBinContent(i+1);
    if(TMath::Abs(x) < edge) {
      sum += TMath::Power((mean - y),2);
      nentries += 1;
    }
  }

  Double_t error = 0.0;
  if(nentries != 0)
    error =  TMath::Sqrt(sum)/nentries;

  cout<<"========================================="<<endl;
  cout<<"Input distribution: "<<hist->GetName()<<endl;
  cout<<"Interval used: -"<<edge<<" -> "<<edge<<endl;
  cout<<"Mean value :"<<mean<<endl;
  cout<<"Error: "<<error<<endl;
  cout<<"========================================="<<endl;

  return 0;
}

//____________________________________________________________________//
Bool_t AliProtonAnalysis::PrintYields(TH1 *hist, Double_t edge) {
  //calculates the (anti)proton yields within the \pm edge
  Double_t sum = 0.0, sumerror = 0.0;
  Double_t error = 0.0;
  for(Int_t i = 0; i < hist->GetXaxis()->GetNbins(); i++) {
    Double_t x = hist->GetBinCenter(i+1);
    Double_t y = hist->GetBinContent(i+1);
    if(TMath::Abs(x) < edge) {
      sum += y;
      sumerror += TMath::Power(hist->GetBinError(i+1),2); 
    }
  }

  error = TMath::Sqrt(sumerror);

  cout<<"========================================="<<endl;
  cout<<"Input distribution: "<<hist->GetName()<<endl;
  cout<<"Interval used: -"<<edge<<" -> "<<edge<<endl;
  cout<<"Yields :"<<sum<<endl;
  cout<<"Error: "<<error<<endl;
  cout<<"========================================="<<endl;

  return 0;
}



