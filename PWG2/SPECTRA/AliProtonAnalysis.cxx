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
#include <AliCFContainer.h>
#include <AliCFEffGrid.h>
#include <AliCFEffGrid.h>

ClassImp(AliProtonAnalysis)

//____________________________________________________________________//
AliProtonAnalysis::AliProtonAnalysis() : 
  TObject(), 
  fNBinsY(0), fMinY(0), fMaxY(0),
  fNBinsPt(0), fMinPt(0), fMaxPt(0),
  fMinTPCClusters(0), fMinITSClusters(0),
  fMaxChi2PerTPCCluster(0), fMaxChi2PerITSCluster(0),
  fMaxCov11(0), fMaxCov22(0), fMaxCov33(0), fMaxCov44(0), fMaxCov55(0),
  fMaxSigmaToVertex(0), fMaxSigmaToVertexTPC(0),
  fMinTPCClustersFlag(kFALSE), fMinITSClustersFlag(kFALSE),
  fMaxChi2PerTPCClusterFlag(kFALSE), fMaxChi2PerITSClusterFlag(kFALSE),
  fMaxCov11Flag(kFALSE), fMaxCov22Flag(kFALSE), fMaxCov33Flag(kFALSE), fMaxCov44Flag(kFALSE), fMaxCov55Flag(kFALSE),
  fMaxSigmaToVertexFlag(kFALSE), fMaxSigmaToVertexTPCFlag(kFALSE),
  fITSRefitFlag(kFALSE), fTPCRefitFlag(kFALSE),
  fESDpidFlag(kFALSE), fTPCpidFlag(kFALSE),
  fQAHistograms(kFALSE),
  fGlobalQAList(0), fQA2DList(0),
  fQAPrimaryProtonsAcceptedList(0),
  fQAPrimaryProtonsRejectedList(0),
  fQASecondaryProtonsAcceptedList(0),
  fQASecondaryProtonsRejectedList(0),
  fQAPrimaryAntiProtonsAcceptedList(0),
  fQAPrimaryAntiProtonsRejectedList(0),
  fQASecondaryAntiProtonsAcceptedList(0),
  fQASecondaryAntiProtonsRejectedList(0),
  fFunctionProbabilityFlag(kFALSE), 
  fElectronFunction(0), fMuonFunction(0),
  fPionFunction(0), fKaonFunction(0), fProtonFunction(0),
  fUseTPCOnly(kFALSE), fHistEvents(0), fHistYPtProtons(0), fHistYPtAntiProtons(0),
  fCorrectionList2D(0), fEfficiencyList1D(0), fCorrectionList1D(0) {
  //Default constructor
  for(Int_t i = 0; i < 5; i++) fPartFrac[i] = 0.0;
  fCorrectionList2D = new TList(); 
  fEfficiencyList1D = new TList(); 
  fCorrectionList1D = new TList();
}

//____________________________________________________________________//
AliProtonAnalysis::AliProtonAnalysis(Int_t nbinsY, Float_t fLowY, Float_t fHighY,Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt) : 
  TObject(),
  fNBinsY(nbinsY), fMinY(fLowY), fMaxY(fHighY),
  fNBinsPt(nbinsPt), fMinPt(fLowPt), fMaxPt(fHighPt),
  fMinTPCClusters(0), fMinITSClusters(0),
  fMaxChi2PerTPCCluster(0), fMaxChi2PerITSCluster(0),
  fMaxCov11(0), fMaxCov22(0), fMaxCov33(0), fMaxCov44(0), fMaxCov55(0),
  fMaxSigmaToVertex(0), fMaxSigmaToVertexTPC(0),
  fMinTPCClustersFlag(kFALSE), fMinITSClustersFlag(kFALSE),
  fMaxChi2PerTPCClusterFlag(kFALSE), fMaxChi2PerITSClusterFlag(kFALSE),
  fMaxCov11Flag(kFALSE), fMaxCov22Flag(kFALSE), fMaxCov33Flag(kFALSE), fMaxCov44Flag(kFALSE), fMaxCov55Flag(kFALSE),
  fMaxSigmaToVertexFlag(kFALSE), fMaxSigmaToVertexTPCFlag(kFALSE),
  fITSRefitFlag(kFALSE), fTPCRefitFlag(kFALSE),
  fESDpidFlag(kFALSE), fTPCpidFlag(kFALSE),
  fQAHistograms(kFALSE), 
  fGlobalQAList(0), fQA2DList(0),
  fQAPrimaryProtonsAcceptedList(0),
  fQAPrimaryProtonsRejectedList(0),
  fQASecondaryProtonsAcceptedList(0),
  fQASecondaryProtonsRejectedList(0),
  fQAPrimaryAntiProtonsAcceptedList(0),
  fQAPrimaryAntiProtonsRejectedList(0),
  fQASecondaryAntiProtonsAcceptedList(0),
  fQASecondaryAntiProtonsRejectedList(0),
  fFunctionProbabilityFlag(kFALSE), 
  fElectronFunction(0), fMuonFunction(0),
  fPionFunction(0), fKaonFunction(0), fProtonFunction(0),
  fUseTPCOnly(kFALSE), fHistEvents(0), fHistYPtProtons(0), fHistYPtAntiProtons(0),
  fCorrectionList2D(0), fEfficiencyList1D(0), fCorrectionList1D(0) {
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
  if(fHistEvents) delete fHistEvents;
  if(fHistYPtProtons) delete fHistYPtProtons;
  if(fHistYPtAntiProtons) delete fHistYPtAntiProtons;
  if(fCorrectionList2D) delete fCorrectionList2D;
  if(fEfficiencyList1D) delete fEfficiencyList1D;
  if(fCorrectionList1D) delete fCorrectionList1D;
  if(fGlobalQAList) delete fGlobalQAList;
  if(fQA2DList) delete fQA2DList;
  if(fQAPrimaryProtonsAcceptedList) delete fQAPrimaryProtonsAcceptedList;
  if(fQAPrimaryProtonsRejectedList) delete fQAPrimaryProtonsRejectedList;
  if(fQASecondaryProtonsAcceptedList) delete fQASecondaryProtonsAcceptedList;
  if(fQASecondaryProtonsRejectedList) delete fQASecondaryProtonsRejectedList;
  if(fQAPrimaryAntiProtonsAcceptedList) delete fQAPrimaryAntiProtonsAcceptedList;
  if(fQAPrimaryAntiProtonsRejectedList) delete fQAPrimaryAntiProtonsRejectedList;
  if(fQASecondaryAntiProtonsAcceptedList) delete fQASecondaryAntiProtonsAcceptedList;
  if(fQASecondaryAntiProtonsRejectedList) delete fQASecondaryAntiProtonsRejectedList; 
}

//____________________________________________________________________//
void AliProtonAnalysis::InitAnalysisHistograms(Int_t nbinsY, Float_t fLowY, Float_t fHighY, 
					       Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt) {
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
  if (nClustersITS!=0)
    chi2PerClusterITS = track->GetITSchi2()/Float_t(nClustersITS);
  Float_t chi2PerClusterTPC = -1;
  if (nClustersTPC!=0)
    chi2PerClusterTPC = track->GetTPCchi2()/Float_t(nClustersTPC);

  Double_t extCov[15];
  track->GetExternalCovariance(extCov);

  if(fMinITSClustersFlag)
    if(nClustersITS < fMinITSClusters) return kFALSE;
  if(fMaxChi2PerITSClusterFlag)
    if(chi2PerClusterITS > fMaxChi2PerITSCluster) return kFALSE; 
  if(fMinTPCClustersFlag)
    if(nClustersTPC < fMinTPCClusters) return kFALSE;
  if(fMaxChi2PerTPCClusterFlag)
    if(chi2PerClusterTPC > fMaxChi2PerTPCCluster) return kFALSE; 
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
  if(fMaxSigmaToVertexTPCFlag)
    if(GetSigmaToVertex(track) > fMaxSigmaToVertexTPC) return kFALSE;
  if(fITSRefitFlag)
    if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) return kFALSE;
  if(fTPCRefitFlag)
    if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) return kFALSE;
  if(fESDpidFlag)
    if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) return kFALSE;
  if(fTPCpidFlag)
    if ((track->GetStatus() & AliESDtrack::kTPCpid) == 0) return kFALSE;

  if((Pt < fMinPt) || (Pt > fMaxPt)) return kFALSE;
  if((Rapidity(Px,Py,Pz) < fMinY) || (Rapidity(Px,Py,Pz) > fMaxY)) return kFALSE;

  return kTRUE;
}

//____________________________________________________________________//
Bool_t AliProtonAnalysis::IsAccepted(AliESDtrack* track, AliStack *stack) {
  // Checks if the track is excluded from the cuts
  Bool_t status = kTRUE;
  Int_t nPrimaries = stack->GetNprimary();
  Int_t label = TMath::Abs(track->GetLabel());

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
  if (nClustersITS!=0)
    chi2PerClusterITS = track->GetITSchi2()/Float_t(nClustersITS);
  Float_t chi2PerClusterTPC = -1;
  if (nClustersTPC!=0)
    chi2PerClusterTPC = track->GetTPCchi2()/Float_t(nClustersTPC);

  Double_t extCov[15];
  track->GetExternalCovariance(extCov);
  
  //protons
  if(track->Charge() > 0) {
    //Primaries
    if(label <= nPrimaries) {
      if(fMinITSClustersFlag) {
	if(nClustersITS < fMinITSClusters) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(0)))->Fill(nClustersITS);
	  status = kFALSE;
	}
	else if(nClustersITS >= fMinITSClusters) 
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(0)))->Fill(nClustersITS);
      }//ITS clusters
      if(fMaxChi2PerITSClusterFlag) {
	if(chi2PerClusterITS > fMaxChi2PerITSCluster) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(1)))->Fill(chi2PerClusterITS);
	  status = kFALSE;
	}
	else if(chi2PerClusterITS <= fMaxChi2PerITSCluster)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(1)))->Fill(chi2PerClusterITS);
      }//chi2 per ITS cluster
      if(fMinTPCClustersFlag) {
	if(nClustersTPC < fMinTPCClusters) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(2)))->Fill(nClustersTPC);
	  status = kFALSE;
	}
	else if(nClustersTPC >= fMinTPCClusters)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(2)))->Fill(nClustersTPC);
      }//TPC clusters
      if(fMaxChi2PerTPCClusterFlag) {
	if(chi2PerClusterTPC > fMaxChi2PerTPCCluster) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(3)))->Fill(chi2PerClusterTPC);
	  status = kFALSE;
	}
	else if(chi2PerClusterTPC <= fMaxChi2PerTPCCluster)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(3)))->Fill(chi2PerClusterTPC);
      }//chi2 per TPC cluster
      if(fMaxCov11Flag) {
	if(extCov[0] > fMaxCov11) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(4)))->Fill(extCov[0]);
	  status = kFALSE;
	}
	else if(extCov[0] <= fMaxCov11)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(4)))->Fill(extCov[0]);
      }//cov11
      if(fMaxCov22Flag) {
	if(extCov[2] > fMaxCov22) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(5)))->Fill(extCov[2]);
	  status = kFALSE;
	}
	else if(extCov[2] <= fMaxCov22)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(5)))->Fill(extCov[2]);
      }//cov11
      if(fMaxCov33Flag) {
	if(extCov[5] > fMaxCov33) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(6)))->Fill(extCov[5]);
	  status = kFALSE;
	}
	else if(extCov[5] <= fMaxCov33)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(6)))->Fill(extCov[5]);
      }//cov11
      if(fMaxCov44Flag) {
	if(extCov[9] > fMaxCov44) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(7)))->Fill(extCov[9]);
	  status = kFALSE;
	}
	else if(extCov[9] <= fMaxCov44)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(7)))->Fill(extCov[9]);
      }//cov11
      if(fMaxCov55Flag) {
	if(extCov[14] > fMaxCov55) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(8)))->Fill(extCov[14]);
	  status = kFALSE;
	}
	else if(extCov[14] <= fMaxCov55)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(8)))->Fill(extCov[14]);
      }//cov55
      if(fMaxSigmaToVertexFlag) {
	if(GetSigmaToVertex(track) > fMaxSigmaToVertex) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(9)))->Fill(GetSigmaToVertex(track));
	  status = kFALSE;
	}
	else if(GetSigmaToVertex(track) <= fMaxSigmaToVertex)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(9)))->Fill(GetSigmaToVertex(track));
      }//sigma to vertex
      if(fMaxSigmaToVertexTPCFlag) {
	if(GetSigmaToVertex(track) > fMaxSigmaToVertexTPC) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(10)))->Fill(GetSigmaToVertex(track));
	  status = kFALSE;
	}
	else if(GetSigmaToVertex(track) <= fMaxSigmaToVertexTPC)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(10)))->Fill(GetSigmaToVertex(track));
      }//sigma to vertex TPC
      if(fITSRefitFlag) {
	if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(11)))->Fill(0);
	status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kITSrefit) != 0)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(11)))->Fill(0);
      }//ITS refit
      if(fTPCRefitFlag) {
	if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(12)))->Fill(0);
	  status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kTPCrefit) != 0)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(12)))->Fill(0);
      }//TPC refit
      if(fESDpidFlag) {
	if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(13)))->Fill(0);
	  status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kESDpid) != 0)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(13)))->Fill(0);
      }//ESD pid
      if(fTPCpidFlag) {
	if ((track->GetStatus() & AliESDtrack::kTPCpid) == 0) {
	  ((TH1F *)(fQAPrimaryProtonsRejectedList->At(13)))->Fill(0);
	  status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kTPCpid) != 0)
	  ((TH1F *)(fQAPrimaryProtonsAcceptedList->At(13)))->Fill(0);
      }//TPC pid
    }//primary particle cut

    //Secondaries
    if(label > nPrimaries) {
      if(fMinITSClustersFlag) {
	if(nClustersITS < fMinITSClusters) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(0)))->Fill(nClustersITS);
	  status = kFALSE;
	}
	else if(nClustersITS >= fMinITSClusters) 
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(0)))->Fill(nClustersITS);
      }//ITS clusters
      if(fMaxChi2PerITSClusterFlag) {
	if(chi2PerClusterITS > fMaxChi2PerITSCluster) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(1)))->Fill(chi2PerClusterITS);
	  status = kFALSE;
	}
	else if(chi2PerClusterITS <= fMaxChi2PerITSCluster)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(1)))->Fill(chi2PerClusterITS);
      }//chi2 per ITS cluster
      if(fMinTPCClustersFlag) {
	if(nClustersTPC < fMinTPCClusters) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(2)))->Fill(nClustersTPC);
	  status = kFALSE;
	}
	else if(nClustersTPC >= fMinTPCClusters)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(2)))->Fill(nClustersTPC);
      }//TPC clusters
      if(fMaxChi2PerTPCClusterFlag) {
	if(chi2PerClusterTPC > fMaxChi2PerTPCCluster) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(3)))->Fill(chi2PerClusterTPC);
	  status = kFALSE;
	}
	else if(chi2PerClusterTPC <= fMaxChi2PerTPCCluster)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(3)))->Fill(chi2PerClusterTPC);
      }//chi2 per TPC cluster
      if(fMaxCov11Flag) {
	if(extCov[0] > fMaxCov11) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(4)))->Fill(extCov[0]);
	  status = kFALSE;
	}
	else if(extCov[0] <= fMaxCov11)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(4)))->Fill(extCov[0]);
      }//cov11
      if(fMaxCov22Flag) {
	if(extCov[2] > fMaxCov22) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(5)))->Fill(extCov[2]);
	  status = kFALSE;
	}
	else if(extCov[2] <= fMaxCov22)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(5)))->Fill(extCov[2]);
      }//cov11
      if(fMaxCov33Flag) {
	if(extCov[5] > fMaxCov33) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(6)))->Fill(extCov[5]);
	  status = kFALSE;
	}
	else if(extCov[5] <= fMaxCov33)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(6)))->Fill(extCov[5]);
      }//cov11
      if(fMaxCov44Flag) {
	if(extCov[9] > fMaxCov44) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(7)))->Fill(extCov[9]);
	  status = kFALSE;
	}
	else if(extCov[9] <= fMaxCov44)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(7)))->Fill(extCov[9]);
      }//cov11
      if(fMaxCov55Flag) {
	if(extCov[14] > fMaxCov55) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(8)))->Fill(extCov[14]);
	  status = kFALSE;
	}
	else if(extCov[14] <= fMaxCov55)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(8)))->Fill(extCov[14]);
      }//cov55
      if(fMaxSigmaToVertexFlag) {
	if(GetSigmaToVertex(track) > fMaxSigmaToVertex) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(9)))->Fill(GetSigmaToVertex(track));
	  status = kFALSE;
	}
	else if(GetSigmaToVertex(track) <= fMaxSigmaToVertex)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(9)))->Fill(GetSigmaToVertex(track));
      }//sigma to vertex
      if(fMaxSigmaToVertexTPCFlag) {
	if(GetSigmaToVertex(track) > fMaxSigmaToVertexTPC) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(10)))->Fill(GetSigmaToVertex(track));
	  status = kFALSE;
	}
	else if(GetSigmaToVertex(track) <= fMaxSigmaToVertexTPC)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(10)))->Fill(GetSigmaToVertex(track));
      }//sigma to vertex TPC
      if(fITSRefitFlag) {
	if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(11)))->Fill(0);
	status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kITSrefit) != 0)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(11)))->Fill(0);
      }//ITS refit
      if(fTPCRefitFlag) {
	if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(12)))->Fill(0);
	  status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kTPCrefit) != 0)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(12)))->Fill(0);
      }//TPC refit
      if(fESDpidFlag) {
	if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(13)))->Fill(0);
	  status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kESDpid) != 0)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(13)))->Fill(0);
      }//ESD pid
      if(fTPCpidFlag) {
	if ((track->GetStatus() & AliESDtrack::kTPCpid) == 0) {
	  ((TH1F *)(fQASecondaryProtonsRejectedList->At(13)))->Fill(0);
	  status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kTPCpid) != 0)
	  ((TH1F *)(fQASecondaryProtonsAcceptedList->At(13)))->Fill(0);
      }//TPC pid
    }//secondary particle cut
  }//protons

  //antiprotons
  if(track->Charge() < 0) {
    //Primaries
    if(label <= nPrimaries) {
      if(fMinITSClustersFlag) {
	if(nClustersITS < fMinITSClusters) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(0)))->Fill(nClustersITS);
	  status = kFALSE;
	}
	else if(nClustersITS >= fMinITSClusters) 
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(0)))->Fill(nClustersITS);
      }//ITS clusters
      if(fMaxChi2PerITSClusterFlag) {
	if(chi2PerClusterITS > fMaxChi2PerITSCluster) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(1)))->Fill(chi2PerClusterITS);
	  status = kFALSE;
	}
	else if(chi2PerClusterITS <= fMaxChi2PerITSCluster)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(1)))->Fill(chi2PerClusterITS);
      }//chi2 per ITS cluster
      if(fMinTPCClustersFlag) {
	if(nClustersTPC < fMinTPCClusters) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(2)))->Fill(nClustersTPC);
	  status = kFALSE;
	}
	else if(nClustersTPC >= fMinTPCClusters)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(2)))->Fill(nClustersTPC);
      }//TPC clusters
      if(fMaxChi2PerTPCClusterFlag) {
	if(chi2PerClusterTPC > fMaxChi2PerTPCCluster) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(3)))->Fill(chi2PerClusterTPC);
	  status = kFALSE;
	}
	else if(chi2PerClusterTPC <= fMaxChi2PerTPCCluster)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(3)))->Fill(chi2PerClusterTPC);
      }//chi2 per TPC cluster
      if(fMaxCov11Flag) {
	if(extCov[0] > fMaxCov11) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(4)))->Fill(extCov[0]);
	  status = kFALSE;
	}
	else if(extCov[0] <= fMaxCov11)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(4)))->Fill(extCov[0]);
      }//cov11
      if(fMaxCov22Flag) {
	if(extCov[2] > fMaxCov22) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(5)))->Fill(extCov[2]);
	  status = kFALSE;
	}
	else if(extCov[2] <= fMaxCov22)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(5)))->Fill(extCov[2]);
      }//cov11
      if(fMaxCov33Flag) {
	if(extCov[5] > fMaxCov33) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(6)))->Fill(extCov[5]);
	  status = kFALSE;
	}
	else if(extCov[5] <= fMaxCov33)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(6)))->Fill(extCov[5]);
      }//cov11
      if(fMaxCov44Flag) {
	if(extCov[9] > fMaxCov44) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(7)))->Fill(extCov[9]);
	  status = kFALSE;
	}
	else if(extCov[9] <= fMaxCov44)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(7)))->Fill(extCov[9]);
      }//cov11
      if(fMaxCov55Flag) {
	if(extCov[14] > fMaxCov55) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(8)))->Fill(extCov[14]);
	  status = kFALSE;
	}
	else if(extCov[14] <= fMaxCov55)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(8)))->Fill(extCov[14]);
      }//cov55
      if(fMaxSigmaToVertexFlag) {
	if(GetSigmaToVertex(track) > fMaxSigmaToVertex) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(9)))->Fill(GetSigmaToVertex(track));
	  status = kFALSE;
	}
	else if(GetSigmaToVertex(track) <= fMaxSigmaToVertex)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(9)))->Fill(GetSigmaToVertex(track));
      }//sigma to vertex
      if(fMaxSigmaToVertexTPCFlag) {
	if(GetSigmaToVertex(track) > fMaxSigmaToVertexTPC) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(10)))->Fill(GetSigmaToVertex(track));
	  status = kFALSE;
	}
	else if(GetSigmaToVertex(track) <= fMaxSigmaToVertexTPC)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(10)))->Fill(GetSigmaToVertex(track));
      }//sigma to vertex TPC
      if(fITSRefitFlag) {
	if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(11)))->Fill(0);
	status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kITSrefit) != 0)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(11)))->Fill(0);
      }//ITS refit
      if(fTPCRefitFlag) {
	if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(12)))->Fill(0);
	  status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kTPCrefit) != 0)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(12)))->Fill(0);
      }//TPC refit
      if(fESDpidFlag) {
	if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(13)))->Fill(0);
	  status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kESDpid) != 0)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(13)))->Fill(0);
      }//ESD pid
      if(fTPCpidFlag) {
	if ((track->GetStatus() & AliESDtrack::kTPCpid) == 0) {
	  ((TH1F *)(fQAPrimaryAntiProtonsRejectedList->At(13)))->Fill(0);
	  status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kTPCpid) != 0)
	  ((TH1F *)(fQAPrimaryAntiProtonsAcceptedList->At(13)))->Fill(0);
      }//TPC pid
    }//primary particle cut

    //Secondaries
    if(label > nPrimaries) {
      if(fMinITSClustersFlag) {
	if(nClustersITS < fMinITSClusters) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(0)))->Fill(nClustersITS);
	  status = kFALSE;
	}
	else if(nClustersITS >= fMinITSClusters) 
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(0)))->Fill(nClustersITS);
      }//ITS clusters
      if(fMaxChi2PerITSClusterFlag) {
	if(chi2PerClusterITS > fMaxChi2PerITSCluster) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(1)))->Fill(chi2PerClusterITS);
	  status = kFALSE;
	}
	else if(chi2PerClusterITS <= fMaxChi2PerITSCluster)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(1)))->Fill(chi2PerClusterITS);
      }//chi2 per ITS cluster
      if(fMinTPCClustersFlag) {
	if(nClustersTPC < fMinTPCClusters) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(2)))->Fill(nClustersTPC);
	  status = kFALSE;
	}
	else if(nClustersTPC >= fMinTPCClusters)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(2)))->Fill(nClustersTPC);
      }//TPC clusters
      if(fMaxChi2PerTPCClusterFlag) {
	if(chi2PerClusterTPC > fMaxChi2PerTPCCluster) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(3)))->Fill(chi2PerClusterTPC);
	  status = kFALSE;
	}
	else if(chi2PerClusterTPC <= fMaxChi2PerTPCCluster)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(3)))->Fill(chi2PerClusterTPC);
      }//chi2 per TPC cluster
      if(fMaxCov11Flag) {
	if(extCov[0] > fMaxCov11) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(4)))->Fill(extCov[0]);
	  status = kFALSE;
	}
	else if(extCov[0] <= fMaxCov11)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(4)))->Fill(extCov[0]);
      }//cov11
      if(fMaxCov22Flag) {
	if(extCov[2] > fMaxCov22) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(5)))->Fill(extCov[2]);
	  status = kFALSE;
	}
	else if(extCov[2] <= fMaxCov22)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(5)))->Fill(extCov[2]);
      }//cov11
      if(fMaxCov33Flag) {
	if(extCov[5] > fMaxCov33) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(6)))->Fill(extCov[5]);
	  status = kFALSE;
	}
	else if(extCov[5] <= fMaxCov33)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(6)))->Fill(extCov[5]);
      }//cov11
      if(fMaxCov44Flag) {
	if(extCov[9] > fMaxCov44) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(7)))->Fill(extCov[9]);
	  status = kFALSE;
	}
	else if(extCov[9] <= fMaxCov44)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(7)))->Fill(extCov[9]);
      }//cov11
      if(fMaxCov55Flag) {
	if(extCov[14] > fMaxCov55) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(8)))->Fill(extCov[14]);
	  status = kFALSE;
	}
	else if(extCov[14] <= fMaxCov55)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(8)))->Fill(extCov[14]);
      }//cov55
      if(fMaxSigmaToVertexFlag) {
	if(GetSigmaToVertex(track) > fMaxSigmaToVertex) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(9)))->Fill(GetSigmaToVertex(track));
	  status = kFALSE;
	}
	else if(GetSigmaToVertex(track) <= fMaxSigmaToVertex)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(9)))->Fill(GetSigmaToVertex(track));
      }//sigma to vertex
      if(fMaxSigmaToVertexTPCFlag) {
	if(GetSigmaToVertex(track) > fMaxSigmaToVertexTPC) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(10)))->Fill(GetSigmaToVertex(track));
	  status = kFALSE;
	}
	else if(GetSigmaToVertex(track) <= fMaxSigmaToVertexTPC)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(10)))->Fill(GetSigmaToVertex(track));
      }//sigma to vertex TPC
      if(fITSRefitFlag) {
	if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(11)))->Fill(0);
	status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kITSrefit) != 0)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(11)))->Fill(0);
      }//ITS refit
      if(fTPCRefitFlag) {
	if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(12)))->Fill(0);
	  status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kTPCrefit) != 0)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(12)))->Fill(0);
      }//TPC refit
      if(fESDpidFlag) {
	if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(13)))->Fill(0);
	  status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kESDpid) != 0)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(13)))->Fill(0);
      }//ESD pid
      if(fTPCpidFlag) {
	if ((track->GetStatus() & AliESDtrack::kTPCpid) == 0) {
	  ((TH1F *)(fQASecondaryAntiProtonsRejectedList->At(13)))->Fill(0);
	  status = kFALSE;
	}
	else if((track->GetStatus() & AliESDtrack::kTPCpid) != 0)
	  ((TH1F *)(fQASecondaryAntiProtonsAcceptedList->At(13)))->Fill(0);
      }//TPC pid
    }//secondary particle cut
  }//antiprotons

  if((Pt < fMinPt) || (Pt > fMaxPt)) status = kFALSE;
  if((Rapidity(Px,Py,Pz) < fMinY) || (Rapidity(Px,Py,Pz) > fMaxY)) status = kFALSE;

  return status;
}

//____________________________________________________________________//
Float_t AliProtonAnalysis::GetSigmaToVertex(AliESDtrack* esdTrack) {
  // Calculates the number of sigma to the vertex.
  
  Float_t b[2];
  Float_t bRes[2];
  Float_t bCov[3];
  if(fUseTPCOnly) 
    esdTrack->GetImpactParametersTPC(b,bCov);
  else
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

//____________________________________________________________________//
Bool_t AliProtonAnalysis::ReadCorrectionContainer(const char* filename) {
  // Reads the outout of the correction framework task
  // Creates the correction maps
  // Puts the results in the different TList objects
  Bool_t status = kTRUE;

  TFile *file = TFile::Open(filename);
  if(!file) {
    cout<<"Could not find the input CORRFW file "<<filename<<endl;
    status = kFALSE;
  }

  AliCFContainer *corrfwContainer = (AliCFContainer*) (file->Get("container"));
  if(!corrfwContainer) {
    cout<<"CORRFW container not found!"<<endl;
    status = kFALSE;
  }
  
  Int_t nSteps = corrfwContainer->GetNStep();
  TH2D *gYPt[4];
  //currently the GRID is formed by the y-pT parameters
  //Add Vz as a next step
  Int_t iRap = 0, iPt = 1;
  for(Int_t iStep = 0; iStep < nSteps; iStep++) {
    gYPt[iStep] = corrfwContainer->ShowProjection(iRap,iPt,iStep);
    //fCorrectionList2D->Add(gYPt[iStep]);
  }

  //construct the efficiency grid from the data container                                           
  TString gTitle = 0;
  AliCFEffGrid *efficiency[3]; //The efficiency array has nStep-1 entries!!!
  TH1D *gEfficiency[2][3]; //efficiency as a function of pT and of y (raws - [2]) 
  TH1D *gCorrection[2][3]; //efficiency as a function of pT and of y (raws - [2]) 

  //Get the 2D efficiency maps
  for(Int_t iStep = 1; iStep < nSteps; iStep++) { 
    gTitle = "Efficiency_Step0_Step"; gTitle += iStep; 
    efficiency[iStep] = new AliCFEffGrid(gTitle.Data(),
					 gTitle.Data(),*corrfwContainer);
    efficiency[iStep]->CalculateEfficiency(iStep,0); //eff= step[i]/step0
    fCorrectionList2D->Add(efficiency[iStep]);  
  }
  //Get the projection of the efficiency maps
  for(Int_t iParameter = 0; iParameter < 2; iParameter++) { 
    for(Int_t iStep = 1; iStep < nSteps; iStep++) { 
      gEfficiency[iParameter][iStep-1] = efficiency[iStep]->Project(iParameter);
      fEfficiencyList1D->Add(gEfficiency[iParameter][iStep-1]);  
      gTitle = "Correction_Parameter"; gTitle += iParameter+1;
      gTitle += "_Step0_Step"; gTitle += iStep; 
      gCorrection[iParameter][iStep-1] = new TH1D(gTitle.Data(),
						   gTitle.Data(),
						   gEfficiency[iParameter][iStep-1]->GetNbinsX(),
						   gEfficiency[iParameter][iStep-1]->GetXaxis()->GetXmin(),
						   gEfficiency[iParameter][iStep-1]->GetXaxis()->GetXmax());
      //initialisation of the correction
      for(Int_t iBin = 1; iBin <= gEfficiency[iParameter][iStep-1]->GetNbinsX(); iBin++)
	gCorrection[iParameter][iStep-1]->SetBinContent(iBin,1.0);
    }//step loop
  }//parameter loop
  //Calculate the 1D correction parameters as a function of y and pT
  for(Int_t iParameter = 0; iParameter < 2; iParameter++) { 
    for(Int_t iStep = 1; iStep < nSteps; iStep++) { 
      gCorrection[iParameter][iStep-1]->Divide(gEfficiency[iParameter][iStep-1]);
      fCorrectionList1D->Add(gCorrection[iParameter][iStep-1]);  
    }
  }

  return status;
}

//____________________________________________________________________//
void AliProtonAnalysis::InitQA() {
  //Initializes the QA histograms and builds the directory structure
  //2D histograms
  TDirectory *dir2D = gDirectory->mkdir("2D");
  fGlobalQAList->Add(dir2D); dir2D->cd();
  TH2F *gHistYPtPrimaryProtonsPass = new TH2F("gHistYPtPrimaryProtonsPass",
					      ";y;P_{T} [GeV/c]",
					      fNBinsY,fMinY,fMaxY,
					      fNBinsPt,fMinPt,fMaxPt);
  gHistYPtPrimaryProtonsPass->SetStats(kTRUE);
  gHistYPtPrimaryProtonsPass->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtPrimaryProtonsPass);
  TH2F *gHistYPtPrimaryAntiProtonsPass = new TH2F("gHistYPtAntiPrimaryProtonsPass",
						  ";y;P_{T} [GeV/c]",
						  fNBinsY,fMinY,fMaxY,
						  fNBinsPt,fMinPt,fMaxPt);
  gHistYPtPrimaryAntiProtonsPass->SetStats(kTRUE);
  gHistYPtPrimaryAntiProtonsPass->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtPrimaryAntiProtonsPass);
  TH2F *gHistYPtSecondaryProtonsPass = new TH2F("gHistYPtSecondaryAntiProtonsPass",
						";y;P_{T} [GeV/c]",
						fNBinsY,fMinY,fMaxY,
						fNBinsPt,fMinPt,fMaxPt);
  gHistYPtSecondaryProtonsPass->SetStats(kTRUE);
  gHistYPtSecondaryProtonsPass->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtSecondaryProtonsPass);
  TH2F *gHistYPtSecondaryAntiAntiProtonsPass = new TH2F("gHistYPtAntiSecondaryAntiProtonsPass",
							";y;P_{T} [GeV/c]",
							fNBinsY,fMinY,fMaxY,
							fNBinsPt,fMinPt,fMaxPt);
  gHistYPtSecondaryAntiAntiProtonsPass->SetStats(kTRUE);
  gHistYPtSecondaryAntiAntiProtonsPass->GetXaxis()->SetTitleColor(1);
  fQA2DList->Add(gHistYPtSecondaryAntiAntiProtonsPass);
  
  gDirectory->cd("../");
  //protons
  TDirectory *dirProtons = gDirectory->mkdir("Protons");
  fGlobalQAList->Add(dirProtons); dirProtons->cd();
  
  //________________________________________________________________//
  TDirectory *dirProtonsPrimary = gDirectory->mkdir("Primaries");
  dirProtonsPrimary->cd();
  TDirectory *dirProtonsPrimaryAccepted = gDirectory->mkdir("Accepted");
  dirProtonsPrimaryAccepted->cd();

  //Accepted primary protons
  TH1F *fPrimaryProtonsITSClustersPass = new TH1F("fPrimaryProtonsITSClustersPass",
					    ";N_{clusters} (ITS);Entries",
					    7,0,7);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsITSClustersPass);
  TH1F *fPrimaryProtonsChi2PerClusterITSPass = new TH1F("fPrimaryProtonsChi2PerClusterITSPass",
						  ";x^{2}/N_{clusters} (ITS);Entries",
						  100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsChi2PerClusterITSPass);
  TH1F *fPrimaryProtonsTPCClustersPass = new TH1F("fPrimaryProtonsTPCClustersPass",
					    ";N_{clusters} (TPC);Entries",
					    100,0,200);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsTPCClustersPass);
  TH1F *fPrimaryProtonsChi2PerClusterTPCPass = new TH1F("fPrimaryProtonsChi2PerClusterTPCPass",
						  ";x^{2}/N_{clusters} (TPC);Entries",
						  100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsChi2PerClusterTPCPass);
  TH1F *fPrimaryProtonsExtCov11Pass = new TH1F("fPrimaryProtonsExtCov11Pass",
					 ";#sigma_{y} [cm];Entries",
					 100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsExtCov11Pass);
  TH1F *fPrimaryProtonsExtCov22Pass = new TH1F("fPrimaryProtonsExtCov22Pass",
					 ";#sigma_{z} [cm];Entries",
					 100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsExtCov22Pass);
  TH1F *fPrimaryProtonsExtCov33Pass = new TH1F("fPrimaryProtonsExtCov33Pass",
					 ";#sigma_{sin(#phi)};Entries",
					 100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsExtCov33Pass);
  TH1F *fPrimaryProtonsExtCov44Pass = new TH1F("fPrimaryProtonsExtCov44Pass",
					 ";#sigma_{tan(#lambda)};Entries",
					 100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsExtCov44Pass);
  TH1F *fPrimaryProtonsExtCov55Pass = new TH1F("fPrimaryProtonsExtCov55Pass",
					 ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
					 100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsExtCov55Pass);
  TH1F *fPrimaryProtonsSigmaToVertexPass = new TH1F("fPrimaryProtonsSigmaToVertexPass",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsSigmaToVertexPass);
  TH1F *fPrimaryProtonsSigmaToVertexTPCPass = new TH1F("fPrimaryProtonsSigmaToVertexTPCPass",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsSigmaToVertexTPCPass);
  TH1F *fPrimaryProtonsITSRefitPass = new TH1F("fPrimaryProtonsITSRefitPass",
					 "",10,-1,1);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsITSRefitPass);
  TH1F *fPrimaryProtonsTPCRefitPass = new TH1F("fPrimaryProtonsTPCRefitPass",
					 "",10,-1,1);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsTPCRefitPass);
  TH1F *fPrimaryProtonsESDpidPass = new TH1F("fPrimaryProtonsESDpidPass",
				       "",10,-1,1);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsESDpidPass);
  TH1F *fPrimaryProtonsTPCpidPass = new TH1F("fPrimaryProtonsTPCpidPass",
				       "",10,-1,1);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsTPCpidPass);

  //Rejected primary protons
  gDirectory->cd("../");
  TDirectory *dirProtonsPrimaryRejected = gDirectory->mkdir("Rejected");
  dirProtonsPrimaryRejected->cd();

  TH1F *fPrimaryProtonsITSClustersReject = new TH1F("fPrimaryProtonsITSClustersReject",
						    ";N_{clusters} (ITS);Entries",
						    7,0,7);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsITSClustersReject);
  TH1F *fPrimaryProtonsChi2PerClusterITSReject = new TH1F("fPrimaryProtonsChi2PerClusterITSReject",
							  ";x^{2}/N_{clusters} (ITS);Entries",
							  100,0,4);
  fQAPrimaryProtonsAcceptedList->Add(fPrimaryProtonsChi2PerClusterITSReject);
  TH1F *fPrimaryProtonsTPCClustersReject = new TH1F("fPrimaryProtonsTPCClustersReject",
					    ";N_{clusters} (TPC);Entries",
					    100,0,200);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsTPCClustersReject);
  TH1F *fPrimaryProtonsChi2PerClusterTPCReject = new TH1F("fPrimaryProtonsChi2PerClusterTPCReject",
						  ";x^{2}/N_{clusters} (TPC);Entries",
						  100,0,4);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsChi2PerClusterTPCReject);
  TH1F *fPrimaryProtonsExtCov11Reject = new TH1F("fPrimaryProtonsExtCov11Reject",
					 ";#sigma_{y} [cm];Entries",
					 100,0,4);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsExtCov11Reject);
  TH1F *fPrimaryProtonsExtCov22Reject = new TH1F("fPrimaryProtonsExtCov22Reject",
					 ";#sigma_{z} [cm];Entries",
					 100,0,4);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsExtCov22Reject);
  TH1F *fPrimaryProtonsExtCov33Reject = new TH1F("fPrimaryProtonsExtCov33Reject",
					 ";#sigma_{sin(#phi)};Entries",
					 100,0,4);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsExtCov33Reject);
  TH1F *fPrimaryProtonsExtCov44Reject = new TH1F("fPrimaryProtonsExtCov44Reject",
					 ";#sigma_{tan(#lambda)};Entries",
					 100,0,4);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsExtCov44Reject);
  TH1F *fPrimaryProtonsExtCov55Reject = new TH1F("fPrimaryProtonsExtCov55Reject",
					 ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
					 100,0,4);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsExtCov55Reject);
  TH1F *fPrimaryProtonsSigmaToVertexReject = new TH1F("fPrimaryProtonsSigmaToVertexReject",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsSigmaToVertexReject);
  TH1F *fPrimaryProtonsSigmaToVertexTPCReject = new TH1F("fPrimaryProtonsSigmaToVertexTPCReject",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsSigmaToVertexTPCReject);
  TH1F *fPrimaryProtonsITSRefitReject = new TH1F("fPrimaryProtonsITSRefitReject",
					 "",10,-1,1);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsITSRefitReject);
  TH1F *fPrimaryProtonsTPCRefitReject = new TH1F("fPrimaryProtonsTPCRefitReject",
					 "",10,-1,1);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsTPCRefitReject);
  TH1F *fPrimaryProtonsESDpidReject = new TH1F("fPrimaryProtonsESDpidReject",
				       "",10,-1,1);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsESDpidReject);
  TH1F *fPrimaryProtonsTPCpidReject = new TH1F("fPrimaryProtonsTPCpidReject",
				       "",10,-1,1);
  fQAPrimaryProtonsRejectedList->Add(fPrimaryProtonsTPCpidReject);

  //________________________________________________________________//
  gDirectory->cd("../../");

  TDirectory *dirProtonsSecondary = gDirectory->mkdir("Secondaries");
  dirProtonsSecondary->cd();
  TDirectory *dirProtonsSecondaryAccepted = gDirectory->mkdir("Accepted");
  dirProtonsSecondaryAccepted->cd();

  //Accepted secondary protons
  TH1F *fSecondaryProtonsITSClustersPass = new TH1F("fSecondaryProtonsITSClustersPass",
						    ";N_{clusters} (ITS);Entries",
						    7,0,7);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsITSClustersPass);
  TH1F *fSecondaryProtonsChi2PerClusterITSPass = new TH1F("fSecondaryProtonsChi2PerClusterITSPass",
							  ";x^{2}/N_{clusters} (ITS);Entries",
							  100,0,4);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsChi2PerClusterITSPass);
  TH1F *fSecondaryProtonsTPCClustersPass = new TH1F("fSecondaryProtonsTPCClustersPass",
					    ";N_{clusters} (TPC);Entries",
					    100,0,200);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsTPCClustersPass);
  TH1F *fSecondaryProtonsChi2PerClusterTPCPass = new TH1F("fSecondaryProtonsChi2PerClusterTPCPass",
						  ";x^{2}/N_{clusters} (TPC);Entries",
						  100,0,4);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsChi2PerClusterTPCPass);
  TH1F *fSecondaryProtonsExtCov11Pass = new TH1F("fSecondaryProtonsExtCov11Pass",
					 ";#sigma_{y} [cm];Entries",
					 100,0,4);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsExtCov11Pass);
  TH1F *fSecondaryProtonsExtCov22Pass = new TH1F("fSecondaryProtonsExtCov22Pass",
					 ";#sigma_{z} [cm];Entries",
					 100,0,4);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsExtCov22Pass);
  TH1F *fSecondaryProtonsExtCov33Pass = new TH1F("fSecondaryProtonsExtCov33Pass",
					 ";#sigma_{sin(#phi)};Entries",
					 100,0,4);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsExtCov33Pass);
  TH1F *fSecondaryProtonsExtCov44Pass = new TH1F("fSecondaryProtonsExtCov44Pass",
					 ";#sigma_{tan(#lambda)};Entries",
					 100,0,4);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsExtCov44Pass);
  TH1F *fSecondaryProtonsExtCov55Pass = new TH1F("fSecondaryProtonsExtCov55Pass",
					 ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
					 100,0,4);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsExtCov55Pass);
  TH1F *fSecondaryProtonsSigmaToVertexPass = new TH1F("fSecondaryProtonsSigmaToVertexPass",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsSigmaToVertexPass);
  TH1F *fSecondaryProtonsSigmaToVertexTPCPass = new TH1F("fSecondaryProtonsSigmaToVertexTPCPass",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsSigmaToVertexTPCPass);
  TH1F *fSecondaryProtonsITSRefitPass = new TH1F("fSecondaryProtonsITSRefitPass",
					 "",10,-1,1);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsITSRefitPass);
  TH1F *fSecondaryProtonsTPCRefitPass = new TH1F("fSecondaryProtonsTPCRefitPass",
					 "",10,-1,1);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsTPCRefitPass);
  TH1F *fSecondaryProtonsESDpidPass = new TH1F("fSecondaryProtonsESDpidPass",
				       "",10,-1,1);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsESDpidPass);
  TH1F *fSecondaryProtonsTPCpidPass = new TH1F("fSecondaryProtonsTPCpidPass",
				       "",10,-1,1);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsTPCpidPass);

  //Rejected secondary protons
  gDirectory->cd("../");
  TDirectory *dirProtonsSecondaryRejected = gDirectory->mkdir("Rejected");
  dirProtonsSecondaryRejected->cd();

  TH1F *fSecondaryProtonsITSClustersReject = new TH1F("fSecondaryProtonsITSClustersReject",
						      ";N_{clusters} (ITS);Entries",
						      7,0,7);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsITSClustersReject);
  TH1F *fSecondaryProtonsChi2PerClusterITSReject = new TH1F("fSecondaryProtonsChi2PerClusterITSReject",
							    ";x^{2}/N_{clusters} (ITS);Entries",
							    100,0,4);
  fQASecondaryProtonsAcceptedList->Add(fSecondaryProtonsChi2PerClusterITSReject);
  TH1F *fSecondaryProtonsTPCClustersReject = new TH1F("fSecondaryProtonsTPCClustersReject",
					    ";N_{clusters} (TPC);Entries",
					    100,0,200);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsTPCClustersReject);
  TH1F *fSecondaryProtonsChi2PerClusterTPCReject = new TH1F("fSecondaryProtonsChi2PerClusterTPCReject",
						  ";x^{2}/N_{clusters} (TPC);Entries",
						  100,0,4);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsChi2PerClusterTPCReject);
  TH1F *fSecondaryProtonsExtCov11Reject = new TH1F("fSecondaryProtonsExtCov11Reject",
					 ";#sigma_{y} [cm];Entries",
					 100,0,4);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsExtCov11Reject);
  TH1F *fSecondaryProtonsExtCov22Reject = new TH1F("fSecondaryProtonsExtCov22Reject",
					 ";#sigma_{z} [cm];Entries",
					 100,0,4);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsExtCov22Reject);
  TH1F *fSecondaryProtonsExtCov33Reject = new TH1F("fSecondaryProtonsExtCov33Reject",
					 ";#sigma_{sin(#phi)};Entries",
					 100,0,4);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsExtCov33Reject);
  TH1F *fSecondaryProtonsExtCov44Reject = new TH1F("fSecondaryProtonsExtCov44Reject",
					 ";#sigma_{tan(#lambda)};Entries",
					 100,0,4);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsExtCov44Reject);
  TH1F *fSecondaryProtonsExtCov55Reject = new TH1F("fSecondaryProtonsExtCov55Reject",
					 ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
					 100,0,4);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsExtCov55Reject);
  TH1F *fSecondaryProtonsSigmaToVertexReject = new TH1F("fSecondaryProtonsSigmaToVertexReject",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsSigmaToVertexReject);
  TH1F *fSecondaryProtonsSigmaToVertexTPCReject = new TH1F("fSecondaryProtonsSigmaToVertexTPCReject",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsSigmaToVertexTPCReject);
  TH1F *fSecondaryProtonsITSRefitReject = new TH1F("fSecondaryProtonsITSRefitReject",
					 "",10,-1,1);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsITSRefitReject);
  TH1F *fSecondaryProtonsTPCRefitReject = new TH1F("fSecondaryProtonsTPCRefitReject",
					 "",10,-1,1);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsTPCRefitReject);
  TH1F *fSecondaryProtonsESDpidReject = new TH1F("fSecondaryProtonsESDpidReject",
				       "",10,-1,1);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsESDpidReject);
  TH1F *fSecondaryProtonsTPCpidReject = new TH1F("fSecondaryProtonsTPCpidReject",
				       "",10,-1,1);
  fQASecondaryProtonsRejectedList->Add(fSecondaryProtonsTPCpidReject);


  gDirectory->cd("../../../");

  //antiprotons
  TDirectory *dirAntiProtons = gDirectory->mkdir("AntiProtons");
  fGlobalQAList->Add(dirAntiProtons); dirAntiProtons->cd();
  
  //________________________________________________________________//
  TDirectory *dirAntiProtonsPrimary = gDirectory->mkdir("Primaries");
  dirAntiProtonsPrimary->cd();
  TDirectory *dirAntiProtonsPrimaryAccepted = gDirectory->mkdir("Accepted");
  dirAntiProtonsPrimaryAccepted->cd();
  
  //Accepted primary antiprotons
  TH1F *fPrimaryAntiProtonsITSClustersPass = new TH1F("fPrimaryAntiProtonsITSClustersPass",
						      ";N_{clusters} (ITS);Entries",
						      7,0,7);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsITSClustersPass);
  TH1F *fPrimaryAntiProtonsChi2PerClusterITSPass = new TH1F("fPrimaryAntiProtonsChi2PerClusterITSPass",
							    ";x^{2}/N_{clusters} (ITS);Entries",
							    100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsChi2PerClusterITSPass);
  TH1F *fPrimaryAntiProtonsTPCClustersPass = new TH1F("fPrimaryAntiProtonsTPCClustersPass",
						      ";N_{clusters} (TPC);Entries",
						      100,0,200);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsTPCClustersPass);
  TH1F *fPrimaryAntiProtonsChi2PerClusterTPCPass = new TH1F("fPrimaryAntiProtonsChi2PerClusterTPCPass",
							    ";x^{2}/N_{clusters} (TPC);Entries",
							    100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsChi2PerClusterTPCPass);
  TH1F *fPrimaryAntiProtonsExtCov11Pass = new TH1F("fPrimaryAntiProtonsExtCov11Pass",
						   ";#sigma_{y} [cm];Entries",
						   100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsExtCov11Pass);
  TH1F *fPrimaryAntiProtonsExtCov22Pass = new TH1F("fPrimaryAntiProtonsExtCov22Pass",
						   ";#sigma_{z} [cm];Entries",
						   100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsExtCov22Pass);
  TH1F *fPrimaryAntiProtonsExtCov33Pass = new TH1F("fPrimaryAntiProtonsExtCov33Pass",
						   ";#sigma_{sin(#phi)};Entries",
						   100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsExtCov33Pass);
  TH1F *fPrimaryAntiProtonsExtCov44Pass = new TH1F("fPrimaryAntiProtonsExtCov44Pass",
						   ";#sigma_{tan(#lambda)};Entries",
						   100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsExtCov44Pass);
  TH1F *fPrimaryAntiProtonsExtCov55Pass = new TH1F("fPrimaryAntiProtonsExtCov55Pass",
						   ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
						   100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsExtCov55Pass);
  TH1F *fPrimaryAntiProtonsSigmaToVertexPass = new TH1F("fPrimaryAntiProtonsSigmaToVertexPass",
							";#sigma_{Vertex};Entries",
							100,0,10);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsSigmaToVertexPass);
  TH1F *fPrimaryAntiProtonsSigmaToVertexTPCPass = new TH1F("fPrimaryAntiProtonsSigmaToVertexTPCPass",
							   ";#sigma_{Vertex};Entries",
							   100,0,10);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsSigmaToVertexTPCPass);
  TH1F *fPrimaryAntiProtonsITSRefitPass = new TH1F("fPrimaryAntiProtonsITSRefitPass",
						   "",10,-1,1);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsITSRefitPass);
  TH1F *fPrimaryAntiProtonsTPCRefitPass = new TH1F("fPrimaryAntiProtonsTPCRefitPass",
						   "",10,-1,1);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsTPCRefitPass);
  TH1F *fPrimaryAntiProtonsESDpidPass = new TH1F("fPrimaryAntiProtonsESDpidPass",
						 "",10,-1,1);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsESDpidPass);
  TH1F *fPrimaryAntiProtonsTPCpidPass = new TH1F("fPrimaryAntiProtonsTPCpidPass",
						 "",10,-1,1);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsTPCpidPass);
  
  //Rejected primary antiprotons
  gDirectory->cd("../");
  TDirectory *dirAntiProtonsPrimaryRejected = gDirectory->mkdir("Rejected");
  dirAntiProtonsPrimaryRejected->cd();
  
  TH1F *fPrimaryAntiProtonsITSClustersReject = new TH1F("fPrimaryAntiProtonsITSClustersReject",
							";N_{clusters} (ITS);Entries",
							7,0,7);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsITSClustersReject);
  TH1F *fPrimaryAntiProtonsChi2PerClusterITSReject = new TH1F("fPrimaryAntiProtonsChi2PerClusterITSReject",
							      ";x^{2}/N_{clusters} (ITS);Entries",
							      100,0,4);
  fQAPrimaryAntiProtonsAcceptedList->Add(fPrimaryAntiProtonsChi2PerClusterITSReject);
  TH1F *fPrimaryAntiProtonsTPCClustersReject = new TH1F("fPrimaryAntiProtonsTPCClustersReject",
							";N_{clusters} (TPC);Entries",
							100,0,200);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsTPCClustersReject);
  TH1F *fPrimaryAntiProtonsChi2PerClusterTPCReject = new TH1F("fPrimaryAntiProtonsChi2PerClusterTPCReject",
							      ";x^{2}/N_{clusters} (TPC);Entries",
							      100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsChi2PerClusterTPCReject);
  TH1F *fPrimaryAntiProtonsExtCov11Reject = new TH1F("fPrimaryAntiProtonsExtCov11Reject",
						     ";#sigma_{y} [cm];Entries",
						     100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsExtCov11Reject);
  TH1F *fPrimaryAntiProtonsExtCov22Reject = new TH1F("fPrimaryAntiProtonsExtCov22Reject",
						     ";#sigma_{z} [cm];Entries",
						     100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsExtCov22Reject);
  TH1F *fPrimaryAntiProtonsExtCov33Reject = new TH1F("fPrimaryAntiProtonsExtCov33Reject",
						     ";#sigma_{sin(#phi)};Entries",
						     100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsExtCov33Reject);
  TH1F *fPrimaryAntiProtonsExtCov44Reject = new TH1F("fPrimaryAntiProtonsExtCov44Reject",
						     ";#sigma_{tan(#lambda)};Entries",
						     100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsExtCov44Reject);
  TH1F *fPrimaryAntiProtonsExtCov55Reject = new TH1F("fPrimaryAntiProtonsExtCov55Reject",
						     ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
						     100,0,4);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsExtCov55Reject);
  TH1F *fPrimaryAntiProtonsSigmaToVertexReject = new TH1F("fPrimaryAntiProtonsSigmaToVertexReject",
							  ";#sigma_{Vertex};Entries",
							  100,0,10);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsSigmaToVertexReject);
  TH1F *fPrimaryAntiProtonsSigmaToVertexTPCReject = new TH1F("fPrimaryAntiProtonsSigmaToVertexTPCReject",
							     ";#sigma_{Vertex};Entries",
							     100,0,10);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsSigmaToVertexTPCReject);
  TH1F *fPrimaryAntiProtonsITSRefitReject = new TH1F("fPrimaryAntiProtonsITSRefitReject",
						     "",10,-1,1);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsITSRefitReject);
  TH1F *fPrimaryAntiProtonsTPCRefitReject = new TH1F("fPrimaryAntiProtonsTPCRefitReject",
						     "",10,-1,1);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsTPCRefitReject);
  TH1F *fPrimaryAntiProtonsESDpidReject = new TH1F("fPrimaryAntiProtonsESDpidReject",
						   "",10,-1,1);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsESDpidReject);
  TH1F *fPrimaryAntiProtonsTPCpidReject = new TH1F("fPrimaryAntiProtonsTPCpidReject",
						   "",10,-1,1);
  fQAPrimaryAntiProtonsRejectedList->Add(fPrimaryAntiProtonsTPCpidReject);
  
  //________________________________________________________________//
  gDirectory->cd("../../");

  TDirectory *dirAntiProtonsSecondary = gDirectory->mkdir("Secondaries");
  dirAntiProtonsSecondary->cd();
  TDirectory *dirAntiProtonsSecondaryAccepted = gDirectory->mkdir("Accepted");
  dirAntiProtonsSecondaryAccepted->cd();

  //Accepted secondary antiprotons
  TH1F *fSecondaryAntiProtonsITSClustersPass = new TH1F("fSecondaryAntiProtonsITSClustersPass",
							";N_{clusters} (ITS);Entries",
							7,0,7);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsITSClustersPass);
  TH1F *fSecondaryAntiProtonsChi2PerClusterITSPass = new TH1F("fSecondaryAntiProtonsChi2PerClusterITSPass",
							      ";x^{2}/N_{clusters} (ITS);Entries",
							      100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsChi2PerClusterITSPass);
  TH1F *fSecondaryAntiProtonsTPCClustersPass = new TH1F("fSecondaryAntiProtonsTPCClustersPass",
							";N_{clusters} (TPC);Entries",
							100,0,200);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsTPCClustersPass);
  TH1F *fSecondaryAntiProtonsChi2PerClusterTPCPass = new TH1F("fSecondaryAntiProtonsChi2PerClusterTPCPass",
							      ";x^{2}/N_{clusters} (TPC);Entries",
							      100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsChi2PerClusterTPCPass);
  TH1F *fSecondaryAntiProtonsExtCov11Pass = new TH1F("fSecondaryAntiProtonsExtCov11Pass",
						     ";#sigma_{y} [cm];Entries",
						     100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsExtCov11Pass);
  TH1F *fSecondaryAntiProtonsExtCov22Pass = new TH1F("fSecondaryAntiProtonsExtCov22Pass",
						     ";#sigma_{z} [cm];Entries",
						     100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsExtCov22Pass);
  TH1F *fSecondaryAntiProtonsExtCov33Pass = new TH1F("fSecondaryAntiProtonsExtCov33Pass",
						     ";#sigma_{sin(#phi)};Entries",
						     100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsExtCov33Pass);
  TH1F *fSecondaryAntiProtonsExtCov44Pass = new TH1F("fSecondaryAntiProtonsExtCov44Pass",
						     ";#sigma_{tan(#lambda)};Entries",
						     100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsExtCov44Pass);
  TH1F *fSecondaryAntiProtonsExtCov55Pass = new TH1F("fSecondaryAntiProtonsExtCov55Pass",
						     ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
						     100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsExtCov55Pass);
  TH1F *fSecondaryAntiProtonsSigmaToVertexPass = new TH1F("fSecondaryAntiProtonsSigmaToVertexPass",
							  ";#sigma_{Vertex};Entries",
							  100,0,10);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsSigmaToVertexPass);
  TH1F *fSecondaryAntiProtonsSigmaToVertexTPCPass = new TH1F("fSecondaryAntiProtonsSigmaToVertexTPCPass",
							     ";#sigma_{Vertex};Entries",
							     100,0,10);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsSigmaToVertexTPCPass);
  TH1F *fSecondaryAntiProtonsITSRefitPass = new TH1F("fSecondaryAntiProtonsITSRefitPass",
						     "",10,-1,1);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsITSRefitPass);
  TH1F *fSecondaryAntiProtonsTPCRefitPass = new TH1F("fSecondaryAntiProtonsTPCRefitPass",
						     "",10,-1,1);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsTPCRefitPass);
  TH1F *fSecondaryAntiProtonsESDpidPass = new TH1F("fSecondaryAntiProtonsESDpidPass",
						   "",10,-1,1);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsESDpidPass);
  TH1F *fSecondaryAntiProtonsTPCpidPass = new TH1F("fSecondaryAntiProtonsTPCpidPass",
						   "",10,-1,1);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsTPCpidPass);
  
  //Rejected secondary antiprotons
  gDirectory->cd("../");
  TDirectory *dirAntiProtonsSecondaryRejected = gDirectory->mkdir("Rejected");
  dirAntiProtonsSecondaryRejected->cd();

  TH1F *fSecondaryAntiProtonsITSClustersReject = new TH1F("fSecondaryAntiProtonsITSClustersReject",
							  ";N_{clusters} (ITS);Entries",
							  7,0,7);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsITSClustersReject);
  TH1F *fSecondaryAntiProtonsChi2PerClusterITSReject = new TH1F("fSecondaryAntiProtonsChi2PerClusterITSReject",
								";x^{2}/N_{clusters} (ITS);Entries",
								100,0,4);
  fQASecondaryAntiProtonsAcceptedList->Add(fSecondaryAntiProtonsChi2PerClusterITSReject);
  TH1F *fSecondaryAntiProtonsTPCClustersReject = new TH1F("fSecondaryAntiProtonsTPCClustersReject",
							  ";N_{clusters} (TPC);Entries",
							  100,0,200);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsTPCClustersReject);
  TH1F *fSecondaryAntiProtonsChi2PerClusterTPCReject = new TH1F("fSecondaryAntiProtonsChi2PerClusterTPCReject",
								";x^{2}/N_{clusters} (TPC);Entries",
								100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsChi2PerClusterTPCReject);
  TH1F *fSecondaryAntiProtonsExtCov11Reject = new TH1F("fSecondaryAntiProtonsExtCov11Reject",
						       ";#sigma_{y} [cm];Entries",
						       100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsExtCov11Reject);
  TH1F *fSecondaryAntiProtonsExtCov22Reject = new TH1F("fSecondaryAntiProtonsExtCov22Reject",
						       ";#sigma_{z} [cm];Entries",
						       100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsExtCov22Reject);
  TH1F *fSecondaryAntiProtonsExtCov33Reject = new TH1F("fSecondaryAntiProtonsExtCov33Reject",
						       ";#sigma_{sin(#phi)};Entries",
						       100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsExtCov33Reject);
  TH1F *fSecondaryAntiProtonsExtCov44Reject = new TH1F("fSecondaryAntiProtonsExtCov44Reject",
						       ";#sigma_{tan(#lambda)};Entries",
						       100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsExtCov44Reject);
  TH1F *fSecondaryAntiProtonsExtCov55Reject = new TH1F("fSecondaryAntiProtonsExtCov55Reject",
						       ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
						       100,0,4);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsExtCov55Reject);
  TH1F *fSecondaryAntiProtonsSigmaToVertexReject = new TH1F("fSecondaryAntiProtonsSigmaToVertexReject",
							    ";#sigma_{Vertex};Entries",
							    100,0,10);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsSigmaToVertexReject);
  TH1F *fSecondaryAntiProtonsSigmaToVertexTPCReject = new TH1F("fSecondaryAntiProtonsSigmaToVertexTPCReject",
							       ";#sigma_{Vertex};Entries",
							       100,0,10);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsSigmaToVertexTPCReject);
  TH1F *fSecondaryAntiProtonsITSRefitReject = new TH1F("fSecondaryAntiProtonsITSRefitReject",
						       "",10,-1,1);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsITSRefitReject);
  TH1F *fSecondaryAntiProtonsTPCRefitReject = new TH1F("fSecondaryAntiProtonsTPCRefitReject",
						       "",10,-1,1);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsTPCRefitReject);
  TH1F *fSecondaryAntiProtonsESDpidReject = new TH1F("fSecondaryAntiProtonsESDpidReject",
						     "",10,-1,1);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsESDpidReject);
  TH1F *fSecondaryAntiProtonsTPCpidReject = new TH1F("fSecondaryAntiProtonsTPCpidReject",
						     "",10,-1,1);
  fQASecondaryAntiProtonsRejectedList->Add(fSecondaryAntiProtonsTPCpidReject);
  
}

//____________________________________________________________________//
void AliProtonAnalysis::RunQA(AliStack *stack, AliESDEvent *fESD) {
  //Runs the QA code
  Int_t nGoodTracks = fESD->GetNumberOfTracks();
  for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    Int_t label = TMath::Abs(track->GetLabel()); 
    Double_t Pt = 0.0, P = 0.0;
    Double_t probability[5];

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
	if(IsAccepted(track, stack)) {
	  if(label <= stack->GetNprimary()) {
            if(track->Charge() > 0)
              ((TH2F *)(fQA2DList->At(0)))->Fill(Rapidity(track->Px(),track->Py(),track->Pz()),Pt);
            else if(track->Charge() < 0)
              ((TH2F *)(fQA2DList->At(1)))->Fill(Rapidity(track->Px(),track->Py(),track->Pz()),Pt);
          }//primary particles
          else if(label > stack->GetNprimary()) {
            if(track->Charge() > 0)
              ((TH2F *)(fQA2DList->At(2)))->Fill(Rapidity(track->Px(),track->Py(),track->Pz()),Pt);
            else if(track->Charge() < 0)
              ((TH2F *)(fQA2DList->At(3)))->Fill(Rapidity(track->Px(),track->Py(),track->Pz()),Pt);
          }//secondary particles
	}//cuts
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
	if(IsAccepted(track, stack)) {
	  if(label <= stack->GetNprimary()) {
	    if(track->Charge() > 0)
	      ((TH2F *)(fQA2DList->At(0)))->Fill(Rapidity(track->Px(),track->Py(),track->Pz()),Pt);
	    else if(track->Charge() < 0)
	      ((TH2F *)(fQA2DList->At(1)))->Fill(Rapidity(track->Px(),track->Py(),track->Pz()),Pt);
	  }//primary particles
	  else if(label > stack->GetNprimary()) {
	    if(track->Charge() > 0)
	      ((TH2F *)(fQA2DList->At(2)))->Fill(Rapidity(track->Px(),track->Py(),track->Pz()),Pt);
	    else if(track->Charge() < 0)
	      ((TH2F *)(fQA2DList->At(3)))->Fill(Rapidity(track->Px(),track->Py(),track->Pz()),Pt);
	  }//secondary particles
	}//cuts
      }//proton check
    }//combined tracking
  }//track loop
    
}









