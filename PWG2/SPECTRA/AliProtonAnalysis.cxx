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
//   Origin: Panos Christakoglou | Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------
#include <Riostream.h>
#include <TFile.h>
#include <TSystem.h>
#include <TF1.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TParticle.h>
#include <TList.h>

#include <AliExternalTrackParam.h>
#include <AliAODEvent.h>
#include <AliESDEvent.h>
//#include <AliLog.h>
#include <AliPID.h>
#include <AliStack.h>
#include <AliCFContainer.h>
#include <AliCFEffGrid.h>
#include <AliCFDataGrid.h>
//#include <AliESDVertex.h>
class AliLog;
class AliESDVertex;

#include "AliProtonAnalysis.h"
#include "AliProtonAnalysisBase.h"

ClassImp(AliProtonAnalysis)

//____________________________________________________________________//
AliProtonAnalysis::AliProtonAnalysis() : 
  TObject(), fProtonAnalysisBase(0),
  fNBinsY(0), fMinY(0), fMaxY(0),
  fNBinsPt(0), fMinPt(0), fMaxPt(0),
  fProtonContainer(0), fAntiProtonContainer(0),
  fHistEvents(0), fHistYPtProtons(0), fHistYPtAntiProtons(0),
  fEffGridListProtons(0), fCorrectionListProtons2D(0), 
  fEfficiencyListProtons1D(0), fCorrectionListProtons1D(0),
  fEffGridListAntiProtons(0), fCorrectionListAntiProtons2D(0), 
  fEfficiencyListAntiProtons1D(0), fCorrectionListAntiProtons1D(0),
  fCorrectProtons(0), fCorrectAntiProtons(0),
  fGlobalQAList(0), fQA2DList(0),
  fQAProtonsAcceptedList(0), fQAProtonsRejectedList(0),
  fQAAntiProtonsAcceptedList(0), fQAAntiProtonsRejectedList(0),
  fInitQAFlag(kFALSE) {
  //Default constructor
}

//____________________________________________________________________//
AliProtonAnalysis::AliProtonAnalysis(Int_t nbinsY, 
				     Float_t fLowY, Float_t fHighY,
				     Int_t nbinsPt, 
				     Float_t fLowPt, Float_t fHighPt) : 
  TObject(), fProtonAnalysisBase(0),
  fNBinsY(nbinsY), fMinY(fLowY), fMaxY(fHighY),
  fNBinsPt(nbinsPt), fMinPt(fLowPt), fMaxPt(fHighPt),
  fProtonContainer(0), fAntiProtonContainer(0),
  fHistEvents(0), fHistYPtProtons(0), fHistYPtAntiProtons(0),
  fEffGridListProtons(0), fCorrectionListProtons2D(0), 
  fEfficiencyListProtons1D(0), fCorrectionListProtons1D(0),
  fEffGridListAntiProtons(0), fCorrectionListAntiProtons2D(0), 
  fEfficiencyListAntiProtons1D(0), fCorrectionListAntiProtons1D(0),
  fCorrectProtons(0), fCorrectAntiProtons(0),
  fGlobalQAList(0), fQA2DList(0),
  fQAProtonsAcceptedList(0), fQAProtonsRejectedList(0),
  fQAAntiProtonsAcceptedList(0), fQAAntiProtonsRejectedList(0),
  fInitQAFlag(kFALSE) {
  //Default constructor
  if(!fInitQAFlag) InitQA();
  fHistEvents = new TH1I("fHistEvents","Analyzed events",1,0,1);

  fHistYPtProtons = new TH2D("fHistYPtProtons","Protons",
			     fNBinsY,fMinY,fMaxY,
			     fNBinsPt,fMinPt,fMaxPt);
  fHistYPtProtons->SetStats(kTRUE);
  fHistYPtProtons->GetYaxis()->SetTitle("P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fHistYPtProtons->GetXaxis()->SetTitle("#eta");
  else
    fHistYPtProtons->GetXaxis()->SetTitle("y");
  fHistYPtProtons->GetXaxis()->SetTitleColor(1);

  fHistYPtAntiProtons = new TH2D("fHistYPtAntiProtons","Antiprotons",
				 fNBinsY,fMinY,fMaxY,
				 fNBinsPt,fMinPt,fMaxPt);
  fHistYPtAntiProtons->SetStats(kTRUE);
  fHistYPtAntiProtons->GetYaxis()->SetTitle("P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fHistYPtAntiProtons->GetXaxis()->SetTitle("#eta");
  else
    fHistYPtAntiProtons->GetXaxis()->SetTitle("y");
  fHistYPtAntiProtons->GetXaxis()->SetTitleColor(1);

  //setting up the containers
  Int_t iBin[2];
  iBin[0] = nbinsY;
  iBin[1] = nbinsPt;
  Double_t *binLimY = new Double_t[nbinsY+1];
  Double_t *binLimPt = new Double_t[nbinsPt+1];
  //values for bin lower bounds
  for(Int_t i = 0; i <= nbinsY; i++) 
    binLimY[i]=(Double_t)fLowY  + (fHighY - fLowY)  /nbinsY*(Double_t)i;
  for(Int_t i = 0; i <= nbinsPt; i++) 
    binLimPt[i]=(Double_t)fLowPt  + (fHighPt - fLowPt)  /nbinsPt*(Double_t)i;

  fProtonContainer = new AliCFContainer("containerProtons",
					"container for protons",
					kNSteps,2,iBin);
  fProtonContainer->SetBinLimits(0,binLimY); //rapidity or eta
  fProtonContainer->SetBinLimits(1,binLimPt); //pT
  fAntiProtonContainer = new AliCFContainer("containerAntiProtons",
					    "container for antiprotons",
					    kNSteps,2,iBin);
  fAntiProtonContainer->SetBinLimits(0,binLimY); //rapidity or eta
  fAntiProtonContainer->SetBinLimits(1,binLimPt); //pT
} 

//____________________________________________________________________//
AliProtonAnalysis::AliProtonAnalysis(Int_t nbinsY, Double_t *gY,
				     Int_t nbinsPt,Double_t *gPt) : 
  TObject(), fProtonAnalysisBase(0),
  fNBinsY(nbinsY), fMinY(gY[0]), fMaxY(gY[nbinsY]),
  fNBinsPt(nbinsPt), fMinPt(gPt[0]), fMaxPt(gPt[nbinsPt]),
  fProtonContainer(0), fAntiProtonContainer(0),
  fHistEvents(0), fHistYPtProtons(0), fHistYPtAntiProtons(0),
  fEffGridListProtons(0), fCorrectionListProtons2D(0), 
  fEfficiencyListProtons1D(0), fCorrectionListProtons1D(0),
  fEffGridListAntiProtons(0), fCorrectionListAntiProtons2D(0), 
  fEfficiencyListAntiProtons1D(0), fCorrectionListAntiProtons1D(0),
  fCorrectProtons(0), fCorrectAntiProtons(0),
  fGlobalQAList(0), fQA2DList(0),
  fQAProtonsAcceptedList(0), fQAProtonsRejectedList(0),
  fQAAntiProtonsAcceptedList(0), fQAAntiProtonsRejectedList(0),
  fInitQAFlag(kFALSE) {
  //Default constructor
  if(!fInitQAFlag) InitQA();
  fHistEvents = new TH1I("fHistEvents","Analyzed events",1,0,1);

  fHistYPtProtons = new TH2D("fHistYPtProtons","Protons",
			     fNBinsY,gY,fNBinsPt,gPt);
  fHistYPtProtons->SetStats(kTRUE);
  fHistYPtProtons->GetYaxis()->SetTitle("P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fHistYPtProtons->GetXaxis()->SetTitle("#eta");
  else
    fHistYPtProtons->GetXaxis()->SetTitle("y");
  fHistYPtProtons->GetXaxis()->SetTitleColor(1);

  fHistYPtAntiProtons = new TH2D("fHistYPtAntiProtons","Antiprotons",
				 fNBinsY,gY,fNBinsPt,gPt);
  fHistYPtAntiProtons->SetStats(kTRUE);
  fHistYPtAntiProtons->GetYaxis()->SetTitle("P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fHistYPtAntiProtons->GetXaxis()->SetTitle("#eta");
  else
    fHistYPtAntiProtons->GetXaxis()->SetTitle("y");
  fHistYPtAntiProtons->GetXaxis()->SetTitleColor(1);

  //setting up the containers
  Int_t iBin[2];
  iBin[0] = nbinsY;
  iBin[1] = nbinsPt;
  fProtonContainer = new AliCFContainer("containerProtons",
					"container for protons",
					kNSteps,2,iBin);
  fProtonContainer->SetBinLimits(0,gY); //rapidity or eta
  fProtonContainer->SetBinLimits(1,gPt); //pT
  fAntiProtonContainer = new AliCFContainer("containerAntiProtons",
					    "container for antiprotons",
					    kNSteps,2,iBin);
  fAntiProtonContainer->SetBinLimits(0,gY); //rapidity or eta
  fAntiProtonContainer->SetBinLimits(1,gPt); //pT
} 

//____________________________________________________________________//
AliProtonAnalysis::~AliProtonAnalysis() {
  //Default destructor
  if(fProtonAnalysisBase) delete fProtonAnalysisBase;

  if(fHistEvents) delete fHistEvents;
  if(fHistYPtProtons) delete fHistYPtProtons;
  if(fHistYPtAntiProtons) delete fHistYPtAntiProtons;
  if(fProtonContainer) delete fProtonContainer;
  if(fAntiProtonContainer) delete fAntiProtonContainer;

  if(fEffGridListProtons) delete fEffGridListProtons;
  if(fCorrectionListProtons2D) delete fCorrectionListProtons2D;
  if(fEfficiencyListProtons1D) delete fEfficiencyListProtons1D;
  if(fCorrectionListProtons1D) delete fCorrectionListProtons1D;
  if(fEffGridListAntiProtons) delete fEffGridListAntiProtons;
  if(fCorrectionListAntiProtons2D) delete fCorrectionListAntiProtons2D;
  if(fEfficiencyListAntiProtons1D) delete fEfficiencyListAntiProtons1D;
  if(fCorrectionListAntiProtons1D) delete fCorrectionListAntiProtons1D;
  if(fCorrectProtons) delete fCorrectProtons;
  if(fCorrectAntiProtons) delete fCorrectAntiProtons;

  //QA lists
  if(fGlobalQAList) delete fGlobalQAList;
  if(fQA2DList) delete fQA2DList;
  if(fQAProtonsAcceptedList) delete fQAProtonsAcceptedList; 
  if(fQAProtonsRejectedList) delete fQAProtonsRejectedList;
  if(fQAAntiProtonsAcceptedList) delete fQAAntiProtonsAcceptedList; 
  if(fQAAntiProtonsRejectedList) delete fQAAntiProtonsRejectedList;
}

//____________________________________________________________________//
void AliProtonAnalysis::InitAnalysisHistograms(Int_t nbinsY, 
					       Float_t fLowY, Float_t fHighY, 
					       Int_t nbinsPt, 
					       Float_t fLowPt, Float_t fHighPt) {
  //Initializes the histograms
  if(!fInitQAFlag) InitQA();
  fNBinsY = nbinsY;
  fMinY = fLowY;
  fMaxY = fHighY;
  fNBinsPt = nbinsPt;
  fMinPt = fLowPt;
  fMaxPt = fHighPt;

  fHistEvents = new TH1I("fHistEvents","Analyzed events",1,0,1);

  fHistYPtProtons = new TH2D("fHistYPtProtons","Protons",
			     fNBinsY,fMinY,fMaxY,
			     fNBinsPt,fMinPt,fMaxPt);
  fHistYPtProtons->SetStats(kTRUE);
  fHistYPtProtons->GetYaxis()->SetTitle("P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fHistYPtProtons->GetXaxis()->SetTitle("#eta");
  else
    fHistYPtProtons->GetXaxis()->SetTitle("y");
  fHistYPtProtons->GetXaxis()->SetTitleColor(1);

  fHistYPtAntiProtons = new TH2D("fHistYPtAntiProtons","Antiprotons",
				 fNBinsY,fMinY,fMaxY,
				 fNBinsPt,fMinPt,fMaxPt);
  fHistYPtAntiProtons->SetStats(kTRUE);
  fHistYPtAntiProtons->GetYaxis()->SetTitle("P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fHistYPtAntiProtons->GetXaxis()->SetTitle("#eta");
  else
    fHistYPtAntiProtons->GetXaxis()->SetTitle("y");
  fHistYPtAntiProtons->GetXaxis()->SetTitleColor(1);

  //setting up the containers
  Int_t iBin[2];
  iBin[0] = nbinsY;
  iBin[1] = nbinsPt;
  Double_t *binLimY = new Double_t[nbinsY+1];
  Double_t *binLimPt = new Double_t[nbinsPt+1];
  //values for bin lower bounds
  for(Int_t i = 0; i <= nbinsY; i++) 
    binLimY[i]=(Double_t)fLowY  + (fHighY - fLowY)  /nbinsY*(Double_t)i;
  for(Int_t i = 0; i <= nbinsPt; i++) 
    binLimPt[i]=(Double_t)fLowPt  + (fHighPt - fLowPt)  /nbinsPt*(Double_t)i;

  fProtonContainer = new AliCFContainer("containerProtons",
					"container for protons",
					kNSteps,2,iBin);
  fProtonContainer->SetBinLimits(0,binLimY); //rapidity
  fProtonContainer->SetBinLimits(1,binLimPt); //pT
  fAntiProtonContainer = new AliCFContainer("containerAntiProtons",
					    "container for antiprotons",
					    kNSteps,2,iBin);
  fAntiProtonContainer->SetBinLimits(0,binLimY); //rapidity
  fAntiProtonContainer->SetBinLimits(1,binLimPt); //pT
}

//____________________________________________________________________//
void AliProtonAnalysis::InitAnalysisHistograms(Int_t nbinsY, Double_t *gY, 
					       Int_t nbinsPt, Double_t *gPt) {
  //Initializes the histograms using asymmetric values - global tracking
  if(!fInitQAFlag) InitQA();
  fNBinsY = nbinsY;
  fMinY = gY[0];
  fMaxY = gY[nbinsY];
  fNBinsPt = nbinsPt;
  fMinPt = gPt[0];
  fMaxPt = gPt[nbinsPt];

  fHistEvents = new TH1I("fHistEvents","Analyzed events",1,0,1);

  fHistYPtProtons = new TH2D("fHistYPtProtons","Protons",
			     fNBinsY,gY,fNBinsPt,gPt);
  fHistYPtProtons->SetStats(kTRUE);
  fHistYPtProtons->GetYaxis()->SetTitle("P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fHistYPtProtons->GetXaxis()->SetTitle("#eta");
  else
    fHistYPtProtons->GetXaxis()->SetTitle("y");
  fHistYPtProtons->GetXaxis()->SetTitleColor(1);

  fHistYPtAntiProtons = new TH2D("fHistYPtAntiProtons","Antiprotons",
				 fNBinsY,gY,fNBinsPt,gPt);
  fHistYPtAntiProtons->SetStats(kTRUE);
  fHistYPtAntiProtons->GetYaxis()->SetTitle("P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fHistYPtAntiProtons->GetXaxis()->SetTitle("#eta");
  else
    fHistYPtAntiProtons->GetXaxis()->SetTitle("y");
  fHistYPtAntiProtons->GetXaxis()->SetTitleColor(1);

  //setting up the containers
  Int_t iBin[2];
  iBin[0] = nbinsY;
  iBin[1] = nbinsPt;

  fProtonContainer = new AliCFContainer("containerProtons",
					"container for protons",
					kNSteps,2,iBin);
  fProtonContainer->SetBinLimits(0,gY); //rapidity
  fProtonContainer->SetBinLimits(1,gPt); //pT
  fAntiProtonContainer = new AliCFContainer("containerAntiProtons",
					    "container for antiprotons",
					    kNSteps,2,iBin);
  fAntiProtonContainer->SetBinLimits(0,gY); //rapidity
  fAntiProtonContainer->SetBinLimits(1,gPt); //pT
}

//____________________________________________________________________//
Bool_t AliProtonAnalysis::ReadFromFile(const char* filename) {
  //Read the containers from the existing file
  Bool_t status = kTRUE;

  TFile *file = TFile::Open(filename);
  if(!file) {
    cout<<"Could not find the input file "<<filename<<endl;
    status = kFALSE;
  }

  TList *list = (TList *)file->Get("outputList");
  if(list) {
    cout<<"Retrieving objects from the list "<<list->GetName()<<"..."<<endl; 
    fHistYPtProtons = (TH2D *)list->At(0);
    fHistYPtAntiProtons = (TH2D *)list->At(1);
    fHistEvents = (TH1I *)list->At(2);
    fProtonContainer = (AliCFContainer *)list->At(3);
    fAntiProtonContainer = (AliCFContainer *)list->At(4);
  }
  else if(!list) {
    cout<<"Retrieving objects from the file... "<<endl;
    fHistYPtProtons = (TH2D *)file->Get("fHistYPtProtons");
    fHistYPtAntiProtons = (TH2D *)file->Get("fHistYPtAntiProtons");
    fHistEvents = (TH1I *)file->Get("fHistEvents");
    fProtonContainer = (AliCFContainer *)file->Get("containerProtons");
    fAntiProtonContainer = (AliCFContainer *)file->Get("containerAntiProtons");
  }
  if((!fHistYPtProtons)||(!fHistYPtAntiProtons)||(!fHistEvents)
     ||(!fProtonContainer)||(!fAntiProtonContainer)) {
    cout<<"Input containers were not found!!!"<<endl;
    status = kFALSE;
  }
  else {
    //fHistYPtProtons = fProtonContainer->ShowProjection(0,1,0);
    //fHistYPtAntiProtons = fAntiProtonContainer->ShowProjection(0,1,0);
    fHistYPtProtons->Sumw2();
    fHistYPtAntiProtons->Sumw2();
  }

  return status;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetProtonYHistogram() {
  //Get the y histogram for protons
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();

  TH1D *fYProtons = (TH1D *)fHistYPtProtons->ProjectionX("fYProtons",0,fHistYPtProtons->GetYaxis()->GetNbins(),"");
  //TH1D *fYProtons = fProtonContainer->ShowProjection(0,0); //variable-step
   
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
  //Get the y histogram for antiprotons
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();
  
  TH1D *fYAntiProtons = (TH1D *)fHistYPtAntiProtons->ProjectionX("fYAntiProtons",0,fHistYPtAntiProtons->GetYaxis()->GetNbins(),"");
  //TH1D *fYAntiProtons = fAntiProtonContainer->ShowProjection(0,0);//variable-step 
 
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
  //Get the Pt histogram for protons
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();
  
  TH1D *fPtProtons = (TH1D *)fHistYPtProtons->ProjectionY("fPtProtons",0,fHistYPtProtons->GetXaxis()->GetNbins(),""); 
  //TH1D *fPtProtons = fProtonContainer->ShowProjection(1,0); //variable-step

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
  //Get the Pt histogram for antiprotons
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();
  
  TH1D *fPtAntiProtons = (TH1D *)fHistYPtAntiProtons->ProjectionY("fPtAntiProtons",0,fHistYPtProtons->GetXaxis()->GetNbins(),""); 
  //TH1D *fPtAntiProtons = fAntiProtonContainer->ShowProjection(1,0); //variable-step

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
TH1D *AliProtonAnalysis::GetProtonCorrectedYHistogram() {
  //Get the corrected y histogram for protons
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();

  TH1D *fYProtons = fCorrectProtons->Project(0); //0: rapidity
   
  fYProtons->SetStats(kFALSE);
  fYProtons->GetYaxis()->SetTitle("(1/N_{events})(dN/dy)");
  fYProtons->GetXaxis()->SetTitle("y");
  fYProtons->SetTitle("dN/dy protons");
  fYProtons->SetMarkerStyle(kFullCircle);
  fYProtons->SetMarkerColor(4);
  if(nAnalyzedEvents > 0)
    fYProtons->Scale(1./nAnalyzedEvents);
  
  return fYProtons;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetAntiProtonCorrectedYHistogram() {
  //Get the corrected y histogram for antiprotons
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();

  TH1D *fYAntiProtons = fCorrectAntiProtons->Project(0); //0: rapidity
   
  fYAntiProtons->SetStats(kFALSE);
  fYAntiProtons->GetYaxis()->SetTitle("(1/N_{events})(dN/dy)");
  fYAntiProtons->GetXaxis()->SetTitle("y");
  fYAntiProtons->SetTitle("dN/dy protons");
  fYAntiProtons->SetMarkerStyle(kFullCircle);
  fYAntiProtons->SetMarkerColor(4);
  if(nAnalyzedEvents > 0)
    fYAntiProtons->Scale(1./nAnalyzedEvents);
  
  return fYAntiProtons;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetProtonCorrectedPtHistogram() {
  //Get the corrected Pt histogram for protons
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();

  TH1D *fPtProtons = fCorrectProtons->Project(0); //0: rapidity
   
  fPtProtons->SetStats(kFALSE);
  fPtProtons->GetYaxis()->SetTitle("(1/N_{events})(dN/dP_{T})");
  fPtProtons->GetXaxis()->SetTitle("P_{T} [GeV/c]");
  fPtProtons->SetTitle("dN/dPt protons");
  fPtProtons->SetMarkerStyle(kFullCircle);
  fPtProtons->SetMarkerColor(4);
  if(nAnalyzedEvents > 0)
    fPtProtons->Scale(1./nAnalyzedEvents);
  
  return fPtProtons;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetAntiProtonCorrectedPtHistogram() {
  //Get the corrected Pt histogram for antiprotons
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();

  TH1D *fPtAntiProtons = fCorrectAntiProtons->Project(0); //0: rapidity
   
  fPtAntiProtons->SetStats(kFALSE);
  fPtAntiProtons->GetYaxis()->SetTitle("(1/N_{events})(dN/dP_{T})");
  fPtAntiProtons->GetXaxis()->SetTitle("P_{T} [GeV/c]");
  fPtAntiProtons->SetTitle("dN/dPt antiprotons");
  fPtAntiProtons->SetMarkerStyle(kFullCircle);
  fPtAntiProtons->SetMarkerColor(4);
  if(nAnalyzedEvents > 0)
    fPtAntiProtons->Scale(1./nAnalyzedEvents);
  
  return fPtAntiProtons;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetYRatioHistogram() {
  //Returns the rapidity dependence of the ratio (uncorrected)
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
TH1D *AliProtonAnalysis::GetYRatioCorrectedHistogram(TH2D *gCorrectionProtons, 
						     TH2D *gCorrectionAntiProtons) {
  //Returns the rapidity dependence of the ratio (corrected)
  fHistYPtProtons->Multiply(gCorrectionProtons);
  TH1D *fYProtons = GetProtonYHistogram();
  fHistYPtAntiProtons->Multiply(gCorrectionAntiProtons);
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
  //Returns the pT dependence of the ratio (uncorrected)
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
TH1D *AliProtonAnalysis::GetPtRatioCorrectedHistogram(TH2D *gCorrectionProtons, 
						      TH2D *gCorrectionAntiProtons) {
  //Returns the Pt dependence of the ratio (corrected)
  fHistYPtProtons->Multiply(gCorrectionProtons);
  TH1D *fPtProtons = GetProtonPtHistogram();
  fHistYPtAntiProtons->Multiply(gCorrectionAntiProtons);
  TH1D *fPtAntiProtons = GetAntiProtonPtHistogram();
  
  TH1D *hRatioPt = new TH1D("hRatioPt","",fPtProtons->GetNbinsX(),fPtProtons->GetXaxis()->GetXmin(),fPtProtons->GetXaxis()->GetXmax());
  hRatioPt->Divide(fPtAntiProtons,fPtProtons,1.0,1.0);
  hRatioPt->SetMarkerStyle(kFullCircle);
  hRatioPt->SetMarkerColor(4);
  hRatioPt->GetYaxis()->SetTitle("#bar{p}/p");
  hRatioPt->GetYaxis()->SetTitleOffset(1.4);
  hRatioPt->GetXaxis()->SetTitle("y");
  hRatioPt->GetXaxis()->SetTitleColor(1);
  hRatioPt->SetStats(kFALSE);

  return hRatioPt;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetYAsymmetryHistogram() {
  //Returns the rapidity dependence of the asymmetry (uncorrected)
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
  //Returns the pT dependence of the asymmetry (uncorrected)
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
void AliProtonAnalysis::Analyze(AliESDEvent* esd,
				const AliESDVertex *vertex) {
  //Main analysis part - ESD
  Int_t nTracks = 0;
  Int_t nIdentifiedProtons = 0, nIdentifiedAntiProtons = 0;
  Int_t nSurvivedProtons = 0, nSurvivedAntiProtons = 0;

  fHistEvents->Fill(0); //number of analyzed events
  Double_t containerInput[2] ;
  Double_t gPt = 0.0, gP = 0.0;
  nTracks = esd->GetNumberOfTracks();
  for(Int_t iTracks = 0; iTracks < nTracks; iTracks++) {
    AliESDtrack* track = esd->GetTrack(iTracks);
    AliESDtrack trackTPC;

    //in case it's a TPC only track relate it to the proper vertex
    /*if(fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kTPC) {
      Float_t p[2],cov[3];
      track->GetImpactParametersTPC(p,cov);
      if (p[0]==0 && p[1]==0)  
	track->RelateToVertexTPC(((AliESDEvent*)esd)->GetPrimaryVertexTPC(),esd->GetMagneticField(),kVeryBig);
      if (!track->FillTPCOnlyTrack(trackTPC)) {
	continue;
      }
      track = &trackTPC ;
      }*/

    if((fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kTPC)||(fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kHybrid)) {
      AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
      if(!tpcTrack) continue;
      gPt = tpcTrack->Pt();
      gP = tpcTrack->P();
      
      if(fProtonAnalysisBase->IsAccepted(esd,vertex,track)) 
	((TH2F *)(fQA2DList->At(0)))->Fill(gP,track->GetTPCsignal());

      if(fProtonAnalysisBase->IsProton(track)) {
	if(fProtonAnalysisBase->IsAccepted(esd,vertex,track)) 
	  ((TH2F *)(fQA2DList->At(1)))->Fill(gP,track->GetTPCsignal());
	FillQA(esd,vertex,track);
	if(tpcTrack->Charge() > 0) {
	  nIdentifiedProtons += 1;
	  if(fProtonAnalysisBase->GetEtaMode())
	    containerInput[0] = tpcTrack->Eta();
	  else
	    containerInput[0] = fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
							      tpcTrack->Py(),
							      tpcTrack->Pz());
	  containerInput[1] = gPt;
	  fProtonContainer->Fill(containerInput,0);   
	  
	  if(!fProtonAnalysisBase->IsAccepted(esd,vertex,track)) continue;//track cuts
	  if(fProtonAnalysisBase->GetEtaMode())
	    containerInput[0] = tpcTrack->Eta();
	  else
	    containerInput[0] = fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
							      tpcTrack->Py(),
							      tpcTrack->Pz());
	  containerInput[1] = gPt;
	  fProtonContainer->Fill(containerInput,1);   
	  
	  if(!fProtonAnalysisBase->IsInPhaseSpace(track)) continue; //track outside the analyzed y-Pt
	  nSurvivedProtons += 1;
	  if(fProtonAnalysisBase->GetEtaMode()) {
	    fHistYPtProtons->Fill(tpcTrack->Eta(),
				  gPt);
	    containerInput[0] = tpcTrack->Eta();
	  }
	  else {
	    fHistYPtProtons->Fill(fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
								tpcTrack->Py(),
								tpcTrack->Pz()),
				  gPt);
	    //fill the container
	    containerInput[0] = fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
							      tpcTrack->Py(),
							      tpcTrack->Pz());
	  }
	  containerInput[1] = gPt;
	  fProtonContainer->Fill(containerInput,2);   
	}//protons
	else if(tpcTrack->Charge() < 0) {
	  nIdentifiedAntiProtons += 1;
	  if(fProtonAnalysisBase->GetEtaMode())
	    containerInput[0] = tpcTrack->Eta();
	  else
	    containerInput[0] = fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
							      tpcTrack->Py(),
							      tpcTrack->Pz());
	  containerInput[1] = gPt;
	  fAntiProtonContainer->Fill(containerInput,0);   
	  
	  if(!fProtonAnalysisBase->IsAccepted(esd,vertex,track)) continue;//track cuts
	  if(fProtonAnalysisBase->GetEtaMode())
	    containerInput[0] = tpcTrack->Eta();
	  else
	    containerInput[0] = fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
							      tpcTrack->Py(),
							      tpcTrack->Pz());
	  containerInput[1] = gPt;
	  fAntiProtonContainer->Fill(containerInput,1);   
	  
	  if(!fProtonAnalysisBase->IsInPhaseSpace(track)) continue; //track outside the analyzed y-Pt
	  nSurvivedAntiProtons += 1;
	  if(fProtonAnalysisBase->GetEtaMode()) {
	    fHistYPtAntiProtons->Fill(tpcTrack->Eta(),
				      gPt);
	    containerInput[0] = tpcTrack->Eta();
	  }
	  else {
	    fHistYPtAntiProtons->Fill(fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
								    tpcTrack->Py(),
								    tpcTrack->Pz()),
				      gPt);
	    //fill the container
	    containerInput[0] = fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
							      tpcTrack->Py(),
							      tpcTrack->Pz());
	  }
	  containerInput[1] = gPt;
	  fAntiProtonContainer->Fill(containerInput,2);
	}//antiprotons   
      }//proton check
    }//TPC only tracks
    else if(fProtonAnalysisBase->GetAnalysisMode() == AliProtonAnalysisBase::kGlobal) {
      gPt = track->Pt();
      gP = track->P();
      
      if(fProtonAnalysisBase->IsAccepted(esd,vertex,track)) 
	((TH2F *)(fQA2DList->At(0)))->Fill(gP,track->GetTPCsignal());
      if(fProtonAnalysisBase->IsProton(track)) {
	if(fProtonAnalysisBase->IsAccepted(esd,vertex,track)) 
	  ((TH2F *)(fQA2DList->At(1)))->Fill(gP,track->GetTPCsignal());
	FillQA(esd,vertex,track);
	if(track->Charge() > 0) {
	  nIdentifiedProtons += 1;
	  if(fProtonAnalysisBase->GetEtaMode())
	    containerInput[0] = track->Eta();
	  else
	    containerInput[0] = fProtonAnalysisBase->Rapidity(track->Px(),
							      track->Py(),
							      track->Pz());
	  containerInput[1] = gPt;
	  fProtonContainer->Fill(containerInput,0);   
	    
	  if(!fProtonAnalysisBase->IsAccepted(esd,vertex,track)) continue;//track cuts
	  if(fProtonAnalysisBase->GetEtaMode())
	    containerInput[0] = track->Eta();
	  else
	    containerInput[0] = fProtonAnalysisBase->Rapidity(track->Px(),
							      track->Py(),
							      track->Pz());
	  containerInput[1] = gPt;
	  fProtonContainer->Fill(containerInput,1);   

	  if(!fProtonAnalysisBase->IsInPhaseSpace(track)) continue; //track outside the analyzed y-Pt
	  nSurvivedProtons += 1;
	  if(fProtonAnalysisBase->GetEtaMode()) {
	    fHistYPtProtons->Fill(track->Eta(),
				  gPt);
	    containerInput[0] = track->Eta();
	  }
	  else {
	    fHistYPtProtons->Fill(fProtonAnalysisBase->Rapidity(track->Px(),
								track->Py(),
								track->Pz()),
				  gPt);
	    //fill the container
	    containerInput[0] = fProtonAnalysisBase->Rapidity(track->Px(),
							      track->Py(),
							      track->Pz());
	  }
	  containerInput[1] = gPt;
	  fProtonContainer->Fill(containerInput,2);   
	}//protons
	else if(track->Charge() < 0) {
	  nIdentifiedAntiProtons += 1;
	  if(fProtonAnalysisBase->GetEtaMode())
	    containerInput[0] = track->Eta();
	  else
	    containerInput[0] = fProtonAnalysisBase->Rapidity(track->Px(),
							      track->Py(),
							      track->Pz());
	  containerInput[1] = gPt;
	  fAntiProtonContainer->Fill(containerInput,0);   

	  if(!fProtonAnalysisBase->IsAccepted(esd,vertex,track)) continue;//track cuts
	  if(fProtonAnalysisBase->GetEtaMode())
	    containerInput[0] = track->Eta();
	  else
	    containerInput[0] = fProtonAnalysisBase->Rapidity(track->Px(),
							      track->Py(),
							      track->Pz());
	  containerInput[1] = gPt;
	  fAntiProtonContainer->Fill(containerInput,1);   
	  
	  if(!fProtonAnalysisBase->IsInPhaseSpace(track)) continue; //track outside the analyzed y-Pt
	  nSurvivedAntiProtons += 1;
	  if(fProtonAnalysisBase->GetEtaMode()) {
	    fHistYPtAntiProtons->Fill(track->Eta(),
				      gPt);
	    containerInput[0] = track->Eta();
	  }
	  else {
	    fHistYPtAntiProtons->Fill(fProtonAnalysisBase->Rapidity(track->Px(),
								    track->Py(),
								    track->Pz()),
				      gPt);
	    //fill the container
	    containerInput[0] = fProtonAnalysisBase->Rapidity(track->Px(),
							      track->Py(),
							      track->Pz());
	  }
	  containerInput[1] = gPt;
	  fAntiProtonContainer->Fill(containerInput,2);   
	}//antiprotons
      }//proton check 
    }//combined tracking
  }//track loop 
  
  if(fProtonAnalysisBase->GetDebugMode())
    Printf("Initial number of tracks: %d | Identified (anti)protons: %d - %d | Survived (anti)protons: %d - %d",nTracks,nIdentifiedProtons,nIdentifiedAntiProtons,nSurvivedProtons,nSurvivedAntiProtons);
}

//____________________________________________________________________//
void AliProtonAnalysis::Analyze(AliAODEvent* const fAOD) {
  //Main analysis part - AOD
  fHistEvents->Fill(0); //number of analyzed events
  Int_t nTracks = fAOD->GetNumberOfTracks();
  for(Int_t iTracks = 0; iTracks < nTracks; iTracks++) {
    AliAODTrack* track = fAOD->GetTrack(iTracks);
    Double_t gPt = track->Pt();
    Double_t gP = track->P();
    
    //pid
    Double_t probability[10];
    track->GetPID(probability);
    Double_t rcc = 0.0;
    for(Int_t i = 0; i < AliPID::kSPECIESN; i++) rcc += probability[i]*fProtonAnalysisBase->GetParticleFraction(i,gP);
    if(rcc == 0.0) continue;
    Double_t w[10];
    for(Int_t i = 0; i < AliPID::kSPECIESN; i++) w[i] = probability[i]*fProtonAnalysisBase->GetParticleFraction(i,gP)/rcc;
    Long64_t fParticleType = TMath::LocMax(AliPID::kSPECIESN,w);
    if(fParticleType == 4) {
      if(track->Charge() > 0) 
	fHistYPtProtons->Fill(track->Y(fParticleType),gPt);
      else if(track->Charge() < 0) 
	fHistYPtAntiProtons->Fill(track->Y(fParticleType),gPt);
    }//proton check
  }//track loop 
}

//____________________________________________________________________//
void AliProtonAnalysis::Analyze(AliStack* const stack, 
				Bool_t iInclusive) {
  //Main analysis part - MC
  fHistEvents->Fill(0); //number of analyzed events

  Int_t nParticles = 0;
  //inclusive protons - 
  if(iInclusive) nParticles = stack->GetNtrack();
  else nParticles = stack->GetNprimary();

  for(Int_t i = 0; i < nParticles; i++) {
    TParticle *particle = stack->Particle(i);
    if(!particle) continue;

    //in case of inclusive protons reject the secondaries from hadronic inter.
    if(particle->GetUniqueID() == 13) continue;

    if(TMath::Abs(particle->Eta()) > 1.0) continue;
    if((particle->Pt() > fMaxPt)||(particle->Pt() < fMinPt)) continue;
    if((fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()) > fMaxY)||(fProtonAnalysisBase->Rapidity(particle->Px(),particle->Py(),particle->Pz()) < fMinY)) continue;

    Int_t pdgcode = particle->GetPdgCode();
    if(pdgcode == 2212) fHistYPtProtons->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
									    particle->Py(),
									    particle->Pz()),
					      particle->Pt());
    if(pdgcode == -2212) fHistYPtAntiProtons->Fill(fProtonAnalysisBase->Rapidity(particle->Px(),
										 particle->Py(),
										 particle->Pz()),
						   particle->Pt());
  }//particle loop
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
void AliProtonAnalysis::Correct(Int_t step) {
  //Applies the correction maps to the initial containers
  fCorrectProtons = new AliCFDataGrid("correctProtons",
				      "corrected data",
				      *fProtonContainer);
  fCorrectProtons->SetMeasured(0);
  fCorrectProtons->ApplyEffCorrection(*(AliCFEffGrid *)fEffGridListProtons->At(step));

  fCorrectAntiProtons = new AliCFDataGrid("correctAntiProtons",
					  "corrected data",
					  *fAntiProtonContainer);
  fCorrectAntiProtons->SetMeasured(0);
  fCorrectAntiProtons->ApplyEffCorrection(*(AliCFEffGrid *)fEffGridListAntiProtons->At(step));
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

  //________________________________________//
  //Protons
  fEffGridListProtons = new TList();
  fCorrectionListProtons2D = new TList(); 
  fEfficiencyListProtons1D = new TList(); 
  fCorrectionListProtons1D = new TList();
  
  AliCFContainer *corrfwContainerProtons = (AliCFContainer*) (file->Get("containerProtons"));
  if(!corrfwContainerProtons) {
    cout<<"CORRFW container for protons not found!"<<endl;
    status = kFALSE;
  }
  
  Int_t nSteps = corrfwContainerProtons->GetNStep();
  TH2D *gYPt[4];
  //currently the GRID is formed by the y-pT parameters
  //Add Vz as a next step
  Int_t iRap = 0, iPt = 1;
  AliCFEffGrid *effProtonsStep0Step1 = new AliCFEffGrid("eff10",
					 "effProtonsStep0Step1",
					 *corrfwContainerProtons);
  effProtonsStep0Step1->CalculateEfficiency(1,0); //eff= step1/step0
  fEffGridListProtons->Add(effProtonsStep0Step1);
  gYPt[0] = effProtonsStep0Step1->Project(iRap,iPt);
  fCorrectionListProtons2D->Add(gYPt[0]);
  
  AliCFEffGrid *effProtonsStep0Step2 = new AliCFEffGrid("eff20",
					 "effProtonsStep0Step2",
					 *corrfwContainerProtons);
  effProtonsStep0Step2->CalculateEfficiency(2,0); //eff= step2/step0
  fEffGridListProtons->Add(effProtonsStep0Step2);
  gYPt[1] = effProtonsStep0Step2->Project(iRap,iPt);
  fCorrectionListProtons2D->Add(gYPt[1]);

  AliCFEffGrid *effProtonsStep0Step3 = new AliCFEffGrid("eff30",
					 "effProtonsStep0Step3",
					 *corrfwContainerProtons);
  effProtonsStep0Step3->CalculateEfficiency(3,0); //eff= step1/step0
  fEffGridListProtons->Add(effProtonsStep0Step3);
  gYPt[2] = effProtonsStep0Step3->Project(iRap,iPt);
  fCorrectionListProtons2D->Add(gYPt[2]);

  TH1D *gEfficiency[2][3]; //efficiency as a function of pT and of y (raws-[2])
  TH1D *gCorrection[2][3]; //efficiency as a function of pT and of y (raws-[2])
  TString gTitle;
  //Get the projection of the efficiency maps
  for(Int_t iParameter = 0; iParameter < 2; iParameter++) {
    gEfficiency[iParameter][0] = effProtonsStep0Step1->Project(iParameter);
    gTitle = "ProtonsEfficiency_Parameter"; gTitle += iParameter+1;
    gTitle += "_Step0_Step1"; 
    gEfficiency[iParameter][0]->SetName(gTitle.Data());
    fEfficiencyListProtons1D->Add(gEfficiency[iParameter][0]);  
    gTitle = "ProtonsCorrection_Parameter"; gTitle += iParameter+1;
    gTitle += "_Step0_Step1"; 
    gCorrection[iParameter][0] = new TH1D(gTitle.Data(),
					  gTitle.Data(),
					  gEfficiency[iParameter][0]->GetNbinsX(),
					  gEfficiency[iParameter][0]->GetXaxis()->GetXmin(),
					  gEfficiency[iParameter][0]->GetXaxis()->GetXmax());
    //initialisation of the correction
    for(Int_t iBin = 1; iBin <= gEfficiency[iParameter][0]->GetNbinsX(); iBin++)
      gCorrection[iParameter][0]->SetBinContent(iBin,1.0);

    gEfficiency[iParameter][1] = effProtonsStep0Step2->Project(iParameter);
    gTitle = "ProtonsEfficiency_Parameter"; gTitle += iParameter+1;
    gTitle += "_Step0_Step2"; 
    gEfficiency[iParameter][1]->SetName(gTitle.Data());
    fEfficiencyListProtons1D->Add(gEfficiency[iParameter][1]);  
    gTitle = "ProtonsCorrection_Parameter"; gTitle += iParameter+1;
    gTitle += "_Step0_Step2"; 
    gCorrection[iParameter][1] = new TH1D(gTitle.Data(),
					  gTitle.Data(),
					  gEfficiency[iParameter][1]->GetNbinsX(),
					  gEfficiency[iParameter][1]->GetXaxis()->GetXmin(),
					  gEfficiency[iParameter][1]->GetXaxis()->GetXmax());
    //initialisation of the correction
    for(Int_t iBin = 1; iBin <= gEfficiency[iParameter][1]->GetNbinsX(); iBin++)
      gCorrection[iParameter][1]->SetBinContent(iBin,1.0);

    gEfficiency[iParameter][2] = effProtonsStep0Step3->Project(iParameter);
    gTitle = "ProtonsEfficiency_Parameter"; gTitle += iParameter+1;
    gTitle += "_Step0_Step3"; 
    gEfficiency[iParameter][2]->SetName(gTitle.Data());
    fEfficiencyListProtons1D->Add(gEfficiency[iParameter][2]);  
    gTitle = "ProtonsCorrection_Parameter"; gTitle += iParameter+1;
    gTitle += "_Step0_Step3"; 
    gCorrection[iParameter][2] = new TH1D(gTitle.Data(),
					  gTitle.Data(),
					  gEfficiency[iParameter][2]->GetNbinsX(),
					  gEfficiency[iParameter][2]->GetXaxis()->GetXmin(),
					  gEfficiency[iParameter][2]->GetXaxis()->GetXmax());
    //initialisation of the correction
    for(Int_t iBin = 1; iBin <= gEfficiency[iParameter][2]->GetNbinsX(); iBin++)
      gCorrection[iParameter][2]->SetBinContent(iBin,1.0);
  }//parameter loop
  //Calculate the 1D correction parameters as a function of y and pT
  for(Int_t iParameter = 0; iParameter < 2; iParameter++) { 
    for(Int_t iStep = 1; iStep < nSteps; iStep++) { 
      gCorrection[iParameter][iStep-1]->Divide(gEfficiency[iParameter][iStep-1]);
      fCorrectionListProtons1D->Add(gCorrection[iParameter][iStep-1]);  
    }
  }

  //________________________________________//
  //AntiProtons
  fEffGridListAntiProtons = new TList();
  fCorrectionListAntiProtons2D = new TList(); 
  fEfficiencyListAntiProtons1D = new TList(); 
  fCorrectionListAntiProtons1D = new TList();
  
  AliCFContainer *corrfwContainerAntiProtons = (AliCFContainer*) (file->Get("containerAntiProtons"));
  if(!corrfwContainerAntiProtons) {
    cout<<"CORRFW container for antiprotons not found!"<<endl;
    status = kFALSE;
  }
  
  nSteps = corrfwContainerAntiProtons->GetNStep();
  //currently the GRID is formed by the y-pT parameters
  //Add Vz as a next step
  AliCFEffGrid *effAntiProtonsStep0Step1 = new AliCFEffGrid("eff10",
					 "effAntiProtonsStep0Step1",
					 *corrfwContainerAntiProtons);
  effAntiProtonsStep0Step1->CalculateEfficiency(1,0); //eff= step1/step0
  fEffGridListAntiProtons->Add(effAntiProtonsStep0Step1);
  gYPt[0] = effAntiProtonsStep0Step1->Project(iRap,iPt);
  fCorrectionListAntiProtons2D->Add(gYPt[0]);
  
  AliCFEffGrid *effAntiProtonsStep0Step2 = new AliCFEffGrid("eff20",
					 "effAntiProtonsStep0Step2",
					 *corrfwContainerAntiProtons);
  effAntiProtonsStep0Step2->CalculateEfficiency(2,0); //eff= step2/step0
  fEffGridListAntiProtons->Add(effAntiProtonsStep0Step2);
  gYPt[1] = effAntiProtonsStep0Step2->Project(iRap,iPt);
  fCorrectionListAntiProtons2D->Add(gYPt[1]);

  AliCFEffGrid *effAntiProtonsStep0Step3 = new AliCFEffGrid("eff30",
					 "effAntiProtonsStep0Step3",
					 *corrfwContainerAntiProtons);
  effAntiProtonsStep0Step3->CalculateEfficiency(3,0); //eff= step1/step0
  fEffGridListAntiProtons->Add(effAntiProtonsStep0Step3);
  gYPt[2] = effAntiProtonsStep0Step3->Project(iRap,iPt);
  fCorrectionListAntiProtons2D->Add(gYPt[2]);

  //Get the projection of the efficiency maps
  for(Int_t iParameter = 0; iParameter < 2; iParameter++) {
    gEfficiency[iParameter][0] = effAntiProtonsStep0Step1->Project(iParameter);
    gTitle = "AntiProtonsEfficiency_Parameter"; gTitle += iParameter+1;
    gTitle += "_Step0_Step1"; 
    gEfficiency[iParameter][0]->SetName(gTitle.Data());
    fEfficiencyListAntiProtons1D->Add(gEfficiency[iParameter][0]);  
    gTitle = "AntiProtonsCorrection_Parameter"; gTitle += iParameter+1;
    gTitle += "_Step0_Step1"; 
    gCorrection[iParameter][0] = new TH1D(gTitle.Data(),
					  gTitle.Data(),
					  gEfficiency[iParameter][0]->GetNbinsX(),
					  gEfficiency[iParameter][0]->GetXaxis()->GetXmin(),
					  gEfficiency[iParameter][0]->GetXaxis()->GetXmax());
    //initialisation of the correction
    for(Int_t iBin = 1; iBin <= gEfficiency[iParameter][0]->GetNbinsX(); iBin++)
      gCorrection[iParameter][0]->SetBinContent(iBin,1.0);

    gEfficiency[iParameter][1] = effAntiProtonsStep0Step2->Project(iParameter);
    gTitle = "AntiProtonsEfficiency_Parameter"; gTitle += iParameter+1;
    gTitle += "_Step0_Step2"; 
    gEfficiency[iParameter][1]->SetName(gTitle.Data());
    fEfficiencyListAntiProtons1D->Add(gEfficiency[iParameter][1]);  
    gTitle = "AntiProtonsCorrection_Parameter"; gTitle += iParameter+1;
    gTitle += "_Step0_Step2"; 
    gCorrection[iParameter][1] = new TH1D(gTitle.Data(),
					  gTitle.Data(),
					  gEfficiency[iParameter][1]->GetNbinsX(),
					  gEfficiency[iParameter][1]->GetXaxis()->GetXmin(),
					  gEfficiency[iParameter][1]->GetXaxis()->GetXmax());
    //initialisation of the correction
    for(Int_t iBin = 1; iBin <= gEfficiency[iParameter][1]->GetNbinsX(); iBin++)
      gCorrection[iParameter][1]->SetBinContent(iBin,1.0);

    gEfficiency[iParameter][2] = effAntiProtonsStep0Step3->Project(iParameter);
    gTitle = "AntiProtonsEfficiency_Parameter"; gTitle += iParameter+1;
    gTitle += "_Step0_Step3"; 
    gEfficiency[iParameter][2]->SetName(gTitle.Data());
    fEfficiencyListAntiProtons1D->Add(gEfficiency[iParameter][2]);  
    gTitle = "AntiProtonsCorrection_Parameter"; gTitle += iParameter+1;
    gTitle += "_Step0_Step3"; 
    gCorrection[iParameter][2] = new TH1D(gTitle.Data(),
					  gTitle.Data(),
					  gEfficiency[iParameter][2]->GetNbinsX(),
					  gEfficiency[iParameter][2]->GetXaxis()->GetXmin(),
					  gEfficiency[iParameter][2]->GetXaxis()->GetXmax());
    //initialisation of the correction
    for(Int_t iBin = 1; iBin <= gEfficiency[iParameter][2]->GetNbinsX(); iBin++)
      gCorrection[iParameter][2]->SetBinContent(iBin,1.0);
  }//parameter loop
  //Calculate the 1D correction parameters as a function of y and pT
  for(Int_t iParameter = 0; iParameter < 2; iParameter++) { 
    for(Int_t iStep = 1; iStep < nSteps; iStep++) { 
      gCorrection[iParameter][iStep-1]->Divide(gEfficiency[iParameter][iStep-1]);
      fCorrectionListAntiProtons1D->Add(gCorrection[iParameter][iStep-1]);  
    }
  }

  return status;
}
 
//____________________________________________________________________//
void AliProtonAnalysis::InitQA() {
  //Applies the correction maps to the initial containers
  fInitQAFlag = kTRUE;
  fGlobalQAList = new TList();
  fGlobalQAList->SetName("fGlobalQAList");

  //========================================================//
  fQA2DList = new TList();
  fQA2DList->SetName("fQA2DList");
  fGlobalQAList->Add(fQA2DList);
  //dEdx plots
  TH2F *gHistdEdxP = new TH2F("gHistdEdxP","dE/dx (TPC); P [GeV/c]; dE/dx [a.u]",100,0.01,10.1,100,0,600);
  fQA2DList->Add(gHistdEdxP);
  TH2F *gHistProtonsdEdxP = new TH2F("gHistProtonsdEdxP","Accepted protons dE/dx (TPC); P [GeV/c]; dE/dx [a.u]",100,0.01,10.1,100,0,600);
  fQA2DList->Add(gHistProtonsdEdxP);
  
  //========================================================//
  fQAProtonsAcceptedList = new TList();
  fQAProtonsAcceptedList->SetName("fQAProtonsAcceptedList");
  fGlobalQAList->Add(fQAProtonsAcceptedList);
  //Accepted protons
  TH1F *gProtonsITSClustersPass = new TH1F("gProtonsITSClustersPass",
					    ";N_{clusters} (ITS);Entries",
					    7,0,7);
  fQAProtonsAcceptedList->Add(gProtonsITSClustersPass);
  TH1F *gProtonsChi2PerClusterITSPass = new TH1F("gProtonsChi2PerClusterITSPass",
						 ";x^{2}/N_{clusters} (ITS);Entries",
						 100,0,4);
  fQAProtonsAcceptedList->Add(gProtonsChi2PerClusterITSPass);
  TH1F *gProtonsTPCClustersPass = new TH1F("gProtonsTPCClustersPass",
					   ";N_{clusters} (TPC);Entries",
					   100,0,200);
  fQAProtonsAcceptedList->Add(gProtonsTPCClustersPass);
  TH1F *gProtonsChi2PerClusterTPCPass = new TH1F("gProtonsChi2PerClusterTPCPass",
						 ";x^{2}/N_{clusters} (TPC);Entries",
						 100,0,4);
  fQAProtonsAcceptedList->Add(gProtonsChi2PerClusterTPCPass);
  TH1F *gProtonsExtCov11Pass = new TH1F("gProtonsExtCov11Pass",
					";#sigma_{y} [cm];Entries",
					100,0,4);
  fQAProtonsAcceptedList->Add(gProtonsExtCov11Pass);
  TH1F *gProtonsExtCov22Pass = new TH1F("gProtonsExtCov22Pass",
					";#sigma_{z} [cm];Entries",
					100,0,4);
  fQAProtonsAcceptedList->Add(gProtonsExtCov22Pass);
  TH1F *gProtonsExtCov33Pass = new TH1F("gProtonsExtCov33Pass",
					";#sigma_{sin(#phi)};Entries",
					100,0,4);
  fQAProtonsAcceptedList->Add(gProtonsExtCov33Pass);
  TH1F *gProtonsExtCov44Pass = new TH1F("gProtonsExtCov44Pass",
					";#sigma_{tan(#lambda)};Entries",
					100,0,4);
  fQAProtonsAcceptedList->Add(gProtonsExtCov44Pass);
  TH1F *gProtonsExtCov55Pass = new TH1F("gProtonsExtCov55Pass",
					";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
					100,0,4);
  fQAProtonsAcceptedList->Add(gProtonsExtCov55Pass);
  TH1F *gProtonsSigmaToVertexPass = new TH1F("gProtonsSigmaToVertexPass",
					     ";#sigma_{Vertex};Entries",
					     100,0,10);
  fQAProtonsAcceptedList->Add(gProtonsSigmaToVertexPass);
  TH1F *gProtonsSigmaToVertexTPCPass = new TH1F("gProtonsSigmaToVertexTPCPass",
						";#sigma_{Vertex};Entries",
						100,0,10);
  fQAProtonsAcceptedList->Add(gProtonsSigmaToVertexTPCPass);
  TH1F *gProtonsDCAXYPass = new TH1F("gProtonsDCAXYPass",
				     ";DCA_{xy} [cm];Entries",
				     100,0,20);
  fQAProtonsAcceptedList->Add(gProtonsDCAXYPass);
  TH1F *gProtonsDCAXYTPCPass = new TH1F("gProtonsDCAXYTPCPass",
					";DCA_{xy} [cm];Entries",
					100,0,20);
  fQAProtonsAcceptedList->Add(gProtonsDCAXYTPCPass);
  TH1F *gProtonsDCAZPass = new TH1F("gProtonsDCAZPass",
				    ";DCA_{z} [cm];Entries",
				    100,0,20);
  fQAProtonsAcceptedList->Add(gProtonsDCAZPass);
  TH1F *gProtonsDCAZTPCPass = new TH1F("gProtonsDCAZTPCPass",
				       ";DCA_{z} [cm];Entries",
				       100,0,20);
  fQAProtonsAcceptedList->Add(gProtonsDCAZTPCPass);
  TH1F *gProtonsConstrainChi2Pass = new TH1F("gProtonsConstrainChi2Pass",
					     ";Log_{10}(#chi^{2});Entries",
					     100,-10,10);
  fQAProtonsAcceptedList->Add(gProtonsConstrainChi2Pass);
  TH1F *gProtonsITSRefitPass = new TH1F("gProtonsITSRefitPass",
					"",10,-1,1);
  fQAProtonsAcceptedList->Add(gProtonsITSRefitPass);
  TH1F *gProtonsTPCRefitPass = new TH1F("gProtonsTPCRefitPass",
					"",10,-1,1);
  fQAProtonsAcceptedList->Add(gProtonsTPCRefitPass);
  TH1F *gProtonsESDpidPass = new TH1F("gProtonsESDpidPass",
				      "",10,-1,1);
  fQAProtonsAcceptedList->Add(gProtonsESDpidPass);
  TH1F *gProtonsTPCpidPass = new TH1F("gProtonsTPCpidPass",
				      "",10,-1,1);
  fQAProtonsAcceptedList->Add(gProtonsTPCpidPass);
  TH1F *gProtonsPointOnITSLayer1Pass = new TH1F("gProtonsPointOnITSLayer1Pass",
						"",10,-1,1);
  fQAProtonsAcceptedList->Add(gProtonsPointOnITSLayer1Pass);
  TH1F *gProtonsPointOnITSLayer2Pass = new TH1F("gProtonsPointOnITSLayer2Pass",
						"",10,-1,1);
  fQAProtonsAcceptedList->Add(gProtonsPointOnITSLayer2Pass);
  TH1F *gProtonsPointOnITSLayer3Pass = new TH1F("gProtonsPointOnITSLayer3Pass",
						"",10,-1,1);
  fQAProtonsAcceptedList->Add(gProtonsPointOnITSLayer3Pass);
  TH1F *gProtonsPointOnITSLayer4Pass = new TH1F("gProtonsPointOnITSLayer4Pass",
						"",10,-1,1);
  fQAProtonsAcceptedList->Add(gProtonsPointOnITSLayer4Pass);
  TH1F *gProtonsPointOnITSLayer5Pass = new TH1F("gProtonsPointOnITSLayer5Pass",
						"",10,-1,1);
  fQAProtonsAcceptedList->Add(gProtonsPointOnITSLayer5Pass);
  TH1F *gProtonsPointOnITSLayer6Pass = new TH1F("gProtonsPointOnITSLayer6Pass",
						"",10,-1,1);
  fQAProtonsAcceptedList->Add(gProtonsPointOnITSLayer6Pass);
  TH1F *gProtonsNumberOfTPCdEdxPointsPass = new TH1F("gProtonsNumberOfTPCdEdxPointsPass","",100,0,200);
  fQAProtonsAcceptedList->Add(gProtonsNumberOfTPCdEdxPointsPass);

  //========================================================//  
  fQAProtonsRejectedList = new TList();
  fQAProtonsRejectedList->SetName("fQAProtonsRejectedList");
  fGlobalQAList->Add(fQAProtonsRejectedList);
  //Rejected protons
  TH1F *gProtonsITSClustersReject = new TH1F("gProtonsITSClustersReject",
					     ";N_{clusters} (ITS);Entries",
					     7,0,7);
  gProtonsITSClustersReject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsITSClustersReject);
  TH1F *gProtonsChi2PerClusterITSReject = new TH1F("gProtonsChi2PerClusterITSReject",
						   ";x^{2}/N_{clusters} (ITS);Entries",
						   100,0,4);
  gProtonsChi2PerClusterITSReject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsChi2PerClusterITSReject);
  TH1F *gProtonsTPCClustersReject = new TH1F("gProtonsTPCClustersReject",
					     ";N_{clusters} (TPC);Entries",
					     100,0,200);
  gProtonsTPCClustersReject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsTPCClustersReject);
  TH1F *gProtonsChi2PerClusterTPCReject = new TH1F("gProtonsChi2PerClusterTPCReject",
						   ";x^{2}/N_{clusters} (TPC);Entries",
						   100,0,4);
  gProtonsChi2PerClusterTPCReject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsChi2PerClusterTPCReject);
  TH1F *gProtonsExtCov11Reject = new TH1F("gProtonsExtCov11Reject",
					  ";#sigma_{y} [cm];Entries",
					  100,0,4);
  gProtonsExtCov11Reject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsExtCov11Reject);
  TH1F *gProtonsExtCov22Reject = new TH1F("gProtonsExtCov22Reject",
					  ";#sigma_{z} [cm];Entries",
					  100,0,4);
  gProtonsExtCov22Reject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsExtCov22Reject);
  TH1F *gProtonsExtCov33Reject = new TH1F("gProtonsExtCov33Reject",
					  ";#sigma_{sin(#phi)};Entries",
					  100,0,4);
  gProtonsExtCov33Reject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsExtCov33Reject);
  TH1F *gProtonsExtCov44Reject = new TH1F("gProtonsExtCov44Reject",
					  ";#sigma_{tan(#lambda)};Entries",
					  100,0,4);
  gProtonsExtCov44Reject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsExtCov44Reject);
  TH1F *gProtonsExtCov55Reject = new TH1F("gProtonsExtCov55Reject",
					  ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
					  100,0,4);
  gProtonsExtCov55Reject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsExtCov55Reject);
  TH1F *gProtonsSigmaToVertexReject = new TH1F("gProtonsSigmaToVertexReject",
					       ";#sigma_{Vertex};Entries",
					       100,0,10);
  gProtonsSigmaToVertexReject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsSigmaToVertexReject);
  TH1F *gProtonsSigmaToVertexTPCReject = new TH1F("gProtonsSigmaToVertexTPCReject",
						  ";#sigma_{Vertex};Entries",
						  100,0,10);
  gProtonsSigmaToVertexTPCReject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsSigmaToVertexTPCReject);
  TH1F *gProtonsDCAXYReject = new TH1F("gProtonsDCAXYReject",
				       ";DCA_{xy} [cm];Entries",
				       100,0,20);
  gProtonsDCAXYReject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsDCAXYReject);
  TH1F *gProtonsDCAXYTPCReject = new TH1F("gProtonsDCAXYTPCReject",
					  ";DCA_{xy} [cm];Entries",
					  100,0,20);
  gProtonsDCAXYTPCReject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsDCAXYTPCReject);
  TH1F *gProtonsDCAZReject = new TH1F("gProtonsDCAZReject",
				      ";DCA_{z} [cm];Entries",
				      100,0,20);
  gProtonsDCAZReject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsDCAZReject);
  TH1F *gProtonsDCAZTPCReject = new TH1F("gProtonsDCAZTPCReject",
					 ";DCA_{z} [cm];Entries",
					 100,0,20);
  gProtonsDCAZTPCReject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsDCAZTPCReject);
  TH1F *gProtonsConstrainChi2Reject = new TH1F("gProtonsConstrainChi2Reject",
					       ";Log_{10}(#chi^{2});Entries",
					       100,-10,10);
  gProtonsConstrainChi2Reject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsConstrainChi2Reject);
  TH1F *gProtonsITSRefitReject = new TH1F("gProtonsITSRefitReject",
					  "",10,-1,1);
  gProtonsITSRefitReject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsITSRefitReject);
  TH1F *gProtonsTPCRefitReject = new TH1F("gProtonsTPCRefitReject",
					  "",10,-1,1);
  gProtonsTPCRefitReject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsTPCRefitReject);
  TH1F *gProtonsESDpidReject = new TH1F("gProtonsESDpidReject",
					"",10,-1,1);
  gProtonsESDpidReject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsESDpidReject);
  TH1F *gProtonsTPCpidReject = new TH1F("gProtonsTPCpidReject",
					"",10,-1,1);
  gProtonsTPCpidReject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsTPCpidReject);
  TH1F *gProtonsPointOnITSLayer1Reject = new TH1F("gProtonsPointOnITSLayer1Reject",
						  "",10,-1,1);
  gProtonsPointOnITSLayer1Reject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsPointOnITSLayer1Reject);
  TH1F *gProtonsPointOnITSLayer2Reject = new TH1F("gProtonsPointOnITSLayer2Reject",
						  "",10,-1,1);
  gProtonsPointOnITSLayer2Reject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsPointOnITSLayer2Reject);
  TH1F *gProtonsPointOnITSLayer3Reject = new TH1F("gProtonsPointOnITSLayer3Reject",
						  "",10,-1,1);
  gProtonsPointOnITSLayer3Reject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsPointOnITSLayer3Reject);
  TH1F *gProtonsPointOnITSLayer4Reject = new TH1F("gProtonsPointOnITSLayer4Reject",
						  "",10,-1,1);
  gProtonsPointOnITSLayer4Reject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsPointOnITSLayer4Reject);
  TH1F *gProtonsPointOnITSLayer5Reject = new TH1F("gProtonsPointOnITSLayer5Reject",
						  "",10,-1,1);
  gProtonsPointOnITSLayer5Reject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsPointOnITSLayer5Reject);
  TH1F *gProtonsPointOnITSLayer6Reject = new TH1F("gProtonsPointOnITSLayer6Reject",
						  "",10,-1,1);
  gProtonsPointOnITSLayer6Reject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsPointOnITSLayer6Reject);
  TH1F *gProtonsNumberOfTPCdEdxPointsReject = new TH1F("gProtonsNumberOfTPCdEdxPointsReject","",100,0,200);
  gProtonsNumberOfTPCdEdxPointsReject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsNumberOfTPCdEdxPointsReject);
    
  //========================================================//
  fQAAntiProtonsAcceptedList = new TList();
  fQAAntiProtonsAcceptedList->SetName("fQAAntiProtonsAcceptedList");
  fGlobalQAList->Add(fQAAntiProtonsAcceptedList);
  //Accepted antiprotons
  TH1F *gAntiProtonsITSClustersPass = new TH1F("gAntiProtonsITSClustersPass",
					       ";N_{clusters} (ITS);Entries",
					       7,0,7);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsITSClustersPass);
  TH1F *gAntiProtonsChi2PerClusterITSPass = new TH1F("gAntiProtonsChi2PerClusterITSPass",
						     ";x^{2}/N_{clusters} (ITS);Entries",
						     100,0,4);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsChi2PerClusterITSPass);
  TH1F *gAntiProtonsTPCClustersPass = new TH1F("gAntiProtonsTPCClustersPass",
					       ";N_{clusters} (TPC);Entries",
					       100,0,200);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsTPCClustersPass);
  TH1F *gAntiProtonsChi2PerClusterTPCPass = new TH1F("gAntiProtonsChi2PerClusterTPCPass",
						     ";x^{2}/N_{clusters} (TPC);Entries",
						     100,0,4);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsChi2PerClusterTPCPass);
  TH1F *gAntiProtonsExtCov11Pass = new TH1F("gAntiProtonsExtCov11Pass",
					    ";#sigma_{y} [cm];Entries",
					    100,0,4);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsExtCov11Pass);
  TH1F *gAntiProtonsExtCov22Pass = new TH1F("gAntiProtonsExtCov22Pass",
					    ";#sigma_{z} [cm];Entries",
					    100,0,4);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsExtCov22Pass);
  TH1F *gAntiProtonsExtCov33Pass = new TH1F("gAntiProtonsExtCov33Pass",
					    ";#sigma_{sin(#phi)};Entries",
					    100,0,4);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsExtCov33Pass);
  TH1F *gAntiProtonsExtCov44Pass = new TH1F("gAntiProtonsExtCov44Pass",
					    ";#sigma_{tan(#lambda)};Entries",
					    100,0,4);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsExtCov44Pass);
  TH1F *gAntiProtonsExtCov55Pass = new TH1F("gAntiProtonsExtCov55Pass",
					    ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
					    100,0,4);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsExtCov55Pass);
  TH1F *gAntiProtonsSigmaToVertexPass = new TH1F("gAntiProtonsSigmaToVertexPass",
						 ";#sigma_{Vertex};Entries",
						 100,0,10);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsSigmaToVertexPass);
  TH1F *gAntiProtonsSigmaToVertexTPCPass = new TH1F("gAntiProtonsSigmaToVertexTPCPass",
						    ";#sigma_{Vertex};Entries",
						    100,0,10);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsSigmaToVertexTPCPass);
  TH1F *gAntiProtonsDCAXYPass = new TH1F("gAntiProtonsDCAXYPass",
					 ";DCA_{xy} [cm];Entries",
					 100,0,20);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsDCAXYPass);
  TH1F *gAntiProtonsDCAXYTPCPass = new TH1F("gAntiProtonsDCAXYTPCPass",
					    ";DCA_{xy} [cm];Entries",
					    100,0,20);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsDCAXYTPCPass);
  TH1F *gAntiProtonsDCAZPass = new TH1F("gAntiProtonsDCAZPass",
					";DCA_{z} [cm];Entries",
					100,0,20);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsDCAZPass);
  TH1F *gAntiProtonsDCAZTPCPass = new TH1F("gAntiProtonsDCAZTPCPass",
					   ";DCA_{z} [cm];Entries",
					   100,0,20);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsDCAZTPCPass);
  TH1F *gAntiProtonsConstrainChi2Pass = new TH1F("gAntiProtonsConstrainChi2Pass",
						 ";Log_{10}(#chi^{2});Entries",
						 100,-10,10);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsConstrainChi2Pass);
  TH1F *gAntiProtonsITSRefitPass = new TH1F("gAntiProtonsITSRefitPass",
					    "",10,-1,1);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsITSRefitPass);
  TH1F *gAntiProtonsTPCRefitPass = new TH1F("gAntiProtonsTPCRefitPass",
					    "",10,-1,1);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsTPCRefitPass);
  TH1F *gAntiProtonsESDpidPass = new TH1F("gAntiProtonsESDpidPass",
					  "",10,-1,1);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsESDpidPass);
  TH1F *gAntiProtonsTPCpidPass = new TH1F("gAntiProtonsTPCpidPass",
					  "",10,-1,1);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsTPCpidPass);
  TH1F *gAntiProtonsPointOnITSLayer1Pass = new TH1F("gAntiProtonsPointOnITSLayer1Pass",
						    "",10,-1,1);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsPointOnITSLayer1Pass);
  TH1F *gAntiProtonsPointOnITSLayer2Pass = new TH1F("gAntiProtonsPointOnITSLayer2Pass",
						    "",10,-1,1);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsPointOnITSLayer2Pass);
  TH1F *gAntiProtonsPointOnITSLayer3Pass = new TH1F("gAntiProtonsPointOnITSLayer3Pass",
						    "",10,-1,1);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsPointOnITSLayer3Pass);
  TH1F *gAntiProtonsPointOnITSLayer4Pass = new TH1F("gAntiProtonsPointOnITSLayer4Pass",
						    "",10,-1,1);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsPointOnITSLayer4Pass);
  TH1F *gAntiProtonsPointOnITSLayer5Pass = new TH1F("gAntiProtonsPointOnITSLayer5Pass",
						    "",10,-1,1);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsPointOnITSLayer5Pass);
  TH1F *gAntiProtonsPointOnITSLayer6Pass = new TH1F("gAntiProtonsPointOnITSLayer6Pass",
						    "",10,-1,1);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsPointOnITSLayer6Pass);
  TH1F *gAntiProtonsNumberOfTPCdEdxPointsPass = new TH1F("gAntiProtonsNumberOfTPCdEdxPointsPass","",100,0,200);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsNumberOfTPCdEdxPointsPass);
  
  //========================================================//
  fQAAntiProtonsRejectedList = new TList();
  fQAAntiProtonsRejectedList->SetName("fQAAntiProtonsRejectedList");
  fGlobalQAList->Add(fQAAntiProtonsRejectedList);
  //Rejected antiprotons
  TH1F *gAntiProtonsITSClustersReject = new TH1F("gAntiProtonsITSClustersReject",
						 ";N_{clusters} (ITS);Entries",
						 7,0,7);
  gAntiProtonsITSClustersReject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsITSClustersReject);
  TH1F *gAntiProtonsChi2PerClusterITSReject = new TH1F("gAntiProtonsChi2PerClusterITSReject",
						       ";x^{2}/N_{clusters} (ITS);Entries",
						       100,0,4);
  gAntiProtonsChi2PerClusterITSReject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsChi2PerClusterITSReject);
  TH1F *gAntiProtonsTPCClustersReject = new TH1F("gAntiProtonsTPCClustersReject",
						 ";N_{clusters} (TPC);Entries",
						 100,0,200);
  gAntiProtonsTPCClustersReject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsTPCClustersReject);
  TH1F *gAntiProtonsChi2PerClusterTPCReject = new TH1F("gAntiProtonsChi2PerClusterTPCReject",
						       ";x^{2}/N_{clusters} (TPC);Entries",
						       100,0,4);
  gAntiProtonsChi2PerClusterTPCReject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsChi2PerClusterTPCReject);
  TH1F *gAntiProtonsExtCov11Reject = new TH1F("gAntiProtonsExtCov11Reject",
					      ";#sigma_{y} [cm];Entries",
					      100,0,4);
  gAntiProtonsExtCov11Reject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsExtCov11Reject);
  TH1F *gAntiProtonsExtCov22Reject = new TH1F("gAntiProtonsExtCov22Reject",
					      ";#sigma_{z} [cm];Entries",
					      100,0,4);
  gAntiProtonsExtCov22Reject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsExtCov22Reject);
  TH1F *gAntiProtonsExtCov33Reject = new TH1F("gAntiProtonsExtCov33Reject",
					      ";#sigma_{sin(#phi)};Entries",
					      100,0,4);
  gAntiProtonsExtCov33Reject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsExtCov33Reject);
  TH1F *gAntiProtonsExtCov44Reject = new TH1F("gAntiProtonsExtCov44Reject",
					      ";#sigma_{tan(#lambda)};Entries",
					      100,0,4);
  gAntiProtonsExtCov44Reject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsExtCov44Reject);
  TH1F *gAntiProtonsExtCov55Reject = new TH1F("gAntiProtonsExtCov55Reject",
					      ";#sigma_{1/P_{T}} [GeV/c]^{-1};Entries",
					      100,0,4);
  gAntiProtonsExtCov55Reject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsExtCov55Reject);
  TH1F *gAntiProtonsSigmaToVertexReject = new TH1F("gAntiProtonsSigmaToVertexReject",
						   ";#sigma_{Vertex};Entries",
						   100,0,10);
  gAntiProtonsSigmaToVertexReject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsSigmaToVertexReject);
  TH1F *gAntiProtonsSigmaToVertexTPCReject = new TH1F("gAntiProtonsSigmaToVertexTPCReject",
						      ";#sigma_{Vertex};Entries",
						      100,0,10);
  gAntiProtonsSigmaToVertexTPCReject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsSigmaToVertexTPCReject);
  TH1F *gAntiProtonsDCAXYReject = new TH1F("gAntiProtonsDCAXYReject",
					   ";DCA_{xy} [cm];Entries",
					   100,0,20);
  gAntiProtonsDCAXYReject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsDCAXYReject);
  TH1F *gAntiProtonsDCAXYTPCReject = new TH1F("gAntiProtonsDCAXYTPCReject",
					      ";DCA_{xy} [cm];Entries",
					      100,0,20);
  gAntiProtonsDCAXYTPCReject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsDCAXYTPCReject);
  TH1F *gAntiProtonsDCAZReject = new TH1F("gAntiProtonsDCAZReject",
					  ";DCA_{z} [cm];Entries",
					  100,0,20);
  gAntiProtonsDCAZReject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsDCAZReject);
  TH1F *gAntiProtonsDCAZTPCReject = new TH1F("gAntiProtonsDCAZTPCReject",
					     ";DCA_{z} [cm];Entries",
					     100,0,20);
  gAntiProtonsDCAZTPCReject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsDCAZTPCReject);
  TH1F *gAntiProtonsConstrainChi2Reject = new TH1F("gAntiProtonsConstrainChi2Reject",
						   ";Log_{10}(#chi^{2});Entries",
						   100,-10,10);
  gAntiProtonsConstrainChi2Reject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsConstrainChi2Reject);
  TH1F *gAntiProtonsITSRefitReject = new TH1F("gAntiProtonsITSRefitReject",
					      "",10,-1,1);
  gAntiProtonsITSRefitReject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsITSRefitReject);
  TH1F *gAntiProtonsTPCRefitReject = new TH1F("gAntiProtonsTPCRefitReject",
					      "",10,-1,1);
  gAntiProtonsTPCRefitReject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsTPCRefitReject);
  TH1F *gAntiProtonsESDpidReject = new TH1F("gAntiProtonsESDpidReject",
					    "",10,-1,1);
  gAntiProtonsESDpidReject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsESDpidReject);
  TH1F *gAntiProtonsTPCpidReject = new TH1F("gAntiProtonsTPCpidReject",
					    "",10,-1,1);
  gAntiProtonsTPCpidReject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsTPCpidReject);
  TH1F *gAntiProtonsPointOnITSLayer1Reject = new TH1F("gAntiProtonsPointOnITSLayer1Reject",
						      "",10,-1,1);
  gAntiProtonsPointOnITSLayer1Reject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsPointOnITSLayer1Reject);
  TH1F *gAntiProtonsPointOnITSLayer2Reject = new TH1F("gAntiProtonsPointOnITSLayer2Reject",
						      "",10,-1,1);
  gAntiProtonsPointOnITSLayer2Reject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsPointOnITSLayer2Reject);
  TH1F *gAntiProtonsPointOnITSLayer3Reject = new TH1F("gAntiProtonsPointOnITSLayer3Reject",
						      "",10,-1,1);
  gAntiProtonsPointOnITSLayer3Reject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsPointOnITSLayer3Reject);
  TH1F *gAntiProtonsPointOnITSLayer4Reject = new TH1F("gAntiProtonsPointOnITSLayer4Reject",
						      "",10,-1,1);
  gAntiProtonsPointOnITSLayer4Reject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsPointOnITSLayer4Reject);
  TH1F *gAntiProtonsPointOnITSLayer5Reject = new TH1F("gAntiProtonsPointOnITSLayer5Reject",
						      "",10,-1,1);
  gAntiProtonsPointOnITSLayer5Reject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsPointOnITSLayer5Reject);
  TH1F *gAntiProtonsPointOnITSLayer6Reject = new TH1F("gAntiProtonsPointOnITSLayer6Reject",
						      "",10,-1,1);
  gAntiProtonsPointOnITSLayer6Reject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsPointOnITSLayer6Reject);
  TH1F *gAntiProtonsNumberOfTPCdEdxPointsReject = new TH1F("gAntiProtonsNumberOfTPCdEdxPointsReject","",100,0,200);
  gAntiProtonsNumberOfTPCdEdxPointsReject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsNumberOfTPCdEdxPointsReject);
}

//____________________________________________________________________//
void AliProtonAnalysis::FillQA(AliESDEvent *esd,
			       const AliESDVertex *vertex, 
			       AliESDtrack* track) {
  //Fills the QA histograms
  Double_t gPt = 0.0, gPx = 0.0, gPy = 0.0, gPz = 0.0;
  Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};  //The impact parameters and their covariance.

  if((fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kTPC)||(fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kHybrid)) {
    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(!tpcTrack) {
      gPt = 0.0; gPx = 0.0; gPy = 0.0; gPz = 0.0;
      dca[0] = -100.; dca[1] = -100.;
      cov[0] = -100.; cov[1] = -100.; cov[2] = -100.;
    }
    else {
      gPt = tpcTrack->Pt();
      gPx = tpcTrack->Px();
      gPy = tpcTrack->Py();
      gPz = tpcTrack->Pz();
      tpcTrack->PropagateToDCA(vertex,
			       esd->GetMagneticField(),
			       100.,dca,cov);
    }
  }
  else{
    gPt = track->Pt();
    gPx = track->Px();
    gPy = track->Py();
    gPz = track->Pz();
    track->PropagateToDCA(vertex,
			  esd->GetMagneticField(),
			  100.,dca,cov);
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
    if(fProtonAnalysisBase->IsUsedMinITSClusters()) {
      if(nClustersITS < fProtonAnalysisBase->GetMinITSClusters()) {
	((TH1F *)(fQAProtonsRejectedList->At(0)))->Fill(nClustersITS);
      }
      else if(nClustersITS >= fProtonAnalysisBase->GetMinITSClusters()) 
	((TH1F *)(fQAProtonsAcceptedList->At(0)))->Fill(nClustersITS);
    }//ITS clusters
    if(fProtonAnalysisBase->IsUsedMaxChi2PerITSCluster()) {
      if(chi2PerClusterITS > fProtonAnalysisBase->GetMaxChi2PerITSCluster()) {
	((TH1F *)(fQAProtonsRejectedList->At(1)))->Fill(chi2PerClusterITS);
      }
      else if(chi2PerClusterITS <= fProtonAnalysisBase->GetMaxChi2PerITSCluster())
	((TH1F *)(fQAProtonsAcceptedList->At(1)))->Fill(chi2PerClusterITS);
    }//chi2 per ITS cluster
    if(fProtonAnalysisBase->IsUsedMinTPCClusters()) {
      if(nClustersTPC < fProtonAnalysisBase->GetMinTPCClusters()) {
	((TH1F *)(fQAProtonsRejectedList->At(2)))->Fill(nClustersTPC);
      }
      else if(nClustersTPC >= fProtonAnalysisBase->GetMinTPCClusters()) {
	((TH1F *)(fQAProtonsAcceptedList->At(2)))->Fill(nClustersTPC);
      }
    }//TPC clusters
    if(fProtonAnalysisBase->IsUsedMaxChi2PerTPCCluster()) {
      if(chi2PerClusterTPC > fProtonAnalysisBase->GetMaxChi2PerTPCCluster()) {
	((TH1F *)(fQAProtonsRejectedList->At(3)))->Fill(chi2PerClusterTPC);
      }
      else if(chi2PerClusterTPC <= fProtonAnalysisBase->GetMaxChi2PerTPCCluster())
	((TH1F *)(fQAProtonsAcceptedList->At(3)))->Fill(chi2PerClusterTPC);
    }//chi2 per TPC cluster
    if(fProtonAnalysisBase->IsUsedMaxCov11()) {
      if(extCov[0] > fProtonAnalysisBase->GetMaxCov11()) {
	((TH1F *)(fQAProtonsRejectedList->At(4)))->Fill(extCov[0]);
      }
      else if(extCov[0] <= fProtonAnalysisBase->GetMaxCov11())
	((TH1F *)(fQAProtonsAcceptedList->At(4)))->Fill(extCov[0]);
    }//cov11
    if(fProtonAnalysisBase->IsUsedMaxCov22()) {
      if(extCov[2] > fProtonAnalysisBase->GetMaxCov22()) {
	((TH1F *)(fQAProtonsRejectedList->At(5)))->Fill(extCov[2]);
      }
      else if(extCov[2] <= fProtonAnalysisBase->GetMaxCov22())
	((TH1F *)(fQAProtonsAcceptedList->At(5)))->Fill(extCov[2]);
    }//cov11
    if(fProtonAnalysisBase->IsUsedMaxCov33()) {
      if(extCov[5] > fProtonAnalysisBase->GetMaxCov33()) {
	((TH1F *)(fQAProtonsRejectedList->At(6)))->Fill(extCov[5]);
      }
      else if(extCov[5] <= fProtonAnalysisBase->GetMaxCov33())
	((TH1F *)(fQAProtonsAcceptedList->At(6)))->Fill(extCov[5]);
    }//cov11
    if(fProtonAnalysisBase->IsUsedMaxCov44()) {
      if(extCov[9] > fProtonAnalysisBase->GetMaxCov44()) {
	((TH1F *)(fQAProtonsRejectedList->At(7)))->Fill(extCov[9]);
      }
      else if(extCov[9] <= fProtonAnalysisBase->GetMaxCov44())
	((TH1F *)(fQAProtonsAcceptedList->At(7)))->Fill(extCov[9]);
    }//cov11
    if(fProtonAnalysisBase->IsUsedMaxCov55()) {
      if(extCov[14] > fProtonAnalysisBase->GetMaxCov55()) {
	((TH1F *)(fQAProtonsRejectedList->At(8)))->Fill(extCov[14]);
      }
      else if(extCov[14] <= fProtonAnalysisBase->GetMaxCov55())
	((TH1F *)(fQAProtonsAcceptedList->At(8)))->Fill(extCov[14]);
    }//cov55
    if(fProtonAnalysisBase->IsUsedMaxSigmaToVertex()) {
      if(fProtonAnalysisBase->GetSigmaToVertex(track) > fProtonAnalysisBase->GetMaxSigmaToVertex()) {
	((TH1F *)(fQAProtonsRejectedList->At(9)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
      }
      else if(fProtonAnalysisBase->GetSigmaToVertex(track) <= fProtonAnalysisBase->GetMaxSigmaToVertex())
	((TH1F *)(fQAProtonsAcceptedList->At(9)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
    }//sigma to vertex
    if(fProtonAnalysisBase->IsUsedMaxSigmaToVertexTPC()) {
      if(fProtonAnalysisBase->GetSigmaToVertex(track) > fProtonAnalysisBase->GetMaxSigmaToVertexTPC()) {
	((TH1F *)(fQAProtonsRejectedList->At(10)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
      }
      else if(fProtonAnalysisBase->GetSigmaToVertex(track) <= fProtonAnalysisBase->GetMaxSigmaToVertexTPC())
	((TH1F *)(fQAProtonsAcceptedList->At(10)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
    }//sigma to vertex TPC
    if(fProtonAnalysisBase->IsUsedMaxDCAXY()) {
      if(TMath::Abs(dca[0]) > fProtonAnalysisBase->GetMaxDCAXY()) {
	((TH1F *)(fQAProtonsRejectedList->At(11)))->Fill(TMath::Abs(dca[0]));
      }
      else if(TMath::Abs(dca[0]) <= fProtonAnalysisBase->GetMaxDCAXY())
	((TH1F *)(fQAProtonsAcceptedList->At(11)))->Fill(TMath::Abs(dca[0]));
    }//DCA xy global tracking
    if(fProtonAnalysisBase->IsUsedMaxDCAXYTPC()) {
      if(TMath::Abs(dca[0]) > fProtonAnalysisBase->GetMaxDCAXYTPC()) {
	((TH1F *)(fQAProtonsRejectedList->At(12)))->Fill(TMath::Abs(dca[0]));
      }
      else if(TMath::Abs(dca[0]) <= fProtonAnalysisBase->GetMaxDCAXYTPC())
	((TH1F *)(fQAProtonsAcceptedList->At(12)))->Fill(TMath::Abs(dca[0]));
    }//DCA xy TPC tracking
    if(fProtonAnalysisBase->IsUsedMaxDCAZ()) {
      if(TMath::Abs(dca[1]) > fProtonAnalysisBase->GetMaxDCAZ()) {
	((TH1F *)(fQAProtonsRejectedList->At(13)))->Fill(TMath::Abs(dca[1]));
      }
      else if(TMath::Abs(dca[1]) <= fProtonAnalysisBase->GetMaxDCAZ())
	((TH1F *)(fQAProtonsAcceptedList->At(13)))->Fill(TMath::Abs(dca[1]));
    }//DCA z global tracking
    if(fProtonAnalysisBase->IsUsedMaxDCAZTPC()) {
      if(TMath::Abs(dca[1]) > fProtonAnalysisBase->GetMaxDCAZTPC()) {
	((TH1F *)(fQAProtonsRejectedList->At(14)))->Fill(TMath::Abs(dca[1]));
      }
      else if(TMath::Abs(dca[1]) <= fProtonAnalysisBase->GetMaxDCAZTPC())
	((TH1F *)(fQAProtonsAcceptedList->At(14)))->Fill(TMath::Abs(dca[1]));
    }//DCA z TPC tracking
    if(fProtonAnalysisBase->IsUsedMaxConstrainChi2()) {
      if(track->GetConstrainedChi2() > 0) {
	if(TMath::Log(track->GetConstrainedChi2()) > fProtonAnalysisBase->GetMaxConstrainChi2()) {
	  ((TH1F *)(fQAProtonsRejectedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
	}
	else if(TMath::Log(track->GetConstrainedChi2()) <= fProtonAnalysisBase->GetMaxConstrainChi2())
	  ((TH1F *)(fQAProtonsAcceptedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
      }
    }//constrain chi2 - vertex
    if(fProtonAnalysisBase->IsUsedITSRefit()) {
      if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) {
	((TH1F *)(fQAProtonsRejectedList->At(16)))->Fill(0);
      }
      else if((track->GetStatus() & AliESDtrack::kITSrefit) != 0)
	((TH1F *)(fQAProtonsAcceptedList->At(16)))->Fill(0);
    }//ITS refit
    if(fProtonAnalysisBase->IsUsedTPCRefit()) {
      if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) {
	((TH1F *)(fQAProtonsRejectedList->At(17)))->Fill(0);
      }
      else if((track->GetStatus() & AliESDtrack::kTPCrefit) != 0)
	((TH1F *)(fQAProtonsAcceptedList->At(17)))->Fill(0);
    }//TPC refit
    if(fProtonAnalysisBase->IsUsedESDpid()) {
      if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) {
	((TH1F *)(fQAProtonsRejectedList->At(18)))->Fill(0);
      }
      else if((track->GetStatus() & AliESDtrack::kESDpid) != 0)
	((TH1F *)(fQAProtonsAcceptedList->At(18)))->Fill(0);
    }//ESD pid
    if(fProtonAnalysisBase->IsUsedTPCpid()) {
      if ((track->GetStatus() & AliESDtrack::kTPCpid) == 0) {
	((TH1F *)(fQAProtonsRejectedList->At(19)))->Fill(0);
      }
      else if((track->GetStatus() & AliESDtrack::kTPCpid) != 0)
	((TH1F *)(fQAProtonsAcceptedList->At(19)))->Fill(0);
    }//TPC pid
    if(fProtonAnalysisBase->IsUsedPointOnITSLayer1()) {
      if(!track->HasPointOnITSLayer(0)) {
	((TH1F *)(fQAProtonsRejectedList->At(20)))->Fill(0);
      }
      else if(track->HasPointOnITSLayer(0))
	((TH1F *)(fQAProtonsAcceptedList->At(20)))->Fill(0);
    }//point on SPD1
    if(fProtonAnalysisBase->IsUsedPointOnITSLayer2()) {
      if(!track->HasPointOnITSLayer(1)) {
	((TH1F *)(fQAProtonsRejectedList->At(21)))->Fill(0);
      }
      else if(track->HasPointOnITSLayer(1))
	((TH1F *)(fQAProtonsAcceptedList->At(21)))->Fill(0);
    }//point on SPD2
    if(fProtonAnalysisBase->IsUsedPointOnITSLayer3()) {
      if(!track->HasPointOnITSLayer(2)) {
	((TH1F *)(fQAProtonsRejectedList->At(22)))->Fill(0);
      }
      else if(track->HasPointOnITSLayer(2))
	((TH1F *)(fQAProtonsAcceptedList->At(22)))->Fill(0);
    }//point on SDD1
    if(fProtonAnalysisBase->IsUsedPointOnITSLayer4()) {
      if(!track->HasPointOnITSLayer(3)) {
	((TH1F *)(fQAProtonsRejectedList->At(23)))->Fill(0);
      }
      else if(track->HasPointOnITSLayer(3))
	((TH1F *)(fQAProtonsAcceptedList->At(23)))->Fill(0);
    }//point on SDD2
    if(fProtonAnalysisBase->IsUsedPointOnITSLayer5()) {
      if(!track->HasPointOnITSLayer(4)) {
	((TH1F *)(fQAProtonsRejectedList->At(24)))->Fill(0);
      }
      else if(track->HasPointOnITSLayer(4))
	((TH1F *)(fQAProtonsAcceptedList->At(24)))->Fill(0);
    }//point on SSD1
    if(fProtonAnalysisBase->IsUsedPointOnITSLayer6()) {
      if(!track->HasPointOnITSLayer(5)) {
	((TH1F *)(fQAProtonsRejectedList->At(25)))->Fill(0);
      }
      else if(track->HasPointOnITSLayer(5))
	((TH1F *)(fQAProtonsAcceptedList->At(25)))->Fill(0);
    }//point on SSD2
    if(fProtonAnalysisBase->IsUsedMinTPCdEdxPoints()) {
      if(track->GetTPCsignalN() < fProtonAnalysisBase->GetMinTPCdEdxPoints()) {
	((TH1F *)(fQAProtonsRejectedList->At(26)))->Fill(track->GetTPCsignalN());
      }
      if(track->GetTPCsignalN() >= fProtonAnalysisBase->GetMinTPCdEdxPoints())
	((TH1F *)(fQAProtonsAcceptedList->At(26)))->Fill(track->GetTPCsignalN());
    }//number of TPC points for the dE/dx
  }//protons

  //antiprotons
  if(track->Charge() < 0) {
    if(fProtonAnalysisBase->IsUsedMinITSClusters()) {
      if(nClustersITS < fProtonAnalysisBase->GetMinITSClusters()) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(0)))->Fill(nClustersITS);
      }
      else if(nClustersITS >= fProtonAnalysisBase->GetMinITSClusters()) 
	((TH1F *)(fQAAntiProtonsAcceptedList->At(0)))->Fill(nClustersITS);
    }//ITS clusters
    if(fProtonAnalysisBase->IsUsedMaxChi2PerITSCluster()) {
      if(chi2PerClusterITS > fProtonAnalysisBase->GetMaxChi2PerITSCluster()) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(1)))->Fill(chi2PerClusterITS);
      }
      else if(chi2PerClusterITS <= fProtonAnalysisBase->GetMaxChi2PerITSCluster())
	((TH1F *)(fQAAntiProtonsAcceptedList->At(1)))->Fill(chi2PerClusterITS);
    }//chi2 per ITS cluster
    if(fProtonAnalysisBase->IsUsedMinTPCClusters()) {
      if(nClustersTPC < fProtonAnalysisBase->GetMinTPCClusters()) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(2)))->Fill(nClustersTPC);
      }
      else if(nClustersTPC >= fProtonAnalysisBase->GetMinTPCClusters()) {
	((TH1F *)(fQAAntiProtonsAcceptedList->At(2)))->Fill(nClustersTPC);
      }
    }//TPC clusters
    if(fProtonAnalysisBase->IsUsedMaxChi2PerTPCCluster()) {
      if(chi2PerClusterTPC > fProtonAnalysisBase->GetMaxChi2PerTPCCluster()) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(3)))->Fill(chi2PerClusterTPC);
      }
      else if(chi2PerClusterTPC <= fProtonAnalysisBase->GetMaxChi2PerTPCCluster())
	((TH1F *)(fQAAntiProtonsAcceptedList->At(3)))->Fill(chi2PerClusterTPC);
    }//chi2 per TPC cluster
    if(fProtonAnalysisBase->IsUsedMaxCov11()) {
      if(extCov[0] > fProtonAnalysisBase->GetMaxCov11()) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(4)))->Fill(extCov[0]);
      }
      else if(extCov[0] <= fProtonAnalysisBase->GetMaxCov11())
	((TH1F *)(fQAAntiProtonsAcceptedList->At(4)))->Fill(extCov[0]);
    }//cov11
    if(fProtonAnalysisBase->IsUsedMaxCov22()) {
      if(extCov[2] > fProtonAnalysisBase->GetMaxCov22()) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(5)))->Fill(extCov[2]);
      }
      else if(extCov[2] <= fProtonAnalysisBase->GetMaxCov22())
	((TH1F *)(fQAAntiProtonsAcceptedList->At(5)))->Fill(extCov[2]);
    }//cov11
    if(fProtonAnalysisBase->IsUsedMaxCov33()) {
      if(extCov[5] > fProtonAnalysisBase->GetMaxCov33()) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(6)))->Fill(extCov[5]);
      }
      else if(extCov[5] <= fProtonAnalysisBase->GetMaxCov33())
	((TH1F *)(fQAAntiProtonsAcceptedList->At(6)))->Fill(extCov[5]);
    }//cov11
    if(fProtonAnalysisBase->IsUsedMaxCov44()) {
      if(extCov[9] > fProtonAnalysisBase->GetMaxCov44()) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(7)))->Fill(extCov[9]);
      }
      else if(extCov[9] <= fProtonAnalysisBase->GetMaxCov44())
	((TH1F *)(fQAAntiProtonsAcceptedList->At(7)))->Fill(extCov[9]);
    }//cov11
    if(fProtonAnalysisBase->IsUsedMaxCov55()) {
      if(extCov[14] > fProtonAnalysisBase->GetMaxCov55()) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(8)))->Fill(extCov[14]);
      }
      else if(extCov[14] <= fProtonAnalysisBase->GetMaxCov55())
	((TH1F *)(fQAAntiProtonsAcceptedList->At(8)))->Fill(extCov[14]);
    }//cov55
    if(fProtonAnalysisBase->IsUsedMaxSigmaToVertex()) {
      if(fProtonAnalysisBase->GetSigmaToVertex(track) > fProtonAnalysisBase->GetMaxSigmaToVertex()) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(9)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
      }
      else if(fProtonAnalysisBase->GetSigmaToVertex(track) <= fProtonAnalysisBase->GetMaxSigmaToVertex())
	((TH1F *)(fQAAntiProtonsAcceptedList->At(9)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
    }//sigma to vertex
    if(fProtonAnalysisBase->IsUsedMaxSigmaToVertexTPC()) {
      if(fProtonAnalysisBase->GetSigmaToVertex(track) > fProtonAnalysisBase->GetMaxSigmaToVertexTPC()) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(10)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
      }
      else if(fProtonAnalysisBase->GetSigmaToVertex(track) <= fProtonAnalysisBase->GetMaxSigmaToVertexTPC())
	((TH1F *)(fQAAntiProtonsAcceptedList->At(10)))->Fill(fProtonAnalysisBase->GetSigmaToVertex(track));
    }//sigma to vertex TPC
    if(fProtonAnalysisBase->IsUsedMaxDCAXY()) {
      if(TMath::Abs(dca[0]) > fProtonAnalysisBase->GetMaxDCAXY()) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(11)))->Fill(TMath::Abs(dca[0]));
      }
      else if(TMath::Abs(dca[0]) <= fProtonAnalysisBase->GetMaxDCAXY())
	((TH1F *)(fQAAntiProtonsAcceptedList->At(11)))->Fill(TMath::Abs(dca[0]));
    }//DCA xy global tracking
    if(fProtonAnalysisBase->IsUsedMaxDCAXYTPC()) {
      if(TMath::Abs(dca[0]) > fProtonAnalysisBase->GetMaxDCAXYTPC()) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(12)))->Fill(TMath::Abs(dca[0]));
      }
      else if(TMath::Abs(dca[0]) <= fProtonAnalysisBase->GetMaxDCAXYTPC())
	((TH1F *)(fQAAntiProtonsAcceptedList->At(12)))->Fill(TMath::Abs(dca[0]));
    }//DCA xy TPC tracking
    if(fProtonAnalysisBase->IsUsedMaxDCAZ()) {
      if(TMath::Abs(dca[1]) > fProtonAnalysisBase->GetMaxDCAZ()) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(13)))->Fill(TMath::Abs(dca[1]));
      }
      else if(TMath::Abs(dca[1]) <= fProtonAnalysisBase->GetMaxDCAZ())
	((TH1F *)(fQAAntiProtonsAcceptedList->At(13)))->Fill(TMath::Abs(dca[1]));
    }//DCA z global tracking
    if(fProtonAnalysisBase->IsUsedMaxDCAZTPC()) {
      if(TMath::Abs(dca[1]) > fProtonAnalysisBase->GetMaxDCAZTPC()) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(14)))->Fill(TMath::Abs(dca[1]));
      }
      else if(TMath::Abs(dca[1]) <= fProtonAnalysisBase->GetMaxDCAZTPC())
	((TH1F *)(fQAAntiProtonsAcceptedList->At(14)))->Fill(TMath::Abs(dca[1]));
    }//DCA z TPC tracking
    if(fProtonAnalysisBase->IsUsedMaxConstrainChi2()) {
      if(track->GetConstrainedChi2() > 0) {
	if(TMath::Log(track->GetConstrainedChi2()) > fProtonAnalysisBase->GetMaxConstrainChi2()) {
	  ((TH1F *)(fQAAntiProtonsRejectedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
	}
	else if(TMath::Log(track->GetConstrainedChi2()) <= fProtonAnalysisBase->GetMaxConstrainChi2())
	  ((TH1F *)(fQAAntiProtonsAcceptedList->At(15)))->Fill(TMath::Log(track->GetConstrainedChi2()));
      }
    }//constrain chi2 - vertex
    if(fProtonAnalysisBase->IsUsedITSRefit()) {
      if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(16)))->Fill(0);
      }
      else if((track->GetStatus() & AliESDtrack::kITSrefit) != 0)
	((TH1F *)(fQAAntiProtonsAcceptedList->At(16)))->Fill(0);
    }//ITS refit
    if(fProtonAnalysisBase->IsUsedTPCRefit()) {
      if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(17)))->Fill(0);
      }
      else if((track->GetStatus() & AliESDtrack::kTPCrefit) != 0)
	((TH1F *)(fQAAntiProtonsAcceptedList->At(17)))->Fill(0);
    }//TPC refit
    if(fProtonAnalysisBase->IsUsedESDpid()) {
      if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(18)))->Fill(0);
      }
      else if((track->GetStatus() & AliESDtrack::kESDpid) != 0)
	((TH1F *)(fQAAntiProtonsAcceptedList->At(18)))->Fill(0);
    }//ESD pid
    if(fProtonAnalysisBase->IsUsedTPCpid()) {
      if ((track->GetStatus() & AliESDtrack::kTPCpid) == 0) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(19)))->Fill(0);
      }
      else if((track->GetStatus() & AliESDtrack::kTPCpid) != 0)
	((TH1F *)(fQAAntiProtonsAcceptedList->At(19)))->Fill(0);
    }//TPC pid
    if(fProtonAnalysisBase->IsUsedPointOnITSLayer1()) {
      if(!track->HasPointOnITSLayer(0)) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(20)))->Fill(0);
      }
      else if(track->HasPointOnITSLayer(0))
	((TH1F *)(fQAAntiProtonsAcceptedList->At(20)))->Fill(0);
    }//point on SPD1
    if(fProtonAnalysisBase->IsUsedPointOnITSLayer2()) {
      if(!track->HasPointOnITSLayer(1)) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(21)))->Fill(0);
      }
      else if(track->HasPointOnITSLayer(1))
	((TH1F *)(fQAAntiProtonsAcceptedList->At(21)))->Fill(0);
    }//point on SPD2
    if(fProtonAnalysisBase->IsUsedPointOnITSLayer3()) {
      if(!track->HasPointOnITSLayer(2)) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(22)))->Fill(0);
      }
      else if(track->HasPointOnITSLayer(2))
	((TH1F *)(fQAAntiProtonsAcceptedList->At(22)))->Fill(0);
    }//point on SDD1
    if(fProtonAnalysisBase->IsUsedPointOnITSLayer4()) {
      if(!track->HasPointOnITSLayer(3)) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(23)))->Fill(0);
      }
      else if(track->HasPointOnITSLayer(3))
	((TH1F *)(fQAAntiProtonsAcceptedList->At(23)))->Fill(0);
    }//point on SDD2
    if(fProtonAnalysisBase->IsUsedPointOnITSLayer5()) {
      if(!track->HasPointOnITSLayer(4)) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(24)))->Fill(0);
      }
      else if(track->HasPointOnITSLayer(4))
	((TH1F *)(fQAAntiProtonsAcceptedList->At(24)))->Fill(0);
    }//point on SSD1
    if(fProtonAnalysisBase->IsUsedPointOnITSLayer6()) {
      if(!track->HasPointOnITSLayer(5)) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(25)))->Fill(0);
      }
      else if(track->HasPointOnITSLayer(5))
	((TH1F *)(fQAAntiProtonsAcceptedList->At(25)))->Fill(0);
    }//point on SSD2
    if(fProtonAnalysisBase->IsUsedMinTPCdEdxPoints()) {
      if(track->GetTPCsignalN() < fProtonAnalysisBase->GetMinTPCdEdxPoints()) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(26)))->Fill(track->GetTPCsignalN());
      }
      if(track->GetTPCsignalN() >= fProtonAnalysisBase->GetMinTPCdEdxPoints())
	((TH1F *)(fQAAntiProtonsAcceptedList->At(26)))->Fill(track->GetTPCsignalN());
    }//number of TPC points for the dE/dx
  }//antiprotons
}
