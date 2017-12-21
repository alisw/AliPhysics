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

/* $Id: AliProtonAnalysis.cxx 42221 2010-07-11 13:13:27Z pchrist $ */

//-----------------------------------------------------------------
//                 AliProtonAnalysis class
//   This is the class to deal with the proton analysis
//   Origin: Panos Christakoglou | Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------
#include <Riostream.h>
#include <TFile.h>
#include <TSystem.h>
#include <TF1.h>
#include <TH3F.h>
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
#include <AliTPCPIDResponse.h>
#include <AliPIDResponse.h>
#include <AliESDpid.h>
#include <AliESDtrackCuts.h>
#include <AliMultiplicity.h>
class AliLog;
class AliESDVertex;

#include "AliProtonAnalysis.h"
#include "AliProtonAnalysisBase.h"

ClassImp(AliProtonAnalysis)

using std::cout;
using std::endl;

//____________________________________________________________________//
AliProtonAnalysis::AliProtonAnalysis() : 
  TObject(), fProtonAnalysisBase(0),
  fNBinsY(0), fMinY(0), fMaxY(0),
  fNBinsPt(0), fMinPt(0), fMaxPt(0),
  fProtonContainer(0), fAntiProtonContainer(0),
  fHistEvents(0),fHistMultiplicity(0), 
  fHistYPtProtonsCorrected(0), fHistYPtAntiProtonsCorrected(0), 
  fHistEventStats(0), fYRatioInPtBinsList(0),
  fHistEfficiencyYPtProtons(0), fHistEfficiencyYPtAntiProtons(0),
  fHistCorrectionForCrossSectionYPtProtons(0),
  fHistCorrectionForCrossSectionYPtAntiProtons(0),
  fHistCorrectionForCrossSectionFlag(kFALSE),
  fHistYPtCorrectionForCutsProtons(0), fHistYPtCorrectionForCutsAntiProtons(0),
  fCorrectForCutsFlag(kFALSE),
  fHistYPtCorrectionForTrackingProtons(0), 
  fHistYPtCorrectionForTrackingAntiProtons(0),
  fCorrectForTrackingFlag(kFALSE),
  fHistYPtCorrectionForFeedDownProtons(0), 
  fHistYPtCorrectionForFeedDownAntiProtons(0),
  fCorrectForFeedDownFlag(kFALSE),
  fHistYPtCorrectionForSecondaries(0), fCorrectForSecondariesFlag(kFALSE),
  fGlobalQAList(0), fQA2DList(0), fSysProtonLow(0),
  fQAProtonsAcceptedList(0), fQAProtonsRejectedList(0),
  fQAAntiProtonsAcceptedList(0), fQAAntiProtonsRejectedList(0) {
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
  fHistEvents(0),fHistMultiplicity(0),
  fHistYPtProtonsCorrected(0), fHistYPtAntiProtonsCorrected(0), 
  fHistEventStats(0), fYRatioInPtBinsList(0),
  fHistEfficiencyYPtProtons(0), fHistEfficiencyYPtAntiProtons(0),
  fHistCorrectionForCrossSectionYPtProtons(0),
  fHistCorrectionForCrossSectionYPtAntiProtons(0),
  fHistCorrectionForCrossSectionFlag(kFALSE),
  fHistYPtCorrectionForCutsProtons(0), fHistYPtCorrectionForCutsAntiProtons(0),
  fCorrectForCutsFlag(kFALSE),
  fHistYPtCorrectionForTrackingProtons(0), 
  fHistYPtCorrectionForTrackingAntiProtons(0),
  fCorrectForTrackingFlag(kFALSE),
  fHistYPtCorrectionForFeedDownProtons(0), 
  fHistYPtCorrectionForFeedDownAntiProtons(0),
  fCorrectForFeedDownFlag(kFALSE),
  fHistYPtCorrectionForSecondaries(0), fCorrectForSecondariesFlag(kFALSE),
  fGlobalQAList(0), fQA2DList(0), fSysProtonLow(0),
  fQAProtonsAcceptedList(0), fQAProtonsRejectedList(0),
  fQAAntiProtonsAcceptedList(0), fQAAntiProtonsRejectedList(0) {
  //Default constructor
  fHistEvents = new TH1I("fHistEvents","Analyzed events",2,0.5,2.5);
  fHistEvents->GetXaxis()->SetBinLabel(1,"Analyzed events");
  fHistEvents->GetXaxis()->SetBinLabel(2,"Events with (anti)protons");
  
  fHistMultiplicity = new TH1D("fHistMultiplicity","Multiplicity",200,0,200); //110
  fHistMultiplicity->GetXaxis()->SetTitle("Multiplicity");
  fHistMultiplicity->GetYaxis()->SetTitle("Counts");

  fHistYPtProtonsCorrected = new TH2D("fHistYPtProtonsCorrected","Protons",
				      fNBinsY,fMinY,fMaxY,
				      fNBinsPt,fMinPt,fMaxPt);
  fHistYPtProtonsCorrected->SetStats(kTRUE);
  fHistYPtProtonsCorrected->GetYaxis()->SetTitle("P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fHistYPtProtonsCorrected->GetXaxis()->SetTitle("#eta");
  else
    fHistYPtProtonsCorrected->GetXaxis()->SetTitle("y");
  fHistYPtProtonsCorrected->GetXaxis()->SetTitleColor(1);

  fHistYPtAntiProtonsCorrected = new TH2D("fHistYPtAntiProtonsCorrected",
					  "Antiprotons",
					  fNBinsY,fMinY,fMaxY,
					  fNBinsPt,fMinPt,fMaxPt);
  fHistYPtAntiProtonsCorrected->SetStats(kTRUE);
  fHistYPtAntiProtonsCorrected->GetYaxis()->SetTitle("P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fHistYPtAntiProtonsCorrected->GetXaxis()->SetTitle("#eta");
  else
    fHistYPtAntiProtonsCorrected->GetXaxis()->SetTitle("y");
  fHistYPtAntiProtonsCorrected->GetXaxis()->SetTitleColor(1);

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
  fProtonContainer->SetVarTitle(1,"P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fProtonContainer->SetVarTitle(0,"#eta");
  else
    fProtonContainer->SetVarTitle(0,"y");
  fAntiProtonContainer = new AliCFContainer("containerAntiProtons",
					    "container for antiprotons",
					    kNSteps,2,iBin);
  fAntiProtonContainer->SetBinLimits(0,binLimY); //rapidity or eta
  fAntiProtonContainer->SetBinLimits(1,binLimPt); //pT
  fAntiProtonContainer->SetVarTitle(1,"P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fAntiProtonContainer->SetVarTitle(0,"#eta");
  else
    fAntiProtonContainer->SetVarTitle(0,"y");

  //Initialize the QA
  if(fProtonAnalysisBase->IsQARun()) InitQA();
  InitSystematicsHistogram();
} 

//____________________________________________________________________//
AliProtonAnalysis::AliProtonAnalysis(Int_t nbinsY, Double_t *gY,
				     Int_t nbinsPt,Double_t *gPt) : 
  TObject(), fProtonAnalysisBase(0),
  fNBinsY(nbinsY), fMinY(gY[0]), fMaxY(gY[nbinsY]),
  fNBinsPt(nbinsPt), fMinPt(gPt[0]), fMaxPt(gPt[nbinsPt]),
  fProtonContainer(0), fAntiProtonContainer(0),
  fHistEvents(0),fHistMultiplicity(0),  
  fHistYPtProtonsCorrected(0), fHistYPtAntiProtonsCorrected(0), 
  fHistEventStats(0), fYRatioInPtBinsList(0),
  fHistEfficiencyYPtProtons(0), fHistEfficiencyYPtAntiProtons(0),
  fHistCorrectionForCrossSectionYPtProtons(0),
  fHistCorrectionForCrossSectionYPtAntiProtons(0),
  fHistCorrectionForCrossSectionFlag(kFALSE),
  fHistYPtCorrectionForCutsProtons(0), fHistYPtCorrectionForCutsAntiProtons(0),
  fCorrectForCutsFlag(kFALSE),
  fHistYPtCorrectionForTrackingProtons(0), 
  fHistYPtCorrectionForTrackingAntiProtons(0),
  fCorrectForTrackingFlag(kFALSE),
  fHistYPtCorrectionForFeedDownProtons(0), 
  fHistYPtCorrectionForFeedDownAntiProtons(0),
  fCorrectForFeedDownFlag(kFALSE),
  fHistYPtCorrectionForSecondaries(0), fCorrectForSecondariesFlag(kFALSE),
  fGlobalQAList(0), fQA2DList(0), fSysProtonLow(0),
  fQAProtonsAcceptedList(0), fQAProtonsRejectedList(0),
  fQAAntiProtonsAcceptedList(0), fQAAntiProtonsRejectedList(0) {
  //Default constructor
  fHistEvents = new TH1I("fHistEvents","Analyzed events",2,0.5,2.5);
  fHistEvents->GetXaxis()->SetBinLabel(1,"Analyzed events");
  fHistEvents->GetXaxis()->SetBinLabel(2,"Events with (anti)protons");
  
  fHistMultiplicity = new TH1D("fHistMultiplicity","Multiplicity",200,0,200);
  fHistMultiplicity->GetXaxis()->SetTitle("Multiplicity");
  fHistMultiplicity->GetYaxis()->SetTitle("Counts");

  fHistYPtProtonsCorrected = new TH2D("fHistYPtProtonsCorrected","Protons",
				      fNBinsY,fMinY,fMaxY,
				      fNBinsPt,fMinPt,fMaxPt);
  fHistYPtProtonsCorrected->SetStats(kTRUE);
  fHistYPtProtonsCorrected->GetYaxis()->SetTitle("P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fHistYPtProtonsCorrected->GetXaxis()->SetTitle("#eta");
  else
    fHistYPtProtonsCorrected->GetXaxis()->SetTitle("y");
  fHistYPtProtonsCorrected->GetXaxis()->SetTitleColor(1);

  fHistYPtAntiProtonsCorrected = new TH2D("fHistYPtAntiProtonsCorrected",
					  "Antiprotons",
					  fNBinsY,fMinY,fMaxY,
					  fNBinsPt,fMinPt,fMaxPt);
  fHistYPtAntiProtonsCorrected->SetStats(kTRUE);
  fHistYPtAntiProtonsCorrected->GetYaxis()->SetTitle("P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fHistYPtAntiProtonsCorrected->GetXaxis()->SetTitle("#eta");
  else
    fHistYPtAntiProtonsCorrected->GetXaxis()->SetTitle("y");
  fHistYPtAntiProtonsCorrected->GetXaxis()->SetTitleColor(1);

  //setting up the containers
  Int_t iBin[2];
  iBin[0] = nbinsY;
  iBin[1] = nbinsPt;
  fProtonContainer = new AliCFContainer("containerProtons",
					"container for protons",
					kNSteps,2,iBin);
  fProtonContainer->SetBinLimits(0,gY); //rapidity or eta
  fProtonContainer->SetBinLimits(1,gPt); //pT
  fProtonContainer->SetVarTitle(1,"P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fProtonContainer->SetVarTitle(0,"#eta");
  else
    fProtonContainer->SetVarTitle(0,"y");
  fAntiProtonContainer = new AliCFContainer("containerAntiProtons",
					    "container for antiprotons",
					    kNSteps,2,iBin);
  fAntiProtonContainer->SetBinLimits(0,gY); //rapidity or eta
  fAntiProtonContainer->SetBinLimits(1,gPt); //pT
  fAntiProtonContainer->SetVarTitle(1,"P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fAntiProtonContainer->SetVarTitle(0,"#eta");
  else
    fAntiProtonContainer->SetVarTitle(0,"y");

  //Initialize the QA
  if(fProtonAnalysisBase->IsQARun()) InitQA();
  InitSystematicsHistogram();
} 

//____________________________________________________________________//
AliProtonAnalysis::~AliProtonAnalysis() {
  //Default destructor
  if(fProtonAnalysisBase) delete fProtonAnalysisBase;
  
  if(fHistMultiplicity) delete fHistMultiplicity;
  if(fHistEvents) delete fHistEvents;
  if(fHistYPtProtonsCorrected) delete fHistYPtProtonsCorrected;
  if(fHistYPtAntiProtonsCorrected) delete fHistYPtAntiProtonsCorrected;
  if(fHistEventStats) delete fHistEventStats;
  if(fYRatioInPtBinsList) delete fYRatioInPtBinsList;

  if(fProtonContainer) delete fProtonContainer;
  if(fAntiProtonContainer) delete fAntiProtonContainer;

  if(fHistEfficiencyYPtProtons) delete fHistEfficiencyYPtProtons;
  if(fHistEfficiencyYPtAntiProtons) delete fHistEfficiencyYPtAntiProtons;
  if(fHistCorrectionForCrossSectionYPtProtons) delete fHistCorrectionForCrossSectionYPtProtons;
  if(fHistCorrectionForCrossSectionYPtAntiProtons) delete fHistCorrectionForCrossSectionYPtAntiProtons;
  if(fHistYPtCorrectionForCutsProtons) delete fHistYPtCorrectionForCutsProtons;
  if(fHistYPtCorrectionForCutsAntiProtons) delete fHistYPtCorrectionForCutsAntiProtons;
  if(fHistYPtCorrectionForTrackingProtons) delete fHistYPtCorrectionForTrackingProtons;
  if(fHistYPtCorrectionForTrackingAntiProtons) delete fHistYPtCorrectionForTrackingAntiProtons;
  if(fHistYPtCorrectionForFeedDownProtons) delete fHistYPtCorrectionForFeedDownProtons;
  if(fHistYPtCorrectionForFeedDownAntiProtons) delete fHistYPtCorrectionForFeedDownAntiProtons;
  if(fHistYPtCorrectionForSecondaries) delete fHistYPtCorrectionForSecondaries;

  //QA lists
  if(fGlobalQAList) delete fGlobalQAList;
  if(fSysProtonLow) delete fSysProtonLow;
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
  fNBinsY = nbinsY;
  fMinY = fLowY;
  fMaxY = fHighY;
  fNBinsPt = nbinsPt;
  fMinPt = fLowPt;
  fMaxPt = fHighPt;

  fHistEvents = new TH1I("fHistEvents","Analyzed events",2,0.5,2.5);
  fHistEvents->GetXaxis()->SetBinLabel(1,"Analyzed events");
  fHistEvents->GetXaxis()->SetBinLabel(2,"Events with (anti)protons");
  
  fHistMultiplicity = new TH1D("fHistMultiplicity","Multiplicity",200,0,200);
  fHistMultiplicity->GetXaxis()->SetTitle("Multiplicity");
  fHistMultiplicity->GetYaxis()->SetTitle("Counts");

  fHistYPtProtonsCorrected = new TH2D("fHistYPtProtonsCorrected","Protons",
				      fNBinsY,fMinY,fMaxY,
				      fNBinsPt,fMinPt,fMaxPt);
  fHistYPtProtonsCorrected->SetStats(kTRUE);
  fHistYPtProtonsCorrected->GetYaxis()->SetTitle("P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fHistYPtProtonsCorrected->GetXaxis()->SetTitle("#eta");
  else
    fHistYPtProtonsCorrected->GetXaxis()->SetTitle("y");
  fHistYPtProtonsCorrected->GetXaxis()->SetTitleColor(1);

  fHistYPtAntiProtonsCorrected = new TH2D("fHistYPtAntiProtonsCorrected",
					  "Antiprotons",
					  fNBinsY,fMinY,fMaxY,
					  fNBinsPt,fMinPt,fMaxPt);
  fHistYPtAntiProtonsCorrected->SetStats(kTRUE);
  fHistYPtAntiProtonsCorrected->GetYaxis()->SetTitle("P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fHistYPtAntiProtonsCorrected->GetXaxis()->SetTitle("#eta");
  else
    fHistYPtAntiProtonsCorrected->GetXaxis()->SetTitle("y");
  fHistYPtAntiProtonsCorrected->GetXaxis()->SetTitleColor(1);

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
  fProtonContainer->SetVarTitle(1,"P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fProtonContainer->SetVarTitle(0,"#eta");
  else
    fProtonContainer->SetVarTitle(0,"y");
  fAntiProtonContainer = new AliCFContainer("containerAntiProtons",
					    "container for antiprotons",
					    kNSteps,2,iBin);
  fAntiProtonContainer->SetBinLimits(0,binLimY); //rapidity
  fAntiProtonContainer->SetBinLimits(1,binLimPt); //pT
  fAntiProtonContainer->SetVarTitle(1,"P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fAntiProtonContainer->SetVarTitle(0,"#eta");
  else
    fAntiProtonContainer->SetVarTitle(0,"y");

  delete [] binLimY;
  delete [] binLimPt;

  //Initialize the QA
  if(fProtonAnalysisBase->IsQARun()) InitQA();
  InitSystematicsHistogram();
}

//____________________________________________________________________//
void AliProtonAnalysis::InitAnalysisHistograms(Int_t nbinsY, Double_t *gY, 
					       Int_t nbinsPt, Double_t *gPt) {
  //Initializes the histograms using asymmetric values - global tracking
  fNBinsY = nbinsY;
  fMinY = gY[0];
  fMaxY = gY[nbinsY];
  fNBinsPt = nbinsPt;
  fMinPt = gPt[0];
  fMaxPt = gPt[nbinsPt];

  fHistEvents = new TH1I("fHistEvents","Analyzed events",2,0.5,2.5);
  fHistEvents->GetXaxis()->SetBinLabel(1,"Analyzed events");
  fHistEvents->GetXaxis()->SetBinLabel(2,"Events with (anti)protons");
  
  fHistMultiplicity = new TH1D("fHistMultiplicity","Multiplicity",200,0,200);
  fHistMultiplicity->GetXaxis()->SetTitle("Multiplicity");
  fHistMultiplicity->GetYaxis()->SetTitle("Counts");

  fHistYPtProtonsCorrected = new TH2D("fHistYPtProtonsCorrected","Protons",
				      fNBinsY,gY,fNBinsPt,gPt);
  fHistYPtProtonsCorrected->SetStats(kTRUE);
  fHistYPtProtonsCorrected->GetYaxis()->SetTitle("P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fHistYPtProtonsCorrected->GetXaxis()->SetTitle("#eta");
  else
    fHistYPtProtonsCorrected->GetXaxis()->SetTitle("y");
  fHistYPtProtonsCorrected->GetXaxis()->SetTitleColor(1);

  fHistYPtAntiProtonsCorrected = new TH2D("fHistYPtAntiProtonsCorrected",
					  "Antiprotons",
					  fNBinsY,gY,fNBinsPt,gPt);
  fHistYPtAntiProtonsCorrected->SetStats(kTRUE);
  fHistYPtAntiProtonsCorrected->GetYaxis()->SetTitle("P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fHistYPtAntiProtonsCorrected->GetXaxis()->SetTitle("#eta");
  else
    fHistYPtAntiProtonsCorrected->GetXaxis()->SetTitle("y");
  fHistYPtAntiProtonsCorrected->GetXaxis()->SetTitleColor(1);

  //setting up the containers
  Int_t iBin[2];
  iBin[0] = nbinsY;
  iBin[1] = nbinsPt;

  fProtonContainer = new AliCFContainer("containerProtons",
					"container for protons",
					kNSteps,2,iBin);
  fProtonContainer->SetBinLimits(0,gY); //rapidity
  fProtonContainer->SetBinLimits(1,gPt); //pT
  fProtonContainer->SetVarTitle(1,"P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fProtonContainer->SetVarTitle(0,"#eta");
  else
    fProtonContainer->SetVarTitle(0,"y");
  fAntiProtonContainer = new AliCFContainer("containerAntiProtons",
					    "container for antiprotons",
					    kNSteps,2,iBin);
  fAntiProtonContainer->SetBinLimits(0,gY); //rapidity
  fAntiProtonContainer->SetBinLimits(1,gPt); //pT
  fAntiProtonContainer->SetVarTitle(1,"P_{T} [GeV/c]");
  if(fProtonAnalysisBase->GetEtaMode())
    fAntiProtonContainer->SetVarTitle(0,"#eta");
  else
    fAntiProtonContainer->SetVarTitle(0,"y");

  //Initialize the QA
  if(fProtonAnalysisBase->IsQARun()) InitQA();
  InitSystematicsHistogram();
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
    fProtonContainer = (AliCFContainer *)list->At(0);
    fAntiProtonContainer = (AliCFContainer *)list->At(1);
	fHistEvents = (TH1I *)list->At(2);
    fHistEventStats = (TH1F *)list->At(3);
    fHistMultiplicity = (TH1D*)list->At(4);
  }
  else if(!list) {
    cout<<"Retrieving objects from the file... "<<endl;
    fHistEvents = (TH1I *)file->Get("fHistEvents");
    fProtonContainer = (AliCFContainer *)file->Get("containerProtons");
    fAntiProtonContainer = (AliCFContainer *)file->Get("containerAntiProtons");
    fHistEventStats = (TH1F *)file->Get("fHistEventStats");
  }
  if((!fHistEvents)||(!fProtonContainer)||(!fAntiProtonContainer)||(!fHistEventStats)) {
    cout<<"Input containers were not found!!!"<<endl;
    status = kFALSE;
}
  return status;
}

//____________________________________________________________________//
TList *AliProtonAnalysis::GetYRatioHistogramsInPtBins() {
  //Returns a TList obkect with the eta (or y) dependent ratios for each 
  //pT bin taken from the 2D histograms (not from the containers)
  fYRatioInPtBinsList = new TList();
  
  TH2D *fHistYPtProtons = GetProtonYPtHistogram();
  TH2D *fHistYPtAntiProtons = GetAntiProtonYPtHistogram();
  
  //Protons
  TH1D *gHistYProtons[100];
  TString title;
  for(Int_t iBin = 1; iBin <= fHistYPtProtons->GetNbinsY(); iBin++) {
    title = "gHistYProtons_PtBin"; title += iBin;
    gHistYProtons[iBin] = (TH1D *)fHistYPtProtons->ProjectionX(title.Data(),
							       iBin,
							       iBin,"b");
    gHistYProtons[iBin]->Sumw2();
  }

  //Antiprotons
  TH1D *gHistYAntiProtons[100];
  for(Int_t iBin = 1; iBin <= fHistYPtAntiProtons->GetNbinsY(); iBin++) {
    title = "gHistYAntiProtons_PtBin"; title += iBin;
    gHistYAntiProtons[iBin] = (TH1D *)fHistYPtAntiProtons->ProjectionX(title.Data(),
								       iBin,
								       iBin,"b");
    gHistYAntiProtons[iBin]->Sumw2();
  }

  Double_t pTmin = fHistYPtProtons->GetYaxis()->GetXmin();
  Double_t pTStep = (fHistYPtProtons->GetYaxis()->GetXmax() - fHistYPtProtons->GetYaxis()->GetXmin())/fHistYPtProtons->GetNbinsY();
  Double_t pTmax = pTmin + pTStep;
  //Ratio
  TH1D *gHistYRatio[100];
  for(Int_t iBin = 1; iBin <= fHistYPtProtons->GetNbinsY(); iBin++) {
    title = "gHistYRatio_PtBin"; title += iBin;
    gHistYRatio[iBin] = new TH1D(title.Data(),"",
				 fHistYPtProtons->GetNbinsX(),
				 fHistYPtProtons->GetXaxis()->GetXmin(),
				 fHistYPtProtons->GetXaxis()->GetXmax());
    title = "Pt: "; title += pTmin; title += " - "; title += pTmax;
    gHistYRatio[iBin]->SetTitle(title.Data());
    gHistYRatio[iBin]->GetYaxis()->SetTitle("#bar{p}/p");
    gHistYRatio[iBin]->GetXaxis()->SetTitle(fHistYPtProtons->GetXaxis()->GetTitle());
    gHistYRatio[iBin]->Divide(gHistYAntiProtons[iBin],
			      gHistYProtons[iBin],1.0,1.0,"b");
    fYRatioInPtBinsList->Add(gHistYRatio[iBin]);
    pTmin += pTStep;
    pTmax += pTStep;
  }
  
  return fYRatioInPtBinsList;
}
//____________________________________________________________________//
TH2D *AliProtonAnalysis::GetProtonYPtHistogram() {
  //Get the y-Pt histogram for protons
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();

  TH2D *fYPtProtons = fProtonContainer->ShowProjection(0,1,kStepInPhaseSpace); //variable-step
   
  fYPtProtons->SetStats(kFALSE);
  fYPtProtons->GetZaxis()->SetTitle("(1/N_{events})(dN/dydPt)");
  fYPtProtons->SetTitle("dN/dydPt protons");
  
  
  return fYPtProtons;
}

//____________________________________________________________________//
TH2D *AliProtonAnalysis::GetAntiProtonYPtHistogram() {
  //Get the y-Pt histogram for AntiProtons
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();

  TH2D *fYPtAntiProtons = fAntiProtonContainer->ShowProjection(0,1,kStepInPhaseSpace); //variable-step
   
  fYPtAntiProtons->SetStats(kFALSE);
  fYPtAntiProtons->GetZaxis()->SetTitle("(1/N_{events})(dN/dydPt)");
  fYPtAntiProtons->SetTitle("dN/dydPt AntiProtons");
  
  
  return fYPtAntiProtons;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetProtonYHistogram() {
  //Get the y histogram for protons
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();

  TH1D *fYProtons = fProtonContainer->ShowProjection(0,kStepInPhaseSpace); //variable-step
   
  fYProtons->SetStats(kFALSE);
  fYProtons->GetYaxis()->SetTitle("(1/N_{events})(dN/dy)");
  fYProtons->SetTitle("dN/dy protons");
  fYProtons->SetMarkerStyle(kFullCircle);
  fYProtons->SetMarkerColor(4);
  
  return fYProtons;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetAntiProtonYHistogram() {
  //Get the y histogram for antiprotons
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();
  
  TH1D *fYAntiProtons = fAntiProtonContainer->ShowProjection(0,kStepInPhaseSpace);//variable-step 
 
  fYAntiProtons->SetStats(kFALSE);
  fYAntiProtons->GetYaxis()->SetTitle("(1/N_{events})(dN/dy)");
  fYAntiProtons->SetTitle("dN/dy antiprotons");
  fYAntiProtons->SetMarkerStyle(kFullCircle);
  fYAntiProtons->SetMarkerColor(4);
  
  return fYAntiProtons;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetProtonPtHistogram() {
  //Get the Pt histogram for protons
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();
  
  TH1D *fPtProtons = fProtonContainer->ShowProjection(1,kStepInPhaseSpace); //variable-step

  fPtProtons->SetStats(kFALSE);
  fPtProtons->GetYaxis()->SetTitle("(1/N_{events})(dN/dP_{T})");
  fPtProtons->SetTitle("dN/dPt protons");
  fPtProtons->SetMarkerStyle(kFullCircle);
  fPtProtons->SetMarkerColor(4);
 
  return fPtProtons;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetAntiProtonPtHistogram() {
  //Get the Pt histogram for antiprotons
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();
  
  TH1D *fPtAntiProtons = fAntiProtonContainer->ShowProjection(1,kStepInPhaseSpace); //variable-step

  fPtAntiProtons->SetStats(kFALSE);
  fPtAntiProtons->GetYaxis()->SetTitle("(1/N_{events})(dN/dP_{T})");
  fPtAntiProtons->SetTitle("dN/dPt antiprotons");
  fPtAntiProtons->SetMarkerStyle(kFullCircle);
  fPtAntiProtons->SetMarkerColor(4);
 
  return fPtAntiProtons;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetProtonCorrectedYHistogram() {
  //Get the corrected y histogram for protons
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();
  
  TH2D *fHistYPtProtons = GetProtonYPtHistogram();
  
  TH1D *fYProtons = (TH1D *)fHistYPtProtonsCorrected->ProjectionX("fYProtons",1,fHistYPtProtons->GetYaxis()->GetNbins(),"e");
   
  fYProtons->SetStats(kFALSE);
  fYProtons->GetYaxis()->SetTitle("(1/N_{events})(dN/dy)");
  fYProtons->GetXaxis()->SetTitle("y");
  fYProtons->SetTitle("dN/dy protons");
  fYProtons->SetMarkerStyle(kFullCircle);
  fYProtons->SetMarkerColor(4);
    
  return fYProtons;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetAntiProtonCorrectedYHistogram() {
  //Get the corrected y histogram for antiprotons
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();
  
  TH2D *fHistYPtAntiProtons = GetAntiProtonYPtHistogram();
  
  TH1D *fYAntiProtons = (TH1D *)fHistYPtAntiProtonsCorrected->ProjectionX("fYAntiProtons",1,fHistYPtAntiProtons->GetYaxis()->GetNbins(),"e");
   
  fYAntiProtons->SetStats(kFALSE);
  fYAntiProtons->GetYaxis()->SetTitle("(1/N_{events})(dN/dy)");
  fYAntiProtons->GetXaxis()->SetTitle("y");
  fYAntiProtons->SetTitle("dN/dy protons");
  fYAntiProtons->SetMarkerStyle(kFullCircle);
  fYAntiProtons->SetMarkerColor(4);
   
  return fYAntiProtons;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetProtonCorrectedPtHistogram(Int_t Low = 1,Int_t High = 100) {
  //Get the corrected Pt histogram for protons
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();

  TH2D *fHistYPtProtons = GetProtonYPtHistogram();
  
  TH1D *fPtProtons = (TH1D *)fHistYPtProtonsCorrected->ProjectionY("fPtProtons",Low,High,"e");
   
  fPtProtons->SetStats(kFALSE);
  fPtProtons->GetYaxis()->SetTitle("(1/N_{events})(dN/dP_{T})");
  fPtProtons->GetXaxis()->SetTitle("P_{T} [GeV/c]");
  fPtProtons->SetTitle("dN/dPt protons");
  fPtProtons->SetMarkerStyle(kFullCircle);
  fPtProtons->SetMarkerColor(4);
   
  return fPtProtons;
}


//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetAntiProtonCorrectedPtHistogram(Int_t Low = 1,Int_t High = 100) {
  //Get the corrected Pt histogram for antiprotons
  Int_t nAnalyzedEvents = GetNumberOfAnalyzedEvents();

  TH2D *fHistYPtAntiProtons = GetAntiProtonYPtHistogram();
  
  TH1D *fPtAntiProtons = (TH1D *)fHistYPtAntiProtonsCorrected->ProjectionY("fPtAntiProtons",Low,High,"e");//fHistYPtAntiProtonsCorrected

   
  fPtAntiProtons->SetStats(kFALSE);
  fPtAntiProtons->GetYaxis()->SetTitle("(1/N_{events})(dN/dP_{T})");
  fPtAntiProtons->GetXaxis()->SetTitle("P_{T} [GeV/c]");
  fPtAntiProtons->SetTitle("dN/dPt antiprotons");
  fPtAntiProtons->SetMarkerStyle(kFullCircle);
  fPtAntiProtons->SetMarkerColor(4);
   
  return fPtAntiProtons;
}

//____________________________________________________________________//
TH1D *AliProtonAnalysis::GetYRatioHistogram() {
  //Returns the rapidity dependence of the ratio (uncorrected)
  TH1D *fYProtons = GetProtonYHistogram();
  TH1D *fYAntiProtons = GetAntiProtonYHistogram();
  
  TH1D *hRatioY = new TH1D("hRatioY","",fYProtons->GetNbinsX(),fYProtons->GetXaxis()->GetXmin(),fYProtons->GetXaxis()->GetXmax());
  hRatioY->Divide(fYAntiProtons,fYProtons,1.0,1.0,"e");
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
TH1D *AliProtonAnalysis::GetYRatioCorrectedHistogram() {
  //Returns the rapidity dependence of the ratio (corrected)
  //fHistYPtProtons->Multiply(gCorrectionProtons);
  TH1D *fYProtons = GetProtonCorrectedYHistogram();
  //fHistYPtAntiProtons->Multiply(gCorrectionAntiProtons);
  TH1D *fYAntiProtons = GetAntiProtonCorrectedYHistogram();
  
  TH1D *hRatioY = new TH1D("hRatioY","",fYProtons->GetNbinsX(),fYProtons->GetXaxis()->GetXmin(),fYProtons->GetXaxis()->GetXmax());
  hRatioY->Divide(fYAntiProtons,fYProtons,1.0,1.0,"e");
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
  hRatioPt->Divide(fPtAntiProtons,fPtProtons,1.0,1.0,"e");
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
TH1D *AliProtonAnalysis::GetPtRatioCorrectedHistogram(Int_t Low = 1,Int_t High = 100,Int_t Rebin = 1) {
  //Returns the Pt dependence of the ratio (corrected)
  //fHistYPtProtons->Multiply(gCorrectionProtons);
  TH1D *fPtProtons = GetProtonCorrectedPtHistogram(Low,High);
  //fHistYPtAntiProtons->Multiply(gCorrectionAntiProtons);
  TH1D *fPtAntiProtons = GetAntiProtonCorrectedPtHistogram(Low,High);
  fPtProtons->Rebin(Rebin);
  fPtAntiProtons->Rebin(Rebin);
  
  TH1D *hRatioPt = new TH1D("hRatioPt","",fPtProtons->GetNbinsX(),fPtProtons->GetXaxis()->GetXmin(),fPtProtons->GetXaxis()->GetXmax());
  hRatioPt->Divide(fPtAntiProtons,fPtProtons,1.0,1.0,"e");
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
  hAsymmetryY->Divide(hdiff,hsum,2.0,1.,"e");
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
TH1D *AliProtonAnalysis::GetYAsymmetryCorrectedHistogram() {
  //Returns the rapidity dependence of the asymmetry (corrected)
  TH1D *fYProtons = GetProtonCorrectedYHistogram();
  TH1D *fYAntiProtons = GetAntiProtonCorrectedYHistogram();
  
  TH1D *hsum = new TH1D("hsumY","",fYProtons->GetNbinsX(),fYProtons->GetXaxis()->GetXmin(),fYProtons->GetXaxis()->GetXmax());
  hsum->Add(fYProtons,fYAntiProtons,1.0,1.0);

  TH1D *hdiff = new TH1D("hdiffY","",fYProtons->GetNbinsX(),fYProtons->GetXaxis()->GetXmin(),fYProtons->GetXaxis()->GetXmax());
  hdiff->Add(fYProtons,fYAntiProtons,1.0,-1.0);

  TH1D *hAsymmetryY = new TH1D("hAsymmetryY","",fYProtons->GetNbinsX(),fYProtons->GetXaxis()->GetXmin(),fYProtons->GetXaxis()->GetXmax());
  hAsymmetryY->Divide(hdiff,hsum,2.0,1.,"e");
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
  hAsymmetryPt->Divide(hdiff,hsum,2.0,1.,"e");
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
TH1D *AliProtonAnalysis::GetPtAsymmetryCorrectedHistogram(Int_t Low = 1,Int_t High = 100) {
  //Returns the pT dependence of the asymmetry (corrected)
  TH1D *fPtProtons = GetProtonCorrectedPtHistogram(Low,High);
  TH1D *fPtAntiProtons = GetAntiProtonCorrectedPtHistogram(Low,High);
  
  TH1D *hsum = new TH1D("hsumPt","",fPtProtons->GetNbinsX(),fPtProtons->GetXaxis()->GetXmin(),fPtProtons->GetXaxis()->GetXmax());
  hsum->Add(fPtProtons,fPtAntiProtons,1.0,1.0);

  TH1D *hdiff = new TH1D("hdiffPt","",fPtProtons->GetNbinsX(),fPtProtons->GetXaxis()->GetXmin(),fPtProtons->GetXaxis()->GetXmax());
  hdiff->Add(fPtProtons,fPtAntiProtons,1.0,-1.0);

  TH1D *hAsymmetryPt = new TH1D("hAsymmetryPt","",fPtProtons->GetNbinsX(),fPtProtons->GetXaxis()->GetXmin(),fPtProtons->GetXaxis()->GetXmax());
  hAsymmetryPt->Divide(hdiff,hsum,2.0,1.,"e");
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

 /* const AliMultiplicity *fMult = esd->GetMultiplicity();
  if (!fMult){
	AliError("Can't get multiplicity object");
	return;
	}
  Int_t nTracklets = fMult->GetNumberOfTracklets();*/

AliCentrality *esdCentrality = esd->GetCentrality();
Float_t nTracklets = esdCentrality->GetCentralityPercentile("V0M");

//if(fProtonAnalysisBase->IsUsedITSSAMultiplicitySelection()){

/*  AliESDtrackCuts *fTrackCuts = new AliESDtrackCuts();
  if(!fTrackCuts){
	AliError("Can't get track cut object");
	}
  Int_t nTracklets = fTrackCuts->GetReferenceMultiplicity(esd,AliESDtrackCuts::kTrackletsITSTPC,0.5);*/
//}


fHistMultiplicity->Fill(nTracklets);
  
  fHistEvents->Fill(1); //number of analyzed events
  Double_t containerInput[2] ;
  Double_t gPt = 0.0;


  nTracks = esd->GetNumberOfTracks();
  for(Int_t iTracks = 0; iTracks < nTracks; iTracks++) {
    AliESDtrack* track = esd->GetTrack(iTracks);

   
if((fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kTPC)||
		(fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kHybrid)|| 
		(fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kFullHybrid)) {

		
      AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
      if(!tpcTrack) continue;


      gPt = tpcTrack->Pt();

//if(TMath::Abs(tpcTrack->Eta()) > 0.5) continue;//acceptance

      if(fProtonAnalysisBase->GetEtaMode())
	containerInput[0] = tpcTrack->Eta();
      if(fProtonAnalysisBase->IsUsedMultiplicitySelection())
        containerInput[0] = nTracklets;
      else
	containerInput[0] = fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
							  tpcTrack->Py(),
							  tpcTrack->Pz());
    containerInput[1] = gPt;

       //Step: kStepSurvived
      if(fProtonAnalysisBase->IsAccepted(track)) {
	
	if(tpcTrack->Charge() > 0) {
	  fProtonContainer->Fill(containerInput,kStepSurvived);   
	}//protons
	else if(tpcTrack->Charge() < 0) {
	  fAntiProtonContainer->Fill(containerInput,kStepSurvived);   
	}//antiprotons

	//Step: kStepIdentified
	if(fProtonAnalysisBase->IsProton(track)) {
	  if(tpcTrack->Charge() > 0)
	    fProtonContainer->Fill(containerInput,kStepIdentified);   
	  else if(tpcTrack->Charge() < 0) 
	    fAntiProtonContainer->Fill(containerInput,kStepIdentified);   
	
	  //Step: kStepIsPrimary 	  
	  if(fProtonAnalysisBase->IsPrimary(esd,vertex,track)) {
	    if(tpcTrack->Charge() > 0) {
	      nIdentifiedProtons += 1;
	      fProtonContainer->Fill(containerInput,kStepIsPrimary);
	    }//protons
	    else if(tpcTrack->Charge() < 0) {
	      nIdentifiedAntiProtons += 1;
	      fAntiProtonContainer->Fill(containerInput,kStepIsPrimary);
	    }//antiprotons
	    
	    //Step: kStepInPhaseSpace
	    if(fProtonAnalysisBase->IsInPhaseSpace(track)) {
	      if(tpcTrack->Charge() > 0) {
		nSurvivedProtons += 1;
		fProtonContainer->Fill(containerInput,kStepInPhaseSpace);   
	      }//protons
	      else if(tpcTrack->Charge() < 0) {
		nSurvivedAntiProtons += 1;
		fAntiProtonContainer->Fill(containerInput,kStepInPhaseSpace);
	      }//antiprotons
	    }//Step: kStepInPhaseSpace
	  }//Step: kStepIsPrimary
	}//Step: kStepIdentified
      }//Step: kStepSurvived
    }//TPC only tracks
    
    else if(fProtonAnalysisBase->GetAnalysisMode() == AliProtonAnalysisBase::kGlobal) {
      gPt = track->Pt();

      if(fProtonAnalysisBase->GetEtaMode())
	containerInput[0] = track->Eta();
      if(fProtonAnalysisBase->IsUsedMultiplicitySelection())
        containerInput[0] = nTracklets;
      else
	containerInput[0] = fProtonAnalysisBase->Rapidity(track->Px(),
							  track->Py(),
							  track->Pz());
    containerInput[1] = gPt;
      
      //Step: kStepSurvived
      if(fProtonAnalysisBase->IsAccepted(track)) {
	if(track->Charge() > 0) {
	  fProtonContainer->Fill(containerInput,kStepSurvived);   
	}//protons
	else if(track->Charge() < 0) {
	  fAntiProtonContainer->Fill(containerInput,kStepSurvived);   
	}//antiprotons
	
	//Step: kStepIsPrimary
	if(fProtonAnalysisBase->IsProton(track)) {
	  if(track->Charge() > 0)
	    fProtonContainer->Fill(containerInput,kStepIdentified);   
	  else if(track->Charge() < 0) 
	    fAntiProtonContainer->Fill(containerInput,kStepIdentified);   
	  
	  //Step: kStepIdentified
	  if(fProtonAnalysisBase->IsPrimary(esd,vertex,track)) {
	    if(track->Charge() > 0) {
	      nIdentifiedProtons += 1;
	      fProtonContainer->Fill(containerInput,kStepIsPrimary);
	    }//protons
	    else if(track->Charge() < 0) {
	      nIdentifiedAntiProtons += 1;
	      fAntiProtonContainer->Fill(containerInput,kStepIsPrimary);
	    }//protons
	    
	    //Step: kStepInPhaseSpace
	    if(fProtonAnalysisBase->IsInPhaseSpace(track)) {
	      if(track->Charge() > 0) {
		nSurvivedProtons += 1;
		fProtonContainer->Fill(containerInput,kStepInPhaseSpace);   
	      }//protons
	      else if(track->Charge() < 0) {
		nSurvivedAntiProtons += 1;
		fAntiProtonContainer->Fill(containerInput,kStepInPhaseSpace);
	      }//antiprotons
	    }//Step: kStepInPhaseSpace
	  }//Step: kStepIsPrimary
	}//Step: kStepIdentified
      }//Step: kStepSurvived
    }//Global tracking
  }//track loop 
  
  if((nIdentifiedProtons > 0)||(nIdentifiedAntiProtons > 0))
    fHistEvents->Fill(2); //number of analyzed events with at least one (anti)proton
  
  if(fProtonAnalysisBase->GetDebugMode())
    Printf("Initial number of tracks: %d | Identified (anti)protons: %d - %d | Survived (anti)protons: %d - %d",nTracks,nIdentifiedProtons,nIdentifiedAntiProtons,nSurvivedProtons,nSurvivedAntiProtons);
}


//____________________________________________________________________//
void AliProtonAnalysis::AnalyzeQA(AliESDEvent* esd,
				const AliESDVertex *vertex) {
  //Main analysis part - ESD
  Int_t nTracks = 0;
  Int_t nIdentifiedProtons = 0, nIdentifiedAntiProtons = 0;
  Int_t nSurvivedProtons = 0, nSurvivedAntiProtons = 0;
  
 /* const AliMultiplicity *fMult = esd->GetMultiplicity();
  if (!fMult){
	AliError("Can't get multiplicity object");
	return;
	}
  Int_t nTracklets = fMult->GetNumberOfTracklets();*/
AliCentrality *esdCentrality = esd->GetCentrality();
Float_t nTracklets = esdCentrality->GetCentralityPercentile("V0M");
//Printf("Centrality in base: %f",nTracklets);

//if(fProtonAnalysisBase->IsUsedITSSAMultiplicitySelection()){ 
//nTracklets = 0;
/*AliESDtrackCuts *fTrackCuts = new AliESDtrackCuts();
  if(!fTrackCuts){
	AliError("Can't get track cut object");
	}
  Int_t nTracklets = fTrackCuts->GetReferenceMultiplicity(esd,AliESDtrackCuts::kTrackletsITSTPC,0.5);*/
//Printf("Multiplicity in base: %i",nTracklets);
//}
  fHistMultiplicity->Fill(nTracklets);
  
  //=========================================//
  AliTPCPIDResponse tpcResponse = dynamic_cast<AliTPCPIDResponse&>(dynamic_cast<AliPIDResponse*>(fProtonAnalysisBase->GetPIDResponse())->GetTPCResponse());
  //=========================================//  

  fHistEvents->Fill(1); //number of analyzed events
  Double_t containerInput[2] ;
  Double_t gPt = 0.0, gP = 0.0;
  Float_t dcaXY = 0.0, dcaZ = 0.0;

  nTracks = esd->GetNumberOfTracks();
  for(Int_t iTracks = 0; iTracks < nTracks; iTracks++) {
    AliESDtrack* track = esd->GetTrack(iTracks);
    AliESDtrack trackTPC;

    Int_t nClustersTPC = track->GetTPCclusters(0x0);
    Int_t npointsTPCdEdx = track->GetTPCsignalN();
    Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};  //The impact parameters and their covariance.
    Double_t dca3D = 0.0;

    if((fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kTPC)||
    		(fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kFullHybrid)|| 
    		(fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kHybrid)) {
      AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
      if(!tpcTrack) continue;
  //    if(TMath::Abs(tpcTrack->Eta()) > 0.5) continue;//acceptance

      gPt = tpcTrack->Pt();
      gP = track->GetInnerParam()->P();
      
      if((fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kTPC)||
      		(fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kHybrid)) {
      tpcTrack->PropagateToDCA(vertex,
			       esd->GetMagneticField(),
			       100.,dca,cov);
      
	}
      else if(fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kFullHybrid){	  
      AliExternalTrackParam cParam;
      track->RelateToVertex(vertex,
			    esd->GetMagneticField(),
			    100.,&cParam);
      track->GetImpactParameters(dcaXY,dcaZ);
      dca[0] = dcaXY; dca[1] = dcaZ;
        }
      
      dca3D = TMath::Sqrt(TMath::Power(dca[0],2) +
			  TMath::Power(dca[1],2));

      if(fProtonAnalysisBase->GetEtaMode())
	containerInput[0] = tpcTrack->Eta();
      if(fProtonAnalysisBase->IsUsedMultiplicitySelection())
        containerInput[0] = nTracklets;
      else
	containerInput[0] = fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
							  tpcTrack->Py(),
							  tpcTrack->Pz());
      containerInput[1] = gPt;
      
      if(fProtonAnalysisBase->IsProton(track)) FillQA(esd,vertex,track);//QA for (anti)protons

       //Step: kStepSurvived
      if(fProtonAnalysisBase->IsAccepted(track)) {
	  ((TH2F *)(fQA2DList->At(1)))->Fill(gP,track->GetTPCsignal());
	  if((track->GetTPCsignal() > 0.0) && (tpcResponse.GetExpectedSignal(gP,AliPID::kProton) > 0.0))((TH2F *)(fQA2DList->At(3)))->Fill(gP,TMath::Log(track->GetTPCsignal()/tpcResponse.GetExpectedSignal(gP,AliPID::kProton)));
	  ((TH3F *)(fQA2DList->At(5)))->Fill(tpcTrack->Eta(),
					     tpcTrack->Phi()*180./TMath::Pi(),
					     npointsTPCdEdx);
	  ((TH3F *)(fQA2DList->At(7)))->Fill(tpcTrack->Eta(),
					     tpcTrack->Phi()*180./TMath::Pi(),
					     nClustersTPC);
	  ((TH3F *)(fQA2DList->At(9)))->Fill(gPt,
					     tpcTrack->Phi()*180./TMath::Pi(),
					     npointsTPCdEdx);
	  ((TH3F *)(fQA2DList->At(11)))->Fill(gPt,
					      tpcTrack->Phi()*180./TMath::Pi(),
					      nClustersTPC);
		
	if(tpcTrack->Charge() > 0) 
	  fProtonContainer->Fill(containerInput,kStepSurvived);    

	else if(tpcTrack->Charge() < 0) 
	  fAntiProtonContainer->Fill(containerInput,kStepSurvived);    

	//Step: kStepIdentified
	if(fProtonAnalysisBase->IsProton(track)) {
	    ((TH2F *)(fQA2DList->At(0)))->Fill(gP,track->GetTPCsignal());
	    if((track->GetTPCsignal() > 0.0) && (tpcResponse.GetExpectedSignal(gP,AliPID::kProton) > 0.0))((TH2F *)(fQA2DList->At(2)))->Fill(gP,TMath::Log(track->GetTPCsignal()/tpcResponse.GetExpectedSignal(gP,AliPID::kProton)));
	    ((TH3F *)(fQA2DList->At(4)))->Fill(tpcTrack->Eta(),
					       tpcTrack->Phi()*180./TMath::Pi(),
					       npointsTPCdEdx);
	    ((TH3F *)(fQA2DList->At(6)))->Fill(tpcTrack->Eta(),
					       tpcTrack->Phi()*180./TMath::Pi(),
					       nClustersTPC);
	    ((TH3F *)(fQA2DList->At(8)))->Fill(gPt,
					       tpcTrack->Phi()*180./TMath::Pi(),
					       npointsTPCdEdx);
	    ((TH3F *)(fQA2DList->At(10)))->Fill(gPt,
						tpcTrack->Phi()*180./TMath::Pi(),
						nClustersTPC);	
	
	  if(tpcTrack->Charge() > 0){
	    fProtonContainer->Fill(containerInput,kStepIdentified);
	    
	      ((TH2F *)(fQA2DList->At(20)))->Fill(gP,track->GetTPCsignal());
	      ((TH2F *)(fQA2DList->At(12)))->Fill(tpcTrack->Eta(),
						tpcTrack->Phi()*180./TMath::Pi());
	    if(fProtonAnalysisBase->GetEtaMode()) {
	      ((TH3F *)(fQA2DList->At(14)))->Fill(tpcTrack->Eta(),
						  tpcTrack->Pt(),
						  dca[0]);
	      ((TH3F *)(fQA2DList->At(15)))->Fill(tpcTrack->Eta(),
						  tpcTrack->Pt(),
						  dca[1]);
	      ((TH3F *)(fQA2DList->At(18)))->Fill(tpcTrack->Eta(),
						  tpcTrack->Pt(),
						  TMath::Abs(dca3D));
	    }
	    if(fProtonAnalysisBase->IsUsedMultiplicitySelection()){
	      ((TH3F *)(fQA2DList->At(14)))->Fill(nTracklets,
						  tpcTrack->Pt(),
						  dca[0]);
	      ((TH3F *)(fQA2DList->At(15)))->Fill(nTracklets,
						  tpcTrack->Pt(),
						  dca[1]);
	      ((TH3F *)(fQA2DList->At(18)))->Fill(nTracklets,
						  tpcTrack->Pt(),
						  TMath::Abs(dca3D));

	    }
	    else {
	      ((TH3F *)(fQA2DList->At(14)))->Fill(fProtonAnalysisBase->Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz()),
						  tpcTrack->Pt(),
						  dca[0]);
	      ((TH3F *)(fQA2DList->At(15)))->Fill(fProtonAnalysisBase->Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz()),
						  tpcTrack->Pt(),
						  dca[1]);
	      ((TH3F *)(fQA2DList->At(18)))->Fill(fProtonAnalysisBase->Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz()),
						  tpcTrack->Pt(),
						  TMath::Abs(dca3D));
	    }
	  }//protons
	       
	  else if(tpcTrack->Charge() < 0){ 
	    fAntiProtonContainer->Fill(containerInput,kStepIdentified);
	      ((TH2F *)(fQA2DList->At(21)))->Fill(gP,track->GetTPCsignal());
	      ((TH2F *)(fQA2DList->At(13)))->Fill(tpcTrack->Eta(),
						tpcTrack->Phi()*180./TMath::Pi());
	    if(fProtonAnalysisBase->GetEtaMode()) {
	      ((TH3F *)(fQA2DList->At(16)))->Fill(tpcTrack->Eta(),
						  tpcTrack->Pt(),
						  dca[0]);
	      ((TH3F *)(fQA2DList->At(17)))->Fill(tpcTrack->Eta(),
						  tpcTrack->Pt(),
						  dca[1]);
	      ((TH3F *)(fQA2DList->At(19)))->Fill(tpcTrack->Eta(),
						  tpcTrack->Pt(),
						  TMath::Abs(dca3D));
	    }
	    if(fProtonAnalysisBase->IsUsedMultiplicitySelection()){
	      ((TH3F *)(fQA2DList->At(16)))->Fill(nTracklets,
						  tpcTrack->Pt(),
						  dca[0]);
	      ((TH3F *)(fQA2DList->At(17)))->Fill(nTracklets,
						  tpcTrack->Pt(),
						  dca[1]);
	      ((TH3F *)(fQA2DList->At(19)))->Fill(nTracklets,
						  tpcTrack->Pt(),
						  TMath::Abs(dca3D));

	    }
	    else {
	      ((TH3F *)(fQA2DList->At(16)))->Fill(fProtonAnalysisBase->Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz()),
						  tpcTrack->Pt(),
						  dca[0]);
	      ((TH3F *)(fQA2DList->At(17)))->Fill(fProtonAnalysisBase->Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz()),
						  tpcTrack->Pt(),
						  dca[1]);
	      ((TH3F *)(fQA2DList->At(19)))->Fill(fProtonAnalysisBase->Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz()),
						  tpcTrack->Pt(),
						  TMath::Abs(dca3D));
	    }  
	  }//antiprotons
	  
	  //Step: kStepIsPrimary	  
	  if(fProtonAnalysisBase->IsPrimary(esd,vertex,track)) {
	    if(tpcTrack->Charge() > 0) {
	      nIdentifiedProtons += 1;
	      fProtonContainer->Fill(containerInput,kStepIsPrimary);
	    }//protons
	    else if(tpcTrack->Charge() < 0) {
	      nIdentifiedAntiProtons += 1;
	      fAntiProtonContainer->Fill(containerInput,kStepIsPrimary);
	    }//protons
	
	    
	    //Step: kStepInPhaseSpace
	    if(fProtonAnalysisBase->IsInPhaseSpace(track)) {
	      if(tpcTrack->Charge() > 0) {
		nSurvivedProtons += 1;
		fProtonContainer->Fill(containerInput,kStepInPhaseSpace);   
	      }//protons
	      else if(tpcTrack->Charge() < 0) {
		nSurvivedAntiProtons += 1;
		fAntiProtonContainer->Fill(containerInput,kStepInPhaseSpace);
	      }//antiprotons
	    }//Step: kStepInPhaseSpace
	  }//Step: kStepIsPrimary
	}//Step: kStepIdentified
      }//Step: kStepSurvived
//**************************SYSTEMATICS******************************
//FillSystematics(esd,vertex,track);
//*******************************************************************
    }//TPC only tracks
    
    else if(fProtonAnalysisBase->GetAnalysisMode() == AliProtonAnalysisBase::kGlobal) {
      if(!track) continue;
      AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
      if(!tpcTrack) continue;
      //
      gPt = track->Pt();
      gP = track->GetInnerParam()->P();

      AliExternalTrackParam cParam;
      track->RelateToVertex(vertex,
			    esd->GetMagneticField(),
			    100.,&cParam);
      track->GetImpactParameters(dcaXY,dcaZ);
      dca[0] = dcaXY; dca[1] = dcaZ;
       
      dca3D = TMath::Sqrt(TMath::Power(dca[0],2) +
			  TMath::Power(dca[1],2));
			  
      if(fProtonAnalysisBase->GetEtaMode())
	containerInput[0] = track->Eta();
      if(fProtonAnalysisBase->IsUsedMultiplicitySelection())
        containerInput[0] = nTracklets;
      else
	containerInput[0] = fProtonAnalysisBase->Rapidity(track->Px(),
							  track->Py(),
							  track->Pz());
      containerInput[1] = gPt;
      
      if(fProtonAnalysisBase->IsProton(track)) FillQA(esd,vertex,track);//QA for (anti)protons
      
      //Step: kStepSurvived
      if(fProtonAnalysisBase->IsAccepted(track)) {
	  ((TH2F *)(fQA2DList->At(1)))->Fill(gP,track->GetTPCsignal());
	  if((track->GetTPCsignal() > 0.0) && (tpcResponse.GetExpectedSignal(gP,AliPID::kProton) > 0.0))((TH2F *)(fQA2DList->At(3)))->Fill(gP,TMath::Log(track->GetTPCsignal()/tpcResponse.GetExpectedSignal(gP,AliPID::kProton)));
	  ((TH3F *)(fQA2DList->At(5)))->Fill(track->Eta(),
					     track->Phi()*180./TMath::Pi(),
					     npointsTPCdEdx);
	  ((TH3F *)(fQA2DList->At(7)))->Fill(track->Eta(),
					     track->Phi()*180./TMath::Pi(),
					     nClustersTPC);
	  ((TH3F *)(fQA2DList->At(9)))->Fill(gPt,
					     track->Phi()*180./TMath::Pi(),
					     npointsTPCdEdx);
	  ((TH3F *)(fQA2DList->At(11)))->Fill(gPt,
					      track->Phi()*180./TMath::Pi(),
					      nClustersTPC);
		
	if(track->Charge() > 0) 
	  fProtonContainer->Fill(containerInput,kStepSurvived);    

	else if(track->Charge() < 0) 
	  fAntiProtonContainer->Fill(containerInput,kStepSurvived);    

	//Step: kStepIdentified
	if(fProtonAnalysisBase->IsProton(track)) {
	    ((TH2F *)(fQA2DList->At(0)))->Fill(gP,track->GetTPCsignal());
	    if((track->GetTPCsignal() > 0.0) && (tpcResponse.GetExpectedSignal(gP,AliPID::kProton) > 0.0)) ((TH2F *)(fQA2DList->At(2)))->Fill(gP,TMath::Log(track->GetTPCsignal()/tpcResponse.GetExpectedSignal(gP,AliPID::kProton)));
	    ((TH3F *)(fQA2DList->At(4)))->Fill(track->Eta(),
					       track->Phi()*180./TMath::Pi(),
					       npointsTPCdEdx);
	    ((TH3F *)(fQA2DList->At(6)))->Fill(track->Eta(),
					       track->Phi()*180./TMath::Pi(),
					       nClustersTPC);
	    ((TH3F *)(fQA2DList->At(8)))->Fill(gPt,
					       track->Phi()*180./TMath::Pi(),
					       npointsTPCdEdx);
	    ((TH3F *)(fQA2DList->At(10)))->Fill(gPt,
						track->Phi()*180./TMath::Pi(),
						nClustersTPC);	
	
	  if(track->Charge() > 0){
	    fProtonContainer->Fill(containerInput,kStepIdentified);
	    
	      ((TH2F *)(fQA2DList->At(20)))->Fill(gP,track->GetTPCsignal());
	      ((TH2F *)(fQA2DList->At(12)))->Fill(track->Eta(),
						track->Phi()*180./TMath::Pi());
	    if(fProtonAnalysisBase->GetEtaMode()) {
	      ((TH3F *)(fQA2DList->At(14)))->Fill(track->Eta(),
						  track->Pt(),
						  dca[0]);
	      ((TH3F *)(fQA2DList->At(15)))->Fill(track->Eta(),
						  track->Pt(),
						  dca[1]);
	      ((TH3F *)(fQA2DList->At(18)))->Fill(track->Eta(),
						  track->Pt(),
						  TMath::Abs(dca3D));
	    }
	    if(fProtonAnalysisBase->IsUsedMultiplicitySelection()){
	      ((TH3F *)(fQA2DList->At(14)))->Fill(nTracklets,
						  tpcTrack->Pt(),
						  dca[0]);
	      ((TH3F *)(fQA2DList->At(15)))->Fill(nTracklets,
						  tpcTrack->Pt(),
						  dca[1]);
	      ((TH3F *)(fQA2DList->At(18)))->Fill(nTracklets,
						  tpcTrack->Pt(),
						  TMath::Abs(dca3D));

	    }
	    else {
	      ((TH3F *)(fQA2DList->At(14)))->Fill(fProtonAnalysisBase->Rapidity(track->Px(),track->Py(),track->Pz()),
						  track->Pt(),
						  dca[0]);
	      ((TH3F *)(fQA2DList->At(15)))->Fill(fProtonAnalysisBase->Rapidity(track->Px(),track->Py(),track->Pz()),
						  track->Pt(),
						  dca[1]);
	      ((TH3F *)(fQA2DList->At(18)))->Fill(fProtonAnalysisBase->Rapidity(track->Px(),track->Py(),track->Pz()),
						  track->Pt(),
						  TMath::Abs(dca3D));
	    }
	  }//protons
	       
	  else if(track->Charge() < 0){ 
	    fAntiProtonContainer->Fill(containerInput,kStepIdentified);
	      ((TH2F *)(fQA2DList->At(21)))->Fill(gP,track->GetTPCsignal());
	      ((TH2F *)(fQA2DList->At(13)))->Fill(track->Eta(),
						track->Phi()*180./TMath::Pi());
	    if(fProtonAnalysisBase->GetEtaMode()) {
	      ((TH3F *)(fQA2DList->At(16)))->Fill(track->Eta(),
						  track->Pt(),
						  dca[0]);
	      ((TH3F *)(fQA2DList->At(17)))->Fill(track->Eta(),
						  track->Pt(),
						  dca[1]);
	      ((TH3F *)(fQA2DList->At(19)))->Fill(track->Eta(),
						  track->Pt(),
						  TMath::Abs(dca3D));
	    }
	    if(fProtonAnalysisBase->IsUsedMultiplicitySelection()){
	      ((TH3F *)(fQA2DList->At(16)))->Fill(nTracklets,
						  tpcTrack->Pt(),
						  dca[0]);
	      ((TH3F *)(fQA2DList->At(17)))->Fill(nTracklets,
						  tpcTrack->Pt(),
						  dca[1]);
	      ((TH3F *)(fQA2DList->At(19)))->Fill(nTracklets,
						  tpcTrack->Pt(),
						  TMath::Abs(dca3D));

	    }
	    else {
	      ((TH3F *)(fQA2DList->At(16)))->Fill(fProtonAnalysisBase->Rapidity(track->Px(),track->Py(),track->Pz()),
						  track->Pt(),
						  dca[0]);
	      ((TH3F *)(fQA2DList->At(17)))->Fill(fProtonAnalysisBase->Rapidity(track->Px(),track->Py(),track->Pz()),
						  track->Pt(),
						  dca[1]);
	      ((TH3F *)(fQA2DList->At(19)))->Fill(fProtonAnalysisBase->Rapidity(track->Px(),track->Py(),track->Pz()),
						  track->Pt(),
						  TMath::Abs(dca3D));
	    }  
	  }//antiprotons
	  
	  //Step: kStepIsPrimary	  
	  if(fProtonAnalysisBase->IsPrimary(esd,vertex,track)) {
	    if(track->Charge() > 0) {
	      nIdentifiedProtons += 1;
	      fProtonContainer->Fill(containerInput,kStepIsPrimary);
	    }//protons
	    else if(track->Charge() < 0) {
	      nIdentifiedAntiProtons += 1;
	      fAntiProtonContainer->Fill(containerInput,kStepIsPrimary);
	    }//protons
	
	    
	    //Step: kStepInPhaseSpace
	    if(fProtonAnalysisBase->IsInPhaseSpace(track)) {
	      if(track->Charge() > 0) {
		nSurvivedProtons += 1;
		fProtonContainer->Fill(containerInput,kStepInPhaseSpace);   
	      }//protons
	      else if(track->Charge() < 0) {
		nSurvivedAntiProtons += 1;
		fAntiProtonContainer->Fill(containerInput,kStepInPhaseSpace);
	      }//antiprotons
	    }//Step: kStepInPhaseSpace
	  }//Step: kStepIsPrimary
	}//Step: kStepIdentified
      }//Step: kStepSurvived  
    }//Global tracking
  }//track loop 
  
  //delete tpcResponse;
  
  if((nIdentifiedProtons > 0)||(nIdentifiedAntiProtons > 0))
    fHistEvents->Fill(2); //number of analyzed events with at least one (anti)proton
  
  if(fProtonAnalysisBase->GetDebugMode())
    Printf("Initial number of tracks: %d | Identified (anti)protons: %d - %d | Survived (anti)protons: %d - %d",nTracks,nIdentifiedProtons,nIdentifiedAntiProtons,nSurvivedProtons,nSurvivedAntiProtons);
}

//_______________________Systematics_____________________________________________//
void AliProtonAnalysis::FillSystematics(AliESDEvent *esd, const AliESDVertex *vertex, AliESDtrack *track){

TF1  *fPtDependentDcaXY = fProtonAnalysisBase->GetPtDependentDcaXY();

Double_t containerInput[2] ;

/*const AliMultiplicity *fMult = esd->GetMultiplicity();
  if (!fMult){
	AliError("Can't get multiplicity object");
	return;
	}
  Int_t nTracklets = fMult->GetNumberOfTracklets();*/
AliCentrality *esdCentrality = esd->GetCentrality();
Float_t nTracklets = esdCentrality->GetCentralityPercentile("V0M");

//if(fProtonAnalysisBase->IsUsedITSSAMultiplicitySelection()){

/* AliESDtrackCuts *fTrackCuts = new AliESDtrackCuts();
  if(!fTrackCuts){
	AliError("Can't get track cut object");
	}
  Int_t nTracklets = fTrackCuts->GetReferenceMultiplicity(esd,AliESDtrackCuts::kTrackletsITSTPC,0.5);*/
//}

AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
      //if(!tpcTrack) continue;

//if(TMath::Abs(tpcTrack->Eta()) > 0.5) continue;//acceptance

if(fProtonAnalysisBase->IsUsedMultiplicitySelection())
        containerInput[0] = nTracklets;
      else
	containerInput[0] = fProtonAnalysisBase->Rapidity(tpcTrack->Px(),
							  tpcTrack->Py(),
							  tpcTrack->Pz());
containerInput[1]=tpcTrack->Pt();

Double_t fPSels[4];

//fPSels[0]=0.2;
fPSels[0]=7*fPtDependentDcaXY->Eval(containerInput[1]); //maximum xy DCA 0.2
//Printf("1 sigma %f",fPSels[0]);
fPSels[1]=80; // minimum TPC Clusters
fPSels[2]=3;  // N sigma
fPSels[3]=10; // z componenta of vertex

Double_t fPSelsMid[4];

//fPSelsMid[0]=0.2;
fPSelsMid[0]=7*fPtDependentDcaXY->Eval(containerInput[1]); //maximum xy DCA 0.2
//Printf("1 sigma Mid %f",fPSelsMid[0]);
fPSelsMid[1]=80; // minimum TPC Clusters 80
fPSelsMid[2]=3;  // N sigma 3
fPSelsMid[3]=10; // z componenta of vertex 10

Double_t fPSelsLow[4];

//fPSelsLow[0]=0.1;
fPSelsLow[0]=6*fPtDependentDcaXY->Eval(containerInput[1]); //maximum xy DCA 0.1 0.1
//Printf("2 sigma Low %f",fPSelsLow[0]);
fPSelsLow[1]=70; // minimum TPC Clusters 70
fPSelsLow[2]=2;  // N sigma 2
fPSelsLow[3]=8; // z componenta of vertex 8

Double_t fPSelsHigh[4];

//fPSelsHigh[0]=0.3;
fPSelsHigh[0]=8*fPtDependentDcaXY->Eval(containerInput[1]); //maximum xy DCA 0.3 0.3
//Printf("3 sigma High %f",fPSelsHigh[0]);
fPSelsHigh[1]=90; // minimum TPC Clusters 90
fPSelsHigh[2]=4;  // N sigma 4
fPSelsHigh[3]=12; // z componenta of vertex 12


for(Int_t j=0; j<4; j++){
	  for(Int_t i=0; i<4; i++){
		  if(i==j) fPSels[i] = fPSelsLow[j];
		  else fPSels[i] = fPSelsMid[i];
	  }


fProtonAnalysisBase->SetMaxDCAXY(fPSels[0]);
fProtonAnalysisBase->SetMinTPCClusters(fPSels[1]);
fProtonAnalysisBase->SetNSigma(fPSels[2]);
fProtonAnalysisBase->SetAcceptedVertexDiamond(1.,1.,fPSels[3]);


const AliESDVertex *vertex2 = fProtonAnalysisBase->GetVertex(esd,fProtonAnalysisBase->GetAnalysisMode(),fProtonAnalysisBase->GetVxMax(),fProtonAnalysisBase->GetVyMax(),fProtonAnalysisBase->GetVzMax());

if(!vertex2) continue;
if(!fProtonAnalysisBase->IsAccepted(track)) continue;
//Printf("IsAccepted");
if(!fProtonAnalysisBase->IsProton(track)) continue;
//Printf("Sigma is %i",fProtonAnalysisBase->GetNSigma());
//Printf("IsProton");
if(!fProtonAnalysisBase->IsPrimary(esd,vertex2,track)) continue;
//Printf("Sigma Low %f",fPSels[0]);
//Printf("IsPrimary");
if(!fProtonAnalysisBase->IsInPhaseSpace(track)) continue;
//Printf("IsInPhaseSpace");

 if(track->Charge() > 0) ((TH2F*)(fSysProtonLow->At(j)))->Fill(containerInput[0],containerInput[1]);
 if(track->Charge() < 0) ((TH2F*)(fSysAntiProtonLow->At(j)))->Fill(containerInput[0],containerInput[1]);
//delete &vertex2;
}
for(Int_t j=0; j<4; j++){
	  for(Int_t i=0; i<4; i++){
		  if(i==j) fPSels[i] = fPSelsHigh[j];
		  else fPSels[i] = fPSelsMid[i];
	  }

fProtonAnalysisBase->SetMaxDCAXY(fPSels[0]);
fProtonAnalysisBase->SetMinTPCClusters(fPSels[1]);
fProtonAnalysisBase->SetNSigma(fPSels[2]);
fProtonAnalysisBase->SetAcceptedVertexDiamond(1.,1.,fPSels[3]);

const AliESDVertex *vertex1 = fProtonAnalysisBase->GetVertex(esd,fProtonAnalysisBase->GetAnalysisMode(),fProtonAnalysisBase->GetVxMax(),fProtonAnalysisBase->GetVyMax(),fProtonAnalysisBase->GetVzMax());

if(!vertex1) continue;
if(!fProtonAnalysisBase->IsAccepted(track)) continue;
//Printf("IsAccepted");
if(!fProtonAnalysisBase->IsProton(track)) continue;
//Printf("Sigma is %i",fProtonAnalysisBase->GetNSigma());
//Printf("IsProton");
if(!fProtonAnalysisBase->IsPrimary(esd,vertex1,track)) continue;
//Printf("Sigma High %f",fPSels[0]);
//Printf("IsPrimary");
if(!fProtonAnalysisBase->IsInPhaseSpace(track)) continue;
//Printf("IsInPhaseSpace");

 if(track->Charge() > 0) ((TH2F*)(fSysProtonHigh->At(j)))->Fill(containerInput[0],containerInput[1]);
 if(track->Charge() < 0) ((TH2F*)(fSysAntiProtonHigh->At(j)))->Fill(containerInput[0],containerInput[1]);
// delete &vertex1;
}

fProtonAnalysisBase->SetMaxDCAXY(fPSelsMid[0]);
fProtonAnalysisBase->SetMinTPCClusters(fPSelsMid[1]);
fProtonAnalysisBase->SetNSigma(fPSelsMid[2]);
fProtonAnalysisBase->SetAcceptedVertexDiamond(1.,1.,fPSelsMid[3]);

}

//____________________________________________________________________//
void AliProtonAnalysis::Analyze(AliAODEvent* const fAOD) {
  //Main analysis part - AOD
/*  fHistEvents->Fill(1); //number of analyzed events
  Double_t containerInput[2] ;
  Int_t nTracks = fAOD->GetNumberOfTracks();
  for(Int_t iTracks = 0; iTracks < nTracks; iTracks++) {
    AliAODTrack* track = fAOD->GetTrack(iTracks);
    Double_t gPt = track->Pt();
    Double_t gP = track->P();
    
    //pid
    Double_t probability[10];
    track->GetPID(probability);
    Double_t rcc = 0.0;
    for(Int_t i = 0; i < AliPID::kSPECIESC; i++) rcc += probability[i]*fProtonAnalysisBase->GetParticleFraction(i,gP);
    if(rcc == 0.0) continue;
    Double_t w[10];
    for(Int_t i = 0; i < AliPID::kSPECIESC; i++) w[i] = probability[i]*fProtonAnalysisBase->GetParticleFraction(i,gP)/rcc;
    Long64_t fParticleType = TMath::LocMax(AliPID::kSPECIESC,w);
    if(fParticleType == 4) {
	containerInput[0] = track->Y(fParticleType);
    containerInput[1] = gPt;
      if(track->Charge() > 0) 
	fProtonContainer->Fill(containerInput,kStepInPhaseSpace);
      else if(track->Charge() < 0) 
	fAntiProtonContainer->Fill(containerInput,kStepInPhaseSpace);
    }//proton check
  }//track loop */
}

//____________________________________________________________________//
void AliProtonAnalysis::Analyze(AliStack* const stack, 
				Bool_t iInclusive) {
  //Main analysis part - MC
  fHistEvents->Fill(1); //number of analyzed events
  Double_t containerInput[2] ;
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

	if(fProtonAnalysisBase->GetEtaMode())
	containerInput[0] = particle->Eta();
      else
	containerInput[0] = fProtonAnalysisBase->Rapidity(particle->Px(),
							  particle->Py(),
							  particle->Pz());
    containerInput[1] = particle->Pt();
	  
    Int_t pdgcode = particle->GetPdgCode();
    if(pdgcode == 2212) fProtonContainer->Fill(containerInput,kStepInPhaseSpace);
    if(pdgcode == -2212) fAntiProtonContainer->Fill(containerInput,kStepInPhaseSpace);
  }//particle loop
}

//____________________________________________________________________//
Bool_t AliProtonAnalysis::PrintMean(TH1 *hist, Double_t edge) {
  //calculates the mean value of the ratio/asymmetry within \pm edge
  Double_t sum = 0.0, sumError = 0.0;
  Int_t nentries = 0;
  //calculate the mean
  for(Int_t i = 1; i <= hist->GetXaxis()->GetNbins(); i++) {
    Double_t x = hist->GetBinCenter(i);
    Double_t y = hist->GetBinContent(i);
    if(TMath::Abs(x) < edge) {
      sum += y;
      sumError += TMath::Power(hist->GetBinError(i),2);
      nentries += 1;
    }
    //Printf("eta: %lf - sum: %lf - sumError: %lf - counter: %d",
    //TMath::Abs(x),sum,sumError,nentries);
  }
  Double_t mean = 0.0;
  Double_t error = 0.0;
  if(nentries != 0) {
    mean = sum/nentries;
    error =  TMath::Sqrt(sumError)/nentries;
  }

  //calculate the error
    /*for(Int_t i = 1; i <= hist->GetXaxis()->GetNbins(); i++) {
    Double_t x = hist->GetBinCenter(i);
    Double_t y = hist->GetBinContent(i);
    if(TMath::Abs(x) < edge) {
      sum += TMath::Power((mean - y),2);
      nentries += 1;
    }
    }*/


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
void AliProtonAnalysis::SetCorrectionMapForCuts(const char* filename) {
  //Reads the file with the correction maps for the cut efficiency
  TFile *gCorrectionForCuts = TFile::Open(filename);
  if(!gCorrectionForCuts) {
    Printf("The TFile object is not valid!!!");
    return;
  }
  if(!gCorrectionForCuts->IsOpen()) {
    Printf("The file is not found!!!");
    return;
  }
  fHistYPtCorrectionForCutsProtons = dynamic_cast<TH2D *>(gCorrectionForCuts->Get("gHistCorrectionForCutsProtons"));
  fHistYPtCorrectionForCutsAntiProtons = dynamic_cast<TH2D *>(gCorrectionForCuts->Get("gHistCorrectionForCutsAntiProtons"));
  fCorrectForCutsFlag = kTRUE;
}

//____________________________________________________________________//
void AliProtonAnalysis::SetCorrectionMapForTracking(const char* filename) {
  //Reads the file with the correction maps for the tracking efficiency
  TFile *gCorrectionForTracking = TFile::Open(filename);
  if(!gCorrectionForTracking) {
    Printf("The TFile object is not valid!!!");
    return;
  }
  if(!gCorrectionForTracking->IsOpen()) {
    Printf("The file is not found!!!");
    return;
  }
  fHistYPtCorrectionForTrackingProtons = dynamic_cast<TH2D *>(gCorrectionForTracking->Get("gHistCorrectionForTrackingProtons"));
  fHistYPtCorrectionForTrackingAntiProtons = dynamic_cast<TH2D *>(gCorrectionForTracking->Get("gHistCorrectionForTrackingAntiProtons"));
  fCorrectForTrackingFlag = kTRUE;
}

//____________________________________________________________________//
void AliProtonAnalysis::SetCorrectionMapForFeedDown(const char* filename) {
  //Reads the file with the correction maps for the feed-down contamination
  TFile *gCorrectionForFeedDown = TFile::Open(filename);
  if(!gCorrectionForFeedDown) {
    Printf("The TFile object is not valid!!!");
    return;
  }
  if(!gCorrectionForFeedDown->IsOpen()) {
    Printf("The file is not found!!!");
    return;
  }
  fHistYPtCorrectionForFeedDownProtons = dynamic_cast<TH2D *>(gCorrectionForFeedDown->Get("gHistCorrectionForFeedDownProtons"));
  fHistYPtCorrectionForFeedDownAntiProtons = dynamic_cast<TH2D *>(gCorrectionForFeedDown->Get("gHistCorrectionForFeedDownAntiProtons"));
  fCorrectForFeedDownFlag = kTRUE;
}

//____________________________________________________________________//
void AliProtonAnalysis::SetCorrectionMapForSecondaries(const char* filename) {
  //Reads the file with the correction maps for the secondaries
  TFile *gCorrectionForSecondaries = TFile::Open(filename);
  if(!gCorrectionForSecondaries) {
    Printf("The TFile object is not valid!!!");
    return;
  }
  if(!gCorrectionForSecondaries->IsOpen()) {
    Printf("The file is not found!!!");
    return;
  }

  fHistYPtCorrectionForSecondaries = dynamic_cast<TH2D *>(gCorrectionForSecondaries->Get("gHistCorrectionForSecondaries"));
  fCorrectForSecondariesFlag = kTRUE;
}

//____________________________________________________________________//
void AliProtonAnalysis::SetCorrectionMapForCrossSection(const char* filename) {
  //Reads the file with the correction maps for the proper x-section
  TFile *gCorrectionForXSection = TFile::Open(filename);
  if(!gCorrectionForXSection) {
    Printf("The TFile object is not valid!!!");
    return;
  }
  if(!gCorrectionForXSection->IsOpen()) {
    Printf("The file is not found!!!");
    return;
  }

  fHistCorrectionForCrossSectionYPtProtons = dynamic_cast<TH2D *>(gCorrectionForXSection->Get("gHistCorrectionForCrossSectionProtons"));
  fHistCorrectionForCrossSectionYPtProtons->Sumw2();
  fHistCorrectionForCrossSectionYPtAntiProtons = dynamic_cast<TH2D *>(gCorrectionForXSection->Get("gHistCorrectionForCrossSectionAntiProtons"));
  fHistCorrectionForCrossSectionYPtAntiProtons->Sumw2();
  fHistCorrectionForCrossSectionFlag = kTRUE;
}

//____________________________________________________________________//
void AliProtonAnalysis::Correct() {
  //Apply the corrections: Fast & dirty way for the absorption corrections
  //Correct the protons for the efficiency
  fHistYPtProtonsCorrected = dynamic_cast<TH2D*>(fProtonContainer->ShowProjection(0,1,kStepInPhaseSpace));
  fHistYPtProtonsCorrected->Divide(fHistEfficiencyYPtProtons);
  //Correct the protons for proper cross-section
  if(fHistCorrectionForCrossSectionFlag)
    fHistYPtProtonsCorrected->Multiply(fHistCorrectionForCrossSectionYPtProtons);
  //Correct the protons for secondaries
  if(fCorrectForSecondariesFlag)
    fHistYPtProtonsCorrected->Divide(fHistYPtCorrectionForSecondaries);
  //Correct the protons for feed-down
  if(fCorrectForFeedDownFlag)
    fHistYPtProtonsCorrected->Divide(fHistYPtCorrectionForFeedDownProtons);
  //Correct the protons for the cut efficiency
  if(fCorrectForCutsFlag)
    fHistYPtProtonsCorrected->Multiply(fHistYPtCorrectionForCutsProtons);
  //Correct the protons for the tracking efficiency
  if(fCorrectForTrackingFlag)
    fHistYPtProtonsCorrected->Multiply(fHistYPtCorrectionForTrackingProtons);
 
  //Correct the antiprotons for the efficiency
  fHistYPtAntiProtonsCorrected = dynamic_cast<TH2D*>(fAntiProtonContainer->ShowProjection(0,1,kStepInPhaseSpace));
  fHistYPtAntiProtonsCorrected->Divide(fHistEfficiencyYPtAntiProtons);
  //Correct the antiprotons for proper cross-section
  if(fHistCorrectionForCrossSectionFlag)
    fHistYPtAntiProtonsCorrected->Multiply(fHistCorrectionForCrossSectionYPtAntiProtons);
  //Correct the antiprotons for feed-down
  if(fCorrectForFeedDownFlag)
    fHistYPtAntiProtonsCorrected->Divide(fHistYPtCorrectionForFeedDownAntiProtons);
  //Correct the antiprotons for the cut efficiency
   if(fCorrectForCutsFlag)
     fHistYPtAntiProtonsCorrected->Multiply(fHistYPtCorrectionForCutsAntiProtons);
  //Correct the antiprotons for the tracking efficiency
   if(fCorrectForTrackingFlag)
     fHistYPtAntiProtonsCorrected->Multiply(fHistYPtCorrectionForTrackingAntiProtons);
}

//____________________________________________________________________//
Bool_t AliProtonAnalysis::ReadCorrectionContainer(const char* filename, Int_t CS) {
  // Reads the outout of the correction framework task
  // Creates the correction maps
  // Puts the results in the different TList objects
  Bool_t status = kTRUE;
  Int_t kProton,kAntiProton;

if (CS==1) {
	kProton=2; kAntiProton=3;
//cout<<"HI HI HI HI HI HI HI HI";
}
else {kProton=0; kAntiProton=1;} 

  TFile *file = TFile::Open(filename);
  if(!file) {
    cout<<"Could not find the input CORRFW file "<<filename<<endl;
    status = kFALSE;
  }
  TList *list = dynamic_cast<TList *>(file->Get("outputList"));
  Int_t iRap = 0, iPt = 1;

  //Calculation of efficiency/correction: Protons
  AliCFContainer *gProtonContainer = dynamic_cast<AliCFContainer *>(list->At(kProton));
  AliCFEffGrid *effProtonsStep0Step2 = new AliCFEffGrid("eff20",
							"effProtonsStep0Step2",
							*gProtonContainer);
  effProtonsStep0Step2->CalculateEfficiency(1,0); 
  fHistEfficiencyYPtProtons = dynamic_cast<TH2D*>(effProtonsStep0Step2->Project(iRap,iPt));
  fHistEfficiencyYPtProtons->Sumw2();

  //Calculation of efficiency/correction: Protons
  AliCFContainer *gAntiProtonContainer = dynamic_cast<AliCFContainer *>(list->At(kAntiProton));
  AliCFEffGrid *effAntiProtonsStep0Step2 = new AliCFEffGrid("eff20",
							    "effAntiProtonsStep0Step2",
							    *gAntiProtonContainer);
  effAntiProtonsStep0Step2->CalculateEfficiency(1,0); 
  fHistEfficiencyYPtAntiProtons = dynamic_cast<TH2D*>(effAntiProtonsStep0Step2->Project(iRap,iPt));
  fHistEfficiencyYPtAntiProtons->Sumw2();

  Correct();

  return status;
}
//____________________________________________________________________//
void AliProtonAnalysis::InitSystematicsHistogram() {
//-----------------------LOW--------------------------------------
fSysProtonLow = new TList();
fSysProtonLow->SetName("fListSysProtonLow");

TH2F *gProtonMaxDCAXYLow = new TH2F("gProtonMaxDCAXYLow","Proton-MaxDCAXY-Low",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
fSysProtonLow->Add(gProtonMaxDCAXYLow);
TH2F *gProtonMinTPCClustersLow = new TH2F("gProtonMinTPCClustersLow","Proton-MinTPCClusters-Low",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
fSysProtonLow->Add(gProtonMinTPCClustersLow);
TH2F *gProtonNsigmaLow = new TH2F("gProtonNSigmaLow","Proton-Nsigma-Low",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
fSysProtonLow->Add(gProtonNsigmaLow);
TH2F *gProtonAcceptedVertexDiamondLow = new TH2F("gProtonAcceptedVertexDiamondLow","Proton-Accepted Vertex Diamond-Low",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
fSysProtonLow->Add(gProtonAcceptedVertexDiamondLow);

fSysAntiProtonLow = new TList();
fSysAntiProtonLow->SetName("fListSysAntiProtonLow");

TH2F *gAntiProtonMaxDCAXYLow = new TH2F("gAntiProtonMaxDCAXYLow","AntiProton-MaxDCAXY-Low",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
fSysAntiProtonLow->Add(gAntiProtonMaxDCAXYLow);
TH2F *gAntiProtonMinTPCClustersLow = new TH2F("gAntiProtonMinTPCClustersLow","AntiProton-MinTPCClusters-Low",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
fSysAntiProtonLow->Add(gAntiProtonMinTPCClustersLow);
TH2F *gAntiProtonNsigmaLow = new TH2F("gAntiProtonNSigmaLow","AntiProton-Nsigma-Low",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
fSysAntiProtonLow->Add(gAntiProtonNsigmaLow);
TH2F *gAntiProtonAcceptedVertexDiamondLow = new TH2F("gAntiProtonAcceptedVertexDiamondLow","AntiProton-Accepted Vertex Diamond-Low",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
fSysAntiProtonLow->Add(gAntiProtonAcceptedVertexDiamondLow);
//----------------------END LOW------------------------------------
fSysProtonHigh = new TList();
fSysProtonHigh->SetName("fListSysProtonHigh");

TH2F *gProtonMaxDCAXYHigh = new TH2F("gProtonMaxDCAXYHigh","Proton-MaxDCAXY-High",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
fSysProtonHigh->Add(gProtonMaxDCAXYHigh);
TH2F *gProtonMinTPCClustersHigh = new TH2F("gProtonMinTPCClustersHigh","Proton-MinTPCClusters-High",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
fSysProtonHigh->Add(gProtonMinTPCClustersHigh);
TH2F *gProtonNsigmaHigh = new TH2F("gProtonNSigmaHigh","Proton-Nsigma-High",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
fSysProtonHigh->Add(gProtonNsigmaHigh);
TH2F *gProtonAcceptedVertexDiamondHigh = new TH2F("gProtonAcceptedVertexDiamondHigh","Proton-Accepted Vertex Diamond-High",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
fSysProtonHigh->Add(gProtonAcceptedVertexDiamondHigh);

fSysAntiProtonHigh = new TList();
fSysAntiProtonHigh->SetName("fListSysAntiProtonHigh");

TH2F *gAntiProtonMaxDCAXYHigh = new TH2F("gAntiProtonMaxDCAXYHigh","AntiProton-MaxDCAXY-High",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
fSysAntiProtonHigh->Add(gAntiProtonMaxDCAXYHigh);
TH2F *gAntiProtonMinTPCClustersHigh = new TH2F("gAntiProtonMinTPCClustersHigh","AntiProton-MinTPCClusters-High",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
fSysAntiProtonHigh->Add(gAntiProtonMinTPCClustersHigh);
TH2F *gAntiProtonNsigmaHigh = new TH2F("gAntiProtonNSigmaHigh","AntiProton-Nsigma-High",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
fSysAntiProtonHigh->Add(gAntiProtonNsigmaHigh);
TH2F *gAntiProtonAcceptedVertexDiamondHigh = new TH2F("gAntiProtonAcceptedVertexDiamondHigh","AntiProton-Accepted Vertex Diamond-Low",fNBinsY,fMinY,fMaxY,fNBinsPt,fMinPt,fMaxPt);
fSysAntiProtonHigh->Add(gAntiProtonAcceptedVertexDiamondHigh);

}

//____________________________________________________________________//
void AliProtonAnalysis::InitQA() {
  //Applies the correction maps to the initial containers
  fGlobalQAList = new TList();
  fGlobalQAList->SetName("fGlobalQAList");

  //========================================================//
  fQA2DList = new TList();
  fQA2DList->SetName("fQA2DList");
  fGlobalQAList->Add(fQA2DList);

  //dEdx plots
  TH2F *gHistdEdxP = new TH2F("gHistdEdxP","dE/dx (TPC); P [GeV/c]; dE/dx [a.u]",1000,0.05,20.05,600,0,600);
  fQA2DList->Add(gHistdEdxP);
  TH2F *gHistProtonsdEdxP = new TH2F("gHistProtonsdEdxP","Accepted protons dE/dx (TPC); P [GeV/c]; dE/dx [a.u]",1000,0.05,20.05,600,0,600);
  fQA2DList->Add(gHistProtonsdEdxP);

  //normalized dEdx plots
  TH2F *gHistZP = new TH2F("gHistZP","Normalized dE/dx (TPC); P [GeV/c]; ln[(dE/dx)_{exp.}/(dE/dx)_{BB}] ",1000,0.05,20.05,100,-2.0,2.0);
  fQA2DList->Add(gHistZP);
  TH2F *gHistProtonsZP = new TH2F("gHistProtonsZP","Normalized dE/dx (TPC); P [GeV/c]; ln[(dE/dx)_{exp.}/(dE/dx)_{BB}] ",1000,0.05,20.05,100,-2.0,2.0);
  fQA2DList->Add(gHistProtonsZP);

  //eta-phi-Npoints(dEdx)
  TH3F *gHistEtaPhiTPCdEdxNPoints = new TH3F("gHistEtaPhiTPCdEdxNPoints",
					     ";#eta;#phi;N_{points}(TPC)",
					     18,-0.9,0.9,
					     180,0,360,
					     100,0,200);
  gHistEtaPhiTPCdEdxNPoints->SetStats(kTRUE);
  fQA2DList->Add(gHistEtaPhiTPCdEdxNPoints);
  TH3F *gHistProtonsEtaPhiTPCdEdxNPoints = new TH3F("gHistProtonsEtaPhiTPCdEdxNPoints",
						    ";#eta;#phi;N_{points}(TPC)",
						    18,-0.9,0.9,
						    180,0,360,
						    100,0,200);
  gHistProtonsEtaPhiTPCdEdxNPoints->SetStats(kTRUE);
  fQA2DList->Add(gHistProtonsEtaPhiTPCdEdxNPoints);

  //eta-phi-Npoints
  TH3F *gHistEtaPhiTPCNPoints = new TH3F("gHistEtaPhiTPCNPoints",
					 ";#eta;#phi;N_{points}(TPC)",
					 18,-0.9,0.9,
					 180,0,360,
					 100,0,200);
  gHistEtaPhiTPCNPoints->SetStats(kTRUE);
  fQA2DList->Add(gHistEtaPhiTPCNPoints);
  TH3F *gHistProtonsEtaPhiTPCNPoints = new TH3F("gHistProtonsEtaPhiTPCNPoints",
						";#eta;#phi;N_{points}(TPC)",
						18,-0.9,0.9,
						180,0,360,
						100,0,200);
  gHistProtonsEtaPhiTPCNPoints->SetStats(kTRUE);
  fQA2DList->Add(gHistProtonsEtaPhiTPCNPoints);

  //pt-phi-Npoints(dEdx)
  TH3F *gHistPtPhiTPCdEdxNPoints = new TH3F("gHistPtPhiTPCdEdxNPoints",
					    ";P_{T} [GeV/c];#phi;N_{points}(TPC)",
					    fNBinsPt,fMinPt,fMaxPt,
					    180,0,360,
					    100,0,200);
  gHistPtPhiTPCdEdxNPoints->SetStats(kTRUE);
  fQA2DList->Add(gHistPtPhiTPCdEdxNPoints);
  TH3F *gHistProtonsPtPhiTPCdEdxNPoints = new TH3F("gHistProtonsPtPhiTPCdEdxNPoints",
						    ";P_{T} [GeV/c];#phi;N_{points}(TPC)",
						    fNBinsPt,fMinPt,fMaxPt,
						    180,0,360,
						    100,0,200);
  gHistProtonsPtPhiTPCdEdxNPoints->SetStats(kTRUE);
  fQA2DList->Add(gHistProtonsPtPhiTPCdEdxNPoints);

  //pt-phi-Npoints
  TH3F *gHistPtPhiTPCNPoints = new TH3F("gHistPtPhiTPCNPoints",
					";P_{T} [GeV/c];#phi;N_{points}(TPC)",
					fNBinsPt,fMinPt,fMaxPt,
					180,0,360,
					100,0,200);
  gHistPtPhiTPCNPoints->SetStats(kTRUE);
  fQA2DList->Add(gHistPtPhiTPCNPoints);
  TH3F *gHistProtonsPtPhiTPCNPoints = new TH3F("gHistProtonsPtPhiTPCNPoints",
					       ";P_{T} [GeV/c];#phi;N_{points}(TPC)",
					       fNBinsPt,fMinPt,fMaxPt,
					       180,0,360,
					       100,0,200);
  gHistProtonsPtPhiTPCNPoints->SetStats(kTRUE);
  fQA2DList->Add(gHistProtonsPtPhiTPCNPoints);

  //eta-phi for protons & antiprotons
  TH2F *gHistProtonsEtaPhi = new TH2F("gHistProtonsEtaPhi",
				      ";#eta;#phi",
				      18,-0.9,0.9,
				      180,0,360);
  gHistProtonsEtaPhi->SetStats(kTRUE);
  fQA2DList->Add(gHistProtonsEtaPhi);
  TH2F *gHistAntiProtonsEtaPhi = new TH2F("gHistAntiProtonsEtaPhi",
					  ";#eta;#phi",
					  18,-0.9,0.9,
					  180,0,360);
  gHistAntiProtonsEtaPhi->SetStats(kTRUE);
  fQA2DList->Add(gHistAntiProtonsEtaPhi);

  const Int_t nBinsdca = 10000;
  Double_t dcamin = -10., dcamax = 10.;
  
  //dca vs pT for protons & antiprotons
  TH3F *gHistProtonsDCAxyEtaPt = new TH3F("gHistProtonsDCAxyEtaPt",
					  ";P_{T} [GeV/c];dca_{xy} [cm]",
					  fNBinsY,fMinY,fMaxY,
					  fNBinsPt,fMinPt,fMaxPt,
					  nBinsdca,dcamin, dcamax);
  if(fProtonAnalysisBase->GetEtaMode())
    gHistProtonsDCAxyEtaPt->GetXaxis()->SetTitle("#eta");
  else
    gHistProtonsDCAxyEtaPt->GetXaxis()->SetTitle("y");
  gHistProtonsDCAxyEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistProtonsDCAxyEtaPt);
  TH3F *gHistProtonsDCAzEtaPt = new TH3F("gHistProtonsDCAzEtaPt",
					 ";P_{T} [GeV/c];dca_{z} [cm]",
					 fNBinsY,fMinY,fMaxY,
					 fNBinsPt,fMinPt,fMaxPt,
					 nBinsdca,dcamin, dcamax);
  if(fProtonAnalysisBase->GetEtaMode())
    gHistProtonsDCAzEtaPt->GetXaxis()->SetTitle("#eta");
  else
    gHistProtonsDCAzEtaPt->GetXaxis()->SetTitle("y");
  gHistProtonsDCAzEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistProtonsDCAzEtaPt);
  TH3F *gHistAntiProtonsDCAxyEtaPt = new TH3F("gHistAntiProtonsDCAxyEtaPt",
					      ";P_{T} [GeV/c];dca_{xy} [cm]",
					      fNBinsY,fMinY,fMaxY,
					      fNBinsPt,fMinPt,fMaxPt,
					      nBinsdca,dcamin, dcamax);
  if(fProtonAnalysisBase->GetEtaMode())
    gHistAntiProtonsDCAxyEtaPt->GetXaxis()->SetTitle("#eta");
  else
    gHistAntiProtonsDCAxyEtaPt->GetXaxis()->SetTitle("y");
  gHistAntiProtonsDCAxyEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistAntiProtonsDCAxyEtaPt);
  TH3F *gHistAntiProtonsDCAzEtaPt = new TH3F("gHistAntiProtonsDCAzEtaPt",
					     ";P_{T} [GeV/c];dca_{z} [cm]",
					     fNBinsY,fMinY,fMaxY,
					     fNBinsPt,fMinPt,fMaxPt,
					     nBinsdca,dcamin, dcamax);
  if(fProtonAnalysisBase->GetEtaMode())
    gHistAntiProtonsDCAzEtaPt->GetXaxis()->SetTitle("#eta");
  else
    gHistAntiProtonsDCAzEtaPt->GetXaxis()->SetTitle("y");
  gHistAntiProtonsDCAzEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistAntiProtonsDCAzEtaPt);

  TH3F *gHistProtonsDCA3DEtaPt = new TH3F("gHistProtonsDCA3DEtaPt",
					  ";P_{T} [GeV/c];dca [cm]",
					  fNBinsY,fMinY,fMaxY,
					  fNBinsPt,fMinPt,fMaxPt,
					  100,0,20);
  if(fProtonAnalysisBase->GetEtaMode())
    gHistProtonsDCA3DEtaPt->GetXaxis()->SetTitle("#eta");
  else
    gHistProtonsDCA3DEtaPt->GetXaxis()->SetTitle("y");
  gHistProtonsDCA3DEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistProtonsDCA3DEtaPt);
  TH3F *gHistAntiProtonsDCA3DEtaPt = new TH3F("gHistAntiProtonsDCA3DEtaPt",
					      ";P_{T} [GeV/c];dca [cm]",
					      fNBinsY,fMinY,fMaxY,
					      fNBinsPt,fMinPt,fMaxPt,
					      100,0,20);
  if(fProtonAnalysisBase->GetEtaMode())
    gHistAntiProtonsDCA3DEtaPt->GetXaxis()->SetTitle("#eta");
  else
    gHistAntiProtonsDCA3DEtaPt->GetXaxis()->SetTitle("y");
  gHistAntiProtonsDCA3DEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistAntiProtonsDCA3DEtaPt);

  TH2F *gHistPosdEdxP = new TH2F("gHistPosdEdxP","dE/dx (TPC); P [GeV/c]; dE/dx [a.u]",1000,0.05,20.05,600,0,600);
  fQA2DList->Add(gHistPosdEdxP);
  TH2F *gHistNegdEdxP = new TH2F("gHistNegdEdxP","dE/dx (TPC); P [GeV/c]; dE/dx [a.u]",1000,0.05,20.05,600,0,600);
  fQA2DList->Add(gHistNegdEdxP);

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
					     100,0.,10.);
  fQAProtonsAcceptedList->Add(gProtonsSigmaToVertexPass);
  TH1F *gProtonsSigmaToVertexTPCPass = new TH1F("gProtonsSigmaToVertexTPCPass",
						";#sigma_{Vertex};Entries",
						100,0.,10.);
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
  TH1F *gProtonsITSClusterMapPass = new TH1F("gProtonsITSClusterMapPass",";ITS Layer;Entries",6,0.5,6.5);
  fQAProtonsAcceptedList->Add(gProtonsITSClusterMapPass);
  TH1F *gProtonsDCA3DPass = new TH1F("gProtonsDCA3DPass",
				     ";dca [cm];Entries",
				     100,0,20);
  fQAProtonsAcceptedList->Add(gProtonsDCA3DPass);

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
  TH1F *gProtonsITSClusterMapReject = new TH1F("gProtonsITSClusterMapReject",";ITS Layer;Entries",6,0.5,6.5);
  gProtonsITSClusterMapReject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsITSClusterMapReject);
  TH1F *gProtonsDCA3DReject = new TH1F("gProtonsDCA3DReject",
				       ";dca [cm];Entries",
				       100,0,20);
  gProtonsDCA3DReject->SetFillColor(kRed-2);
  fQAProtonsRejectedList->Add(gProtonsDCA3DReject);
    
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
  TH1F *gAntiProtonsITSClusterMapPass = new TH1F("gAntiProtonsITSClusterMapPass",";ITS Layer;Entries",6,0.5,6.5);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsITSClusterMapPass);
  TH1F *gAntiProtonsDCA3DPass = new TH1F("gAntiProtonsDCA3DPass",
					 ";dca [cm];Entries",
					 100,0,20);
  fQAAntiProtonsAcceptedList->Add(gAntiProtonsDCA3DPass);

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
  TH1F *gAntiProtonsITSClusterMapReject = new TH1F("gAntiProtonsITSClusterMapReject",";ITS Layer;Entries",6,0.5,6.5);
  gAntiProtonsITSClusterMapReject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsITSClusterMapReject);
  TH1F *gAntiProtonsDCA3DReject = new TH1F("gAntiProtonsDCA3DReject",
					   ";dca [cm];Entries",
					   100,0,20);
  gAntiProtonsDCA3DReject->SetFillColor(kRed-2);
  fQAAntiProtonsRejectedList->Add(gAntiProtonsDCA3DReject);
}

//____________________________________________________________________//
void AliProtonAnalysis::FillQA(AliESDEvent *esd,
			       const AliESDVertex *vertex, 
			       AliESDtrack* track) {
  //Fills the QA histograms
  Double_t gPt = 0.0, gPx = 0.0, gPy = 0.0, gPz = 0.0;
  Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};  //The impact parameters and their covariance.
  Double_t dca3D = 0.0;
  Float_t dcaXY = 0.0, dcaZ = 0.0;

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
  if(fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kFullHybrid) {
    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    AliExternalTrackParam cParam;
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
      track->RelateToVertex(vertex,
			    esd->GetMagneticField(),
			    100.,&cParam);
      track->GetImpactParameters(dcaXY,dcaZ);
      dca[0] = dcaXY; dca[1] = dcaZ;
    }
  }
  else{
    gPt = track->Pt();
    gPx = track->Px();
    gPy = track->Py();
    gPz = track->Pz();
    AliExternalTrackParam cParam;
    track->RelateToVertex(vertex,
			    esd->GetMagneticField(),
			    100.,&cParam);
    track->GetImpactParameters(dcaXY,dcaZ);
    dca[0] = dcaXY; dca[1] = dcaZ;
  }

  dca3D = TMath::Sqrt(TMath::Power(dca[0],2) +
		      TMath::Power(dca[1],2));

  
  Int_t nClustersITS = track->GetITSclusters(0x0);
  Int_t nClustersTPC = track->GetTPCclusters(0x0);

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
    if(fProtonAnalysisBase->IsUsedPointOnSPDLayer()) {
      if((!track->HasPointOnITSLayer(0))&&(!track->HasPointOnITSLayer(1))) {
	for(Int_t iLayer = 0; iLayer < 6; iLayer++) {
	  if(track->HasPointOnITSLayer(iLayer))
	    ((TH1F *)(fQAProtonsRejectedList->At(27)))->Fill(iLayer+1);
	}
      }
      else if((track->HasPointOnITSLayer(0))||(track->HasPointOnITSLayer(1))) {
	for(Int_t iLayer = 0; iLayer < 6; iLayer++) {
	  if(track->HasPointOnITSLayer(iLayer))
	    ((TH1F *)(fQAProtonsAcceptedList->At(27)))->Fill(iLayer+1);
	}
      }
    }//point on either SPD layers
    if(fProtonAnalysisBase->IsUsedMaxDCA3D()) {
      if(dca3D > fProtonAnalysisBase->GetMaxDCA3D()) {
	((TH1F *)(fQAProtonsRejectedList->At(28)))->Fill(dca3D);
      }
      if(dca3D < fProtonAnalysisBase->GetMaxDCA3D()) 
	((TH1F *)(fQAProtonsAcceptedList->At(28)))->Fill(dca3D);
    }//dca3D
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
    if(fProtonAnalysisBase->IsUsedPointOnSPDLayer()) {
      if((!track->HasPointOnITSLayer(0))&&(!track->HasPointOnITSLayer(1))) {
	for(Int_t iLayer = 0; iLayer < 6; iLayer++) {
	  if(track->HasPointOnITSLayer(iLayer))
	    ((TH1F *)(fQAAntiProtonsRejectedList->At(27)))->Fill(iLayer+1);
	}
      }
      else if((track->HasPointOnITSLayer(0))||(track->HasPointOnITSLayer(1))) {
	for(Int_t iLayer = 0; iLayer < 6; iLayer++) {
	  if(track->HasPointOnITSLayer(iLayer))
	    ((TH1F *)(fQAAntiProtonsAcceptedList->At(27)))->Fill(iLayer+1);
	}
      }
    }//point on either SPD layers
    if(fProtonAnalysisBase->IsUsedMaxDCA3D()) {
      if(dca3D > fProtonAnalysisBase->GetMaxDCA3D()) {
	((TH1F *)(fQAAntiProtonsRejectedList->At(28)))->Fill(dca3D);
      }
      if(dca3D < fProtonAnalysisBase->GetMaxDCA3D()) 
	((TH1F *)(fQAAntiProtonsAcceptedList->At(28)))->Fill(dca3D);
    }//dca3D
  }//antiprotons
}
