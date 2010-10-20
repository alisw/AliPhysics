// **************************************
// Task used for the correction of determiantion of reconstructed jet spectra
// Compares input (gen) and output (rec) jets   
// *******************************************


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

 
#include <TROOT.h>
#include <TRandom.h>
#include <TString.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TKey.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include "TDatabasePDG.h"

#include "AliAnalysisTaskJetServices.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisManager.h"
#include "AliJetFinder.h"
#include "AliJetHeader.h"
#include "AliJetReader.h"
#include "AliJetReaderHeader.h"
#include "AliUA1JetHeaderV1.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODTrack.h"
#include "AliAODJet.h"
#include "AliAODMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenPythiaEventHeader.h"
#include "AliJetKineReaderHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliPhysicsSelection.h"
#include "AliTriggerAnalysis.h"

#include "AliAnalysisHelperJetTasks.h"

ClassImp(AliAnalysisTaskJetServices)

AliAnalysisTaskJetServices::AliAnalysisTaskJetServices(): AliAnalysisTaskSE(),
  fUseAODInput(kFALSE),
  fUsePhysicsSelection(kFALSE),
  fMC(kFALSE),
  fSelectionInfoESD(0),
  fEventCutInfoESD(0),
  fAvgTrials(1),
  fVtxXMean(0),
  fVtxYMean(0),
  fVtxZMean(0),
  fVtxRCut(1.),
  fVtxZCut(8.),
  fPtMinCosmic(5.),
  fRIsolMinCosmic(3.),
  fMaxCosmicAngle(0.01),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardTrials(0x0),
  fh1SelectionInfoESD(0x0),
  fh1EventCutInfoESD(0),
  fh2TriggerCount(0x0),
  fh2ESDTriggerCount(0x0),
  fh2TriggerVtx(0x0),
  fh2ESDTriggerVtx(0x0),
  fh2ESDTriggerRun(0x0),
  fh2VtxXY(0x0),
  fh1NCosmicsPerEvent(0x0),
  fTriggerAnalysis(0x0),
  fHistList(0x0)  
{
  fRunRange[0] = fRunRange[1] = 0; 
}

AliAnalysisTaskJetServices::AliAnalysisTaskJetServices(const char* name):
  AliAnalysisTaskSE(name),
  fUseAODInput(kFALSE),
  fUsePhysicsSelection(kFALSE),
  fMC(kFALSE),
  fSelectionInfoESD(0),
  fEventCutInfoESD(0),
  fAvgTrials(1),
  fVtxXMean(0),
  fVtxYMean(0),
  fVtxZMean(0),
  fVtxRCut(1.),
  fVtxZCut(8.),
  fPtMinCosmic(5.),
  fRIsolMinCosmic(3.),
  fMaxCosmicAngle(0.01),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardTrials(0x0),
  fh1SelectionInfoESD(0x0),
  fh1EventCutInfoESD(0),
  fh2TriggerCount(0x0),
  fh2ESDTriggerCount(0x0),
  fh2TriggerVtx(0x0),
  fh2ESDTriggerVtx(0x0),
  fh2ESDTriggerRun(0x0),
  fh2VtxXY(0x0),
  fh1NCosmicsPerEvent(0x0),
  fTriggerAnalysis(0x0),
  fHistList(0x0)  
{
  fRunRange[0] = fRunRange[1] = 0; 
  DefineOutput(1,TList::Class());
}



Bool_t AliAnalysisTaskJetServices::Notify()
{
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // 

  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  Float_t xsection = 0;
  Float_t ftrials  = 1;

  fAvgTrials = 1;
  if(tree){
    TFile *curfile = tree->GetCurrentFile();
    if (!curfile) {
      Error("Notify","No current file");
      return kFALSE;
    }
    if(!fh1Xsec||!fh1Trials){
      Printf("%s%d No Histogram fh1Xsec",(char*)__FILE__,__LINE__);
      return kFALSE;
    }
    AliAnalysisHelperJetTasks::PythiaInfoFromFile(curfile->GetName(),xsection,ftrials);
    fh1Xsec->Fill("<#sigma>",xsection);
    // construct a poor man average trials 
    Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
    if(ftrials>=nEntries && nEntries>0.)fAvgTrials = ftrials/nEntries;
  }  
  return kTRUE;
}

void AliAnalysisTaskJetServices::UserCreateOutputObjects()
{

  //
  // Create the output container
  //


  if (fDebug > 1) printf("AnalysisTaskJetServices::UserCreateOutputObjects() \n");

  OpenFile(1);
  if(!fHistList)fHistList = new TList();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  fh1Xsec = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
  fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
  fHistList->Add(fh1Xsec);

  fh1Trials = new TH1F("fh1Trials","trials root file",1,0,1);
  fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
  fHistList->Add(fh1Trials);

  const Int_t nBinPt = 100;
  Double_t binLimitsPt[nBinPt+1];
  for(Int_t iPt = 0;iPt <= nBinPt;iPt++){
    if(iPt == 0){
      binLimitsPt[iPt] = 0.0;
    }
    else {// 1.0
      binLimitsPt[iPt] =  binLimitsPt[iPt-1] + 2.5;
    }
  }
  
  fh2TriggerCount = new TH2F("fh2TriggerCount",";Trigger No.;constrained;Count",AliAnalysisHelperJetTasks::kTrigger,-0.5,AliAnalysisHelperJetTasks::kTrigger-0.5,kConstraints,-0.5,kConstraints-0.5); 
  fHistList->Add(fh2TriggerCount);

  fh2ESDTriggerCount = new TH2F("fh2ESDTriggerCount",";Trigger No.;constrained;Count",AliAnalysisHelperJetTasks::kTrigger,-0.5,AliAnalysisHelperJetTasks::kTrigger-0.5,kConstraints,-0.5,kConstraints-0.5); 
  fHistList->Add(fh2ESDTriggerCount);
  const Int_t nBins = AliAnalysisHelperJetTasks::kTrigger*kConstraints;
  fh2TriggerVtx = new TH2F("fh2TriggerVtx",";Constraint No. * (trig no+1);Vtx (cm);Count",nBins,-0.5,nBins-0.5,400,-20.,20.); 
  fHistList->Add(fh2TriggerVtx);

  fh2ESDTriggerVtx = new TH2F("fh2ESDTriggerVtx",";Constraint No.* (trg no+1);Vtx (cm);Count",nBins,-0.5,nBins-0.5,400,-20.,20.); 
  fHistList->Add(fh2ESDTriggerVtx);
  

  fh1PtHard = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",nBinPt,binLimitsPt);
  fHistList->Add(fh1PtHard);
  fh1PtHardTrials = new TH1F("fh1PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard}",nBinPt,binLimitsPt);
  fHistList->Add(fh1PtHardTrials);
  fh1SelectionInfoESD = new TH1F("fh1SelectionInfoESD","Bit Masks that satisfy the selection info",
				 AliAnalysisHelperJetTasks::kTotalSelections,0.5,AliAnalysisHelperJetTasks::kTotalSelections+0.5);
  fHistList->Add(fh1SelectionInfoESD);

  fh1EventCutInfoESD = new TH1F("fh1EventCutInfoESD","Bit Masks that for the events after each step of cuts",
				kTotalEventCuts,0.5,kTotalEventCuts+0.5);
  fHistList->Add(fh1EventCutInfoESD);

  // 3 decisions, 0 trigger X, X + SPD vertex, X + SPD vertex in range  
  // 3 triggers BB BE/EB EE

  fh2ESDTriggerRun = new TH2F("fh2ESDTriggerRun","Trigger vs run number:run;trigger",(Int_t)(1+fRunRange[1]-fRunRange[0]),fRunRange[0]-0.5,fRunRange[1]+0.5,10,-0.5,9.5);
  fHistList->Add(fh2ESDTriggerRun);

  fh2VtxXY = new TH2F("fh2VtxXY","Beam Spot all INT triggered events;x (cm);y (cm)",160,-10,10,160,-10,10);
  fHistList->Add(fh2VtxXY);
  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fHistList->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fHistList->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fHistList->At(i));
    if(hn)hn->Sumw2();
  }

  fh1NCosmicsPerEvent = new TH1F("fh1NCosmicsPerEvent","Number of cosmic candidates per event",10,0.,10.);
  fHistList->Add(fh1NCosmicsPerEvent),


  TH1::AddDirectory(oldStatus);
}

void AliAnalysisTaskJetServices::Init()
{
  //
  // Initialization
  //
  if (fDebug > 1) printf("AnalysisTaskJetServices::Init() \n");

}

void AliAnalysisTaskJetServices::UserExec(Option_t */*option*/)
{

  //
  // Execute analysis for current event
  //
 
  AliAODEvent *aod = 0;
  AliESDEvent *esd = 0;
  
  AliAnalysisHelperJetTasks::Selected(kTRUE,kFALSE); // set slection to false
  fSelectionInfoESD = 0; // reset
  fEventCutInfoESD = 0; // reset
  AliAnalysisHelperJetTasks::SelectInfo(kTRUE,fSelectionInfoESD); // set slection to false


  if(fUseAODInput){    
    aod = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!aod){
      Printf("%s:%d AODEvent not found in Input Manager %d",(char*)__FILE__,__LINE__,fUseAODInput);
      return;
    }
    // fethc the header
  }
  else{
    //  assume that the AOD is in the general output...
    aod  = AODEvent();
    if(!aod){
      Printf("%s:%d AODEvent not found in the Output",(char*)__FILE__,__LINE__);
      return;
    }
    esd = dynamic_cast<AliESDEvent*>(InputEvent());
  }
  
  fSelectionInfoESD |= kNoEventCut;
  fEventCutInfoESD |= kNoEventCut;

  Bool_t esdVtxValid = false;
  Bool_t esdVtxIn = false;
  Bool_t aodVtxValid = false;
  Bool_t aodVtxIn = false;

  if(esd){
    // trigger analyisis
    if(!fTriggerAnalysis){
      fTriggerAnalysis = new AliTriggerAnalysis;
      fTriggerAnalysis->SetAnalyzeMC(fMC);
      fTriggerAnalysis->SetSPDGFOThreshhold(1);
    }
    //    fTriggerAnalysis->FillTriggerClasses(esd);
    Bool_t v0A       = fTriggerAnalysis->IsOfflineTriggerFired(esd, AliTriggerAnalysis::kV0A);
    Bool_t v0C       = fTriggerAnalysis->IsOfflineTriggerFired(esd, AliTriggerAnalysis::kV0C);
    Bool_t v0ABG = fTriggerAnalysis->IsOfflineTriggerFired(esd, AliTriggerAnalysis::kV0ABG);
    Bool_t v0CBG = fTriggerAnalysis->IsOfflineTriggerFired(esd, AliTriggerAnalysis::kV0CBG);
    Bool_t spdFO      = fTriggerAnalysis->SPDFiredChips(esd, 0);;
    if(v0A)fSelectionInfoESD |=  AliAnalysisHelperJetTasks::kV0A;
    if(v0C)fSelectionInfoESD |=  AliAnalysisHelperJetTasks::kV0C;
    if(!(v0ABG||v0CBG))fSelectionInfoESD |=  AliAnalysisHelperJetTasks::kNoV0BG;
    if(spdFO)fSelectionInfoESD |=  AliAnalysisHelperJetTasks::kSPDFO;

    Float_t run = (Float_t)esd->GetRunNumber();
    const AliESDVertex *vtxESD = esd->GetPrimaryVertex();
    esdVtxValid = 
    esdVtxIn = IsVertexIn(vtxESD);
    Float_t zvtx = -999;
    Float_t xvtx = -999;
    Float_t yvtx = -999;

    if(esdVtxValid){
      zvtx = vtxESD->GetZ();
      yvtx = vtxESD->GetY();
      xvtx = vtxESD->GetX();
    }

    // CKB this can be cleaned up a bit...
    
    Int_t iTrig = -1;
    if(esd->GetFiredTriggerClasses().Contains("CINT1B")
       ||esd->GetFiredTriggerClasses().Contains("CSMBB")
       ||esd->GetFiredTriggerClasses().Contains("MB1")
       ||esd->GetFiredTriggerClasses().Contains("CINT6B")){
      iTrig = 0;
      fSelectionInfoESD |=  AliAnalysisHelperJetTasks::kBunchBunch;
    }
    else if(esd->GetFiredTriggerClasses().Contains("CINT1A")
	    ||esd->GetFiredTriggerClasses().Contains("CSMBA")
	    ||esd->GetFiredTriggerClasses().Contains("CINT6A")
	    ||esd->GetFiredTriggerClasses().Contains("CINT1C")
	    ||esd->GetFiredTriggerClasses().Contains("CSMBC")
	    ||esd->GetFiredTriggerClasses().Contains("CINT6C")){
      // empty bunch or bunch empty
      iTrig = 1;
      fSelectionInfoESD |=  AliAnalysisHelperJetTasks::kBunchEmpty;
    }
    else if(esd->GetFiredTriggerClasses().Contains("CINT1-E")
       ||esd->GetFiredTriggerClasses().Contains("CINT6-E")){
      iTrig = 2;
      fSelectionInfoESD |=  AliAnalysisHelperJetTasks::kEmptyEmpty;
    }

    
    if(iTrig>=0){
      iTrig *= 3;
      fh2ESDTriggerRun->Fill(run,iTrig+1);
      if(vtxESD->GetNContributors()>2){
	fh2ESDTriggerRun->Fill(run,iTrig+2);
	fh2VtxXY->Fill(xvtx,yvtx);
      }
      xvtx -= fVtxXMean; 
      yvtx -= fVtxYMean; 
      zvtx -= fVtxZMean; 
      Float_t r2 = xvtx *xvtx + yvtx *yvtx; 
      if(TMath::Abs(zvtx)<fVtxZCut&&r2<(fVtxRCut*fVtxRCut))fh2ESDTriggerRun->Fill(run,iTrig+3);
    }
    else{
      fh2ESDTriggerRun->Fill(run,0);
    }
    // BKC
  }
  

  // Apply additional constraints
  Bool_t esdEventSelected = IsEventSelected(esd);
  Bool_t esdEventPileUp = IsEventPileUp(esd);
  Bool_t esdEventCosmic = IsEventCosmic(esd);

  Bool_t aodEventSelected = IsEventSelected(aod);

  Bool_t physicsSelection = ((fInputHandler->IsEventSelected())&AliVEvent::kMB);
  fEventCutInfoESD |= kPhysicsSelectionCut; // other alreay set via IsEventSelected
  fh1EventCutInfoESD->Fill(fEventCutInfoESD);

  if(esdEventSelected) fSelectionInfoESD |=  AliAnalysisHelperJetTasks::kVertexIn;
  if(esdEventPileUp)   fSelectionInfoESD |=  AliAnalysisHelperJetTasks::kIsPileUp;
  if(esdEventCosmic)   fSelectionInfoESD |=  AliAnalysisHelperJetTasks::kIsCosmic;
  if(physicsSelection) fSelectionInfoESD |=  AliAnalysisHelperJetTasks::kPhysicsSelection;


  // here we have all selection information, fill histogram
  for(unsigned int i = 1;i<(UInt_t)fh1SelectionInfoESD->GetNbinsX();i++){
    if((i&fSelectionInfoESD)==i)fh1SelectionInfoESD->Fill(i);
  }
  AliAnalysisHelperJetTasks::SelectInfo(kTRUE,fSelectionInfoESD); 

  if(esd&&aod&&fDebug){
    if(esdEventSelected&&!aodEventSelected){
      Printf("%s:%d Different Selection for ESD and AOD",(char*)__FILE__,__LINE__);
      const AliESDVertex *vtxESD = esd->GetPrimaryVertex();
      const AliAODVertex *vtxAOD = aod->GetPrimaryVertex();
      Printf("ESD Vtx %s %s %d",vtxESD->GetName(),vtxESD->GetTitle(),vtxESD->GetNContributors());
      vtxESD->Print();
      Printf("AOD Vtx %s %s %d",vtxAOD->GetName(),vtxAOD->GetTitle(),vtxAOD->GetNContributors());
      vtxAOD->Print();
    }
  }

  // loop over all possible triggers for esd 

  if(esd){
    const AliESDVertex *vtxESD = esd->GetPrimaryVertex();
    //      Printf(">> ESDvtx %s %s",vtxESD->GetName(),vtxESD->GetTitle());vtxESD->Print();
    TString vtxName(vtxESD->GetName());
    Float_t zvtx = vtxESD->GetZ();
    for(int it = AliAnalysisHelperJetTasks::kAcceptAll;it < AliAnalysisHelperJetTasks::kTrigger;it++){
      Bool_t esdTrig = kFALSE;
      esdTrig = AliAnalysisHelperJetTasks::IsTriggerFired(esd,(AliAnalysisHelperJetTasks::Trigger)it);
      if(esdTrig)fh2ESDTriggerCount->Fill(it,kAllTriggered);
      Bool_t cand = physicsSelection;
      if(cand){
	fh2ESDTriggerCount->Fill(it,kSelectedALICE); 
	fh2ESDTriggerVtx->Fill(kSelectedALICE*(it+1),zvtx);
      }
      if(!fUsePhysicsSelection)cand =  AliAnalysisHelperJetTasks::IsTriggerFired(esd,AliAnalysisHelperJetTasks::kMB1);
      if(vtxESD->GetNContributors()>2&&!vtxName.Contains("TPCVertex")){
	if(esdTrig){
	  fh2ESDTriggerCount->Fill(it,kTriggeredVertex);
	  fh2ESDTriggerVtx->Fill(kTriggeredVertex*(it+1),zvtx);
	}
	if(esdEventSelected&&esdTrig){
	  fh2ESDTriggerCount->Fill(it,kTriggeredVertexIn);
	  fh2ESDTriggerVtx->Fill(kTriggeredVertexIn*(it+1),zvtx);
	}
	if(cand){
	  fh2ESDTriggerCount->Fill(it,kSelectedALICEVertexValid);
	  fh2ESDTriggerVtx->Fill(kSelectedALICEVertexValid*(it+1),zvtx);
	}
      }
      if(cand&&esdEventSelected){
	fh2ESDTriggerCount->Fill(it,kSelectedALICEVertexIn);
	fh2ESDTriggerVtx->Fill(kSelectedALICEVertexIn*(it+1),zvtx);
	fh2ESDTriggerVtx->Fill(kSelected*(it+1),zvtx);
	fh2ESDTriggerCount->Fill(it,kSelected);
	AliAnalysisHelperJetTasks::Selected(kTRUE,kTRUE);// select this event
      }
    }
  }

  if(aod){
    const AliAODVertex *vtxAOD = aod->GetPrimaryVertex();
    aodVtxValid = IsVertexValid(vtxAOD);
    aodVtxIn = IsVertexIn(vtxAOD);

    for(int it = AliAnalysisHelperJetTasks::kAcceptAll;it < AliAnalysisHelperJetTasks::kTrigger;it++){
      Bool_t aodTrig = kFALSE;
      aodTrig = AliAnalysisHelperJetTasks::IsTriggerFired(aod,(AliAnalysisHelperJetTasks::Trigger)it);
      Bool_t cand = aod->GetHeader()->GetOfflineTrigger()&AliVEvent::kMB;
      if(aodTrig)fh2TriggerCount->Fill(it,kAllTriggered);
      if(aodVtxValid){
	Float_t zvtx = vtxAOD->GetZ();
	if(aodTrig){
	  fh2TriggerCount->Fill(it,kTriggeredVertex);
	  fh2TriggerVtx->Fill(kTriggeredVertex*(it+1),zvtx);
	}
	if(cand&&aodEventSelected){
	  fh2TriggerCount->Fill(it,kSelected);
	}
	if(aodTrig&&aodEventSelected){
	  fh2TriggerVtx->Fill(kTriggeredVertexIn*(it+1),zvtx);
	  fh2TriggerCount->Fill(it,kTriggeredVertexIn);
	}
	if(cand){
	  fh2TriggerCount->Fill(it,kSelectedALICEVertexValid);
	  fh2TriggerVtx->Fill(kSelectedALICEVertexValid*(it+1),zvtx);
	  if(aodEventSelected){
	    fh2TriggerCount->Fill(it,kSelectedALICEVertexIn);
	    fh2TriggerVtx->Fill(kSelectedALICEVertexIn*(it+1),zvtx);
	  }
	}
      }
    }
  }

  if (fDebug > 1)printf(" AliAnalysisTaskJetServices: Analysing event # %5d\n", (Int_t) fEntry);

  
  Double_t ptHard = 0; 
  Double_t nTrials = 1; // Trials for MC trigger 

  fh1Trials->Fill("#sum{ntrials}",fAvgTrials); 
  AliMCEvent* mcEvent = MCEvent();
  //    AliStack *pStack = 0; 
  if(mcEvent){
    AliGenPythiaEventHeader*  pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(mcEvent);
    if(pythiaGenHeader){
      nTrials = pythiaGenHeader->Trials();
      ptHard  = pythiaGenHeader->GetPtHard();
      int iProcessType = pythiaGenHeader->ProcessType();
      // 11 f+f -> f+f
      // 12 f+barf -> f+barf
      // 13 f+barf -> g+g
      // 28 f+g -> f+g
      // 53 g+g -> f+barf
      // 68 g+g -> g+g
      if (fDebug > 10)Printf("%d iProcessType %d",__LINE__, iProcessType);
      if(fDebug>20)AliAnalysisHelperJetTasks::PrintStack(mcEvent);
      fh1PtHard->Fill(ptHard);
      fh1PtHardTrials->Fill(ptHard,nTrials);

    }// if pythia gen header
  }

  // trigger selection
  

  PostData(1, fHistList);
}

Bool_t AliAnalysisTaskJetServices::IsEventSelected(const AliESDEvent* esd){
  if(!esd)return kFALSE;
  const AliESDVertex *vtx = esd->GetPrimaryVertex();
  return IsVertexIn(vtx); // vertex in calls vertex valid
}


Bool_t AliAnalysisTaskJetServices::IsEventSelected(const AliAODEvent* aod) const {
  if(!aod)return kFALSE;
  const AliAODVertex *vtx = aod->GetPrimaryVertex();
  return IsVertexIn(vtx); // VertexIn calls VertexValid
}

Bool_t  AliAnalysisTaskJetServices::IsVertexValid ( const AliESDVertex* vtx) {

  // We can check the type of the vertex by its name and title
  // if vertexer failed title is empty (default c'tor) and ncontributors is 0
  // vtx       name                  title
  // Tracks    PrimaryVertex         VertexerTracksNoConstraint
  // Tracks    PrimaryVertex         VertexerTracksWithConstraint
  // TPC       TPCVertex             VertexerTracksNoConstraint
  // TPC       TPCVertex             VertexerTracksWithConstraint
  // SPD       SPDVertex             vertexer: 3D
  // SPD       SPDVertex             vertexer: Z
  
  Int_t nCont = vtx->GetNContributors();

  if(nCont>=1){
    fEventCutInfoESD |= kContributorsCut1;    
    if(nCont>=2){
      fEventCutInfoESD |= kContributorsCut2;    
      if(nCont>=3){
	fEventCutInfoESD |= kContributorsCut3;    
      }
    }
  }
  
  if(nCont<3)return kFALSE;

  // do not want tpc only primary vertex
  TString vtxName(vtx->GetName());
  if(vtxName.Contains("TPCVertex")){
    fEventCutInfoESD |= kVertexTPC;
    return kFALSE;
  }
  if(vtxName.Contains("SPDVertex"))fEventCutInfoESD |= kVertexSPD;
  if(vtxName.Contains("PrimaryVertex"))fEventCutInfoESD |= kVertexGlobal;


  TString vtxTitle(vtx->GetTitle());
  if(vtxTitle.Contains("vertexer: Z")){
    if(vtx->GetDispersion()>0.02)return kFALSE;   
  }
  fEventCutInfoESD |= kSPDDispersionCut;
  return kTRUE;
}


Bool_t  AliAnalysisTaskJetServices::IsVertexValid ( const AliAODVertex* vtx) const {

  // We can check the type of the vertex by its name and title
  // if vertexer failed title is empty (default c'tor) and ncontributors is 0
  // vtx       name                  title
  // Tracks    PrimaryVertex         VertexerTracksNoConstraint
  // TPC       TPCVertex             VertexerTracksNoConstraint
  // SPD       SPDVertex             vertexer: 3D
  // SPD       SPDVertex             vertexer: Z
  
  if(vtx->GetNContributors()<3)return kFALSE;
  // do not want tpc only primary vertex
  TString vtxName(vtx->GetName());
  if(vtxName.Contains("TPCVertex"))return kFALSE;

  // no dispersion yet...
  /* 
  TString vtxTitle(vtx->GetTitle());
  if(vtxTitle.Contains("vertexer: Z")){
    if(vtx->GetDispersion()>0.02)return kFALSE;
  }
  */
  return kTRUE;
}


Bool_t  AliAnalysisTaskJetServices::IsVertexIn (const AliESDVertex* vtx) {

  if(!IsVertexValid(vtx))return kFALSE;

  Float_t zvtx = vtx->GetZ();
  Float_t yvtx = vtx->GetY();
  Float_t xvtx = vtx->GetX();

  xvtx -= fVtxXMean;
  yvtx -= fVtxYMean;
  zvtx -= fVtxZMean;



  if(TMath::Abs(zvtx)>fVtxZCut){
    return kFALSE;
  }
  fEventCutInfoESD |= kVertexZCut;  
  Float_t r2   = yvtx*yvtx+xvtx*xvtx;      
  if(r2>(fVtxRCut*fVtxRCut)){
    return kFALSE;
  }
  fEventCutInfoESD |= kVertexRCut;  
  return kTRUE;
}


Bool_t  AliAnalysisTaskJetServices::IsVertexIn (const AliAODVertex* vtx) const {

  if(!IsVertexValid(vtx))return kFALSE;

  Float_t zvtx = vtx->GetZ();
  Float_t yvtx = vtx->GetY();
  Float_t xvtx = vtx->GetX();

  xvtx -= fVtxXMean;
  yvtx -= fVtxYMean;
  zvtx -= fVtxZMean;

  Float_t r2   = yvtx*yvtx+xvtx*xvtx;      

  Bool_t vertexIn = TMath::Abs(zvtx)<fVtxZCut&&r2<(fVtxRCut*fVtxRCut);
  return vertexIn;
}

Bool_t AliAnalysisTaskJetServices::IsEventPileUp(const AliESDEvent* esd) const{
  if(!esd)return kFALSE;
  return esd->IsPileupFromSPD();
}

Bool_t AliAnalysisTaskJetServices::IsEventCosmic(const AliESDEvent* esd) const {
  if(!esd)return kFALSE;
  // add track cuts for which we look for cosmics...

  Bool_t isCosmic = kFALSE;
  Int_t nTracks = esd->GetNumberOfTracks();
  Int_t nCosmicCandidates = 0;

  for (Int_t iTrack1 = 0; iTrack1 < nTracks; iTrack1++) {
    AliESDtrack* track1 = (AliESDtrack*)esd->GetTrack(iTrack1);
    if (!track1)  continue;
    UInt_t status1 = track1->GetStatus();
    //If track is ITS stand alone track, skip the track
    if (((status1 & AliESDtrack::kITSin) == 0 || (status1 & AliESDtrack::kTPCin))) continue;
    if(track1->Pt()<fPtMinCosmic) continue;
    //Start 2nd track loop to look for correlations
    for (Int_t iTrack2 = iTrack1+1; iTrack2 < nTracks; iTrack2++) {
      AliESDtrack* track2 = (AliESDtrack*)esd->GetTrack(iTrack2);
      if(!track2) continue;
      UInt_t status2 = track2->GetStatus();
      //If track is ITS stand alone track, skip the track
      if (((status2 & AliESDtrack::kITSin) == 0 || (status2 & AliESDtrack::kTPCin))) continue;
      if(track2->Pt()<fPtMinCosmic) continue;
      //Check if back-to-back
      Double_t mom1[3],mom2[3];
      track1->GetPxPyPz(mom1);
      track2->GetPxPyPz(mom2);
      TVector3 momv1(mom1[0],mom1[1],mom1[2]);
      TVector3 momv2(mom2[0],mom2[1],mom2[2]);
      Float_t theta = (float)(momv1.Phi()-momv2.Phi());
      if(theta<-0.5*TMath::Pi()) theta+=2.*TMath::Pi();

      Float_t deltaPhi = track1->Phi()-track2->Phi();
      if(deltaPhi<-0.5*TMath::Pi()) deltaPhi+=2.*TMath::Pi();

      Float_t rIsol = (float)(TMath::Sqrt( deltaPhi*deltaPhi+(track1->Eta()-track2->Eta())*(track1->Eta()-track2->Eta()) ));
      if(rIsol<fRIsolMinCosmic) continue;

      if(TMath::Abs(TMath::Pi()-theta)<fMaxCosmicAngle) {
	nCosmicCandidates+=1;
	isCosmic = kTRUE;
      }
      
    }
  }

  fh1NCosmicsPerEvent->Fill((float)nCosmicCandidates);

  return isCosmic;
}




void AliAnalysisTaskJetServices::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
}
