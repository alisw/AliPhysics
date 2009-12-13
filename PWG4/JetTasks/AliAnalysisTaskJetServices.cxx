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
#include  "TDatabasePDG.h"

#include "AliAnalysisTaskJetServices.h"
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


#include "AliAnalysisHelperJetTasks.h"

ClassImp(AliAnalysisTaskJetServices)

AliAnalysisTaskJetServices::AliAnalysisTaskJetServices(): AliAnalysisTaskSE(),
  fUseAODInput(kFALSE),
  fAvgTrials(1),
  fZVtxCut(10.),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardTrials(0x0),
  fh2TriggerCount(0x0),
  fh2ESDTriggerCount(0x0),
  fh2TriggerVtx(0x0),
  fh2ESDTriggerVtx(0x0),
  fHistList(0x0)  
{

}

AliAnalysisTaskJetServices::AliAnalysisTaskJetServices(const char* name):
  AliAnalysisTaskSE(name),
  fUseAODInput(kFALSE),
  fAvgTrials(1),
  fZVtxCut(8),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardTrials(0x0),
  fh2TriggerCount(0x0),
  fh2ESDTriggerCount(0x0),
  fh2TriggerVtx(0x0),
  fh2ESDTriggerVtx(0x0),
  fHistList(0x0)  
{
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


  // Connect the AOD

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

  fh2TriggerVtx = new TH2F("fh2TriggerVtx",";Trigger No.;Vtx (cm);Count",AliAnalysisHelperJetTasks::kTrigger,-0.5,AliAnalysisHelperJetTasks::kTrigger-0.5,400,-20.,20.); 
  fHistList->Add(fh2TriggerVtx);

  fh2ESDTriggerVtx = new TH2F("fh2ESDTriggerVtx",";Trigger No.;Vtx (cm);Count",AliAnalysisHelperJetTasks::kTrigger,-0.5,AliAnalysisHelperJetTasks::kTrigger-0.5,400,-20.,20.); 
  fHistList->Add(fh2ESDTriggerVtx);
  

  fh1PtHard = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",nBinPt,binLimitsPt);
  fHistList->Add(fh1PtHard);
  fh1PtHardTrials = new TH1F("fh1PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard}",nBinPt,binLimitsPt);
  fHistList->Add(fh1PtHardTrials);

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


  TH1::AddDirectory(oldStatus);
}

void AliAnalysisTaskJetServices::Init()
{
  //
  // Initialization
  //

  Printf(">>> AnalysisTaskJetServices::Init() debug level %d\n",fDebug);
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
  
  
  // loop over all possible trigger and 
  for(int it = AliAnalysisHelperJetTasks::kAcceptAll;it < AliAnalysisHelperJetTasks::kTrigger;it++){
    Bool_t esdTrig = kFALSE;
    if(esd){
      esdTrig = AliAnalysisHelperJetTasks::IsTriggerFired(esd,(AliAnalysisHelperJetTasks::Trigger)it);
      if(esdTrig)fh2ESDTriggerCount->Fill(it,kAllTriggered);
    }
    Bool_t aodTrig = kFALSE;
    if(aod){
      aodTrig = AliAnalysisHelperJetTasks::IsTriggerFired(aod,(AliAnalysisHelperJetTasks::Trigger)it);
      if(aodTrig)fh2TriggerCount->Fill(it,kAllTriggered);
    }

    // Check wether we have also an SPD vertex
    
    if(aod){
      const AliAODVertex *vtxAOD = aod->GetPrimaryVertex();
      //      Printf(">> AODvtx %s %s",vtxAOD->GetName(),vtxAOD->GetTitle());vtxAOD->Print();
      
      if(vtxAOD->GetNContributors()>0){
	if(aodTrig)fh2TriggerCount->Fill(it,kTriggeredSPDVertex);
	Float_t zvtx = vtxAOD->GetZ();
	Float_t yvtx = vtxAOD->GetY();
	Float_t xvtx = vtxAOD->GetX();
	fh2TriggerVtx->Fill(it,zvtx);
	if(TMath::Abs(zvtx)<fZVtxCut&&aodTrig&&TMath::Abs(xvtx)<0.5&&TMath::Abs(yvtx)<0.5){
	  fh2TriggerCount->Fill(it,kTriggeredVertexIn);
	}
      }
    }
    if(esd){
      const AliESDVertex *vtxESD = esd->GetPrimaryVertex();
      //      Printf(">> ESDvtx %s %s",vtxESD->GetName(),vtxESD->GetTitle());vtxESD->Print();
      if(vtxESD->GetNContributors()>0){
	if(esdTrig)fh2ESDTriggerCount->Fill(it,kTriggeredSPDVertex);
	Float_t zvtx = vtxESD->GetZ();
    	Float_t yvtx = vtxESD->GetY();
    	Float_t xvtx = vtxESD->GetX();
	fh2ESDTriggerVtx->Fill(it,zvtx);
	if(TMath::Abs(zvtx)<fZVtxCut&&esdTrig&&TMath::Abs(xvtx)<0.5&&TMath::Abs(yvtx)<0.5){
	  fh2ESDTriggerCount->Fill(it,kTriggeredVertexIn);
	  // here we select based on ESD info...
	  fh2ESDTriggerCount->Fill(it,kSelected);
	  AliAnalysisHelperJetTasks::Selected(kTRUE,kTRUE);// select this event
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


void AliAnalysisTaskJetServices::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if (fDebug > 1) printf("AnalysisJetServices: Terminate() \n");
}
