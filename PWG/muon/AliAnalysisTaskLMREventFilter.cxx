/**************************************************************************
 * Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
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

//========================================================================
//
//     Contact author: boris.teyssier@cern.ch | antonio.uras@cern.ch
//
//=========================================================================


#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDVZERO.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliVParticle.h"
#include "AliCentrality.h"
#include "AliLog.h"
#include "TBits.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGRPObject.h"
#include "AliTriggerInput.h"
#include "AliTriggerConfiguration.h"
#include "AliInputEventHandler.h"

#include "AliAnalysisMuonUtility.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliLMRMuon.h"
#include "AliLMREvent.h"
#include "AliAnalysisTaskLMREventFilter.h"

ClassImp(AliAnalysisTaskLMREventFilter)

//====================================================================================================================================================

AliAnalysisTaskLMREventFilter::AliAnalysisTaskLMREventFilter() : 
  AliAnalysisTaskSE(), 
  fMuonTrackCuts(0x0),
  fEventTree(0x0),
  fOutputList(0x0),
  fAliLMREvent(0),  
  fhTriggers(0),
  fhBeamType(0),
  fhL0TriggerInputMLL(0),
  fhL0TriggerInputMUL(0),
  fhL0TriggerInputMSL(0),
  fhL0TriggerInputTVX(0),
  fhNMu(0)
{
  //
  // Default constructor
  //
  fNTrigClass=13;
  fTriggerClasses[0]="-B-";
  fTriggerClasses[1]="CMSL7";
  fTriggerClasses[2]="CMSH7";
  fTriggerClasses[3]="CMUL7";
  fTriggerClasses[4]="CMLL7";
  fTriggerClasses[5]="CMSL8";
  fTriggerClasses[6]="CMSH8";
  fTriggerClasses[7]="CMUL8";
  fTriggerClasses[8]="CMLL8";
  fTriggerClasses[9]="C0TVX";
  fTriggerClasses[10]="CINT7-B-NOPF-MUFAST";
  fTriggerClasses[11]="CINT7-B-NOPF-CENTNOTRD";
  fTriggerClasses[12]="CINT7-B-NOPF-CENT";

  fL0TriggerInputMLL = 20; // reference to MLL L0 trigger
  fL0TriggerInputMUL = 21; // reference to MUL L0 trigger
  fL0TriggerInputMSL = 18; // reference to MSL L0 trigger
  fL0TriggerInputTVX =  3; // reference to TVX L0 trigger
  fminContributorsPileUp = 3;
}


//====================================================================================================================================================

AliAnalysisTaskLMREventFilter::AliAnalysisTaskLMREventFilter(const Char_t *name, AliMuonTrackCuts *cuts) : 
  AliAnalysisTaskSE(name), 
  fMuonTrackCuts(cuts),
  fEventTree(0x0),
  fOutputList(0x0),
  fAliLMREvent(0),  
  fhTriggers(0),
  fhBeamType(0),
  fhL0TriggerInputMLL(0),
  fhL0TriggerInputMUL(0),
  fhL0TriggerInputMSL(0),
  fhL0TriggerInputTVX(0),
  fhNMu(0)
{
  // Constructor
  fNTrigClass=13;
  fTriggerClasses[0]="-B-";
  fTriggerClasses[1]="CMSL7";
  fTriggerClasses[2]="CMSH7";
  fTriggerClasses[3]="CMUL7";
  fTriggerClasses[4]="CMLL7";
  fTriggerClasses[5]="CMSL8";
  fTriggerClasses[6]="CMSH8";
  fTriggerClasses[7]="CMUL8";
  fTriggerClasses[8]="CMLL8";
  fTriggerClasses[9]="C0TVX";
  fTriggerClasses[10]="CINT7-B-NOPF-MUFAST";
  fTriggerClasses[11]="CINT7-B-NOPF-CENTNOTRD";
  fTriggerClasses[12]="CINT7-B-NOPF-CENT";

  fL0TriggerInputMLL = 20; // reference to MLL L0 trigger
  fL0TriggerInputMUL = 21; // reference to MUL L0 trigger
  fL0TriggerInputMSL = 18; // reference to MSL L0 trigger
  fL0TriggerInputTVX =  3; // reference to TVX L0 trigger
  fminContributorsPileUp = 3;

  // Define input and output slots here
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());

}

//====================================================================================================================================================

AliAnalysisTaskLMREventFilter::~AliAnalysisTaskLMREventFilter() 
{
  delete fMuonTrackCuts; 
  fMuonTrackCuts=NULL;
  delete fEventTree;
  fEventTree=NULL;
  delete fOutputList;
  fOutputList=NULL;
  delete fAliLMREvent;
  fAliLMREvent=NULL;
  delete fhTriggers;
  delete fhBeamType;
  delete fhL0TriggerInputMLL;
  delete fhL0TriggerInputMUL;
  delete fhL0TriggerInputMSL;
  delete fhL0TriggerInputTVX;
  delete fhNMu;
  fhTriggers=NULL;
  fhNMu=NULL;
}

//====================================================================================================================================================

void AliAnalysisTaskLMREventFilter::UserCreateOutputObjects()
 {
  // Called once
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca); 
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);

  fEventTree = new TTree("Data","Data");
  fAliLMREvent = new AliLMREvent();
  fEventTree->Branch("fAliLMREvent", &fAliLMREvent);
  
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
  
  fhTriggers = new TH1D("hTriggers","L2 Triggers",46,0,fNTrigClass);
  fOutputList->Add(fhTriggers);	  
  fhTriggers->Sumw2();

  fhNMu = new TH2D("hNMu","Number of Muon",20,0,20,46,0,fNTrigClass);
  fOutputList->Add(fhNMu);
  fhNMu->Sumw2();

  fhL0TriggerInputMLL = new TH1D("fhL0TriggerInputMLL","",120,-0.5,119.5);
  fOutputList->Add(fhL0TriggerInputMLL);
  fhL0TriggerInputMLL->Sumw2();

  fhL0TriggerInputMUL = new TH1D("fhL0TriggerInputMUL","",120,-0.5,119.5);
  fOutputList->Add(fhL0TriggerInputMUL);
  fhL0TriggerInputMUL->Sumw2();

  fhL0TriggerInputMSL = new TH1D("fhL0TriggerInputMSL","",120,-0.5,119.5);
  fOutputList->Add(fhL0TriggerInputMSL);
  fhL0TriggerInputMSL->Sumw2();

  fhL0TriggerInputTVX = new TH1D("fhL0TriggerInputTVX","",120,-0.5,119.5);
  fOutputList->Add(fhL0TriggerInputTVX);
  fhL0TriggerInputTVX->Sumw2();

  fhBeamType = new TH1D("fhBeamType","",4,-0.5,3.5);
  fhBeamType->Sumw2();
  fhBeamType->GetXaxis()->SetBinLabel(1,"p-p");
  fhBeamType->GetXaxis()->SetBinLabel(2,"p-A");
  fhBeamType->GetXaxis()->SetBinLabel(3,"A-p");
  fhBeamType->GetXaxis()->SetBinLabel(4,"A-A");
  fOutputList->Add(fhBeamType);

  fhTriggers->GetXaxis()->SetBinLabel(1,fTriggerClasses[0]);
  fhNMu->GetYaxis()->SetBinLabel(1,fTriggerClasses[0]);
  for (Int_t i=1;i<fNTrigClass;i++)
    {
      fhTriggers->GetXaxis()->SetBinLabel(i+1,fTriggerClasses[i]);
      fhNMu->GetYaxis()->SetBinLabel(i+1,fTriggerClasses[i]);
    }

  Int_t CINT7Shift=fNTrigClass;
  Int_t CINT8Shift=CINT7Shift+4;
  for(Int_t i=1;i<5;i++) 
    {
      fhTriggers->GetXaxis()->SetBinLabel(i+CINT7Shift,Form("%s (PS)",fTriggerClasses[i].Data()));
      fhNMu->GetYaxis()->SetBinLabel(i+CINT7Shift,Form("%s (PS)",fTriggerClasses[i].Data()));
      
      fhTriggers->GetXaxis()->SetBinLabel(i+CINT8Shift,Form("%s (PS)",fTriggerClasses[i+4].Data()));
      fhNMu->GetYaxis()->SetBinLabel(i+CINT8Shift,Form("%s (PS)",fTriggerClasses[i+4].Data()));
    }
  
  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+5,Form("%s (PS)",fTriggerClasses[9].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+5,Form("%s (PS)",fTriggerClasses[9].Data()));

  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+6,Form("%s (PS)",fTriggerClasses[10].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+6,Form("%s (PS)",fTriggerClasses[10].Data()));

  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+7,Form("%s (PS)",fTriggerClasses[11].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+7,Form("%s (PS)",fTriggerClasses[11].Data()));

  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+8,Form("%s (PS)",fTriggerClasses[12].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+8,Form("%s (PS)",fTriggerClasses[12].Data()));

  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+9,Form("%s &0MUL (PS)",fTriggerClasses[1].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+9,Form("%s &0MUL (PS)",fTriggerClasses[1].Data()));
  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+10,Form("%s &0MLL (PS)",fTriggerClasses[1].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+10,Form("%s &0MLL (PS)",fTriggerClasses[1].Data()));

  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+11,Form("%s &0MUL (PS)",fTriggerClasses[5].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+11,Form("%s &0MUL (PS)",fTriggerClasses[5].Data()));
  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+12,Form("%s &0MLL (PS)",fTriggerClasses[5].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+12,Form("%s &0MLL (PS)",fTriggerClasses[5].Data()));
  
  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+13,Form("%s &0MUL (PS)",fTriggerClasses[9].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+13,Form("%s &0MUL (PS)",fTriggerClasses[9].Data()));
  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+14,Form("%s &0MSL (PS)",fTriggerClasses[9].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+14,Form("%s &0MSL (PS)",fTriggerClasses[9].Data()));
  
  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+15,Form("%s &0MUL (PS)",fTriggerClasses[10].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+15,Form("%s &0MUL (PS)",fTriggerClasses[10].Data()));
  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+16,Form("%s &0MSL (PS)",fTriggerClasses[10].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+16,Form("%s &0MSL (PS)",fTriggerClasses[10].Data()));

  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+17,Form("%s &0MUL (PS)",fTriggerClasses[11].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+17,Form("%s &0MUL (PS)",fTriggerClasses[11].Data()));
  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+18,Form("%s &0MSL (PS)",fTriggerClasses[11].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+18,Form("%s &0MSL (PS)",fTriggerClasses[11].Data()));
    
  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+19,Form("%s &0MUL (PS)",fTriggerClasses[12].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+19,Form("%s &0MUL (PS)",fTriggerClasses[12].Data()));
  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+20,Form("%s &0MSL (PS)",fTriggerClasses[12].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+20,Form("%s &0MSL (PS)",fTriggerClasses[12].Data()));
  
  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+21,Form("%s &0TVX (PS)",fTriggerClasses[10].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+21,Form("%s &0TVX &0MUL (PS)",fTriggerClasses[10].Data()));
  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+22,Form("%s &0TVX &0MUL (PS)",fTriggerClasses[10].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+22,Form("%s &0TVX &0MUL (PS)",fTriggerClasses[10].Data()));
  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+23,Form("%s &0TVX &0MSL (PS)",fTriggerClasses[10].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+23,Form("%s &0TVX &0MSL (PS)",fTriggerClasses[10].Data()));

  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+24,Form("%s &0TVX (PS)",fTriggerClasses[11].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+24,Form("%s &0TVX &0MUL (PS)",fTriggerClasses[11].Data()));
  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+25,Form("%s &0TVX &0MUL (PS)",fTriggerClasses[11].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+25,Form("%s &0TVX &0MUL (PS)",fTriggerClasses[11].Data()));
  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+26,Form("%s &0TVX &0MSL (PS)",fTriggerClasses[11].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+26,Form("%s &0TVX &0MSL (PS)",fTriggerClasses[11].Data()));
    
  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+27,Form("%s &0TVX (PS)",fTriggerClasses[12].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+27,Form("%s &0TVX &0MUL (PS)",fTriggerClasses[12].Data()));
  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+28,Form("%s &0TVX &0MUL (PS)",fTriggerClasses[12].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+28,Form("%s &0TVX &0MUL (PS)",fTriggerClasses[12].Data()));
  fhTriggers->GetXaxis()->SetBinLabel(CINT8Shift+29,Form("%s &0TVX &0MSL (PS)",fTriggerClasses[12].Data()));
  fhNMu->GetYaxis()->SetBinLabel(CINT8Shift+29,Form("%s &0TVX &0MSL (PS)",fTriggerClasses[12].Data()));
    
  PostData(1, fOutputList);
  PostData(2, fEventTree);
  printf("End of create Output\n");
}


void AliAnalysisTaskLMREventFilter::NotifyRun()
{
  fMuonTrackCuts->SetRun(fInputHandler);
  AliCDBManager *man = AliCDBManager::Instance();
  man->Init();
  man->SetDefaultStorage("raw://"); 
  man->SetRun(fInputHandler->GetEvent()->GetRunNumber()); 

  AliGRPObject* fGRPData = (AliGRPObject*) man->Get("GRP/GRP/Data")->GetObject();
  if(fGRPData->GetBeamType().Contains("p-p"))
    {
      fhBeamType->Fill("p-p",1);
      fminContributorsPileUp = 3;
    }
  else if(fGRPData->GetBeamType().Contains("p-A")) 
    {
      fhBeamType->Fill("p-A",1);
      fminContributorsPileUp = 5;
    }
  else if(fGRPData->GetBeamType().Contains("A-p"))
    {
      fhBeamType->Fill("A-p",1);
      fminContributorsPileUp = 5;
    }
  else if(fGRPData->GetBeamType().Contains("A-A"))
    {
      fhBeamType->Fill("A-A",1);
      fminContributorsPileUp = 5;
    }
  AliTriggerConfiguration *cfg=(AliTriggerConfiguration*)man->Get("GRP/CTP/Config")->GetObject(); 
  TObjArray  inputs = cfg->GetInputs(); 
  for(Int_t i=0;i<inputs.GetEntriesFast();i++)
    {
      AliTriggerInput* inp =  (AliTriggerInput*) inputs[i];
      TString name=inp->GetName();
      if(name.Contains("0MUL"))
	fL0TriggerInputMUL = inp->GetIndexCTP(); // reference to MUL L0 trigger
      if(name.Contains("0MLL"))
	fL0TriggerInputMLL = inp->GetIndexCTP(); // reference to MLL L0 trigger
      if(name.Contains("0MSL"))
	fL0TriggerInputMSL = inp->GetIndexCTP(); // reference to MSL L0 trigger
      if(name.Contains("0TVX"))
	fL0TriggerInputTVX = inp->GetIndexCTP(); // reference to TVX L0 trigger
    }
  fhL0TriggerInputMLL->Fill(fL0TriggerInputMLL);
  fhL0TriggerInputMUL->Fill(fL0TriggerInputMUL);
  fhL0TriggerInputMSL->Fill(fL0TriggerInputMSL);
  fhL0TriggerInputTVX->Fill(fL0TriggerInputTVX);
}
//====================================================================================================================================================

void AliAnalysisTaskLMREventFilter::UserExec(Option_t *)
 {
  //   Main loop
  //   Called for each event
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent *>(InputEvent());  
  if (!fAOD) 
    return;

  UShort_t physicsSelectionMask=1<<0;
  UShort_t L0TriggerInput=1<<0;
  
  if (!IsSelectedTrigger(fAOD, kTRUE,physicsSelectionMask,L0TriggerInput)) 
    return; 
  TString triggerWord(((AliAODHeader*) fAOD->GetHeader())->GetFiredTriggerClasses());
  // ---

  
  AliMultSelection *MultSelection = (AliMultSelection*)fAOD-> FindListObject("MultSelection");

  // All multiplicity are initialized at 166 and is used for error code of multselection non actived
  Float_t Multiplicity_V0M          = 166.;
  Float_t Multiplicity_ADM          = 166.;
  Float_t Multiplicity_SPDTracklets = 166.;
  Float_t Multiplicity_SPDClusters  = 166.;
  Float_t Multiplicity_RefMult05    = 166.;
  Float_t Multiplicity_RefMult08    = 166.;
  Float_t Multiplicity_V0A          = 166.;
  Float_t Multiplicity_V0C          = 166.;
  Float_t Multiplicity_V0EqA        = 166.;
  Float_t Multiplicity_V0EqC        = 166.;
  Float_t Multiplicity_V0EqM        = 166.;
  Float_t Multiplicity_ZNA          = 166.;
  Float_t Multiplicity_ZNC          = 166.;

  if (MultSelection) 
    {
      Multiplicity_V0M          = MultSelection->GetMultiplicityPercentile("V0M");
      Multiplicity_ADM          = MultSelection->GetMultiplicityPercentile("ADM");
      Multiplicity_SPDTracklets = MultSelection->GetMultiplicityPercentile("SPDTracklets");
      Multiplicity_SPDClusters  = MultSelection->GetMultiplicityPercentile("SPDClusters");
      Multiplicity_RefMult05    = MultSelection->GetMultiplicityPercentile("RefMult05");
      Multiplicity_RefMult08    = MultSelection->GetMultiplicityPercentile("RefMult08");
      Multiplicity_RefMult05    = MultSelection->GetMultiplicityPercentile("RefMult05");
      Multiplicity_V0A          = MultSelection->GetMultiplicityPercentile("V0A");
      Multiplicity_V0C          = MultSelection->GetMultiplicityPercentile("V0C");
      Multiplicity_V0EqA        = MultSelection->GetMultiplicityPercentile("V0EqA");
      Multiplicity_V0EqC        = MultSelection->GetMultiplicityPercentile("V0EqC");
      Multiplicity_V0EqM        = MultSelection->GetMultiplicityPercentile("V0EqM");
      Multiplicity_ZNA          = MultSelection->GetMultiplicityPercentile("ZNA");
      Multiplicity_ZNC          = MultSelection->GetMultiplicityPercentile("ZNC");
    }
  

  AliAODVertex *vert = fAOD->GetPrimaryVertex();
  if (!vert) {
    printf ("No vertex found\n");
    return;
  }
  Double_t xvert  = vert->GetX();
  Double_t yvert  = vert->GetY();
  Double_t zvert  = vert->GetZ();
  Int_t vtxcontrib= vert->GetNContributors();
  Double_t evtPlane = fAOD->GetEventplane()->GetEventplane("V0",fAOD,2);
    
  Int_t runNumber = ((AliAODHeader*) fAOD->GetHeader())->GetRunNumber();
  fAliLMREvent->SetRunNumber(runNumber);
  fAliLMREvent->SetEventPlane(evtPlane);
  fAliLMREvent->SetXVertex(xvert);
  fAliLMREvent->SetYVertex(yvert);
  fAliLMREvent->SetZVertex(zvert);
  fAliLMREvent->SetVtxContributors(vtxcontrib);
  fAliLMREvent->SetMultiplicity("V0M",Multiplicity_V0M);
  fAliLMREvent->SetMultiplicity("ADM",Multiplicity_ADM);
  fAliLMREvent->SetMultiplicity("SPDTracklets",Multiplicity_SPDTracklets);
  fAliLMREvent->SetMultiplicity("SPDClusters",Multiplicity_SPDClusters);
  fAliLMREvent->SetMultiplicity("RefMult05",Multiplicity_RefMult05);
  fAliLMREvent->SetMultiplicity("RefMult08",Multiplicity_RefMult08);
  fAliLMREvent->SetMultiplicity("V0A",Multiplicity_V0A);
  fAliLMREvent->SetMultiplicity("V0C",Multiplicity_V0C);
  fAliLMREvent->SetMultiplicity("V0EqA",Multiplicity_V0EqA);
  fAliLMREvent->SetMultiplicity("V0EqC",Multiplicity_V0EqC);
  fAliLMREvent->SetMultiplicity("V0EqM",Multiplicity_V0EqM);
  fAliLMREvent->SetMultiplicity("ZNA",Multiplicity_ZNA);
  fAliLMREvent->SetMultiplicity("ZNC",Multiplicity_ZNC);
  fAliLMREvent->SetTriggerString(triggerWord);
  fAliLMREvent->SetL0TriggerInput(L0TriggerInput);
  fAliLMREvent->SetPhysicsSelectionMask(physicsSelectionMask);
  fAliLMREvent->SetIsPileupFromSPD(fAOD->IsPileupFromSPD(fminContributorsPileUp));
  Int_t nmu=0;
  if(fAOD->GetNumberOfTracks())
     nmu= fAOD->GetNumberOfMuonTracks();
  if (nmu>0)
    {
      Int_t ntotTr = fAOD->GetNumberOfTracks(); 
      AliLMRMuon *trk = NULL;
      for (Int_t itr=0; itr<ntotTr; itr++) 
	{ 
	  AliAODTrack *track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(itr));
	  if(!track->IsMuonTrack()) continue;
    
	  Double_t p[3];
    
	  Short_t charge = track->Charge(); 
	  track->GetPxPyPz(p);
	  Double_t chi2Match = track->GetChi2MatchTrigger();
	  Int_t match = track->GetMatchTrigger();	
	  Double_t chi2 = track->Chi2perNDF();
	  Double_t rAbs = track->GetRAtAbsorberEnd();
    
	  TVector3 dcaAtVz  = fMuonTrackCuts->GetCorrectedDCA(track);
	  Double_t pTotMean = fMuonTrackCuts->GetAverageMomentum(track);
	  Double_t pDca = pTotMean * dcaAtVz.Mag();
	  // Create new Muon
	  trk=fAliLMREvent->AddMuon();
	  trk->SetMomentum(p[0],p[1],p[2]);
	  trk->SetCharge(charge);
	  trk->SetChi2Match(chi2Match);
	  trk->SetChi2(chi2);
	  trk->SetRabs(rAbs);
	  trk->SetpDCA(pDca);
	  trk->SetTriggerMatch(match);
	  trk->SetSelectionMask(fMuonTrackCuts->GetSelectionMask(track));
	  trk->SetLocalBoard((UShort_t)AliAnalysisMuonUtility::GetLoCircuit(track));
	  }
    }
  
 
 fEventTree->Fill();
 fAliLMREvent->Clear("");

 PostData(1, fOutputList); 
 PostData(2, fEventTree);
 // fAOD=NULL;

}

//====================================================================================================================================================

void AliAnalysisTaskLMREventFilter::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  
  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }  
}

//====================================================================================================================================================

Bool_t AliAnalysisTaskLMREventFilter::IsSelectedTrigger(AliAODEvent *fAOD, Bool_t fillHisto,UShort_t &physicsSelectionMask,UShort_t &L0TriggerInput)
{
  Bool_t evtToBeProcessed = kFALSE;
  Bool_t goodTrig = kFALSE;
  Int_t nmu=0;
   if (fAOD->GetNumberOfTracks())
     nmu= fAOD->GetNumberOfMuonTracks();
  // -- Check Trigger 
  TString trigStr(((AliAODHeader*) fAOD->GetHeader())->GetFiredTriggerClasses());
  TObjArray * tokens = trigStr.Tokenize(" ");
  Int_t ntokens = tokens->GetEntriesFast();
  for (Int_t itoken = 0; itoken < ntokens; ++itoken)
    {
      if ((((TObjString*)tokens->At(itoken))->String()).Contains(fTriggerClasses[0])) //for -B-
	{
	  goodTrig = kTRUE;
	  fhNMu->Fill(nmu,fTriggerClasses[0],1);
	  break;
	}
    }
  if(!goodTrig)
    {
      delete tokens;
      return evtToBeProcessed;
    }

  if (fillHisto) 
    fhTriggers->Fill(fTriggerClasses[0].Data(),1);

  
  for(Int_t i=1;i<fNTrigClass;i++)
    {
      goodTrig = kFALSE;
      for (Int_t itoken = 0; itoken < ntokens; ++itoken)
	{
	  if ((((TObjString*)tokens->At(itoken))->String()).Contains(fTriggerClasses[i])) 
	    {
	      goodTrig = kTRUE;
	      fhNMu->Fill(nmu,fTriggerClasses[i],1);
	      evtToBeProcessed = kTRUE;
	      break;
	    }
	}
      if(goodTrig)
	if (fillHisto) 
	  fhTriggers->Fill(fTriggerClasses[i].Data(),1);
    }

  UInt_t fSelectedMaskMuonINT7=AliVEvent::kMuonSingleLowPt7|AliVEvent::kMuonSingleHighPt7|AliVEvent::kMuonLikeLowPt7|AliVEvent::kMuonUnlikeLowPt7;
  UInt_t fSelectedMaskMuonINT8=AliVEvent::kMuonSingleLowPt8|AliVEvent::kMuonSingleHighPt8|AliVEvent::kMuonLikeLowPt8|AliVEvent::kMuonUnlikeLowPt8;

  UInt_t fSelectMask = fInputHandler->IsEventSelected(); 
  Bool_t isMuonINT7selected  = fSelectMask&fSelectedMaskMuonINT7; 
  Bool_t isMuonC0TVXselected = fSelectMask&fSelectedMaskMuonINT8; 

  if(fSelectMask&(AliVEvent::kINT7|AliVEvent::kINT7inMUON))
    physicsSelectionMask|=1<<1;
  if(fSelectMask&AliVEvent::kINT8)
    physicsSelectionMask|=1<<2;
  if(isMuonINT7selected)
    physicsSelectionMask|=1<<3;
  if(isMuonC0TVXselected)
    physicsSelectionMask|=1<<4;

  UInt_t inpmask = fAOD->GetHeader()->GetL0TriggerInputs();
  Int_t is0TVXfired = (inpmask & (1<<(fL0TriggerInputTVX-1)));
  Int_t is0MSLfired = (inpmask & (1<<(fL0TriggerInputMSL-1)));
  Int_t is0MULfired = (inpmask & (1<<(fL0TriggerInputMUL-1)));
  Int_t is0MLLfired = (inpmask & (1<<(fL0TriggerInputMLL-1)));

  if(is0TVXfired)
    L0TriggerInput|=1<<1;
  if(is0MSLfired)
    L0TriggerInput|=1<<2;
  if(is0MULfired)
    L0TriggerInput|=1<<3;
  if(is0MLLfired)
    L0TriggerInput|=1<<4;

  if(fSelectMask&(AliVEvent::kINT7|AliVEvent::kINT7inMUON))
    {
      for(Int_t i=0;i<3;i++)
	{
	  if(trigStr.Contains(fTriggerClasses[10+i].Data()))
	    {
	      fhTriggers->Fill(Form("%s (PS)",fTriggerClasses[10+i].Data()),1);
	      fhNMu ->Fill(nmu,Form("%s (PS)",fTriggerClasses[10+i].Data()),1);
	      if(is0TVXfired)
		{
		  fhTriggers->Fill(Form("%s &0TVX (PS)",fTriggerClasses[10+i].Data()),1);
		  fhNMu ->Fill(nmu,Form("%s &0TVX (PS)",fTriggerClasses[10+i].Data()),1);
		}
	      if(is0MSLfired)
		{	 
		  fhTriggers->Fill(Form("%s &0MSL (PS)",fTriggerClasses[10+i].Data()),1);
		  fhNMu ->Fill(nmu,Form("%s &0MSL (PS)",fTriggerClasses[10+i].Data()),1);
		  if(is0TVXfired)
		    {
		      fhTriggers->Fill(Form("%s &0TVX &0MSL (PS)",fTriggerClasses[10+i].Data()),1);
		      fhNMu ->Fill(nmu,Form("%s &0TVX &0MSL (PS)",fTriggerClasses[10+i].Data()),1);
		    }
		}
	      if(is0MULfired)
		{	 
		  fhTriggers->Fill(Form("%s &0MUL (PS)",fTriggerClasses[10+i].Data()),1);
		  fhNMu ->Fill(nmu,Form("%s &0MUL (PS)",fTriggerClasses[10+i].Data()),1);
		  if(is0TVXfired)
		    {
		      fhTriggers->Fill(Form("%s &0TVX &0MUL (PS)",fTriggerClasses[10+i].Data()),1);
		      fhNMu ->Fill(nmu,Form("%s &0TVX &0MUL (PS)",fTriggerClasses[10+i].Data()),1);
		    }
		}
	    }
	}
    }
  if(fSelectMask&AliVEvent::kINT8)
    {
      if(trigStr.Contains(fTriggerClasses[9].Data()))
	{
	  fhTriggers->Fill(Form("%s (PS)",fTriggerClasses[9].Data()),1);
	  fhNMu ->Fill(nmu,Form("%s (PS)",fTriggerClasses[9].Data()),1);
	  if(is0MSLfired)
	    {
	      fhTriggers->Fill(Form("%s &0MSL (PS)",fTriggerClasses[9].Data()),1);
	      fhNMu ->Fill(nmu,Form("%s &0MSL (PS)",fTriggerClasses[9].Data()),1);
	    }
	  if(is0MULfired)
	    {
	      fhTriggers->Fill(Form("%s &0MUL (PS)",fTriggerClasses[9].Data()),1);
	      fhNMu ->Fill(nmu,Form("%s &0MUL (PS)",fTriggerClasses[9].Data()),1);
	    }
	}
    }
  if(isMuonINT7selected)
    {
      for(Int_t i=1;i<5;i++) 
	{
	  if(trigStr.Contains(fTriggerClasses[i].Data()))
	    {
	      fhTriggers->Fill(Form("%s (PS)",fTriggerClasses[i].Data()),1);
	      fhNMu ->Fill(nmu,Form("%s (PS)",fTriggerClasses[i].Data()),1);
	    }
	 }
      if(trigStr.Contains("CMSL7"))
	{
	  if(is0MULfired)
	    {
	      fhTriggers->Fill(Form("%s &0MUL (PS)",fTriggerClasses[1].Data()),1);
	      fhNMu ->Fill(nmu,Form("%s &0MUL (PS)",fTriggerClasses[1].Data()),1);
	    }
	  if(is0MLLfired)
	    {
	      fhTriggers->Fill(Form("%s &0MLL (PS)",fTriggerClasses[1].Data()),1);
	      fhNMu ->Fill(nmu,Form("%s &0MLL (PS)",fTriggerClasses[1].Data()),1);
	    }
	}
    }
  if(isMuonC0TVXselected)
    {
      for(Int_t i=5;i<9;i++) 
	{
	  if(trigStr.Contains(fTriggerClasses[i].Data()))
	    {
	      fhTriggers->Fill(Form("%s (PS)",fTriggerClasses[i].Data()),1);
	      fhNMu ->Fill(nmu,Form("%s (PS)",fTriggerClasses[i].Data()),1);
	    }
	}
      if(trigStr.Contains("CMSL8"))
	{
	  if(is0MULfired)
	    {
	      fhTriggers->Fill(Form("%s &0MUL (PS)",fTriggerClasses[5].Data()),1);
	      fhNMu ->Fill(nmu,Form("%s &0MUL (PS)",fTriggerClasses[5].Data()),1);
	    }
	  if(is0MLLfired)
	    {
	      fhTriggers->Fill(Form("%s &0MLL (PS)",fTriggerClasses[5].Data()),1);
	      fhNMu ->Fill(nmu,Form("%s &0MLL (PS)",fTriggerClasses[5].Data()),1);
	    }
	}
    }
  delete tokens;

  return evtToBeProcessed;
}
