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
#include <TRandom3.h>
#include <TClonesArray.h>
#include "TDatabasePDG.h"

#include "AliAnalysisTaskJetServices.h"
#include "AliCentrality.h"
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

AliAODHeader*  AliAnalysisTaskJetServices::fgAODHeader = NULL;
TClonesArray*   AliAnalysisTaskJetServices::fgAODVertices = NULL;

AliAnalysisTaskJetServices::AliAnalysisTaskJetServices(): 
  AliAnalysisTaskSE(),
  fUseAODInput(kFALSE),
  fUsePhysicsSelection(kFALSE),
  fMC(kFALSE),
  fFilterAODCollisions(kFALSE),
  fPhysicsSelectionFlag(AliVEvent::kMB),
  fSelectionInfoESD(0),
  fEventCutInfoESD(0),
  fFilterMask(0),
  fRPSubeventMethod(0),
  fAvgTrials(1),
  fVtxXMean(0),
  fVtxYMean(0),
  fVtxZMean(0),
  fVtxRCut(1.),
  fVtxZCut(8.),
  fPtMinCosmic(5.),
  fRIsolMinCosmic(3.),
  fMaxCosmicAngle(0.01),
  fCentrality(101),
  fTrackRecEtaWindow(0.9),
  fMinTrackPt(0.15),
  fRPAngle(0),
  fRandomizer(0),
  fNonStdFile(""),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardTrials(0x0),
  fh1SelectionInfoESD(0x0),
  fh1EventCutInfoESD(0),
  fh1CentralityESD(0),
  fh1Centrality(0),
  fh1RP(0),
  fh2TriggerCount(0x0),
  fh2ESDTriggerCount(0x0),
  fh2TriggerVtx(0x0),
  fh2ESDTriggerVtx(0x0),
  fh2ESDTriggerRun(0x0),
  fh2VtxXY(0x0),
  fh1NCosmicsPerEvent(0x0),
  fh2RPSubevents(0x0),
  fh2RPCentrality(0x0),
  fh2RPDeltaRP(0x0),
  fh2RPQxQy(0x0),
  fh2RPCosDeltaRP(0x0),
  fh3PhiWeights(0x0),
  fh3RPPhiTracks(0x0),
  fTriggerAnalysis(0x0),
  fHistList(0x0)  
{
  fRunRange[0] = fRunRange[1] = 0; 
  fFlatA[0] =   fFlatA[1] = 0;
  fFlatB[0] =   fFlatB[1] = 0;
  fDeltaQxy[0] =   fDeltaQxy[1] = 0; 

}

AliAnalysisTaskJetServices::AliAnalysisTaskJetServices(const char* name):
  AliAnalysisTaskSE(name),
  fUseAODInput(kFALSE),
  fUsePhysicsSelection(kFALSE),
  fMC(kFALSE),
  fFilterAODCollisions(kFALSE),
  fPhysicsSelectionFlag(AliVEvent::kMB),
  fSelectionInfoESD(0),
  fEventCutInfoESD(0),
  fFilterMask(0),
  fRPSubeventMethod(0),
  fAvgTrials(1),
  fVtxXMean(0),
  fVtxYMean(0),
  fVtxZMean(0),
  fVtxRCut(1.),
  fVtxZCut(8.),
  fPtMinCosmic(5.),
  fRIsolMinCosmic(3.),
  fMaxCosmicAngle(0.01),
  fCentrality(101),
  fTrackRecEtaWindow(0.9),
  fMinTrackPt(0.15),
  fRPAngle(0),
  fRandomizer(0),
  fNonStdFile(""),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardTrials(0x0),
  fh1SelectionInfoESD(0x0),
  fh1EventCutInfoESD(0),
  fh1CentralityESD(0),
  fh1Centrality(0),
  fh1RP(0),
  fh2TriggerCount(0x0),
  fh2ESDTriggerCount(0x0),
  fh2TriggerVtx(0x0),
  fh2ESDTriggerVtx(0x0),
  fh2ESDTriggerRun(0x0),
  fh2VtxXY(0x0),
  fh1NCosmicsPerEvent(0x0),
  fh2RPSubevents(0x0),
  fh2RPCentrality(0x0),
  fh2RPDeltaRP(0x0),
  fh2RPQxQy(0x0),
  fh2RPCosDeltaRP(0x0),
  fh3PhiWeights(0x0),
  fh3RPPhiTracks(0x0),

  fTriggerAnalysis(0x0),
  fHistList(0x0)  
{
  fRunRange[0] = fRunRange[1] = 0; 
  fFlatA[0] =   fFlatA[1] = 0;
  fFlatB[0] =   fFlatB[1] = 0;
  fDeltaQxy[0] =   fDeltaQxy[1] = 0; 
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

  fRandomizer = new TRandom3(0);
  if (fDebug > 1) printf("AnalysisTaskJetServices::UserCreateOutputObjects() \n");

  OpenFile(1);
  if(!fHistList)fHistList = new TList();
  fHistList->SetOwner();
  PostData(1, fHistList); // post data in any case once

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
  

  fh1CentralityESD = new TH1F("fh1CentralityESD","cent",102,-0.5,101.5);
  fHistList->Add(fh1CentralityESD);
  
  fh1Centrality = new TH1F("fh1Centrality","cent",102,-0.5,101.5);
  fHistList->Add(fh1Centrality);

  fh1RP = new TH1F("fh1RP","RP;#Psi",440, -1.*TMath::Pi(), 2.*TMath::Pi());
  fHistList->Add(fh1RP);

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

  fh2ESDTriggerRun = new TH2F("fh2ESDTriggerRun","Eventclass vs run number:run;trigger",(Int_t)(1+fRunRange[1]-fRunRange[0]),fRunRange[0]-0.5,fRunRange[1]+0.5,10,-0.5,9.5);
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

  
  fh2RPSubevents = new TH2F("fh2RPSubevents" ,"Reaction Plane Angle" , 180, 0, TMath::Pi(), 180, 0, TMath::Pi());
  fHistList->Add( fh2RPSubevents);

  fh2RPCentrality = new TH2F("fh2RPCentrality" ,"Reaction Plane Angle" , 20, 0.,100., 180, 0, TMath::Pi());
  fHistList->Add(fh2RPCentrality);

  fh2RPDeltaRP   = new TH2F("fh2DeltaRP" ,"Delta Reaction Plane Angle" , 100, -TMath::Pi()/2, TMath::Pi()/2,20,0.,100.0);
  fHistList->Add(fh2RPDeltaRP);

  fh2RPQxQy      = new TH2F("fh2RPQxQy" ,"" , 100, -100,100,100,-100,100);
  fHistList->Add(fh2RPQxQy);

  fh2RPCosDeltaRP = new TH2F("fh2RPCosDeltaRP" ,"" , 20, 0.001,100.001,100,-1,1);
  fHistList->Add(fh2RPCosDeltaRP);

  fh3RPPhiTracks = new TH3F("fh3RPPhiTracks","Phi Tracks Pt Centrality", 10, 0.,100.,20,-5,5,180, 0, 2*TMath::Pi());
  fHistList->Add(fh3RPPhiTracks);
  

  fHistList->Add(fh1NCosmicsPerEvent),


  TH1::AddDirectory(oldStatus);

  // Add an AOD branch for replication
  if(fNonStdFile.Length()){
     if (fDebug > 1) AliInfo("Replicating header");
     fgAODHeader = new AliAODHeader;
     AddAODBranch("AliAODHeader",&fgAODHeader,fNonStdFile.Data());
     if (fDebug > 1) AliInfo("Replicating primary vertices");
     fgAODVertices = new TClonesArray("AliAODVertex",3);
     fgAODVertices->SetName("vertices");
     AddAODBranch("TClonesArray",&fgAODVertices,fNonStdFile.Data());
  }
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
  AliAnalysisHelperJetTasks::EventClass(kTRUE,0);
  AliAnalysisHelperJetTasks::ReactionPlane(kTRUE,0); // set slection to false
  fSelectionInfoESD = 0; // reset
  fEventCutInfoESD = 0; // reset
  AliAnalysisHelperJetTasks::SelectInfo(kTRUE,fSelectionInfoESD); // set slection to false


  static AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
   


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
  if(aod&&fDebug>2){
    aod->Print();
    Printf("Vertices %d",aod->GetNumberOfVertices());
    Printf("tracks %d",aod->GetNumberOfTracks());
    Printf("jets %d",aod->GetNJets());
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
  }
 
  // Apply additional constraints
  Bool_t esdEventSelected = IsEventSelected(esd);
  Bool_t esdEventPileUp = IsEventPileUp(esd);
  Bool_t esdEventCosmic = IsEventCosmic(esd);

  Bool_t aodEventSelected = IsEventSelected(aod);

  Bool_t physicsSelection = ((fInputHandler->IsEventSelected())&fPhysicsSelectionFlag);


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

  Float_t cent = 100;
  if(aod)cent = aod->GetHeader()->GetCentrality();
  if(cent<0)cent = 101;
  fCentrality = cent;
  fRPAngle = 0;

  if(esd){
    const AliESDVertex *vtxESD = esd->GetPrimaryVertex();
    esdVtxValid = IsVertexValid(vtxESD);
    esdVtxIn = IsVertexIn(vtxESD);
    if(aodH&&physicsSelection&&fFilterAODCollisions&&aod){
      if(fDebug)Printf("%s:%d Centrality %3.3f vtxin %d",(char*)__FILE__,__LINE__,cent,esdVtxIn);
      if(cent<=80&&esdVtxIn){
	aodH->SetFillAOD(kTRUE);
	aodH->SetFillExtension(kTRUE);
      }
    }


    Float_t zvtx = vtxESD->GetZ();
    Int_t  iCl = GetEventClass(esd);
    AliAnalysisHelperJetTasks::EventClass(kTRUE,iCl);
    Bool_t cand = physicsSelection;

    if(fDebug)Printf("%s:%d %d %d %d Icl %d",(char*)__FILE__,__LINE__,esdVtxValid,esdVtxIn,cand,iCl);
    fh2ESDTriggerCount->Fill(0.,kAllTriggered); 
    fh2ESDTriggerCount->Fill(iCl,kAllTriggered); 
    if(cand){
      fh2ESDTriggerCount->Fill(0.,kSelectedALICE); 
      fh2ESDTriggerCount->Fill(iCl,kSelectedALICE); 
      fh2ESDTriggerVtx->Fill(kSelectedALICE*(iCl+1),zvtx);
    }
    //    if(!fUsePhysicsSelection)cand =  AliAnalysisHelperJetTasks::IsTriggerFired(esd,AliAnalysisHelperJetTasks::kMB1);
    if(esdVtxValid){
      fh2ESDTriggerCount->Fill(0.,kTriggeredVertex);
      fh2ESDTriggerCount->Fill(iCl,kTriggeredVertex);
      fh2ESDTriggerVtx->Fill(iCl,zvtx);
      if(esdVtxIn){
	fh2ESDTriggerCount->Fill(0.,kTriggeredVertexIn);
	fh2ESDTriggerCount->Fill(iCl,kTriggeredVertexIn);
	fh2ESDTriggerVtx->Fill(kTriggeredVertexIn*(iCl+1),zvtx);
      }
      if(cand){
	fh2ESDTriggerCount->Fill(0.,kSelectedALICEVertexValid);
	fh2ESDTriggerCount->Fill(iCl,kSelectedALICEVertexValid);
	fh2ESDTriggerVtx->Fill(kSelectedALICEVertexValid*(iCl+1),zvtx);
      }
    }

    if(cand&&esdVtxIn&&iCl<5){
      fh2ESDTriggerCount->Fill(0.,kSelectedALICEVertexIn);
      fh2ESDTriggerCount->Fill(iCl,kSelectedALICEVertexIn);
      fh2ESDTriggerVtx->Fill(kSelectedALICEVertexIn*(iCl+1),zvtx);
      fh2ESDTriggerVtx->Fill(kSelected*(iCl+1),zvtx);
      fh2ESDTriggerCount->Fill(iCl,kSelected);
      fh2ESDTriggerCount->Fill(0.,kSelected);
      AliAnalysisHelperJetTasks::Selected(kTRUE,kTRUE);// select this event
      if(esd->GetCentrality()){
	Float_t tmpCent = 100;
	tmpCent = esd->GetCentrality()->GetCentralityPercentile("V0M");
	if(tmpCent<0)tmpCent = 101;
	fh1CentralityESD->Fill(tmpCent);
      }
    }
  }



  if(aod){
    const AliAODVertex *vtxAOD = aod->GetPrimaryVertex();
    aodVtxValid = IsVertexValid(vtxAOD);
    aodVtxIn = IsVertexIn(vtxAOD);
    Float_t zvtx = vtxAOD->GetZ();
    Int_t  iCl = GetEventClass(aod);
    AliAnalysisHelperJetTasks::EventClass(kTRUE,iCl);
    Bool_t cand = aod->GetHeader()->GetOfflineTrigger()&fPhysicsSelectionFlag;
    if(fDebug)Printf("%s:%d AOD selection %d %d",(char*)__FILE__,__LINE__,cand,aod->GetHeader()->GetOfflineTrigger());
    fh2TriggerCount->Fill(0.,kAllTriggered); 
    fh2TriggerCount->Fill(iCl,kAllTriggered); 
    if(cand){
      fh2TriggerCount->Fill(0.,kSelectedALICE); 
      fh2TriggerCount->Fill(iCl,kSelectedALICE); 
      fh2TriggerVtx->Fill(kSelectedALICE*(iCl+1),zvtx);
    }
    if(aodVtxValid){
      fh2TriggerCount->Fill(0.,kTriggeredVertex);
      fh2TriggerCount->Fill(iCl,kTriggeredVertex);
      fh2TriggerVtx->Fill(iCl,zvtx);
      if(aodVtxIn){
	fh2TriggerCount->Fill(0.,kTriggeredVertexIn);
	fh2TriggerCount->Fill(iCl,kTriggeredVertexIn);
	fh2TriggerVtx->Fill(kTriggeredVertexIn*(iCl+1),zvtx);
      }
      if(cand){
	fh2TriggerCount->Fill(0.,kSelectedALICEVertexValid);
	fh2TriggerCount->Fill(iCl,kSelectedALICEVertexValid);
	fh2TriggerVtx->Fill(kSelectedALICEVertexValid*(iCl+1),zvtx);
      }
    }
    if(cand&&aodVtxIn&&iCl<5){
      fh2TriggerCount->Fill(0.,kSelectedALICEVertexIn);
      fh2TriggerCount->Fill(iCl,kSelectedALICEVertexIn);
      fh2TriggerVtx->Fill(kSelectedALICEVertexIn*(iCl+1),zvtx);
      fh2TriggerVtx->Fill(kSelected*(iCl+1),zvtx);
      fh2TriggerCount->Fill(iCl,kSelected);
      fh2TriggerCount->Fill(0.,kSelected);
      fh1Centrality->Fill(cent);
      AliAnalysisHelperJetTasks::Selected(kTRUE,kTRUE);// select this event
      TList recTracks;
      GetListOfTracks(&recTracks);
      CalculateReactionPlaneAngle(&recTracks);
      AliAnalysisHelperJetTasks::ReactionPlane(kTRUE,fRPAngle); // set slection to false
      if(fUseAODInput&&cent<=80){
	if(fFilterAODCollisions&&aod){
	  aodH->SetFillAOD(kTRUE);
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
  
  // replication of 
  if(fNonStdFile.Length()&&aod){
    if (fgAODHeader){
      *fgAODHeader =  *(dynamic_cast<AliAODHeader*>(aod->GetHeader()));
      Double_t q[2] = {fRPAngle,fRPAngle};
      fgAODHeader->SetQTheta(q,2);
    }
    if(fgAODVertices){
      fgAODVertices->Delete();
      TClonesArray &vertices = *fgAODVertices;
      const AliAODVertex *vtxAOD = aod->GetPrimaryVertex();
      // we only use some basic information, 
      

      Double_t pos[3];
      Double_t covVtx[6];

      vtxAOD->GetXYZ(pos); // position                                                                
      vtxAOD->GetCovMatrix(covVtx); //covariance matrix                                                      
      Int_t jVertices = 0;
      AliAODVertex * vtx = new(vertices[jVertices++])
        AliAODVertex(pos, covVtx, vtxAOD->GetChi2perNDF(), NULL, -1, AliAODVertex::kPrimary);
      vtx->SetName(vtxAOD->GetName());
      vtx->SetTitle(vtxAOD->GetTitle());
      TString vtitle = vtxAOD->GetTitle();
      vtx->SetNContributors(vtxAOD->GetNContributors());

      // Add SPD "main" vertex                                                                    
      const AliAODVertex *vtxS = aod->GetPrimaryVertexSPD();
      vtxS->GetXYZ(pos); // position
      vtxS->GetCovMatrix(covVtx); //covariance matrix                                                      
      AliAODVertex * mVSPD = new(vertices[jVertices++])
	AliAODVertex(pos, covVtx, vtxS->GetChi2perNDF(), NULL, -1, AliAODVertex::kMainSPD);
      mVSPD->SetName(vtxS->GetName());
      mVSPD->SetTitle(vtxS->GetTitle());
      mVSPD->SetNContributors(vtxS->GetNContributors());
      
      // Add tpc only vertex
      if(esd){
	const AliESDVertex *vtxT =  esd->GetPrimaryVertexTPC();
	vtxT->GetXYZ(pos); // position                                                    
	vtxT->GetCovMatrix(covVtx); //covariance matrix                                               
	AliAODVertex * mVTPC = new(vertices[jVertices++])
	  AliAODVertex(pos, covVtx, vtxT->GetChi2toNDF(), NULL, -1, AliAODVertex::kMainTPC);
	mVTPC->SetName(vtxT->GetName());
	mVTPC->SetTitle(vtxT->GetTitle());
	mVTPC->SetNContributors(vtxT->GetNContributors());
      }
    }
  }
  
  PostData(1, fHistList);
}

Bool_t AliAnalysisTaskJetServices::IsEventSelected(const AliESDEvent* esd){
  if(!esd)return kFALSE;
  const AliESDVertex *vtx = esd->GetPrimaryVertex();
  return IsVertexIn(vtx); // vertex in calls vertex valid
}

AliAnalysisTaskJetServices::~AliAnalysisTaskJetServices(){
  if(fgAODVertices){
    fgAODVertices->Delete();
    delete fgAODVertices;
  }
  delete fRandomizer;
  if(fgAODHeader)delete fgAODHeader;
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

  if(fDebug){
    Printf(" n contrib %d",vtx->GetNContributors());
    vtx->Print();
  }
  
  //  if(vtx->GetNContributors()<3)return kFALSE;
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


Int_t AliAnalysisTaskJetServices::GetEventClass(AliESDEvent *esd){

  Float_t cent = 999;
  if(esd->GetCentrality()){
    cent = esd->GetCentrality()->GetCentralityPercentile("V0M");
  }
  if(cent>80||cent<0)return 5;
  if(cent>50)return 4;
  if(cent>30)return 3;
  if(cent>10)return 2;
  return 1;

}


Int_t AliAnalysisTaskJetServices::GetEventClass(AliAODEvent *aod){

  Float_t cent = aod->GetHeader()->GetCentrality();
  if(cent>80||cent<0)return 5;
  if(cent>50)return 4;
  if(cent>30)return 3;
  if(cent>10)return 2;
  return 1;

}


void AliAnalysisTaskJetServices::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
}

Bool_t AliAnalysisTaskJetServices::CalculateReactionPlaneAngle(const TList *trackList)
{

  if(!trackList)return kFALSE;
  fRPAngle=0;

  // need to get this info from elsewhere??

  Double_t fPsiRP =0,fDeltaPsiRP = 0;
   
   
    
  TVector2 mQ,mQ1,mQ2;
  Float_t mQx= fDeltaQxy[0], mQy=fDeltaQxy[1];
  
  Float_t mQx1=fDeltaQxy[0], mQy1=fDeltaQxy[1];
  Float_t mQx2=fDeltaQxy[0], mQy2=fDeltaQxy[1];
  
  AliVParticle *track=0x0;
  Int_t count[3]={0,0,0};
  

  for (Int_t iter=0;iter<trackList->GetEntries();iter++){

    track=(AliVParticle*)trackList->At(iter);
    
    //cuts already applied before
    // Comment DCA not correctly implemented yet for AOD tracks
    
    Double_t momentum;
    if(track->Charge()>0){momentum=track->Pt();}
    else{momentum=-track->Pt();}

       

    // For Weighting
    fh3RPPhiTracks->Fill(fCentrality,momentum,track->Phi());
    count[0]++;

    Double_t phiweight=GetPhiWeight(track->Phi(),momentum);
    //    Double_t phiweight=1; 
    Double_t weight=2;
    if(track->Pt()<2){weight=track->Pt();}
    

    mQx += (cos(2*track->Phi()))*weight*phiweight;
    mQy += (sin(2*track->Phi()))*weight*phiweight;

    // Make random Subevents

    if(fRPSubeventMethod==0){
      if(fRandomizer->Binomial(1,0.5)){
	mQx1 += (cos(2*track->Phi()))*weight*phiweight;
	mQy1 += (sin(2*track->Phi()))*weight*phiweight;
	count[1]++;}
      else{
	mQx2 += (cos(2*track->Phi()))*weight*phiweight;
	mQy2 += (sin(2*track->Phi()))*weight*phiweight;
	count[2]++;}
    }
    else if(fRPSubeventMethod==1){
      // Make eta dependent subevents
      if(track->Eta()>0){
	mQx1 += (cos(2*track->Phi()))*weight*phiweight;
	mQy1 += (sin(2*track->Phi()))*weight*phiweight;
	count[1]++;}
      else{
	mQx2 += (cos(2*track->Phi()))*weight*phiweight;
	mQy2 += (sin(2*track->Phi()))*weight*phiweight;
	count[2]++;}
    }

  }



  //If no track passes the cuts, the ,Q.Phi() will return Pi and a peak at Pi/2 in the RP Angular Distribution will appear
  if(count[0]==0||count[1]==0||count[2]==0){
    return kFALSE;
  }

  mQ.Set(mQx,mQy);
  mQ1.Set(mQx1,mQy1);
  mQ2.Set(mQx2,mQy2);

  // cout<<"MQ"<<mQx<<" " <<mQy<<" psi"<<endl;

  fPsiRP=mQ.Phi()/2;
    
  //Correction
  fPsiRP+=fFlatA[0]*TMath::Cos(2*fPsiRP)+fFlatB[0]*TMath::Sin(2*fPsiRP)+fFlatA[1]*TMath::Cos(4*fPsiRP)+fFlatB[1]*TMath::Sin(4*fPsiRP);

  Double_t fPsiRP1=mQ1.Phi()/2;
  fPsiRP1+=fFlatA[0]*TMath::Cos(2*fPsiRP1)+fFlatB[0]*TMath::Sin(2*fPsiRP1)+fFlatA[1]*TMath::Cos(4*fPsiRP1)+fFlatB[1]*TMath::Sin(4*fPsiRP1);
  Double_t fPsiRP2=mQ2.Phi()/2;
  fPsiRP2+=fFlatA[0]*TMath::Cos(2*fPsiRP2)+fFlatB[0]*TMath::Sin(2*fPsiRP2)+fFlatA[1]*TMath::Cos(4*fPsiRP2)+fFlatB[1]*TMath::Sin(4*fPsiRP2);
  fDeltaPsiRP=fPsiRP1-fPsiRP2;
  
  if(fPsiRP>TMath::Pi()){fPsiRP-=TMath::Pi();}
  if(fPsiRP<0){fPsiRP+=TMath::Pi();}
  
  // reactionplaneangle + Pi() is the same angle
  if(TMath::Abs(fDeltaPsiRP)>TMath::Pi()/2){
    if(fDeltaPsiRP>0)fDeltaPsiRP-=TMath::Pi();
    else fDeltaPsiRP+=TMath::Pi();
  }
  
  Double_t cos2deltaRP=TMath::Cos(2*fDeltaPsiRP);
  
  // FillHistograms
  fh2RPSubevents->Fill(fPsiRP1,fPsiRP2);
  fh1RP->Fill(fPsiRP);
  fh2RPCentrality->Fill(fCentrality,fPsiRP);
  fh2RPDeltaRP->Fill(fDeltaPsiRP,fCentrality);
  fh2RPQxQy->Fill(mQx,mQy);
  fh2RPCosDeltaRP->Fill(fCentrality,cos2deltaRP);
  
  fRPAngle=fPsiRP;  
  return kTRUE;
}

Double_t AliAnalysisTaskJetServices::GetPhiWeight(Double_t phi,Double_t signedpt){
  if(!fh3PhiWeights)return 1;
  else return fh3PhiWeights->GetBinContent(fh3PhiWeights->GetXaxis()->FindBin(fCentrality),fh3PhiWeights->GetYaxis()->FindBin(signedpt),fh3PhiWeights->GetZaxis()->FindBin(phi));
}

 //________________________________________________________________________

Int_t  AliAnalysisTaskJetServices::GetListOfTracks(TList *list){
  Int_t iCount = 0;
  AliAODEvent *aod = 0;
  if(fUseAODInput)aod = dynamic_cast<AliAODEvent*>(InputEvent());
  else aod = AODEvent();
  if(!aod){
    return iCount;
  }
  for(int it = 0;it < aod->GetNumberOfTracks();++it){
    AliAODTrack *tr = aod->GetTrack(it);
    if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask)))continue;
    if(TMath::Abs(tr->Eta())>fTrackRecEtaWindow)continue;
    if(tr->Pt()<fMinTrackPt)continue;
    list->Add(tr);
    iCount++;
  }
  list->Sort();
  return iCount;

}

