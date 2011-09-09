
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
#include "AliAODVZERO.h"
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

AliAODVZERO*  AliAnalysisTaskJetServices::fgAODVZERO = NULL;
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
  fRPMethod(0),
  fCollisionType(kPbPb),
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
  fPsiVZEROA(0),
  fPsiVZEROC(0),
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
  fp1RPXA(0x0),
  fp1RPYA(0x0),
  fp1RPXC(0x0),
  fp1RPYC(0x0),
  fp1CalibRPXA(0x0),
  fp1CalibRPYA(0x0),
  fp1CalibRPXC(0x0),
  fp1CalibRPYC(0x0),
  fh2RPAC(0x0),
  fh2RPAT(0x0),
  fh2RPCT(0x0),
  fh2XYA(0x0),
  fh2XYC(0x0),
  fh2RPCentrality(0x0),
  fh2RPACentrality(0x0),
  fh2RPCCentrality(0x0),
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
  fFilterAODCollisions(kFALSE),
  fPhysicsSelectionFlag(AliVEvent::kMB),
  fSelectionInfoESD(0),
  fEventCutInfoESD(0),
  fFilterMask(0),
  fRPMethod(0),
  fCollisionType(kPbPb),
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
  fPsiVZEROA(0),
  fPsiVZEROC(0),
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
  fp1RPXA(0x0),
  fp1RPYA(0x0),
  fp1RPXC(0x0),
  fp1RPYC(0x0),
  fp1CalibRPXA(0x0),
  fp1CalibRPYA(0x0),
  fp1CalibRPXC(0x0),
  fp1CalibRPYC(0x0),
  fh2RPAC(0x0),
  fh2RPAT(0x0),
  fh2RPCT(0x0),
  fh2XYA(0x0),
  fh2XYC(0x0),
  fh2RPCentrality(0x0),
  fh2RPACentrality(0x0),
  fh2RPCCentrality(0x0),
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
  

  fh1CentralityESD = new TH1F("fh1CentralityESD","cent",103,-1,102);
  fHistList->Add(fh1CentralityESD);
  
  fh1Centrality = new TH1F("fh1Centrality","cent",103,-1,102);
  fHistList->Add(fh1Centrality);

  fh1RP = new TH1F("fh1RP","RP;#Psi",440, -1.*TMath::Pi(), 2.*TMath::Pi());
  fHistList->Add(fh1RP);

  fh2TriggerCount = new TH2F("fh2TriggerCount",";Trigger No.;constrained;Count",6,-0.5,5.5,kConstraints,-0.5,kConstraints-0.5); 
  fHistList->Add(fh2TriggerCount);

  fh2ESDTriggerCount = new TH2F("fh2ESDTriggerCount",";Trigger No.;constrained;Count",6,-0.5,5.5,kConstraints,-0.5,kConstraints-0.5); 
  fHistList->Add(fh2ESDTriggerCount);
  const Int_t nBins = 6*kConstraints;
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
  fHistList->Add(fh1NCosmicsPerEvent),


  fh2RPCentrality = new TH2F("fh2RPCentrality" ,"Reaction Plane Angle from tracks" , 20, 0.,100., 180, 0, TMath::Pi());
  fHistList->Add(fh2RPCentrality);


  fh2RPACentrality = new TH2F("fh2RPACentrality" ,"Reaction Plane Angle from vzero A" , 20, 0.,100., 180, 0, TMath::Pi());
  fHistList->Add(fh2RPACentrality);

  fh2RPCCentrality = new TH2F("fh2RPCCentrality" ,"Reaction Plane Angle from vzero C" , 20, 0.,100., 180, 0, TMath::Pi());
  fHistList->Add(fh2RPCCentrality);

  fh2RPAC = new TH2F("fh2RPAC" ,"Reaction Plane Angle vzero a vs c" , 180, 0, TMath::Pi(), 180, 0, TMath::Pi());
  fHistList->Add( fh2RPAC);

  fh2RPAT = new TH2F("fh2RPAT" ,"Reaction Plane Angle vzero a vs tracks" , 180, 0, TMath::Pi(), 180, 0, TMath::Pi());
  fHistList->Add( fh2RPAT);

  fh2RPCT = new TH2F("fh2RPCT" ,"Reaction Plane Angle vzero c vs tracks" , 180, 0, TMath::Pi(), 180, 0, TMath::Pi());
  fHistList->Add( fh2RPCT);

  fh2XYA = new TH2F("fh2XYA" ,"XY vzeroa ;X;Y;" ,100,-0.3,0.3,100,-0.3,0.3);
  fHistList->Add(fh2XYA);
  fh2XYC = new TH2F("fh2XYC" ,"XY vzeroc ;X;Y;" ,100,-0.3,0.3,100,-0.3,0.3);
  fHistList->Add(fh2XYC);

  // profiles for mean 
  fp1RPXA = new TProfile("fp1RPXA","mean vzeroa x vs run number;run;x",(Int_t)(1+fRunRange[1]-fRunRange[0]),fRunRange[0]-0.5,fRunRange[1]+0.5);
  fHistList->Add(fp1RPXA);

  fp1RPYA = new TProfile("fp1RPYA","mean vzeroa y vs run number;run;y",(Int_t)(1+fRunRange[1]-fRunRange[0]),fRunRange[0]-0.5,fRunRange[1]+0.5);
  fHistList->Add(fp1RPYA);


  fp1RPXC = new TProfile("fp1RPXC","mean vzeroc x vs run number;run;x",(Int_t)(1+fRunRange[1]-fRunRange[0]),fRunRange[0]-0.5,fRunRange[1]+0.5);
  fHistList->Add(fp1RPXC);

  fp1RPYC = new TProfile("fp1RPYC","mean vzeroa y vs run number;run;y",(Int_t)(1+fRunRange[1]-fRunRange[0]),fRunRange[0]-0.5,fRunRange[1]+0.5);
  fHistList->Add(fp1RPYC);


  TH1::AddDirectory(oldStatus);

  // Add an AOD branch for replication
  if(fNonStdFile.Length()){
     if (fDebug > 1) AliInfo("Replicating header");
     fgAODHeader = new AliAODHeader;
     AddAODBranch("AliAODHeader",&fgAODHeader,fNonStdFile.Data());
     if (fDebug > 1) AliInfo("Replicating vzeros");
     fgAODVZERO = new AliAODVZERO;
     AddAODBranch("AliAODVZERO",&fgAODVZERO,fNonStdFile.Data());
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

  Float_t cent = 0;
  if(fCollisionType==kPbPb){
    if(aod)cent = aod->GetHeader()->GetCentrality();
    if(fDebug)Printf("%s:%d %3.3f",(char*)__FILE__,__LINE__,cent);
    if(cent<0)cent = 101;
  }
  fCentrality = cent;
  fRPAngle = 0;
  fPsiVZEROA = 0;
  fPsiVZEROC = 0;


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
      if(aodH&&cand&&fFilterAODCollisions&&!esd){
	if(fCentrality<=80&&aodVtxIn){
	  aodH->SetFillAOD(kTRUE);
	  aodH->SetFillExtension(kTRUE);
	}
      }

      TList recTracks;
      GetListOfTracks(&recTracks);
      CalculateReactionPlaneAngleVZERO(aod);
      fRPAngle = aod->GetHeader()->GetEventplane();
      fh1RP->Fill(fRPAngle);
      fh2RPCentrality->Fill(fCentrality,fRPAngle);
      fh2RPACentrality->Fill(fCentrality,fPsiVZEROA);
      fh2RPCCentrality->Fill(fCentrality,fPsiVZEROC);
      fh2RPAC->Fill(fPsiVZEROA,fPsiVZEROC);
      fh2RPAT->Fill(fPsiVZEROA,fRPAngle);
      fh2RPCT->Fill(fPsiVZEROC,fRPAngle);
      if(fRPMethod==kRPTracks)AliAnalysisHelperJetTasks::ReactionPlane(kTRUE,fRPAngle); // set slection to false
      else if(fRPMethod==kRPVZEROA)AliAnalysisHelperJetTasks::ReactionPlane(kTRUE,fPsiVZEROA); // set slection to false
      else if(fRPMethod==kRPVZEROC)AliAnalysisHelperJetTasks::ReactionPlane(kTRUE,fPsiVZEROA); // set slection to false

      if(fUseAODInput&&fCentrality<=80){
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
      Double_t q[kRPMethods] = {fRPAngle,fPsiVZEROA,fPsiVZEROC};
      fgAODHeader->SetQTheta(q,kRPMethods);
    }
    if (fgAODVZERO){
      *fgAODVZERO =  *(dynamic_cast<AliAODVZERO*>(aod->GetVZEROData()));
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
  if(fgAODVZERO)delete fgAODVZERO;
  delete fp1CalibRPXA;
  delete fp1CalibRPYA;
  delete fp1CalibRPXC;
  delete fp1CalibRPYC;

}


void AliAnalysisTaskJetServices::SetV0Centroids(TProfile *xa,TProfile *ya,TProfile *xc, TProfile *yc){

  if(xa){
    if(fp1CalibRPXA)delete fp1CalibRPXA;
    fp1CalibRPXA =  (TProfile*)xa->Clone(Form("%sCalib",xa->GetName()));
  }
  else{
    Printf("%s:%d centroid histogram is 0x0",(char*)__FILE__,__LINE__);
  }
  if(ya){
    if(fp1CalibRPYA)delete fp1CalibRPYA;
    fp1CalibRPYA =  (TProfile*)ya->Clone(Form("%sCalib",ya->GetName()));
  }
  else{
    Printf("%s:%d centroid histogram is 0x0",(char*)__FILE__,__LINE__);
  }
  if(xc){
    if(fp1CalibRPXC)delete fp1CalibRPXC;
    fp1CalibRPXC =  (TProfile*)xc->Clone(Form("%sCalib",xc->GetName()));
  }
  else{
    Printf("%s:%d centroid histogram is 0x0",(char*)__FILE__,__LINE__);
  }
  if(ya){
    if(fp1CalibRPYC)delete fp1CalibRPYC;
    fp1CalibRPYC =  (TProfile*)yc->Clone(Form("%sCalib",yc->GetName()));
  }
  else{
    Printf("%s:%d centroid histogram is 0x0",(char*)__FILE__,__LINE__);
  }
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
  
  if(!vtx){
    Printf("%s:%d No ESD vertex found",(char*)__FILE__,__LINE__);
    return kFALSE;
  }
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

  if(!vtx){
    Printf("%s:%d No AOD vertex found",(char*)__FILE__,__LINE__);
    return kFALSE;
  }


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

  Float_t cent = 0;
  if(fCollisionType==kPbPb){
    if(esd->GetCentrality()){
      cent = esd->GetCentrality()->GetCentralityPercentile("V0M");
    }
    if(cent<0)cent = 101;
    if(cent>80||cent<0)return 5;
    if(cent>50)return 4;
    if(cent>30)return 3;
    if(cent>10)return 2;
    return 1;
  }
  return 1;
}


Int_t AliAnalysisTaskJetServices::GetEventClass(AliAODEvent *aod){

  if(fCollisionType==kPbPb){
    Float_t cent = aod->GetHeader()->GetCentrality();
    if(cent>80||cent<0)return 5;
    if(cent>50)return 4;
    if(cent>30)return 3;
    if(cent>10)return 2;
    return 1;
  }
  return 1;

}


void AliAnalysisTaskJetServices::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
}

Bool_t AliAnalysisTaskJetServices::CalculateReactionPlaneAngleVZERO(AliAODEvent *aod){

  //  const Double_t arr_eta[11]={-3.7, -3.2, -2.7, -2.2, -1.7,0, 2.8, 3.4, 3.9, 4.5,5.1};
  if(!aod)return kFALSE;
  static bool bFirst = true;
  static Double_t v0phi[64] = {0,};

  if(bFirst){
    int is=0;
    for(int iArm = 0; iArm<2; iArm++){
      for(int iRing = 0; iRing<4 ; iRing++){
	for(int iSec = 0; iSec<8 ; iSec++){
	  v0phi[is] = 22.5 + 45. * iSec;
	  v0phi[is] *= TMath::Pi()/180; 
	  // cout<< is <<" "<< v0phi[is]<<endl;
	  is++;
	}
      }
    }
    bFirst = false;
  }

  // 
  const AliAODVZERO *aodVZERO = aod->GetVZEROData();
  Double_t numYZNA = 0,numXZNA = 0,sumZNA = 0;
  Double_t meanXA = 0,meanYA = 0;

  Double_t numYZNC = 0,numXZNC = 0,sumZNC = 0;
  Double_t meanXC = 0,meanYC = 0;



  static Int_t iOldRun = -1;
  static Int_t iFoundBin = -1;

  if(aod->GetRunNumber()!=iOldRun&&(fp1CalibRPYA)){
    // search only or the bin in case of new runs
    iFoundBin = -1;
    Int_t ib = fp1CalibRPYA->FindBin(aod->GetRunNumber());
    Float_t err = fp1CalibRPYA->GetBinError(ib);
    if(err>0){// value can be zero...
      iFoundBin = ib;
    }
    else{
      Int_t ibLo = ib-1;
      Int_t ibUp = ib+1;
      while(iFoundBin<0&&(ibLo>0||ibUp<=fp1CalibRPYA->GetNbinsX())){
	err = fp1CalibRPYA->GetBinError(ibLo);
	if(err>0){
	  iFoundBin = ibLo;
	}
	else{
	  err = fp1CalibRPYA->GetBinError(ibUp);
	  if(err>0)iFoundBin = ibUp;
	}
	ibUp++;
	ibLo--;
      }
    }
    iOldRun = aod->GetRunNumber();
  }

  if(fDebug)Printf("%s:%d iFoundBin %d",(char*)__FILE__,__LINE__,iFoundBin);

  if(iFoundBin>0&&(fp1CalibRPYA)){
    meanXA = fp1CalibRPXA->GetBinContent(iFoundBin);
    meanYA = fp1CalibRPYA->GetBinContent(iFoundBin);
    meanXC = fp1CalibRPXC->GetBinContent(iFoundBin);
    meanYC = fp1CalibRPYC->GetBinContent(iFoundBin);
  }

  if(fDebug)Printf("%s:%d iFoundBin %1.3E %1.3E %1.3E %1.3E",(char*)__FILE__,__LINE__,meanXA,meanYA,meanXC,meanYC);

  for (int i=0; i<64; i++) {  
    Double_t mult = aodVZERO->GetMultiplicity(i);
    Double_t phi = v0phi[i];
    if (mult>0) {
      if (i<32) { //C-side
        Double_t wZNC= mult;
	numYZNC += sin(2.*phi)*wZNC; 
        numXZNC += cos(2.*phi)*wZNC;
	sumZNC+=wZNC;
      }
      else if(i>31){ //A-side
	Double_t wZNA=mult;
	numYZNA += sin(2.*phi)*wZNA; 
        numXZNA += cos(2.*phi)*wZNA;
	sumZNA+=wZNA; 
      } 
    }// mult>0
  }// 64 sectors

  Double_t   XC = numXZNC/sumZNC; 
  Double_t   YC = numYZNC/sumZNC; 
  
  Double_t   XA = numXZNA/sumZNA;
  Double_t   YA = numYZNA/sumZNA;
  

  fPsiVZEROA = 0.5*TMath::ATan2(YA-meanYA, XA-meanXA);
  if(fPsiVZEROA>TMath::Pi()){fPsiVZEROA-=TMath::Pi();}
  if(fPsiVZEROA<0){fPsiVZEROA+=TMath::Pi();}

  fPsiVZEROC = 0.5*TMath::ATan2(YC-meanYC, XA-meanXC);
  if(fPsiVZEROC>TMath::Pi()){fPsiVZEROC-=TMath::Pi();}
  if(fPsiVZEROC<0){fPsiVZEROC+=TMath::Pi();}
  
  fh2XYA->Fill(XA-meanXA,YA-meanYA); // control
  fp1RPXA->Fill(aod->GetRunNumber(),XA);
  fp1RPYA->Fill(aod->GetRunNumber(),YA);
  fh2XYC->Fill(XC-meanXC,YC-meanYC); // control
  fp1RPXC->Fill(aod->GetRunNumber(),XC);
  fp1RPYC->Fill(aod->GetRunNumber(),YC);
  return kTRUE;

}

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

