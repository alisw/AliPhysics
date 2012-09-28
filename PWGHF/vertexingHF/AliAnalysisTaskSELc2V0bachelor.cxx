/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appeuear in the supporting documentation. The authors make no claims   *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//
//
//               Base class for Lc2V0 Analysis
//
//
//  The Lc spectra study is done 2D histograms:
//   cascade_invariantMass VS cascade_pT
//
//  Cuts have been centralized in AliRDHFCutsLctoV0 class
//
//-------------------------------------------------------------------------
//
//                 Authors: A.De Caro(a,b), P. Pagano(b)
//  (a) Centro 'E.Fermi' - Roma
//  (b) INFN and University of Salerno
//
//  Contatcs: decaro@sa.infn.it
//            paola.pagano@sa.infn.it
//-------------------------------------------------------------------------

#include <TSystem.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TH1F.h>
#include <TH2F.h>
#include "TROOT.h"
#include <TDatabasePDG.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSELc2V0bachelor.h"
#include "AliNormalizationCounter.h"
#include "AliAODPidHF.h"
#include "AliPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliInputEventHandler.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskSELc2V0bachelor)

//__________________________________________________________________________
AliAnalysisTaskSELc2V0bachelor::AliAnalysisTaskSELc2V0bachelor():  
  AliAnalysisTaskSE(),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fOutputPIDBach(0),
  fCEvents(0),
  fPIDResponse(0),
  fIsK0sAnalysis(kFALSE),
  fCounter(0),
  fProdCuts(0),
  fAnalCuts(0),
  fListCuts(0),
  fUseOnTheFlyV0(kFALSE),
  fIsEventSelected(kFALSE)
{
  //
  // Default ctor
  //
}
//___________________________________________________________________________
AliAnalysisTaskSELc2V0bachelor::AliAnalysisTaskSELc2V0bachelor(const Char_t* name, AliRDHFCutsLctoV0* prodCuts,
							       AliRDHFCutsLctoV0* analCuts, Bool_t useOnTheFly) :
  AliAnalysisTaskSE(name),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fOutputPIDBach(0),
  fCEvents(0),
  fPIDResponse(0),
  fIsK0sAnalysis(kFALSE),
  fCounter(0),
  fProdCuts(prodCuts),
  fAnalCuts(analCuts),
  fListCuts(0),
  fUseOnTheFlyV0(useOnTheFly),
  fIsEventSelected(kFALSE)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSELc2V0bachelor","Calling Constructor");

  DefineOutput(1,TList::Class());  //conters
  DefineOutput(2,TList::Class());  //All Entries output
  DefineOutput(3,TList::Class());  //3sigma PID output
  DefineOutput(4,AliNormalizationCounter::Class());
  DefineOutput(5,TList::Class());

}

//___________________________________________________________________________
AliAnalysisTaskSELc2V0bachelor::~AliAnalysisTaskSELc2V0bachelor() {
  //
  // destructor
  //
  Info("~AliAnalysisTaskSELc2V0bachelor","Calling Destructor");
  
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }

  if (fOutputAll) {
    delete fOutputAll;
    fOutputAll = 0;
  }

  if (fOutputPIDBach) {
    delete fOutputPIDBach;
    fOutputPIDBach = 0;
  }

  if (fPIDResponse) {
    delete  fPIDResponse;
  }

  if (fCounter) {
    delete fCounter;
    fCounter = 0;
  }

  if (fProdCuts) {
    delete fProdCuts;
    fProdCuts = 0;
  }

  if (fAnalCuts) {
    delete fAnalCuts;
    fAnalCuts = 0;
  }

  if (fListCuts) {
    delete fListCuts;
    fListCuts = 0;
  }

}
//_________________________________________________
void AliAnalysisTaskSELc2V0bachelor::Init() {
  //
  // Initialization
  //

  fIsEventSelected=kFALSE;

  if (fDebug > 1) AliInfo("Init");

  fListCuts = new TList();
  fListCuts->SetOwner();
  fListCuts->Add(new AliRDHFCutsLctoV0(*fProdCuts));
  fListCuts->Add(new AliRDHFCutsLctoV0(*fAnalCuts));
  PostData(5,fListCuts);

  return;
}

//_________________________________________________
void AliAnalysisTaskSELc2V0bachelor::UserExec(Option_t *)
{
  // user exec
  if (!fInputEvent) {
    AliError("NO EVENT FOUND!");
    return;
  }

  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
  TClonesArray *arrayLctopKos=0;

  if (!aodEvent && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
   
    if (aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      arrayLctopKos=(TClonesArray*)aodFromExt->GetList()->FindObject("CascadesHF");
    }
  } else {
    arrayLctopKos=(TClonesArray*)aodEvent->GetList()->FindObject("CascadesHF");
  }

  fCEvents->Fill(1);
  fCounter->StoreEvent(aodEvent,fAnalCuts,fUseMCInfo);

  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();
  if (!vtx1) return;

  // fix for temporary bug in ESDfilter 
  if (TMath::Abs(aodEvent->GetMagneticField())<0.001) return;
  fCEvents->Fill(2);

  Float_t zVertex = vtx1->GetZ();
  ((TH1F*)(fOutput->FindObject("hZ2")))->Fill(zVertex);

  if (vtx1->GetNContributors()<1) return;
  fCEvents->Fill(3);

  ((TH1F*)(fOutput->FindObject("hZ3")))->Fill(zVertex);

  if (!arrayLctopKos) {
    AliInfo("Could not find array of HF cascades, skipping the event");
    return;
  } else {
    if (arrayLctopKos->GetEntriesFast()) {
      AliInfo(Form("Found %d cascades",arrayLctopKos->GetEntriesFast()));
    }
  }
  fCEvents->Fill(4);
  ((TH1F*)(fOutput->FindObject("hZ4")))->Fill(zVertex);

  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD
  /*
  TString trigclass = aodEvent->GetFiredTriggerClasses();
  if (trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD")||trigclass.Contains("C0SMH-B-NOPF-ALL"))
    fCEvents->Fill(5); // in case of RealData events
  */

  ///////////////////////
  Bool_t check1 = kFALSE;
  TString firedTriggerClasses = aodEvent->GetFiredTriggerClasses(); // trigger class
  if ( !fUseMCInfo && // don't do for MC...
       (aodEvent->GetRunNumber()<136851 || aodEvent->GetRunNumber()>139517) ) { // ...and for PbPb 2010 data
    if ( !(firedTriggerClasses.Contains("CINT1")) ) {
      AliInfo(Form(" ======================== firedTriggerClasses.Data() = %s",firedTriggerClasses.Data()));
      fCEvents->Fill(8);
      ((TH1F*)(fOutput->FindObject("hZ8")))->Fill(zVertex);
      check1 = kTRUE;
    }
  }

  ULong64_t fTriggerMask=AliVEvent::kAnyINT;
  Bool_t isSelectedAAA = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);
  if (!isSelectedAAA) {
    fCEvents->Fill(9);
    ((TH1F*)(fOutput->FindObject("hZ9")))->Fill(zVertex);
  }

  if (!isSelectedAAA || check1) {
    fCEvents->Fill(16);
    ((TH1F*)(fOutput->FindObject("hZ16")))->Fill(zVertex);
  }

  fTriggerMask=AliVEvent::kAny;
  Bool_t isSelectedBBB = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);
  if (!isSelectedBBB) {
    fCEvents->Fill(10);
    ((TH1F*)(fOutput->FindObject("hZ10")))->Fill(zVertex);
  }

  TString title=vtx1->GetTitle();
  if (title.Contains("Z")) {
    fCEvents->Fill(11);
    ((TH1F*)(fOutput->FindObject("hZ11")))->Fill(zVertex);
  }
  else if (title.Contains("3D")) {
    fCEvents->Fill(12);
    ((TH1F*)(fOutput->FindObject("hZ12")))->Fill(zVertex);
  } else {
    fCEvents->Fill(13);
    ((TH1F*)(fOutput->FindObject("hZ13")))->Fill(zVertex);
  }

  if (TMath::Abs(zVertex)<=fAnalCuts->GetMaxVtxZ()) {
    fCEvents->Fill(14);
    ((TH1F*)(fOutput->FindObject("hZ14")))->Fill(zVertex);
  }

  fIsEventSelected = fAnalCuts->IsEventSelected(aodEvent);
  if ( fIsEventSelected ) {
    fCEvents->Fill(7);
    ((TH1F*)(fOutput->FindObject("hZ7")))->Fill(zVertex);
  } else {
    fCEvents->Fill(15);
    ((TH1F*)(fOutput->FindObject("hZ15")))->Fill(zVertex);
    return;
  }
  ///////////////////////



  // mc analysis 
  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader=0;

  if (fUseMCInfo) {
    // MC array need for maching
    mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!mcArray) {
      AliError("Could not find Monte-Carlo in AOD");
      return;
    }
    fCEvents->Fill(5); // in case of MC events
    ((TH1F*)(fOutput->FindObject("hZ5")))->Fill(zVertex);

    // load MC header
    mcHeader = (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!mcHeader) {
      AliError("AliAnalysisTaskSELc2V0bachelor::UserExec: MC header branch not found!\n");
      return;
    }
    fCEvents->Fill(6);
    ((TH1F*)(fOutput->FindObject("hZ6")))->Fill(zVertex);

    //AliInfo("~~~~~~~~~~Sono dentro fUseMCInfo 2");

    // check on MC Lc Daughter
    SearchLcDaughter(mcArray);

  }

  //AliInfo("~~~~~~~~~~Sono prima di isEvSelA");
  Int_t nSelectedProd = 0;
  Int_t nSelectedAnal = 0;
  if (fIsK0sAnalysis) {
    MakeAnalysisForLc2prK0S(vtx1,arrayLctopKos,mcArray,
			    nSelectedProd, fProdCuts, nSelectedAnal, fAnalCuts);

    if (nSelectedAnal) {

      ((TH1I*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(4);
      ((TH1F*)(fOutput->FindObject("hZ4a")))->Fill(zVertex);

      //TString firedTriggerClasses = aodEvent->GetFiredTriggerClasses(); // trigger class
      Bool_t check1a = kFALSE;
      if ( !fUseMCInfo && // don't do for MC...
	   (aodEvent->GetRunNumber()<136851 || aodEvent->GetRunNumber()>139517) ) { // ...and for PbPb 2010 data
	if ( !(firedTriggerClasses.Contains("CINT1")) ) {
	  AliInfo(Form(" ======================== firedTriggerClasses.Data() = %s",firedTriggerClasses.Data()));
	  ((TH1I*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(8);
	  ((TH1F*)(fOutput->FindObject("hZ8a")))->Fill(zVertex);
	  check1a = kTRUE;
	}
      }

      //ULong64_t fTriggerMask=AliVEvent::kAnyINT;
      fTriggerMask=AliVEvent::kAnyINT;
      Bool_t isSelectedAAAa = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);
      if (!isSelectedAAAa) {
	((TH1I*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(9);
	((TH1F*)(fOutput->FindObject("hZ9a")))->Fill(zVertex);
      }

      if (!isSelectedAAAa || check1a) {
	((TH1I*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(16);
	((TH1F*)(fOutput->FindObject("hZ16a")))->Fill(zVertex);
      }

      fTriggerMask=AliVEvent::kAny;
      Bool_t isSelectedBBBa = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);
      if (!isSelectedBBBa) {
	((TH1I*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(10);
	((TH1F*)(fOutput->FindObject("hZ10a")))->Fill(zVertex);
      }

      //TString title=vtx1->GetTitle();
      if (title.Contains("Z")) {
	((TH1I*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(11);
	((TH1F*)(fOutput->FindObject("hZ11a")))->Fill(zVertex);
      }
      else if (title.Contains("3D")) {
	((TH1I*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(12);
	((TH1F*)(fOutput->FindObject("hZ12a")))->Fill(zVertex);
      } else {
	((TH1I*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(13);
	((TH1F*)(fOutput->FindObject("hZ13a")))->Fill(zVertex);
      }

      if (TMath::Abs(zVertex)<=fAnalCuts->GetMaxVtxZ()) {
	((TH1I*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(14);
	((TH1F*)(fOutput->FindObject("hZ14a")))->Fill(zVertex);
      }

      //Bool_t eventSelected = fAnalCuts->IsEventSelected(aodEvent);
      if ( fIsEventSelected ) {
	((TH1I*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(7);
	((TH1F*)(fOutput->FindObject("hZ7a")))->Fill(zVertex);
      } else {
	((TH1I*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(15);
	((TH1F*)(fOutput->FindObject("hZ15a")))->Fill(zVertex);
      }
    }

  }


  fCounter->StoreCandidates(aodEvent,nSelectedProd,kTRUE);
  fCounter->StoreCandidates(aodEvent,nSelectedAnal,kFALSE);

  PostData(1,fOutput);
  PostData(2,fOutputAll);
  PostData(3,fOutputPIDBach);
  PostData(4,fCounter);

}
//________________________________________ terminate ___________________________
void AliAnalysisTaskSELc2V0bachelor::Terminate(Option_t*)
{    
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  //AliInfo("Terminate","");
  AliAnalysisTaskSE::Terminate();
  
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    AliError("fOutput not available");
    return;
  }
  
  //fCEvents = dynamic_cast<TH1I*>(fOutput->FindObject("fCEvents"));

  fOutputAll = dynamic_cast<TList*> (GetOutputData(2));
  if (!fOutputAll) {
    AliError("fOutputAll not available");
    return;
  }

  fOutputPIDBach = dynamic_cast<TList*> (GetOutputData(3));
  if (!fOutputPIDBach) {
    AliError("fOutputPIDBach not available");
    return;
  }

  return;
}
//___________________________________________________________________________
void AliAnalysisTaskSELc2V0bachelor::UserCreateOutputObjects() { 
  // output
  AliInfo(Form("CreateOutputObjects of task %s\n", GetName()));
  
  //slot #1  
  //OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("chist0");

  fOutputAll = new TList();
  fOutputAll->SetOwner();
  fOutputAll->SetName("listAll");

  fOutputPIDBach = new TList();
  fOutputPIDBach->SetOwner();
  fOutputPIDBach->SetName("listPIDBach");

  // define histograms
  DefineHistograms();
  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  if (fProdCuts->GetIsUsePID()){
    fProdCuts->GetPidHF()->SetPidResponse(fPIDResponse);
    fProdCuts->GetPidV0pos()->SetPidResponse(fPIDResponse);
    fProdCuts->GetPidV0neg()->SetPidResponse(fPIDResponse);
    fProdCuts->GetPidHF()->SetOldPid(kFALSE);
    fProdCuts->GetPidV0pos()->SetOldPid(kFALSE);
    fProdCuts->GetPidV0neg()->SetOldPid(kFALSE);
  }
  if (fAnalCuts->GetIsUsePID()){
    fAnalCuts->GetPidHF()->SetPidResponse(fPIDResponse);
    fAnalCuts->GetPidV0pos()->SetPidResponse(fPIDResponse);
    fAnalCuts->GetPidV0neg()->SetPidResponse(fPIDResponse);
    fAnalCuts->GetPidHF()->SetOldPid(kFALSE);
    fAnalCuts->GetPidV0pos()->SetOldPid(kFALSE);
    fAnalCuts->GetPidV0neg()->SetOldPid(kFALSE);
  }

  PostData(1,fOutput);
  PostData(2,fOutputAll);
  PostData(3,fOutputPIDBach);

  fCounter = new AliNormalizationCounter("NormalizationCounter");
  fCounter->Init();
  PostData(4,fCounter);
  
  return;
}
//___________________________________ hiostograms _______________________________________
void  AliAnalysisTaskSELc2V0bachelor::DefineHistograms() {

  fCEvents = new TH1I("fCEvents","conter",17,0,17);
  fCEvents->SetStats(kTRUE);
  fCEvents->GetXaxis()->SetBinLabel(2,"Analyzed events");
  fCEvents->GetXaxis()->SetBinLabel(3,"AliAODVertex exists");
  fCEvents->GetXaxis()->SetBinLabel(4,"GetNContributors()>0");
  fCEvents->GetXaxis()->SetBinLabel(5,"CascadesHF exists");
  fCEvents->GetXaxis()->SetBinLabel(9,"triggerClass!=CINT1");
  fCEvents->GetXaxis()->SetBinLabel(10,"triggerMask!=kAnyINT");
  fCEvents->GetXaxis()->SetBinLabel(11,"triggerMask!=kAny");
  fCEvents->GetXaxis()->SetBinLabel(12,"vtxTitle.Contains(Z)");
  fCEvents->GetXaxis()->SetBinLabel(13,"vtxTitle.Contains(3D)");
  fCEvents->GetXaxis()->SetBinLabel(14,"vtxTitle.Doesn'tContain(Z-3D)");
  fCEvents->GetXaxis()->SetBinLabel(15,"zVtx<=10cm");
  fCEvents->GetXaxis()->SetBinLabel(16,"!IsEventSelected");
  fCEvents->GetXaxis()->SetBinLabel(17,"triggerMask!=kAnyINT || triggerClass!=CINT1");
  fCEvents->GetXaxis()->SetBinLabel(8,"IsEventSelected");
  //fCEvents->GetXaxis()->SetTitle("");
  fCEvents->GetYaxis()->SetTitle("counts");

  fOutput->Add(fCEvents);
  TString fillthis="";

  if (fUseMCInfo) {
    fillthis="histMcStatLc";
    TH1F* mcStatisticLc = new TH1F(fillthis.Data(),"#Lambda_{C} generated and their decays",19,-9.5,9.5);
    fOutput->Add(mcStatisticLc);
  }

  fillthis="histopionV0SigmaVspTOF";
  TH2F *hpionV0SigmaVspTOF=new TH2F(fillthis.Data(),fillthis.Data(),100,0.,5.0,100,-10.0,10.0);
  fillthis="histoprotonBachSigmaVspTOF";
  TH2F *hprotonBachSigmaVspTOF=new TH2F(fillthis.Data(),fillthis.Data(),100,0.,5.0,100,-10.0,10.0);

  fOutput->Add(hpionV0SigmaVspTOF);
  fOutput->Add(hprotonBachSigmaVspTOF);

  fillthis="histopionV0SigmaVspTPC";
  TH2F *hpionV0SigmaVspTPC=new TH2F(fillthis.Data(),fillthis.Data(),100,0.,5.0,100,-10.0,10.0);
  fillthis="histoprotonBachSigmaVspTPC";
  TH2F *hprotonBachSigmaVspTPC=new TH2F(fillthis.Data(),fillthis.Data(),100,0.,5.0,100,-10.0,10.0);

  fOutput->Add(hpionV0SigmaVspTPC);
  fOutput->Add(hprotonBachSigmaVspTPC);

  TH1F *hZ2 = new TH1F("hZ2","",100,-50.,50.);
  fOutput->Add(hZ2);
  TH1F *hZ3 = new TH1F("hZ3","",100,-50.,50.);
  fOutput->Add(hZ3);
  TH1F *hZ4 = new TH1F("hZ4","",100,-50.,50.);
  fOutput->Add(hZ4);
  TH1F *hZ5 = new TH1F("hZ5","",100,-50.,50.);
  fOutput->Add(hZ5);
  TH1F *hZ6 = new TH1F("hZ6","",100,-50.,50.);
  fOutput->Add(hZ6);
  TH1F *hZ7 = new TH1F("hZ7","",100,-50.,50.);
  fOutput->Add(hZ7);
  TH1F *hZ8 = new TH1F("hZ8","",100,-50.,50.);
  fOutput->Add(hZ8);
  TH1F *hZ9 = new TH1F("hZ9","",100,-50.,50.);
  fOutput->Add(hZ9);
  TH1F *hZ10 = new TH1F("hZ10","",100,-50.,50.);
  fOutput->Add(hZ10);
  TH1F *hZ11 = new TH1F("hZ11","",100,-50.,50.);
  fOutput->Add(hZ11);
  TH1F *hZ12 = new TH1F("hZ12","",100,-50.,50.);
  fOutput->Add(hZ12);
  TH1F *hZ13 = new TH1F("hZ13","",100,-50.,50.);
  fOutput->Add(hZ13);
  TH1F *hZ14 = new TH1F("hZ14","",100,-50.,50.);
  fOutput->Add(hZ14);
  TH1F *hZ15 = new TH1F("hZ15","",100,-50.,50.);
  fOutput->Add(hZ15);
  TH1F *hZ16 = new TH1F("hZ16","",100,-50.,50.);
  fOutput->Add(hZ16);

  TH1I *hCandidateSelection = new TH1I("hCandidateSelection","",12,0,12);
  hCandidateSelection->GetXaxis()->SetBinLabel(1,"IsEventSelected");
  hCandidateSelection->GetXaxis()->SetBinLabel(2,"IsSecondaryVtx");
  hCandidateSelection->GetXaxis()->SetBinLabel(3,"V0toPosNeg");
  hCandidateSelection->GetXaxis()->SetBinLabel(4,"prodCuts::kTracks");
  hCandidateSelection->GetXaxis()->SetBinLabel(5,"prodCuts::kCandidate");
  hCandidateSelection->GetXaxis()->SetBinLabel(6,"prodCuts::kPID");
  hCandidateSelection->GetXaxis()->SetBinLabel(7,"prodCuts::kAll");
  hCandidateSelection->GetXaxis()->SetBinLabel(8,"analCuts::kTracks");
  hCandidateSelection->GetXaxis()->SetBinLabel(9,"analCuts::kCandidate");
  hCandidateSelection->GetXaxis()->SetBinLabel(10,"analCuts::kPID");
  hCandidateSelection->GetXaxis()->SetBinLabel(11,"analCuts::kAll");
  fOutput->Add(hCandidateSelection);

  TH1I *hEventsWithCandidates = new TH1I("hEventsWithCandidates","conter",17,0,17);
  fOutput->Add(hEventsWithCandidates);

  TH1F *hZ2a = new TH1F("hZ2a","",100,-50.,50.);
  fOutput->Add(hZ2a);
  TH1F *hZ3a = new TH1F("hZ3a","",100,-50.,50.);
  fOutput->Add(hZ3a);
  TH1F *hZ4a = new TH1F("hZ4a","",100,-50.,50.);
  fOutput->Add(hZ4a);
  TH1F *hZ5a = new TH1F("hZ5a","",100,-50.,50.);
  fOutput->Add(hZ5a);
  TH1F *hZ6a = new TH1F("hZ6a","",100,-50.,50.);
  fOutput->Add(hZ6a);
  TH1F *hZ7a = new TH1F("hZ7a","",100,-50.,50.);
  fOutput->Add(hZ7a);
  TH1F *hZ8a = new TH1F("hZ8a","",100,-50.,50.);
  fOutput->Add(hZ8a);
  TH1F *hZ9a = new TH1F("hZ9a","",100,-50.,50.);
  fOutput->Add(hZ9a);
  TH1F *hZ10a = new TH1F("hZ10a","",100,-50.,50.);
  fOutput->Add(hZ10a);
  TH1F *hZ11a = new TH1F("hZ11a","",100,-50.,50.);
  fOutput->Add(hZ11a);
  TH1F *hZ12a = new TH1F("hZ12a","",100,-50.,50.);
  fOutput->Add(hZ12a);
  TH1F *hZ13a = new TH1F("hZ13a","",100,-50.,50.);
  fOutput->Add(hZ13a);
  TH1F *hZ14a = new TH1F("hZ14a","",100,-50.,50.);
  fOutput->Add(hZ14a);
  TH1F *hZ15a = new TH1F("hZ15a","",100,-50.,50.);
  fOutput->Add(hZ15a);
  TH1F *hZ16a = new TH1F("hZ16a","",100,-50.,50.);
  fOutput->Add(hZ16a);

  TH1I *hSwitchOnCandidates0 = new TH1I("hSwitchOnCandidates0","",8,0,8);
  fOutput->Add(hSwitchOnCandidates0);
  TH1I *hSwitchOnCandidates1 = new TH1I("hSwitchOnCandidates1","",8,0,8);
  fOutput->Add(hSwitchOnCandidates1);
  TH1I *hSwitchOnCandidates2 = new TH1I("hSwitchOnCandidates2","",8,0,8);
  fOutput->Add(hSwitchOnCandidates2);
  TH1I *hSwitchOnCandidates3 = new TH1I("hSwitchOnCandidates3","",8,0,8);
  fOutput->Add(hSwitchOnCandidates3);
  TH1I *hSwitchOnCandidates4 = new TH1I("hSwitchOnCandidates4","",8,0,8);
  fOutput->Add(hSwitchOnCandidates4);

  if (fIsK0sAnalysis) DefineK0SHistos();// hK0S histos declarations

  return;
}
//________________________________________________________________________
void AliAnalysisTaskSELc2V0bachelor::FillLc2pK0Sspectrum(AliAODRecoCascadeHF *part,
						   Int_t isLc,
						   Int_t &nSelectedProd,
						   AliRDHFCutsLctoV0 *cutsProd,
						   Int_t &nSelectedAnal,
						   AliRDHFCutsLctoV0 *cutsAnal)
{
  //
  // Fill histos for Lc -> K0S+proton
  //
  
  if ( ( ( (cutsProd->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) )
    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(4);
  if ( ( ( (cutsProd->IsSelected(part,AliRDHFCuts::kPID))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) )
    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(5);
  if ( ( ( (cutsProd->IsSelected(part,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) ) {
    nSelectedProd++;
    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(6);
  }
  /*
  if ( ( (cutsProd->IsSelectedPID(part))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)  && // to be used since fUsePID is kFALSE for cutsProd
       ( ( (cutsProd->IsSelectedSingleCut(part,AliRDHFCuts::kCandidate,2))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) ) nSelectedProd++;
  */


  if ( (((cutsAnal->IsSelected(part,AliRDHFCuts::kTracks))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) )
    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(7);

  //Bool_t usePID = cutsAnal->GetIsUsePID();
  //cutsAnal->SetUsePID(kFALSE);
  if ( (((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) )
    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(8);
  //cutsAnal->SetUsePID(usePID);
  if ( (((cutsAnal->IsSelected(part,AliRDHFCuts::kPID))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) )
    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(9);

  TString fillthis="";

  Double_t invmassLc = part->InvMassLctoK0sP();
  Double_t lambdacpt=part->Pt();
  Double_t cosOAK0Sp = part->PxProng(0)*part->PxProng(1)+
    part->PyProng(0)*part->PyProng(1)+
    part->PzProng(0)*part->PzProng(1);
  cosOAK0Sp /= (part->PProng(0)*part->PProng(1));

  AliAODv0 * v0part = (AliAODv0*)part->Getv0();
  Bool_t onFlyV0 = v0part->GetOnFlyStatus(); // on-the-flight V0s
  Double_t momK0s  = TMath::Sqrt(v0part->Ptot2V0());
  Double_t ptK0s = TMath::Sqrt(v0part->Pt2V0());
  Double_t dcaV0 = v0part->DcaV0ToPrimVertex();
  Double_t invmassK0s = v0part->MassK0Short();
  Bool_t isInV0window = (((cutsAnal->IsSelectedSingleCut(part,AliRDHFCuts::kCandidate,2))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)); // cut on V0 invMass

  AliAODTrack *bachelor = (AliAODTrack*)part->GetBachelor();
  Double_t momBach  = bachelor->P();
  Bool_t isBachelorID = (((cutsAnal->IsSelected(part,AliRDHFCuts::kPID))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)); // ID x bachelor

  //if (isBachelorID && isInV0window) nSelectedAnal++;
  if ( (((cutsAnal->IsSelected(part,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) ) {
    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(10);
    if (!onFlyV0) {
      nSelectedAnal++;
      ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(11);
    }
  }

  if ( !onFlyV0 ) {

    ((TH1I*)(fOutput->FindObject("hSwitchOnCandidates0")))->Fill( cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate) );

    if ( (((cutsAnal->IsSelected(part,AliRDHFCuts::kTracks))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) )
      ((TH1I*)(fOutput->FindObject("hSwitchOnCandidates1")))->Fill( cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate) );

    if ( (((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) )
      ((TH1I*)(fOutput->FindObject("hSwitchOnCandidates2")))->Fill( cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate) );

    if ( (((cutsAnal->IsSelected(part,AliRDHFCuts::kPID))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) )
      ((TH1I*)(fOutput->FindObject("hSwitchOnCandidates3")))->Fill( cutsAnal->IsSelected(part,AliRDHFCuts::kPID) );

    if ( (((cutsAnal->IsSelected(part,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) )
      ((TH1I*)(fOutput->FindObject("hSwitchOnCandidates4")))->Fill( cutsAnal->IsSelected(part,AliRDHFCuts::kAll) );

  }

  if (onFlyV0 && fUseOnTheFlyV0) {  

    fillthis="histK0SMass";
    ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassK0s,ptK0s);
    if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassK0s,ptK0s);

    if (isInV0window) {

      fillthis="histpK0Svsp";
      ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momBach,momK0s);
      if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momBach,momK0s);

      fillthis="histcosOAK0Spvsp";
      ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(cosOAK0Sp,lambdacpt);
      if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(cosOAK0Sp,lambdacpt);

      fillthis="histoDCAtoPVvsinvmassK0s";
      ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassK0s,dcaV0);
      if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassK0s,dcaV0);

      fillthis="histArmPodK0s";
      FillArmPodDistribution(v0part,fillthis,fOutputAll);
      if (isBachelorID) FillArmPodDistribution(v0part,fillthis,fOutputPIDBach);

      fillthis="histLcMassByK0S";
      ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);
      if (isBachelorID)((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);
 
    }
  }
  else if (!onFlyV0) {

    fillthis="histK0SMassOffline";
    ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassK0s,ptK0s);
    if (isBachelorID) ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassK0s,ptK0s);

    if (isInV0window) {

      fillthis="histpK0SvspOffline";
      ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momBach,momK0s);
      if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momBach,momK0s);

      fillthis="histcosOAK0SpvspOffline";
      ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(cosOAK0Sp,lambdacpt);
      if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(cosOAK0Sp,lambdacpt);
  
      fillthis="histoDCAtoPVvsinvmassK0sOffline";
      ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassK0s,dcaV0);
      if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassK0s,dcaV0);

      fillthis="histOfflineArmPodK0s";
      FillArmPodDistribution(v0part,fillthis,fOutputAll);
      if (isBachelorID) FillArmPodDistribution(v0part,fillthis,fOutputPIDBach);

      fillthis="histLcMassOfflineByK0S";
      ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);
      if (isBachelorID)((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);

      if (fIsEventSelected) {
	fillthis="hist1LcMassOfflineByK0S";
	if (isBachelorID) ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);
      } else {
	fillthis="hist0LcMassOfflineByK0S";
	if (isBachelorID) ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);
      }
 
    }
  }


  if (fUseMCInfo) {
    if (isLc==1) {
      if (onFlyV0 && fUseOnTheFlyV0) {

	fillthis="histK0SMassSgn";
	((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassK0s,ptK0s);
	if (isBachelorID) ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassK0s,ptK0s);

	if (isInV0window) {

	  fillthis="histpK0SvspSgn";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momBach,momK0s);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momBach,momK0s);

	  fillthis="histcosOAK0SpSgnvsp";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(cosOAK0Sp,lambdacpt);
	  if (isBachelorID)     ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(cosOAK0Sp,lambdacpt);

	  fillthis="histoDCAtoPVvsinvmassK0sSgn";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassK0s,dcaV0);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassK0s,dcaV0);

	  fillthis="histArmPodK0sSgn";
	  FillArmPodDistribution(v0part,fillthis,fOutputAll);
	  if (isBachelorID) FillArmPodDistribution(v0part,fillthis,fOutputPIDBach);

	  fillthis="histLcMassByK0SSgn";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);

	}
    
      }     
      else if (!onFlyV0) {  

	fillthis="histK0SMassOfflineSgn";
	((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassK0s,ptK0s);
	if (isBachelorID) ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassK0s,ptK0s);

 	if (isInV0window) {

	  fillthis="histpK0SvspOfflineSgn";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momBach,momK0s);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momBach,momK0s);
  
	  fillthis="histcosOAK0SpSgnvspOffline";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(cosOAK0Sp,lambdacpt);
	  if (isBachelorID)     ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(cosOAK0Sp,lambdacpt);

	  fillthis="histoDCAtoPVvsinvmassK0sOfflineSgn";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassK0s,dcaV0);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassK0s,dcaV0);
   
	  fillthis="histOfflineArmPodK0sSgn";
	  FillArmPodDistribution(v0part,fillthis,fOutputAll);
	  if (isBachelorID) FillArmPodDistribution(v0part,fillthis,fOutputPIDBach);

	  fillthis="histLcMassOfflineByK0SSgn";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);
	  
	}

      }

    }// sgn
    else { // bkg
      if (onFlyV0 && fUseOnTheFlyV0) {

	fillthis="histK0SMassBkg";
	((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassK0s,ptK0s);
	if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassK0s,ptK0s);

	if (isInV0window) {

	  fillthis="histpK0SvspBkg";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momBach,momK0s);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momBach,momK0s);

	  fillthis="histcosOAK0SpBkgvsp";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(cosOAK0Sp,lambdacpt);
	  if (isBachelorID)   ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(cosOAK0Sp,lambdacpt);

	  fillthis="histoDCAtoPVvsinvmassK0sBkg";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassK0s,dcaV0);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassK0s,dcaV0);

	  fillthis="histArmPodK0sBkg";
	  FillArmPodDistribution(v0part,fillthis,fOutputAll);
	  if (isBachelorID) FillArmPodDistribution(v0part,fillthis,fOutputPIDBach);

	  fillthis="histLcMassByK0SBkg";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);
 
	}

      }
      else if (!onFlyV0) {

	fillthis="histK0SMassOfflineBkg";
	((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassK0s,ptK0s);
	if (isBachelorID) ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassK0s,ptK0s);
  
	if (isInV0window) {

	  fillthis="histpK0SvspOfflineBkg";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momBach,momK0s);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momBach,momK0s);
   
	  fillthis="histcosOAK0SpBkgvspOffline";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(cosOAK0Sp,lambdacpt);
	  if (isBachelorID)   ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(cosOAK0Sp,lambdacpt);
    
	  fillthis="histoDCAtoPVvsinvmassK0sOfflineBkg";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassK0s,dcaV0);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassK0s,dcaV0);

	  fillthis="histOfflineArmPodK0sBkg";
	  FillArmPodDistribution(v0part,fillthis,fOutputAll);
	  if (isBachelorID) FillArmPodDistribution(v0part,fillthis,fOutputPIDBach);

	  fillthis="histLcMassOfflineByK0SBkg";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);
 
	}

      }

    }
  } // if fUseMCInfo
 
  return;
}
//-------------------------------------------------------------------------------
Int_t AliAnalysisTaskSELc2V0bachelor::CheckOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const {		
  //
  // checking whether the mother of the particles come from a charm or a bottom quark
  //
	
  Int_t pdgGranma = 0;
  Int_t abspdgGranma =0;
  Bool_t isFromB=kFALSE;
  Bool_t isQuarkFound=kFALSE;
  Int_t mother = mcPartCandidate->GetMother();
  Int_t istep = 0;
  while (mother>0) {
    istep++;
    AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcGranma) {
      pdgGranma = mcGranma->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ( (abspdgGranma > 500  && abspdgGranma < 600 ) ||
	   (abspdgGranma > 5000 && abspdgGranma < 6000) ) isFromB=kTRUE;
      else if (abspdgGranma==4 || abspdgGranma==5) isQuarkFound=kTRUE;
      mother = mcGranma->GetMother();
    } else {
      AliError("Failed casting the mother particle!");
      break;
    }
  }
  
  if (isFromB) return 5;
  else return 4;

}

//-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2V0bachelor::MakeAnalysisForLc2prK0S(AliAODVertex */*vtx1*/,
						       TClonesArray *arrayLctopKos,
						       TClonesArray *mcArray,
						       Int_t &nSelectedProd,
						       AliRDHFCutsLctoV0 *cutsProd,
						       Int_t &nSelectedAnal,
						       AliRDHFCutsLctoV0 *cutsAnal)
{

  // counters for efficiencies
  Int_t icountReco = 0;

  //Lc prong needed to MatchToMC method

  Int_t pdgCand = 4122;
  Int_t pdgDgLctoV0bachelorOld[2]={2212,310};
  Int_t pdgDgLctoV0bachelor[2]={310,2212};
  Int_t pdgDgV0toDaughters[2]={211,211};

  // loop over cascades to search for candidates Lc->p+K0S
  Int_t nCascades= arrayLctopKos->GetEntriesFast();
  if (nCascades==0) {
    AliInfo("Could not find cascades, skipping the event");
    return;
  }
  for (Int_t iLctopK0s = 0; iLctopK0s<nCascades; iLctopK0s++) {

    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(0);

    // Lc candidates and K0s from Lc
    AliAODRecoCascadeHF* lcK0spr = (AliAODRecoCascadeHF*)arrayLctopKos->At(iLctopK0s);
    if (!lcK0spr->GetSecondaryVtx()) {
      AliInfo("No secondary vertex");
      continue;
    }

    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(1);

    AliAODTrack * v0Pos = lcK0spr->Getv0PositiveTrack();
    AliAODTrack * v0Neg = lcK0spr->Getv0NegativeTrack();
    if (v0Pos->Charge() ==  v0Neg->Charge()) continue;
  
    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(2);

    Int_t isLc = 0;

    if (fUseMCInfo) {

      Bool_t isPrimary=kTRUE;
  
      Int_t pdgCode=-2;

      // find associated MC particle for Lc -> p+K0 and K0S->pi+pi
      Int_t mcLabelOld = MatchToMC(lcK0spr,pdgDgLctoV0bachelorOld,pdgDgV0toDaughters,mcArray);
      Int_t mcLabel = lcK0spr->MatchToMC(pdgCand,pdgDgLctoV0bachelor[0],pdgDgLctoV0bachelor,pdgDgV0toDaughters,mcArray,kTRUE);
      if (mcLabelOld!=mcLabel) AliInfo(Form(" Changed MC label: oldONE=%d wrt rightONE=%d",mcLabelOld,mcLabel));
      if (mcLabel>=0) {
	AliInfo(Form(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cascade numero %d di %d", iLctopK0s,nCascades));

	AliAODMCParticle *partLc = (AliAODMCParticle*)mcArray->At(mcLabel);
	Int_t checkOrigin = CheckOrigin(mcArray,partLc);
	if (checkOrigin==5) isPrimary=kFALSE;

	pdgCode = partLc->GetPdgCode();
	if (pdgCode<0) AliInfo(Form(" ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ %d ~~~~~~~~~~ %d", mcLabel, pdgCode));
	pdgCode = TMath::Abs(pdgCode);
	isLc = 1;
      } else {
	pdgCode=-1;
      }
    }

    if ( !( ( (cutsProd->IsSelected(lcK0spr,AliRDHFCuts::kTracks))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) ) continue;
    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(3);

    FillLc2pK0Sspectrum(lcK0spr, isLc,
			nSelectedProd, cutsProd,
			nSelectedAnal, cutsAnal);

  }
  
  AliDebug(2, Form("Found %d Reco particles that are Lc!!", icountReco));

}
//-------------------------------------------------------------------------------
Int_t AliAnalysisTaskSELc2V0bachelor::MatchToMC(AliAODRecoCascadeHF *lc2bacV0,
					  Int_t *pdgDgLc2bacV0, Int_t *pdgDgV0,
					  TClonesArray *mcArray)
{

  // bachelor
  AliAODTrack *bachelor = (AliAODTrack*)lc2bacV0->GetBachelor();
  Int_t labBachelor = bachelor->GetLabel();
  if (labBachelor<0) return -1;
  AliAODMCParticle *partBachelor = (AliAODMCParticle*)mcArray->At(labBachelor);
  if (TMath::Abs(partBachelor->GetPdgCode())!=pdgDgLc2bacV0[0]) return -1;

  Int_t labBacMother = partBachelor->GetMother();
  if (labBacMother<0) return -1;
  AliAODMCParticle *partBacMother = (AliAODMCParticle*)mcArray->At(labBacMother);
  if (TMath::Abs(partBacMother->GetPdgCode())!=4122) return -1;

  // V0
  AliAODTrack *posV0Daugh = (AliAODTrack*)lc2bacV0->Getv0PositiveTrack();
  AliAODTrack *negV0Daugh = (AliAODTrack*)lc2bacV0->Getv0NegativeTrack();
  Int_t labV0pos = posV0Daugh->GetLabel();
  Int_t labV0neg = negV0Daugh->GetLabel();

  if (labV0pos<0 || labV0neg<0) return -1;
  AliAODMCParticle *partV0pos = (AliAODMCParticle*)mcArray->At(labV0neg);
  AliAODMCParticle *partV0neg = (AliAODMCParticle*)mcArray->At(labV0pos);

  if ( ! ( (TMath::Abs(partV0pos->GetPdgCode())==pdgDgV0[0] &&
	    TMath::Abs(partV0neg->GetPdgCode())==pdgDgV0[1]) ||
	   (TMath::Abs(partV0pos->GetPdgCode())==pdgDgV0[1] &&
	    TMath::Abs(partV0neg->GetPdgCode())==pdgDgV0[0]) ) ) return -1;
  Int_t labV0posMother = partV0pos->GetMother();
  Int_t labV0negMother = partV0neg->GetMother();

  if (labV0posMother<0 || labV0negMother<0) return -1;
  if (labV0posMother!=labV0negMother) return -1;

  AliAODMCParticle *motherV0 = (AliAODMCParticle*)mcArray->At(labV0posMother);

  if (TMath::Abs(motherV0->GetPdgCode())!=pdgDgLc2bacV0[1]) return -1;
  Int_t labV0mother = motherV0->GetMother();
  if (labV0mother<0) return -1;
  AliAODMCParticle *gMotherV0 = (AliAODMCParticle*)mcArray->At(labV0mother);

  if ( !(pdgDgLc2bacV0[1]==310 && TMath::Abs(gMotherV0->GetPdgCode())==311) &&
       !(pdgDgLc2bacV0[1]==3122 && TMath::Abs(motherV0->GetPdgCode())==3122) ) return -1;

  if ( (pdgDgLc2bacV0[1]==310 && TMath::Abs(gMotherV0->GetPdgCode())==311) ) {
    Int_t labV0GMother = gMotherV0->GetMother();
    if (labV0GMother<0) return -1;
    AliAODMCParticle *ggMotherV0 = (AliAODMCParticle*)mcArray->At(labV0GMother);

    if (TMath::Abs(ggMotherV0->GetPdgCode())!=4122) return -1;
    gMotherV0 = (AliAODMCParticle*)ggMotherV0;
    labV0mother=labV0GMother;
  }
  else if (pdgDgLc2bacV0[1]==3122 && TMath::Abs(motherV0->GetPdgCode())==3122) {
    if (TMath::Abs(gMotherV0->GetPdgCode())!=4122) return -1;
  }

  if (labBacMother!=labV0mother) {
    AliInfo(Form("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++A - %d %d",
		 labBacMother, labV0mother));
    return -1;
  }

  AliInfo(Form(" V0: %d (%d) ->  %d (%d) + %d (%d) --- GM %d (PDG=%d)",
	       motherV0->GetPdgCode(), labV0posMother,
	       partV0pos->GetPdgCode(), labV0pos,
	       partV0neg->GetPdgCode(), labV0neg,
	       labV0mother, gMotherV0->GetPdgCode()));

  AliInfo(Form(" Bachelor: %d (%d) --- M %d (%d)",
	       partBachelor->GetPdgCode(), labBachelor,
	       labBacMother, partBacMother->GetPdgCode()));
  return labBacMother;//labV0mother;//

}

//-----------------------
void AliAnalysisTaskSELc2V0bachelor::SearchLcDaughter(TClonesArray *arrayMC)
{
 
  AliAODMCParticle *searchLc=0;
  AliAODMCParticle *daugh=0;
  AliAODMCParticle *daugh1=0;
  AliAODMCParticle *daugh2=0;
  AliAODMCParticle *daughK0=0;
  AliAODMCParticle *daughK0s1=0;
  AliAODMCParticle *daughK0s2=0;
  AliAODMCParticle *daughL1=0;
  AliAODMCParticle *daughL2=0;

  Int_t nDaughLc=0;
  Int_t nDaughK0=0;
  Int_t nDaughK0s=0;
  Int_t searchLcpdg=0;
  Int_t daughPdg1=0;
  Int_t daughPdg2=0;
  Int_t daughK0Pdg=0;
  Int_t nDaughL=0;
  Int_t daughK0s1pdg;
  Int_t daughK0s2pdg;
  Int_t daughL1pdg=0;
  Int_t daughL2pdg=0;

  TString fillthis="";
  fillthis="histMcStatLc";

  for (Int_t iii=0; iii<arrayMC->GetEntries(); iii++) {
    searchLc = (AliAODMCParticle*)arrayMC->At(iii);
    searchLcpdg =  searchLc->GetPdgCode();
    if (TMath::Abs(searchLcpdg) == 4122) {
      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(0);
      nDaughLc= searchLc->GetNDaughters();

      if (searchLcpdg == 4122) { // It is Lc+
	((TH1F*)(fOutput->FindObject(fillthis)))->Fill(1);
	if (nDaughLc!=2) continue;
	daugh1 = (AliAODMCParticle*)arrayMC->At(searchLc->GetDaughter(0));
	daughPdg1=daugh1->GetPdgCode();
	daugh2 = (AliAODMCParticle*)arrayMC->At(searchLc->GetDaughter(1));
	daughPdg2=daugh2->GetPdgCode();
	if ( (daughPdg1==2212 && daughPdg2==-311) ||
	     (daughPdg2==2212 && daughPdg1==-311) ) { // Lc+ -> p K0bar
	  AliInfo(Form(" Ecco %d(%d) -> %d(%d) %d(%d)",
		       iii,searchLcpdg,
		       searchLc->GetDaughter(0),daughPdg1,
		       searchLc->GetDaughter(1),daughPdg2));
	  ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(2);
	  if (daughPdg1==-311) {
	    nDaughK0=daugh1->GetNDaughters();
	    if (nDaughK0!=1) {
	      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(4);
	      continue;
	    } else {
	      daughK0 = (AliAODMCParticle*)arrayMC->At(daugh1->GetDaughter(0)); // K0S
	    }
	  } else {
	    nDaughK0=daugh2->GetNDaughters();
	    if (nDaughK0!=1) {
	      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(4);
	      continue;
	    } else {
	      daughK0 = (AliAODMCParticle*)arrayMC->At(daugh2->GetDaughter(0)); // K0S
	    }
	  }
	  if (!daughK0) {
	    //((TH1F*)(fOutput->FindObject(fillthis)))->Fill(4);
	    continue;
	  }
	  //cout << " positive daughK0 = " << daughK0 << endl;
	  AliInfo(" Found positive daughK0 ");
	  daughK0Pdg=daughK0->GetPdgCode();
	  if (daughK0Pdg!=310) {
	    ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(4);
	    continue;
	  } else {
	    ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(3);
	    nDaughK0s=daughK0->GetNDaughters();
	    if (nDaughK0s!=2) {
	      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(6);
	      continue;
	    }
	    daughK0s1= (AliAODMCParticle*)arrayMC->At(daughK0->GetDaughter(0));
	    daughK0s2= (AliAODMCParticle*)arrayMC->At(daughK0->GetDaughter(1));
	    daughK0s1pdg=daughK0s1->GetPdgCode();
	    daughK0s2pdg=daughK0s2->GetPdgCode();

	    if ( ((daughK0s1pdg==211) && (daughK0s2pdg==-211)) ||
		 ((daughK0s2pdg==211) && (daughK0s1pdg==-211)) ) ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(5);
	    else ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(6);
	  } // else if (daughK0Pdg==310)
	}//if ((daughPdg1==2212 && daughPdg2==-311)||(daughPdg2==2212 && daughPdg1==-311))

	else if ( (daughPdg1==3122 && daughPdg2==211) ||
		  (daughPdg2==3122 && daughPdg1==211) ) { // Lc+ -> pi+ Lambda
	  AliInfo(Form(" Ecco %d(%d) -> %d(%d) %d(%d)",
		       iii,searchLcpdg,
		       searchLc->GetDaughter(0),daughPdg1,
		       searchLc->GetDaughter(1),daughPdg2));
	  ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(7);
	  if (daughPdg1==3122)
	    daugh = (AliAODMCParticle*)arrayMC->At(searchLc->GetDaughter(0)); // Lambda
	  else
	    daugh = (AliAODMCParticle*)arrayMC->At(searchLc->GetDaughter(1)); // Lambda
	  if (!daugh) continue;
	  //cout << " positive daughL = " << daugh << endl;
	  AliInfo(" Found positive daughL ");
	  nDaughL=daugh->GetNDaughters();
	  if (nDaughL==2) {
	    daughL1= (AliAODMCParticle*)arrayMC->At(daugh->GetDaughter(0));
	    daughL2= (AliAODMCParticle*)arrayMC->At(daugh->GetDaughter(1));
	    daughL1pdg=daughL1->GetPdgCode();
	    daughL2pdg=daughL2->GetPdgCode();
	    if ( ((daughL1pdg==-211) && (daughL2pdg==2212)) ||
		 ((daughL2pdg==-211) && (daughL1pdg==2212)) ) ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(8);
	    else ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(9);
	  }//if (nDaughL==2) 
	  else ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(9);
	}//else if ((daughPdg1==3122 && daughPdg2==211)||(daughPdg2==3122 && daughPdg1==211))



      }//if (searchLcpdg == 4122)

      if (searchLcpdg == -4122) { // It is Lc+
	((TH1F*)(fOutput->FindObject(fillthis)))->Fill(-1);
	if (nDaughLc!=2) continue;
	daugh1 = (AliAODMCParticle*)arrayMC->At(searchLc->GetDaughter(0));
	daughPdg1=daugh1->GetPdgCode();
	daugh2 = (AliAODMCParticle*)arrayMC->At(searchLc->GetDaughter(1));
	daughPdg2=daugh2->GetPdgCode();
	if ( (daughPdg1==-2212 && daughPdg2==311) ||
	     (daughPdg2==-2212 && daughPdg1==311) ) { // Lc- -> pbar K0
	  ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(-2);
	  AliInfo(Form(" Ecco %d(%d) -> %d(%d) %d(%d)",
		       iii,searchLcpdg,
		       searchLc->GetDaughter(0),daughPdg1,
		       searchLc->GetDaughter(1),daughPdg2));
	  if (daughPdg1==311) {
	    nDaughK0=daugh1->GetNDaughters();
	    if (nDaughK0!=1) {
	      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(-4);
	      continue;
	    } else {
	      daughK0 = (AliAODMCParticle*)arrayMC->At(daugh1->GetDaughter(0));
	    }
	  } else {
	    nDaughK0=daugh2->GetNDaughters();
	    if (nDaughK0!=1) {
	      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(-4);
	      continue;
	    } else {
	      daughK0 = (AliAODMCParticle*)arrayMC->At(daugh2->GetDaughter(0));
	    }
	  }
	  if (!daughK0) {
            //((TH1F*)(fOutput->FindObject(fillthis)))->Fill(4);
	    continue;
	  }
	  //cout << " negative daughK0 = " << daughK0 << endl;
	  AliInfo(" Found negative daughK0 ");
	  daughK0Pdg=daughK0->GetPdgCode();
	  if (daughK0Pdg!=310) {
	    ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(-4);
	    continue;
	  } else {
	    ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(-3);
	    nDaughK0s=daughK0->GetNDaughters();
	    if (nDaughK0s!=2) {
	      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(-6);
	      continue;
	    }
	    daughK0s1= (AliAODMCParticle*)arrayMC->At(daughK0->GetDaughter(0));
	    daughK0s2= (AliAODMCParticle*)arrayMC->At(daughK0->GetDaughter(1));
	    daughK0s1pdg=daughK0s1->GetPdgCode();
	    daughK0s2pdg=daughK0s2->GetPdgCode();
	    if ( ((daughK0s1pdg==211) && (daughK0s2pdg==-211)) ||
		 ((daughK0s2pdg==211) && (daughK0s1pdg==-211)) ) ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(-5);
	    else ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(-6);
	  } // else if (daughK0Pdg==310)
	}//if ((daughPdg1==-2212 && daughPdg2==-311)||(daughPdg2==-2212 && daughPdg1==-311))

	else if ( (daughPdg1==-3122 && daughPdg2==-211) ||
		  (daughPdg2==-3122 && daughPdg1==-211) ) {
	  AliInfo(Form(" Ecco %d(%d) -> %d(%d) %d(%d)",
		       iii,searchLcpdg,
		       searchLc->GetDaughter(0),daughPdg1,
		       searchLc->GetDaughter(1),daughPdg2));
	  ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(-7);
	  if (daughPdg1==-3122)
	    daugh = (AliAODMCParticle*)arrayMC->At(searchLc->GetDaughter(0));
	  else
	    daugh = (AliAODMCParticle*)arrayMC->At(searchLc->GetDaughter(1));
	  if (!daugh) continue;
	  //cout << " negative daughL = " << daugh << endl;
	  AliInfo(" Found negative daughL ");
	  nDaughL=daugh->GetNDaughters();
	  if (nDaughL==2) {
	    daughL1 = (AliAODMCParticle*)arrayMC->At(daugh->GetDaughter(0));
	    daughL2 = (AliAODMCParticle*)arrayMC->At(daugh->GetDaughter(1));
	    daughL1pdg=daughL1->GetPdgCode();
	    daughL2pdg= daughL2->GetPdgCode();
	    if ( ((daughL1pdg==211) && (daughL2pdg==-2212)) ||
		 ((daughL2pdg==211) && (daughL1pdg==-2212)) ) ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(-8);
	    else ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(-9);
	  }//if (nDaughL==2) 
	  else ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(-9);
	}//else if ((daughPdg1==-3122 && daughPdg2==-211)||(daughPdg2==-3122 && daughPdg1==-211))


      }
    }// if (TMath::Abs(searchLcpdg) == 4122)
  }// for (Int_t iii=0; iii<arrayMC->GetEntries(); iii++)
  
}
//----------------------------------------------------

void AliAnalysisTaskSELc2V0bachelor::FillArmPodDistribution(AliAODv0 *vZero,
						      TString histoTitle,
						      TList *histoList) {

  Double_t alpha = vZero->AlphaV0();
  Double_t qT    = vZero->PtArmV0();

  ((TH2F*)(histoList->FindObject(histoTitle)))->Fill(alpha,qT);

}
//----------------------------------------------------

void AliAnalysisTaskSELc2V0bachelor::DefineK0SHistos()
{ 

  TString nameMass=" ", nameSgn=" ", nameBkg=" ";

  ///---------------- START  K0S HISTOS DECLARATIONS -------------------///

  if (fUseOnTheFlyV0) {

    // V0 invariant masses (on-the-fly)
    nameMass="histK0SMass";
    TH2F* spectrumK0SMass = new TH2F(nameMass.Data(),"K^{0}_{S} invariant mass VS p_{T}; M(K^{0}_{S}) [GeV/c^{2}]; p_{T}(K^{0}_{S}) [GeV/c]; Entries",
				     520,0.43,0.56,200,0.,20.);
    spectrumK0SMass->Sumw2();
    spectrumK0SMass->SetLineColor(6);
    spectrumK0SMass->SetMarkerStyle(20);
    spectrumK0SMass->SetMarkerSize(0.6);
    spectrumK0SMass->SetMarkerColor(6);

    // Lc invariant masses (x K0S on-the-fly)
    nameMass="histLcMassByK0S";
    TH2F* spectrumLcMassByK0S = new TH2F(nameMass.Data(),"#Lambda_{C} invariant mass (by K^{0}_{S}) vs p_{T} ; M(#Lambda_{C}) [GeV/c^{2}]; p_{T}",
					    1000,2.,2.5,200,0.,20.);
    spectrumLcMassByK0S->Sumw2();
    spectrumLcMassByK0S->SetLineColor(6);
    spectrumLcMassByK0S->SetMarkerStyle(20);
    spectrumLcMassByK0S->SetMarkerSize(0.6);
    spectrumLcMassByK0S->SetMarkerColor(6);

    nameMass="histcosOAK0Spvsp";
    TH2F* cosOpeningAngleK0Spvsp = new TH2F(nameMass.Data(),"#Lambda_{C}  : K^{0}_{S} - p Opening Angle vs p; Cos(Opening Angle)  ; p ",
					   200,-1.,1.,200,0.,20.);
 
    nameMass="histpK0Svsp";
    TH2F* momentumDistributionK0Svsp = new TH2F(nameMass.Data(),"#Lambda_{C}  : K^{0}_{S} vs p Total Momentum Distribution;  p_{p}; p_{K^{0}_{S}}  ",
					       200,0.,20.,200,0.,20.);

    nameMass="histArmPodK0s";
    TH2F* armenterosPodK0s = new TH2F(nameMass.Data(),"K^{0}_{S}  Armenteros-Podolanski distribution; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]",
				      200,-1.,1.,300,0.,0.3);
 
    nameMass="histoDCAtoPVvsinvmassK0s";
    TH2F *dcatoPVvspK0s = new TH2F(nameMass.Data(),"K^{0}_{S}: DCA to Primary Vertex vs  K^{0}_{S} invariant mass ; M(K^{0}_{S}) [GeV/c^{2}]; DCA to Primary Vertex []; Entries",
				   520,0.43,0.56,100,0.,10.0);

    TH2F* allspectrumK0SMass = (TH2F*)spectrumK0SMass->Clone(); 
    TH2F* allspectrumLcMassByK0S    = (TH2F*)spectrumLcMassByK0S->Clone(); 
    TH2F* allcosOpeningAngleK0Spvsp= (TH2F*)cosOpeningAngleK0Spvsp->Clone(); 
    TH2F* allmomentumDistributionK0Svsp= (TH2F*)momentumDistributionK0Svsp->Clone(); 
    TH2F* alldcatoPVvspK0s=(TH2F*)dcatoPVvspK0s->Clone(); 

    TH2F* pidBachspectrumK0SMass = (TH2F*)spectrumK0SMass->Clone(); 
    TH2F* pidBachspectrumLcMassByK0S    = (TH2F*)spectrumLcMassByK0S->Clone(); 
    TH2F* pidBachcosOpeningAngleK0Spvsp= (TH2F*)cosOpeningAngleK0Spvsp->Clone(); 
    TH2F* pidBachmomentumDistributionK0Svsp= (TH2F*)momentumDistributionK0Svsp->Clone(); 
    TH2F* pidBachdcatoPVvspK0s=(TH2F*)dcatoPVvspK0s->Clone(); 

    TH2F* allArmenterosPodK0s = (TH2F*)armenterosPodK0s->Clone();
    TH2F* pidBachArmenterosPodK0s = (TH2F*)armenterosPodK0s->Clone();

    fOutputAll->Add(allspectrumK0SMass);
    fOutputAll->Add(allspectrumLcMassByK0S);
    fOutputAll->Add(allcosOpeningAngleK0Spvsp); 
    fOutputAll->Add(allmomentumDistributionK0Svsp); 
    fOutputAll->Add(allArmenterosPodK0s);
    fOutputAll->Add(alldcatoPVvspK0s);

    fOutputPIDBach->Add(pidBachspectrumK0SMass);
    fOutputPIDBach->Add(pidBachspectrumLcMassByK0S);
    fOutputPIDBach->Add(pidBachcosOpeningAngleK0Spvsp); 
    fOutputPIDBach->Add(pidBachmomentumDistributionK0Svsp); 
    fOutputPIDBach->Add(pidBachArmenterosPodK0s);
    fOutputPIDBach->Add(pidBachdcatoPVvspK0s);
 
  }

  // V0 invariant masses (offline)
  nameMass="histK0SMassOffline";
  TH2F* spectrumK0SMassOffline = new TH2F(nameMass.Data(),"K^{0}_{S} invariant mass VS p_{T}; M(K^{0}_{S}) [GeV/c^{2}]; p_{T}(K^{0}_{S}) [GeV/c]; Entries",
					  520,0.43,0.56,200,0.,20.);
  spectrumK0SMassOffline->Sumw2();
  spectrumK0SMassOffline->SetLineColor(6);
  spectrumK0SMassOffline->SetMarkerStyle(20);
  spectrumK0SMassOffline->SetMarkerSize(0.6);
  spectrumK0SMassOffline->SetMarkerColor(6);


  // Lc invariant masses (x K0S offline)
  nameMass="histLcMassOfflineByK0S";
  TH2F* spectrumLcMassOfflineByK0S = new TH2F(nameMass.Data(),"#Lambda_{C} invariant mass (by K^{0}_{S})  vs p_{T}; M(#Lambda_{C}) [GeV/c^{2}]; p_{T}",
					      1000,2.,2.5,200,0.,20.);
  spectrumLcMassOfflineByK0S->Sumw2();
  spectrumLcMassOfflineByK0S->SetLineColor(6);
  spectrumLcMassOfflineByK0S->SetMarkerStyle(20);
  spectrumLcMassOfflineByK0S->SetMarkerSize(0.6);
  spectrumLcMassOfflineByK0S->SetMarkerColor(6);


  nameMass="histcosOAK0SpvspOffline";
  TH2F* cosOpeningAngleK0SpvspOffline = new TH2F(nameMass.Data(),"#Lambda_{C}  : K^{0}_{S} - p Opening Angle vs p  -  Offline ; Cos(Opening Angle)  ; p ",
						200,-1.,1.,200,0.,20.);

  nameMass="histpK0SvspOffline";
  TH2F* momentumDistributionK0SvspOffline = new TH2F(nameMass.Data(),"#Lambda_{C}  : K^{0}_{S} vs p Total Momentum Distribution - Offline ;  p_{p}; p_{K^{0}_{S}}  ",
						    200,0.,20.,200,0.,20.);

  nameMass="histOfflineArmPodK0s";
  TH2F* armenterosPodK0sOff = new TH2F(nameMass.Data(),"K^{0}_{S}  Armenteros-Podolanski distribution -offline-; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]",
				       200,-1.,1.,300,0.,0.3);

  nameMass="histoDCAtoPVvsinvmassK0sOffline";
  TH2F *dcatoPVvspK0sOffline = new TH2F(nameMass.Data(),"K^{0}_{S}: DCA to Primary Vertex vs  K^{0}_{S} invariant mass -offline -; M(K^{0}_{S}) [GeV/c^{2}]; DCA to Primary Vertex []; Entries",
					520,0.43,0.56,100,0.,10.0);



  TH2F* allspectrumK0SMassOffline = (TH2F*)spectrumK0SMassOffline->Clone(); 
  TH2F* allspectrumLcMassOfflineByK0S    = (TH2F*)spectrumLcMassOfflineByK0S->Clone(); 
  TH2F* allcosOpeningAngleK0SpvspOffline= (TH2F*)cosOpeningAngleK0SpvspOffline->Clone(); 
  TH2F* allmomentumDistributionK0SvspOffline= (TH2F*)momentumDistributionK0SvspOffline->Clone(); 
  TH2F* alldcatoPVvspK0sOffline=(TH2F*)dcatoPVvspK0sOffline->Clone(); 

  TH2F* pidBachspectrumK0SMassOffline = (TH2F*)spectrumK0SMassOffline->Clone(); 
  TH2F* pidBachspectrumLcMassOfflineByK0S    = (TH2F*)spectrumLcMassOfflineByK0S->Clone(); 
  TH2F* pidBachcosOpeningAngleK0SpvspOffline= (TH2F*)cosOpeningAngleK0SpvspOffline->Clone(); 
  TH2F* pidBachmomentumDistributionK0SvspOffline= (TH2F*)momentumDistributionK0SvspOffline->Clone(); 
  TH2F* pidBachdcatoPVvspK0sOffline=(TH2F*)dcatoPVvspK0sOffline->Clone(); 

  TH2F* allArmenterosPodK0sOff = (TH2F*)armenterosPodK0sOff->Clone();
  TH2F* pidBachArmenterosPodK0sOff = (TH2F*)armenterosPodK0sOff->Clone();


  fOutputAll->Add(allspectrumK0SMassOffline);
  fOutputAll->Add(allspectrumLcMassOfflineByK0S);
  fOutputAll->Add(allcosOpeningAngleK0SpvspOffline); 
  fOutputAll->Add(allmomentumDistributionK0SvspOffline); 
  fOutputAll->Add(allArmenterosPodK0sOff);
  fOutputAll->Add(alldcatoPVvspK0sOffline);

  fOutputPIDBach->Add(pidBachspectrumK0SMassOffline);
  fOutputPIDBach->Add(pidBachspectrumLcMassOfflineByK0S);
  fOutputPIDBach->Add(pidBachcosOpeningAngleK0SpvspOffline); 
  fOutputPIDBach->Add(pidBachmomentumDistributionK0SvspOffline); 
  fOutputPIDBach->Add(pidBachArmenterosPodK0sOff);
  fOutputPIDBach->Add(pidBachdcatoPVvspK0sOffline);


  nameMass="hist1LcMassOfflineByK0S";
  TH2D* h1 = new TH2D(nameMass.Data(),"IsEventSelected; M(#Lambda_{C}) [GeV/c^{2}]; p_{T}",1000,2.,2.5,200,0.,20.);
  fOutput->Add(h1);
  nameMass="hist0LcMassOfflineByK0S";
  TH2D* h0 = new TH2D(nameMass.Data(),"!IsEventSelected; M(#Lambda_{C}) [GeV/c^{2}]; p_{T}",1000,2.,2.5,200,0.,20.);
  fOutput->Add(h0);

  if (fUseMCInfo) {

    if (fUseOnTheFlyV0) {

      nameSgn="histK0SMassSgn";
      nameBkg="histK0SMassBkg";
      TH2F* spectrumK0SMassSgn = new TH2F(nameSgn.Data(), "K^{0}_{S} Signal invariant mass VS p_{T} - MC; M(K^{0}_{S}) [GeV/c^{2}]; p_{T}(K^{0}_{S}) [GeV/c]; Entries", 520,0.43,0.56,200,0.,20.);
      TH2F* spectrumK0SMassBkg = new TH2F(nameBkg.Data(), "K^{0}_{S} Background invariant mass VS p_{T} - MC; M(K^{0}_{S}) [GeV/c^{2}]; p_{T}(K^{0}_{S}) [GeV/c]; Entries",  520,0.43,0.56,200,0.,20.);
      spectrumK0SMassSgn->Sumw2();
      spectrumK0SMassBkg->Sumw2();
      spectrumK0SMassSgn->SetLineColor(2);
      spectrumK0SMassBkg->SetLineColor(4);
      spectrumK0SMassSgn->SetMarkerStyle(20);
      spectrumK0SMassBkg->SetMarkerStyle(20);
      spectrumK0SMassSgn->SetMarkerSize(0.6);
      spectrumK0SMassBkg->SetMarkerSize(0.6);
      spectrumK0SMassSgn->SetMarkerColor(2);
      spectrumK0SMassBkg->SetMarkerColor(4);

      TH2F* allspectrumK0SMassSgn = (TH2F*)spectrumK0SMassSgn->Clone(); 
      TH2F* allspectrumK0SMassBkg = (TH2F*) spectrumK0SMassBkg->Clone();  
      TH2F* pidBachspectrumK0SMassSgn = (TH2F*)spectrumK0SMassSgn->Clone(); 
      TH2F* pidBachspectrumK0SMassBkg = (TH2F*) spectrumK0SMassBkg->Clone();  
  

      fOutputAll->Add(allspectrumK0SMassSgn);
      fOutputAll->Add(allspectrumK0SMassBkg);
      fOutputPIDBach->Add(pidBachspectrumK0SMassSgn);
      fOutputPIDBach->Add(pidBachspectrumK0SMassBkg);


      nameSgn="histLcMassByK0SSgn";
      nameBkg="histLcMassByK0SBkg";
      TH2F* spectrumLcMassByK0SSgn = new TH2F(nameSgn.Data(), "#Lambda_{C} Signal invariant mass (by K^{0}_{S}) vs p_{T}  - MC; M(#Lambda_{C}) [GeV/c^{2}];  p_{T}",
					      1000,2.,2.5,200,0.,20.);
      TH2F* spectrumLcMassByK0SBkg = new TH2F(nameBkg.Data(), "#Lambda_{C} Background invariant mass (by K^{0}_{S}) vs p_{T}  - MC; M(#Lambda_{C}) [GeV/c^{2}]; p_{T}",
					      1000,2.,2.5,200,0.,20.);
      spectrumLcMassByK0SSgn->Sumw2();
      spectrumLcMassByK0SBkg->Sumw2();
      spectrumLcMassByK0SSgn->SetLineColor(2);
      spectrumLcMassByK0SBkg->SetLineColor(4);
      spectrumLcMassByK0SSgn->SetMarkerStyle(20);
      spectrumLcMassByK0SBkg->SetMarkerStyle(20);
      spectrumLcMassByK0SSgn->SetMarkerSize(0.6);
      spectrumLcMassByK0SBkg->SetMarkerSize(0.6);
      spectrumLcMassByK0SSgn->SetMarkerColor(2);
      spectrumLcMassByK0SBkg->SetMarkerColor(4);

      TH2F* allspectrumLcMassByK0SSgn = (TH2F*)spectrumLcMassByK0SSgn->Clone(); 
      TH2F* allspectrumLcMassByK0SBkg = (TH2F*) spectrumLcMassByK0SBkg->Clone();  
      TH2F* pidBachspectrumLcMassByK0SSgn = (TH2F*)spectrumLcMassByK0SSgn->Clone(); 
      TH2F* pidBachspectrumLcMassByK0SBkg = (TH2F*) spectrumLcMassByK0SBkg->Clone();  
      fOutputAll->Add(allspectrumLcMassByK0SSgn);
      fOutputAll->Add(allspectrumLcMassByK0SBkg);
      fOutputPIDBach->Add(pidBachspectrumLcMassByK0SSgn);
      fOutputPIDBach->Add(pidBachspectrumLcMassByK0SBkg);


      nameSgn="histcosOAK0SpSgnvsp";
      nameBkg="histcosOAK0SpBkgvsp";
      TH2F* cosOpeningAngleK0SpSgnvsp= new TH2F(nameSgn.Data(),"#Lambda_{C} SGN : K^{0}_{S} - p Opening Angle  vs p - MC ; Cos(Opening Angle);  p ",
						200,-1.,1.,200,0.,20.);
      TH2F* cosOpeningAngleK0SpBkgvsp= new TH2F(nameBkg.Data(),"#Lambda_{C} BKG : K^{0}_{S} - p Opening Angle  vs p - MC;  Cos(Opening Angle);  p ",
						200,-1.,1.,200,0.,20.);

      TH2F* allcosOpeningAngleK0SpSgnvsp= (TH2F*)cosOpeningAngleK0SpSgnvsp->Clone(); 
      TH2F* allcosOpeningAngleK0SpBkgvsp= (TH2F*)cosOpeningAngleK0SpBkgvsp->Clone(); 
      TH2F* pidBachcosOpeningAngleK0SpSgnvsp= (TH2F*)cosOpeningAngleK0SpSgnvsp->Clone(); 
      TH2F* pidBachcosOpeningAngleK0SpBkgvsp= (TH2F*)cosOpeningAngleK0SpBkgvsp->Clone(); 
      fOutputAll->Add(allcosOpeningAngleK0SpSgnvsp); 
      fOutputAll->Add(allcosOpeningAngleK0SpBkgvsp); 
      fOutputPIDBach->Add(pidBachcosOpeningAngleK0SpSgnvsp); 
      fOutputPIDBach->Add(pidBachcosOpeningAngleK0SpBkgvsp); 

      nameSgn="histpK0SvspSgn";
      nameBkg="histpK0SvspBkg";
      TH2F* momentumDistributionK0SvspSgn= new TH2F(nameSgn.Data(),"#Lambda_{C} SGN : K^{0}_{S} vs p Total Momentum Distribution - MC; p_{p}; p_{K^{0}_{S}}",
						    200,0.,20.,200,0.,20.);
      TH2F* momentumDistributionK0SvspBkg= new TH2F(nameBkg.Data(),"#Lambda_{C} BKG : K^{0}_{S} vs p Total Momentum Distribution - MC; p_{p}; p_{K^{0}_{S}}",
						    200,0.,20.,200,0.,20.);

      TH2F* allmomentumDistributionK0SvspSgn= (TH2F*)momentumDistributionK0SvspSgn->Clone(); 
      TH2F* allmomentumDistributionK0SvspBkg= (TH2F*)momentumDistributionK0SvspBkg->Clone(); 
      TH2F* pidBachmomentumDistributionK0SvspSgn= (TH2F*)momentumDistributionK0SvspSgn->Clone(); 
      TH2F* pidBachmomentumDistributionK0SvspBkg= (TH2F*)momentumDistributionK0SvspBkg->Clone(); 
      fOutputAll->Add(allmomentumDistributionK0SvspSgn); 
      fOutputAll->Add(allmomentumDistributionK0SvspBkg); 
      fOutputPIDBach->Add(pidBachmomentumDistributionK0SvspSgn); 
      fOutputPIDBach->Add(pidBachmomentumDistributionK0SvspBkg); 


      // armenteros-podolanski plots K0S
      nameSgn="histArmPodK0sSgn";
      nameBkg="histArmPodK0sBkg";
      TH2F* armenterosPodK0sSgn = new TH2F(nameSgn.Data(),"K^{0}_{S}  Armenteros-Podolanski distribution (sgn); #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]",
					   200,-1.,1.,300,0.,0.3);
      TH2F* armenterosPodK0sBkg = new TH2F(nameBkg.Data(),"K^{0}_{S}  Armenteros-Podolanski distribution (bkg); #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]",
					   200,-1.,1.,300,0.,0.3);

      TH2F* allArmenterosPodK0sSgn = (TH2F*)armenterosPodK0sSgn->Clone();
      TH2F* allArmenterosPodK0sBkg = (TH2F*)armenterosPodK0sBkg->Clone();
      TH2F* pidBachArmenterosPodK0sSgn = (TH2F*)armenterosPodK0sSgn->Clone();
      TH2F* pidBachArmenterosPodK0sBkg = (TH2F*)armenterosPodK0sBkg->Clone();

      fOutputAll->Add(allArmenterosPodK0sSgn);
      fOutputAll->Add(allArmenterosPodK0sBkg);

      fOutputPIDBach->Add(pidBachArmenterosPodK0sSgn);
      fOutputPIDBach->Add(pidBachArmenterosPodK0sBkg);


      nameSgn="histoDCAtoPVvsinvmassK0sSgn";
      nameBkg="histoDCAtoPVvsinvmassK0sBkg";
      TH2F *dcatoPVvspK0sSgn=new TH2F(nameSgn.Data(),"K^{0}_{S}: DCA to Primary Vertex vs  K^{0}_{S} invariant mass (sgn); M(K^{0}_{S}) [GeV/c^{2}]; DCA to Primary Vertex []; Entries",520,0.43,0.56,100,0.,10.0);
      TH2F *dcatoPVvspK0sBkg=new TH2F(nameBkg.Data(),"K^{0}_{S}: DCA to Primary Vertex vs  K^{0}_{S} invariant mass (bkg); M(K^{0}_{S}) [GeV/c^{2}]; DCA to Primary Vertex []; Entries",520,0.43,0.56,100,0.,10.0);

      TH2F* alldcatoPVvspK0sSgn= (TH2F*)dcatoPVvspK0sSgn->Clone();
      TH2F* alldcatoPVvspK0sBkg= (TH2F*)dcatoPVvspK0sBkg->Clone();
      TH2F* pidBachdcatoPVvspK0sSgn= (TH2F*)dcatoPVvspK0sSgn->Clone();
      TH2F* pidBachdcatoPVvspK0sBkg= (TH2F*)dcatoPVvspK0sBkg->Clone();

      fOutputAll->Add(alldcatoPVvspK0sSgn);
      fOutputPIDBach->Add(pidBachdcatoPVvspK0sSgn);
      fOutputAll->Add(alldcatoPVvspK0sBkg);
      fOutputPIDBach->Add(pidBachdcatoPVvspK0sBkg);

    }


    nameSgn="histK0SMassOfflineSgn";
    nameBkg="histK0SMassOfflineBkg";
    TH2F* spectrumK0SMassOfflineSgn = new TH2F(nameSgn.Data(), "K^{0}_{S} Signal invariant mass VS p_{T} - MC; M(K^{0}_{S}) [GeV/c^{2}]; p_{T}(K^{0}_{S}) [GeV/c]; Entries",
					       520,0.43,0.56,200,0.,20.);
    TH2F* spectrumK0SMassOfflineBkg = new TH2F(nameBkg.Data(), "K^{0}_{S} Background invariant mass VS p_{T} - MC; M(K^{0}_{S}) [GeV/c^{2}]; p_{T}(K^{0}_{S}) [GeV/c]; Entries",
					       520,0.43,0.56,200,0.,20.);
    spectrumK0SMassOfflineSgn->Sumw2();
    spectrumK0SMassOfflineBkg->Sumw2();
    spectrumK0SMassOfflineSgn->SetLineColor(2);
    spectrumK0SMassOfflineBkg->SetLineColor(4);
    spectrumK0SMassOfflineSgn->SetMarkerStyle(20);
    spectrumK0SMassOfflineBkg->SetMarkerStyle(20);
    spectrumK0SMassOfflineSgn->SetMarkerSize(0.6);
    spectrumK0SMassOfflineBkg->SetMarkerSize(0.6);
    spectrumK0SMassOfflineSgn->SetMarkerColor(2);
    spectrumK0SMassOfflineBkg->SetMarkerColor(4);

    TH2F* allspectrumK0SMassOfflineSgn = (TH2F*)spectrumK0SMassOfflineSgn->Clone(); 
    TH2F* allspectrumK0SMassOfflineBkg = (TH2F*) spectrumK0SMassOfflineBkg->Clone();  
    fOutputAll->Add(allspectrumK0SMassOfflineSgn);
    fOutputAll->Add(allspectrumK0SMassOfflineBkg);



    TH2F* pidBachspectrumK0SMassOfflineSgn = (TH2F*)spectrumK0SMassOfflineSgn->Clone(); 
    TH2F* pidBachspectrumK0SMassOfflineBkg = (TH2F*) spectrumK0SMassOfflineBkg->Clone();
    fOutputPIDBach->Add(pidBachspectrumK0SMassOfflineSgn);
    fOutputPIDBach->Add(pidBachspectrumK0SMassOfflineBkg);


    nameSgn="histLcMassOfflineByK0SSgn";
    nameBkg="histLcMassOfflineByK0SBkg";
    TH2F* spectrumLcMassOfflineByK0SSgn = new TH2F(nameSgn.Data(), "#Lambda_{C} Signal invariant mass (by K^{0}_{S})  vs p_{T} - MC; M(#Lambda_{C}) [GeV/c^{2}]; p_{T}",
						   1000,2.,2.5,200,0.,20.);
    TH2F* spectrumLcMassOfflineByK0SBkg = new TH2F(nameBkg.Data(), "#Lambda_{C} Background invariant mass (by K^{0}_{S})  vs p_{T} - MC; M(#Lambda_{C}) [GeV/c^{2}]; p_{T}",
						   1000,2.,2.5,200,0.,20.);
    spectrumLcMassOfflineByK0SSgn->Sumw2();
    spectrumLcMassOfflineByK0SBkg->Sumw2();
    spectrumLcMassOfflineByK0SSgn->SetLineColor(2);
    spectrumLcMassOfflineByK0SBkg->SetLineColor(4);
    spectrumLcMassOfflineByK0SSgn->SetMarkerStyle(20);
    spectrumLcMassOfflineByK0SBkg->SetMarkerStyle(20);
    spectrumLcMassOfflineByK0SSgn->SetMarkerSize(0.6);
    spectrumLcMassOfflineByK0SBkg->SetMarkerSize(0.6);
    spectrumLcMassOfflineByK0SSgn->SetMarkerColor(2);
    spectrumLcMassOfflineByK0SBkg->SetMarkerColor(4);


    TH2F* allspectrumLcMassOfflineByK0SSgn = (TH2F*)spectrumLcMassOfflineByK0SSgn->Clone(); 
    TH2F* allspectrumLcMassOfflineByK0SBkg = (TH2F*) spectrumLcMassOfflineByK0SBkg->Clone();  
    TH2F* pidBachspectrumLcMassOfflineByK0SSgn = (TH2F*)spectrumLcMassOfflineByK0SSgn->Clone(); 
    TH2F* pidBachspectrumLcMassOfflineByK0SBkg = (TH2F*) spectrumLcMassOfflineByK0SBkg->Clone();  
    fOutputAll->Add(allspectrumLcMassOfflineByK0SSgn);
    fOutputAll->Add(allspectrumLcMassOfflineByK0SBkg);
    fOutputPIDBach->Add(pidBachspectrumLcMassOfflineByK0SSgn);
    fOutputPIDBach->Add(pidBachspectrumLcMassOfflineByK0SBkg);
  
 


    nameSgn="histcosOAK0SpSgnvspOffline";
    nameBkg="histcosOAK0SpBkgvspOffline";
    TH2F* cosOpeningAngleK0SpSgnvspOffline= new TH2F(nameSgn.Data(),"#Lambda_{C} SGN : K^{0}_{S} - p Opening Angle  vs p  -  Offline  - MC ; Cos(Opening Angle);  p ",
						     200,-1.,1.,200,0.,20.);
    TH2F* cosOpeningAngleK0SpBkgvspOffline= new TH2F(nameBkg.Data(),"#Lambda_{C} BKG : K^{0}_{S} - p Opening Angle  vs p  -  Offline  - MC;  Cos(Opening Angle);  p ",
						     200,-1.,1.,200,0.,20.);
    

    TH2F* allcosOpeningAngleK0SpSgnvspOffline= (TH2F*)cosOpeningAngleK0SpSgnvspOffline->Clone(); 
    TH2F* allcosOpeningAngleK0SpBkgvspOffline= (TH2F*)cosOpeningAngleK0SpBkgvspOffline->Clone(); 
    TH2F* pidBachcosOpeningAngleK0SpSgnvspOffline= (TH2F*)cosOpeningAngleK0SpSgnvspOffline->Clone(); 
    TH2F* pidBachcosOpeningAngleK0SpBkgvspOffline= (TH2F*)cosOpeningAngleK0SpBkgvspOffline->Clone(); 
    fOutputAll->Add(allcosOpeningAngleK0SpSgnvspOffline); 
    fOutputAll->Add(allcosOpeningAngleK0SpBkgvspOffline); 
    fOutputPIDBach->Add(pidBachcosOpeningAngleK0SpSgnvspOffline); 
    fOutputPIDBach->Add(pidBachcosOpeningAngleK0SpBkgvspOffline); 




    nameSgn="histpK0SvspOfflineSgn";
    nameBkg="histpK0SvspOfflineBkg";
    TH2F* momentumDistributionK0SvspOfflineSgn= new TH2F(nameSgn.Data(),"#Lambda_{C} SGN : K^{0}_{S} vs p Total Momentum Distribution - Offline  - MC; p_{p};  p_{K^{0}_{S}}",
							 200,0.,20.,200,0.,20.);
    TH2F* momentumDistributionK0SvspOfflineBkg= new TH2F(nameBkg.Data(),"#Lambda_{C} BKG : K^{0}_{S} vs p Total Momentum Distribution - Offline  - MC;  p_{p}; p_{K^{0}_{S}}",
							 200,0.,20.,200,0.,20.);


    TH2F* allmomentumDistributionK0SvspOfflineSgn= (TH2F*)momentumDistributionK0SvspOfflineSgn->Clone(); 
    TH2F* allmomentumDistributionK0SvspOfflineBkg= (TH2F*)momentumDistributionK0SvspOfflineBkg->Clone(); 
    TH2F* pidBachmomentumDistributionK0SvspOfflineSgn= (TH2F*)momentumDistributionK0SvspOfflineSgn->Clone(); 
    TH2F* pidBachmomentumDistributionK0SvspOfflineBkg= (TH2F*)momentumDistributionK0SvspOfflineBkg->Clone(); 
    fOutputAll->Add(allmomentumDistributionK0SvspOfflineSgn); 
    fOutputAll->Add(allmomentumDistributionK0SvspOfflineBkg); 
    fOutputPIDBach->Add(pidBachmomentumDistributionK0SvspOfflineSgn); 
    fOutputPIDBach->Add(pidBachmomentumDistributionK0SvspOfflineBkg); 





    // armenteros-podolanski plots K0S (offline)
    nameSgn="histOfflineArmPodK0sSgn";
    nameBkg="histOfflineArmPodK0sBkg";
    TH2F* armenterosPodK0sOffSgn = new TH2F(nameSgn.Data(),"K^{0}_{S}  Armenteros-Podolanski distribution (sgn) -offline-; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]",200,-1.,1.,300,0.,0.3);
    TH2F* armenterosPodK0sOffBkg = new TH2F(nameBkg.Data(),"K^{0}_{S}  Armenteros-Podolanski distribution (bkg) -offline-; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]",200,-1.,1.,300,0.,0.3);
  

    TH2F* allArmenterosPodK0sOffSgn = (TH2F*)armenterosPodK0sOffSgn->Clone();
    TH2F* allArmenterosPodK0sOffBkg = (TH2F*)armenterosPodK0sOffBkg->Clone();
    TH2F* pidBachArmenterosPodK0sOffSgn = (TH2F*)armenterosPodK0sOffSgn->Clone();
    TH2F* pidBachArmenterosPodK0sOffBkg = (TH2F*)armenterosPodK0sOffBkg->Clone();


    fOutputAll->Add(allArmenterosPodK0sOffSgn);
    fOutputAll->Add(allArmenterosPodK0sOffBkg);
    fOutputPIDBach->Add(pidBachArmenterosPodK0sOffSgn);
    fOutputPIDBach->Add(pidBachArmenterosPodK0sOffBkg);


    nameSgn="histoDCAtoPVvsinvmassK0sOfflineSgn";
    nameBkg="histoDCAtoPVvsinvmassK0sOfflineBkg";
    TH2F *dcatoPVvspK0sOfflineSgn=new TH2F(nameSgn.Data(),"K^{0}_{S}: DCA to Primary Vertex vs  K^{0}_{S} invariant mass  -offline - (sgn); M(K^{0}_{S}) [GeV/c^{2}]; DCA to Primary Vertex []; Entries",520,0.43,0.56,100,0.,10.0);
    TH2F *dcatoPVvspK0sOfflineBkg=new TH2F(nameBkg.Data(),"K^{0}_{S}: DCA to Primary Vertex vs  K^{0}_{S} invariant mass  -offline - (bkg); M(K^{0}_{S}) [GeV/c^{2}]; DCA to Primary Vertex []; Entries",520,0.43,0.56,100,0.,10.0);
    

    TH2F* alldcatoPVvspK0sOfflineSgn= (TH2F*)dcatoPVvspK0sOfflineSgn->Clone();
    TH2F* pidBachdcatoPVvspK0sOfflineSgn= (TH2F*)dcatoPVvspK0sOfflineSgn->Clone();
    TH2F* alldcatoPVvspK0sOfflineBkg= (TH2F*)dcatoPVvspK0sOfflineBkg->Clone();
    TH2F* pidBachdcatoPVvspK0sOfflineBkg= (TH2F*)dcatoPVvspK0sOfflineBkg->Clone();



    fOutputAll->Add(alldcatoPVvspK0sOfflineSgn);
    fOutputPIDBach->Add(pidBachdcatoPVvspK0sOfflineSgn);
    fOutputAll->Add(alldcatoPVvspK0sOfflineBkg);
    fOutputPIDBach->Add(pidBachdcatoPVvspK0sOfflineBkg);

  }


  ///---------------- END  K0S HISTOS DECLARATIONS -------------------///
}
