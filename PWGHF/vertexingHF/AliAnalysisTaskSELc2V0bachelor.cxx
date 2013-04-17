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
#include <TH1I.h>
#include <TH2F.h>
#include <TTree.h>
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
  fIsEventSelected(kFALSE),
  fWriteVariableTree(kFALSE),
  fVariablesTree(0),
  fCandidateVariables(),
  fVtx1(0),
  fBzkG(0)
{
  //
  // Default ctor
  //
}
//___________________________________________________________________________
AliAnalysisTaskSELc2V0bachelor::AliAnalysisTaskSELc2V0bachelor(const Char_t* name, AliRDHFCutsLctoV0* prodCuts,
							       AliRDHFCutsLctoV0* analCuts, Bool_t useOnTheFly,
							       Bool_t writeVariableTree) :
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
  fIsEventSelected(kFALSE),
  fWriteVariableTree(writeVariableTree),
  fVariablesTree(0),
  fCandidateVariables(),
  fVtx1(0),
  fBzkG(0)
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

  // Output slot #6 keeps a tree of the candidate variables after track selection
  if (fWriteVariableTree) DefineOutput(6,TTree::Class());  //My private output

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

  if(fVariablesTree){
    delete fVariablesTree;
    fVariablesTree = 0;
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
  fVtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();
  if (!fVtx1) return;

  // fix for temporary bug in ESDfilter 
  fBzkG = (Double_t)aodEvent->GetMagneticField(); 
  if (TMath::Abs(fBzkG)<0.001) return;
  fCEvents->Fill(2);

  Float_t zVertex = fVtx1->GetZ();
  ((TH1F*)(fOutput->FindObject("hZ2")))->Fill(zVertex);

  if (fVtx1->GetNContributors()<1) return;
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

  TString titleVtx=fVtx1->GetTitle();
  if (titleVtx.Contains("Z")) {
    fCEvents->Fill(11);
    ((TH1F*)(fOutput->FindObject("hZ11")))->Fill(zVertex);
  }
  else if (titleVtx.Contains("3D")) {
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
  }
  ///////////////////////

  if ( !fIsEventSelected ) return; // don't take into account not selected events 

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

    // check on MC Lc Daughter
    for (Int_t iii=0; iii<mcArray->GetEntries(); iii++) {
      SearchLcDaughter(mcArray,iii);
    }

  }

  Int_t nSelectedProd = 0;
  Int_t nSelectedAnal = 0;
  if (fIsK0sAnalysis) {
    MakeAnalysisForLc2prK0S(arrayLctopKos,mcArray,
			    nSelectedProd, fProdCuts, nSelectedAnal, fAnalCuts);

    if (nSelectedAnal) {

      ((TH1I*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(4);
      ((TH1F*)(fOutput->FindObject("hZ4a")))->Fill(zVertex);

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

      if (titleVtx.Contains("Z")) {
	((TH1I*)(fOutput->FindObject("hEventsWithCandidates")))->Fill(11);
	((TH1F*)(fOutput->FindObject("hZ11a")))->Fill(zVertex);
      }
      else if (titleVtx.Contains("3D")) {
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
  if (fWriteVariableTree) PostData(6,fVariablesTree);

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
  
  //fCEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fCEvents"));

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

  // Output slot 6: tree of the candidate variables
  if (fWriteVariableTree) {
    const char* nameoutput = GetOutputSlot(6)->GetContainer()->GetName();
    fVariablesTree = new TTree(nameoutput,"Candidates variables tree");
    Int_t nVar = 80;
    fCandidateVariables = new Float_t [nVar];
    TString * fCandidateVariableNames = new TString[nVar];
    fCandidateVariableNames[0]="isLcByMC";
    fCandidateVariableNames[1]="isV0ByMC";
    fCandidateVariableNames[2]="massLc2K0Sp";
    fCandidateVariableNames[3]="massLc2Lambdapi";
    fCandidateVariableNames[4]="massK0S";
    fCandidateVariableNames[5]="massLambda";
    fCandidateVariableNames[6]="massLambdaBar";
    fCandidateVariableNames[7]="cosPAK0S";
    fCandidateVariableNames[8]="dcaV0ptp";
    fCandidateVariableNames[9]="tImpParBach";
    fCandidateVariableNames[10]="tImpParV0";
    fCandidateVariableNames[11]="nSigmaTPCpr";
    fCandidateVariableNames[12]="nSigmaTPCpi";
    fCandidateVariableNames[13]="nSigmaTPCka";
    fCandidateVariableNames[14]="nSigmaTOFpr";
    fCandidateVariableNames[15]="nSigmaTOFpi";
    fCandidateVariableNames[16]="nSigmaTOFka";
    fCandidateVariableNames[17]="bachelorPx";
    fCandidateVariableNames[18]="bachelorPy";
    fCandidateVariableNames[19]="bachelorPz";
    fCandidateVariableNames[20]="V0positivePx";
    fCandidateVariableNames[21]="V0positivePy";
    fCandidateVariableNames[22]="V0positivePz";
    fCandidateVariableNames[23]="V0negativePx";
    fCandidateVariableNames[24]="V0negativePy";
    fCandidateVariableNames[25]="V0negativePz";
    fCandidateVariableNames[26]="qtLc";
    fCandidateVariableNames[27]="alphaLc";
    fCandidateVariableNames[28]="dcaV0postoPV";
    fCandidateVariableNames[29]="dcaV0negtoPV";
    fCandidateVariableNames[30]="bachelorPxDCA";
    fCandidateVariableNames[31]="bachelorPyDCA";
    fCandidateVariableNames[32]="bachelorPzDCA";
    fCandidateVariableNames[33]="v0PxDCA";
    fCandidateVariableNames[34]="v0PyDCA";
    fCandidateVariableNames[35]="v0PzDCA";
    fCandidateVariableNames[36]="V0positivePxDCA";
    fCandidateVariableNames[37]="V0positivePyDCA";
    fCandidateVariableNames[38]="V0positivePzDCA";
    fCandidateVariableNames[39]="V0negativePxDCA";
    fCandidateVariableNames[40]="V0negativePyDCA";
    fCandidateVariableNames[41]="V0negativePzDCA";
    fCandidateVariableNames[42]="flagToCheckCandidate";
    fCandidateVariableNames[43]="massGamma";

    fCandidateVariableNames[44]="bachelorP"; // @ prim vtx
    fCandidateVariableNames[45]="bachelorPt"; // @ prim vtx
    fCandidateVariableNames[46]="V0positiveP"; // @ prim vtx
    fCandidateVariableNames[47]="V0positivePt"; // @ prim vtx
    fCandidateVariableNames[48]="V0negativeP"; // @ prim vtx
    fCandidateVariableNames[49]="V0negativePt"; // @ prim vtx
    fCandidateVariableNames[50]="bachelorPDCA"; // @ DCA
    fCandidateVariableNames[51]="bachelorPtDCA"; // @ DCA
    fCandidateVariableNames[52]="v0PDCA"; // @ DCA
    fCandidateVariableNames[53]="v0PtDCA"; // @ DCA
    fCandidateVariableNames[54]="V0positivePDCA"; // @ DCA
    fCandidateVariableNames[55]="V0positivePtDCA"; // @ DCA
    fCandidateVariableNames[56]="V0negativePDCA"; // @ DCA
    fCandidateVariableNames[57]="V0negativePtDCA"; // @ DCA
    fCandidateVariableNames[58]="LcP"; // @ DCA
    fCandidateVariableNames[59]="LcPt"; // @ DCA
    fCandidateVariableNames[60]="v0P"; // @ V0 DCA
    fCandidateVariableNames[61]="v0Pt"; // @ V0 DCA

    fCandidateVariableNames[62]="cosPALc";
    fCandidateVariableNames[63]="decayLengthLc";
    fCandidateVariableNames[64]="decayLengthV0";

    fCandidateVariableNames[65]="yLc";

    fCandidateVariableNames[66]="massD2K0Spi"; // D+ -> pi+ K0S
    fCandidateVariableNames[67]="massDS2K0SK"; // D+S -> K+ K0S

    fCandidateVariableNames[68]="nSigmaITSpi"; // nSigmaITSpi
    fCandidateVariableNames[69]="nSigmaITSka"; // nSigmaITSka
    fCandidateVariableNames[70]="nSigmaITSpr"; // nSigmaITSpr

    fCandidateVariableNames[71]="dcaLcptp"; // DCA Lc prong-to-prong

    fCandidateVariableNames[72]="cosPAV0XY"; // cosPA XY x V0
    fCandidateVariableNames[73]="cosPALcXY"; // cosPA XY x V0

    fCandidateVariableNames[74]="decayLengthV0XY"; // decay length XY x V0
    fCandidateVariableNames[75]="decayLengthLcXY"; // decay length XY x V0

    fCandidateVariableNames[76]="normalizedDecayLengthV0"; // normalized decay length x V0
    fCandidateVariableNames[77]="normalizedDecayLengthLc"; // normalized decay length x Lc

    fCandidateVariableNames[78]="normalizedDecayLengthXYV0"; // normalized decay length XY x V0
    fCandidateVariableNames[79]="normalizedDecayLengthXYLc"; // normalized decay length XY x Lc

    for(Int_t ivar=0; ivar<nVar; ivar++){
      fVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
    }
    PostData(6,fVariablesTree);
  }

  return;
}
//___________________________________ hiostograms _______________________________________
void  AliAnalysisTaskSELc2V0bachelor::DefineHistograms() {

  fCEvents = new TH1F("fCEvents","conter",17,0,17);
  fCEvents->SetStats(kTRUE);
  fCEvents->GetXaxis()->SetBinLabel(1,"X1");
  fCEvents->GetXaxis()->SetBinLabel(2,"Analyzed events");
  fCEvents->GetXaxis()->SetBinLabel(3,"AliAODVertex exists");
  fCEvents->GetXaxis()->SetBinLabel(4,"GetNContributors()>0");
  fCEvents->GetXaxis()->SetBinLabel(5,"MCarray exists");
  fCEvents->GetXaxis()->SetBinLabel(6,"CascadesHF exists");
  fCEvents->GetXaxis()->SetBinLabel(7,"MCheader exists");
  fCEvents->GetXaxis()->SetBinLabel(8,"IsEventSelected");
  fCEvents->GetXaxis()->SetBinLabel(9,"triggerClass!=CINT1");
  fCEvents->GetXaxis()->SetBinLabel(10,"triggerMask!=kAnyINT");
  fCEvents->GetXaxis()->SetBinLabel(11,"triggerMask!=kAny");
  fCEvents->GetXaxis()->SetBinLabel(12,"vtxTitle.Contains(Z)");
  fCEvents->GetXaxis()->SetBinLabel(13,"vtxTitle.Contains(3D)");
  fCEvents->GetXaxis()->SetBinLabel(14,"vtxTitle.Doesn'tContain(Z-3D)");
  fCEvents->GetXaxis()->SetBinLabel(15,Form("zVtx<=%2.0fcm",fAnalCuts->GetMaxVtxZ()));
  fCEvents->GetXaxis()->SetBinLabel(16,"!IsEventSelected");
  fCEvents->GetXaxis()->SetBinLabel(17,"triggerMask!=kAnyINT || triggerClass!=CINT1");
  //fCEvents->GetXaxis()->SetTitle("");
  fCEvents->GetYaxis()->SetTitle("counts");

  fOutput->Add(fCEvents);
  TString fillthis="";

  if (fUseMCInfo) {
    fillthis="histMcStatLc";
    TH1F* mcStatisticLc = new TH1F(fillthis.Data(),"#Lambda_{C} generated and their decays",19,-9.5,9.5);
    fOutput->Add(mcStatisticLc);
  }

  //fillthis="histopionV0SigmaVspTOF";
  //TH2F *hpionV0SigmaVspTOF=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);
  fillthis="histoprotonBachSigmaVspTOF";
  TH2F *hprotonBachSigmaVspTOF=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);

  //fOutput->Add(hpionV0SigmaVspTOF);
  fOutput->Add(hprotonBachSigmaVspTOF);

  //fillthis="histopionV0SigmaVspTPC";
  //TH2F *hpionV0SigmaVspTPC=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);
  fillthis="histoprotonBachSigmaVspTPC";
  TH2F *hprotonBachSigmaVspTPC=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);

  //fOutput->Add(hpionV0SigmaVspTPC);
  fOutput->Add(hprotonBachSigmaVspTPC);

  if (fUseMCInfo) {

    //fillthis="histopionV0SigmaVspTOFsgn";
    //TH2F *hpionV0SigmaVspTOFsgn=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);
    fillthis="histoprotonBachSigmaVspTOFsgn";
    TH2F *hprotonBachSigmaVspTOFsgn=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);

    //fOutput->Add(hpionV0SigmaVspTOFsgn);
    fOutput->Add(hprotonBachSigmaVspTOFsgn);

    //fillthis="histopionV0SigmaVspTPCsgn";
    //TH2F *hpionV0SigmaVspTPCsgn=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);
    fillthis="histoprotonBachSigmaVspTPCsgn";
    TH2F *hprotonBachSigmaVspTPCsgn=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);

    //fOutput->Add(hpionV0SigmaVspTPCsgn);
    fOutput->Add(hprotonBachSigmaVspTPCsgn);


    //fillthis="histopionV0SigmaVspTOFbkg";
    //TH2F *hpionV0SigmaVspTOFbkg=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);
    fillthis="histoprotonBachSigmaVspTOFbkg";
    TH2F *hprotonBachSigmaVspTOFbkg=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);

    //fOutput->Add(hpionV0SigmaVspTOFbkg);
    fOutput->Add(hprotonBachSigmaVspTOFbkg);

    //fillthis="histopionV0SigmaVspTPCbkg";
    //TH2F *hpionV0SigmaVspTPCbkg=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);
    fillthis="histoprotonBachSigmaVspTPCbkg";
    TH2F *hprotonBachSigmaVspTPCbkg=new TH2F(fillthis.Data(),fillthis.Data(),300,0.,30.,100,-5.,5.);

    //fOutput->Add(hpionV0SigmaVspTPCbkg);
    fOutput->Add(hprotonBachSigmaVspTPCbkg);


  }


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

  TH1I *hCandidateSelection = new TH1I("hCandidateSelection","",14,0,14);
  hCandidateSelection->GetXaxis()->SetBinLabel(1,"IsEventSelected");
  hCandidateSelection->GetXaxis()->SetBinLabel(2,"IsSecondaryVtx");
  hCandidateSelection->GetXaxis()->SetBinLabel(3,"V0toPosNeg");
  hCandidateSelection->GetXaxis()->SetBinLabel(4,"offlineV0");
  hCandidateSelection->GetXaxis()->SetBinLabel(5,"prodCuts::kTracks");
  hCandidateSelection->GetXaxis()->SetBinLabel(6,"prodCuts::kCandidate");
  hCandidateSelection->GetXaxis()->SetBinLabel(7,"prodCuts::kPID");
  hCandidateSelection->GetXaxis()->SetBinLabel(8,"prodCuts::kAll");
  hCandidateSelection->GetXaxis()->SetBinLabel(9,"offlineV0");
  hCandidateSelection->GetXaxis()->SetBinLabel(10,"isInFiducialAcceptance");
  hCandidateSelection->GetXaxis()->SetBinLabel(11,"analCuts::kTracks");
  hCandidateSelection->GetXaxis()->SetBinLabel(12,"analCuts::kCandidate");
  hCandidateSelection->GetXaxis()->SetBinLabel(13,"analCuts::kPID");
  hCandidateSelection->GetXaxis()->SetBinLabel(14,"analCuts::kAll");
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

  TH1I *hSwitchOnCandidates1 = new TH1I("hSwitchOnCandidates1","",16,-8,8);
  fOutput->Add(hSwitchOnCandidates1);
  TH1I *hSwitchOnCandidates2 = new TH1I("hSwitchOnCandidates2","",16,-8,8);
  fOutput->Add(hSwitchOnCandidates2);
  TH1I *hSwitchOnCandidates3 = new TH1I("hSwitchOnCandidates3","",16,-8,8);
  fOutput->Add(hSwitchOnCandidates3);
  TH1I *hSwitchOnCandidates4 = new TH1I("hSwitchOnCandidates4","",16,-8,8);
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
							 AliRDHFCutsLctoV0 *cutsAnal,
							 TClonesArray *mcArray)
{
  //
  // Fill histos for Lc -> K0S+proton
  //

  TString fillthis="";

  Double_t invmassLc = part->InvMassLctoK0sP();
  Double_t lambdacpt = part->Pt();

  AliAODv0 * v0part = (AliAODv0*)part->Getv0();
  Bool_t onFlyV0 = v0part->GetOnFlyStatus(); // on-the-flight V0s
  Double_t momK0s = TMath::Sqrt(v0part->Ptot2V0());
  Double_t ptK0s = TMath::Sqrt(v0part->Pt2V0());
  Double_t dcaV0ptp = v0part->GetDCA();
  Double_t invmassK0s = v0part->MassK0Short();
  Bool_t isInV0windowProd = (((cutsProd->IsSelectedSingleCut(part,AliRDHFCuts::kCandidate,2))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)); // cut on V0 invMass
  Bool_t isInCascadeWindow = (((cutsAnal->IsSelectedSingleCut(part,AliRDHFCuts::kCandidate,0))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)); // cut on Lc->p+K0S invMass
  Bool_t isCandidateSelectedCuts = (((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)); // kinematic/topological cuts

  AliAODTrack *bachelor = (AliAODTrack*)part->GetBachelor();
  Double_t momBach  = bachelor->P();
  Bool_t isBachelorID = (((cutsAnal->IsSelected(part,AliRDHFCuts::kPID))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)); // ID x bachelor

  /*
  if (fIsEventSelected) {
    if ( ( !onFlyV0 ||
	   (onFlyV0 && fUseOnTheFlyV0) ) &&
	 isCandidateSelectedCuts && isBachelorID) {
      fillthis="hist1LcMassOfflineByK0S";
      if (isBachelorID) ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);
    }
  } else {
    if ( ( !onFlyV0 ||
	   (onFlyV0 && fUseOnTheFlyV0) ) &&
	 isCandidateSelectedCuts && isBachelorID) {
      fillthis="hist0LcMassOfflineByK0S";
      if (isBachelorID) ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);
    }
    return; // don't take into account not selected events
  }
  */

  if ( !onFlyV0 )
    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(3); // it counts number of candidates coming from offline V0s

  if ( !( !onFlyV0 || (onFlyV0 && fUseOnTheFlyV0) ) ) return;

  if ( !( ( (cutsProd->IsSelected(part,AliRDHFCuts::kTracks))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) ) return;
  ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(4);
  if ( ( ( (cutsProd->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) )
    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(5);
  if ( ( ( (cutsProd->IsSelected(part,AliRDHFCuts::kPID))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) )
    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(6);
  if ( ( ( (cutsProd->IsSelected(part,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) ) {
    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(7);
    nSelectedProd++;
  }


  if ( !onFlyV0 )
    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(8); // it counts number of candidates coming from offline V0s

  if ( cutsAnal->IsInFiducialAcceptance(part->Pt(),part->Y(4122)) )
    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(9);
  if ( (((cutsAnal->IsSelected(part,AliRDHFCuts::kTracks))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) )
    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(10);
  if ( (((cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) )
    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(11);
  if ( (((cutsAnal->IsSelected(part,AliRDHFCuts::kPID))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) )
    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(12);
  if ( (((cutsAnal->IsSelected(part,AliRDHFCuts::kAll))&(AliRDHFCutsLctoV0::kLcToK0Spr))==(AliRDHFCutsLctoV0::kLcToK0Spr)) ) {
    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(13);
    nSelectedAnal++;
  }

  if ( !(cutsAnal->IsInFiducialAcceptance(part->Pt(),part->Y(4122))) ) return;

  if ( !( ( (cutsAnal->IsSelected(part,AliRDHFCuts::kTracks))&(AliRDHFCutsLctoV0::kLcToK0Spr)) == (AliRDHFCutsLctoV0::kLcToK0Spr) ) ) return;

  Int_t aaa = cutsAnal->IsSelected(part,AliRDHFCuts::kTracks);
  if ( (aaa&AliRDHFCutsLctoV0::kLcToK0Spr)==AliRDHFCutsLctoV0::kLcToK0Spr ) {
    if ( ( (aaa&AliRDHFCutsLctoV0::kLcToLpi)==AliRDHFCutsLctoV0::kLcToLpi && bachelor->Charge()==-1)  ||
	 ( (aaa&AliRDHFCutsLctoV0::kLcToLBarpi)==AliRDHFCutsLctoV0::kLcToLBarpi && bachelor->Charge()==+1) )
      ((TH1I*)(fOutput->FindObject("hSwitchOnCandidates1")))->Fill( -aaa );
    else
      ((TH1I*)(fOutput->FindObject("hSwitchOnCandidates1")))->Fill( aaa );
  }

  aaa = cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate);
  if ((aaa&AliRDHFCutsLctoV0::kLcToK0Spr)==AliRDHFCutsLctoV0::kLcToK0Spr) {
    if ( ( (aaa&AliRDHFCutsLctoV0::kLcToLpi)==AliRDHFCutsLctoV0::kLcToLpi && bachelor->Charge()==-1) ||
	 ( (aaa&AliRDHFCutsLctoV0::kLcToLBarpi)==AliRDHFCutsLctoV0::kLcToLBarpi && bachelor->Charge()==+1) )
      ((TH1I*)(fOutput->FindObject("hSwitchOnCandidates2")))->Fill( -aaa );
    else
      ((TH1I*)(fOutput->FindObject("hSwitchOnCandidates2")))->Fill( aaa );
  }

  aaa = cutsAnal->IsSelected(part,AliRDHFCuts::kPID);
  if ((aaa&AliRDHFCutsLctoV0::kLcToK0Spr)==AliRDHFCutsLctoV0::kLcToK0Spr) {
    if ( ( (aaa&AliRDHFCutsLctoV0::kLcToLpi)==AliRDHFCutsLctoV0::kLcToLpi && bachelor->Charge()==-1) ||
	 ( (aaa&AliRDHFCutsLctoV0::kLcToLBarpi)==AliRDHFCutsLctoV0::kLcToLBarpi && bachelor->Charge()==+1) )
      ((TH1I*)(fOutput->FindObject("hSwitchOnCandidates3")))->Fill( -aaa );
    else
      ((TH1I*)(fOutput->FindObject("hSwitchOnCandidates3")))->Fill( aaa );
  }

  aaa = cutsAnal->IsSelected(part,AliRDHFCuts::kAll);
  if ((aaa&AliRDHFCutsLctoV0::kLcToK0Spr)==AliRDHFCutsLctoV0::kLcToK0Spr) {
    if ( ( (aaa&AliRDHFCutsLctoV0::kLcToLpi)==AliRDHFCutsLctoV0::kLcToLpi && bachelor->Charge()==-1) ||
	 ( (aaa&AliRDHFCutsLctoV0::kLcToLBarpi)==AliRDHFCutsLctoV0::kLcToLBarpi && bachelor->Charge()==+1) )
      ((TH1I*)(fOutput->FindObject("hSwitchOnCandidates4")))->Fill( -aaa );
    else
      ((TH1I*)(fOutput->FindObject("hSwitchOnCandidates4")))->Fill( aaa );
  }




  aaa = cutsAnal->IsSelected(part,AliRDHFCuts::kCandidate);
  Int_t flagToCheckCandidate = 0;
  if ( (aaa&AliRDHFCutsLctoV0::kLcToK0Spr)==AliRDHFCutsLctoV0::kLcToK0Spr ) {
    if ( aaa==AliRDHFCutsLctoV0::kLcToK0Spr ) {
      flagToCheckCandidate = aaa; // Lc->K0S+p OK
    } else {
      if ( (aaa&AliRDHFCutsLctoV0::kLcToLpi)==AliRDHFCutsLctoV0::kLcToLpi ) {
	if (bachelor->Charge()==+1)
	  flagToCheckCandidate = aaa; // Lc->Lambda+pi+
	else if (bachelor->Charge()==-1)
	  flagToCheckCandidate =-aaa;//+(AliRDHFCutsLctoV0::kLcToK0Spr); // Lambda+pi- AS Lc->K0S+p candidate
      }
      if ( (aaa&AliRDHFCutsLctoV0::kLcToLBarpi)==AliRDHFCutsLctoV0::kLcToLBarpi ) {
	if (bachelor->Charge()==-1)
	  flagToCheckCandidate = aaa; // Lc->LambdaBar+pi-
	else if (bachelor->Charge()==+1)
	  flagToCheckCandidate =-aaa;//+(AliRDHFCutsLctoV0::kLcToK0Spr); // LambdaBar+pi+ AS Lc->K0S+p candidate
      }
    }
  } else {
    if ( aaa==AliRDHFCutsLctoV0::kLcToK0Spr ) {
      flagToCheckCandidate = -10-(AliRDHFCutsLctoV0::kLcToK0Spr); // NEVER
    } else {
      if ( (aaa&AliRDHFCutsLctoV0::kLcToLpi)==AliRDHFCutsLctoV0::kLcToLpi ) {
	if (bachelor->Charge()==+1)
	  flagToCheckCandidate = aaa; // Lc->Lambda+pi+ OK
	else if (bachelor->Charge()==-1)
	  flagToCheckCandidate =-aaa;//+(AliRDHFCutsLctoV0::kLcToK0Spr); // Lambda+pi- AS Lc->Lambda+pi+ candidate
      }
      if ( (aaa&AliRDHFCutsLctoV0::kLcToLBarpi)==AliRDHFCutsLctoV0::kLcToLBarpi ) {
	if (bachelor->Charge()==-1)
	  flagToCheckCandidate = aaa; // Lc->LambdaBar+pi- OK
	else if (bachelor->Charge()==+1)
	  flagToCheckCandidate =-aaa;//+(AliRDHFCutsLctoV0::kLcToK0Spr); // LambdaBar+pi+ AS Lc->LambdaBar+pi- candidate
      }
    }
  }



  Int_t pdgCand = 4122;
  Int_t pdgDgLctoV0bachelor[2]={3122,211};
  Int_t pdgDgV0toDaughters[2]={2212,211};
  Int_t isLc2LBarpi=0, isLc2Lpi=0;
  Int_t mcLabel = 0;
  if (fUseMCInfo) {
    mcLabel = part->MatchToMC(pdgCand,pdgDgLctoV0bachelor[0],pdgDgLctoV0bachelor,pdgDgV0toDaughters,mcArray,kTRUE);
    if (mcLabel>=0) {
      if (bachelor->Charge()==-1) isLc2LBarpi=1;
      if (bachelor->Charge()==+1) isLc2Lpi=1;
    }
  }

  Int_t pdgDg2prong[2] = {211, 211};
  Int_t labelK0S = 0;
  Int_t isK0S = 0;
  if (fUseMCInfo) {
    labelK0S = v0part->MatchToMC(310,mcArray,2,pdgDg2prong);
    if (labelK0S>=0) isK0S = 1;
  }

  pdgDg2prong[0] = 211;
  pdgDg2prong[1] = 2212;
  Int_t isLambda = 0;
  Int_t isLambdaBar = 0;
  Int_t lambdaLabel = 0;
  if (fUseMCInfo) {
    lambdaLabel = v0part->MatchToMC(3122,mcArray,2,pdgDg2prong);
    if (lambdaLabel>=0) {
      AliAODMCParticle *lambdaTrack = (AliAODMCParticle*)mcArray->At(lambdaLabel);
      if (lambdaTrack->GetPdgCode()==3122) isLambda = 1;
      else if (lambdaTrack->GetPdgCode()==-3122) isLambdaBar = 1;
    }
  }

  pdgDg2prong[0] = 11;
  pdgDg2prong[1] = 11;
  Int_t isGamma = 0;
  Int_t gammaLabel = 0;
  if (fUseMCInfo) {
    gammaLabel = v0part->MatchToMC(22,mcArray,2,pdgDg2prong);
    if (gammaLabel>=0) {
      AliAODMCParticle *gammaTrack = (AliAODMCParticle*)mcArray->At(gammaLabel);
      if (gammaTrack->GetPdgCode()==22) isGamma = 1;
    }
  }

  Double_t invmassLc2Lpi = part->InvMassLctoLambdaPi();
  Double_t invmassLambda = v0part->MassLambda();
  Double_t invmassLambdaBar = v0part->MassAntiLambda();

  Int_t isLcByMC = isLc+isLc2LBarpi*2+isLc2Lpi*4;
  Int_t isV0ByMC = isK0S+isLambdaBar*2+isLambda*4+isGamma*8;

  Double_t nSigmaITSpr=-999.;
  cutsAnal->GetPidHF()->GetnSigmaITS(bachelor,4,nSigmaITSpr);
  Double_t nSigmaTPCpr=-999.;
  cutsAnal->GetPidHF()->GetnSigmaTPC(bachelor,4,nSigmaTPCpr);
  Double_t nSigmaTOFpr=-999.;
  cutsAnal->GetPidHF()->GetnSigmaTOF(bachelor,4,nSigmaTOFpr);

  Double_t nSigmaITSpi=-999.;
  cutsAnal->GetPidHF()->GetnSigmaITS(bachelor,2,nSigmaITSpi);
  Double_t nSigmaTPCpi=-999.;
  cutsAnal->GetPidHF()->GetnSigmaTPC(bachelor,2,nSigmaTPCpi);
  Double_t nSigmaTOFpi=-999.;
  cutsAnal->GetPidHF()->GetnSigmaTOF(bachelor,2,nSigmaTOFpi);

  Double_t nSigmaITSka=-999.;
  cutsAnal->GetPidHF()->GetnSigmaITS(bachelor,3,nSigmaITSka);
  Double_t nSigmaTPCka=-999.;
  cutsAnal->GetPidHF()->GetnSigmaTPC(bachelor,3,nSigmaTPCka);
  Double_t nSigmaTOFka=-999.;
  cutsAnal->GetPidHF()->GetnSigmaTOF(bachelor,3,nSigmaTOFka);


  // Fill candidate variable Tree (track selection, V0 invMass selection)
  if ( fWriteVariableTree && !onFlyV0 && isInV0windowProd && isInCascadeWindow && part->CosV0PointingAngle()>0.99) {

    fCandidateVariables[0] = fUseMCInfo+isLcByMC; // 0: real data; 1: bkg; 2: Lc->K0S+p; 3: Lc->LambdaBar+pbar; 5: Lc->Lambda+p
    fCandidateVariables[1] = fUseMCInfo+isV0ByMC; // 0: real data; 1: bkg; 2: K0S->pi+pi; 3: LambdaBar->pbar+pi+; 5: Lambda->p+pi-
    fCandidateVariables[2] = invmassLc;
    fCandidateVariables[3] = invmassLc2Lpi;
    fCandidateVariables[4] = invmassK0s;
    fCandidateVariables[5] = invmassLambda;
    fCandidateVariables[6] = invmassLambdaBar;
    fCandidateVariables[7] = part->CosV0PointingAngle();
    fCandidateVariables[8] = dcaV0ptp;
    fCandidateVariables[9] = part->Getd0Prong(0);
    fCandidateVariables[10] = part->Getd0Prong(1);
    fCandidateVariables[11] = nSigmaTPCpr;
    fCandidateVariables[12] = nSigmaTPCpi;
    fCandidateVariables[13] = nSigmaTPCka;
    fCandidateVariables[14] = nSigmaTOFpr;
    fCandidateVariables[15] = nSigmaTOFpi;
    fCandidateVariables[16] = nSigmaTOFka;
    fCandidateVariables[17] = bachelor->Px();
    fCandidateVariables[18] = bachelor->Py();
    fCandidateVariables[19] = bachelor->Pz();
    AliAODTrack *v0neg = (AliAODTrack*)part->Getv0NegativeTrack(); 
    fCandidateVariables[20] = v0neg->Px();
    fCandidateVariables[21] = v0neg->Py();
    fCandidateVariables[22] = v0neg->Pz();
    AliAODTrack *v0pos = (AliAODTrack*)part->Getv0PositiveTrack();
    fCandidateVariables[23] = v0pos->Px();
    fCandidateVariables[24] = v0pos->Py();
    fCandidateVariables[25] = v0pos->Pz();
    fCandidateVariables[26] = part->QtProng(0);
    fCandidateVariables[27] = part->Alpha();
    fCandidateVariables[28] = v0part->Getd0Prong(0);
    fCandidateVariables[29] = v0part->Getd0Prong(1);
    fCandidateVariables[30] = part->PxProng(0);
    fCandidateVariables[31] = part->PyProng(0);
    fCandidateVariables[32] = part->PzProng(0);
    fCandidateVariables[33] = part->PxProng(1);
    fCandidateVariables[34] = part->PyProng(1);
    fCandidateVariables[35] = part->PzProng(1);
    fCandidateVariables[36] = v0part->PxProng(0);
    fCandidateVariables[37] = v0part->PyProng(0);
    fCandidateVariables[38] = v0part->PzProng(0);
    fCandidateVariables[39] = v0part->PxProng(1);
    fCandidateVariables[40] = v0part->PyProng(1);
    fCandidateVariables[41] = v0part->PzProng(1);
    fCandidateVariables[42] = flagToCheckCandidate;
    fCandidateVariables[43] = v0part->InvMass2Prongs(0,1,11,11);

    fCandidateVariables[44] = bachelor->P();
    fCandidateVariables[45] = bachelor->Pt();
    fCandidateVariables[46] = v0pos->P();
    fCandidateVariables[47] = v0pos->Pt();
    fCandidateVariables[48] = v0neg->P();
    fCandidateVariables[49] = v0neg->Pt();
    fCandidateVariables[50] = part->PProng(0);
    fCandidateVariables[51] = part->PtProng(0);
    fCandidateVariables[52] = part->PProng(1);
    fCandidateVariables[53] = part->PtProng(1);
    fCandidateVariables[54] = v0part->PProng(0);
    fCandidateVariables[55] = v0part->PtProng(0);
    fCandidateVariables[56] = v0part->PProng(1);
    fCandidateVariables[57] = v0part->PtProng(1);
    fCandidateVariables[58] = part->P();
    fCandidateVariables[59] = part->Pt();
    fCandidateVariables[60] = v0part->P();
    fCandidateVariables[61] = v0part->Pt();

    fCandidateVariables[62] = part->CosPointingAngle();
    fCandidateVariables[63] = part->DecayLength();
    fCandidateVariables[64] = part->DecayLengthV0();

    fCandidateVariables[65] = part->Y(4122);

    fCandidateVariables[66] = part->InvMass2Prongs(0,1,211,310); // D+ -> pi+ K0S
    fCandidateVariables[67] = part->InvMass2Prongs(0,1,321,310); // D+S -> K+ K0S

    fCandidateVariables[68] = nSigmaITSpr;
    fCandidateVariables[69] = nSigmaITSpi;
    fCandidateVariables[70] = nSigmaITSka;

    fCandidateVariables[71] = part->GetDCA();

    fCandidateVariables[72] = part->CosV0PointingAngleXY();
    fCandidateVariables[73] = part->CosPointingAngleXY();

    fCandidateVariables[74] = part->DecayLengthXYV0();
    fCandidateVariables[75] = part->DecayLengthXY();

    fCandidateVariables[76] = part->NormalizedV0DecayLength();
    fCandidateVariables[77] = part->NormalizedDecayLength();

    fCandidateVariables[78] = part->NormalizedV0DecayLengthXY();
    fCandidateVariables[79] = part->NormalizedDecayLengthXY();

    Double_t v0Momentum = (v0part->PxProng(0)+v0part->PxProng(1))*(v0part->PxProng(0)+v0part->PxProng(1));
    v0Momentum += (v0part->PyProng(0)+v0part->PyProng(1))*(v0part->PyProng(0)+v0part->PyProng(1));
    v0Momentum += (v0part->PzProng(0)+v0part->PzProng(1))*(v0part->PzProng(0)+v0part->PzProng(1));
    v0Momentum = TMath::Sqrt(v0Momentum);

    Double_t lcMomentum = (part->PxProng(0)+part->PxProng(1))*(part->PxProng(0)+part->PxProng(1));
    lcMomentum += (part->PyProng(0)+part->PyProng(1))*(part->PyProng(0)+part->PyProng(1));
    lcMomentum += (part->PzProng(0)+part->PzProng(1))*(part->PzProng(0)+part->PzProng(1));
    lcMomentum = TMath::Sqrt(lcMomentum);

    fVariablesTree->Fill();
  }



  if (onFlyV0 && fUseOnTheFlyV0) {  

    fillthis="histK0SMass";
    ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassK0s,ptK0s);
    if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassK0s,ptK0s);

    if (isCandidateSelectedCuts) {

      fillthis="histpK0Svsp";
      ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momBach,momK0s);
      if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momBach,momK0s);

      fillthis="histDCAtoPVvspK0S";
      ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momK0s,dcaV0ptp);
      if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momK0s,dcaV0ptp);

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

    if (isCandidateSelectedCuts) {

      fillthis="histoprotonBachSigmaVspTOF";
      ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(momBach,nSigmaTOFpr);
      fillthis="histoprotonBachSigmaVspTPC";
      ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(momBach,nSigmaTPCpr);


      fillthis="histpK0SvspOffline";
      ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momBach,momK0s);
      if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momBach,momK0s);

      fillthis="histDCAtoPVvspK0SOffline";
      ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momK0s,dcaV0ptp);
      if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momK0s,dcaV0ptp);

      fillthis="histOfflineArmPodK0s";
      FillArmPodDistribution(v0part,fillthis,fOutputAll);
      if (isBachelorID) FillArmPodDistribution(v0part,fillthis,fOutputPIDBach);

      fillthis="histLcMassOfflineByK0S";
      ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassLc,lambdacpt);
      if (isBachelorID)((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassLc,lambdacpt); // main histogram

    }
  }


  if (fUseMCInfo) {
    if (isLc==1) {
      if (onFlyV0 && fUseOnTheFlyV0) {

	fillthis="histK0SMassSgn";
	((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(invmassK0s,ptK0s);
	if (isBachelorID) ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(invmassK0s,ptK0s);

	if (isCandidateSelectedCuts) {

	  fillthis="histpK0SvspSgn";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momBach,momK0s);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momBach,momK0s);

	  fillthis="histDCAtoPVvspK0SSgn";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momK0s,dcaV0ptp);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momK0s,dcaV0ptp);

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

 	if (isCandidateSelectedCuts) {

	  fillthis="histoprotonBachSigmaVspTOFsgn";
	  ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(momBach,nSigmaTOFpr);
	  fillthis="histoprotonBachSigmaVspTPCsgn";
	  ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(momBach,nSigmaTPCpr);


	  fillthis="histpK0SvspOfflineSgn";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momBach,momK0s);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momBach,momK0s);
  
	  fillthis="histDCAtoPVvspK0SOfflineSgn";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momK0s,dcaV0ptp);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momK0s,dcaV0ptp);
   
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

	if (isCandidateSelectedCuts) {

	  fillthis="histpK0SvspBkg";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momBach,momK0s);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momBach,momK0s);

	  fillthis="histDCAtoPVvspK0SBkg";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momK0s,dcaV0ptp);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momK0s,dcaV0ptp);

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
  
	if (isCandidateSelectedCuts) {

	  fillthis="histoprotonBachSigmaVspTOFbkg";
	  ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(momBach,nSigmaTOFpr);
	  fillthis="histoprotonBachSigmaVspTPCbkg";
	  ((TH2F*)(fOutput->FindObject(fillthis)))->Fill(momBach,nSigmaTPCpr);


	  fillthis="histpK0SvspOfflineBkg";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momBach,momK0s);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momBach,momK0s);
   
	  fillthis="histDCAtoPVvspK0SOfflineBkg";
	  ((TH2F*)(fOutputAll->FindObject(fillthis)))->Fill(momK0s,dcaV0ptp);
	  if (isBachelorID)  ((TH2F*)(fOutputPIDBach->FindObject(fillthis)))->Fill(momK0s,dcaV0ptp);

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
void AliAnalysisTaskSELc2V0bachelor::MakeAnalysisForLc2prK0S(TClonesArray *arrayLctopKos,
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
    AliAODRecoCascadeHF* lcK0spr = dynamic_cast<AliAODRecoCascadeHF*>(arrayLctopKos->At(iLctopK0s));
    if (!lcK0spr) {
      AliDebug(2,Form("Cascade %d doens't exist, skipping",iLctopK0s));
      continue;
    }

    if (!lcK0spr->GetSecondaryVtx()) {
      AliInfo("No secondary vertex");
      continue;
    }

    if (lcK0spr->GetNDaughters()!=2) {
      AliDebug(2,Form("Cascade %d has not 2 daughters (nDaughters=%d)",iLctopK0s,lcK0spr->GetNDaughters()));
      continue;
    }

    AliAODv0 * v0part = dynamic_cast<AliAODv0*>(lcK0spr->Getv0());
    AliAODTrack * bachPart = dynamic_cast<AliAODTrack*>(lcK0spr->GetBachelor());
    if (!v0part || !bachPart) {
      AliDebug(2,Form("Cascade %d has no V0 or no bachelor object",iLctopK0s));
      continue;
    }


    if (!v0part->GetSecondaryVtx()) {
      AliDebug(2,Form("No secondary vertex for V0 by cascade %d",iLctopK0s));
      continue;
    }

    if (v0part->GetNDaughters()!=2) {
      AliDebug(2,Form("current V0 has not 2 daughters (onTheFly=%d, nDaughters=%d)",v0part->GetOnFlyStatus(),v0part->GetNDaughters()));
      continue;
    }

    AliAODTrack * v0Pos = dynamic_cast<AliAODTrack*>(lcK0spr->Getv0PositiveTrack());
    AliAODTrack * v0Neg = dynamic_cast<AliAODTrack*>(lcK0spr->Getv0NegativeTrack());
    if (!v0Neg || !v0Neg) {
      AliDebug(2,Form("V0 by cascade %d has no V0positive of V0negative object",iLctopK0s));
      continue;
    }

    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(1);

    if (v0Pos->Charge() == v0Neg->Charge()) continue;

    ((TH1I*)(fOutput->FindObject("hCandidateSelection")))->Fill(2);

    Int_t isLc = 0;

    if (fUseMCInfo) {

      Int_t pdgCode=-2;

      // find associated MC particle for Lc -> p+K0 and K0S->pi+pi
      Int_t mcLabelOld = MatchToMC(lcK0spr,pdgDgLctoV0bachelorOld,pdgDgV0toDaughters,mcArray);
      Int_t mcLabel = lcK0spr->MatchToMC(pdgCand,pdgDgLctoV0bachelor[0],pdgDgLctoV0bachelor,pdgDgV0toDaughters,mcArray,kTRUE);
      if (mcLabelOld!=mcLabel) AliDebug(2,Form(" Changed MC label: oldONE=%d wrt rightONE=%d",mcLabelOld,mcLabel));
      if (mcLabel>=0) {
	AliDebug(2,Form(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cascade number %d (total cascade number = %d)", iLctopK0s,nCascades));

	AliAODMCParticle *partLc = dynamic_cast<AliAODMCParticle*>(mcArray->At(mcLabel));
	pdgCode = partLc->GetPdgCode();
	if (pdgCode<0) AliDebug(2,Form(" ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ MClabel=%d ~~~~~~~~~~ pdgCode=%d", mcLabel, pdgCode));
	pdgCode = TMath::Abs(pdgCode);
	isLc = 1;
      } else {
	pdgCode=-1;
      }
    }

    FillLc2pK0Sspectrum(lcK0spr, isLc,
			nSelectedProd, cutsProd,
			nSelectedAnal, cutsAnal,
			mcArray);

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
    //AliInfo(Form(" bachelor mother label=%d - V0 mother label=%d",labBacMother, labV0mother));
    return -1;
  }

  //AliInfo(Form(" V0: %d (label=%d) ->  %d (label=%d) + %d (label=%d) --- GM %d (PDG=%d)",
  //motherV0->GetPdgCode(), labV0posMother,
  //partV0pos->GetPdgCode(), labV0pos,
  //partV0neg->GetPdgCode(), labV0neg,
  //labV0mother, gMotherV0->GetPdgCode()));

  //AliInfo(Form(" Bachelor: %d (label=%d) --- M %d (label=%d)",
  //partBachelor->GetPdgCode(), labBachelor,
  //labBacMother, partBacMother->GetPdgCode()));
  return labBacMother;//labV0mother;//

}

//________________________________________________________________
Int_t AliAnalysisTaskSELc2V0bachelor::SearchLcDaughter(TClonesArray *arrayMC, Int_t iii)
{

  Int_t indexToBeReturned=-999;

  Int_t pdgLc=4122;
  Int_t pdgLambda=3122;
  Int_t pdgV0=310;
  Int_t pdgK0=311;
  Int_t pdgBachelor=2212;
  Int_t pdgBachelorPi=211;

  TString fillthis="";
  fillthis="histMcStatLc";

  AliAODMCParticle *searchLc = dynamic_cast<AliAODMCParticle*>(arrayMC->At(iii));
  if(!searchLc) return -999;
  if (TMath::Abs(searchLc->GetPdgCode()) != pdgLc) return -999;

  ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(0);
  indexToBeReturned = 0;

  ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*1);
  indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*1;

  Int_t nDaughLc = searchLc->GetNDaughters();
  if (nDaughLc!=2) {
    ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*10);
    indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*10;
    return indexToBeReturned;
  }

  Int_t index1=searchLc->GetDaughter(0);
  Int_t index2=searchLc->GetDaughter(1);
  if (index1<=0 || index2<=0) {
    return -999;
  }

  AliAODMCParticle *daugh1 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(index1));
  AliAODMCParticle *daugh2 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(index2));
  if (!daugh1 || !daugh2) return -999;

  Int_t daughPdg1 = TMath::Abs(daugh1->GetPdgCode());
  Int_t daughPdg2 = TMath::Abs(daugh2->GetPdgCode());
  if ( !( (daughPdg1==pdgBachelor && daughPdg2==pdgK0) ||
	  (daughPdg2==pdgBachelor && daughPdg1==pdgK0) ||
	  (daughPdg1==pdgLambda && daughPdg2==pdgBachelorPi) ||
	  (daughPdg2==pdgLambda && daughPdg1==pdgBachelorPi) ) ) {
    ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*10);
    indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*10;
    return indexToBeReturned;
  }

  if (daughPdg1==pdgK0 || daughPdg1==pdgLambda) {
    index1=searchLc->GetDaughter(1);
    index2=searchLc->GetDaughter(0);
  }
  daugh1 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(index1));
  daugh2 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(index2));
  if (!daugh1 || !daugh2) return -999;

  daughPdg1=TMath::Abs(daugh1->GetPdgCode());
  daughPdg2=TMath::Abs(daugh2->GetPdgCode());

  if ( daughPdg1==pdgBachelor && daughPdg2==pdgK0 ) { // Lc+ -> p K0bar

    ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*2);
    indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*2;

    Int_t nDaughK0 = daugh2->GetNDaughters();
    if (nDaughK0!=1) return -999;

    Int_t indexK0daugh=daugh2->GetDaughter(0);
    if (indexK0daugh<=0) return -999;

    AliAODMCParticle *daughK0 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(indexK0daugh));
    if (!daughK0) return -999;

    Int_t daughK0Pdg=TMath::Abs(daughK0->GetPdgCode());
    if (daughK0Pdg!=pdgV0) {
      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*4); // K0L
      indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*4;
      return indexToBeReturned;
    }
    ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*3); // K0S
    indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*3;

    Int_t nDaughK0S = daughK0->GetNDaughters();
    if (nDaughK0S!=2) {
      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*5); // other decays for K0S
      indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*5;
      return indexToBeReturned;
    }

    index1=daughK0->GetDaughter(0);
    index2=daughK0->GetDaughter(1);
    if(index1<=0 || index2<=0) {
      return -999;
    }

    AliAODMCParticle *daughK0s1 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(index1));
    AliAODMCParticle *daughK0s2 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(index2));
    if (!daughK0s1 || !daughK0s2) return -999;

    Int_t daughK0s1pdg=TMath::Abs(daughK0s1->GetPdgCode());
    Int_t daughK0s2pdg=TMath::Abs(daughK0s2->GetPdgCode());

    if ( daughK0s1pdg==211 && daughK0s2pdg==211 ) {
      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*6); // K0S -> pi+ pi-
      indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*6;
    } else {
      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*5); // other decays for K0S
    indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*5;
    }

  } //if (daughPdg1==pdgBachelor && daughPdg2==pdgK0)
  else if ( daughPdg1==pdgBachelorPi && daughPdg2==pdgLambda ) { // Lc+ -> pi+ Lambda

    ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*7);
    indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*7;

    Int_t nDaughL = daugh2->GetNDaughters();
    if (nDaughL!=2) {
      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*8);
      indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*8;
      return indexToBeReturned;
    }

    index1=daugh2->GetDaughter(0);
    index2=daugh2->GetDaughter(1);
    if(index1<=0 || index2<=0) {
      return -999;
    }

    AliAODMCParticle *daughL1 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(index1));
    AliAODMCParticle *daughL2 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(index2));
    if (!daughL1 || !daughL2) return -999;

    Int_t daughL1pdg=TMath::Abs(daughL1->GetPdgCode());
    Int_t daughL2pdg=TMath::Abs(daughL2->GetPdgCode());
    if ( (daughL1pdg==211 && daughL2pdg==2212) ||
	 (daughL2pdg==211 && daughL1pdg==2212) ) {
      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*9);
      indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*9;
    } else {
      ((TH1F*)(fOutput->FindObject(fillthis)))->Fill(TMath::Nint(searchLc->Charge()/3.)*8);
      indexToBeReturned = TMath::Nint(searchLc->Charge()/3.)*8;
    }

  } //else if (daughPdg1==pdgBachelorPi && daughPdg2==pdgLambda)

  return indexToBeReturned;

}
//________________________________________________________________
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
    TH2F* spectrumK0SMass = new TH2F(nameMass.Data(),"K^{0}_{S} invariant mass VS p_{T}; M(#pi^{+}#pi^{-}) [GeV/c^{2}]; p_{T}(K^{0}_{S}) [GeV/c]; Entries",
				     520,0.43,0.56,200,0.,20.);

    // Lc invariant masses (x K0S on-the-fly)
    nameMass="histLcMassByK0S";
    TH2F* spectrumLcMassByK0S = new TH2F(nameMass.Data(),"#Lambda_{C} invariant mass (by K^{0}_{S}) vs p_{T} ; M(#Lambda_{C}) [GeV/c^{2}]; p_{T} [GeV/c]",
					 1200,2.,2.6,200,0.,20.);

    //nameMass="histcosOAK0Spvsp";
    //TH2F* cosOpeningAngleK0Spvsp = new TH2F(nameMass.Data(),"#Lambda_{C}: cosine of K^{0}_{S} - p opening angle vs #Lambda_{c} momentum; Cos(Opening Angle)  ; p [GeV/c]",
    //200,-1.,1.,200,0.,20.);

    nameMass="histpK0Svsp";
    TH2F* momentumDistributionK0Svsp = new TH2F(nameMass.Data(),"#Lambda_{C}: p(K^{0}_{S}) vs p(p);  p_{p}; p_{K^{0}_{S}}  ",
						200,0.,20.,200,0.,20.);

    nameMass="histArmPodK0s";
    TH2F* armenterosPodK0s = new TH2F(nameMass.Data(),"K^{0}_{S}: Armenteros-Podolanski distribution; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]",
				      200,-1.,1.,300,0.,0.3);
 
    nameMass="histDCAtoPVvspK0S";
    TH2F *dcatoPVvspK0s = new TH2F(nameMass.Data(),"K^{0}_{S}: DCA to Primary Vertex vs K^{0}_{S} momentum ; p(K^{0}_{S}) [GeV/c]; DCA to Primary Vertex []; Entries",
				   200,0.,20.,100,0.,10.);

    nameMass="histK0ScosPAwrtPVvspK0s";
    TH2F *cosPAwrtPVvspK0s = new TH2F(nameMass.Data(),"K^{0}_{S}: cosine of pointing angle wrt primary vertex vs K^{0}_{S} momentum ; p(K^{0}_{S}) [GeV/c]; cosine; Entries",
				      200,0.,20.,100,0.99,1.);

    TH2F* allspectrumK0SMass = (TH2F*)spectrumK0SMass->Clone(); 
    TH2F* allspectrumLcMassByK0S    = (TH2F*)spectrumLcMassByK0S->Clone(); 
    //TH2F* allcosOpeningAngleK0Spvsp= (TH2F*)cosOpeningAngleK0Spvsp->Clone(); 
    TH2F* allmomentumDistributionK0Svsp= (TH2F*)momentumDistributionK0Svsp->Clone(); 
    TH2F* alldcatoPVvspK0s=(TH2F*)dcatoPVvspK0s->Clone(); 
    TH2F* allcosV0PAwrtPVvspK0s=(TH2F*)cosPAwrtPVvspK0s->Clone(); 

    TH2F* pidBachspectrumK0SMass = (TH2F*)spectrumK0SMass->Clone(); 
    TH2F* pidBachspectrumLcMassByK0S    = (TH2F*)spectrumLcMassByK0S->Clone(); 
    //TH2F* pidBachcosOpeningAngleK0Spvsp= (TH2F*)cosOpeningAngleK0Spvsp->Clone(); 
    TH2F* pidBachmomentumDistributionK0Svsp= (TH2F*)momentumDistributionK0Svsp->Clone(); 
    TH2F* pidBachdcatoPVvspK0s=(TH2F*)dcatoPVvspK0s->Clone(); 
    TH2F* pidBachcosV0PAwrtPVvspK0s=(TH2F*)cosPAwrtPVvspK0s->Clone(); 

    TH2F* allArmenterosPodK0s = (TH2F*)armenterosPodK0s->Clone();
    TH2F* pidBachArmenterosPodK0s = (TH2F*)armenterosPodK0s->Clone();

    fOutputAll->Add(allspectrumK0SMass);
    fOutputAll->Add(allspectrumLcMassByK0S);
    //fOutputAll->Add(allcosOpeningAngleK0Spvsp); 
    fOutputAll->Add(allmomentumDistributionK0Svsp); 
    fOutputAll->Add(allArmenterosPodK0s);
    fOutputAll->Add(alldcatoPVvspK0s);
    fOutputAll->Add(allcosV0PAwrtPVvspK0s);

    fOutputPIDBach->Add(pidBachspectrumK0SMass);
    fOutputPIDBach->Add(pidBachspectrumLcMassByK0S);
    //fOutputPIDBach->Add(pidBachcosOpeningAngleK0Spvsp); 
    fOutputPIDBach->Add(pidBachmomentumDistributionK0Svsp); 
    fOutputPIDBach->Add(pidBachArmenterosPodK0s);
    fOutputPIDBach->Add(pidBachdcatoPVvspK0s);
    fOutputPIDBach->Add(pidBachcosV0PAwrtPVvspK0s);
 
  }

  // V0 invariant masses (offline)
  nameMass="histK0SMassOffline";
  TH2F* spectrumK0SMassOffline = new TH2F(nameMass.Data(),"K^{0}_{S} invariant mass VS p_{T}; M(#pi^{+}#pi^{-}) [GeV/c^{2}]; p_{T}(K^{0}_{S}) [GeV/c]; Entries",
					  520,0.43,0.56,200,0.,20.);

  // Lc invariant masses (x K0S offline)
  nameMass="histLcMassOfflineByK0S";
  TH2F* spectrumLcMassOfflineByK0S = new TH2F(nameMass.Data(),"#Lambda_{C} invariant mass (by K^{0}_{S}) vs p_{T}; M(K^{0}_{S}p) [GeV/c^{2}]; p_{T} [GeV/c]",
					      1200,2.,2.6,200,0.,20.);


  //nameMass="histcosOAK0SpvspOffline";
  //TH2F* cosOpeningAngleK0SpvspOffline = new TH2F(nameMass.Data(),"#Lambda_{C}: K^{0}_{S} - p opening angle vs p  -  Offline ; Cos(Opening Angle)  ; p [GeV/c]",
  //200,-1.,1.,200,0.,20.);

  nameMass="histpK0SvspOffline";
  TH2F* momentumDistributionK0SvspOffline = new TH2F(nameMass.Data(),"#Lambda_{C}: p(K^{0}_{S}) vs p(p) - Offline ;  p_{p} [GeV/c]; p_{K^{0}_{S}} [GeV/c]",
						    200,0.,20.,200,0.,20.);

  nameMass="histOfflineArmPodK0s";
  TH2F* armenterosPodK0sOff = new TH2F(nameMass.Data(),"K^{0}_{S}  Armenteros-Podolanski distribution - Offline; #frac{p_{L}^{+}-p_{L}^{-}}{p_{L}^{+}+p_{L}^{-}}; p_{T}^{+} [GeV/c]",
				       200,-1.,1.,300,0.,0.3);

  nameMass="histDCAtoPVvspK0SOffline";
  TH2F *dcatoPVvspK0sOffline = new TH2F(nameMass.Data(),"K^{0}_{S}: DCA to Primary Vertex vs  K^{0}_{S} invariant mass - Offline; p(K^{0}_{S}) [GeV/c]; DCA to Primary Vertex []; Entries",
					200,0.,20.,100,0.,10.);

  nameMass="histK0ScosPAwrtPVvspK0sOffline";
  TH2F *cosPAwrtPVvspK0sOffline = new TH2F(nameMass.Data(),"K^{0}_{S}: cosine of pointing angle wrt primary vertex vs K^{0}_{S} momentum - Offline; p(K^{0}_{S}) [GeV/c]; cosine; Entries",
					   200,0.,20.,100,0.99,1.);



  TH2F* allspectrumK0SMassOffline = (TH2F*)spectrumK0SMassOffline->Clone(); 
  TH2F* allspectrumLcMassOfflineByK0S    = (TH2F*)spectrumLcMassOfflineByK0S->Clone(); 
  //TH2F* allcosOpeningAngleK0SpvspOffline= (TH2F*)cosOpeningAngleK0SpvspOffline->Clone(); 
  TH2F* allmomentumDistributionK0SvspOffline= (TH2F*)momentumDistributionK0SvspOffline->Clone(); 
  TH2F* alldcatoPVvspK0sOffline=(TH2F*)dcatoPVvspK0sOffline->Clone(); 
  TH2F* allcosPAwrtPVvspK0sOffline=(TH2F*)cosPAwrtPVvspK0sOffline->Clone(); 

  TH2F* pidBachspectrumK0SMassOffline = (TH2F*)spectrumK0SMassOffline->Clone(); 
  TH2F* pidBachspectrumLcMassOfflineByK0S    = (TH2F*)spectrumLcMassOfflineByK0S->Clone(); 
  //TH2F* pidBachcosOpeningAngleK0SpvspOffline= (TH2F*)cosOpeningAngleK0SpvspOffline->Clone(); 
  TH2F* pidBachmomentumDistributionK0SvspOffline= (TH2F*)momentumDistributionK0SvspOffline->Clone(); 
  TH2F* pidBachdcatoPVvspK0sOffline=(TH2F*)dcatoPVvspK0sOffline->Clone(); 
  TH2F* pidBachcosPAwrtPVvspK0sOffline=(TH2F*)cosPAwrtPVvspK0sOffline->Clone(); 

  TH2F* allArmenterosPodK0sOff = (TH2F*)armenterosPodK0sOff->Clone();
  TH2F* pidBachArmenterosPodK0sOff = (TH2F*)armenterosPodK0sOff->Clone();


  fOutputAll->Add(allspectrumK0SMassOffline);
  fOutputAll->Add(allspectrumLcMassOfflineByK0S);
  //fOutputAll->Add(allcosOpeningAngleK0SpvspOffline); 
  fOutputAll->Add(allmomentumDistributionK0SvspOffline); 
  fOutputAll->Add(allArmenterosPodK0sOff);
  fOutputAll->Add(alldcatoPVvspK0sOffline);
  fOutputAll->Add(allcosPAwrtPVvspK0sOffline);

  fOutputPIDBach->Add(pidBachspectrumK0SMassOffline);
  fOutputPIDBach->Add(pidBachspectrumLcMassOfflineByK0S);
  //fOutputPIDBach->Add(pidBachcosOpeningAngleK0SpvspOffline); 
  fOutputPIDBach->Add(pidBachmomentumDistributionK0SvspOffline); 
  fOutputPIDBach->Add(pidBachArmenterosPodK0sOff);
  fOutputPIDBach->Add(pidBachdcatoPVvspK0sOffline);
  fOutputPIDBach->Add(pidBachcosPAwrtPVvspK0sOffline);

  /*
  nameMass="hist1LcMassOfflineByK0S";
  TH2D* h1 = new TH2D(nameMass.Data(),"IsEventSelected; M(#Lambda_{C}) [GeV/c^{2}]; p_{T}",1200,2.,2.6,200,0.,20.);
  fOutput->Add(h1);
  nameMass="hist0LcMassOfflineByK0S";
  TH2D* h0 = new TH2D(nameMass.Data(),"!IsEventSelected; M(#Lambda_{C}) [GeV/c^{2}]; p_{T}",1200,2.,2.6,200,0.,20.);
  fOutput->Add(h0);
  */

  if (fUseMCInfo) {

    if (fUseOnTheFlyV0) {

      nameSgn="histK0SMassSgn";
      nameBkg="histK0SMassBkg";
      TH2F* spectrumK0SMassSgn = new TH2F(nameSgn.Data(), "K^{0}_{S} Signal invariant mass VS p_{T} - MC; M(K^{0}_{S}) [GeV/c^{2}]; p_{T}(K^{0}_{S}) [GeV/c]; Entries", 520,0.43,0.56,200,0.,20.);
      TH2F* spectrumK0SMassBkg = new TH2F(nameBkg.Data(), "K^{0}_{S} Background invariant mass VS p_{T} - MC; M(K^{0}_{S}) [GeV/c^{2}]; p_{T}(K^{0}_{S}) [GeV/c]; Entries",  520,0.43,0.56,200,0.,20.);

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
					      1200,2.,2.6,200,0.,20.);
      TH2F* spectrumLcMassByK0SBkg = new TH2F(nameBkg.Data(), "#Lambda_{C} Background invariant mass (by K^{0}_{S}) vs p_{T}  - MC; M(#Lambda_{C}) [GeV/c^{2}]; p_{T}",
					      1200,2.,2.6,200,0.,20.);

      TH2F* allspectrumLcMassByK0SSgn = (TH2F*)spectrumLcMassByK0SSgn->Clone(); 
      TH2F* allspectrumLcMassByK0SBkg = (TH2F*) spectrumLcMassByK0SBkg->Clone();  
      TH2F* pidBachspectrumLcMassByK0SSgn = (TH2F*)spectrumLcMassByK0SSgn->Clone(); 
      TH2F* pidBachspectrumLcMassByK0SBkg = (TH2F*) spectrumLcMassByK0SBkg->Clone();  
      fOutputAll->Add(allspectrumLcMassByK0SSgn);
      fOutputAll->Add(allspectrumLcMassByK0SBkg);
      fOutputPIDBach->Add(pidBachspectrumLcMassByK0SSgn);
      fOutputPIDBach->Add(pidBachspectrumLcMassByK0SBkg);

      /*
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
      */

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


      nameSgn="histDCAtoPVvspK0SSgn";
      nameBkg="histDCAtoPVvspK0SBkg";
      TH2F *dcatoPVvspK0sSgn=new TH2F(nameSgn.Data(),"K^{0}_{S}: DCA to Primary Vertex vs  K^{0}_{S} invariant mass (sgn); p(K^{0}_{S}) [GeV/c]; DCA to Primary Vertex []; Entries",200,0.,20.,100,0.,10.0);
      TH2F *dcatoPVvspK0sBkg=new TH2F(nameBkg.Data(),"K^{0}_{S}: DCA to Primary Vertex vs  K^{0}_{S} invariant mass (bkg); p(K^{0}_{S}) [GeV/c]; DCA to Primary Vertex []; Entries",200,0.,20.,100,0.,10.0);

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
						   1200,2.,2.6,200,0.,20.);
    TH2F* spectrumLcMassOfflineByK0SBkg = new TH2F(nameBkg.Data(), "#Lambda_{C} Background invariant mass (by K^{0}_{S})  vs p_{T} - MC; M(#Lambda_{C}) [GeV/c^{2}]; p_{T}",
						   1200,2.,2.6,200,0.,20.);


    TH2F* allspectrumLcMassOfflineByK0SSgn = (TH2F*)spectrumLcMassOfflineByK0SSgn->Clone(); 
    TH2F* allspectrumLcMassOfflineByK0SBkg = (TH2F*) spectrumLcMassOfflineByK0SBkg->Clone();  
    TH2F* pidBachspectrumLcMassOfflineByK0SSgn = (TH2F*)spectrumLcMassOfflineByK0SSgn->Clone(); 
    TH2F* pidBachspectrumLcMassOfflineByK0SBkg = (TH2F*) spectrumLcMassOfflineByK0SBkg->Clone();  
    fOutputAll->Add(allspectrumLcMassOfflineByK0SSgn);
    fOutputAll->Add(allspectrumLcMassOfflineByK0SBkg);
    fOutputPIDBach->Add(pidBachspectrumLcMassOfflineByK0SSgn);
    fOutputPIDBach->Add(pidBachspectrumLcMassOfflineByK0SBkg);
  
 

    /*
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
    */



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


    nameSgn="histDCAtoPVvspK0SOfflineSgn";
    nameBkg="histDCAtoPVvspK0SOfflineBkg";
    TH2F *dcatoPVvspK0sOfflineSgn=new TH2F(nameSgn.Data(),"K^{0}_{S}: DCA to Primary Vertex vs  K^{0}_{S} invariant mass  -offline - (sgn); p(K^{0}_{S}) [GeV/c]; DCA to Primary Vertex []; Entries",200,0.,20.,100,0.,10.0);
    TH2F *dcatoPVvspK0sOfflineBkg=new TH2F(nameBkg.Data(),"K^{0}_{S}: DCA to Primary Vertex vs  K^{0}_{S} invariant mass  -offline - (bkg); p(K^{0}_{S}) [GeV/c]; DCA to Primary Vertex []; Entries",200,0.,20.,100,0.,10.0);
    

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
