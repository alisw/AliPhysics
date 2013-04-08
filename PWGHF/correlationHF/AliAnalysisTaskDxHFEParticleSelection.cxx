// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE Project            * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Sedat Altinpinar <Sedat.Altinpinar@cern.ch>           *
//*                  Hege Erdal       <hege.erdal@gmail.com>               *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliAnalysisTaskDxHFEParticleSelection.cxx
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-03-19
/// @brief  AnalysisTask electron selection for D0 - HFE correlation
///

#include "AliAnalysisTaskDxHFEParticleSelection.h"
#include "AliDxHFEParticleSelection.h"
#include "AliDxHFEParticleSelectionD0.h"
#include "AliDxHFEParticleSelectionMCD0.h"
#include "AliDxHFEParticleSelectionEl.h"
#include "AliDxHFEParticleSelectionMCEl.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisCuts.h"
#include "AliLog.h"
#include "TH1F.h"
#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliHFEtools.h"
#include "TChain.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "AliAODMCParticle.h"
#include "TFile.h"
#include "TList.h"
#include "TObjArray.h"
#include <memory>

using namespace std;

/// ROOT macro for the implementation of ROOT specific class methods
ClassImp(AliAnalysisTaskDxHFEParticleSelection)

AliAnalysisTaskDxHFEParticleSelection::AliAnalysisTaskDxHFEParticleSelection(const char* opt)
  : AliAnalysisTaskSE("AliAnalysisTaskDxHFEParticleSelection")
  , fOutput(0)
  , fOption(opt)
  , fCutList(NULL)
  , fCutsD0(NULL)
  , fSelector(NULL)
  , fUseMC(kFALSE)
  , fFillOnlyD0D0bar(0)
  , fSelectedTracks(NULL)
  , fMCArray(NULL)
  , fParticleType(kD0)
  , fSystem(0)
  , fSelectionParticleOptions()
  , fUseKine(kFALSE)
{
  // constructor
  //
  //

  DefineSlots();
}

int AliAnalysisTaskDxHFEParticleSelection::DefineSlots()
{
  // define the data slots
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  return 0;
}

AliAnalysisTaskDxHFEParticleSelection::~AliAnalysisTaskDxHFEParticleSelection()
{
  // destructor
  //
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor

  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
    fOutput = 0;
  }

  if (fSelector) {
    fSelector->Clear();
    delete fSelector;
    fSelector=NULL;
  }
  if(fSelectedTracks){
    delete fSelectedTracks;
    fSelectedTracks=NULL;
  }
  if(fMCArray) delete fMCArray;
  fMCArray=NULL;
  // external object, do not delete
  fCutsD0=NULL;

}

void AliAnalysisTaskDxHFEParticleSelection::UserCreateOutputObjects()
{
  // create result objects and add to output list

  Int_t iResult=0;

  fOutput = new TList;
  fOutput->SetOwner();

  ParseArguments(fOption.Data());

  // REVIEW: this has only effect if the list is deleted, it should
  // be done where the list is created, in any case, all objects of the
  // list are presumably external objects, they should not be deleted in
  // this class
  fCutList->SetOwner(); // NOt sure if needed

  // setting up for D0s
  if(fParticleType==kD0){

    if(fUseMC) fSelector=new AliDxHFEParticleSelectionMCD0(fOption);
    else fSelector=new AliDxHFEParticleSelectionD0(fOption);
    fSelector->SetCuts(fCutList,AliDxHFEParticleSelectionD0::kCutList);
    iResult=fSelector->Init();

    TObject *obj=NULL;
    TIter iter(fCutList);
    while((obj = iter())){
      fCutsD0=dynamic_cast<AliRDHFCuts*>(obj);
      if (!fCutsD0) {
	AliError(Form("Cut object is not of required type AliRDHFCuts but %s", obj->ClassName()));
      }
    }
  }
  

  // Setting up for electrons
  if(fParticleType==kElectron){
    if(fUseMC) fSelector=new AliDxHFEParticleSelectionMCEl(fOption);
    else fSelector=new AliDxHFEParticleSelectionEl(fOption);
    fSelector->SetCuts(fCutList, AliDxHFEParticleSelectionEl::kCutList);
    iResult=fSelector->Init();

    // If running on electrons, for now use RDHFCuts for event selection. 
    // Set up default cut object:
    AliRDHFCutsD0toKpi* cuts=new AliRDHFCutsD0toKpi();
    cuts->SetStandardCutsPP2010();
    fCutsD0=cuts;
  }

  if (iResult<0) {
    AliFatal(Form("initialization of worker class instance fElectrons failed with error %d", iResult));
  }


  // Retrieving the list containing histos and THnSparse
  // and storing them instead of fSelector
  // Fix to be able to merge
  TList *list =(TList*)fSelector->GetControlObjects();
  TObject *obj=NULL;

  TIter next(list);
  while((obj = next())){
    fOutput->Add(obj);
  }

  // all tasks must post data once for all outputs
  PostData(1, fOutput);
  PostData(2, fCutList);
}

void AliAnalysisTaskDxHFEParticleSelection::UserExec(Option_t* /*option*/)
{
  // process the event

  // TODO: implement correct input, this is likely not to be the
  // ESD
  TObject* pInput=InputEvent();
  if (!pInput) {
    AliError("failed to get input");
    return;
  }

  // check if input is an ESD
  AliVEvent *pEvent = dynamic_cast<AliVEvent*>(pInput);
  TClonesArray *inputArray=0;

  if(!pEvent && AODEvent() && IsStandardAOD()) { //Not sure if this is needed.. Keep it for now. 
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    pEvent = AODEvent();
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent* aodFromExt = ext->GetAOD();
      inputArray=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
    }
  } else if(pEvent) {
    inputArray=(TClonesArray*)pEvent->GetList()->FindObject("D0toKpi");
  }
  if(!inputArray || !pEvent) {
    AliError("Input branch not found!\n");
    return;
  }
  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!pEvent->GetPrimaryVertex() || TMath::Abs(pEvent->GetMagneticField())<0.001){
    AliDebug(2,"Rejected at GetPrimaryvertex");
    return;
  }

  fSelector->HistogramEventProperties(AliDxHFEParticleSelection::kEventsAll);

  //At the moment: Use AliRDHFCuts for event selection
  AliRDHFCuts* cutsd0=dynamic_cast<AliRDHFCuts*>(fCutsD0);
  if (!cutsd0) return; // Fatal thrown already in initialization

  if(!cutsd0->IsEventSelected(pEvent)) {
    AliDebug(2,"rejected at IsEventSelected");
    return;
  }

  fSelector->HistogramEventProperties(AliDxHFEParticleSelection::kEventsSel);

  if(fUseMC && fUseKine){
    fMCArray = dynamic_cast<TObjArray*>(pEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!fMCArray){
      AliError("Array of MC particles not found");
      return;
    }
  }

  if (fSelectedTracks) delete fSelectedTracks;

  if(fParticleType==kElectron){
    // Gets the PID response from the analysis manager
    AliPIDResponse *pidResponse = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler())->GetPIDResponse();
    if(!pidResponse){
      // TODO: consider issuing fatal instead of debug in case pidresponse not available
      AliDebug(1, "Using default PID Response");
      pidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class()); 
    }
  
    // Fetching the PID objects from the list, checks if the objects are AliHFEpids
    // If so, checks if they are initialized and also sets the pidresponse
    TObject *obj=NULL;
    TIter next(fCutList);
    while((obj = next())){
      if(strcmp(obj->ClassName(),"AliHFEpid")==0){
	AliHFEpid* pidObj=(AliHFEpid*)obj;
	if(!pidObj->IsInitialized()){
	  AliWarning("PID not initialised, get from Run no");
	  pidObj->InitializePID(pEvent->GetRunNumber());
	}
	pidObj->SetPIDResponse(pidResponse);
      }
    }

    // Also sends the pidresponse to the particle selection class for electron
    fSelector->SetPIDResponse(pidResponse); 

    if(fUseKine) fSelectedTracks=(fSelector->Select(fMCArray, pEvent));
    else fSelectedTracks=(fSelector->Select(pEvent));
  }
  else{
    if(fUseKine)  fSelectedTracks=(fSelector->Select(fMCArray,pEvent));
    else fSelectedTracks=(fSelector->Select(inputArray,pEvent));
  }

  Int_t nSelected = fSelectedTracks->GetEntriesFast();
  if(nSelected>0)
    fSelector->HistogramEventProperties(AliDxHFEParticleSelection::kEventsWithParticle);

  PostData(1, fOutput);
}

int AliAnalysisTaskDxHFEParticleSelection::ParseArguments(const char* arguments)
{
  // parse arguments and set internal flags
  TString strArguments(arguments);
  auto_ptr<TObjArray> tokens(strArguments.Tokenize(" "));
  if (!tokens.get()) return 0;

  TIter next(tokens.get());
  TObject* token;
  while ((token=next())) {
    TString argument=token->GetName();
   
    if (argument.BeginsWith("mc")) {
      fUseMC=true;
      continue;
    }
    if (argument.BeginsWith("system=")) {
      argument.ReplaceAll("system=", "");
      if (argument.CompareTo("pp")==0) {fSystem=0; }
      else if (argument.CompareTo("Pb-Pb")==0){ fSystem=1;}
      else if (argument.CompareTo("p-Pb")==0){ fSystem=2;}
      else {
	AliWarning(Form("can not set collision system, unknown parameter '%s'", argument.Data()));
	// TODO: check what makes sense
	fSystem=0;
      }
      continue;
    }
    if (argument.BeginsWith("usekine") || argument.BeginsWith("kine")) {
      fUseKine=true;
      AliInfo("Running on MC stack");
      continue;
    }
    if (argument.BeginsWith("fillD0scheme=")){
      argument.ReplaceAll("fillD0scheme=", "");
      if (argument.CompareTo("both")==0){ fFillOnlyD0D0bar=0;}
      else if (argument.CompareTo("D0")==0){ fFillOnlyD0D0bar=1;}
      else if (argument.CompareTo("D0bar")==0){ fFillOnlyD0D0bar=2;}
      else {
	AliWarning(Form("can not set D0 filling scheme, unknown parameter '%s'", argument.Data()));
	fFillOnlyD0D0bar=0;
      }
      continue;
    }
    if (argument.BeginsWith("particle=")){
      argument.ReplaceAll("particle=","");     
      if (argument.CompareTo("D0")==0){ fParticleType=kD0; AliInfo("Running over D0s ");}
      else if (argument.CompareTo("electron")==0){ fParticleType=kElectron; AliInfo("Running over Electrons ");}
      else{
	AliWarning(Form("can not set particle, unknown parameter '%s'. Particle set to D0", argument.Data()));
	fParticleType=kD0;
      }
      continue;
    }
    AliWarning(Form("unknown argument '%s'", argument.Data()));
      
  }

  return 0;
}

void AliAnalysisTaskDxHFEParticleSelection::FinishTaskOutput()
{
  // end of the processing
}

void AliAnalysisTaskDxHFEParticleSelection::Terminate(Option_t *)
{
  // last action on the client
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    AliFatal("failed to get output container");
    return;
  }

}
