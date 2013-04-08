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

/// @file   AliAnalysisTaskDxHFECorrelation.cxx
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-03-19
/// @brief  AnalysisTask D0 - HFE correlation
///

#include "AliAnalysisTaskDxHFECorrelation.h"
#include "AliDxHFECorrelation.h"
#include "AliDxHFECorrelationMC.h"
#include "AliDxHFEParticleSelectionD0.h"
#include "AliDxHFEParticleSelectionMCD0.h"
#include "AliDxHFEParticleSelectionEl.h"
#include "AliDxHFEParticleSelectionMCEl.h"
#include "AliDxHFEParticleSelection.h"
#include "AliHFCorrelator.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliVertexerTracks.h"
#include "AliAODHandler.h"
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliHFEcontainer.h"
#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEcuts.h"
#include "AliHFEtools.h"
#include "TObject.h"
#include "TChain.h"
#include "TSystem.h"
#include "AliReducedParticle.h"
#include "AliHFAssociatedTrackCuts.h" // initialization of event pool
#include "TFile.h"
#include <memory>

using namespace std;

/// ROOT macro for the implementation of ROOT specific class methods
ClassImp(AliAnalysisTaskDxHFECorrelation)

AliAnalysisTaskDxHFECorrelation::AliAnalysisTaskDxHFECorrelation(const char* opt)
  : AliAnalysisTaskSE("AliAnalysisTaskDxHFECorrelation")
  , fOutput(0)
  , fOption(opt)
  , fCorrelation(NULL)
  , fD0s(NULL)
  , fElectrons(NULL)
  , fCutsD0(NULL)
  , fCuts(NULL)
  , fUseMC(kFALSE)
  , fUseEventMixing(kFALSE)
  , fSystem(0)
  , fSelectedD0s(NULL)
  , fSelectedElectrons(NULL)
  , fListHFE(NULL)
  , fTriggerParticle(AliDxHFECorrelation::kD)
  , fUseKine(kFALSE)
  , fMCArray(NULL)
  , fCorrelationArguments("")
{
  // constructor
  //
  //
  DefineSlots();
}

int AliAnalysisTaskDxHFECorrelation::DefineSlots()
{
  // define the data slots
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2,AliRDHFCutsD0toKpi::Class());
  DefineOutput(3,TList::Class());
  DefineOutput(4,AliHFAssociatedTrackCuts::Class());
  return 0;
}

AliAnalysisTaskDxHFECorrelation::~AliAnalysisTaskDxHFECorrelation()
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
  if (fD0s) delete fD0s;
  fD0s=NULL;
  if (fElectrons) delete fElectrons;
  fElectrons=NULL;
  if (fCorrelation) delete fCorrelation;
  fCorrelation=NULL;
  // external object, do not delete
  fCutsD0=NULL;
  // external object, do not delete
  // TODO: also delete it here? this class is set as owner...
  fListHFE=NULL;
  if(fSelectedElectrons) delete fSelectedElectrons;
  fSelectedElectrons=NULL;
  if(fSelectedD0s) delete fSelectedD0s;
  fSelectedD0s=NULL;
  if(fListHFE) delete fListHFE;
  fListHFE=NULL;
  if(fMCArray) delete fMCArray;
  fMCArray=NULL;


}

void AliAnalysisTaskDxHFECorrelation::UserCreateOutputObjects()
{
  // create result objects and add to output list
  int iResult=0;
  // ParseArguments will also define the strings that are used as input
  // for the particle selection classes and the correlation class

  ParseArguments(fOption.Data());

  fOutput = new TList;
  fOutput->SetOwner();

  // D0s ===============================================
  if(fUseMC) fD0s=new AliDxHFEParticleSelectionMCD0(fOption);
  else fD0s=new AliDxHFEParticleSelectionD0(fOption);
  fD0s->SetCuts(fCutsD0,AliDxHFEParticleSelectionD0::kCutD0);
  iResult=fD0s->Init();
  if (iResult<0) {
    AliFatal(Form("initialization of worker class instance fD0s failed with error %d", iResult));
  }

  //Electrons ============================================
  fListHFE->SetOwner(); // Not sure if needed
  if(fUseMC) fElectrons=new AliDxHFEParticleSelectionMCEl(fOption);
  else fElectrons=new AliDxHFEParticleSelectionEl(fOption);
  fElectrons->SetCuts(fListHFE, AliDxHFEParticleSelectionEl::kCutList);
  iResult=fElectrons->Init();
  if (iResult<0) {
    AliFatal(Form("initialization of worker class instance fElectrons failed with error %d", iResult));
  }

  //Correlation ===========================================
  if(fUseMC) fCorrelation=new AliDxHFECorrelationMC;
  else fCorrelation=new AliDxHFECorrelation;
  fCorrelation->SetCuts(fCuts);
  iResult=fCorrelation->Init(fOption);
  if (iResult<0) {
    AliFatal(Form("initialization of worker class instance fCorrelation failed with error %d", iResult));
  }

  // Fix for merging:
  // Retrieving the individual objects created
  // and storing them instead of fD0s, fElectrons etc.. 
  TList *list =(TList*)fCorrelation->GetControlObjects();
  TObject *obj=NULL;

  TIter next(list);
  while((obj = next())){
    fOutput->Add(obj);
  }

  list=(TList*)fD0s->GetControlObjects();
  next=TIter(list);
  while((obj= next())){
    fOutput->Add(obj);
  }

  list=(TList*)fElectrons->GetControlObjects();
  next=TIter(list);
  while((obj = next()))
    fOutput->Add(obj);

  if (!fCutsD0) {
    AliFatal(Form("cut object for D0 missing"));
    return;
  }
 
  if (!dynamic_cast<AliRDHFCutsD0toKpi*>(fCutsD0)) {
    AliFatal(Form("cut object %s is of incorrect type %s, expecting AliRDHFCutsD0toKpi", fCutsD0->GetName(), fCutsD0->ClassName()));
    return;
  }
  // that's the copy for the output stream
  AliRDHFCutsD0toKpi* copyfCuts=new AliRDHFCutsD0toKpi(dynamic_cast<AliRDHFCutsD0toKpi&>(*fCutsD0));
  const char* nameoutput=GetOutputSlot(2)->GetContainer()->GetName();
  copyfCuts->SetName(nameoutput);

  // all tasks must post data once for all outputs
  PostData(1, fOutput);
  PostData(2,copyfCuts);
  PostData(3,fListHFE);
  PostData(4,fCuts);

}

void AliAnalysisTaskDxHFECorrelation::UserExec(Option_t* /*option*/)
{
  // process the event
  TObject* pInput=InputEvent();
  if (!pInput) {
    AliError("failed to get input");
    return;
  }
  AliVEvent *pEvent = dynamic_cast<AliVEvent*>(pInput);
  TClonesArray *inputArray=0;

  fCorrelation->HistogramEventProperties(AliDxHFECorrelation::kEventsAll);

  if(!pEvent && AODEvent() && IsStandardAOD()) { //Not sure if this is needed.. Keep it for now. 
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    pEvent = dynamic_cast<AliAODEvent*> (AODEvent());
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

  AliRDHFCuts* cutsd0=dynamic_cast<AliRDHFCuts*>(fCutsD0);
  if (!cutsd0) return; // Fatal thrown already in initialization

  if(!cutsd0->IsEventSelected(pEvent)) {
    // TODO: Fill histograms based on why the event is rejected
    AliDebug(2,"rejected at IsEventSelected");
    return;
  }

  if(fSystem==1){
    // Not really used anywhere, just as a test of centralityselection
    // (Could also be used in histo, so have not removed it)
    AliCentrality *centralityObj = 0;
    Double_t MultipOrCent = -1;
    AliAODEvent* aodEvent=dynamic_cast<AliAODEvent*>(pEvent);
    centralityObj = aodEvent->GetHeader()->GetCentralityP();
    MultipOrCent = centralityObj->GetCentralityPercentileUnchecked("V0M");
    AliInfo(Form("Centrality is %f", MultipOrCent));
  }

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
  TIter next(fListHFE);
  while((obj = next())){
    AliHFEpid* pidObj=dynamic_cast<AliHFEpid*>(obj);
    if(pidObj){
      if(!pidObj->IsInitialized()){
	AliWarning("PID not initialised, get from Run no");
	pidObj->InitializePID(pEvent->GetRunNumber());
      }
      pidObj->SetPIDResponse(pidResponse);
    }
  }

  // Also sends the pidresponse to the particle selection class for electron
  fElectrons->SetPIDResponse(pidResponse); 

  // Retrieving process from the AODMCHeader. 
  // TODO: Move it somewhere else? (keep it here for the moment since only need to read once pr event)
  if(fUseMC){
    AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(pEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    
    if (!mcHeader) {
      AliError("Could not find MC Header in AOD");
      return;
    }
    Int_t eventType = mcHeader->GetEventType();
    fCorrelation->SetEventType(eventType);
    if(fUseKine){
      fMCArray = dynamic_cast<TObjArray*>(pEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if(fUseMC && !fMCArray){
	AliError("Array of MC particles not found");
	return;
      }
    }
  }

  Int_t nInD0toKpi = inputArray->GetEntriesFast();

  fCorrelation->HistogramEventProperties(AliDxHFECorrelation::kEventsSel);

  Int_t nElSelected=0;
  // Run electron selection first if trigger particle is an electron
  if(fTriggerParticle==AliDxHFECorrelation::kElectron){

    if (fSelectedElectrons) delete fSelectedElectrons;
    // If run on kinematical level, send in MCarray instead of reconstructed tracks
    if(fUseKine) fSelectedElectrons=(fElectrons->Select(fMCArray, pEvent));
    else fSelectedElectrons=(fElectrons->Select(pEvent));
  
    if(! fSelectedElectrons) {
      return;
    }

    nElSelected =  fSelectedElectrons->GetEntriesFast();

    // No need to go further if no electrons are found, except for event mixing. Will here anyway correlate D0s with electrons from previous events
    if(!fUseEventMixing && nElSelected==0){
      AliDebug(4,"No electrons found in this event");
      return;
    }
    // Fill bin with triggered events if electrons are the trigger (only for events with nr electrons >0)
    if(nElSelected>0) fCorrelation->HistogramEventProperties(AliDxHFECorrelation::kEventsTriggered);

  }


  // D0 selection 
  if(fSelectedD0s) delete fSelectedD0s;
  if(fUseKine)  fSelectedD0s=(fD0s->Select(fMCArray,pEvent));
  else fSelectedD0s=(fD0s->Select(inputArray,pEvent));
  
  if(! fSelectedD0s) {
    return;
  }
  Int_t nD0Selected = fSelectedD0s->GetEntriesFast();

  // When not using EventMixing, no need to go further if no D0s are found.
  // For Event Mixing, need to store all found electrons in the pool
  if(!fUseEventMixing && nD0Selected==0){
    AliDebug(4,"No D0s found in this event");
    return;
  }
  //  Run electron selection second if trigger particle is D0
  if(fTriggerParticle!=AliDxHFECorrelation::kElectron){

    // Fill bin with triggered events here if D0 are the trigger (only for events with nr D0 >0)
    if(nD0Selected>0) fCorrelation->HistogramEventProperties(AliDxHFECorrelation::kEventsTriggered);

    if (fSelectedElectrons) delete fSelectedElectrons;

    // If run on kinematical level, send in MCarray instead of reconstructed tracks
    if(fUseKine) fSelectedElectrons=(fElectrons->Select(fMCArray, pEvent));
    else fSelectedElectrons=(fElectrons->Select(pEvent));
    if(! fSelectedElectrons) {
      return;
    }

    nElSelected =  fSelectedElectrons->GetEntriesFast();

    // No need to go further if no electrons are found, except for event mixing. Will here anyway correlate D0s with electrons from previous events
    if(!fUseEventMixing && nElSelected==0){
      AliDebug(4,"No electrons found in this event");
      return;
    }
  }

  // Should not be necessary:
  if(!fUseEventMixing && (nD0Selected==0 && nElSelected==0)){
    AliDebug(4,"Neither D0 nor electrons in this event");
    return;
  }
  
  AliDebug(4,Form("Number of D0->Kpi Start: %d , End: %d    Electrons Selected: %d\n", nInD0toKpi, nD0Selected, nElSelected));

  if(nD0Selected >0 && nElSelected>0) fCorrelation->HistogramEventProperties(AliDxHFECorrelation::kEventsCorrelated);

  int iResult=0;
  if(fTriggerParticle==AliDxHFECorrelation::kD) 
    fCorrelation->Fill(fSelectedD0s, fSelectedElectrons, pEvent);
  else 
    fCorrelation->Fill(fSelectedElectrons, fSelectedD0s, pEvent);

  if (iResult<0) {
    AliFatal(Form("%s processing failed with error %d", fCorrelation->GetName(), iResult));
  }

  PostData(1, fOutput);
  return;

}

int AliAnalysisTaskDxHFECorrelation::ParseArguments(const char* arguments)
{
  // parse arguments and set internal flags
  TString strArguments(arguments);
  auto_ptr<TObjArray> tokens(strArguments.Tokenize(" "));
  if (!tokens.get()) return 0;
  AliInfo(strArguments);
  TIter next(tokens.get());
  TObject* token;
  while ((token=next())) {
    TString argument=token->GetName();
   
    if (argument.BeginsWith("event-mixing")) {
      fUseEventMixing=true;
      AliInfo("Running with Event mixing");
      continue;
    }
      
    if (argument.BeginsWith("mc")) {
      fUseMC=true;
      AliInfo("Running on MC data");
      continue;
    }
    if (argument.BeginsWith("usekine") || argument.BeginsWith("kine")) {
      fUseKine=true;
      AliInfo("Running on MC stack");
      continue;
    }
    if (argument.BeginsWith("system=")) {
      argument.ReplaceAll("system=", "");
      if (argument.CompareTo("pp")==0) {fSystem=0;}
      else if (argument.CompareTo("Pb-Pb")==0){ fSystem=1;}
      else if (argument.CompareTo("p-Pb")==0){ fSystem=2;}
      else {
	AliWarning(Form("can not set collision system, unknown parameter '%s'", argument.Data()));
	// TODO: check what makes sense
	fSystem=0;
      }
      continue;
    }
    if (argument.BeginsWith("trigger=")){
      argument.ReplaceAll("trigger=", "");
      if (argument.CompareTo("D0")==0) {fTriggerParticle=AliDxHFECorrelation::kD; AliInfo("CorrTask: trigger on D0");}
      else if (argument.CompareTo("D")==0){ fTriggerParticle=AliDxHFECorrelation::kD; AliInfo("CorrTask: trigger on D0");}
      else if (argument.CompareTo("electron")==0) { fTriggerParticle=AliDxHFECorrelation::kElectron; AliInfo("CorrTask: trigger on electron");}
      continue;
    }
    AliWarning(Form("unknown argument '%s'", argument.Data()));
      
  }

  return 0;
}

void AliAnalysisTaskDxHFECorrelation::FinishTaskOutput()
{
  // end of the processing
}

void AliAnalysisTaskDxHFECorrelation::Terminate(Option_t *)
{
  // last action on the client
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    // looks like that is a valid condition if the task is run
    // in mode "terminate"
    return;
  }
}
