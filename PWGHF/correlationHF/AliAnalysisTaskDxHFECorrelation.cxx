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
#include "TProfile.h"
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
  , fQASelection(0)
  , fOption(opt)
  , fCorrelation(NULL)
  , fCorrelationCharm(NULL)
  , fCorrelationBeauty(NULL)
  , fCorrelationNonHF(NULL)
  , fCorrelationHadron(NULL)
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
  , fStoreSeparateOrigins(kFALSE)
  , fReqD0InEvent(kFALSE)
  , fpidResponse(NULL)
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
  DefineOutput(5, TList::Class());

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
  if (fQASelection && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fQASelection;
    fQASelection = 0;
  }
  if (fD0s) delete fD0s;
  fD0s=NULL;
  if (fElectrons) delete fElectrons;
  fElectrons=NULL;
  if (fCorrelation) delete fCorrelation;
  fCorrelation=NULL;
  if (fCorrelationCharm) delete fCorrelationCharm;
  fCorrelationCharm=NULL;
  if (fCorrelationBeauty) delete fCorrelationBeauty;
  fCorrelationBeauty=NULL;
  if (fCorrelationNonHF) delete fCorrelationNonHF;
  fCorrelationNonHF=NULL;
  if (fCorrelationHadron) delete fCorrelationHadron;
  fCorrelationHadron=NULL;
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

  if(fpidResponse) delete fpidResponse;
  fpidResponse=NULL;

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

  fQASelection = new TList;
  fQASelection->SetName("QA histos");
  fQASelection->SetOwner();
  // Gets the PID response from the input handler
  fpidResponse = fInputHandler->GetPIDResponse();
  if(!fpidResponse){
    // TODO: consider issuing fatal instead of debug in case pidresponse not available
    AliDebug(1, "Using default PID Response");
    fpidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class()); 
  }

  // D0s ===============================================
  if(fUseMC) fD0s=new AliDxHFEParticleSelectionMCD0(fOption);
  else fD0s=new AliDxHFEParticleSelectionD0(fOption);
  fD0s->SetCuts(fCutsD0,AliDxHFEParticleSelectionD0::kCutList);
  iResult=fD0s->Init();
  if (iResult<0) {
    AliFatal(Form("initialization of worker class instance fD0s failed with error %d", iResult));
  }

  //Electrons ============================================
  fListHFE->SetOwner(); // Not sure if needed
  if(fUseMC) fElectrons=new AliDxHFEParticleSelectionMCEl(fOption);
  else fElectrons=new AliDxHFEParticleSelectionEl(fOption);
  fElectrons->SetCuts(fListHFE, AliDxHFEParticleSelectionEl::kCutList);
  fElectrons->SetPIDResponse(fpidResponse);
  iResult=fElectrons->Init();
  if (iResult<0) {
    AliFatal(Form("initialization of worker class instance fElectrons failed with error %d", iResult));
  }

  //Correlation ===========================================
  if(fUseMC) fCorrelation=new AliDxHFECorrelationMC;
  else fCorrelation=new AliDxHFECorrelation;
  fCorrelation->SetCuts(fCuts);
  fCorrelation->SetCutsD0(fCutsD0);
  iResult=fCorrelation->Init(fOption);
  if (iResult<0) {
    AliFatal(Form("initialization of worker class instance fCorrelation failed with error %d", iResult));
  }

  if(fUseMC && fStoreSeparateOrigins){
    TString option;
    
    //Correlation only charm===========================================
    fCorrelationCharm=new AliDxHFECorrelationMC("AliDxHFECorrelationMCCharm");
    fCorrelationCharm->SetCuts(fCuts);
    fCorrelationCharm->SetCutsD0(fCutsD0);
    option=fOption+" storeoriginD=charm storeoriginEl=charm";
    iResult=fCorrelationCharm->Init(option);
    if (iResult<0) {
      AliFatal(Form("initialization of worker class instance fCorrelation failed with error %d", iResult));
    }

    //Correlation only beauty==========================================
    fCorrelationBeauty=new AliDxHFECorrelationMC("AliDxHFECorrelationMCBeauty");
    fCorrelationBeauty->SetCuts(fCuts);
    fCorrelationBeauty->SetCutsD0(fCutsD0);
    option=fOption+" storeoriginD=beauty storeoriginEl=beauty";
    iResult=fCorrelationBeauty->Init(option);
    if (iResult<0) {
      AliFatal(Form("initialization of worker class instance fCorrelation failed with error %d", iResult));
    }

    //Correlation nonHFE with c+b Ds====================================
    fCorrelationNonHF=new AliDxHFECorrelationMC("AliDxHFECorrelationMCNonHF");
    fCorrelationNonHF->SetCuts(fCuts);
    fCorrelationNonHF->SetCutsD0(fCutsD0);
    option=fOption+" storeoriginD=HF storeoriginEl=nonHF";
    iResult=fCorrelationNonHF->Init(option);
    if (iResult<0) {
      AliFatal(Form("initialization of worker class instance fCorrelation failed with error %d", iResult));
    }

    //Correlation hadrons with c+b Ds====================================
    fCorrelationHadron=new AliDxHFECorrelationMC("AliDxHFECorrelationMCHadrons");
    fCorrelationHadron->SetCuts(fCuts);
    fCorrelationHadron->SetCutsD0(fCutsD0);
    option=fOption+" storeoriginD=HF storeoriginEl=hadrons";
    iResult=fCorrelationHadron->Init(option);
    if (iResult<0) {
      AliFatal(Form("initialization of worker class instance fCorrelation failed with error %d", iResult));
    }
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
  if(fUseMC && fStoreSeparateOrigins){
    
    list=(TList*)fCorrelationCharm->GetControlObjects();
    next=TIter(list);
    while((obj= next())){
      fOutput->Add(obj);
    }

    list=(TList*)fCorrelationBeauty->GetControlObjects();
    next=TIter(list);
    while((obj= next())){
      fOutput->Add(obj);
    }

    list=(TList*)fCorrelationNonHF->GetControlObjects();
    next=TIter(list);
    while((obj= next())){
      fOutput->Add(obj);
    }

    list=(TList*)fCorrelationHadron->GetControlObjects();
    next=TIter(list);
    while((obj= next())){
      fOutput->Add(obj);
    }
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
    AliFatal(Form("cut object list for D0 missing"));
    return;
  }

  //[FIXME]Under development. At the moment the run ranges are set to the LHC13b and c ranges
  TProfile * numD0EvtRnAll = new TProfile("numD0EvtRnAll", "Number of D0 per Event vs Run Number, All Events", 400, 195300, 195700, 0, 1000);
  TProfile * numD0EvtRnCan = new TProfile("numD0EvtRnCan", "Number of D0 per Event vs Run Number, All D0toKpi", 400, 195300, 195700, 0, 1000);
  TProfile * numD0EvtRnSel = new TProfile("numD0EvtRnSel", "Number of D0 per Event vs Run Number, Selected D0", 400, 195300, 195700, 0, 1000);
  TProfile * numD0EvtRnCorr = new TProfile("numD0EvtRnCorr", "Number of D0 per Event vs Run Number, Correlated D0", 400, 195300, 195700, 0, 1000);

  TProfile * numElEvtRnSel = new TProfile("numElEvtRnSel", "Number of El per Event vs Run Number, Selected El", 400, 195300, 195700, 0, 1000);
  TProfile * numElEvtRnTrig = new TProfile("numElEvtRnTrig", "Number of El per Event vs Run Number, Triggered El", 400, 195300, 195700, 0, 1000);
  TProfile * numElEvtRnCorr = new TProfile("numElEvtRnCorr", "Number of El per Event vs Run Number, Correlated El", 400, 195300, 195700, 0, 1000);

  fQASelection->Add(numD0EvtRnAll);
  fQASelection->Add(numD0EvtRnCan);
  fQASelection->Add(numD0EvtRnSel);
  fQASelection->Add(numD0EvtRnCorr);

  fQASelection->Add(numElEvtRnSel);
  fQASelection->Add(numElEvtRnTrig);
  fQASelection->Add(numElEvtRnCorr);

  // Fetching out the RDHF-cut objects for D0 to store in output stream
  TObject *obj2=NULL;
  TIter next2(fCutsD0);
  obj2 = next2();
  AliRDHFCutsD0toKpi* copyfCuts=new AliRDHFCutsD0toKpi(dynamic_cast<AliRDHFCutsD0toKpi&>(*obj2));
  if(!copyfCuts){
    AliFatal(Form("cut object is of incorrect type %s, expecting AliRDHFCutsD0toKpi", obj->ClassName()));
  }

  const char* nameoutput=GetOutputSlot(2)->GetContainer()->GetName();
  copyfCuts->SetName(nameoutput);
  // all tasks must post data once for all outputs
  PostData(1, fOutput);
  PostData(2,copyfCuts);
  PostData(3,fListHFE);
  PostData(4,fCuts);
  PostData(5,fQASelection);
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
  Int_t runNumber=pEvent->GetRunNumber();
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
    ((TProfile*)fQASelection->FindObject("numD0EvtRnAll"))->Fill(runNumber, inputArray->GetEntriesFast());
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
  TObject *obj2=NULL;
  TIter next2(fCutsD0);
  obj2 = next2();
  AliRDHFCuts* cutsd0=dynamic_cast<AliRDHFCuts*>(obj2);
  //  AliRDHFCutsD0toKpi* cutsd0=dynamic_cast<AliRDHFCutsD0toKpi&>(*obj2);
  if (!cutsd0) return; // Fatal thrown already in initialization

  // Retrieving process from the AODMCHeader. 
  // TODO: Move it somewhere else? (keep it here for the moment since only need to read once pr event)
  if(fUseMC){
    AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(pEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    AliAODEvent* aod=dynamic_cast<AliAODEvent*>(pEvent);
    if (!mcHeader) {
      AliError("Could not find MC Header in AOD");
      return;
    }
    Int_t eventType = mcHeader->GetEventType();
    fCorrelation->SetEventType(eventType);
    if(fUseKine || fReqD0InEvent){
      fMCArray = dynamic_cast<TObjArray*>(pEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if(fUseMC && !fMCArray){
	AliError("Array of MC particles not found");
	return;
      }
    }

    if(!fUseKine && !cutsd0->IsEventSelected(pEvent)) { //Using different approach for fUseKine
      // TODO: Fill histograms based on why the event is rejected
      AliDebug(2,"rejected at IsEventSelected");
      return;
    }
    //On Kine, instead of IsEventSelected just select on zVtx and trigger mask in pPb
    if(fUseKine){
      Double_t zVtxMC = mcHeader->GetVtxZ();
      if(TMath::Abs(zVtxMC)>10) return;
      if(aod->GetTriggerMask()==0 && (aod->GetRunNumber()>=195344 && aod->GetRunNumber()<=195677)) return;
    }
  }else{
    if(!cutsd0->IsEventSelected(pEvent)) { //Using different approach for fUseKine
      // TODO: Fill histograms based on why the event is rejected
      AliDebug(2,"rejected at IsEventSelected");
      return;
    }
  }

  if(fSystem==1){
    // Not really used anywhere, just as a test of centralityselection
    // (Could also be used in histo, so have not removed it)
    AliCentrality *centralityObj = 0;
    Double_t MultipOrCent = -1;
    AliAODEvent* aodEvent=dynamic_cast<AliAODEvent*>(pEvent);
    if (aodEvent) {
      centralityObj = ((AliVAODHeader*)aodEvent->GetHeader())->GetCentralityP();
      if (centralityObj) {
	MultipOrCent = centralityObj->GetCentralityPercentileUnchecked("V0M");
      }
    }
    AliInfo(Form("Centrality is %f", MultipOrCent));
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
      pidObj->SetPIDResponse(fpidResponse);
    }
  }
  
  Int_t nInD0toKpi = inputArray->GetEntriesFast();
  ((TProfile*)fQASelection->FindObject("numD0EvtRnCan"))->Fill(runNumber, nInD0toKpi);

  //For studies reco/kine ratio -> select only events where there is a MC truth D0
  if(fUseMC && fReqD0InEvent){
    TIter itr(fMCArray);
    TObject* otr=NULL;
    Bool_t selectEvent=kFALSE;
    while ((otr=itr())!=NULL) {
      AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(otr);
      if(!mcPart) continue;
      Int_t PDG =TMath::Abs(mcPart->PdgCode()); 
      if(PDG==421) {selectEvent=kTRUE;}
    }
    if(!selectEvent) return;

  }

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
    ((TProfile*)fQASelection->FindObject("numElEvtRnTrig"))->Fill(runNumber, nElSelected);

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
  ((TProfile*)fQASelection->FindObject("numD0EvtRnSel"))->Fill(runNumber, nD0Selected);

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
    ((TProfile*)fQASelection->FindObject("numElEvtRnSel"))->Fill(runNumber, nElSelected);

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
  if(nD0Selected >0 && nElSelected>0){
    fCorrelation->HistogramEventProperties(AliDxHFECorrelation::kEventsCorrelated);
    ((TProfile*)fQASelection->FindObject("numD0EvtRnCorr"))->Fill(runNumber, nD0Selected);
    ((TProfile*)fQASelection->FindObject("numElEvtRnCorr"))->Fill(runNumber, nElSelected);
  }
  int iResult=0;
  if(fTriggerParticle==AliDxHFECorrelation::kD){ 
    fCorrelation->Fill(fSelectedD0s, fSelectedElectrons, pEvent);
    if(fUseMC && fStoreSeparateOrigins){
      fCorrelationCharm->Fill(fSelectedD0s, fSelectedElectrons, pEvent);
      fCorrelationBeauty->Fill(fSelectedD0s, fSelectedElectrons, pEvent);
      fCorrelationNonHF->Fill(fSelectedD0s, fSelectedElectrons, pEvent);
      fCorrelationHadron->Fill(fSelectedD0s, fSelectedElectrons, pEvent);
    }
  }
  else {
    fCorrelation->Fill(fSelectedElectrons, fSelectedD0s, pEvent);
    if(fUseMC && fStoreSeparateOrigins){
      fCorrelationCharm->Fill(fSelectedElectrons, fSelectedD0s, pEvent);
      fCorrelationBeauty->Fill(fSelectedElectrons, fSelectedD0s, pEvent);
      fCorrelationNonHF->Fill(fSelectedElectrons, fSelectedD0s, pEvent);
      fCorrelationHadron->Fill(fSelectedElectrons, fSelectedD0s, pEvent);
    }
  }
  if (iResult<0) {
    AliFatal(Form("%s processing failed with error %d", fCorrelation->GetName(), iResult));
  }

  PostData(1, fOutput);
  PostData(5, fQASelection);
  return;

}

int AliAnalysisTaskDxHFECorrelation::ParseArguments(const char* arguments)
{
  // parse arguments and set internal flags
  TString strArguments(arguments);
  unique_ptr<TObjArray> tokens(strArguments.Tokenize(" "));
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
    if (argument.BeginsWith("storeseparateorigins")) {
      fStoreSeparateOrigins=true;
      AliInfo("Store separate origins");
      continue;
    }
    if (argument.BeginsWith("usekine") || argument.BeginsWith("kine")) {
      fUseKine=true;
      AliInfo("Running on MC stack");
      continue;
    }
    if (argument.BeginsWith("reqD0inevent") || argument.BeginsWith("ReqD0InEvent")) {
      fReqD0InEvent=true;
      AliInfo("Requiring MC truth D0 in event");
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
