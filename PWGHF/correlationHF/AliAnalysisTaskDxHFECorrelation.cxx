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
#include "AliDxHFEParticleSelectionD0.h"
#include "AliDxHFEParticleSelectionMCD0.h"
#include "AliDxHFEParticleSelectionEl.h"
#include "AliDxHFEParticleSelectionMCEl.h"
#include "AliDxHFEParticleSelection.h"
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
  , fCutsHFE(NULL)
  , fPID(NULL)
  , fFillOnlyD0D0bar(0)
  , fUseMC(kFALSE)
{
  // constructor
  //
  //
  DefineSlots();
  fPID = new AliHFEpid("hfePid");

}

int AliAnalysisTaskDxHFECorrelation::DefineSlots()
{
  // define the data slots
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2,AliRDHFCutsD0toKpi::Class());
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
  fCutsHFE=NULL;
  if(fPID) delete fPID;
  fPID=NULL;


}

void AliAnalysisTaskDxHFECorrelation::UserCreateOutputObjects()
{
  // create result objects and add to output list

  //Initialize PID for electron selection
  if(!fPID->GetNumberOfPIDdetectors()) { 
    fPID->AddDetector("TOF",0);
    fPID->AddDetector("TPC",1);
  }
  fPID->InitializePID();

  fOutput = new TList;
  fOutput->SetOwner();

  // setting up for D0s
  TString selectionD0Options;
  switch (fFillOnlyD0D0bar) {
  case 1: selectionD0Options+="FillOnlyD0 "; break;
  case 2: selectionD0Options+="FillOnlyD0bar "; break;
  default: selectionD0Options+="FillD0D0bar ";
  }

  if(fUseMC) fD0s=new AliDxHFEParticleSelectionMCD0(selectionD0Options);
  else fD0s=new AliDxHFEParticleSelectionD0(selectionD0Options);
  fD0s->SetCuts(fCutsD0);
  fD0s->Init();

  //Electrons
  if(fUseMC) fElectrons=new AliDxHFEParticleSelectionMCEl;
  else fElectrons=new AliDxHFEParticleSelectionEl;
  fElectrons->SetCuts(fPID, AliDxHFEParticleSelectionEl::kCutPID);
  fElectrons->SetCuts(fCutsHFE, AliDxHFEParticleSelectionEl::kCutHFE);
  fElectrons->Init();

  //Correlation
  fCorrelation=new AliDxHFECorrelation;
  fCorrelation->SetCuts(dynamic_cast<AliRDHFCutsD0toKpi*>(fCutsD0));
  // TODO: check if we can get rid of the mc flag in the correlation analysis class
  fCorrelation->SetUseMC(fUseMC);
  fCorrelation->Init();

  // Fix for merging:
  // Retrieving the individual objects created
  // and storing them instead of fD0s, fElectrons etc.. 
  TList *list =(TList*)fD0s->GetControlObjects();
  TObject *obj=NULL;

  TIter next(list);
  while((obj = next())){
    fOutput->Add(obj);
  }

  list=(TList*)fCorrelation->GetControlObjects();
  next=TIter(list);
  while((obj= next())){
    fOutput->Add(obj);
  }

  list=(TList*)fElectrons->GetControlObjects();
  next=TIter(list);
  while((obj = next()))
    fOutput->Add(obj);

  PostData(1, fOutput);
  if (fCutsD0) {
    AliRDHFCutsD0toKpi* cuts=dynamic_cast<AliRDHFCutsD0toKpi*>(fCutsD0);
    if (!cuts) {
      AliFatal(Form("cut object %s is of incorrect type %s, expecting AliRDHFCutsD0toKpi", fCutsD0->GetName(), fCutsD0->ClassName()));
      return;
    }
    // TODO: why copy? cleanup?
    AliRDHFCutsD0toKpi* copyfCuts=new AliRDHFCutsD0toKpi(*cuts);
    const char* nameoutput=GetOutputSlot(2)->GetContainer()->GetName();
    copyfCuts->SetName(nameoutput);

    PostData(2,copyfCuts);
  }
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

  fCorrelation->HistogramEventProperties(AliDxHFECorrelation::kEventsAll);
  AliRDHFCuts* cutsd0=dynamic_cast<AliRDHFCuts*>(fCutsD0);
  if (!cutsd0) return; // Fatal thrown already in initialization

  if(!cutsd0->IsEventSelected(pEvent)) {
    AliDebug(2,"rejected at IsEventSelected");
    return;
  }

  if(!fPID->IsInitialized()){ 
    // Initialize PID with the given run number
    AliWarning("PID not initialised, get from Run no");
    fPID->InitializePID(pEvent->GetRunNumber());
  }

  AliPIDResponse *pidResponse = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler())->GetPIDResponse();
  if(!pidResponse){
    AliDebug(1, "Using default PID Response");
    pidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class()); 
  }
 
  fPID->SetPIDResponse(pidResponse);


  Int_t nInD0toKpi = inputArray->GetEntriesFast();

  fCorrelation->HistogramEventProperties(AliDxHFECorrelation::kEventsSel);

  //TODO: pSelectedD0s will now contain both D0 and D0bar. Probably add argument to Select
  // which determines which ones to return. 
  std::auto_ptr<TObjArray> pSelectedD0s(fD0s->Select(inputArray,pEvent));
  if(! pSelectedD0s.get()) {
    return;
  }
  Int_t nD0Selected = pSelectedD0s->GetEntriesFast();

  // No need to go further if no D0s are found. Not sure if this is the best way though..
  // At the moment this means no control histos can be implemented in AliDxHFECorrelation
  if(nD0Selected==0){
    //AliInfo("No D0s found in this event");
    return;
  }
  fCorrelation->HistogramEventProperties(AliDxHFECorrelation::kEventsD0);


  std::auto_ptr<TObjArray> pSelectedElectrons(fElectrons->Select(pEvent));

  // TODO: use the array of selected track for something, right now
  // only the control histograms of the selection class are filled
  // note: the pointer is deleted automatically once the scope is left
  // if the array should be published, the auto pointer must be released
  // first, however some other cleanup will be necessary in that case
  // probably a clone with a reduced AliVParticle implementation is
  // appropriate.


  if(! pSelectedElectrons.get()) {
    return;
  }

  Int_t nElSelected =  pSelectedElectrons->GetEntriesFast();

  // No need to go further if no electrons are found. Not sure if this is the best way though..
  // At the moment this means no control histos can be implemented in AliDxHFECorrelation
  if(nElSelected==0){
    //AliInfo("No electrons found in this event");
    return;
  }

  AliDebug(4,Form("Number of D0->Kpi Start: %d , End: %d    Electrons Selected: %d\n", nInD0toKpi, nD0Selected, nElSelected));

  fCorrelation->HistogramEventProperties(AliDxHFECorrelation::kEventsD0e);

  //This is only called if there are electrons and D0s present.
  int iResult=fCorrelation->Fill(pSelectedD0s.get(), pSelectedElectrons.get());

  if (iResult<0) {
    AliError(Form("%s processing failed with error %d", fCorrelation->GetName(), iResult));
  }

  PostData(1, fOutput);
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
    AliFatal("failed to get output container");
    return;
  }

}
