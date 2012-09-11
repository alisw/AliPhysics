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
#include "AliAnalysisManager.h"
#include "AliAnalysisCuts.h"
#include "AliLog.h"
#include "TH1F.h"
#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliRDHFCutsD0toKpi.h"
#include "TChain.h"
#include "TSystem.h"
#include "TFile.h"
#include "TObjArray.h"
#include <memory>

using namespace std;

/// ROOT macro for the implementation of ROOT specific class methods
ClassImp(AliAnalysisTaskDxHFEParticleSelection)

AliAnalysisTaskDxHFEParticleSelection::AliAnalysisTaskDxHFEParticleSelection(const char* opt)
  : AliAnalysisTaskSE("AliAnalysisTaskDxHFEParticleSelection")
  , fOutput(0)
  , fOption(opt)
  , fCuts(NULL)
  , fSelector(NULL)
  , fUseMC(kFALSE)
  , fFillOnlyD0D0bar(0)
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

  if (fCuts) {
    delete fCuts;
    fCuts=NULL;
  }
}

void AliAnalysisTaskDxHFEParticleSelection::UserCreateOutputObjects()
{
  // create result objects and add to output list

  fOutput = new TList;
  fOutput->SetOwner();

  // setting up for D0s
  TString selectionD0Options;
  switch (fFillOnlyD0D0bar) {
  case 1: selectionD0Options+="FillOnlyD0 "; break;
  case 2: selectionD0Options+="FillOnlyD0bar "; break;
  default: selectionD0Options+="FillD0D0bar ";
  }

  if(fUseMC) fSelector=new AliDxHFEParticleSelectionMCD0(selectionD0Options);
  else fSelector=new AliDxHFEParticleSelectionD0(selectionD0Options);
  fSelector->SetCuts(fCuts);
  fSelector->Init();

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
  AliRDHFCuts* cuts=dynamic_cast<AliRDHFCuts*>(fCuts);
  if (!cuts) return; // Fatal thrown already in initialization

  if(!cuts->IsEventSelected(pEvent)) {
    AliDebug(2,"rejected at IsEventSelected");
    return;
  }

  //  Int_t nInD0toKpi = inputArray->GetEntriesFast();

  fSelector->HistogramEventProperties(AliDxHFEParticleSelection::kEventsSel);
             
  std::auto_ptr<TObjArray> selectedTracks(fSelector->Select(inputArray,pEvent));
  // TODO: use the array of selected track for something, right now
  // only the control histograms of the selection class are filled
  // note: the pointer is deleted automatically once the scope is left
  // if the array should be published, the auto pointer must be released
  // first, however some other cleanup will be necessary in that case
  // probably a clone with a reduced AliVParticle implementation is
  // appropriate.

  if(! selectedTracks.get()) {
    cout << "No selected D0s in this event" << endl; 
    return;
  }

  //Test to see if I have read in D0s and retrieved them after selection
  Int_t nD0Selected = selectedTracks->GetEntriesFast();

  fSelector->HistogramEventProperties(AliDxHFEParticleSelection::kEventsWithParticle);

  for(Int_t iD0toKpi = 0; iD0toKpi < nD0Selected; iD0toKpi++) {
    AliAODRecoDecayHF2Prong *particle = (AliAODRecoDecayHF2Prong*)selectedTracks->UncheckedAt(iD0toKpi);
    if (!particle) continue;
  }

  PostData(1, fOutput);
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
