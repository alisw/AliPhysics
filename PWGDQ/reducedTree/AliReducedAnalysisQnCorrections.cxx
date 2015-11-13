/*
 Creation date: 2015/10/01
 Author: Jacobus Onderwaater, Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
 */


#include "AliReducedAnalysisQnCorrections.h"

#include <iostream>
using std::cout;
using std::endl;

#include <TClonesArray.h>

#include "AliQnCorrectionsManager.h"
#include "AliReducedQnFillEvent.h"

#include "AliReducedVarManager.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedEventInfo.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedCaloClusterInfo.h"
#include "AliReducedPairInfo.h"
#include "AliHistogramManager.h"

ClassImp(AliReducedAnalysisQnCorrections);


#ifdef ALIREDUCEDVARMANAGER_H
#define VAR AliReducedVarManager
#endif

//___________________________________________________________________________
AliReducedAnalysisQnCorrections::AliReducedAnalysisQnCorrections() :
  AliReducedAnalysisTaskSE()
{
  //
  // default constructor
  //
}


//___________________________________________________________________________
AliReducedAnalysisQnCorrections::AliReducedAnalysisQnCorrections(const Char_t* name, const Char_t* title) :
  AliReducedAnalysisTaskSE(name,title),
  fQnFillEvent(0x0),
  fQnManager(0x0)
{
  //
  // named constructor
  //
}


//___________________________________________________________________________
AliReducedAnalysisQnCorrections::~AliReducedAnalysisQnCorrections() 
{
  //
  // destructor
  //
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisQnCorrections::IsEventSelected(AliReducedBaseEvent* event) {
  //
  // apply event cuts
  //
  if(fEventCuts->GetEntries()==0) return kTRUE;
  // loop over all the cuts and make a logical and between all cuts in the list
  // TODO: Cut masks or more complicated cut configurations can be implemented here
  for(Int_t i=0; i<fEventCuts->GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fEventCuts->At(i);
    if(!cut->IsSelected(event)) return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisQnCorrections::IsTrackSelected(AliReducedBaseTrack* track) {
  //
  // apply event cuts
  //
  if(fTrackCuts->GetEntries()==0) return kTRUE;
  // loop over all the cuts and make a logical and between all cuts in the list
  // TODO: Cut masks or more complicated cut configurations can be implemented here
  for(Int_t i=0; i<fTrackCuts->GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fTrackCuts->At(i);
    if(!cut->IsSelected(track)) return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisQnCorrections::IsPairSelected(AliReducedBaseTrack* pair) {
  //
  // apply event cuts
  //
  return kTRUE;
}

//___________________________________________________________________________
void AliReducedAnalysisQnCorrections::Init() {
  //
  // initialize stuff
  //
  fHistosManager = new AliHistogramManager("Histogram Manager", VAR::kNVars);
  fHistosManager->SetUseDefaultVariableNames(kTRUE);
  fHistosManager->SetDefaultVarNames(AliReducedVarManager::fgVariableNames,AliReducedVarManager::fgVariableUnits);

  fQnFillEvent = new AliReducedQnFillEvent();
  fQnFillEvent->SetQnCorrectionsManager(fQnManager);
  fQnFillEvent->SetHistogramManager(fHistosManager);

}


//___________________________________________________________________________
void AliReducedAnalysisQnCorrections::Process() {
  //
  // process the current event
  //  
  if(fEvent->IsA()!=AliReducedEventInfo::Class()) return;

  AliReducedEventInfo* event = (AliReducedEventInfo*)fEvent;
  AliReducedVarManager::SetEvent(event);
  AliReducedVarManager::FillEventInfo(event, fValues);
  fHistosManager->FillHistClass("Event_NoCuts", fValues);
  for(UShort_t ibit=0; ibit<64; ++ibit) {
    AliReducedVarManager::FillEventOnlineTriggers(ibit, fValues);
    fHistosManager->FillHistClass("OnlineTriggers_NoCuts", fValues);
  }
  if(!IsEventSelected(event)) return;
  
    
  fHistosManager->FillHistClass("Event_AfterCuts", fValues);


  fQnManager->ClearEvent();
  fQnFillEvent->Process(event, fValues);
  for(Int_t i=0; i<VAR::kNVars; i++) fQnManager->SetDataContainer(i, fValues[i]);

  fQnManager->SetCalibrationFileDirectoryName(Form("%d",(Int_t) (fValues[AliReducedVarManager::kRunNo]+10E-6)));

  fQnManager->Process();

}


//___________________________________________________________________________
void AliReducedAnalysisQnCorrections::Finish() {
  //
  // run stuff after the event loop
  //
  fQnManager->Finalize();
}




