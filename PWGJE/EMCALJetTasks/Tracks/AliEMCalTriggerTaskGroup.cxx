/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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
/*
 * Group of analysis components with the same event selection
 * Analysis components are initialised via the Initialise function, and executed, if
 * the event is selected, via the function process
 *
 *   Author: Markus Fasel
 */
#include <THashList.h>
#include <TList.h>
#include "AliEMCalTriggerTracksAnalysisComponent.h"
#include "AliEMCalTriggerEventSelection.h"

#include "AliEMCalTriggerTaskGroup.h"

ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup)

namespace EMCalTriggerPtAnalysis {

//______________________________________________________________________________
AliEMCalTriggerTaskGroup::AliEMCalTriggerTaskGroup() :
    TNamed(),
    fAnalysisComponents(NULL),
    fEventSelection(NULL),
    fBinning(NULL),
    fKineCuts(NULL),
    fWeightHandler(NULL)
{
  /*
   * Dummy constructor, not to be used
   */
}

//______________________________________________________________________________
AliEMCalTriggerTaskGroup::AliEMCalTriggerTaskGroup(const char* name) :
    TNamed(name, ""),
    fAnalysisComponents(NULL),
    fEventSelection(NULL),
    fBinning(NULL),
    fKineCuts(NULL),
    fWeightHandler(NULL)
{
  /*
   * Main constructor: to be used by the users
   */
  fAnalysisComponents = new TObjArray();
  fAnalysisComponents->SetOwner();
}

//______________________________________________________________________________
AliEMCalTriggerTaskGroup::~AliEMCalTriggerTaskGroup() {
  /*
   * Destructor
   */
  if(fEventSelection) delete fEventSelection;
  if(fAnalysisComponents) delete fAnalysisComponents;
}

//______________________________________________________________________________
TList *AliEMCalTriggerTaskGroup::InitialiseAnalysisComponents() {
  /*
   * Initialise all analysis components. Build a global histlist for the full group
   *
   * @return: the global histogram list
   */
  TIter compIter(fAnalysisComponents);
  AliEMCalTriggerTracksAnalysisComponent *ana(NULL);
  // Build a global histogram list
  TList *histlist = new TList;
  TObject *htmp(NULL);
  while((ana = dynamic_cast<AliEMCalTriggerTracksAnalysisComponent *>(compIter()))){
    ana->SetBinning(fBinning);
    ana->SetKineCuts(fKineCuts);
    if(fWeightHandler) ana->SetWeightHandler(fWeightHandler);
    ana->CreateHistos();
    TList *ltmp = ana->GetHistList();
    TIter hiter(ltmp);
    while((htmp = hiter())) histlist->Add(htmp);
  }
  return histlist;
}

//______________________________________________________________________________
void AliEMCalTriggerTaskGroup::Process(const AliEMCalTriggerEventData* const event) {
  /*
   * Run analysis of the different groups. Apply event selection if requested;
   *
   * @param event: The combined event data
   */
  if(fEventSelection && !fEventSelection->IsEventSelected(event)) return;
  TIter compIter(fAnalysisComponents);
  AliEMCalTriggerTracksAnalysisComponent *ana(NULL);
  while((ana = dynamic_cast<AliEMCalTriggerTracksAnalysisComponent *>(compIter())))
    ana->Process(event);
}

//______________________________________________________________________________
void AliEMCalTriggerTaskGroup::AddAnalysisComponent(AliEMCalTriggerTracksAnalysisComponent * const analysis) {
  /*
   * Add new analysis component to the task group
   *
   * @param analysis: the analysis component to be added
   */
  fAnalysisComponents->Add(analysis);
}

//______________________________________________________________________________
void EMCalTriggerPtAnalysis::AliEMCalTriggerTaskGroup::SetTriggerDecision(
    const AliEMCalTriggerAnaTriggerDecision* trigger) {
  /*
   * Forward trigger decision to the analysis components
   *
   * @param trigger: the trigger decision
   */
  AliEMCalTriggerTracksAnalysisComponent *myana(NULL);
  TIter compIter(fAnalysisComponents);
  while((myana = dynamic_cast<AliEMCalTriggerTracksAnalysisComponent *>(compIter())))
    myana->SetTriggerDecision(trigger);
}


} /* namespace EMCalTriggerPtAnalysis */
