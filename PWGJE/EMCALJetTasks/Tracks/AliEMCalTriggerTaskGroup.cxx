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
#include <THashList.h>
#include <TList.h>
#include "AliEMCalTriggerTracksAnalysisComponent.h"
#include "AliEMCalTriggerEventSelection.h"

#include "AliEMCalTriggerTaskGroup.h"

ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerTaskGroup)

using namespace PWGJE::EMCALJetTasks;

/**
 * Dummy constructor, not to be used
 */
AliEMCalTriggerTaskGroup::AliEMCalTriggerTaskGroup() :
    TNamed(),
    fAnalysisComponents(NULL),
    fEventSelection(NULL),
    fBinning(NULL),
    fKineCuts(NULL),
    fWeightHandler(NULL)
{
}

/**
 * Main constructor: to be used by the users
 */
AliEMCalTriggerTaskGroup::AliEMCalTriggerTaskGroup(const char* name) :
    TNamed(name, ""),
    fAnalysisComponents(NULL),
    fEventSelection(NULL),
    fBinning(NULL),
    fKineCuts(NULL),
    fWeightHandler(NULL)
{
  fAnalysisComponents = new TObjArray();
  fAnalysisComponents->SetOwner();
}

/**
 * Copy constructor. Copy process is delegated to the copy function handling both copy constructor
 * and assignment operator.
 * \param ref Reference for the copy
 */
AliEMCalTriggerTaskGroup::AliEMCalTriggerTaskGroup(const AliEMCalTriggerTaskGroup& ref):
    TNamed(ref),
    fAnalysisComponents(NULL),
    fEventSelection(NULL),
    fBinning(NULL),
    fKineCuts(NULL),
    fWeightHandler(NULL)
{
  ref.Copy(*this);
}

/**
 * Copy constructor. Copy process is delegated to the copy function handling both copy constructor
 * and assignment operator.
 * \param ref Reference for the copy
 */
AliEMCalTriggerTaskGroup& AliEMCalTriggerTaskGroup::operator=(const AliEMCalTriggerTaskGroup& ref) {
  TNamed::operator=(ref);
  if(this != &ref){
    ref.Copy(*this);
  }
  return *this;
}

/**
 * Destructor, cleaning up objects belonging to the task group (analysis components and event selection)
 */
AliEMCalTriggerTaskGroup::~AliEMCalTriggerTaskGroup() {
  if(fEventSelection) delete fEventSelection;
  if(fAnalysisComponents) delete fAnalysisComponents;
}

/**
 * Initialise all analysis components. Build a global histlist for the full group
 *
 * \param classmgr The trigger class manager
 * \return: the global histogram list
 */
TList *AliEMCalTriggerTaskGroup::InitialiseAnalysisComponents(const AliEMCalTriggerAnaClassManager *classmgr) {
  TIter compIter(fAnalysisComponents);
  AliEMCalTriggerTracksAnalysisComponent *ana(NULL);
  // Build a global histogram list
  TList *histlist = new TList;
  TObject *htmp(NULL);
  while((ana = dynamic_cast<AliEMCalTriggerTracksAnalysisComponent *>(compIter()))){
    ana->SetBinning(fBinning);
    ana->SetKineCuts(fKineCuts);
    ana->SetTriggerClassManager(classmgr);
    if(fWeightHandler) ana->SetWeightHandler(fWeightHandler);
    ana->CreateHistos();
    TList *ltmp = ana->GetHistList();
    TIter hiter(ltmp);
    while((htmp = hiter())) histlist->Add(htmp);
  }
  return histlist;
}

/**
 * Run analysis of the different groups. Apply event selection if requested;
 *
 * \param event The combined event data
 */
void AliEMCalTriggerTaskGroup::Process(const AliEMCalTriggerEventData* const event) {
  if(fEventSelection && !fEventSelection->IsEventSelected(event)) return;
  TIter compIter(fAnalysisComponents);
  AliEMCalTriggerTracksAnalysisComponent *ana(NULL);
  while((ana = dynamic_cast<AliEMCalTriggerTracksAnalysisComponent *>(compIter())))
    ana->Process(event);
}

/**
 * Add new analysis component to the task group
 *
 * \param analysis: the analysis component to be added
 */
void AliEMCalTriggerTaskGroup::AddAnalysisComponent(AliEMCalTriggerTracksAnalysisComponent * const analysis) {
  if(fWeightHandler) analysis->SetWeightHandler(fWeightHandler);
  fAnalysisComponents->Add(analysis);
}

/**
 * Set the weight handler to the group. If components are registered, publish the weight handler
 * to them
 * \param handler The weight handler for the group
 */
void AliEMCalTriggerTaskGroup::SetWeightHandler(const AliEMCalTriggerWeightHandler* handler) {
  fWeightHandler = handler;
  for(TIter compIter = TIter(fAnalysisComponents).Begin(); compIter != TIter::End(); ++compIter){
    AliEMCalTriggerTracksAnalysisComponent *mycomp = static_cast<AliEMCalTriggerTracksAnalysisComponent *>(*compIter);
    mycomp->SetWeightHandler(fWeightHandler);
  }
}

/**
 * Copy function. Does a flat copy of all elements which are not under control of the task group. Will
 * make a flat copy of all analysis components as well for the moment.
 * \param other Target where the information from this object is set to.
 */
void AliEMCalTriggerTaskGroup::Copy(TObject& other) const {
  AliEMCalTriggerTaskGroup &target = static_cast<AliEMCalTriggerTaskGroup &>(other);
  target.fKineCuts = fKineCuts;
  target.fBinning = fBinning;
  target.fWeightHandler = fWeightHandler;
  target.fEventSelection = new AliEMCalTriggerEventSelection(*fEventSelection);
  target.fAnalysisComponents = new TObjArray();
  target.fAnalysisComponents->SetOwner(false);
  for(TIter compIter = TIter(fAnalysisComponents).Begin(); compIter != TIter::End(); ++compIter){
    target.fAnalysisComponents->Add(*compIter);
  }
}
