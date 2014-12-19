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
 * Re-structured analysis task of the pt analysis on EMCal-triggered events:
 * Analysis steps are moved to analysis components, which are grouped by a common
 * event selection. The analysis task steers the event builder, runs each group,
 * and collects the output of all groups.
 *
 *   Author: Markus Fasel
 */
#include "AliInputEventHandler.h"
#include "AliParticleContainer.h"
#include "AliJetContainer.h"

#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerBinningFactory.h"
#include "AliEMCalTriggerEventData.h"
#include "AliEMCalTriggerTaskGroup.h"
#include "AliEMCalTriggerAnaTriggerDecision.h"
#include "AliEMCalHistoContainer.h"
#include "AliAnalysisTaskPtEMCalTriggerV1.h"

ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskPtEMCalTriggerV1)

namespace EMCalTriggerPtAnalysis {

//______________________________________________________________________________
AliAnalysisTaskPtEMCalTriggerV1::AliAnalysisTaskPtEMCalTriggerV1() :
    AliAnalysisTaskEmcalJet(),
    fTaskGroups(NULL),
    fBinning(NULL),
    fMCJetContainer(),
    fDataJetContainer(),
    fSwapTriggerThresholds(kFALSE)
{
  /*
   * Dummy constructor
   */
}

//______________________________________________________________________________
AliAnalysisTaskPtEMCalTriggerV1::AliAnalysisTaskPtEMCalTriggerV1(const char* name) :
    AliAnalysisTaskEmcalJet(name, kTRUE),
    fTaskGroups(NULL),
    fBinning(NULL),
    fMCJetContainer(),
    fDataJetContainer(),
    fSwapTriggerThresholds(kFALSE)
{
  /*
   * Main Constructor
   */
  fTaskGroups = new TObjArray;
  fTaskGroups->SetOwner();
  fBinning = new AliEMCalTriggerBinningComponent();
  SetMakeGeneralHistograms(kTRUE);
}

//______________________________________________________________________________
AliAnalysisTaskPtEMCalTriggerV1::~AliAnalysisTaskPtEMCalTriggerV1() {
  /*
   * Destructor
   */
  delete fTaskGroups;
  delete fBinning;
}

//______________________________________________________________________________
void AliAnalysisTaskPtEMCalTriggerV1::UserCreateOutputObjects() {
  /*
   * Initialise all analysis components
   */
  AliAnalysisTaskEmcal::UserCreateOutputObjects();
  SetCaloTriggerPatchInfoName("EmcalTriggers");

  AliEMCalTriggerBinningFactory binmaker;
  binmaker.Create(fBinning);

  TIter groupIter(fTaskGroups);
  AliEMCalTriggerTaskGroup *mygroup(NULL);
  TList *outputList = new TList;
  outputList->SetName(Form("histos%s", GetName()));
  while((mygroup = dynamic_cast<AliEMCalTriggerTaskGroup *>(groupIter()))){
    mygroup->SetGlobalBinning(fBinning);
    TList *ltmp = mygroup->InitialiseAnalysisComponents();
    // Collect output list and append it to the global output list
    TIter listIter(ltmp);
    TObject *hist(NULL);
    while((hist = listIter())) outputList->Add(hist);
  }
  fOutput->Add(outputList);
  PostData(1, fOutput);

}

//______________________________________________________________________________
Bool_t AliAnalysisTaskPtEMCalTriggerV1::Run() {
  /*
   * Run the analysis:
   * 1st build the event data shared among the tasks
   * 2nd loop over task groups and run them
   */
  AliEMCalTriggerEventData *event = BuildEvent();
  AliEMCalTriggerAnaTriggerDecision triggerDecision;
  if(fSwapTriggerThresholds) triggerDecision.SetSwapThresholds();
  triggerDecision.Create(event);
  triggerDecision.SetIsMinBias(fInputHandler->IsEventSelected() & AliVEvent::kINT7);
  TIter groupIter(fTaskGroups);
  AliEMCalTriggerTaskGroup *mygroup(NULL);
  while((mygroup = dynamic_cast<AliEMCalTriggerTaskGroup *>(groupIter()))){
    mygroup->SetTriggerDecision(&triggerDecision);
    mygroup->Process(event);
  }

  delete event;

  PostData(1, fOutput);
  return kTRUE;
}

//______________________________________________________________________________
void AliAnalysisTaskPtEMCalTriggerV1::SetBinning(const char* dimname, int nbins, double* binning) {
  /*
   * Set binning for a give dimension
   *
   * @param dimname: name of the axis
   * @param nbins: number of bins
   * @param binning: the bin limits
   */
  fBinning->SetBinning(dimname, nbins, binning);
}

//______________________________________________________________________________
void AliAnalysisTaskPtEMCalTriggerV1::SetBinning(const char* dimname, const TArrayD &binning) {
  /*
   * Set binning for a give dimension
   *
   * @param binning: the bin limits
   */
  fBinning->SetBinning(dimname, binning);
}

//______________________________________________________________________________
AliEMCalTriggerEventData* AliAnalysisTaskPtEMCalTriggerV1::BuildEvent() const {
  /*
   * Build event structure. Take the information about the different containers
   * from the base analysis task.
   *
   * @return: the resulting event structure
   */
  AliEMCalTriggerEventData *eventstruct = new AliEMCalTriggerEventData;
  eventstruct->SetRecEvent(fInputEvent);
  eventstruct->SetMCEvent(fMCEvent);
  eventstruct->SetTriggerPatchContainer(fTriggerPatchInfo);
  eventstruct->SetClusterContainer(fCaloClusters);
  eventstruct->SetTrackContainer(fTracks);
  if(fMCJetContainer.Length()){
    AliJetContainer *jcmc = dynamic_cast<AliJetContainer *>(fJetCollArray.FindObject(fMCJetContainer.Data()));
    eventstruct->SetParticleContainer(jcmc->GetParticleContainer()->GetArray());
    eventstruct->SetMCJetContainer(jcmc);
  }
  if(fDataJetContainer.Length()){
    AliJetContainer *jcdat = dynamic_cast<AliJetContainer *>(fJetCollArray.FindObject(fDataJetContainer.Data()));
    eventstruct->SetDataJetContainer(jcdat);
  }
  return eventstruct;
}

//______________________________________________________________________________
void AliAnalysisTaskPtEMCalTriggerV1::AddAnalysisGroup(AliEMCalTriggerTaskGroup *taskGroup) {
  /*
   * Add group of analysis components to the task
   *
   * @param taskGroup: Group of analysis components
   */
  fTaskGroups->Add(taskGroup);
}

} /* namespace EMCalTriggerPtAnalysis */
