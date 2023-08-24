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
#include "AliInputEventHandler.h"
#include "AliParticleContainer.h"
#include "AliJetContainer.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliPicoTrack.h"

#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerBinningFactory.h"
#include "AliEMCalTriggerEventData.h"
#include "AliEMCalTriggerTaskGroup.h"
#include "AliEMCalTriggerAnaTriggerDecision.h"
#include "AliAnalysisTaskPtEMCalTriggerV1.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskPtEMCalTriggerV1)
using namespace PWGJE::EMCALJetTasks;

/**
 * Dummy (I/O) constructor, not to be used.
 */
AliAnalysisTaskPtEMCalTriggerV1::AliAnalysisTaskPtEMCalTriggerV1() :
    AliAnalysisTaskEmcalJet(),
    fTaskGroups(NULL),
    fBinning(NULL),
    fTriggerDecisionConfig(NULL),
    fTriggerClassManager(NULL),
    fMCJetContainer(),
    fDataJetContainer(),
    fSwapTriggerThresholds(kFALSE),
    fDoTriggerDebug(kFALSE)
{
}

/**
 * Main Constructor: Initialises all values with default values. Generating also output container.
 * \param name Name of the component
 */
AliAnalysisTaskPtEMCalTriggerV1::AliAnalysisTaskPtEMCalTriggerV1(const char* name) :
    AliAnalysisTaskEmcalJet(name, kTRUE),
    fTaskGroups(NULL),
    fBinning(NULL),
    fTriggerDecisionConfig(NULL),
    fTriggerClassManager(NULL),
    fMCJetContainer(),
    fDataJetContainer(),
    fSwapTriggerThresholds(kFALSE),
    fDoTriggerDebug(kFALSE)
{
  fTriggerClassManager = new AliEMCalTriggerAnaClassManager("triggermanager");
  fTaskGroups = new TObjArray;
  fTaskGroups->SetOwner();
  fBinning = new AliEMCalTriggerBinningComponent();
  SetMakeGeneralHistograms(kTRUE);
  SetCaloTriggerPatchInfoName("EmcalTriggers");   // Default settings here, to be able to override it in the wagon configuration
}

/**
 * Destructor
 */
AliAnalysisTaskPtEMCalTriggerV1::~AliAnalysisTaskPtEMCalTriggerV1() {
  delete fTriggerClassManager;
  delete fTaskGroups;
  delete fBinning;
}

  /**
   * Initialise all analysis components. Collect histograms from all components
   * and store them in a combined list.
   */
void AliAnalysisTaskPtEMCalTriggerV1::UserCreateOutputObjects() {
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  AliEMCalTriggerBinningFactory binmaker;
  binmaker.Create(fBinning);

  TIter groupIter(fTaskGroups);
  AliEMCalTriggerTaskGroup *mygroup(NULL);
  TList *outputList = new TList;
  outputList->SetName(Form("histos%s", GetName()));
  while((mygroup = dynamic_cast<AliEMCalTriggerTaskGroup *>(groupIter()))){
    mygroup->SetGlobalBinning(fBinning);
    TList *ltmp = mygroup->InitialiseAnalysisComponents(fTriggerClassManager);
    // Collect output list and append it to the global output list
    TIter listIter(ltmp);
    TObject *hist(NULL);
    while((hist = listIter())) outputList->Add(hist);
  }
  fOutput->Add(outputList);
  PostData(1, fOutput);

}

/**
 * Run the analysis:
 *  -# Build the event data shared among the tasks
 *  -# Create the trigger decision and forward it to all components
 *  -# Loop over task groups and execute all components connected to the task group
 */
Bool_t AliAnalysisTaskPtEMCalTriggerV1::Run() {
  AliEMCalTriggerEventData *event = BuildEvent();
  AliEMCalTriggerAnaTriggerDecision triggerDecision;
  if(fDoTriggerDebug) triggerDecision.SetDebugMode();
  if(fTriggerDecisionConfig) triggerDecision.ConfigureTriggerDecision(*fTriggerDecisionConfig);
  triggerDecision.Create(event);
  fTriggerClassManager->SetTriggerDecision(&triggerDecision);
  fTriggerClassManager->PerformEventSelection(event);
  TIter groupIter(fTaskGroups);
  AliEMCalTriggerTaskGroup *mygroup(NULL);
  while((mygroup = dynamic_cast<AliEMCalTriggerTaskGroup *>(groupIter()))){
    mygroup->Process(event);
  }

  delete event;

  PostData(1, fOutput);
  return kTRUE;
}

/**
 * Set binning for a give dimension. Binning is handed over to the binning handler.
 *
 * \param dimname name of the axis
 * \param nbins number of bins
 * \param binning the bin limits
 */
void AliAnalysisTaskPtEMCalTriggerV1::SetBinning(const char* dimname, int nbins, double* binning) {
  fBinning->SetBinning(dimname, nbins, binning);
}

/**
 * Set binning for a give dimension. Binning is handed over to the binning handler.
 *
 * \param binning the bin limits
 */
void AliAnalysisTaskPtEMCalTriggerV1::SetBinning(const char* dimname, const TArrayD &binning) {
  fBinning->SetBinning(dimname, binning);
}

/**
 * Build event structure. Take the information about the different containers
 * from the base analysis task. Also checks whether the track has the event pointer set,
 * and in case not sets it back to this event.
 *
 * \return the resulting event structure
 */
AliEMCalTriggerEventData* AliAnalysisTaskPtEMCalTriggerV1::BuildEvent() {
  AliEMCalTriggerEventData *eventstruct = new AliEMCalTriggerEventData;
  eventstruct->SetRecEvent(fInputEvent);
  // Check whether all tracks have the back pointer to the event itself set
  for(int itrk = 0; itrk < fInputEvent->GetNumberOfTracks(); itrk++){
    FixTrackInputEvent(static_cast<AliVTrack *>(fInputEvent->GetTrack(itrk)));
  }
  eventstruct->SetTriggerBitSelection(fInputHandler->IsEventSelected());
  eventstruct->SetMCEvent(fMCEvent);
  eventstruct->SetTriggerPatchContainer(fTriggerPatchInfo);
  eventstruct->SetClusterContainer(fCaloClusters);
  eventstruct->SetTrackContainer(fTracks);
  for(TIter trackIter = TIter(fTracks).Begin(); trackIter != TIter::End(); ++trackIter)
    FixTrackInputEvent(static_cast<AliVTrack *>(*trackIter));
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

/**
 * Add group of analysis components to the task
 * \param taskGroup Group of analysis components
 */
void AliAnalysisTaskPtEMCalTriggerV1::AddAnalysisGroup(AliEMCalTriggerTaskGroup *taskGroup) {
  fTaskGroups->Add(taskGroup);
}

/**
 * Set the corresponding pointer to the original event to the track in a transparent way for
 * ESD, AOD and pico tracks
 * \param trk The track to handle
 */
void AliAnalysisTaskPtEMCalTriggerV1::FixTrackInputEvent(AliVTrack* trk) {
  if(!trk->GetEvent()){
    if(trk->IsA() == AliESDtrack::Class())
      (static_cast<AliESDtrack *>(trk))->SetESDEvent(static_cast<AliESDEvent *>(fInputEvent));
    else if(trk->IsA() == AliAODTrack::Class())
      (static_cast<AliAODTrack *>(trk))->SetAODEvent(static_cast<AliAODEvent *>(fInputEvent));
    else if(trk->IsA() == AliPicoTrack::Class()){
      AliPicoTrack *mytrk = static_cast<AliPicoTrack *>(trk);
      if(!mytrk->GetEvent()){
        if(mytrk->GetTrack()->IsA() == AliESDtrack::Class())
          (static_cast<AliESDtrack *>(mytrk->GetTrack()))->SetESDEvent(static_cast<AliESDEvent *>(fInputEvent));
        else
          (static_cast<AliAODTrack *>(mytrk->GetTrack()))->SetAODEvent(static_cast<AliAODEvent *>(fInputEvent));
      }
    }
  }
}
