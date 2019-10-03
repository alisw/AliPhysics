/**
 * \file AddTaskHighPtReducedEventCreator.C
 * \brief Add macro for the reduced event creator task
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Apr 16, 2015
 */
# if !defined __CINT__ || defined __MAKECINT__
#include <TList.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAODInputHandler.h"
#include "AliReducedHighPtEventCreator.h"
#endif

/**
 * Configure reduced event creator task and add it to the analysis manager
 * \param isMC Flag for MC analysis
 * \param isPythia Flag for pythia hard analysis
 * \param trackContainer Name of the track container
 * \param clusterContainer Name of the cluster container
 * \param triggerContainer Name of the trigger container
 * \return Fully configured tree selector task
 */
HighPtTracks::AliReducedHighPtEventCreator *AddTaskHighPtReducedEventCreator(
    Bool_t isMC = kFALSE,
    Bool_t isPythia = kFALSE,
    const char *trackContainer = "",
    const char *clusterContainer = "",
    const char *triggerContainer = "",
    const char *outputfilename = ""
    ){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  // Create task, make event specific settings
  TString taskname = "ReducedEventCreator";
  TString outputname = outputfilename;
  if(outputname.Length()){
    taskname = outputname(0, outputname.Index(".root"));
  }
  HighPtTracks::AliReducedHighPtEventCreator *reducedEventCreator = new HighPtTracks::AliReducedHighPtEventCreator(taskname.Data());
  if(isPythia) reducedEventCreator->SetIsPythia(kTRUE);
  if(isMC) reducedEventCreator->SetSwapTriggerThresholds(kTRUE);
  reducedEventCreator->AddParticleContainer(trackContainer);
  reducedEventCreator->AddClusterContainer(clusterContainer);
  reducedEventCreator->SetCaloTriggerPatchInfoName(triggerContainer);

  // Set kinematic range for tracks and particles
  reducedEventCreator->SetEtaRange(-0.8, 0.8);
  reducedEventCreator->SetPtRange(0.1,100);
  reducedEventCreator->SetClusterEnergyCut(2., 100.);

  // Handle track selection, outsourced to a different macro
  gROOT->LoadMacro(Form("%s/PWGJE/EMCALJetTasks/macros/ReducedEventCutVariation.C", gSystem->Getenv("ALICE_PHYSICS")));
  ReducedEventCutVariation(reducedEventCreator, mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class());

  TString commonoutput = TString(mgr->GetCommonFileName()) + ":" + taskname;
  TString treeoutputfile = "ReducedHighPtEvent.root";
  if(outputname.Length()) treeoutputfile = outputname;
  mgr->ConnectInput(reducedEventCreator, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(reducedEventCreator, 1, mgr->CreateContainer(reducedEventCreator->GetName(), TList::Class(), AliAnalysisManager::kOutputContainer, commonoutput.Data()));
  mgr->ConnectOutput(reducedEventCreator, 2, mgr->CreateContainer("ReducedEvent", TTree::Class(), AliAnalysisManager::kOutputContainer, treeoutputfile.Data()));

  return reducedEventCreator;
}
