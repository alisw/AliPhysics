//=============================================================================
//
// *** AddTRDdigitsFilter
//
// This macro sets up the TRD digits filter to write the digits for a
// PID reference sample of electrons and pions from V0 decays.
//
//=============================================================================

#if ! defined (__CINT__) || defined (__MAKECINT__)
#include <AliAnalysisManager.h>
#include <AliAnalysisDataContainer.h>
#include <AliAnalysisTask.h>
#include <AliESDv0KineCuts.h>
#include <AliTRDdigitsFilter.h>
#endif

AliAnalysisTask  *AddTRDdigitsFilter(TString cfg)
{
  //gSystem->Load("libTRDrec");
  // pointer to the analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTRDdigitsFilter", "No analysis manager to connect to.");
    return NULL;
  }

  // check the input handler
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask", "This task requires an input event handler");
    return NULL;
  }

  /////////////////////////
  // The TRD filter Task
  /////////////////////////
  AliTRDdigitsFilter *filterTask = new AliTRDdigitsFilter();

  //if (runNumber >= 295274 && runNumber <= 29762) {
  if ( cfg == "PbPb-2018" ) {

    // LHC18q,r  -  Pb-Pb 2018

    filterTask->GetV0KineCuts()->SetMode(AliESDv0KineCuts::kPurity,
                                     AliESDv0KineCuts::kPbPb);

    filterTask->AcceptTracks("v0elec", AliTRDdigitsFilter::kPidV0Electron,
                                2.0, 99999999., 1.0);

    filterTask->AcceptTracks("v0pilo", AliTRDdigitsFilter::kPidV0Pion,
                                2.0, 2.5, 0.3);

    filterTask->AcceptTracks("v0pihi", AliTRDdigitsFilter::kPidV0Pion,
                                2.5, 99999999., 1.0);

    filterTask->AcceptTracks("v0prot", AliTRDdigitsFilter::kPidV0Proton,
                                2.0, 99999999., 1.0);

    filterTask->AcceptEvents("cent", 0.0, 2.0, 1.0e-2);
  }

  cout << "filter task set up" << endl;
  filterTask->PrintSettings();

  mgr->AddTask(filterTask);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  if (!cinput) cinput = mgr->CreateContainer("cchain",TChain::Class(),
                                      AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *coutput =mgr->CreateContainer("TRDdigitsFilter",TList::Class(), AliAnalysisManager::kOutputContainer, "DigitsFilter.root");


  mgr->ConnectInput(filterTask,0,cinput);
  mgr->ConnectOutput(filterTask,1,coutput);
  return filterTask;

}
