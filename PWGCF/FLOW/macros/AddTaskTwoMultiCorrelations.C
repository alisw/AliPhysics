AliAnalysisTaskTwoMultiCorrelations* AddTaskTwoMultiCorrelations(
    const char* trigger = "MB",
    const char* weightFile = "",
    const char* mainTask = "LHC10h-SCs-0070",
    const char* suffix = "")
{
// Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    Error("AddTaskTwoMultiCorrelations()", "No analysis manager to connect to.");
    return NULL;
  } // End: if (!mgr).

// Get the input event handler, again via a static method. 
// This handler is part of the managing system and feeds events to your task.
  if (!mgr->GetInputEventHandler())
  {
    Error("AddTaskTwoMultiCorrelations()", "No input event handler found.");
    return NULL;
  }

// Define the input container.
  AliAnalysisDataContainer *cInput = mgr->GetCommonInputContainer();

// ====================== Default configuration of the analysis task. ====================== //
  TString taskName;
  taskName.Form("%s_%s", mainTask, suffix);
  AliAnalysisTaskTwoMultiCorrelations* myTask = new AliAnalysisTaskTwoMultiCorrelations(taskName,
    kFALSE);
  if (!myTask) {return 0x0;}

// Select minimum-bias events from AliMultSelection in runTwoMultiCorrelations.C
  if (trigger == "MB") {myTask->SelectCollisionCandidates(AliVEvent::kMB);}
// Analysis parameters: Full analysis (= Not kine), Get SCs (= Not ACs), Full writing
//    (= Not minimal output), AC(1,2,3) 
  myTask->SetAnalysisParameters(kFALSE, kFALSE, kFALSE, 1, 2, 3);
// Centrality: 0-70%, 8 bins, 1st estimator: CL1, 2nd estimrator: V0M, No correlation TH2F.
  myTask->SetCentrality(1, 0., 70., 8, "CL1", "V0M", kFALSE);

// PV_x: No selection.
  myTask->SetPVxSelection(kFALSE, -44., 44.);
// PV_y: No selection.
  myTask->SetPVySelection(kFALSE, -44., 44.);
// |PV_z| < 10cm.
  myTask->SetPVzSelection(kTRUE, -10., 10.);
// Multiplicity: 6 tracks/event, 1st FB: TPConly, 2nd FB: global hybrids,
//    Default cuts for HMOS for TPConly, No correlations TH2F.
  myTask->SetHMOsSelection(6, 128, 256, kTRUE, 1.54, -65., 2.3, 90., kFALSE);

// pT in [0.2, 5.0] GeV/c.
  myTask->SetPtSelection(kTRUE, 0.2, 5.);
// eta < |0.8|.
  myTask->SetEtaSelection(kTRUE, -0.8, 0.8);
// Number of TPC clusters > 70.
  myTask->SetNTPCSelection(kTRUE, 70);
// chi^2/NDF in [0.1, 4.0], with personal selection (1).
  myTask->SetChiSelection(1, kTRUE, 0.1, 4.);
// Number of ITS clusters: no selection (we use TPConly).
  myTask->SetNITSSelection(kFALSE, 2);
// DCAxy < 2.4cm.
  myTask->SetDCAxySelection(kTRUE, 2.4);
// DCAz < 3.2cm.
  myTask->SetDCAzSelection(kTRUE, 3.2);
// Electric charge: no charge.
  myTask->SetChargeSelection(kFALSE, kFALSE);

// Keep both primaries and secondaries (not only primaries) at kine (MC only).
  myTask->SetKineSpecifics(kFALSE);
// Without table; with particle weight; with pT; no Phi, no Eta.
  myTask->SetParticleWeights(kFALSE, kTRUE, kTRUE, kFALSE, kFALSE);
// Data period for the list of runs.
  myTask->SetListOfRuns("LHC10h");
// Path to the weight.root file given as argument of the macros.
  TString inputWeightFile;
  inputWeightFile.Form("%s", weightFile);
  myTask->SetInputParticleWeights(inputWeightFile);
// Don't use JEfficiency; 0: index for TPConly.
  myTask->SetJWeights(kFALSE, 0);

// Default power of the reduced Q-vector = 8; with eta gaps
  myTask->SetReducedQvectors(8, 8, 0, kTRUE);
// Default binning of the event histograms.
  myTask->SetBinningEvents(30000, 1000, 20);
// Default binning of the pT, eta and phi distributions.
  myTask->SetBinningTracks(1000, 1000, 720);

// Add the task to the manager.
  mgr->AddTask(myTask);

// ================================ Prepare the containers. ================================ //
// Create the output container.
  char* outputFileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *cOutput = mgr->CreateContainer(taskName,
    TList::Class(), AliAnalysisManager::kOutputContainer,
    Form("%s:%s", outputFileName, "outputAnalysis") );
  /*AliAnalysisDataContainer *cOutput = mgr->CreateContainer(taskName.Data(),
    TList::Class(), AliAnalysisManager::kOutputContainer,
    Form("%s:%s", outputFileName.Data(), "outputAnalysis") );*/

// Connect the manager to the input and output containers.
  mgr->ConnectInput(myTask, 0, cInput);
  mgr->ConnectOutput(myTask, 1, cOutput);

// Return a pointer to the analysis task. Useful when running with the train.
  return myTask;

} // End: AliAnalysisTaskTwoMultiCorrelations *AddTaskTwoMultiCorrelations().
// End of file.
