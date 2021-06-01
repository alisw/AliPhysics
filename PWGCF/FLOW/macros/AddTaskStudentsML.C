AliAnalysisTaskStudentsML* AddTaskStudentsML(
    const char* trigger = "MB",
    const char* mainTask = "LHC10h-SPCs-0080",
    const char* WeightFile = "",
    const char* suffix = "")
{
// Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    Error("AddTaskStudentsML()", "No analysis manager to connect to.");
    return NULL;
  } // End: if (!mgr).

// Get the input event handler, again via a static method. 
// This handler is part of the managing system and feeds events to your task.
  if (!mgr->GetInputEventHandler())
  {
    Error("AddTaskStudentsML()", "No input event handler found.");
    return NULL;
  }

// Define the input container.
  AliAnalysisDataContainer *cInput = mgr->GetCommonInputContainer();

// ====================== Default configuration of the analysis task. ====================== //
  TString taskName;
  taskName.Form("%s_%s", mainTask, suffix);
  AliAnalysisTaskStudentsML* myTask = new AliAnalysisTaskStudentsML(taskName,
    kFALSE);
  if (!myTask) {return 0x0;}

// Select minimum-bias events from AliMultSelection in runStudentsML.C
  if (trigger == "MB") {myTask->SelectCollisionCandidates(AliVEvent::kMB);}

// Analysis parameters: Physics (SetDoAnalysis: kTRUE), Minimal Writing (SetSaveAllQA: kTRUE to save all QA plots)
  myTask->SetDoAnalysis(kTRUE); 
  myTask->SetSaveAllQA(kTRUE); 

// Centrality: Set edges of the centrality bin, negative value for breaking point. Ascending order important. First two values have to be valid centrality values
  myTask->SetCentrality(0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,-10.,-10.,-10.,-10.,-10.,-10.,-10.);
  myTask->SetInitializeCentralityArray();
  myTask->SetCentralityEstimator("CL1");

// PV_x: No selection.
  myTask->SetVertexX(kFALSE,-44.,-44.);
// PV_y: No selection.
  myTask->SetVertexY(kFALSE,-44.,-44.);
// |PV_z| < 10cm.
  myTask->SetVertexZ(kTRUE,-10.,10.);
// Multiplicity: 14 tracks/event (SetMinNuPar), 1st FB: TPConly, 2nd FB: global hybrids,
//    Default cuts for HMOS for TPConly
  myTask->SetMinNuPar(14.);
  myTask->SetFilter(128);
  myTask->SetBoolMultCut(kTRUE, 256);
  myTask->SetUpperLineCut(2.3, 90.);
  myTask->SetLowerLineCut(1.54,-65.);

// pT in [0.2, 5.0] GeV/c.
  myTask->SetPtCut(kTRUE,0.2,5.0);
// eta < |0.8|.
  myTask->SetEtaCut(kTRUE,-0.8,0.8);
// Number of TPC clusters > 70.
  myTask->SetNumberTPCClusters(kTRUE, 70.);
// chi^2/NDF in [0.1, 4.0]. Chi^2 selection method chosen by the second parameter (here 1)
  myTask->SetChiSquareTPC(kTRUE, 1, 0.1, 4.0);
// Number of ITS clusters: no selection (we use TPConly).
  myTask->SetNumberITSClusters(kFALSE, 2.);
// DCAxy < 2.4cm.
  myTask->SetDCAxy(kTRUE, 2.4);
// DCAz < 3.2cm.
  myTask->SetDCAz(kTRUE, 3.2);
// Charge Cut. First parameter: if kTRUE: Cut on charge. Second parameter: if kTRUE: Select only positive.
  myTask->SetChargeCut(kFALSE,kTRUE);
//Fisher Yates Randomizing. First Parameter: Do FY or not. Second Parameter: How much of the originally passed tracks should be kept
  myTask->SetFisherYates(kFALSE, 1.); 
// Use (or not use) Reco-Kine-Table (1. parameter to kTRUE), keep weak secondaries (2. parameter kTRUE)
  myTask->SetKineSpecifics(kFALSE,kTRUE);
// Use weights or not (1. parameter), use Phi (2. parameter), pt (3. parameter), eta (4. parameter)
  myTask->SetUseWeights(kFALSE,kFALSE,kTRUE,kFALSE);
// Data period for the list of runs.
  myTask->SetListOfRuns("LHC10h");
// Path to the weight.root file.

  TString InputWeightFile;
  InputWeightFile.Form("%s",WeightFile);

  myTask->SetInputParticleWeights(InputWeightFile);

// Configure Multiparticle Correlators
  myTask->SetCorrSet1(4., 2., -2., 4, -4., 0., 0.,0.);
  myTask->SetCorrSet2(2., 2., -2., 0., 0., 0., 0.,0.);
  myTask->SetCorrSet3(2., 4., -4., 0., 0., 0., 0.,0.);
  myTask->SetCorrSet4(4., 2., -2., 3.,-3., 0., 0.,0.);
  myTask->SetCorrSet5(2., 3., -3., 0., 0., 0., 0.,0.);
  myTask->SetCorrSet6(0., 0., 0., 0., 0., 0., 0.,0.);
  myTask->SetCorrSet7(0., 0., 0., 0., 0., 0., 0.,0.);
  myTask->SetCorrSet8(0., 0., 0., 0., 0., 0., 0.,0.);
  myTask->SetMixed(kFALSE,2., kFALSE, kTRUE);


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

} // End: AliAnalysisTaskStudentsML *AddTaskStudentsML().
// End of file.
