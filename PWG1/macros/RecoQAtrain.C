// RS: This is macro to build an analysis train which will run from the reco loop
//     Uses special AliRecoInputHandler 
// Should be attached in rec.C as
// AliReconstruction rec; rec.SetAnalysisMacro(thisMacroPath);

void LoadLibraries();
void AddAnalysisTasks(); 


void RecoQAtrain()
{
  //
  TGrid::Connect("alien://");
  if (!gGrid || !gGrid->IsConnected()) {
    ::Error("QAtrain", "No grid connection");
    return;
  }   
  // Set temporary merging directory to current one
  gSystem->Setenv("TMPDIR", gSystem->pwd());
  // Set temporary compilation directory to current one
  gSystem->SetBuildDir(gSystem->pwd(), kTRUE);
  // Load libraries
  LoadLibraries();
  printf("Include path: %s\n", gSystem->GetIncludePath());
  // Create manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis", "Reco QA train");
  // Input handler
  AliRecoInputHandler *esdHandler = new AliRecoInputHandler();
  esdHandler->SetReadFriends(kTRUE);
  esdHandler->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdHandler);
  //  
  // AnalysisTasks
  AddAnalysisTasks();
  //
  if (!mgr->InitAnalysis()) printf("Failed to initialize Reco Analysis");
}

void LoadLibraries()
{
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/PWG1 -I$ALICE_ROOT/PWG1/TRD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libTENDER");
  gSystem->Load("libPWG0base.so");
  gSystem->Load("libPWG0dep.so");
  gSystem->Load("libPWG0selectors.so");
  gSystem->Load("libPWG1.so");
}

void AddAnalysisTasks()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  mgr->SetCommonFileName("RecoQAresults.root");
  //
  // Event Statistics (Jan Fiete)
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kFALSE /*MC*/);
  //
  // add standard ITSAlignQA task with only SDD calibration histos activated
  gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskITSAlign.C");
  AliAnalysisTaskITSAlignQA *itsAlign = AddTaskITSAlign(0,2011,kTRUE);
  itsAlign->SetLoadGeometryFromOCDB(kFALSE);
  itsAlign->SetDoSPDResiduals(kFALSE);
  itsAlign->SetDoSDDResiduals(kFALSE);
  itsAlign->SetDoSSDResiduals(kFALSE);
  itsAlign->SetDoSDDDriftTime(kFALSE);
  itsAlign->SetMinMaxMult(20.,1070.);
  //
}
