void anaTask()
{
  //Macro to run pi0 calibration on local ESD files.
  //Author: Boris Polishchuk.
  

  //Uncomment the line below if your ESD files are from old productions
  //so there are no PHOS geometry matrices written in the ESD.
  //You can find an ideal geometry here: $ALICE_ROOT/test/QA/geometry.root
  //and copy to your working directory.

  //AliGeomManager::LoadGeometry("geometry.root");

  //You can apply misalignment by following
  //(remember that your local OCDB in ./PHOS should contain a copy of the
  //$ALICE_ROOT/OCDB/PHOS/Align directory while ./PHOS/Calib/EmcGainPedestals
  //should contain the corrections to the calibration coefficients (~1) 
  //instead of real CC (~0.005)! ):

//   AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
//   AliCDBManager::Instance()->SetSpecificStorage("PHOS/*","local://./");

//   AliCDBManager::Instance()->SetRun(0);
//   AliCDBEntry* e = AliCDBManager::Instance()->Get("PHOS/Align/Data");
//   TClonesArray *array = (TClonesArray*) e->GetObject();
//   AliGeomManager::ApplyAlignObjsToGeom(*array);


  // load analysis framework
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice"); //AliAnalysisTaskSE
  gSystem->Load("libPHOSpi0Calib");

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("list.txt", 300);

  // for includes use either global setting in $HOME/.rootrc
  // ACLiC.IncludePaths: -I$ALICE_ROOT/include
  // or in each macro
  // gSystem->AddIncludePath("$ALICE_ROOT/include");

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("Pi0Calib");



  AliVEventHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  // Create task

  AliAnalysisTaskPi0CalibSelection *task = new AliAnalysisTaskPi0CalibSelection("Pi0CalibSelection");
  //task->SetClusterMinEnergy(0.4); 

  // Add task
  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("histos", TList::Class(),  AliAnalysisManager::kOutputContainer, "histos.root");

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

  // Enable debug printouts
  mgr->SetDebugLevel(10);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();

  mgr->StartAnalysis("local", chain);
}
