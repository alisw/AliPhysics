void anaTask()
{
  //Macro to run pi0 calibration on local AOD files.
  //Author: Boris Polishchuk, adapted to AOD by Gustavo Conesa
	
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
  gSystem->Load("libPWG4CaloCalib");
  
  //Analysis settings
 Bool_t bCopy = kTRUE; // True if ESD filter not called before.
 TString eventType = "ESD";

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("list.txt", 15);
  
  // for includes use either global setting in $HOME/.rootrc
  // ACLiC.IncludePaths: -I$ALICE_ROOT/include
  // or in each macro
  // gSystem->AddIncludePath("$ALICE_ROOT/include");
  
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("PHOSPi0Calib");
  
  //Input event handler
  if(eventType =="ESD"){
    AliVEventHandler* esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);
  }
  else if   if(eventType =="AOD"){
    AliVEventHandler* aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);
  }
  
  // ESD filter task
  //gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskESDfilter.C");
  //AliAnalysisTaskESDfilter *esdfilter = AddTaskESDfilter(kFALSE);
  
  // Calibration task 
  //
  AliAnalysisTaskPHOSPi0CalibSelection *task = new AliAnalysisTaskPHOSPi0CalibSelection("PHOSPi0CalibSelection");
  task->SetClusterMinEnergy(0.4); 
  task->CopyAOD(bCopy);
  mgr->AddTask(task);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("histos", TList::Class(),  AliAnalysisManager::kOutputContainer, "PHOShistos.root");
  
  // Connect input/output
  mgr->ConnectInput   (task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput2);
  
  if(bCopy){

    AliAODHandler* aodoutHandler   = new AliAODHandler();
    aodoutHandler->SetOutputFileName("AliAODs_PHOS.root");
    //aodoutHandler->SetCreateNonStandardAOD();
    mgr->SetOutputEventHandler(aodoutHandler);

    AliAnalysisDataContainer *coutput1 =    mgr->GetCommonOutputContainer();
    //mgr->CreateContainer("tree", TTree::Class(),
    //		   AliAnalysisManager::kOutputContainer, "AliAODs_PHOS.root");
    mgr->ConnectOutput (task,     0, coutput1 );

  }


  // Enable debug printouts
  //mgr->SetDebugLevel(10);
  
  if (!mgr->InitAnalysis())
    return;
  
  mgr->PrintStatus();
  
  mgr->StartAnalysis("local", chain);
}
