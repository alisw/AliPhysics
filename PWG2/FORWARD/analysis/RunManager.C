void RunManager()
{
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  
  gSystem->Load("libSTEERBase");
  
  gSystem->Load("libESD") ;
  gSystem->Load("libAOD") ;
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  
  gSystem->Load("libPhysics");
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG0dep");
  gSystem->Load("libPWG2forward");
  
  //You can expand this chain if you have more data :-)
  TChain* chain = new TChain("esdTree");
  
  chain->Add("AliESDs.root");
  
  //Creating the manager and handlers
  AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", "A test setup for the analysis train");
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  
    
  esdHandler->SetInactiveBranches("AliESDACORDE AliRawDataErrorLogs CaloClusters Cascades EMCALCells EMCALTrigger ESDfriend Kinks Kinks Cascades ALIESDACORDE MuonTracks TrdTracks CaloClusters");

       
  mgr->SetInputEventHandler(esdHandler);      
       
       
  // Monte Carlo handler
    
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);
  
  //mcHandler->SetReadTR(readTR);    
  
  // AOD output handler
  AliAODHandler* aodHandler   = new AliAODHandler();
  mgr->SetOutputEventHandler(aodHandler);
  aodHandler->SetOutputFileName("AliAODs.root");
    
           
  
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis/AddTaskFMD.C");
  AliFMDAnalysisTaskSE *fmdtask = AddTaskFMD();
  
  // Run the analysis
    
  TStopwatch t;
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    t.Start();
    mgr->StartAnalysis("local", chain);
    t.Stop();
    t.Print();
  }  
}
