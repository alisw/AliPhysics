void AnalyzeDiff()
{
    
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  
  //load analysis framework
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice"); //AliAnalysisTaskSE
  
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT/PHOS");

  // A task can be compiled dynamically with AClic
  gROOT->LoadMacro("AliCaloPhoton.cxx+g");
  gROOT->LoadMacro("AliAnalysisTaskPi0DiffEfficiency.cxx+g");
  
  // Create the chain
  TChain* chain = new TChain("aodTree");
  chain->AddFile("AliAODout.root") ;

  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("Pi0EmbeddingManager");
  
  // ESD input handler
  AliAODInputHandler* aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);

  // Output
  
  // Debug level
  mgr->SetDebugLevel(0);

  AliAnalysisTaskPi0DiffEfficiency * task = new AliAnalysisTaskPi0DiffEfficiency("Pi0Efficiency") ;
  
  TFile *fBadMap = TFile::Open("BadMap_LHC10h.root");
  if(fBadMap->IsOpen()){
    printf("\n\n...Adding PHOS bad channel map \n") ;
    gROOT->cd();
    char key[55] ;
    for(Int_t mod=1;mod<4; mod++){
      sprintf(key,"PHOS_BadMap_mod%d",mod) ;
      TH2I * h = (TH2I*)fBadMap->Get(key) ;
      if(h)
	task->SetPHOSBadMap(mod,h) ;
    }
    fBadMap->Close() ;
  }


  mgr->AddTask(task);


  // Create containers for input/output
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer(); 
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("histESD",TList::Class(),AliAnalysisManager::kOutputContainer,"histosDiff.root");
  
  // Connect input/output
  mgr->ConnectInput(task , 0, cinput);
  mgr->ConnectOutput(task, 1,coutput);

  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
  }
  
}
