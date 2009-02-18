void runInteractive() {
  TStopwatch timer;
  timer.Start();
  gSystem->Load("libProofPlayer.so");

  printf("*** Connect to AliEn ***\n");
  TGrid::Connect("alien://");
 
  //____________________________________________________//
  //_____________Setting up STEERBase.par_______________//
  //____________________________________________________//
  setupPar("STEERBase");
  gSystem->Load("libSTEERBase.so");

  //____________________________________________________//
  //_____________Setting up ESD.par_____________________//
  //____________________________________________________//
  setupPar("ESD");
  gSystem->Load("libVMC.so");
  gSystem->Load("libESD.so");

  //____________________________________________________//
  //_____________Setting up AOD.par_____________________//
  //____________________________________________________//
  setupPar("AOD");
  gSystem->Load("libAOD.so");

  //_____________________________________________________________//
  //_____________Setting up ANALYSIS_NEW.par_____________________//
  //_____________________________________________________________//
  setupPar("ANALYSIS");
  gSystem->Load("libANALYSIS.so");

  gROOT->LoadMacro("AliAnalysisTaskPt.cxx+");
  
  //____________________________________________//
  AliTagAnalysis *TagAna = new AliTagAnalysis("ESD");
 
  AliRunTagCuts *runCuts = new AliRunTagCuts();
  AliLHCTagCuts *lhcCuts = new AliLHCTagCuts();
  AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
  AliEventTagCuts *evCuts = new AliEventTagCuts();
  evCuts->SetMultiplicityRange(11,12);

  //grid tags
  TAlienCollection* coll = TAlienCollection::Open("tag.xml");
  TGridResult* TagResult = coll->GetGridResult("",0,0);
  TagAna->ChainGridTags(TagResult);
  TChain* chain = 0x0;
  chain = TagAna->QueryTags(runCuts,lhcCuts,detCuts,evCuts);
  chain->SetBranchStatus("*Calo*",0);

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);  
  //____________________________________________//
  // 1st Pt task
  AliAnalysisTaskPt *task1 = new AliAnalysisTaskPt("TaskPt");
  mgr->AddTask(task1);
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist1", TH1::Class(),AliAnalysisManager::kOutputContainer,"Pt.ESD.root");
  
  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();
}

Int_t setupPar(const char* pararchivename) {
  ///////////////////
  // Setup PAR File//
  ///////////////////
  if (pararchivename) {
    char processline[1024];
    sprintf(processline,".! tar xvzf %s.par",pararchivename);
    gROOT->ProcessLine(processline);
    const char* ocwd = gSystem->WorkingDirectory();
    gSystem->ChangeDirectory(pararchivename);

    // check for BUILD.sh and execute
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      printf("*******************************\n");
      printf("*** Building PAR archive    ***\n");
      printf("*******************************\n");

      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
        Error("runAnalysis","Cannot Build the PAR Archive! - Abort!");
        return -1;
      }
    }
    // check for SETUP.C and execute
    if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
      printf("*******************************\n");
      printf("*** Setup PAR archive       ***\n");
      printf("*******************************\n");
      gROOT->Macro("PROOF-INF/SETUP.C");
    }
    
    gSystem->ChangeDirectory("../");
  }

  return 1;
}
