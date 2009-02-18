void runBatch() {
  TStopwatch timer;
  timer.Start();

  printf("*** Connect to AliEn ***\n");
  TGrid::Connect("alien://");
  gSystem->Load("libProofPlayer.so");

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

  //CreateXML();
  
  //ANALYSIS PART
  gROOT->LoadMacro("AliAnalysisTaskPt.cxx+");
  const char *collectionfile = "wn.xml";

  //____________________________________________//
  //Usage of event tags
  AliTagAnalysis *analysis = new AliTagAnalysis();
  TChain *chain = 0x0;
  chain = analysis->GetChainFromCollection(collectionfile,"esdTree");
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

void CreateXML() {
  // Create A tag analysis object and impose some selection criteria
  AliTagAnalysis *tagAna = new AliTagAnalysis(); 

  //Case where the tag files are stored locally
  //TagAna->ChainLocalTags(".");

  //Case where the tag files are stored in the file catalog
  //pp.xml is the xml collection of tag files that was produced 
  //by querying the file catalog.
  TAlienCollection* coll = TAlienCollection::Open("tag.xml");
  TGridResult* tagResult = coll->GetGridResult("",0,0);
  tagAna->ChainGridTags(tagResult);

  //__________________________//
  //Usage of string statements//
  //__________________________//
  /*const char* fRunCuts = "fAliceRunId == 340";
  const char* fEventCuts = "(fEventTag.fTopPtMin >= 1.0)&&(fEventTag.fNumberOfTracks >= 11)&&(fEven
tTag.fNumberOfTracks <= 12)";
  tagAna->CreateXMLCollection("global",fRunCuts,fEventCuts);*/

  //________________________________________________//
  //          Usage of Ali*TagCuts classes          //
  //________________________________________________//
  // create a RunTagCut object
  AliRunTagCuts *runCuts = new AliRunTagCuts();
  // create a LHCTagCut object
  AliLHCTagCuts *lhcCuts = new AliLHCTagCuts();
  // create a DetectorTagCut object
  AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
  // create an EventTagCut object
  AliEventTagCuts *evCuts = new AliEventTagCuts();
  evCuts->SetMultiplicityRange(11,12);
  tagAna->CreateXMLCollection("global",runCuts,lhcCuts,detCuts,evCuts);
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
        Error("runProcess","Cannot Build the PAR Archive! - Abort!");
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
