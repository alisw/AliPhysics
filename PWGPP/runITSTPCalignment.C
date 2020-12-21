void runITSTPCalignment( UInt_t saveinterval = 1000000 )
{
  TStopwatch timer;
  timer.Start();

  //runProof("/COMMON/COMMON/LHC09a4_run8100X#esdTree");
  runLocal("AliESDs.root");

  timer.Stop();
  timer.Print();
}

//_________________________________________________//
void runLocal(const char* filenamestr = "AliESDs.root" ) {

  TString inputFilename(filenamestr);
  TString outputArrayFilename = "ITSTPCalignmentArray.root";
  TString outputHistFilename = "ITSTPCalignmentHist.root";

  //____________________________________________________//
  //_____________Setting up the par files_______________//
  //____________________________________________________//
  setupPar("STEERBase");
  gSystem->Load("libSTEERBase");
  setupPar("ESD");
  gSystem->Load("libVMC");
  gSystem->Load("libESD");
  setupPar("AOD");
  gSystem->Load("libAOD");
  setupPar("ANALYSIS");
  gSystem->Load("libANALYSIS");
  setupPar("ANALYSISalice");
  gSystem->Load("libANALYSISalice");
  //____________________________________________________//  

  //add input files from dirs
  TChain* chain = new TChain("esdTree");
  chain->SetBranchStatus("*ESDfriend*",0);
  TString workingDir = gSystem->pwd();
  void* dirhandle = gSystem->OpenDirectory(workingDir.Data());
  if (!dirhandle) return;
  const char* filenamestr;
  while ((filenamestr = gSystem->GetDirEntry(dirhandle)))
  {
    TString filename(filenamestr);
    if (filename=="." || filename=="..") continue;
    if (gSystem->cd(filename.Data()))
    {
      if (!gSystem->AccessPathName(inputFilename.Data()))//return value is here inverted!
      {
        TString inputESDfile(workingDir);
        inputESDfile += "/";
        inputESDfile += filename;
        inputESDfile += "/";
        inputESDfile += inputFilename;
        chain->Add(inputESDfile.Data());
        printf("found file: %s\n", inputESDfile.Data());
      }
      gSystem->cd(workingDir.Data());
    }
  }

  //____________________________________________//
  gROOT->LoadMacro("AliRelAlignerKalman.cxx++");
  gROOT->LoadMacro("AliRelAlignerKalmanArray.cxx++");
  gROOT->LoadMacro("AliAnalysisTaskITSTPCalignment.cxx++");

  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("ITSTPCalignmentAnalysisManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);

  //create the task
  AliAnalysisTaskITSTPCalignment *taskITSTPCalignment = 
                        new AliAnalysisTaskITSTPCalignment("TaskITSTPCalignment");
  taskITSTPCalignment->SetDoQA(kTRUE);

  mgr->AddTask(taskITSTPCalignment);

  AliAnalysisDataContainer* cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer* coutput0 = mgr->CreateContainer("outputArray",
                                       AliRelAlignerKalmanArray::Class(),
                                       AliAnalysisManager::kOutputContainer,
                                       outputArrayFilename.Data());
  AliAnalysisDataContainer* coutput1 = mgr->CreateContainer("outputList",
                                       TList::Class(),
                                       AliAnalysisManager::kOutputContainer,
                                       outputHistFilename.Data());

  mgr->ConnectInput(taskITSTPCalignment,0,cinput0);
  mgr->ConnectOutput(taskITSTPCalignment,0,coutput0);
  mgr->ConnectOutput(taskITSTPCalignment,1,coutput1);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();

  mgr->StartAnalysis("local",chain);
}

//_________________________________________________//
void runInteractive(const char* collectionName = "tag.xml") {
  gSystem->Load("libProofPlayer");

  TString outputArrayFilename = "ITSTPCalignmentArray.root";
  TString outputHistFilename = "ITSTPCalignmentHist.root";

  printf("*** Connect to AliEn ***\n");
  TGrid::Connect("alien://");
 
  //____________________________________________________//
  //_____________Setting up the par files_______________//
  //____________________________________________________//
  setupPar("STEERBase");
  gSystem->Load("libSTEERBase");
  setupPar("ESD");
  gSystem->Load("libVMC");
  gSystem->Load("libESD");
  setupPar("AOD");
  gSystem->Load("libAOD");
  setupPar("ANALYSIS");
  gSystem->Load("libANALYSIS");
  setupPar("ANALYSISalice");
  gSystem->Load("libANALYSISalice");
  //____________________________________________________//  
  
  //____________________________________________//
  AliTagAnalysis *tagAnalysis = new AliTagAnalysis("ESD");
 
  AliRunTagCuts *runCuts = new AliRunTagCuts();
  AliLHCTagCuts *lhcCuts = new AliLHCTagCuts();
  AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
  AliEventTagCuts *evCuts = new AliEventTagCuts();
 
  //grid tags
  TGridCollection* coll = gGrid->OpenCollection(collectionName);
  TGridResult* TagResult = coll->GetGridResult("",0,0);
  tagAnalysis->ChainGridTags(TagResult);
  TChain* chain = 0x0;
  chain = tagAnalysis->QueryTags(runCuts,lhcCuts,detCuts,evCuts);
  chain->SetBranchStatus("*Calo*",0);
  
  //____________________________________________//
  gROOT->LoadMacro("AliRelAlignerKalman.cxx++");
  gROOT->LoadMacro("AliRelAlignerKalmanArray.cxx++");
  gROOT->LoadMacro("AliAnalysisTaskITSTPCalignment.cxx++");

  // analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("ITSTPCalignmentAnalysisManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);

  //create the task
  AliAnalysisTaskITSTPCalignment *taskITSTPCalignment = 
                        new AliAnalysisTaskITSTPCalignment("TaskITSTPCalignment");
  taskITSTPCalignment->SetDoQA(kTRUE);

  mgr->AddTask(taskITSTPCalignment);

  AliAnalysisDataContainer* cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer* coutput0 = mgr->CreateContainer("outputArray",
                                       AliRelAlignerKalmanArray::Class(),
                                       AliAnalysisManager::kOutputContainer,
                                       outputArrayFilename.Data());
  AliAnalysisDataContainer* coutput1 = mgr->CreateContainer("outputList",
                                       TList::Class(),
                                       AliAnalysisManager::kOutputContainer,
                                       outputHistFilename.Data());

  mgr->ConnectInput(taskITSTPCalignment,0,cinput0);
  mgr->ConnectOutput(taskITSTPCalignment,0,coutput0);
  mgr->ConnectOutput(taskITSTPCalignment,1,coutput1);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();

  mgr->StartAnalysis("local",chain);
}

//______________________________________________________________________________
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
