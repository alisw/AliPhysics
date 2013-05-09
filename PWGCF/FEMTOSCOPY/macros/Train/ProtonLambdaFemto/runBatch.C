//_________________________________
// Femto analysis
//_________________________________

void runBatch() {
  TStopwatch timer;
  timer.Start();


//______ connect to Alien
  printf("*** Connect to AliEn ***\n");
  TGrid::Connect("alien://");


  gSystem->Load("libProofPlayer.so");
  gSystem->Load("libVMC.so");


//_______ Load libraries
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");

   Bool_t useParFiles   = kFALSE;  
   Bool_t useTenderPars = kFALSE;
   Bool_t usePWGCFParFiles = kFALSE;
   
   Bool_t useMC = kFALSE;

   TString format = "aod";
   //format = "aod";
   format.ToLower();

//______ Dynamic Path ______
   if (usePWGCFParFiles) {
   TString dynpath(".:");
   dynpath += gSystem->GetDynamicPath();
   gSystem->SetDynamicPath(dynpath);
    //char dynpath[10000];
    //sprintf(dynpath, ".:%s", gSystem->GetDynamicPath());
    //gSystem->SetDynamicPath(dynpath);
  }

  // Use precompiled libraries for the analysis framework
  if (useParFiles)
    setupPar("STEERBase");
	  if (gSystem->Load("libSTEERBase.so")<0) {
    cout << "Cannot load libSTEERBase.so . Exiting" << endl;
    exit(0);
  }
  gSystem->Load("libVMC.so");

  if (useParFiles)
    setupPar("ESD");
  if (gSystem->Load("libESD.so")<0) {
    cout << "Cannot load libESD.so . Exiting" << endl;
    exit(0);
  }

  if (useParFiles)
    setupPar("AOD");
  if (gSystem->Load("libAOD.so")<0) {
    cout << "Cannot load libAOD.so . Exiting" << endl;
    exit(0);
  }

  if (useParFiles)
    setupPar("ANALYSIS");
  if (gSystem->Load("libANALYSIS.so")<0) {
    cout << "Cannot load libANALYSIS.so . Exiting" << endl;
    exit(0);
  }

  if (useParFiles)
    setupPar("ANALYSISalice");
  if (gSystem->Load("libANALYSISalice.so")<0) {
    cout << "Cannot load libANALYSISalice.so . Exiting" << endl;
    exit(0);
  }
  

  //____________________________________________________//
  //_____________Setting up PWGCFAOD.par_________________//

  // if (usePWGCFParFiles)
  //   setupPar("PWGCFAOD");
  // if (gSystem->Load("libPWGCFAOD.so")<0) {
  //   cout << "Cannot load libPWGCFAOD.so . Exiting" << endl;
  //   exit(0);
  // }
  //____________________________________________________//
  //_____________Setting up PWGCFfemtoscopy.par__________//
  if (usePWGCFParFiles)
    setupPar("PWGCFfemtoscopy");
  if (gSystem->Load("libPWGCFfemtoscopy.so")<0) {
    cout << "Cannot load libPWGCFfemtoscopy.so . Exiting" << endl;
    exit(0);
  }
  //____________________________________________________//
  //_____________Setting up PWGCFfemtoscopyUser.par______//
    if (usePWGCFParFiles)
    setupPar("PWGCFfemtoscopyUser");
  if (gSystem->Load("libPWGCFfemtoscopyUser.so")<0) {
    cout << "Cannot load libPWGCFfemtoscopyUser.so . Exiting" << endl;
    exit(0);
    }



   cout <<"_____GetDynamicPath______\n " <<gSystem->GetDynamicPath() <<endl;  


//____________ include path
   gSystem->AddIncludePath(Form("-I\"%s/include\"", gSystem->Getenv("ALICE_ROOT")));
   gROOT->ProcessLine(Form(".include %s/include", gSystem->ExpandPathName("$ALICE_ROOT")));




  
  //_________________________________________________ 
  //_______Create chain for Alien data collection __________
  
   const char *collectionfile="wn.xml";
  
  //____________________________________________//
  //Usage of event tags
  // AliTagAnalysis *analysis = new AliTagAnalysis();
  // TChain *chain = 0x0;
  // chain = analysis->CreateChainFromCollection(collectionfile,"esdTree");

  // gROOT->LoadMacro("./CreateESDChain.C");
  // const char* chainlistfile = "./list.ESD.txt";
  // chain = CreateESDChain(chainlistfile,500);

  TChain *chain = new TChain("aodTree");
  //chain->Add("../Data/LHC10h/AOD/139465/AliAOD.root");
  ifstream *istr = new ifstream(collectionfile);

  char fname[2000];
  char pname[2000];
  while (!istr->eof()) {
    fname[0] = '\0';
    (*istr) >> fname;
    if (strlen(fname) > 10) {
      sprintf(pname, "alien://%s", fname);
      chain->Add(pname);
    }
    }
  

/*
//__________________________________________________
//___________Create chain for Local data files ____________

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("files.txt", 2);
*/



  //___________ Analysis  __________________________

  //______ Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");

  //______ Set ESD InputHandler
  AliAODInputHandler* aodH = new AliAODInputHandler;
  mgr->SetInputEventHandler(aodH);

  //AliESDInputHandler* esdH = new AliESDInputHandler;
  //mgr->SetInputEventHandler(esdH);
  
  //______ Set MC EventHandler
  //AliMCEventHandler *mcH = new AliMCEventHandler;
  //mgr->SetMCtruthEventHandler(mcH);

  //______ Set Print Debug Level
  mgr->SetDebugLevel(0);  //0, 1, 2, 3 ...
  
  

  //AddTaskPIDResponse
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskSE *pidresponse = AddTaskPIDResponse(kTRUE,kFALSE);

                                                                           
   //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
   //AliCentralitySelectionTask *centrality = AddTaskCentrality();

  //________AddTaskFemto_______________
    gROOT->LoadMacro("AddTaskFemto.C");
  AliAnalysisTaskFemto *taskfemto = AddTaskFemto("ConfigFemtoAnalysis.C");
  taskfemto->SelectCollisionCandidates(AliVEvent::kCentral|AliVEvent::kSemiCentral|AliVEvent::kMB);
  //taskfemto->SelectCollisionCandidates(AliVEvent::kMB);

  //____________________________________________//
  // Run the analysis
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();
}


//*********************************************
//____________________________________________
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
