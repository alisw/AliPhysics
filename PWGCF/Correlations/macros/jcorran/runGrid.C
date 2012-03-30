const Bool_t kMC = kFALSE; //With real data kMC = kFALSE
const TString kInputData = "AOD";
const TString kJCORRANInputFormat = "AOD"; // ESD, AODout, AODin

//_____________________________________________________________________
void runGrid(){
  // Load Custom Configuration and parameters
  // override values with parameters

  //==== Load common libraries
  gSystem->Load("libCore.so");  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libXMLIO.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS"); 
  gSystem->Load("libANALYSISalice");   
  gSystem->Load("libEMCALUtils");
  gSystem->Load("libPHOSUtils");
  gSystem->Load("libPWG4PartCorrBase");
  gSystem->Load("libPWG4PartCorrDep");
  gSystem->Load("libPWGCFJCORRAN");

  //==== Load Ours
// SetupPar("PWG4JCORRAN"); 

  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/EMCAL");
	gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/ANALYSIS ");

  // Create and configure the alien handler plugin
  gROOT->LoadMacro("CreateAlienHandler.C");
  AliAnalysisGrid *alienHandler = CreateAlienHandler();  
  if (!alienHandler) return;
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("JCHII");

  // Connect plug-in to the analysis manager
  mgr->SetGridHandler(alienHandler);

  // MC handler
  if(kMC){
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mcHandler->SetReadTR(kFALSE);//Do not search TrackRef file
    mgr->SetMCtruthEventHandler(mcHandler);
  }
 
  if(kInputData == "ESD"){
    // ESD handler
    AliESDInputHandler *esdHandler = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdHandler);
  }
  if(kInputData == "AOD"){
    // AOD handler
    AliAODInputHandler *aodHandler = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodHandler);
  }

  //================================================
  // TASKS
  //================================================

  //==== PID
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AddTaskPIDResponse(kMC);  
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
  AddTaskPIDqa();
 
  //==== Statistics
  mgr->AddStatisticsTask();

  //==== CENTRALITY
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentrality = AddTaskCentrality(); 
  //taskCentrality->SetPass(2);

  //==== Physics Selection
  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(kMC, kTRUE);
  physSelTask->SelectCollisionCandidates(AliVEvent::kMB);

  //============================
  //   JCORRANTask
  //===========================
  //==== Basic Track Cut
  AliESDtrackCuts* fEsdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  //http://aliceinfo.cern.ch/static/aliroot-new/html/roothtml/src/AliESDtrackCuts.cxx.html#DY4BnB

  //==== Additional Track Filter ;
  gROOT->ProcessLine(".L ./AddESDFilter.C");
  AliAnalysisFilter *filter = AddESDFilter();

  //==== JCORRAN TASK
  AliJCORRANTask *jctask = new AliJCORRANTask("PWG4JCORRANTask",kJCORRANInputFormat);
  jctask->SetESDtrackCuts(fEsdTrackCuts);
  jctask->SetRealOrMC(kMC);
  jctask->SetOutputAODName("jcorran.root");
  jctask->SetDebugLevel(0);
  jctask->SetESDFilter( filter );
  jctask->SetRunType("LHC11h");
  jctask->SetStoreEventPlaneSource(true);
  jctask->SetStoreTPCTrack(true);

  //==event selection
  jctask->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kHighMult);  //Apply offline trigger selection by AliPhysicsSelectionTask

  mgr->AddTask((AliAnalysisTask*) jctask);

  //==== Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAODHandler* aodoutHandler   = new AliAODHandler();
  aodoutHandler->SetCreateNonStandardAOD();
  mgr->SetOutputEventHandler(aodoutHandler);

  // Connect input/output
  AliAnalysisDataContainer *runinfoOutput = mgr->CreateContainer("RunInfo",  TList::Class(), AliAnalysisManager::kOutputContainer, "jcorran.root");
  mgr->ConnectInput(jctask, 0, cinput);
  mgr->ConnectOutput(jctask, 1, runinfoOutput );


  //===============================
  // START MANAGER
  //===============================

  // Enable debug printouts
  mgr->SetDebugLevel(0);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();
  // Start analysis in grid.
  printf("********** LD_LIBRARY_PATH *******\n");
  printf("%s\n", gSystem->Getenv("LD_LIBRARY_PATH"));
  printf("**********************************\n");


//  printf("%s\n", gSystem->Getenv("LD_LIBRARY_PATH"));

  mgr->StartAnalysis("grid");
};




//========================================================================
// Setup Par Files
//========================================================================
void SetupPar(char* pararchivename)
{
  //Load par files, create analysis libraries
  //For testing, if par file already decompressed and modified
  //classes then do not decompress.

  TString cdir(Form("%s", gSystem->WorkingDirectory() )) ;
  TString parpar(Form("%s.par", pararchivename)) ;
  // create par file if it not exist
  if ( gSystem->AccessPathName(parpar.Data()) ) {
    gSystem->ChangeDirectory(gSystem->Getenv("ALICE_ROOT")) ;
    TString processline(Form(".! make %s", parpar.Data())) ;
    gROOT->ProcessLine(processline.Data()) ;
    gSystem->ChangeDirectory(cdir) ;
    processline = Form(".! mv /tmp/%s .", parpar.Data()) ;
    gROOT->ProcessLine(processline.Data()) ;
  }
  // decompres par file
  if ( gSystem->AccessPathName(pararchivename) ) {
    TString processline = Form(".! tar xvzf %s",parpar.Data()) ;
    gROOT->ProcessLine(processline.Data());
  }

  TString ocwd = gSystem->WorkingDirectory();
  gSystem->ChangeDirectory(pararchivename);

  // check for BUILD.sh and execute
  if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
    printf("*******************************\n");
    printf("*** Building PAR archive    ***\n");
    cout<<pararchivename<<endl;
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
    cout<<pararchivename<<endl;
    printf("*******************************\n");
    gROOT->Macro("PROOF-INF/SETUP.C");
  }

  gSystem->ChangeDirectory(ocwd.Data());
  printf("Current dir: %s\n", ocwd.Data());

}


void LoadConf( TString filename ){
  ifstream in(filename.Data());
  TString str;
  str.ReadFile(in);
  str.ReplaceAll("\n",  "");
  cout<<str<<endl;

  gROOT->ProcessLine( str.Data() );
}

