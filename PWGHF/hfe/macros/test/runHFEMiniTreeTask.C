//****************************************
//Nirbhay K. Behera (nbehera@cern.ch)
//***************************************
//READEME
//This is the runGrid macro to run the HFE minitree task-
//Following parameters setting must be checked with care before running

//Task BASE Name---
TString ctaskname = "MiniTree";//User choice 

//Set these parameters before running job ---
//1) ESD or AOD : isAOD-----------------------
Bool_t  isAOD     = 1; //0 for ESD, 1 for AOD

//2) MC or data : isMC-----------------------
Bool_t  isMC      = 1; //1 for MC events, 0 for data

//3) Collision system - "pp" or "AA" (PbPb, pPb)------
TString collisionSystem = "pp"; //change the collision system

//4) Grid or local---------------(for local option, set the list file, for grid Check AliPhysics version) 
const char *RunType = "grid"; //"local" or "grid"

//5)Test run or full job submission or terminate--
const char *RunMode = "test"; //"full", "terminate"

//6)Merge via JDL---
Bool_t mergeOption = 1;// mergeViaJDL = 1; 0 for no

//7) Optional task name--
const char *TaskName = "Test"; //user choice

//==================================================
//Task specific cuts----change with care-------
Double_t TPCchi2 = 4.;
Int_t MinTPCNcluster = 100;
Int_t MinTPCclusterPID = 80;
Double_t TPCclusterRatio = 0.6;
Int_t MinNclusterITS = 3;
Bool_t checkITSLayerstatus = kFALSE;
Double_t eta = 0.8;
Double_t ptMin = 0.5;
Double_t ptMax = 100.;
Double_t Vz = 10.;
Double_t dcaxy = 1.;
Double_t dcaz = 2.;
Double_t prodVz = 0.5;
Double_t spdResolution = 0.25;
Double_t nsigmaTPClow = -1.;
Double_t nsigmaTPChigh = 3.;
Double_t nsigmaTOF = 3.;

//===================================================

//______________________________________________________________________________
void runHFEMiniTreeTask(  const char *runtype     = RunType, 
			  const char *gridmode    = RunMode,
			  Bool_t      mergeviajdl = mergeOption,
			  const char *taskname    = TaskName) {
  
  TStopwatch gtimer;
  gtimer.Start() ;
  
  ctaskname += taskname;
  
  cout << "Output Dir = : " << ctaskname.Data() << endl;
  
  
  if(runtype != "local" && runtype != "proof" && runtype != "grid"){
    Printf("\n\tIncorrect run option, check first argument of run macro");
    Printf("\tint runtype = local, proof or grid\n");
    return;
  }
  Printf("%s analysis chosen",runtype);
  
  
  gSystem->Load("libCore");
  gSystem->Load("libEG");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libGui");
  gSystem->Load("libMinuit");
  gSystem->Load("libPhysics");
  gSystem->Load("libProof");
  gSystem->Load("libMLP");
  gSystem->Load("libTreePlayer");
  gSystem->Load("libXMLIO");
  gSystem->Load("libXMLParser");
  gSystem->Load("libXrdClient");
  
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libANALYSISaliceBase");
  gSystem->Load("libOADB");
  gSystem->Load("libpythia6");
  gSystem->Load("libAliPythia6");
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");
  gSystem->Load("libTender");
  gSystem->Load("libTOFbase");
  gSystem->Load("libPWGTRD");
  gSystem->Load("libPWGHFbase");
  gSystem->Load("libPWGHFhfe");
  
  
  gROOT->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_PHYSICS")));
  gROOT->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_ROOT")));
  
  printf("Include path: %s\n\n", gSystem->GetIncludePath());
  
  
  AliAnalysisManager* mgr = new AliAnalysisManager(taskname);
  
  if( runtype == "grid"){
    AliAnalysisGrid *plugin = CreateAlienHandler(taskname, gridmode, isAOD, mergeviajdl); 
    mgr->SetGridHandler(plugin);
    gROOT->LoadMacro("AliHFEminiTrack.cxx+g");
    gROOT->LoadMacro("AliHFEminiEvent.cxx+g");
    gROOT->LoadMacro("AliHFEminiEventCreator.cxx+g");
  }
  else if(runtype == "local"){
    gROOT->LoadMacro("AliHFEminiTrack.cxx+g");
    gROOT->LoadMacro("AliHFEminiEvent.cxx+g");
    gROOT->LoadMacro("AliHFEminiEventCreator.cxx+g");
    if(isAOD){
      TChain* chain = new TChain("aodTree");
      ifstream file_collect("PbPbaodData.list");
    }
    else{
      TChain* chain = new TChain("esdTree");
      ifstream file_collect("mcesdfile.list");
    }
    
    
    TString line;
    while (line.ReadLine(file_collect) ) {
      chain->Add(line.Data());
    }
    
    chain->GetListOfFiles()->Print();
    
  }
  
  if(isAOD) {
    AliVEventHandler* aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);
  }
  else {
    AliVEventHandler* esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);
    if(isMC) { // galice.root Kinematics.root AliESDs.root  TrackRefs.root AliESDfriends.root
      AliMCEventHandler *mc = new AliMCEventHandler();
      mc->SetReadTR(kFALSE);
      mgr->SetMCtruthEventHandler(mc);
    } 
  }
  
  if(collisionSystem != "pp"){
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AddTaskMultSelection(kFALSE);
  }
  
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  
  if(isMC){
    AliPhysicsSelectionTask* physicsSelTask = AddTaskPhysicsSelection(isMC);
  }
  else AliPhysicsSelectionTask* physicsSelTask = AddTaskPhysicsSelection(isMC,kTRUE);  
  
  
  //PID response
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  
  if(isMC) AddTaskPIDResponse(kTRUE);   
  else AddTaskPIDResponse();

  //Not implemented now
  /*
  if(runImproverTask){  
    gROOT->LoadMacro("AddTaskImproveITS.C"); 
    AddTaskImproveITS();
  }
  */
  
  gROOT->LoadMacro("AddTaskHFEminiEventCreator.C");
  AddTaskHFEminiEventCreator(isMC, TPCchi2, MinTPCNcluster, MinTPCclusterPID, TPCclusterRatio, MinNclusterITS, checkITSLayerstatus, eta, ptMin, ptMax, Vz, dcaxy, dcaz, prodVz, spdResolution, nsigmaTPClow, nsigmaTPChigh, nsigmaTOF, collisionSystem, TaskName);

  //AddTaskHFEminiEventCreator(4., 100, 80, 0.6, 3, kFALSE, 0.8, 0.5, 100., 10., 1., 2., 0.5, 0.25, -1, 3., 3., "pp", "TestRun");
  
  // mgr->Dump();
  mgr->SetDebugLevel(1);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  Printf("Starting Analysis....");
  if(runtype == "grid") mgr->StartAnalysis(runtype);
  else mgr->StartAnalysis("local",chain);
  
  cout << "---------------------------------------------------------" << endl << endl;
  gtimer.Stop();
  gtimer.Print();
  cout << "---------------------------------------------------------" << endl << endl;
  
}
//______________________________________________________________________________
//______________________________________________________________________________
AliAnalysisGrid* CreateAlienHandler(const char *taskname, const char *gridmode, Bool_t isAOD, Bool_t mergeviajdl) {
    

  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetOverwriteMode();
  plugin->SetRunMode(gridmode);
  
  plugin->SetAPIVersion("V1.1x");
  plugin->SetAliPhysicsVersion("vAN-20170801-1");

  //for test purpose only--
  
  if(isAOD){//---------AOD data set---------------------------------------------------------->>>
    if(isMC){
      //AOD MC-----
      if(collisionSystem =="pp"){
	//for pp MC
	plugin->SetGridDataDir("/alice/sim/2017/LHC17c3b2/");
	plugin->SetDataPattern("/AOD/*/AliAOD.root");
	Int_t runNum[1] = { 259888 };
	for(Int_t i = 0; i < 1; i++) plugin->AddRunNumber(runNum[i]);
      }
      else{//AOD MC PbPb-----
	plugin->SetGridDataDir("/alice/sim/2016/LHC16i3a/");
	plugin->SetDataPattern("/AOD/*/AliAOD.root");
	Int_t runNum[1] = { 246994 };
	for(Int_t i = 0; i < 1; i++) plugin->AddRunNumber(runNum[i]);
      }
    }
    else {
      //AOD Data-----------
      if(collisionSystem =="pp"){
	//for AOD pp data
	plugin->SetGridDataDir("/alice/data/2016/LHC16l/");
	plugin->SetDataPattern("/pass1/AOD/*/AliAOD.root");
	plugin->SetRunPrefix("000");
	Int_t runNum[1] = { 259860 };
	for(Int_t i = 0; i < 1; i++) plugin->AddRunNumber(runNum[i]);
      }
      else{
	//AOD PbPb data--- 
	plugin->SetGridDataDir("/alice/data/2015/LHC15o/");
	plugin->SetDataPattern("/pass1_pidfix/AOD186/*/AliAOD.root");
	plugin->SetRunPrefix("000");
	Int_t runNum[1] = { 245554 };
	for(Int_t i = 0; i < 1; i++) plugin->AddRunNumber(runNum[i]);
      }
    }
  }
  else{//--------Set the ESD dataset------------------------------------------------------------->>>>
    
    if(isMC){
      if(collisionSystem =="pp"){//ESD MC pp--
	plugin->SetGridDataDir("/alice/sim/2017/LHC17c3b2/");
	plugin->SetDataPattern("*AliESDs.root");
	Int_t runNum[1] = { 259888 };
	for(Int_t i = 0; i < 1; i++) plugin->AddRunNumber(runNum[i]);
      }
      else{//ESD MC PbPb----
	plugin->SetGridDataDir("/alice/sim/2016/LHC16i3a/");
	plugin->SetDataPattern("*AliESDs.root");
	Int_t runNum[1] = { 246994 };
	for(Int_t i = 0; i < 1; i++) plugin->AddRunNumber(runNum[i]);
      }
    }
    else{//ESD Data-------
      if(collisionSystem =="pp"){
	//for ESD pp data
	plugin->SetGridDataDir("/alice/data/2016/LHC16l/");
	plugin->SetDataPattern("/pass1/*AliESDs.root");
	plugin->SetRunPrefix("000");
	Int_t runNum[1] = { 259860 };
	for(Int_t i = 0; i < 1; i++) plugin->AddRunNumber(runNum[i]);
      }
      else{
	//ESD PbPb data---- 
	plugin->SetGridDataDir("/alice/data/2015/LHC15o/");
	plugin->SetDataPattern("/pass1_pidfix/*AliESDs.root");
	plugin->SetRunPrefix("000");
	Int_t runNum[1] = { 245554 };
	for(Int_t i = 0; i < 1; i++) plugin->AddRunNumber(runNum[i]);
	
      }	  
    }
  }//else ESD finish-
  
  
  plugin->SetNtestFiles(1);
  
  //---------------------------------------
  plugin->SetGridWorkingDir(ctaskname.Data());
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  
  plugin->SetOutputToRunNo();  
  plugin->SetAnalysisSource("AliHFEminiTrack.cxx AliHFEminiEvent.cxx AliHFEminiEventCreator.cxx");
  plugin->SetAdditionalLibs("libSTEERBase.so libOADB.so AliHFEminiTrack.cxx AliHFEminiEvent.cxx AliHFEminiEventCreator.cxx AliHFEminiTrack.h AliHFEminiEvent.h AliHFEminiEventCreator.h");
  
  
  TString includes_str = "-Wno-deprecated -I$. -I$CGAL_DIR/include -I$FASTJET/include -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include";
  plugin->AddIncludePath(includes_str.Data()); // for grid running
  
  //  plugin->SetCheckCopy(kFALSE);
  plugin->SetDefaultOutputs(kFALSE);
  plugin->SetOutputFiles("AnalysisResults.root");
  plugin->SetExecutable(Form("%s.sh",taskname));
  plugin->SetAnalysisMacro(Form("%s.C",taskname));
  plugin->SetTTL(30000);
  plugin->SetInputFormat("xml-single");
  plugin->SetJDLName(Form("%s.jdl",taskname));
  plugin->SetPrice(1);
  plugin->SetSplitMode("se");
  plugin->SetMergeViaJDL(mergeviajdl);
  if(isAOD) plugin->SetSplitMaxInputFileNumber(50);
  else plugin->SetSplitMaxInputFileNumber(100);
  plugin->SetNrunsPerMaster(1);
  plugin->SetMaxMergeFiles(30);
  // Give the local files
  plugin->SetFileForTestMode("file.txt"); 
  //  else plugin->SetFileForTestMode("esd.txt"); 
  
  // Proof cluster
  plugin->SetProofCluster(" ");
  // Dataset to be used   
  //  plugin->SetProofDataSet("/default/sjena/test");
  // May need to reset proof. Supported modes: 0-no reset, 1-soft, 2-hard
  plugin->SetProofReset(0);
  // May limit number of workers
  // plugin->SetNproofWorkers(0);
  // May limit the number of workers per slave
  //  plugin->SetNproofWorkersPerSlave(1);
  // May use a specific version of root installed in proof
  //  plugin->SetRootVersionForProof("current");
  // May set the aliroot mode. Check http://aaf.cern.ch/node/83 
  //  plugin->SetAliRootMode("default"); // Loads AF libs by default
  // May request ClearPackages (individual ClearPackage not supported)
  plugin->SetClearPackages(kFALSE);
  // Plugin test mode works only providing a file containing test file locations, used in "local" mode also
  //  plugin->SetFileForTestMode("files.txt"); // file should contain path name to a local directory containg *ESDs.root etc
  // Request connection to alien upon connection to grid
  plugin->SetProofConnectGrid(kFALSE);
  // Other PROOF specific parameters
  plugin->SetProofParameter("PROOF_UseMergers","-1");
  printf("Using: PROOF_UseMergers   : %s\n", plugin->GetProofParameter("PROOF_UseMergers"));
  return plugin;
}



