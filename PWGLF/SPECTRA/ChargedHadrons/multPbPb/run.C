// TODO:
// 1. Check cuts for 2010 (Jochen?)
// 2. Run with many centrality bins at once
#include <string.h>

enum { kMyRunModeLocal = 0, kMyRunModeGRID = 2};

TList * listToLoad = new TList();

TChain * GetAnalysisChain(const char * incollection);

void run(Char_t* data, Long64_t nev = -1, Long64_t offset = 0, Bool_t debug = kFALSE, Int_t runMode = 0, Bool_t isMC = 0, 
	 Int_t centrBin = 0, const char * centrEstimator = "VOM", Int_t useOtherCentralityCut = 0, Int_t trackMin=0, Int_t trackMax=10000, 
	 const char* option = "",TString customSuffix = "", Int_t workers = -1, Bool_t useSingleBin=kTRUE, const char * runList = 0)
{
  // runMode:
  //
  // 0 local 
  // 1 proof

  if (nev < 0)
    nev = 1234567890;

  InitAndLoadLibs(runMode,workers,debug);

  // Create the analysis manager
  mgr = new AliAnalysisManager;

  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler;
  // Do I need any of this? 
  esdH->SetInactiveBranches("AliESDACORDE FMD ALIESDTZERO ALIESDZDC AliRawDataErrorLogs CaloClusters Cascades EMCALCells EMCALTrigger ESDfriend Kinks AliESDTZERO ALIESDACORDE MuonTracks TrdTracks");
  mgr->SetInputEventHandler(esdH);

  if(isMC) {
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetPreReadMode(AliMCEventHandler::kLmPreRead);
    mgr->SetMCtruthEventHandler(handler);
  }



  // Parse option strings
  TString optionStr(option);
  
  // remove SAVE option if set
  // This  is copied from a macro by Jan. The reason I kept it is that I may want to pass textual options to the new task at some point
  Bool_t doSave = kFALSE;
  TString optionStr(option);
  if (optionStr.Contains("SAVE"))
    {
      optionStr = optionStr(0,optionStr.Index("SAVE")) + optionStr(optionStr.Index("SAVE")+4, optionStr.Length());
      doSave = kTRUE;
    }

  AliESDtrackCuts * cuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE);  
  TString pathsuffix = "";

  if(!useSingleBin) pathsuffix += "_AllCentr";

  if (optionStr.Contains("DCA")) {
    delete cuts;
    //    cuts = AliESDtrackCuts::GetStandardITSPureSATrackCuts2009();
    cout << ">>>> USING DCA cut" << endl;
    cuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();  
    cuts->SetMaxChi2TPCConstrainedGlobal(); // remove golden cut
    pathsuffix+="_DCAcut";
  }

  if (optionStr.Contains("ITSsa")) {
    delete cuts;
    cuts = AliESDtrackCuts::GetStandardITSPureSATrackCuts2009();
    cout << ">>>> USING ITS sa tracks" << endl;
    pathsuffix+="_ITSsa";
  }

  if (optionStr.Contains("TPC")) {
    delete cuts;
    cuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    cout << ">>>> USING TPC only tracks" << endl;
    pathsuffix+="_TPC";
  }

  if(optionStr.Contains("NoElectrons")) pathsuffix +="_NoElectrons";



  Bool_t useMCKinematics = isMC;
  if (optionStr.Contains("NOMCKIN")) {
    cout << ">>>> Ignoring MC kinematics" << endl;
    useMCKinematics=kFALSE;
    pathsuffix+="_NOMCKIN";
  }
  
  // If we are running on grid, we need the alien handler
  if (runMode == kMyRunModeGRID) {
    // Create and configure the alien handler plugin
    gROOT->LoadMacro("CreateAlienHandler.C");
    AliAnalysisGrid *alienHandler = CreateAlienHandler(data, runList, pathsuffix.Data(), listToLoad, "full", isMC);  // full
    if (!alienHandler) {
      cout << "Cannot create alien handler" << endl;    
      exit(1);
    }
    mgr->SetGridHandler(alienHandler);  
  }



  // Add tasks
  // physics selection
  gROOT->ProcessLine(".L $ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
  physicsSelectionTask = AddTaskPhysicsSelection(isMC);

  //PID
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AddTaskPIDResponse(isMC); 

  // // PID QA:
  // gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
  // AddTaskPIDqa();


  // Centrality
  gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentr = AddTaskCentrality();
  // OBSOLETE
  const char * file1 = "$ALICE_ROOT/ANALYSIS/macros/test_AliCentralityBy1D.root";
  const char * file2 = "$ALICE_ROOT/ANALYSIS/macros/test_AliCentralityByFunction.root";
  // END of OBSOLETE
  if(isMC) taskCentr-> SetMCInput();
  taskCentr->SetPass(2);

  // Create my own centrality selector
  AliAnalysisMultPbCentralitySelector * centrSelector = new AliAnalysisMultPbCentralitySelector();
  centrSelector->SetIsMC(isMC);
  centrSelector->SetCentrTaskFiles(file1,file2); // for bookkeping only
  centrSelector->SetCentralityBin(centrBin);
  if (!useSingleBin) centrSelector->SetCentralityBin(0); // FIXME: ok?
  centrSelector->SetCentralityEstimator(centrEstimator);

  if(useOtherCentralityCut == 1){
    cout << "Setting centrality by MULT" << endl;
    centrSelector->SetUseMultRange();
    centrSelector->SetMultRange(trackMin,trackMax);
  }
  if(useOtherCentralityCut == 2){
    cout << "Setting centrality by V0" << endl;
    
    centrSelector->SetUseV0Range();
    centrSelector->SetMultRange(trackMin,trackMax);
  }
  if(useOtherCentralityCut == 3){
    cout << "Setting centrality by SPD outer" << endl;    
    centrSelector->SetUseSPDOuterRange();
    centrSelector->SetMultRange(trackMin,trackMax);
  }

  // load my task
  if (useSingleBin) {
    gROOT->ProcessLine(".L $ALICE_ROOT/PWG0/multPbPb/AddTaskMultPbPbTracks.C");
    AliAnalysisTaskMultPbTracks * task = AddTaskMultPbPbTracks("multPbPbtracks.root", cuts, centrSelector); 
    task->SetIsMC(useMCKinematics);
    task->SetOfflineTrigger(AliVEvent::kMB);
    if(optionStr.Contains("TPC")) task->SetTPCOnly();
    if(optionStr.Contains("NoElectrons")) task->RejectElectrons(kTRUE);
    if(useMCKinematics) task->GetHistoManager()->SetSuffix("MC");
    if(customSuffix!=""){
      cout << "Setting custom suffix: " << customSuffix << endl;    
      task->GetHistoManager()->SetSuffix(customSuffix);
    }
  } else {
    gROOT->ProcessLine(".L $ALICE_ROOT/PWG0/multPbPb/AddTaskMultPbPbTracksAllCentrality.C");
    centrSelector->SetUseV0Range(kTRUE);
    Int_t ncentr = 11;
   
    const Float_t minCentr[] = {0 ,79 ,239,559 ,1165,2135,3555,5525,8213 ,12191,15079};
    const Float_t maxCentr[] = {79,239,559,1165,2135,3555,5525,8213,12191,15079,21000};
    AliAnalysisTaskMultPbTracks ** tasks = AddTaskMultPbPbTracksAllCentrality("multPbPbtracks.root", cuts, centrSelector, ncentr,minCentr,maxCentr); 
    for(Int_t icentr = 0; icentr < ncentr; icentr++){
      tasks[icentr]->Print();
      cout << "MC KINEMATICS:" << useMCKinematics << endl;
      
      tasks[icentr]->SetIsMC(useMCKinematics);
      tasks[icentr]->SetOfflineTrigger(AliVEvent::kMB);
      if(optionStr.Contains("TPC")) tasks[icentr]->SetTPCOnly();
      if(optionStr.Contains("NoElectrons")) task[icentr]->RejectElectrons(kTRUE);
      if(useMCKinematics) tasks[icentr]->GetHistoManager()->SetSuffix("MC");
      if(customSuffix!=""){
	cout << "Setting custom suffix: " << customSuffix+long(icentr) << endl;    
	tasks[icentr]->GetHistoManager()->SetSuffix(customSuffix+long(icentr));
      }	
    }    
  }
  // Init and run the analy
  if (!mgr->InitAnalysis()) return;

  mgr->PrintStatus();
  
  if (runMode == kMyRunModeLocal ) {
    // If running in local mode, create chain of ESD files
    cout << "RUNNING LOCAL, CHAIN" << endl;    
    TChain * chain = GetAnalysisChain(data);
    chain->Print();
    mgr->StartAnalysis("local",chain,nev);
  } else if (runMode == kMyRunModeGRID) {
    mgr->StartAnalysis("grid");
  } else {
    cout << "ERROR: unknown run mode" << endl;        
  }

  if (!useOtherCentralityCut) {
    pathsuffix = pathsuffix + "_" + centrEstimator + "_bin_"+long(centrBin);
  } else if(useOtherCentralityCut==1){
    pathsuffix = pathsuffix + "_TrackRange_" + long(trackMin) + "_" + long(trackMax);
  } else if(useOtherCentralityCut==2){
    pathsuffix = pathsuffix + "_V0Range_" + long(trackMin) + "_" + long(trackMax);
  } else if(useOtherCentralityCut==3){
    pathsuffix = pathsuffix + "_SPDOutRange_" + long(trackMin) + "_" + long(trackMax);
  }
  pathsuffix += customSuffix;

  if (doSave) MoveOutput(data, pathsuffix.Data());

  

  
}


void MoveOutput(const char * data, const char * suffix = ""){

  TString path("output/");
  path = path + TString(data).Tokenize("/")->Last()->GetName() + suffix;
  
  TString fileName = "multPbPbtracks.root";
  gSystem->mkdir(path, kTRUE);
  gSystem->Rename(fileName, path + "/" + fileName);
  for(Int_t ibin = 0; ibin < 20; ibin++){
    TString fileBin = fileName;
    fileBin.ReplaceAll(".root",Form("_%2.2d.root",ibin));
    gSystem->Rename(fileBin, path + "/" + fileBin);    
  }
  
  gSystem->Rename("event_stat.root", path + "/event_stat.root");      
  Printf(">>>>> Moved files to %s", path.Data());
}  



TChain * GetAnalysisChain(const char * incollection){
  // Builds a chain of esd files
  // incollection can be
  // - a single root file
  // - an xml collection of files on alien
  // - a ASCII containing a list of local root files
  TChain* analysisChain = 0;
  // chain
  analysisChain = new TChain("esdTree");
  if (TString(incollection).Contains(".root")){
    analysisChain->Add(incollection);
  }
  else if (TString(incollection).Contains("xml")){
    TGrid::Connect("alien://");
    TGridCollection * coll = gGrid->OpenCollection(incollection);
    while(coll->Next()){
      analysisChain->Add(TString("alien://")+coll->GetLFN());
    }
  } else {
    ifstream file_collect(incollection);
    TString line;
    while (line.ReadLine(file_collect) ) {
      analysisChain->Add(line.Data());
    }
  }
  analysisChain->GetListOfFiles()->Print();

  return analysisChain;
}


void InitAndLoadLibs(Int_t runMode=kMyRunModeLocal, Int_t workers=0,Bool_t debug=0) {
  // Loads libs and par files + custom task and classes

  // Custom stuff to be loaded
  //  listToLoad->Add(new TObjString("$ALICE_ROOT/ANALYSIS/AliCentralitySelectionTask.cxx+"));
  listToLoad->Add(new TObjString("/Users/mfloris/Work/AliSoft/AliRoot-trunk/PWGPP/background/AliHistoListWrapper.cxx+")); // FIXME
  listToLoad->Add(new TObjString("AliAnalysisMultPbTrackHistoManager.cxx+"));
  listToLoad->Add(new TObjString("AliAnalysisMultPbCentralitySelector.cxx+"));
  listToLoad->Add(new TObjString("AliAnalysisTaskMultPbTracks.cxx+"));
  // listToLoad->Add(new TObjString("$ALICE_ROOT/PWG0/multPbPb/AliAnalysisMultPbTrackHistoManager.cxx+"));
  // listToLoad->Add(new TObjString("$ALICE_ROOT/PWG0/multPbPb/AliAnalysisMultPbCentralitySelector.cxx+"));
  // listToLoad->Add(new TObjString("$ALICE_ROOT/PWG0/multPbPb/AliAnalysisTaskMultPbTracks.cxx+"));


  {
      cout << "Init in Local or Grid mode" << endl;
      gSystem->Load("libCore");
      gSystem->Load("libTree");
      gSystem->Load("libGeom");
      gSystem->Load("libVMC");
      gSystem->Load("libPhysics");
      gSystem->Load("libMinuit");
      gSystem->Load("libSTEERBase");
      gSystem->Load("libESD");
      gSystem->Load("libAOD");
      gSystem->Load("libANALYSIS");
      gSystem->Load("libOADB");
      gSystem->Load("libANALYSISalice");   
      // Use AliRoot includes to compile our task
      gROOT->ProcessLine(".include $ALICE_ROOT/include");

      // gSystem->Load("libVMC");
      // gSystem->Load("libTree");
      // gSystem->Load("libSTEERBase");
      // gSystem->Load("libESD");
      // gSystem->Load("libAOD");
      // gSystem->Load("libANALYSIS");
      // gSystem->Load("libANALYSISalice");
      // gSystem->Load("libPWG0base");
    
      gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG0/multPb"));
      gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWGPP/background"));
      //    gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWGPP/background/"));
    }
  // Load helper classes
  TIterator * iter = listToLoad->MakeIterator();
  TObjString * name = 0;
  while (name = (TObjString *)iter->Next()) {
    gSystem->ExpandPathName(name->String());
    cout << name->String().Data();
    {
      gROOT->LoadMacro(name->String()+(debug?"+g":""));   
    }
  }

}
