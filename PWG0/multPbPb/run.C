// TODO:
// 1. Check cuts for 2010 (Jochen?)
// 2. Run with many centrality bins at once
#include <string.h>

enum { kMyRunModeLocal = 0, kMyRunModeCAF, kMyRunModeGRID};

TList * listToLoad = new TList();

TChain * GetAnalysisChain(const char * incollection);

void run(Char_t* data, Long64_t nev = -1, Long64_t offset = 0, Bool_t debug = kFALSE, Int_t runMode = 0, Bool_t isMC = 0, Int_t centrBin = 0, const char * centrEstimator = "VOM", const char* option = "",TString customSuffix = "", Int_t workers = -1)
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
  //  esdH->SetInactiveBranches("AliESDACORDE FMD ALIESDTZERO ALIESDZDC AliRawDataErrorLogs CaloClusters Cascades EMCALCells EMCALTrigger ESDfriend Kinks AliESDTZERO ALIESDACORDE MuonTracks TrdTracks");
  mgr->SetInputEventHandler(esdH);

  if(isMC) {
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetPreReadMode(AliMCEventHandler::kLmPreRead);
    mgr->SetMCtruthEventHandler(handler);
  }


  // If we are running on grid, we need the alien handler
  if (runMode == kMyRunModeGRID) {
    // Create and configure the alien handler plugin
    gROOT->LoadMacro("CreateAlienHandler.C");
    AliAnalysisGrid *alienHandler = CreateAlienHandler(data, listToLoad, "full", isMC);  
    if (!alienHandler) {
      cout << "Cannot create alien handler" << endl;    
      exit(1);
    }
    mgr->SetGridHandler(alienHandler);  
  }



  // physics selection
  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  physicsSelectionTask = AddTaskPhysicsSelection(isMC);

  // Centrality
  AliCentralitySelectionTask *taskCentr = new AliCentralitySelectionTask("CentralitySelection");
  taskCentr->SetPercentileFile("$ALICE_ROOT/ANALYSIS/macros/test_AliCentralityBy1D.root");
  taskCentr->SetPercentileFile2("$ALICE_ROOT/ANALYSIS/macros/test_AliCentralityByFunction.root");
  mgr->AddTask(taskCentr);
  mgr->ConnectInput (taskCentr,0, mgr->GetCommonInputContainer());


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

  AliESDtrackCuts * cuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
  TString pathsuffix = "";
  // cuts->SetPtRange(0.15,0.2);// FIXME pt cut
  // const char * pathsuffix = "_pt_015_020_nofakes";

  if (optionStr.Contains("ITSsa")) {
    delete cuts;
    cuts = AliESDtrackCuts::GetStandardITSPureSATrackCuts2009();
    cout << ">>>> USING ITS sa tracks" << endl;
    pathsuffix="ITSsa";
  }

  if (optionStr.Contains("TPC")) {
    delete cuts;
    cuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    cout << ">>>> USING TPC only tracks" << endl;
    pathsuffix="TPC";
  }

  Bool_t useMCKinematics = isMC;
  if (optionStr.Contains("NOMCKIN")) {
    cout << ">>>> Ignoring MC kinematics" << endl;
    useMCKinematics=kFALSE;
  }
  
  
  // load my task
  gROOT->ProcessLine(".L $ALICE_ROOT/PWG0/multPbPb/AddTaskMultPbPbTracks.C");
  AliAnalysisTaskMultPbTracks * task = AddTaskMultPbPbTracks("multPbPbtracks.root", cuts); // kTRUE enables DCA cut
  task->SetIsMC(useMCKinematics);
  if(useMCKinematics) task->GetHistoManager()->SetSuffix("MC");
  if(customSuffix!=""){
    cout << "Setting custom suffix: " << customSuffix << endl;    
    task->GetHistoManager()->SetSuffix(customSuffix);
  }
  task->SetCentralityBin(centrBin);
  task->SetCentralityEstimator(centrEstimator);
  
  if (!mgr->InitAnalysis()) return;
	
  mgr->PrintStatus();
  
  if (runMode == kMyRunModeLocal ) {
    // If running in local mode, create chain of ESD files
    cout << "RUNNING LOCAL, CHAIN" << endl;    
    TChain * chain = GetAnalysisChain(data);
    chain->Print();
    mgr->StartAnalysis("local",chain,nev);
  } else if (runMode == kMyRunModeCAF) {
    mgr->StartAnalysis("proof",TString(data)+"#esdTree",nev);
  } else if (runMode == kMyRunModeGRID) {
    mgr->StartAnalysis("grid");
  } else {
    cout << "ERROR: unknown run mode" << endl;        
  }

  pathsuffix = pathsuffix + "_" + centrEstimator + "_bin_"+long(centrBin);
  if (doSave) MoveOutput(data, pathsuffix.Data());

  
}


void MoveOutput(const char * data, const char * suffix = ""){

  TString path("output/");
  path = path + TString(data).Tokenize("/")->Last()->GetName() + suffix;
  
  TString fileName = "multPbPbtracks.root";
  gSystem->mkdir(path, kTRUE);
  gSystem->Rename(fileName, path + "/" + fileName);
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
    TAlienCollection * coll = TAlienCollection::Open (incollection);
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
  listToLoad->Add(new TObjString("$ALICE_ROOT/ANALYSIS/AliCentralitySelectionTask.cxx+"));
  listToLoad->Add(new TObjString("$ALICE_ROOT/PWG1/background/AliHistoListWrapper.cxx+"));
  listToLoad->Add(new TObjString("$ALICE_ROOT/PWG0/multPbPb/AliAnalysisMultPbTrackHistoManager.cxx+"));
  listToLoad->Add(new TObjString("$ALICE_ROOT/PWG0/multPbPb/AliAnalysisTaskMultPbTracks.cxx+"));


  if (runMode == kMyRunModeCAF)
  {
    cout << "Init in CAF mode" << endl;
    
    gEnv->SetValue("XSec.GSI.DelegProxy", "2");
    TProof::Open("alice-caf.cern.ch", workers>0 ? Form("workers=%d",workers) : "");
    
    // Enable the needed package
    gProof->UploadPackage("$ALICE_ROOT/STEERBase");
    gProof->EnablePackage("$ALICE_ROOT/STEERBase");
    gProof->UploadPackage("$ALICE_ROOT/ESD");
    gProof->EnablePackage("$ALICE_ROOT/ESD");
    gProof->UploadPackage("$ALICE_ROOT/AOD");
    gProof->EnablePackage("$ALICE_ROOT/AOD");
    gProof->UploadPackage("$ALICE_ROOT/ANALYSIS");
    gProof->EnablePackage("$ALICE_ROOT/ANALYSIS");
    gProof->UploadPackage("$ALICE_ROOT/ANALYSISalice");
    gProof->EnablePackage("$ALICE_ROOT/ANALYSISalice");
    gProof->UploadPackage("$ALICE_ROOT/PWG0base");
    gProof->EnablePackage("$ALICE_ROOT/PWG0base");
    gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG0/multPb"));
    gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG1/background"));
  }
  else
  {
    cout << "Init in Local or Grid mode" << endl;
    gSystem->Load("libCore.so");  
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
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
    gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG1/background"));
    //    gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG1/background/"));
  }
  // Load helper classes
  TIterator * iter = listToLoad->MakeIterator();
  TObjString * name = 0;
  while (name = (TObjString *)iter->Next()) {
    gSystem->ExpandPathName(name->String());
    cout << name->String().Data();
    if (runMode == kMyRunModeCAF) {
      gProof->Load(name->String()+(debug?"+g":""));   
    } else {
      gROOT->LoadMacro(name->String()+(debug?"+g":""));   
    }
  }

}
