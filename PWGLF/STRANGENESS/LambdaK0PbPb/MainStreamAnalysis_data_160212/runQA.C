// #include <iostream>
// #include "AliAnalysisManager.h"
// #include "AliESDInputHandler.h"
// #include "AliMCEventHandler.h"
// #include "AliAnalysisGrid.h"
// #include "AliCentralitySelectionTask.h"
// #include "AliAnalysisCentralitySelector.h"
// #include "AliAnalysisTaskPerformanceStrange.h"
// #include "TString.h"
// #include "TChain.h"
// #include "TAlienCollection.h"
// #include <fstream>
// #include "TObjString.h"
// #include "TIterator.h"
// #include "TGrid.h"
// #include "TROOT.h"

// #include "CreateAlienHandler.C"
// #include 

using namespace std;

enum { kMyRunModeLocal = 0, kMyRunModeCAF, kMyRunModeGRID};

TList * listToLoad = new TList(); // Additional classes to be loaded, see InitAndLoadLibs

TChain * GetAnalysisChain(const char * incollection);
void InitAndLoadLibs(Int_t runMode=kMyRunModeLocal, Int_t workers=0,Bool_t debug=0) ;

void run(const char * data, const char * passOrPath, Long64_t nev = -1, Long64_t offset = 0, Bool_t debug = kFALSE, Int_t runMode = 0, Bool_t isMC = 0, Bool_t usePID = kTRUE, const char* option = "",TString customSuffix = "", Int_t workers = -1, const char * gridMode="full", Int_t binMin=0, Int_t binMax = 6)
{
  // runMode:
  //
  // 0 local 
  // 1 proof
  // 2 grid

  if (nev < 0)
    nev = 1234567890;
  InitAndLoadLibs(runMode,workers,debug);

  // Create the analysis manager
  AliAnalysisManager * mgr = new AliAnalysisManager;

  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler;
  // Do I need any of this? 
  // esdH->SetInactiveBranches("AliESDACORDE FMD ALIESDTZERO ALIESDZDC AliRawDataErrorLogs CaloClusters Cascades EMCALCells EMCALTrigger ESDfriend Kinks AliESDTZERO ALIESDACORDE MuonTracks TrdTracks");
  mgr->SetInputEventHandler(esdH);

  if(isMC) {
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetPreReadMode(AliMCEventHandler::kLmPreRead);
    mgr->SetMCtruthEventHandler(handler);
  }


  // If we are running on grid, we need the alien handler
  if (runMode == kMyRunModeGRID) {
    // Create and configure the alien handler plugin
    TGrid::Connect("alien://");// Why do I need this? Without a get a bus error...
        gROOT->LoadMacro("CreateAlienHandler.C");
    AliAnalysisGrid *alienHandler = CreateAlienHandler(data, listToLoad, gridMode, isMC);  
    if (!alienHandler) {
      cout << "Cannot create alien handler" << endl;    
      exit(1);
    }
    mgr->SetGridHandler(alienHandler);  
  }
  

  // Physics selection
  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask * physicsSelectionTask = AddTaskPhysicsSelection(isMC);

  // Centrality
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
  taskCentrality->SetPass(2);
  if(isMC) taskCentrality->SetMCInput();

  AliAnalysisCentralitySelector * centrSelector = new AliAnalysisCentralitySelector();
  centrSelector->SetIsMC(isMC);
  centrSelector->SetCentralityEstimator("V0M"); // Todo: add parameter to macro?

  // Parse option strings
  TString optionStr(option);
  
  // remove SAVE option if set
  Bool_t doSave = kFALSE;

  if (optionStr.Contains("SAVE"))
    {
      optionStr = optionStr(0,optionStr.Index("SAVE")) + optionStr(optionStr.Index("SAVE")+4, optionStr.Length());
      doSave = kTRUE;
    }
  TString pathsuffix = "";
  // Not used, but may be useful
  Bool_t useMCKinematics = isMC;
  if (optionStr.Contains("NOMCKIN")) {
    cout << ">>>> Ignoring MC kinematics" << endl;
    useMCKinematics=kFALSE;
    pathsuffix+="_NOMCKIN";
  }
  
  gROOT->ProcessLine(".L AddTaskLambdaK0PbPb.C");
  Int_t nbin = 0; // will contain the number of centrality bins
  AliAnalysisTaskPerformanceStrange ** task = AddTaskLambdaK0PbPb("lambdak0.root", centrSelector, nbin, binMin, binMax); // FIXME also pass cuts, centrality bin type selection(5,10% percentiles, ranges...)
  // configure task
  //  else if (iAODanalysis) task->SetAnalysisType("AOD");
  // FIXME: add options to macro
  // FIXME: put physics selection inside the task
  cout << nbin << endl;
  
  for(Int_t ibin = 0; ibin < nbin; ibin++){
    cout << "ibin " << ibin << "  "<< endl;//task[ibin] << endl;
    
    task[ibin]->SetAnalysisType("ESD");
    cout << "1" << endl;
    task[ibin]->SetAnalysisMC(isMC); // 0 or 1
    cout << "2" << endl;
    task[ibin]->SetCollidingSystems(1); // 0 =pp, 1=AA
    cout << "3" << endl;
    task[ibin]->SetAnalysisCut("no");
    cout << "4" << endl;
    if(usePID) 
      task[ibin]->SetUsePID("withPID"); // withPID or withoutPID
    else
      task[ibin]->SetUsePID("withoutPID"); // withPID or withoutPID
    cout << "5" << endl;
  }

  // Init and run the analy
  if (!mgr->InitAnalysis()) return;

  mgr->PrintStatus();
  
  if (runMode == kMyRunModeLocal ) {
    // If running in local mode, create chain of ESD files
    cout << "RUNNING LOCAL, CHAIN" << endl;    
    TChain * chain = GetAnalysisChain(data);
    //    chain->Print();
    mgr->StartAnalysis("local",chain,nev);
  } else if (runMode == kMyRunModeCAF) {
    mgr->StartAnalysis("proof",TString(passOrPath)+data+"#esdTree",nev);
  } else if (runMode == kMyRunModeGRID) {
    mgr->StartAnalysis("grid");
  } else {
    cout << "ERROR: unknown run mode" << endl;        
  }

  pathsuffix += customSuffix;

  if (doSave) MoveOutput(data, pathsuffix.Data());

  
}

void MoveOutput(const char * data, const char * suffix = ""){

  TString path("output/");
  path = path + TString(data).Tokenize("/")->Last()->GetName() + suffix;
  
  TString fileName = "lambdak0.root";
  gSystem->mkdir(path, kTRUE);
  gSystem->Rename(fileName, path + "/" + fileName);
  for(Int_t ibin = 0; ibin < 20; ibin++){
    TString fileBin = fileName;
    fileBin.ReplaceAll(".root",Form("_%2.2d.root",ibin));
    gSystem->Rename(fileBin, path + "/" + fileBin);    
  }
  
  gSystem->Rename("event_stat.root", path + "/event_stat.root");      
  gSystem->Rename("EventStat_temp.root", path + "/EventStat_temp.root");      
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
    TGridCollection * coll = TAlienCollection::Open (incollection);
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


void InitAndLoadLibs(Int_t runMode, Int_t workers,Bool_t debug) {
  // Loads libs and par files + custom task and classes (the order is important)
  // listToLoad->Add(new TObjString("$ALICE_ROOT/STEER/AliCentrality.cxx")); // FIXME: why do I have to load it?!?
  listToLoad->Add(new TObjString("AliAnalysisCentralitySelector.cxx"));
  listToLoad->Add(new TObjString("AliAnalysisTaskPerformanceStrange.cxx"));

  if (runMode == kMyRunModeCAF)
    {
      cout << "Init in CAF mode" << endl;
    
      gEnv->SetValue("XSec.GSI.DelegProxy", "2");
      Char_t* alienuser = gSystem->Getenv("alien_API_USER");
      TProof * p = TProof::Open(alienuser!=0 ? Form("%s@alice-caf.cern.ch",alienuser) : "alice-caf.cern.ch", workers>0 ? Form("workers=%d",workers) : "");
      //TProof * p = TProof::Open("skaf.saske.sk", workers>0 ? Form("workers=%d",workers) : "");    
      p->Exec("TObject *o = gEnv->GetTable()->FindObject(\"Proof.UseMergers\"); gEnv->GetTable()->Remove(o);", kTRUE); // avoid submerging
      //gProof->EnablePackage("VO_ALICE@AliRoot::v4-21-18-AN");
      gProof->EnablePackage("VO_ALICE@AliRoot::v4-21-20-AN");


      // Enable the needed package
      // FIXME: what if I don't want to use par files?
      gSystem->AddIncludePath("-I${ALICE_ROOT}/include/");
      gSystem->AddIncludePath("-I${ALICE_ROOT}/STEER/");
      // gProof->UploadPackage("$ALICE_ROOT/obj/STEERBase");
      // gProof->EnablePackage("$ALICE_ROOT/obj/STEERBase");
      // gProof->UploadPackage("$ALICE_ROOT/obj/ESD");
      // gProof->EnablePackage("$ALICE_ROOT/obj/ESD");
      // gProof->UploadPackage("$ALICE_ROOT/obj/AOD");
      // gProof->EnablePackage("$ALICE_ROOT/obj/AOD");
      // gProof->UploadPackage("$ALICE_ROOT/obj/ANALYSIS");
      // gProof->EnablePackage("$ALICE_ROOT/obj/ANALYSIS");
      // gProof->UploadPackage("$ALICE_ROOT/obj/ANALYSISalice");
      // gProof->EnablePackage("$ALICE_ROOT/obj/ANALYSISalice");
      // gProof->UploadPackage("$ALICE_ROOT/obj/CORRFW");
      // gProof->EnablePackage("$ALICE_ROOT/obj/CORRFW");
      // gProof->UploadPackage("~/Desktop/OADB");//FIXME
      // gProof->EnablePackage("~/Desktop/OADB");//FIXME
      
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

      // gSystem->Load("libVMC");
      gROOT->ProcessLine(".include $ALICE_ROOT/include");
      gROOT->ProcessLine(".include $ALICE_ROOT/STEER");
    }
  // Load helper classes
  TIterator * iter = listToLoad->MakeIterator();
  TObjString * name = 0;
  while ((name = (TObjString *)iter->Next())) {
    gSystem->ExpandPathName(name->String());
    cout << name->String().Data() << endl;
    if (runMode == kMyRunModeCAF) {
      gProof->Load(name->String()+(debug?"++g":"+"));   
    } else {
      gROOT->LoadMacro(name->String()+(debug?"++g":"+"));   
    }
  }

}
