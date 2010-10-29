// TODO:
// 1. Check cuts for 2010 (Jochen?)


enum { kMyRunModeLocal = 0, kMyRunModeCAF};

TChain * GetAnalysisChain(const char * incollection);

void run(Char_t* data, Long64_t nev = -1, Long64_t offset = 0, Bool_t debug = kFALSE, Int_t runMode = 0, Bool_t isMC = 0, const char* option = "",Int_t workers = -1)
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
    mgr->SetMCtruthEventHandler(handler);
  }

  // physics selection
  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  physicsSelectionTask = AddTaskPhysicsSelection(isMC);



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

  AliESDtrackCuts * cuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2009(kTRUE);
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
  
  
  // load my task
  gROOT->ProcessLine(".L $ALICE_ROOT/PWG0/multPbPb/AddTaskMultPbPbTracks.C");
  AliAnalysisTaskMultPbTracks * task = AddTaskMultPbPbTracks("multPbPbtracks.root", cuts); // kTRUE enables DCA cut
  task->SetIsMC(isMC);
   if(isMC) task->GetHistoManager()->SetSuffix("MC");
  
  if (!mgr->InitAnalysis()) return;
	
  mgr->PrintStatus();
  
  if (runMode == kMyRunModeLocal ) {
    // If running in local mode, create chain of ESD files
    TChain * chain = GetAnalysisChain(data);
    chain->Print();
    mgr->StartAnalysis("local",chain,nev);
  } else if (runMode == kMyRunModeCAF) {
    mgr->StartAnalysis("proof",TString(data)+"#esdTree",nev);
  } else {
    cout << "ERROR: unknown run mode" << endl;        
  }

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
  }
  else
  {
    cout << "Init in Local mode" << endl;

    gSystem->Load("libVMC");
    gSystem->Load("libTree");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libPWG0base");
    
    gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG0"));
  }
  // Load helper classes
  // TODO: replace this by a list of TOBJStrings
  TString taskName("AliAnalysisTaskMultPbTracks.cxx+");
  TString histoManName("AliAnalysisMultPbTrackHistoManager.cxx+");
  TString centrName("AliAnalysisMultPbCentralitySelector.cxx+");
  TString listName("AliHistoListWrapper.cxx+");

  // Create, add task
  if (runMode == kMyRunModeCAF) {
    gProof->Load(listName+(debug?"+g":""));   
    gProof->Load(histoManName+(debug?"+g":""));
    gProof->Load(centrName+(debug?"+g":""));   
    gProof->Load(taskName+(debug?"+g":""));
  } else {
    gROOT->LoadMacro(listName+(debug?"+g":""));   
    gROOT->LoadMacro(histoManName+(debug?"+g":""));
    gROOT->LoadMacro(centrName+(debug?"+g":""));   
    gROOT->LoadMacro(taskName+(debug?"+g":""));    
  }


}
