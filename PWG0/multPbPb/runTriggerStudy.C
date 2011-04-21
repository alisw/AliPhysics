enum { kMyRunModeLocal = 0, kMyRunModeCAF, kMyRunModeGRID};

TChain * GetAnalysisChain(const char * incollection);

void runTriggerStudy(Char_t* data, Long64_t nev = -1, Long64_t offset = 0, Bool_t debug = kFALSE, Int_t runMode = 0, Bool_t isMC = 0, Int_t ntrackletsKine = 100, Bool_t rejectBGV0Trigger = kFALSE, const char* option = "", Int_t workers = -1)
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
    gROOT->LoadMacro("CreateAlienHandlerTrigger.C");
    AliAnalysisGrid *alienHandler = CreateAlienHandlerTrigger(data,"pass1",isMC);  
    if (!alienHandler) {
      cout << "Cannot create alien handler" << endl;    
      exit(1);
    }
    mgr->SetGridHandler(alienHandler);  
  }

  // Add physics selection
  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  physicsSelectionTask = AddTaskPhysicsSelection(isMC,0);//FIXME
  physicsSelectionTask->GetPhysicsSelection()->SetSkipZDCTime(1);// Skip ZDC - applyied later



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

  
  
  // load my task
  AliAnalysisTaskTriggerStudy *task = new AliAnalysisTaskTriggerStudy("TaskOfflineTrigger");
  mgr->AddTask(task);
  // Set I/O
  AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("cTrigStudy",
							    AliHistoListWrapper::Class(),
							    AliAnalysisManager::kOutputContainer,
							    "Trig_Temp.root");
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);



  task->SetIsMC(isMC);
  task->SetNTrackletsCutKine(ntrackletsKine);
  task->SetRejectBGWithV0(rejectBGV0Trigger);

  if (!mgr->InitAnalysis()) return;
	
  mgr->PrintStatus();
  
  if (runMode == kMyRunModeLocal ) {
    // If running in local mode, create chain of ESD files
    cout << "RUNNING LOCAL, CHAIN" << endl;    
    TChain * chain = GetAnalysisChain(data);
    //    chain->Print();
    mgr->StartAnalysis("local",chain,nev);
  } else if (runMode == kMyRunModeCAF) {
    mgr->StartAnalysis("proof",TString(data)+"#esdTree",nev);
  } else if (runMode == kMyRunModeGRID) {
    mgr->StartAnalysis("grid");
  }else {
    cout << "ERROR: unknown run mode" << endl;        
  }

  if (doSave) MoveOutput(data, Form("_TrkCut_%d_V0BGCUT_%d",ntrackletsKine,rejectBGV0Trigger));

  
}


void MoveOutput(const char * data, const char * suffix = ""){

  TString path("outTrigger/");
  path = path + TString(data).Tokenize("/")->Last()->GetName() + suffix;
  
  TString fileName = "trigger_study.root";
  gSystem->mkdir(path, kTRUE);
  gSystem->Rename(fileName, path + "/" + fileName);
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
    TProof * p = TProof::Open("alice-caf.cern.ch", workers>0 ? Form("workers=%d",workers) : "");
    p->Exec("TObject *o = gEnv->GetTable()->FindObject(\"Proof.UseMergers\"); gEnv->GetTable()->Remove(o);", kTRUE);
    //TProof::Open("skaf.saske.sk", workers>0 ? Form("workers=%d",workers) : "");
    
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

    gSystem->Load("libCore");  
    gSystem->Load("libGeom");
    gSystem->Load("libPhysics");
    gSystem->Load("libVMC");
    gSystem->Load("libTree");
    gSystem->Load("libProof");
    gSystem->Load("libTree");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libPWG0base");
    
    gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG0/multPb"));
    gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG1/background"));
    //    gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG1/background/"));
  }
  // Load helper classes
  // TODO: replace this by a list of TOBJStrings
  TString taskName("$ALICE_ROOT/PWG0/multPbPb/AliAnalysisTaskTriggerStudy.cxx+");
  TString listName("$ALICE_ROOT/PWG1/background/AliHistoListWrapper.cxx+");

  gSystem->ExpandPathName(taskName);
  gSystem->ExpandPathName(listName);



  // Create, add task
  if (runMode == kMyRunModeCAF) {
    gProof->Load(listName+(debug?"+g":""));   
    gProof->Load(taskName+(debug?"+g":""));
  } else {
    gROOT->LoadMacro(listName+(debug?"+g":""));   
    gROOT->LoadMacro(taskName+(debug?"+g":""));    
  }


}
