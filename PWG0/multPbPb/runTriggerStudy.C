enum { kMyRunModeLocal = 0, kMyRunModeCAF, kMyRunModeGRID};
//#define TENDER
TChain * GetAnalysisChain(const char * incollection);

void runTriggerStudy(Char_t* data, Long64_t nev = -1, Long64_t offset = 0, Bool_t debug = kFALSE, Int_t runMode = 0, Bool_t isMC = 0, Int_t ntrackletsKine = 100, Bool_t rejectBGV0Trigger = kFALSE, const char* option = "", Int_t workers = -1)
{
  // runMode:
  //
  // 0 local 
  // 1 proof
#ifdef TENDER
  TGrid::Connect("alien://");
#endif
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

  // Add tender
#ifdef TENDER
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/TenderSupplies/AddTaskTender.C");
  AliAnalysisTask* tender=0x0;
  if(!isMC)
    {
      tender = AddTaskTender(kTRUE);
      // tender->SetDebugLevel(10);
    }
  else
    {
      tender = AddTaskTender(kFALSE);
      // tender->SetDebugLevel(10);
    }
#endif
  
  // Add physics selection
  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  physicsSelectionTask = AddTaskPhysicsSelection(isMC,1,!isMC);//FIXME
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
  gSystem->Rename("event_stat.root", path + "/" + "event_stat.root");
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
    // cout << workers>0 ? Form("workers=%d",workers) : "workers=1x" << endl;
    // exit(1);
    TProof * p = TProof::Open("alice-caf.cern.ch", workers>0 ? Form("workers=%d",workers) : "workers=1x");
    p->Exec("TObject *o = gEnv->GetTable()->FindObject(\"Proof.UseMergers\"); gEnv->GetTable()->Remove(o);", kTRUE);
    //TProof::Open("skaf.saske.sk", workers>0 ? Form("workers=%d",workers) : "");
    
    // Enable the needed package (par fileS)
    // gProof->UploadPackage("$ALICE_ROOT/obj/STEERBase");
    // gProof->EnablePackage("$ALICE_ROOT/obj/STEERBase");
    // gProof->UploadPackage("$ALICE_ROOT/obj/ESD");
    // gProof->EnablePackage("$ALICE_ROOT/obj/ESD");
    // gProof->UploadPackage("$ALICE_ROOT/obj/AOD");
    // gProof->EnablePackage("$ALICE_ROOT/obj/AOD");
    // gProof->UploadPackage("$ALICE_ROOT/obj/ANALYSIS");
    // gProof->EnablePackage("$ALICE_ROOT/obj/ANALYSIS");
    // gProof->UploadPackage("$ALICE_ROOT/obj/OADB");
    // gProof->EnablePackage("$ALICE_ROOT/obj/OADB");
    // gProof->UploadPackage("$ALICE_ROOT/obj/ANALYSISalice");
    // gProof->EnablePackage("$ALICE_ROOT/obj/ANALYSISalice");
    // gProof->UploadPackage("$ALICE_ROOT/PWG0base");
    // gProof->EnablePackage("$ALICE_ROOT/PWG0base");
    gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG0/multPb"));
    gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG1/background"));
    gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/include "));
    gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/TOF "));

    // Use a precompiled tag
    TString alirootMode="";    // STEERBase,ESD,AOD,ANALYSIS,ANALYSISalice (default aliroot mode)
    //alirootMode="ALIROOT";     // $ALICE_ROOT/macros/loadlibs.C
    //  alirootMode="REC";     // $ALICE_ROOT/macros/loadlibsrec.C
    //  alirootMode="SIM";     // $ALICE_ROOT/macros/loadlibssim.C
    //  alirootMode="TRAIN";   // $ALICE_ROOT/macros/loadlibstrain.C (not working yet)
    //  alirootMode="CUSTOM";  // nothing is loaded, but aliroot variables are set (not working yet)
 
    TString extraLibs;
    extraLibs= ""; // not needed in default aliroot mode
    extraLibs+="CDB:RAWDatabase:STEER:TENDER:TRDbase:STAT:TRDrec:VZERObase:VZEROsim:VZEROrec:RAWDatarec:TPCbase:TPCrec:TPCcalib:TENDERSupplies:RAWDatabase:RAWDatarec:RAWDatasim:TOFbase:TOFrec";
    TList *list = new TList();
    // sets $ALIROOT_MODE on each worker to let proof to know to run in special mode
    list->Add(new TNamed("ALIROOT_MODE", alirootMode.Data()));
    // sets $ALIROOT_EXTRA_LIBS on each worker to let proof to know to load extra libs
    list->Add(new TNamed("ALIROOT_EXTRA_LIBS", extraLibs.Data()));
#ifdef TENDER
    list->Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));
#endif
    // connect to proof
    gProof->EnablePackage("VO_ALICE@AliRoot::v4-21-22-AN", list);
    //    gProof->Exec("TGrid::Connect(\"alien://\");");
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
    gSystem->Load("libOADB");
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
