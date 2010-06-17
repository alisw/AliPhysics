// Steer TRD QA train for Reconstruction (Clusterizer, Tracking and PID).
// 
// Usage:
//   run.C(tasks, files)
//   tasks : "ALL" or one/more of the following:
//     "EFF"  : TRD Tracking Efficiency 
//     "EFFC" : TRD Tracking Efficiency Combined (barrel + stand alone) - only in case of simulations
//     "MULT"  : TRD single track selection
//     "RES"  : TRD tracking Resolution
//     "CLRES": clusters Resolution
//     "CAL"  : TRD calibration
//     "ALGN" : TRD alignment
//     "PID"  : TRD PID - pion efficiency 
//     "PIDR" : TRD PID - reference data
//     "DET"  : Basic TRD Detector checks
//     "NOFR" : Data set does not have AliESDfriends.root 
//     "NOMC" : Data set does not have Monte Carlo Informations (real data), so all tasks which rely
//              on MC information are switched off
//
// In compiled mode : 
// Don't forget to load first the libraries
// gSystem->Load("libMemStat.so")
// gSystem->Load("libMemStatGui.so")
// gSystem->Load("libANALYSIS.so")
// gSystem->Load("libANALYSISalice.so")
// gSystem->Load("libPWG1.so");
// gSystem->Load("libNetx.so") ;
// gSystem->Load("libRAliEn.so");
//
// Authors:
//   Alex Bercuci (A.Bercuci@gsi.de) 
//   Markus Fasel (m.Fasel@gsi.de) 

#if ! defined (__CINT__) || defined (__MAKECINT__)
//#ifndef __CINT__
#include <Riostream.h>

#include "TStopwatch.h"
#include "TMemStat.h"
#include "TMemStatViewerGUI.h"

#include "TROOT.h"
#include "TClass.h"
#include "TSystem.h"
#include "TError.h"
#include "TChain.h"
#include "TGrid.h"
#include "TAlienCollection.h"
#include "TGridCollection.h"
#include "TGridResult.h"
#include "TGeoGlobalMagField.h"

#include "AliMagF.h"
#include "AliTracker.h"
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliGeomManager.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h"

#include "TRD/AliTRDtrackerV1.h"
#include "TRD/AliTRDcalibDB.h"

#include "PWG1/TRD/macros/AliTRDperformanceTrain.h"
#include "PWG1/TRD/macros/AddTRDcheckESD.C"
#include "PWG1/TRD/macros/AddTRDinfoGen.C"
#include "PWG1/TRD/macros/AddTRDcheckDET.C"
#include "PWG1/TRD/macros/AddTRDefficiency.C"
#include "PWG1/TRD/macros/AddTRDresolution.C"
#include "PWG1/TRD/macros/AddTRDcheckPID.C"

#endif

#include "macros/AliTRDperformanceTrain.h"


Bool_t MEM = kFALSE;

TChain* MakeChainLST(const char* filename = 0x0);
TChain* MakeChainXML(const char* filename = 0x0);
void run(Char_t *optList="ALL", const Char_t *files=0x0, Long64_t nev=1234567890, Long64_t first = 0, Int_t runNo=0, const Char_t *ocdb_uri="local://$ALICE_ROOT/OCDB", const Char_t *grp_uri=Form("local://%s", gSystem->pwd()))
{
  TMemStat *mem = 0x0;
  if(MEM){ 
    if(gSystem->Load("libMemStat.so")<0) return;
    if(gSystem->Load("libMemStatGui.so")<0) return;
    mem = new TMemStat("new, gnubuildin");
    mem->AddStamp("Start");
  }
  TStopwatch timer;
  timer.Start();

  // VERY GENERAL SETTINGS
  AliLog::SetGlobalLogLevel(AliLog::kError);
  if(gSystem->Load("libANALYSIS.so")<0) return;
  if(gSystem->Load("libANALYSISalice.so")<0) return;
  if(gSystem->Load("libTENDER.so")<0) return;
  if(gSystem->Load("libPWG1.so")<0) return;

  Bool_t fHasMCdata =  HasReadMCData(optList);
  Bool_t fHasFriends = HasReadFriendData(optList);
  
  // INITIALIZATION OF RUNNING ENVIRONMENT
  // initialize OCDB manager
  AliCDBManager *cdbManager = AliCDBManager::Instance();
  cdbManager->SetDefaultStorage(ocdb_uri);
  if(!cdbManager->IsDefaultStorageSet()){
    Error("run.C", "Error setting OCDB.");
    return;
  }
  cdbManager->SetRun(runNo);
  cdbManager->SetSpecificStorage("GRP/GRP/Data", grp_uri);
  cdbManager->SetCacheFlag(kFALSE);
  // initialize magnetic field from the GRP manager.
  AliGRPManager grpMan;
  if(!grpMan.ReadGRPEntry()) return;
  if(!grpMan.SetMagField()) return;
  //AliRunInfo *runInfo = grpMan.GetRunInfo();
  AliGeomManager::LoadGeometry();

  // DEFINE DATA CHAIN
  TChain *chain = 0x0;
  if(!files) chain = MakeChainLST();
  else{
    TString fn(files);
    if(fn.EndsWith("xml")) chain = MakeChainXML(files);
    else chain = MakeChainLST(files);
  }
  if(!chain) return;
  chain->Lookup();
  chain->GetListOfFiles()->Print();
  Int_t nfound=(Int_t)chain->GetEntries();
  printf("\tENTRIES FOUND [%d] REQUESTED [%d]\n", nfound, nev>nfound?nfound:nev);


  // BUILD ANALYSIS MANAGER
  AliAnalysisManager *mgr = new AliAnalysisManager("TRD Reconstruction Performance & Calibration");
  AliESDInputHandlerRP *esdH(NULL);
  mgr->SetInputEventHandler(esdH = new AliESDInputHandlerRP);
  esdH->SetReadFriends(kTRUE);
  esdH->SetActiveBranches("ESDfriend");
  AliMCEventHandler *mcH(NULL);
  if(fHasMCdata) mgr->SetMCtruthEventHandler(mcH = new AliMCEventHandler());
  //mgr->SetDebugLevel(10);

  gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTrainPerformanceTRD.C");
  if(!AddTrainPerformanceTRD(optList)) {
    Error("run.C", "Error loading TRD train.");
    return;
  }

  if (!mgr->InitAnalysis()) return;
  // verbosity
  printf("\tRUNNING TRAIN FOR TASKS:\n");
  TObjArray *taskList=mgr->GetTasks();
  for(Int_t itask=0; itask<taskList->GetEntries(); itask++){ 
    AliAnalysisTask *task=(AliAnalysisTask*)taskList->At(itask);
    printf(" %s [%s]\n", task->GetName(), task->GetTitle());
  }
  //mgr->PrintStatus();
  mgr->StartAnalysis("local", chain, nev, first);

  timer.Stop();
  timer.Print();  

  delete cdbManager;

  // verbosity
  printf("\tCLEANING TASK LIST:\n");
  mgr->GetTasks()->Delete();

  if(mcH) delete mcH;
  delete esdH;
  delete chain;
  if(MEM) delete mem;
  if(MEM) TMemStatViewerGUI::ShowGUI();
}

//____________________________________________
TChain* MakeChainLST(const char* filename)
{
  // Create the chain
  TChain* chain = new TChain("esdTree");

  if(!filename){
    chain->Add(Form("%s/AliESDs.root", gSystem->pwd()));
    return chain;
  }


  // read ESD files from the input list.
  ifstream in;
  in.open(filename);
  TString esdfile;
  while(in.good()) {
    in >> esdfile;
    if (!esdfile.Contains("root")) continue; // protection
    chain->Add(esdfile.Data());
  }

  in.close();

  return chain;
}

//____________________________________________
TChain* MakeChainXML(const char* xmlfile)
{
  if (!TFile::Open(xmlfile)) {
    Error("MakeChainXML", Form("No file %s was found", xmlfile));
    return 0x0;
  }

  if(gSystem->Load("libNetx.so")<0) return 0x0;
  if(gSystem->Load("libRAliEn.so")<0) return 0x0;
  TGrid::Connect("alien://") ;

  TGridCollection *collection = (TGridCollection*) TAlienCollection::Open(xmlfile);
  if (!collection) {
    Error("MakeChainXML", Form("No collection found in %s", xmlfile)) ; 
    return 0x0; 
  }
  //collection->CheckIfOnline();

  TGridResult* result = collection->GetGridResult("",0 ,0);
  if(!result->GetEntries()){
    Error("MakeChainXML", Form("No entries found in %s", xmlfile)) ; 
    return 0x0; 
  }
  // Makes the ESD chain 
  TChain* chain = new TChain("esdTree");
  for (Int_t idx = 0; idx < result->GetEntries(); idx++) {
    chain->Add(result->GetKey(idx, "turl")); 
  }
  return chain;
}
