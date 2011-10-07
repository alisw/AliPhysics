// Steer TRD QA train for Reconstruction (Clusterizer, Tracking and PID).
// 
// Usage:
//   run.C(optList, files, nev, first, runNo, ocdb_uri, grp_uri)
//
//   optList : "ALL" [default] or one/more of the following:
//     "EFF"  : TRD Tracking Efficiency 
//     "EFFC" : TRD Tracking Efficiency Combined (barrel + stand alone) - only in case of simulations
//     "MULT"  : TRD single track selection
//     "RES"  : TRD tracking Resolution
//     "CLRES": clusters Resolution
//     "CAL"  : TRD calibration
//     "ALGN" : TRD alignment
//     "PID"  : TRD PID - pion efficiency 
//     "PIDR" : TRD PID - reference data
//     "V0"   : monitor V0 performance for use in TRD PID calibration
//     "DET"  : Basic TRD Detector checks
//      ****** SPECIAL OPTIONS **********
//     "NOFR" : Data set does not have AliESDfriends.root 
//     "NOMC" : Data set does not have Monte Carlo Informations (default have MC), 
//
//     files : the list of ESD files to be processed [default AliESds.root from cwd]
//     nev   : number of events to be processed [default all]
//     first : first event to process [default 0]
//     runNo : run number [default 0]
//     ocdb_uri : OCDB location [default local, $ALICE_ROOT]. In case of AliEn the syntax should be of the form
//                 alien://folder=/alice/data/2010/OCDB?user=username?cacheF=yourDirectory?cacheS=yourSize
//     grp_uri  : GRP/GRP/Data location [default cwd]. In case of AliEn the syntax should be of the form
//                 alien://folder=/alice/data/2010/OCDB?user=username?cacheF=yourDirectory?cacheS=yourSize
// In compiled mode : 
// Don't forget to load first the libraries
// gSystem->Load("libMemStat.so")
// gSystem->Load("libMemStatGui.so")
// gSystem->Load("libANALYSIS.so")
// gSystem->Load("libANALYSISalice.so")
// gSystem->Load("libTENDER.so");
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

#include "PWG1/TRD/macros/AddTRDcheckESD.C"
#include "PWG1/TRD/macros/AddTRDinfoGen.C"
#include "PWG1/TRD/macros/AddTRDcheckDET.C"
#include "PWG1/TRD/macros/AddTRDefficiency.C"
#include "PWG1/TRD/macros/AddTRDresolution.C"
#include "PWG1/TRD/macros/AddTRDcheckPID.C"
#endif

Bool_t MEM = kFALSE;

TChain* MakeChainLST(const char* filename = NULL);
TChain* MakeChainXML(const char* filename = NULL);
Bool_t UseMC(Char_t *opt);
Bool_t UseFriends(Char_t *opt);
void run(Char_t *optList="ALL", Int_t run, const Char_t *files=NULL, Long64_t nev=1234567890, Long64_t first = 0)
{
  TMemStat *mem = NULL;
  if(MEM){ 
    if(gSystem->Load("libMemStat.so")<0) return;
    if(gSystem->Load("libMemStatGui.so")<0) return;
    mem = new TMemStat("new, gnubuildin");
    mem->AddStamp("Start");
  }
  TStopwatch timer;
  timer.Start();

  // VERY GENERAL SETTINGS
  //AliLog::SetGlobalLogLevel(AliLog::kError);
  if(gSystem->Load("libANALYSIS.so")<0) return;
  if(gSystem->Load("libANALYSISalice.so")<0) return;
  if(gSystem->Load("libTENDER.so")<0) return;
  if(gSystem->Load("libPWG1.so")<0) return;

  Bool_t fHasMCdata = UseMC(optList);
  Bool_t fHasFriends = UseFriends(optList);
  
  // DEFINE DATA CHAIN
  TChain *chain = NULL;
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
  if(fHasFriends){
    esdH->SetReadFriends(kTRUE);
    esdH->SetActiveBranches("ESDfriend");
  }
  AliMCEventHandler *mcH(NULL);
  if(fHasMCdata) mgr->SetMCtruthEventHandler(mcH = new AliMCEventHandler());
  //mgr->SetDebugLevel(10);
  mgr->SetSkipTerminate(kTRUE);

  // add CDB task
  gROOT->LoadMacro("$ALICE_ROOT/PWG1/PilotTrain/AddTaskCDBconnect.C");
  AliTaskCDBconnect *taskCDB = AddTaskCDBconnect();
  if (!taskCDB) return;
  taskCDB->SetRunNumber(run);

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
    return NULL;
  }

  if(gSystem->Load("libNetx.so")<0) return NULL;
  if(gSystem->Load("libRAliEn.so")<0) return NULL;
  TGrid::Connect("alien://") ;

  TGridCollection *collection = (TGridCollection*) TAlienCollection::Open(xmlfile);
  if (!collection) {
    Error("MakeChainXML", Form("No collection found in %s", xmlfile)) ; 
    return NULL;
  }
  //collection->CheckIfOnline();

  TGridResult* result = collection->GetGridResult("",0 ,0);
  if(!result->GetEntries()){
    Error("MakeChainXML", Form("No entries found in %s", xmlfile)) ; 
    return NULL;
  }
  // Makes the ESD chain 
  TChain* chain = new TChain("esdTree");
  for (Int_t idx = 0; idx < result->GetEntries(); idx++) {
    chain->Add(result->GetKey(idx, "turl")); 
  }
  return chain;
}

//______________________________________________________
Bool_t UseMC(Char_t *opt){
  return !(Bool_t)strstr(opt, "NOMC");
}

//____________________________________________
Bool_t UseFriends(Char_t *opt){
  return !(Bool_t)strstr(opt, "NOFR");
}
