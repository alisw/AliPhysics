// Steer TRD QA train for Reconstruction (Clusterizer, Tracking and PID).
// 
// Usage:Char_t *optList="ALL", Int_t run, const Char_t *files=NULL, Long64_t nev=1234567890, Long64_t first = 0
//   run.C(optList, run, files, nev, first)
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
//     run   : run number [default 0]; if negative value is specified it should be the -YEAR of the run and the OCDB connection is left to AddTRDinfoGen.C.
//             A local cache should be performed first via cacheOCDB.C macro 
//     files : the list of ESD files to be processed [default AliESds.root from cwd]
//     nev   : number of events to be processed [default all]
//     first : first event to process [default 0]
//
// In compiled mode :
// Don't forget to load first the libraries
// gSystem->Load("libMemStat")
// gSystem->Load("libMemStatGui")
// gSystem->Load("libANALYSIS")
// gSystem->Load("libANALYSISalice")
// gSystem->Load("libTender");
// gSystem->Load("libCORRFW");
// gSystem->Load("libPWGPP");
// gSystem->Load("libPWGmuon");
// gSystem->Load("libNetx") ;
// gSystem->Load("libRAliEn");
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

#include "TGrid.h"
#include "TROOT.h"
#include "TClass.h"
#include "TSystem.h"
#include "TError.h"
#include "TChain.h"
#include "TGrid.h"
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

#include "PWGPP/TRD/macros/AddTRDcheckESD.C"
#include "PWGPP/TRD/macros/AddTRDinfoGen.C"
#include "PWGPP/TRD/macros/AddTRDcheckDET.C"
#include "PWGPP/TRD/macros/AddTRDefficiency.C"
#include "PWGPP/TRD/macros/AddTRDresolution.C"
#include "PWGPP/TRD/macros/AddTRDcheckPID.C"
#endif

Bool_t MEM = kFALSE;

TChain* MakeChainLST(const char* filename = NULL);
TChain* MakeChainXML(const char* filename = NULL);
Bool_t UseMC(Char_t *opt);
Bool_t UseFriends(Char_t *opt);
void run(Char_t *optList="ALL", Int_t run=0, const Char_t *files=NULL, Long64_t nev=1234567890, Long64_t first = 0, const Char_t *macroDir=0)
{
  TMemStat *mem = NULL;
  if(MEM){ 
    if(gSystem->Load("libMemStat")<0) return;
    if(gSystem->Load("libMemStatGui")<0) return;
    mem = new TMemStat("new, gnubuildin");
    mem->AddStamp("Start");
  }
  TStopwatch timer;
  timer.Start();

  // VERY GENERAL SETTINGS
  //AliLog::SetGlobalLogLevel(AliLog::kError);
  if(gSystem->Load("libANALYSIS")<0) return;
  //if(gSystem->Load("libANALYSISalice")<0) return;
  if(gSystem->Load("libTender")<0) return;
  if(gSystem->Load("libCORRFW")<0) return;
  if(gSystem->Load("libPWGPP")<0) return;
  //if(gSystem->Load("libPWGmuon")<0) return;

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
  if(run>=0){
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
    AliTaskCDBconnect *taskCDB = AddTaskCDBconnect("cvmfs://", run);
    if (!taskCDB) return;
    //taskCDB->SetRunNumber(run);
  } else {
    Warning("run.C", "OCDB connection via AliTRDinfoGen.");
    AliTRDpwgppHelper::SetRunYear(-run);
  }
  if(!AliTRDpwgppHelper::AddTrainPerformanceTRD(optList, macroDir?macroDir:"$ALICE_PHYSICS/PWGPP/TRD/macros")) {
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
  FILE *fp(NULL);
  if(!(fp = fopen(filename, "rt"))){
    Error("run::MakeChainLST()", Form("Input list \"%s\" not readable", filename));
    return NULL;
  }
  TString esdFile;
  while(esdFile.Gets(fp)){
    if (!esdFile.Contains("root")) continue; // protection
    if(esdFile.BeginsWith("alien://") && !gGrid){
      if(gSystem->Load("libNetx")<0) return NULL;
      if(gSystem->Load("libRAliEn")<0) return NULL;
      TGrid::Connect("alien://");
    }
    chain->Add(esdFile.Data());
  }
  fclose(fp);

  return chain;
}

//____________________________________________
TChain* MakeChainXML(const char* xmlfile)
{
  if (!TFile::Open(xmlfile)) {
    Error("MakeChainXML", Form("No file %s was found", xmlfile));
    return NULL;
  }

  if(gSystem->Load("libNetx")<0) return NULL;
  if(gSystem->Load("libRAliEn")<0) return NULL;
  TGrid::Connect("alien://") ;

  TGridCollection *collection = gGrid->OpenCollection(xmlfile);
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
