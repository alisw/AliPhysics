#if ! defined (__CINT__) || defined (__MAKECINT__)
#include <Riostream.h>

#include <TStopwatch.h>

#include <TROOT.h>
#include <TClass.h>
#include <TSystem.h>
#include <TError.h>
#include <TChain.h>
#include <TGrid.h>
#include <TList.h>
#include <TMethodCall.h>
#include <TGridCollection.h>
#include <TGridResult.h>
#include <TGeoGlobalMagField.h>

#include <AliLog.h>
#include <AliCDBManager.h>
#include <AliGRPManager.h>
#include <AliGeomManager.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisDataContainer.h>
#include <AliMCEventHandler.h>
#include <AliESDInputHandler.h>

#include <AliMESbaseTask.h>
#include <AliMEStender.h>
#endif

TChain* MakeChainLST(const char* filename = NULL);
void run(const Char_t *files=NULL, Bool_t mc=kFALSE, Bool_t tpid=kTRUE,  Bool_t tchg=kTRUE,  Bool_t tpp=kTRUE, Long64_t nev=1234567890, Long64_t first = 0)
{
  TStopwatch timer;
  timer.Start();

  // VERY GENERAL SETTINGS
  //AliLog::SetGlobalLogLevel(AliLog::kError);
  if(gSystem->Load("libANALYSIS.so")<0) return;
  if(gSystem->Load("libANALYSISalice.so")<0) return;
  if(gSystem->Load("libTender.so")<0) return;
  if(gSystem->Load("libTenderSupplies.so")<0) return;
//   if(gSystem->Load("libMES.so")<0) return;
    if(gSystem->Load("libPWGLFspectra.so")<0) return;

  // DEFINE DATA CHAIN
  TChain *chain = NULL;
  if(!files) chain = MakeChainLST();
  else chain = MakeChainLST(files);

  if(!chain) return;
  chain->Lookup();
  chain->GetListOfFiles()->Print();
  Long64_t nfound=chain->GetEntries();
  printf("\tENTRIES FOUND [%lli] REQUESTED [%lli]\n", nfound, nev>nfound?nfound:nev);

  // BUILD ANALYSIS MANAGER
  AliAnalysisManager *mgr = new AliAnalysisManager("Multiplicity and Event Shape");
  AliESDInputHandler *esdH = new AliESDInputHandler();
  AliMCEventHandler *mcH(NULL);
  mgr->SetInputEventHandler(esdH);
  if(mc) mgr->SetMCtruthEventHandler(mcH = new AliMCEventHandler());
  //mgr->SetDebugLevel(10);
  mgr->SetSkipTerminate(kTRUE);

  // LOAD tasks
  // *******************  PID response  ******************
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  if(!mc) AddTaskPIDResponse();
  // else AddTaskPIDResponse(kTRUE,kTRUE,kTRUE,2);
  else AddTaskPIDResponse(kTRUE);

  // *******************  Tenders  ***********************
  AliTender *aliTender(NULL);
  gROOT->LoadMacro("$ALICE_PHYSICS/TENDER/TenderSupplies/AddTaskTender.C");
  if(!mc){    // for DATA
    aliTender = (AliTender*)AddTaskTender(!mc, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kTRUE, kFALSE, kFALSE);
       // (useV0, useTPC,  !!! useTOF=kFALSE for MC !!!, useTRD, usePID, useVTX, useT0, useEmc, usePtFix)
  } else {  // for MC
    aliTender = (AliTender*)AddTaskTender(!mc, kTRUE, kFALSE, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kFALSE);  // (useV0, useTPC,  !!! useTOF=kFALSE for MC !!!, useTRD, usePID, useVTX, useT0, useEmc, usePtFix)
  }
  aliTender->SetHandleOCDB(kTRUE);
  // aliTender->SetDefaultCDBStorage(Form("alien://folder=/alice/data/2010/OCDB?cacheFolder=%s/local", gSystem->ExpandPathName("$HOME")));
  // aliTender->SetDefaultCDBStorage(Form("local://%s/local/alice/data/2010/OCDB", gSystem->ExpandPathName("$HOME")));

// *******************  Physics Selection  *************
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(mc); // 0 = real data; 1 = MC

  // *******************  MES Tender  ******************
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/SPECTRA/MultEvShape/AddMEStender.C");
  AddMEStender(mc);

  // *******************  MES PID task  ******************
  if(tpid){
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/SPECTRA/MultEvShape/AddMESpidTask.C");
	AddMESpidTask(mc);
  }

  // *******************  MES CHG task  ******************
  if(tchg){
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/SPECTRA/MultEvShape/AddMESchgTask.C");
    AddMESchgTask(mc);
  }

  // *******************  MES ppCol task  ******************
  if(tpp){
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/SPECTRA/MultEvShape/AddMESppColTask.C");
    AddMESppColTask(mc);
  }


  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local", chain, nev, first);
  timer.Stop();
  timer.Print();
  // verbosity
  printf("\tCLEANING TASK LIST:\n");
  mgr->GetTasks()->Delete();

  if(mcH) delete mcH;
  delete esdH;
  delete chain;
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
      if(gSystem->Load("libNetx.so")<0) return NULL;
      if(gSystem->Load("libRAliEn.so")<0) return NULL;
      TGrid::Connect("alien://");
    }
    chain->Add(esdFile.Data());
  }
  fclose(fp);

  return chain;
}
