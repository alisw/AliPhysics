// Steer TRD QA train for Reconstruction (Clusterizer, Tracking and PID).
// 
// Usage:
//   run.C(tasks, files, entries)
//   tasks : "ALL" or one/more of the following:
//     "EFF"  : TRD Tracking Efficiency 
//     "EFFC" : TRD Tracking Efficiency Combined (barrel + stand alone) - only in case of simulations
//     "RES"  : TRD tracking Resolution
//     "CAL"  : TRD calibration
//     "PID"  : TRD PID - pion efficiency 
//     "PIDR" : TRD PID - reference data
//     "NOMC" : Data set does not have Monte Carlo Informations (real data), so all tasks which rely
//              on MC information are switched off
//
// In compiled mode : 
// Don't forget to load first the libraries
// gSystem->Load("libMemStat.so")
// gSystem->Load("libMemStatGui.so")
// gSystem->Load("libANALYSIS.so")
// gSystem->Load("libTRDqaRec.so")
//
// Authors:
//   Alex Bercuci (A.Bercuci@gsi.de) 
//   Markus Fasel (m.Fasel@gsi.de) 

#ifndef __CINT__
#include <Riostream.h>

#include "TStopwatch.h"
#include "TMemStat.h"
#include "TMemStatViewerGUI.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TError.h"
#include "TChain.h"

#include "AliMagFMaps.h"
#include "AliTracker.h"
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h"

#include "TRD/AliTRDtrackerV1.h"
#include "TRD/AliTRDcalibDB.h"
#include "TRD/qaRec/AliTRDtrackInfoGen.h"
#include "TRD/qaRec/AliTRDtrackingEfficiency.h"
#include "TRD/qaRec/AliTRDtrackingEfficiencyCombined.h"
#include "TRD/qaRec/AliTRDtrackingResolution.h"
#include "TRD/qaRec/AliTRDcalibration.h"
#endif

#define BIT(n)        (1 << (n))
#define SETBIT(n,i)   ((n) |= BIT(i))
#define TESTBIT(n,i)  ((Bool_t)(((n) & BIT(i)) != 0))
#define CLEARBIT(n,i) ((n) &= ~BIT(i))

Bool_t MEM = kFALSE;
const Int_t fknTasks = 4;
Char_t *fTaskName[fknTasks] = {"Barrel Tracking Effiency", "Combined Tracking Efficiency", "Tracking Resolution", "Calibration"};
enum AliTRDrecoTasks{
  kTrackingEfficiency = 0
  ,kTrackingCombinedEfficiency = 1
  ,kTrackingResolution = 2
  ,kCalibration = 3
};

TChain* CreateESDChain(const char* filename = 0x0, Int_t nfiles=-1 );
void run(const Char_t *files=0x0, Char_t *tasks="ALL", Int_t nmax=-1)
{
  TMemStat *mem = 0x0;
  if(MEM){ 
    gSystem->Load("libMemStat.so");
    gSystem->Load("libMemStatGui.so");
    mem = new TMemStat("new, gnubuildin");
    mem->AddStamp("Start");
  }

  TStopwatch timer;
  timer.Start();

  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libTRDqaRec.so");
  
  Int_t fSteerTask = 0; 
  Bool_t fHasMCdata = kTRUE;
  TObjArray *task = TString(tasks).Tokenize(" ");
  for(Int_t isel = 0; isel < task->GetEntriesFast(); isel++){
    TString s = (dynamic_cast<TObjString *>(task->UncheckedAt(isel)))->String();
    if(s.CompareTo("ALL") == 0){
      for(Int_t itask = 0; itask < fknTasks; itask++) SETBIT(fSteerTask, itask);
      continue;
    } else if(s.CompareTo("EFF") == 0){
      SETBIT(fSteerTask, kTrackingEfficiency);
      continue;
    } else if(s.CompareTo("EFFC") == 0){
      SETBIT(fSteerTask, kTrackingCombinedEfficiency);
      continue;
    } else if(s.CompareTo("RES") == 0){
      SETBIT(fSteerTask, kTrackingResolution);
      continue;
    } else if(s.CompareTo("CAL" ) == 0){
      SETBIT(fSteerTask, kCalibration);
      continue;
    } else if(s.CompareTo("NOMC") == 0){
    	CLEARBIT(fSteerTask, kTrackingEfficiency);
    	CLEARBIT(fSteerTask, kTrackingCombinedEfficiency);
      fHasMCdata = kFALSE;
    } else{
      Info("run.C", Form("Task %s not implemented (yet).", s.Data()));
      continue;
    }
  }
  printf("\n\tRUNNING TRAIN FOR TASKS:\n");
  for(Int_t itask = 0; itask < fknTasks; itask++){
    if(TESTBIT(fSteerTask, itask)) printf("\t%s\n", fTaskName[itask]);
  }

  // define task list pointers;
  AliAnalysisTask *taskPtr[fknTasks];
  memset(taskPtr, 0, fknTasks*sizeof(AliAnalysisTask*));
  Int_t jtask = 0;

  //____________________________________________//
  //gROOT->LoadMacro(Form("%s/TRD/qaRec/CreateESDChain.C", gSystem->ExpandPathName("$ALICE_ROOT")));
  TChain *chain = CreateESDChain(files, nmax);
  //chain->SetBranchStatus("*", 0);
  chain->SetBranchStatus("*FMD*",0);
  chain->SetBranchStatus("*Calo*",0);
  chain->SetBranchStatus("Tracks", 1);
  chain->SetBranchStatus("ESDfriend*",1);
  chain->Lookup();
  chain->GetListOfFiles()->Print();
  printf("\n ----> CHAIN HAS %d ENTRIES <----\n\n", (Int_t)chain->GetEntries());
  
  AliLog::SetGlobalLogLevel(AliLog::kError);

  //____________________________________________
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TRD Track Info Manager");
  //mgr->SetSpecialOutputLocation(source); // To Be Changed
  AliVEventHandler *esdH = 0x0, *mcH = 0x0;
  mgr->SetInputEventHandler(esdH = new AliESDInputHandler);
  if(fHasMCdata) mgr->SetMCtruthEventHandler(mcH = new AliMCEventHandler());
  //mgr->SetDebugLevel(10);

  //____________________________________________
  // TRD track summary generator
  AliTRDtrackInfoGen *task1 = new AliTRDtrackInfoGen();
  task1->SetDebugLevel(1);
  task1->SetMCdata(fHasMCdata);
  mgr->AddTask(task1);
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("data", TChain::Class(), AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("TrackInfoList", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  mgr->ConnectInput( task1, 0, cinput1);
  mgr->ConnectOutput(task1, 0, coutput1);

  //____________________________________________
  // TRD barrel tracking efficiency
  if(TESTBIT(fSteerTask, kTrackingEfficiency)){
    AliTRDtrackingEfficiency *task2 = new AliTRDtrackingEfficiency();
    task2->SetDebugLevel(1);
    mgr->AddTask(task2);
    //Create containers for input/output
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("TrackingEff", TList::Class(), AliAnalysisManager::kOutputContainer, "TRD.TaskTrackingEff.root");
    mgr->ConnectInput( task2, 0, coutput1);
    mgr->ConnectOutput(task2, 0, coutput2);
    taskPtr[jtask++] = task2;
  }

  //____________________________________________
  // TRD combined tracking efficiency
  if(TESTBIT(fSteerTask, kTrackingCombinedEfficiency)){
    AliTRDtrackingEfficiencyCombined *task3 = new AliTRDtrackingEfficiencyCombined();
    task3->SetDebugLevel(0);
    mgr->AddTask(task3);
    // Create containers for input/output
    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("TrackingEffMC", TObjArray::Class(), AliAnalysisManager::kOutputContainer, "TRD.TaskTrackingEffMC.root");
    mgr->ConnectInput( task3, 0, coutput1);
    mgr->ConnectOutput(task3, 0, coutput3);
    taskPtr[jtask++] = task3;
  }

  //____________________________________________
  // TRD tracking resolution
  if(TESTBIT(fSteerTask, kTrackingResolution)){
    AliTRDtrackingResolution *task4 = new AliTRDtrackingResolution();
    task4->SetMCdata(fHasMCdata);
    task4->SetDebugLevel(1);
    mgr->AddTask(task4);
    // Create containers for input/output
    AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("Resolution", TList::Class(), AliAnalysisManager::kOutputContainer, "TRD.TaskResolution.root");
    mgr->ConnectInput( task4, 0, coutput1);
    mgr->ConnectOutput(task4, 0, coutput4);
    taskPtr[jtask++] = task4;
  }

  //____________________________________________
  // TRD calibration
  if(TESTBIT(fSteerTask, kCalibration)){
    AliTRDcalibration *task5 = new AliTRDcalibration();
    task5->SetLow(0);
    task5->SetHigh(30);
    task5->SetDebugLevel(0);
    task5->SetFillZero(kFALSE);
    mgr->AddTask(task5);
    // Create containers for input/output
    AliAnalysisDataContainer *coutput5 = mgr->CreateContainer("Calibration", TList::Class(), AliAnalysisManager::kOutputContainer, "TRD.TaskCalibration.root");
    mgr->ConnectInput(task5, 0, coutput1);
    mgr->ConnectOutput(task5, 0, coutput5);
    taskPtr[jtask++] = task5;
  }
  
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();


  AliCDBManager *cdbManager = AliCDBManager::Instance();
  cdbManager->SetDefaultStorage("local://$ALICE_ROOT");
  //cdbManager->SetSpecificStorage("TRD/Calib/FEE","local:///u/bailhach/aliroot/database30head/database");
  cdbManager->SetRun(0);
  cdbManager->SetCacheFlag(kFALSE);
 
  // initialize TRD settings
  AliMagFMaps *field = new AliMagFMaps("Maps","Maps", 2, 1., 10., AliMagFMaps::k5kG);
  AliTracker::SetFieldMap(field, kTRUE);
  AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
  AliTRDtrackerV1::SetNTimeBins(cal->GetNumberOfTimeBins());


  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();  

  cal->Terminate();
  delete field;
  delete cdbManager;
  for(Int_t it=jtask-1; it>=0; it--) delete taskPtr[it];
  delete task1;
  if(mcH) delete mcH;
  delete esdH;
  delete mgr;
  delete chain;
  if(MEM) delete mem;
  if(MEM) TMemStatViewerGUI::ShowGUI();
}


TChain* CreateESDChain(const char* filename, Int_t nfiles)
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
  while(in.good() && (nfiles--) ) {
    in >> esdfile;
    if (!esdfile.Contains("root")) continue; // protection
    chain->Add(esdfile.Data());
  }

  in.close();

  return chain;
}
