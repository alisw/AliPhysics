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
#include "TRD/qaRec/AliTRDpidChecker.h"
#endif

#define BIT(n)        (1 << (n))
#define SETBIT(n,i)   ((n) |= BIT(i))
#define TESTBIT(n,i)  ((Bool_t)(((n) & BIT(i)) != 0))
#define CLEARBIT(n,i) ((n) &= ~BIT(i))

#include "run.h"

Bool_t MEM = kFALSE;

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
  TObjArray *tasksArray = TString(tasks).Tokenize(" ");
  for(Int_t isel = 0; isel < tasksArray->GetEntriesFast(); isel++){
    TString s = (dynamic_cast<TObjString *>(tasksArray->UncheckedAt(isel)))->String();
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
    } else if(s.CompareTo("PID" ) == 0){
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

  // define task list pointers;
  AliTRDrecoTask *taskPtr[fknTasks], *task = 0x0;
  memset(taskPtr, 0, fknTasks*sizeof(AliAnalysisTask*));

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
  AliAnalysisManager *mgr = new AliAnalysisManager("TRD QA Reconstruction Manager");
  //mgr->SetSpecialOutputLocation(source); // To Be Changed
  AliVEventHandler *esdH = 0x0, *mcH = 0x0;
  mgr->SetInputEventHandler(esdH = new AliESDInputHandler);
  if(fHasMCdata) mgr->SetMCtruthEventHandler(mcH = new AliMCEventHandler());
  //mgr->SetDebugLevel(10);

  //____________________________________________
  // TRD track summary generator
  mgr->AddTask(task = new AliTRDtrackInfoGen());
  taskPtr[(Int_t)kInfoGen] = task;
  task->SetDebugLevel(0);
  task->SetMCdata(fHasMCdata);
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("data", TChain::Class(), AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("TrackInfoList", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  mgr->ConnectInput( task, 0, cinput1);
  mgr->ConnectOutput(task, 0, coutput1);

  //____________________________________________
  // TRD barrel tracking efficiency
  if(TESTBIT(fSteerTask, kTrackingEfficiency)){
    mgr->AddTask(task = new AliTRDtrackingEfficiency());
    taskPtr[(Int_t)kTrackingEfficiency] = task;
    task->SetDebugLevel(0);

    //Create containers for input/output
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName()));
    mgr->ConnectInput( task, 0, coutput1);
    mgr->ConnectOutput(task, 0, coutput2);
  }

  //____________________________________________
  // TRD combined tracking efficiency
  if(TESTBIT(fSteerTask, kTrackingCombinedEfficiency)){
    mgr->AddTask(task = new AliTRDtrackingEfficiencyCombined());
    taskPtr[(Int_t)kTrackingCombinedEfficiency] = task;
    task->SetDebugLevel(0);

    // Create containers for input/output
    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName()));
    mgr->ConnectInput( task, 0, coutput1);
    mgr->ConnectOutput(task, 0, coutput3);
  }

  //____________________________________________
  // TRD tracking resolution
  if(TESTBIT(fSteerTask, kTrackingResolution)){
    mgr->AddTask(task = new AliTRDtrackingResolution());
    taskPtr[(Int_t)kTrackingResolution] = task;
    task->SetMCdata(fHasMCdata);
    task->SetDebugLevel(0);
    
    // Create containers for input/output
    AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName()));
    mgr->ConnectInput( task, 0, coutput1);
    mgr->ConnectOutput(task, 0, coutput4);
  }

  //____________________________________________
  // TRD calibration
  if(TESTBIT(fSteerTask, kCalibration)){
    mgr->AddTask(task = new AliTRDcalibration());
    taskPtr[(Int_t)kCalibration] = task;
    ((AliTRDcalibration*)task)->SetLow(0);
    ((AliTRDcalibration*)task)->SetHigh(30);
    ((AliTRDcalibration*)task)->SetFillZero(kFALSE);
    task->SetDebugLevel(0);

    // Create containers for input/output
    AliAnalysisDataContainer *coutput5 = mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName()));
    mgr->ConnectInput(task, 0, coutput1);
    mgr->ConnectOutput(task, 0, coutput5);
  }
  
  //____________________________________________
  // TRD pid checker
  if(TESTBIT(fSteerTask, kPIDChecker)){
    mgr->AddTask(task = new AliTRDpidChecker());
    taskPtr[(Int_t)kPIDChecker] = task;
    task->SetDebugLevel(0);
    
    // Create containers for input/output
    AliAnalysisDataContainer *coutput6 = mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName()));
    mgr->ConnectInput( task, 0, coutput1);
    mgr->ConnectOutput(task, 0, coutput6);
  }

  if (!mgr->InitAnalysis()) return;
  printf("\n\tRUNNING TRAIN FOR TASKS:\n");
  for(Int_t itask = 1; itask < fknTasks; itask++){
    if(TESTBIT(fSteerTask, itask)) printf("\t   %s [%s]\n", taskPtr[itask]->GetTitle(), taskPtr[itask]->GetName());
  }
  printf("\n\n");
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
  for(Int_t it=fknTasks-1; it>=0; it--) if(taskPtr[it]) delete taskPtr[it];
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
