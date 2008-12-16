// Steer TRD QA train for Reconstruction (Clusterizer, Tracking and PID).
// 
// Usage:
//   run.C(tasks, files, entries)
//   tasks : "ALL" or one/more of the following:
//     "EFF"  : TRD Tracking Efficiency 
//     "EFFC" : TRD Tracking Efficiency Combined (barrel + stand alone) - only in case of simulations
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
#include "TRD/qaRec/AliTRDtrackInfo/AliTRDeventInfo.h"
#include "TRD/qaRec/AliTRDtrackInfoGen.h"
#include "TRD/qaRec/AliTRDtrackingEfficiency.h"
#include "TRD/qaRec/AliTRDtrackingEfficiencyCombined.h"
#include "TRD/qaRec/AliTRDtrackingResolution.h"
#include "TRD/qaRec/AliTRDcalibration.h"
#include "TRD/qaRec/AliTRDpidChecker.h"
#include "TRD/qaRec/AliTRDpidRefMaker.h"
#include "TRD/qaRec/AliTRDcheckDetector.h"
#include "TRD/qaRec/AliTRDclusterResolution.h"
#endif

#include "run.h"

Bool_t MEM = kFALSE;

TChain* CreateESDChain(const char* filename = 0x0, Int_t nfiles=-1 );
void run(Char_t *tasks="ALL", const Char_t *files=0x0, Int_t nmax=-1)
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

  if(gSystem->Load("libANALYSIS.so")<0) return;
  if(gSystem->Load("libTRDqaRec.so")<0) return;
  
  Bool_t fHasMCdata = kTRUE;
  Bool_t fHasFriends = kTRUE;
  TObjArray *tasksArray = TString(tasks).Tokenize(" ");

  Int_t fSteerTask = 0; SETBIT(fSteerTask, kInfoGen);
  for(Int_t isel = 0; isel < tasksArray->GetEntriesFast(); isel++){
    TString s = (dynamic_cast<TObjString *>(tasksArray->UncheckedAt(isel)))->String();
    if(s.CompareTo("ALL") == 0){
      for(Int_t itask = 1; itask <= NTRDTASKS; itask++) SETBIT(fSteerTask, itask);
      continue;
    } else if(s.CompareTo("NOFR") == 0){ 
      fHasFriends = kFALSE;
    } else if(s.CompareTo("NOMC") == 0){ 
      fHasMCdata = kFALSE;
    } else { 
      Bool_t foundOpt = kFALSE;  
      for(Int_t itask = 1; itask <= NTRDTASKS; itask++){
        if(s.CompareTo(fgkTRDtaskOpt[itask]) != 0) continue;
        SETBIT(fSteerTask, itask);
        foundOpt = kTRUE;
        break;
      }
      if(!foundOpt) Info("run.C", Form("Task %s not implemented (yet).", s.Data()));
    }
  }
  // extra rules for calibration tasks
  if(TSTBIT(fSteerTask, kClusterErrorParam)) SETBIT(fSteerTask, kTrackingResolution);

  // define task list pointers;
  AliTRDrecoTask *taskPtr[2*NTRDTASKS], *task = 0x0;
  memset(taskPtr, 0, 2*NTRDTASKS*sizeof(AliAnalysisTask*));

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
  task->SetDebugLevel(1);
  task->SetMCdata(fHasMCdata);
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("data", TChain::Class(), AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("trackInfo", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  AliAnalysisDataContainer *coutput1a = mgr->CreateContainer("eventInfo", AliTRDeventInfo::Class(), AliAnalysisManager::kExchangeContainer);
  mgr->ConnectInput( task, 0, cinput1);
  mgr->ConnectOutput(task, 0, coutput1);
  mgr->ConnectOutput(task, 1, coutput1a);

  //____________________________________________
  // TRD detector checker
	if(TSTBIT(fSteerTask, kCheckDetector)){
    mgr->AddTask(task = new AliTRDcheckDetector());
    taskPtr[(Int_t)kCheckDetector] = task;
    task->SetDebugLevel(4);
    task->SetMCdata(fHasMCdata);
    
    // Create containers for input/output
    mgr->ConnectInput( task, 0, coutput1);
    mgr->ConnectInput( task, 1, coutput1a);
    mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
  }

  //____________________________________________
  // TRD barrel tracking efficiency
  if(fHasMCdata && TSTBIT(fSteerTask, kTrackingEfficiency)){
    mgr->AddTask(task = new AliTRDtrackingEfficiency());
    taskPtr[(Int_t)kTrackingEfficiency] = task;
    task->SetDebugLevel(0);

    //Create containers for input/output
    mgr->ConnectInput( task, 0, coutput1);
    mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
  }

  //____________________________________________
  // TRD combined tracking efficiency
  if(fHasMCdata && TSTBIT(fSteerTask, kTrackingCombinedEfficiency)){
    mgr->AddTask(task = new AliTRDtrackingEfficiencyCombined());
    taskPtr[(Int_t)kTrackingCombinedEfficiency] = task;
    task->SetDebugLevel(0);

    // Create containers for input/output
    mgr->ConnectInput( task, 0, coutput1);
    mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
  }

  //____________________________________________
  // TRD tracking resolution
  if(TSTBIT(fSteerTask, kTrackingResolution)){
    mgr->AddTask(task = new AliTRDtrackingResolution());
    taskPtr[(Int_t)kTrackingResolution] = task;
    task->SetMCdata(fHasMCdata);
    task->SetPostProcess(kFALSE);
    task->SetDebugLevel(1);
    
    // Create containers for input/output
    mgr->ConnectInput( task, 0, coutput1);
    mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));

    // Create output containers for calibration tasks
    const Int_t nc = 5;
    const Char_t *cn[nc] = {"ClRez", "TrkltRez", "TrkltPhiRez", "ClRes", "TrkltRes"}; 
    AliAnalysisDataContainer *co[nc]; 
    for(Int_t ic = 0; ic<nc; ic++){
      co[ic] = mgr->CreateContainer(Form("%s%s", task->GetName(), cn[ic]), TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
      mgr->ConnectOutput(task, 1+ic, co[ic]);
    }
    
    // test reconstruction calibration plugin
    if(TSTBIT(fSteerTask, kClusterErrorParam)){
      mgr->AddTask(task = new AliTRDclusterResolution());
      taskPtr[(Int_t)kClusterErrorParam] = task;
      mgr->ConnectInput(task, 0, co[0]);
      mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
  
      mgr->AddTask(task = new AliTRDclusterResolution());
      taskPtr[(Int_t)kClusterErrorParam+1] = task;
      mgr->ConnectInput(task, 0, co[1]);
      mgr->ConnectOutput(task, 0, mgr->CreateContainer(Form("%sMC", task->GetName()), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%sMC.root", task->GetName())));
    }
  }

  //____________________________________________
  // TRD calibration
  if(TSTBIT(fSteerTask, kCalibration)){
    mgr->AddTask(task = new AliTRDcalibration());
    taskPtr[(Int_t)kCalibration] = task;
    ((AliTRDcalibration*)task)->SetLow(0);
    ((AliTRDcalibration*)task)->SetHigh(30);
    ((AliTRDcalibration*)task)->SetFillZero(kFALSE);
    task->SetDebugLevel(0);

    // Create containers for input/output
    mgr->ConnectInput(task, 0, coutput1);
    mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
  }
  
  //____________________________________________
  // TRD alignment
  if(TSTBIT(fSteerTask, kAlignment)){
    mgr->AddTask(task = new AliTRDalignmentTask());
    taskPtr[(Int_t)kAlignment] = task;
    task->SetDebugLevel(0);

    // Create containers for input/output
    mgr->ConnectInput(task, 0, coutput1);
    mgr->ConnectOutput(task, 0, mgr->CreateContainer(Form("h%s", task->GetName()), TObjArray::Class(), AliAnalysisManager::kExchangeContainer));

    mgr->ConnectOutput(task, 1, mgr->CreateContainer(task->GetName(), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
  }
  
  //____________________________________________
  // TRD pid checker
  if(TSTBIT(fSteerTask, kPIDChecker)){
    mgr->AddTask(task = new AliTRDpidChecker());
    taskPtr[(Int_t)kPIDChecker] = task;
    task->SetDebugLevel(0);
    task->SetMCdata(fHasMCdata);
    
    // Create containers for input/output
    mgr->ConnectInput( task, 0, coutput1);
    mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
  }


  //____________________________________________
  // TRD pid reference 
  if(TSTBIT(fSteerTask, kPIDRefMaker)){
    mgr->AddTask(task = new AliTRDpidRefMaker());
    taskPtr[(Int_t)kPIDRefMaker] = task;
    task->SetDebugLevel(0);
    task->SetMCdata(fHasMCdata);
    
    // Create containers for input/output
    mgr->ConnectInput( task, 0, coutput1);
    mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
    mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("%sNN", task->GetName()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%sNN.root", task->GetName())));
    mgr->ConnectOutput(task, 2, mgr->CreateContainer(Form("%sLQ", task->GetName()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%sLQ.root", task->GetName())));
  }


  if (!mgr->InitAnalysis()) return;
  printf("\n\tRUNNING TRAIN FOR TASKS:\n");
  for(Int_t itask = 1; itask <= NTRDTASKS; itask++){
    if(TSTBIT(fSteerTask, itask)) printf("\t   %s [%s]\n", taskPtr[itask]->GetTitle(), taskPtr[itask]->GetName());
  }
  printf("\n\n");
  //mgr->PrintStatus();


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
  AliGeomManager::LoadGeometry();


  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();  

  cal->Terminate();
  delete field;
  delete cdbManager;
  for(Int_t it=2*NTRDTASKS; it--; ){ 
    if(taskPtr[it]) delete taskPtr[it];
  }
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
