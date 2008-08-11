// Steer TRD QA train for Reconstruction (Clusterizer, Tracking and PID).
// 
// Usage:
//   run.C(tasks, files, entries)
//   tasks : "ALL" or one/more of the following:
//     "EFF"  : TRD Tracking Efficiency 
//     "EFFC" : TRD Tracking Efficiency Combined (barrel + stand alone) - only in case of simulations
//     "RES"  : TRD tracking Resolution
//     "PID"  : TRD PID - pion efficiency 
//     "PIDR" : TRD PID - reference data
// 
// Authors:
//   Alex Bercuci (A.Bercuci@gsi.de) 
//   Markus Fasel (m.Fasel@gsi.de) 

#define BIT(n)       (1 << (n))
#define SETBIT(n,i)  ((n) |= BIT(i))
#define TESTBIT(n,i) ((Bool_t)(((n) & BIT(i)) != 0))

const Int_t fknTasks = 3;
Char_t *fTaskName[fknTasks] = {"Barrel Tracking Effiency", "Combined Tracking Efficiency", "Tracking Resolution"};
enum AliTRDrecoTasks{
  kTrackingEfficiency = 0
  ,kTrackingCombinedEfficiency = 1
  ,kTrackingResolution = 2
};
void run(Char_t *tasks="ALL", const Char_t *files=0x0, Int_t nmax=-1)
{


  TStopwatch timer;
  timer.Start();

  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libTRDqaRec.so");
  
  Int_t fSteerTask = 0; 
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
    } else if(s.CompareTo("RES" ) == 0){
      SETBIT(fSteerTask, kTrackingResolution);
      continue;
    } else{
      Info("run.C", Form("Task %s not implemented (yet).", s.Data()));
      continue;
    }
  }
  printf("\n\tRUNNING TRAIN FOR TASKS:\n");
  for(itask = 0; itask < fknTasks; itask++){
    if(TESTBIT(fSteerTask, itask)) printf("\t%s\n", fTaskName[itask]);
  }

  //____________________________________________//
  gROOT->LoadMacro(Form("%s/TRD/qaRec/CreateESDChain.C", gSystem->ExpandPathName("$ALICE_ROOT")));
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
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);  
  AliMCEventHandler *mc = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mc);
  //mgr->SetDebugLevel(10);

  //____________________________________________
  // TRD track summary generator
  AliTRDtrackInfoGen *task1 = new AliTRDtrackInfoGen();
  task1->SetDebugLevel(1);
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
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("TrackingEfficiency", TList::Class(), AliAnalysisManager::kOutputContainer, "TRD.TrackingEfficiency.root");
    mgr->ConnectInput( task2, 0, coutput1);
    mgr->ConnectOutput(task2, 0, coutput2);
  }

  //____________________________________________
  // TRD combined tracking efficiency
  if(TESTBIT(fSteerTask, kTrackingCombinedEfficiency)){
    AliTRDtrackingEfficiencyCombined *task3 = new AliTRDtrackingEfficiencyCombined();
    task3->SetDebugLevel(0);
    mgr->AddTask(task3);
    // Create containers for input/output
    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("TrackingEfficiencyCombined", TObjArray::Class(), AliAnalysisManager::kOutputContainer, "TRD.TrackingEfficiencyCombined.root");
    mgr->ConnectInput( task3, 0, coutput1);
    mgr->ConnectOutput(task3, 0, coutput3);
  }

  //____________________________________________
  // TRD combined tracking efficiency
  if(TESTBIT(fSteerTask, kTrackingResolution)){
    AliTRDtrackingResolution *task4 = new AliTRDtrackingResolution();
    task4->SetDebugLevel(1);
    mgr->AddTask(task4);
    // Create containers for input/output
    AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("Tracking Resolution", TList::Class(), AliAnalysisManager::kOutputContainer, "TRD.TrackingResolution.root");
    mgr->ConnectInput( task4, 0, coutput1);
    mgr->ConnectOutput(task4, 0, coutput4);
  }

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();
}
