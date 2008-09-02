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
// Authors:
//   Alex Bercuci (A.Bercuci@gsi.de) 
//   Markus Fasel (m.Fasel@gsi.de) 

#define BIT(n)        (1 << (n))
#define SETBIT(n,i)   ((n) |= BIT(i))
#define TESTBIT(n,i)  ((Bool_t)(((n) & BIT(i)) != 0))
#define CLEARBIT(n,i) ((n) &= ~BIT(i))

const Int_t fknTasks = 4;
Char_t *fTaskName[fknTasks] = {"Barrel Tracking Effiency", "Combined Tracking Efficiency", "Tracking Resolution", "Calibration"};
enum AliTRDrecoTasks{
  kTrackingEfficiency = 0
  ,kTrackingCombinedEfficiency = 1
  ,kTrackingResolution = 2
  ,kCalibration = 3
};
void run(const Char_t *files=0x0, Char_t *tasks="ALL", Int_t nmax=-1)
{


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
  if(fHasMCdata) mgr->SetMCtruthEventHandler(new AliMCEventHandler());
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
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("Efficiency", TList::Class(), AliAnalysisManager::kOutputContainer, "TRD.TrackingEfficiency.root");
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
    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("Efficiency2", TObjArray::Class(), AliAnalysisManager::kOutputContainer, "TRD.TrackingEfficiencyCombined.root");
    mgr->ConnectInput( task3, 0, coutput1);
    mgr->ConnectOutput(task3, 0, coutput3);
  }

  //____________________________________________
  // TRD tracking resolution
  if(TESTBIT(fSteerTask, kTrackingResolution)){
    AliTRDtrackingResolution *task4 = new AliTRDtrackingResolution();
    task4->SetMCdata(fHasMCdata);
    task4->SetDebugLevel(1);
    mgr->AddTask(task4);
    // Create containers for input/output
    AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("Resolution", TList::Class(), AliAnalysisManager::kOutputContainer, "TRD.TrackingResolution.root");
    mgr->ConnectInput( task4, 0, coutput1);
    mgr->ConnectOutput(task4, 0, coutput4);
  }

  //____________________________________________
  // TRD calibration
  if(TESTBIT(fSteerTask, kCalibration)){
    AliTRDcalib *task5 = new AliTRDcalib();
    task5->SetLow(0);
    task5->SetHigh(30);
    task5->SetDebugLevel(0);
    task5->SetFillZero(kFALSE);
    mgr->AddTask(task5);
    // Create containers for input/output
    AliAnalysisDataContainer *coutput5 = mgr->CreateContainer("Calibration", TList::Class(), AliAnalysisManager::kOutputContainer, "TRD.Calibration.root");
    mgr->ConnectInput(task5,0,cinput1);
    mgr->ConnectOutput(task5,0,coutput5);
  }
  
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  
  // initialize TRD settings
  AliTRDtrackerV1::SetNTimeBins(24);
  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();
}
