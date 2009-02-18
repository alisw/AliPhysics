// Steer TRD QA train for Reconstruction (Clusterizer, Tracking and PID).
// 
// Usage:
//   run.C(tasks, files)
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
// gSystem->Load("libNetx.so") ;
// gSystem->Load("libRAliEn.so");
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
#include "TGrid.h"
#include "TAlienCollection.h"
#include "TGridCollection.h"
#include "TGridResult.h"
#include "TGeoGlobalMagField.h"

#include "AliMagF.h"
#include "AliTracker.h"
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h"

#include "TRD/AliTRDtrackerV1.h"
#include "TRD/AliTRDcalibDB.h"
#include "TRD/qaRec/AliTRDtrackInfo/AliTRDeventInfo.h"
#include "TRD/qaRec/AliTRDcheckESD.h"
#include "TRD/qaRec/AliTRDtrackInfoGen.h"
#include "TRD/qaRec/AliTRDtrackingEfficiency.h"
#include "TRD/qaRec/AliTRDtrackingEfficiencyCombined.h"
#include "TRD/qaRec/AliTRDtrackingResolution.h"
#include "TRD/qaRec/AliTRDcalibration.h"
#include "TRD/qaRec/AliTRDalignmentTask.h"
#include "TRD/qaRec/AliTRDpidChecker.h"
#include "TRD/qaRec/AliTRDpidRefMaker.h"
#include "TRD/qaRec/AliTRDcheckDetector.h"
#include "TRD/qaRec/AliTRDclusterResolution.h"
#endif

#include "run.h"

Bool_t MEM = kFALSE;

TChain* MakeChainLST(const char* filename = 0x0);
TChain* MakeChainXML(const char* filename = 0x0);
void run(Char_t *tasks="ALL", const Char_t *files=0x0)
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
  
  // DB INITIALIZATION
  //TODO We should use the GRP if available similar to AliReconstruction::InitGRP()!
  // initialize OCDB manager
  AliCDBManager *cdbManager = AliCDBManager::Instance();
  cdbManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdbManager->SetRun(0);
  cdbManager->SetCacheFlag(kFALSE);
  // initialize magnetic field.
  AliMagF *field = 0x0;
  field = new AliMagF("Maps","Maps", 2, 1., 10., AliMagF::k5kG);
  //field = new AliMagF("Maps","Maps", 2, 0., 10., AliMagF::k2kG);
  TGeoGlobalMagField::Instance()->SetField(field);

  // initialize TRD settings
  AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
  AliTRDtrackerV1::SetNTimeBins(cal->GetNumberOfTimeBins());
  AliGeomManager::LoadGeometry();

  Bool_t fHasMCdata = kTRUE;
  Bool_t fHasFriends = kTRUE;
  TObjArray *tasksArray = TString(tasks).Tokenize(" ");

  Int_t fSteerTask = 0; SETBIT(fSteerTask, kInfoGen);
  for(Int_t isel = 0; isel < tasksArray->GetEntriesFast(); isel++){
    TString s = (dynamic_cast<TObjString *>(tasksArray->UncheckedAt(isel)))->String();
    if(s.CompareTo("ALL") == 0){
      for(Int_t itask = 1; itask < NQATASKS; itask++) SETBIT(fSteerTask, itask);
      continue;
    } else if(s.CompareTo("NOFR") == 0){ 
      fHasFriends = kFALSE;
    } else if(s.CompareTo("NOMC") == 0){ 
      fHasMCdata = kFALSE;
    } else { 
      Bool_t foundOpt = kFALSE;  
      for(Int_t itask = 1; itask < NTRDTASKS; itask++){
        if(s.CompareTo(fgkTRDtaskOpt[itask]) != 0) continue;
        SETBIT(fSteerTask, itask);
        foundOpt = kTRUE;
        break;
      }
      if(!foundOpt) Info("run.C", Form("Task %s not implemented (yet).", s.Data()));
    }
  }
  // extra rules for calibration tasks
  if(TSTBIT(fSteerTask, kClErrParam)) SETBIT(fSteerTask, kResolution);
  if(TSTBIT(fSteerTask, kPIDRefMaker)) SETBIT(fSteerTask, kPIDChecker);
  if(TSTBIT(fSteerTask, kAlignment)) SETBIT(fSteerTask, kResolution);

  // define task list pointers;
  AliTRDrecoTask *taskPtr[NTRDTASKS], *task = 0x0;
  memset(taskPtr, 0, NTRDTASKS*sizeof(AliAnalysisTask*));

  //____________________________________________//
  // DEFINE DATA CHAIN
  TChain *chain = 0x0;
  if(!files) chain = MakeChainLST();
  else{
    TString fn(files);
    if(fn.EndsWith("xml")) chain = MakeChainXML(files);
    else chain = MakeChainLST(files);
  }
  if(!chain) return;

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
  AliAnalysisManager *mgr = new AliAnalysisManager("TRD Reconstruction QA");
  //mgr->SetSpecialOutputLocation(source); // To Be Changed
  AliVEventHandler *esdH = 0x0, *mcH = 0x0;
  mgr->SetInputEventHandler(esdH = new AliESDInputHandler);
  if(fHasMCdata) mgr->SetMCtruthEventHandler(mcH = new AliMCEventHandler());
  //mgr->SetDebugLevel(10);

  //____________________________________________
  // TRD check ESD
  AliTRDcheckESD *checkESD = new AliTRDcheckESD();
  mgr->AddTask(checkESD);
  checkESD->SetMC(fHasMCdata);
  mgr->ConnectInput(checkESD, 0, mgr->GetCommonInputContainer());  mgr->ConnectOutput(checkESD, 0, mgr->CreateContainer(checkESD->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", checkESD->GetName())));

  //____________________________________________
  // TRD track summary generator
  mgr->AddTask(task = new AliTRDtrackInfoGen());
  taskPtr[(Int_t)kInfoGen] = task;
  task->SetDebugLevel(0);
  task->SetMCdata(fHasMCdata);
  mgr->ConnectInput( task, 0, mgr->GetCommonInputContainer());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("trackInfo", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  AliAnalysisDataContainer *coutput1a = mgr->CreateContainer("eventInfo", AliTRDeventInfo::Class(), AliAnalysisManager::kExchangeContainer);
  mgr->ConnectOutput(task, 0, coutput1);
  mgr->ConnectOutput(task, 1, coutput1a);

  //____________________________________________
  // TRD detector checker
	if(TSTBIT(fSteerTask, kCheckDetector)){
    mgr->AddTask(task = new AliTRDcheckDetector());
    taskPtr[(Int_t)kCheckDetector] = task;
    task->SetDebugLevel(0);
    task->SetMCdata(fHasMCdata);
    
    // Create containers for input/output
    mgr->ConnectInput( task, 0, coutput1);
    mgr->ConnectInput( task, 1, coutput1a);
    mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
  }

  //____________________________________________
  // TRD barrel tracking efficiency
  if(fHasMCdata && TSTBIT(fSteerTask, kTrackingEff)){
    mgr->AddTask(task = new AliTRDtrackingEfficiency());
    taskPtr[(Int_t)kTrackingEff] = task;
    task->SetDebugLevel(0);

    //Create containers for input/output
    mgr->ConnectInput( task, 0, coutput1);
    mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
  }

  //____________________________________________
  // TRD combined tracking efficiency
  if(fHasMCdata && TSTBIT(fSteerTask, kTrackingEffMC)){
    mgr->AddTask(task = new AliTRDtrackingEfficiencyCombined());
    taskPtr[(Int_t)kTrackingEffMC] = task;
    task->SetDebugLevel(0);

    // Create containers for input/output
    mgr->ConnectInput( task, 0, coutput1);
    mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
  }

  //____________________________________________
  // TRD tracking resolution
  if(TSTBIT(fSteerTask, kResolution)){
    mgr->AddTask(task = new AliTRDtrackingResolution());
    taskPtr[(Int_t)kResolution] = task;
    task->SetMCdata(fHasMCdata);
    task->SetPostProcess(kFALSE);
    task->SetDebugLevel(0);
    
    // Create containers for input/output
    mgr->ConnectInput( task, 0, coutput1);
    mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));

    // Create output containers for calibration tasks
    const Int_t nc = 4;
    const Char_t *cn[nc] = {"Cl", "Trklt", "MC_Cl", "MC_Trklt"}; 
    AliAnalysisDataContainer *co[nc]; 
    for(Int_t ic = 0; ic<nc; ic++){
      co[ic] = mgr->CreateContainer(Form("%s%s", task->GetName(), cn[ic]), TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
      mgr->ConnectOutput(task, 1+ic, co[ic]);
    }
    
    // test reconstruction calibration plugin
    if(TSTBIT(fSteerTask, kClErrParam)){
      mgr->AddTask(task = new AliTRDclusterResolution());
      taskPtr[(Int_t)kClErrParam] = task;
      ((AliTRDclusterResolution*)task)->SetExB();
      mgr->ConnectInput(task, 0, co[0]);
      mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
  
      mgr->AddTask(task = new AliTRDclusterResolution("ClErrParamMC"));
      taskPtr[(Int_t)kClErrParam+1] = task;
      ((AliTRDclusterResolution*)task)->SetExB();
      mgr->ConnectInput(task, 0, co[2]);
      mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
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
  // TRD PID
  if(TSTBIT(fSteerTask, kPIDChecker)){
    mgr->AddTask(task = new AliTRDpidChecker());
    taskPtr[(Int_t)kPIDChecker] = task;
    task->SetDebugLevel(0);
    task->SetMCdata(fHasMCdata);
    
    // Create containers for input/output
    mgr->ConnectInput( task, 0, coutput1);
    mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));

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
  }


  if (!mgr->InitAnalysis()) return;
  printf("\n\tRUNNING TRAIN FOR TASKS:\n");
  for(Int_t itask = 1; itask < NTRDTASKS; itask++){
    if(TSTBIT(fSteerTask, itask)) printf("\t   %s [%s]\n",  taskPtr[itask]->GetName(), taskPtr[itask]->GetTitle());
  }
  printf("\n\n");
  //mgr->PrintStatus();

  mgr->StartAnalysis("local", chain);

  timer.Stop();
  timer.Print();  

  cal->Terminate();
  TGeoGlobalMagField::Instance()->SetField(NULL);
  delete cdbManager;
  for(Int_t it=NTRDTASKS; it--; ){ 
    if(taskPtr[it]){ 
      printf("Cleaning %s [%s] ...\n", fgkTRDtaskClassName[it], taskPtr[it]->GetTitle());
      delete taskPtr[it];
    }
  }
  if(mcH) delete mcH;
  delete esdH;
  delete mgr;
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
