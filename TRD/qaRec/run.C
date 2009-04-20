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
// gSystem->Load("libTRDqaRec.so")
// gSystem->Load("libPWG1.so");
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
#include "TRD/qaRec/AliTRDresolution.h"
#include "TRD/qaRec/AliTRDcalibration.h"
#include "TRD/qaRec/AliTRDalignmentTask.h"
#include "TRD/qaRec/AliTRDpidChecker.h"
#include "TRD/qaRec/AliTRDpidRefMaker.h"
#include "TRD/qaRec/AliTRDcheckDetector.h"
#include "TRD/qaRec/AliTRDclusterResolution.h"
#include "TRD/qaRec/AliTRDmultiplicity.h"


#include "PWG1/AliPerformanceTask.h"
#include "PWG1/AliPerformanceEff.h"
#include "PWG1/AliMCInfoCuts.h"
#include "PWG1/AliRecInfoCuts.h"
#endif

#include "run.h"

Bool_t MEM = kFALSE;
Bool_t fHasMCdata = kTRUE;
Bool_t fHasFriends = kTRUE;
const Int_t kTPCmode = 0; //
const Int_t kTPChpt = 0; //

TChain* MakeChainLST(const char* filename = 0x0);
TChain* MakeChainXML(const char* filename = 0x0);
Int_t   ParseTRD(Char_t *opt);
Int_t   ParseTPC(Char_t *opt);
void run(Char_t *trd="ALL", Char_t *tpc="ALL", const Char_t *files=0x0, Long64_t nev=1234567890, Long64_t first = 0)
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

  AliLog::SetGlobalLogLevel(AliLog::kError);
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

  // Parse TRD options
  Int_t fSteerTRD = ParseTRD(trd);
  Int_t fSteerTPC = ParseTPC(tpc);



  // DEFINE DATA CHAIN
  TChain *chain = 0x0;
  if(!files) chain = MakeChainLST();
  else{
    TString fn(files);
    if(fn.EndsWith("xml")) chain = MakeChainXML(files);
    else chain = MakeChainLST(files);
  }
  if(!chain) return;
  chain->SetBranchStatus("*FMD*",0);
  chain->SetBranchStatus("*Calo*",0);
  chain->SetBranchStatus("Tracks", 1);
  chain->SetBranchStatus("ESDfriend*",1);
  chain->Lookup();
  chain->GetListOfFiles()->Print();
  printf("\n ----> CHAIN HAS %d ENTRIES <----\n\n", (Int_t)chain->GetEntries());
  

  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("Post Reconstruction Calibration/QA");
  AliVEventHandler *esdH = 0x0, *mcH = 0x0;
  mgr->SetInputEventHandler(esdH = new AliESDInputHandler);
  if(fHasMCdata) mgr->SetMCtruthEventHandler(mcH = new AliMCEventHandler());
  //mgr->SetDebugLevel(10);


///////////////////////////////////////////////////////////
///////////////         TRD                     ///////////
///////////////////////////////////////////////////////////
  // define task list pointers for TRD
  AliTRDrecoTask *taskTRD[NTRDTASKS], *task = 0x0;
  memset(taskTRD, 0, NTRDTASKS*sizeof(AliTRDrecoTask*));

  //____________________________________________
  // TRD check ESD
  AliTRDcheckESD *checkESD = new AliTRDcheckESD();
  checkESD->SetMC(fHasMCdata);
  mgr->AddTask(checkESD);
  mgr->ConnectInput(checkESD, 0, mgr->GetCommonInputContainer());  mgr->ConnectOutput(checkESD, 0, mgr->CreateContainer(checkESD->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", checkESD->GetName())));

  //____________________________________________
  // TRD track summary generator
  AliAnalysisDataContainer *coutput1 = 0x0, *coutput1a = 0x0;
	if(TSTBIT(fSteerTRD, kInfoGen)){
    mgr->AddTask(task = new AliTRDtrackInfoGen());
    taskTRD[(Int_t)kInfoGen] = task;
    task->SetDebugLevel(0);
    task->SetMCdata(fHasMCdata);
    mgr->ConnectInput( task, 0, mgr->GetCommonInputContainer());
    coutput1 = mgr->CreateContainer("trackInfo", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
    coutput1a = mgr->CreateContainer("eventInfo", AliTRDeventInfo::Class(), AliAnalysisManager::kExchangeContainer);
    mgr->ConnectOutput(task, 0, coutput1);
    mgr->ConnectOutput(task, 1, coutput1a);
  }

  //____________________________________________
  // TRD detector checker
	if(TSTBIT(fSteerTRD, kCheckDetector)){
    mgr->AddTask(task = new AliTRDcheckDetector());
    taskTRD[(Int_t)kCheckDetector] = task;
    task->SetDebugLevel(0);
    task->SetMCdata(fHasMCdata);
    
    // Create containers for input/output
    mgr->ConnectInput( task, 0, coutput1);
    mgr->ConnectInput( task, 1, coutput1a);
    mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
  }

  //____________________________________________
  // TRD barrel tracking efficiency
  if(fHasMCdata && TSTBIT(fSteerTRD, kTrackingEff)){
    mgr->AddTask(task = new AliTRDtrackingEfficiency());
    taskTRD[(Int_t)kTrackingEff] = task;
    task->SetDebugLevel(0);

    //Create containers for input/output
    mgr->ConnectInput( task, 0, coutput1);
    mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
    
    // TRD single track selection
    if(TSTBIT(fSteerTRD, kMultiplicity)){
      mgr->AddTask(task = new AliTRDmultiplicity());
      taskTRD[(Int_t)kMultiplicity] = task;
      task->SetDebugLevel(0);
      // Create containers for input/output
      mgr->ConnectInput( task, 0, coutput1);
      mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
    }
  }

  //____________________________________________
  // TRD combined tracking efficiency
  if(fHasMCdata && TSTBIT(fSteerTRD, kTrackingEffMC)){
    mgr->AddTask(task = new AliTRDtrackingEfficiencyCombined());
    taskTRD[(Int_t)kTrackingEffMC] = task;
    task->SetDebugLevel(0);

    // Create containers for input/output
    mgr->ConnectInput( task, 0, coutput1);
    mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
  }

  //____________________________________________
  // TRD tracking resolution
  if(TSTBIT(fSteerTRD, kResolution)){
    mgr->AddTask(task = new AliTRDresolution());
    taskTRD[(Int_t)kResolution] = task;
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
    
    // Cluster Error Parameterization
    if(TSTBIT(fSteerTRD, kClErrParam)){
      mgr->AddTask(task = new AliTRDclusterResolution());
      taskTRD[(Int_t)kClErrParam] = task;
      ((AliTRDclusterResolution*)task)->SetExB();
      mgr->ConnectInput(task, 0, co[0]);
      mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
  
      mgr->AddTask(task = new AliTRDclusterResolution("ClErrParamMC"));
      taskTRD[(Int_t)kClErrParam+1] = task;
      ((AliTRDclusterResolution*)task)->SetExB();
      mgr->ConnectInput(task, 0, co[2]);
      mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
    }

    // TRD alignment
    if(TSTBIT(fSteerTRD, kAlignment)){
      mgr->AddTask(task = new AliTRDalignmentTask());
      taskTRD[(Int_t)kAlignment] = task;
      task->SetDebugLevel(0);
  
      // Create containers for input/output
      mgr->ConnectInput(task, 0, coutput1);
      mgr->ConnectOutput(task, 0, mgr->CreateContainer(Form("h%s", task->GetName()), TObjArray::Class(), AliAnalysisManager::kExchangeContainer));
  
      mgr->ConnectOutput(task, 1, mgr->CreateContainer(task->GetName(), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
    }
  }

  //____________________________________________
  // TRD calibration
  if(TSTBIT(fSteerTRD, kCalibration)){
    mgr->AddTask(task = new AliTRDcalibration());
    taskTRD[(Int_t)kCalibration] = task;
    ((AliTRDcalibration*)task)->SetLow(0);
    ((AliTRDcalibration*)task)->SetHigh(30);
    ((AliTRDcalibration*)task)->SetFillZero(kFALSE);
    task->SetDebugLevel(0);

    // Create containers for input/output
    mgr->ConnectInput(task, 0, coutput1);
    mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
  }
  
  
  // TRD PID
  if(TSTBIT(fSteerTRD, kPIDChecker)){
    mgr->AddTask(task = new AliTRDpidChecker());
    taskTRD[(Int_t)kPIDChecker] = task;
    task->SetDebugLevel(0);
    task->SetMCdata(fHasMCdata);
    
    // Create containers for input/output
    mgr->ConnectInput( task, 0, coutput1);
    mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));

    // TRD pid reference 
    if(TSTBIT(fSteerTRD, kPIDRefMaker)){
      mgr->AddTask(task = new AliTRDpidRefMaker());
      taskTRD[(Int_t)kPIDRefMaker] = task;
      task->SetDebugLevel(0);
      task->SetMCdata(fHasMCdata);
      
      // Create containers for input/output
      mgr->ConnectInput( task, 0, coutput1);
      mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", task->GetName())));
      mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("%sNN", task->GetName()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%sNN.root", task->GetName())));
      mgr->ConnectOutput(task, 2, mgr->CreateContainer(Form("%sLQ", task->GetName()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%sLQ.root", task->GetName())));
    }
  }
  printf("\n\tRUNNING TRAIN FOR TASKS:\n");
  for(Int_t itask = 1; itask < NTRDTASKS; itask++){
    if(TSTBIT(fSteerTRD, itask)) printf("\t   %s [%s]\n",  taskTRD[itask]->GetName(), taskTRD[itask]->GetTitle());
  }


/////////////////////////////////////////////////////////
/////////////////     TPC PERFORMANCE      //////////////
/////////////////////////////////////////////////////////
  if(gSystem->Load("libANALYSISalice.so")<0) return;
  if(gSystem->Load("libPWG1.so")<0) return;


  // Create ESD track reconstruction cuts
  AliRecInfoCuts *pRecInfoCuts = new AliRecInfoCuts(); 
  pRecInfoCuts->SetPtRange(0.20,200.0);
  //pRecInfoCuts->SetEtaRange(-0.9,0.9);
  pRecInfoCuts->SetMaxDCAToVertexXY(3.0);
  pRecInfoCuts->SetMaxDCAToVertexZ(3.0);
  pRecInfoCuts->SetMinNClustersTPC(50);
  pRecInfoCuts->SetMinNClustersITS(2);
  pRecInfoCuts->SetMinTPCsignalN(50);
  pRecInfoCuts->SetHistogramsOn(kFALSE); 

  // Create MC track reconstruction cuts
  AliMCInfoCuts  *pMCInfoCuts = new AliMCInfoCuts();
  pMCInfoCuts->SetMinRowsWithDigits(50);
  pMCInfoCuts->SetMaxR(0.025); // from diamond xy size (pp@10TeV) 
  pMCInfoCuts->SetMaxVz(15.);  // from diamond z size  (pp@10TeV)
  pMCInfoCuts->SetRangeTPCSignal(0.5,1.4); 
  pMCInfoCuts->SetMinTrackLength(70);
  
  AliPerformanceObject *taskTPC[NTPCTASKS]; 
  memset(taskTPC, 0, NTPCTASKS*sizeof(AliPerformanceObject*));

  // Create TPC steering task
  AliPerformanceTask *TPC = 0x0;
  if(TSTBIT(fSteerTPC, 0)){
    TPC = new AliPerformanceTask("Performance");
    TPC->SetUseMCInfo(fHasMCdata);
    mgr->AddTask(TPC);
    mgr->ConnectInput(TPC, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(TPC, 0, mgr->CreateContainer("coutput", TList::Class(), AliAnalysisManager::kOutputContainer, Form("TPC.%s.root", TPC->GetName())));
  }

  TClass *ctask = 0x0;
  for(Int_t icomp=1; icomp<NTPCTASKS; icomp++){
    if(!TSTBIT(fSteerTPC, icomp)) continue;
    if(!ctask) ctask = new TClass;
    new(ctask) TClass(fgkTPCtaskClassName[icomp]);
    taskTPC[icomp] = (AliPerformanceObject*)ctask->New();
    taskTPC[icomp]->SetAnalysisMode(kTPCmode);
    taskTPC[icomp]->SetHptGenerator(kTPChpt);
    taskTPC[icomp]->SetAliRecInfoCuts(pRecInfoCuts);
    taskTPC[icomp]->SetAliMCInfoCuts(pMCInfoCuts);
    TPC->AddPerformanceObject(taskTPC[icomp]);
  }
  // verbosity 
  for(Int_t ic = 0; ic < NTPCTASKS; ic++){
    if(taskTPC[ic]) printf("\t   %s [%s]\n",  taskTPC[ic]->GetName(), taskTPC[ic]->GetTitle());
  }


  if (!mgr->InitAnalysis()) return;
  //mgr->PrintStatus();

  mgr->StartAnalysis("local", chain, nev, first);

  timer.Stop();
  timer.Print();  

  cal->Terminate();
  TGeoGlobalMagField::Instance()->SetField(NULL);
  delete cdbManager;
  for(Int_t it=NTRDTASKS; it--; ){ 
    if(!taskTRD[it]) continue;
    printf("Cleaning %s [%s] ...\n", fgkTRDtaskClassName[it], taskTRD[it]->GetTitle());
    delete taskTRD[it];
  }
  delete checkESD;

  for(Int_t it=NTPCTASKS; it--; ){ 
    if(!taskTPC[it]) continue;
    printf("Cleaning %s [%s] ...\n", fgkTPCtaskClassName[it], taskTPC[it]->GetTitle());
    delete taskTPC[it];
  }
  if(TPC) delete TPC;

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


//____________________________________________
Int_t ParseTRD(Char_t *trd)
{
  Int_t fSteerTask = 0;
  TObjArray *tasksArray = TString(trd).Tokenize(" ");
  for(Int_t isel = 0; isel < tasksArray->GetEntriesFast(); isel++){
    TString s = (dynamic_cast<TObjString *>(tasksArray->UncheckedAt(isel)))->String();
    if(s.CompareTo("ALL") == 0){
      for(Int_t itask = 0; itask < NQATASKS; itask++) SETBIT(fSteerTask, itask);
      continue;
    } else if(s.CompareTo("NOFR") == 0){ 
      fHasFriends = kFALSE;
    } else if(s.CompareTo("NOMC") == 0){ 
      fHasMCdata = kFALSE;
    } else { 
      Bool_t foundOpt = kFALSE;  
      for(Int_t itask = 0; itask < NTRDTASKS; itask++){
        if(s.CompareTo(fgkTRDtaskOpt[itask]) != 0) continue;
        SETBIT(fSteerTask, itask); SETBIT(fSteerTask, 0);
        foundOpt = kTRUE;
        break;
      }
      if(!foundOpt) Info("run.C", Form("TRD task %s not implemented (yet).", s.Data()));
    }
  }
  // extra rules for calibration tasks
  if(TSTBIT(fSteerTask, kClErrParam)) SETBIT(fSteerTask, kResolution);
  if(TSTBIT(fSteerTask, kAlignment)) SETBIT(fSteerTask, kResolution);
  if(TSTBIT(fSteerTask, kMultiplicity)) SETBIT(fSteerTask, kTrackingEff);
  if(TSTBIT(fSteerTask, kPIDRefMaker)) SETBIT(fSteerTask, kPIDChecker);

  return fSteerTask;
}



//____________________________________________
Int_t ParseTPC(Char_t *tpc)
{
  Int_t fSteerTask = 0;
  TObjArray *tasksArray = TString(tpc).Tokenize(" ");
  for(Int_t isel = 0; isel < tasksArray->GetEntriesFast(); isel++){
    TString s = (dynamic_cast<TObjString *>(tasksArray->UncheckedAt(isel)))->String();
    if(s.CompareTo("ALL") == 0){
      for(Int_t itask = 0; itask < NTPCPERFORMANCE; itask++) SETBIT(fSteerTask, itask);
      continue;
    } else if(s.CompareTo("NOFR") == 0){ 
      fHasFriends = kFALSE;
    } else if(s.CompareTo("NOMC") == 0){ 
      fHasMCdata = kFALSE;
    } else { 
      Bool_t foundOpt = kFALSE;  
      for(Int_t itask = 0; itask < NTPCTASKS; itask++){
        if(s.CompareTo(fgkTPCtaskOpt[itask]) != 0) continue;
        SETBIT(fSteerTask, itask); SETBIT(fSteerTask, 0);
        foundOpt = kTRUE;
        break;
      }
      if(!foundOpt) Info("run.C", Form("TPC task %s not implemented (yet).", s.Data()));
    }
  }
  return fSteerTask;
}


