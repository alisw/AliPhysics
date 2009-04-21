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
#include "TRD/qaRec/info/AliTRDeventInfo.h"
#include "TRD/qaRec/AliTRDcheckESD.h"
#include "TRD/qaRec/AliTRDinfoGen.h"
#include "TRD/qaRec/AliTRDefficiency.h"
#include "TRD/qaRec/AliTRDefficiencyMC.h"
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
#include "PWG1/AliPerformanceDEdx.h"
#include "PWG1/AliPerformanceTPC.h"
#include "PWG1/AliPerformanceDCA.h"
#include "PWG1/AliPerformanceRes.h"
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



  // VERY GENERAL SETTINGS
  AliLog::SetGlobalLogLevel(AliLog::kError);
  if(gSystem->Load("libANALYSIS.so")<0) return;
  if(gSystem->Load("libANALYSISalice.so")<0) return;


  
  // INITIALIZATION OF RUNNING ENVIRONMENT
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
  


  // BUILD ANALYSIS MANAGER
  AliAnalysisManager *mgr = new AliAnalysisManager("Post Reconstruction Calibration/QA");
  AliVEventHandler *esdH = 0x0, *mcH = 0x0;
  mgr->SetInputEventHandler(esdH = new AliESDInputHandler);
  if(fHasMCdata) mgr->SetMCtruthEventHandler(mcH = new AliMCEventHandler());
  //mgr->SetDebugLevel(10);



///////////////////////////////////////////////////////////
///////////////         TRD                     ///////////
///////////////////////////////////////////////////////////
  // TRD specific library
  if(gSystem->Load("libTRDqaRec.so")<0) return;
  // Parse TRD options
  Int_t fSteerTRD = ParseTRD(trd);
  // TRD data containers
  AliAnalysisDataContainer *ci[] = {0x0, 0x0};
  AliAnalysisDataContainer *co[] = {0x0, 0x0, 0x0, 0x0};

  // plug (set of) TRD wagons in the train
  for(Int_t it=0; it<NQATASKS; it++){
    if(!(TSTBIT(fSteerTRD, it))) continue;
    if(gROOT->LoadMacro(Form("$ALICE_ROOT/TRD/qaRec/macros/Add%s.C+", TString(fgkTRDtaskClassName[it])(3,20).Data()))) {
      Error("run.C", Form("Error loading %s task.", fgkTRDtaskClassName[it]));
      return;
    } 

    switch(it){
    case kCheckESD:
      AddTRDcheckESD(mgr); break;
    case kInfoGen:
      AddTRDinfoGen(mgr, 0x0, ci); break;
    case kCheckDetector:
      AddTRDcheckDetector(mgr, ci, co, fSteerTRD); break;
    case kEfficiency:
      AddTRDefficiency(mgr, ci, co, fSteerTRD); break;
    case kResolution:
      AddTRDresolution(mgr, ci, co, fSteerTRD); break;
    case kPID:
      AddTRDpidChecker(mgr, ci, co, fSteerTRD); break;
    default:
      Warning("run.C", Form("No performance task registered at slot %d.", it)); 
    }
  }


///////////////////////////////////////////////////////////
///////////////         TPC                     ///////////
///////////////////////////////////////////////////////////
  if(gSystem->Load("libPWG1.so")<0) return;
  // Parse TPC options
  Int_t fSteerTPC = ParseTPC(tpc);
  // Create TPC-ESD track reconstruction cuts
  AliRecInfoCuts *pRecInfoCuts = new AliRecInfoCuts(); 
  pRecInfoCuts->SetPtRange(0.20,200.0);
  //pRecInfoCuts->SetEtaRange(-0.9,0.9);
  pRecInfoCuts->SetMaxDCAToVertexXY(3.0);
  pRecInfoCuts->SetMaxDCAToVertexZ(3.0);
  pRecInfoCuts->SetMinNClustersTPC(50);
  pRecInfoCuts->SetMinNClustersITS(2);
  pRecInfoCuts->SetMinTPCsignalN(50);
  pRecInfoCuts->SetHistogramsOn(kFALSE); 
  // Create TPC-MC track reconstruction cuts
  AliMCInfoCuts  *pMCInfoCuts = new AliMCInfoCuts();
  pMCInfoCuts->SetMinRowsWithDigits(50);
  pMCInfoCuts->SetMaxR(0.025); // from diamond xy size (pp@10TeV) 
  pMCInfoCuts->SetMaxVz(15.);  // from diamond z size  (pp@10TeV)
  pMCInfoCuts->SetRangeTPCSignal(0.5,1.4); 
  pMCInfoCuts->SetMinTrackLength(70);
  
  // BUILD STEERING TASK FOR TPC
  if(fSteerTPC){
    if(gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddPerformanceTask.C+")) {
      Error("run.C", "Error loading AliPerformanceTask task.");
      return;
    } 
    AliPerformanceTask *TPC = AddPerformanceTask(mgr);

    // plug (set of) TPC wagons in the train
    TClass ctask; AliPerformanceObject *perf = 0x0;
    for(Int_t icomp=1; icomp<NTPCTASKS; icomp++){
      if(!(TSTBIT(fSteerTPC, icomp))) continue;
      new(&ctask) TClass(fgkTPCtaskClassName[icomp]);
      TPC->AddPerformanceObject((perf = (AliPerformanceObject*)ctask.New()));
      perf->SetAnalysisMode(kTPCmode);
      perf->SetHptGenerator(kTPChpt);
      perf->SetAliRecInfoCuts(pRecInfoCuts);
      perf->SetAliMCInfoCuts(pMCInfoCuts);
    }
  }

  if (!mgr->InitAnalysis()) return;
  // verbosity
  printf("\n\tRUNNING TRAIN FOR TASKS:\n");
  mgr->GetTasks()->ls();
  //mgr->PrintStatus();

  mgr->StartAnalysis("local", chain, nev, first);

  timer.Stop();
  timer.Print();  

  cal->Terminate();
  TGeoGlobalMagField::Instance()->SetField(NULL);
  delete cdbManager;

  // verbosity
  printf("\n\tCLEANING UP TRAIN:\n");
  mgr->GetTasks()->Delete();
//   for(Int_t it=tt->GetEntriesFast(); it--;){
//     if(!(task = (AliAnalysisTask*)tt->At(it))) continue;
//     printf("Cleaning up %s [%s] ...\n", task->GetName(), task->GetTitle());
//     delete task;
//   }

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
  Int_t fSteerTask = 1;
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
  if(TSTBIT(fSteerTask, kCalibration)) SETBIT(fSteerTask, kCheckDetector);
  if(TSTBIT(fSteerTask, kMultiplicity)) SETBIT(fSteerTask, kEfficiency);
  if(TSTBIT(fSteerTask, kEfficiencyMC)) SETBIT(fSteerTask, kEfficiency);
  if(TSTBIT(fSteerTask, kClErrParam)) SETBIT(fSteerTask, kResolution);
  if(TSTBIT(fSteerTask, kAlignment)) SETBIT(fSteerTask, kResolution);
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


