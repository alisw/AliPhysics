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

#include "TRD/qaRec/macros/AddTRDcheckESD.C"
#include "TRD/qaRec/macros/AddTRDinfoGen.C"
#include "TRD/qaRec/macros/AddTRDcheckDET.C"
#include "TRD/qaRec/macros/AddTRDefficiency.C"
#include "TRD/qaRec/macros/AddTRDresolution.C"
#include "TRD/qaRec/macros/AddTRDcheckPID.C"

#include "../../PWG1/macros/AddPerformanceTask.C"
#endif

#include "macros/AliTRDperformanceTrain.h"
#include "../../PWG1/macros/AddPerformanceTask.h"


Bool_t MEM = kFALSE;
Bool_t fHasMCdata = kTRUE;
Bool_t fHasFriends = kTRUE;

TChain* MakeChainLST(const char* filename = 0x0);
TChain* MakeChainXML(const char* filename = 0x0);
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

/*    } else if(s.CompareTo("NOFR") == 0){ 
      fHasFriends = kFALSE;
    } else if(s.CompareTo("NOMC") == 0){ 
      fHasMCdata = kFALSE;
*/
  
  // INITIALIZATION OF RUNNING ENVIRONMENT
  //TODO We should use the GRP if available similar to AliReconstruction::InitGRP()!
  // initialize OCDB manager
  AliCDBManager *cdbManager = AliCDBManager::Instance();
  cdbManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdbManager->SetRun(0);
  cdbManager->SetCacheFlag(kFALSE);
  // initialize magnetic field.
  if(!TGeoGlobalMagField::Instance()->GetField()){
    TGeoGlobalMagField::Instance()->SetField(
      new AliMagF("Maps","Maps", 2, 1., 10., AliMagF::k5kG)
      //AliMagF("Maps","Maps", 2, 0., 10., AliMagF::k2kG);
    );
    AliGeomManager::LoadGeometry();
  }


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
  // TRD data containers
  AliAnalysisDataContainer *ci[] = {0x0, 0x0};


  // initialize TRD settings
  AliTRDcalibDB *cal = AliTRDcalibDB::Instance();
  AliTRDtrackerV1::SetNTimeBins(cal->GetNumberOfTimeBins());

  // plug (set of) TRD wagons in the train
  if(trd){
    for(Int_t it=0; it<NTRDQATASKS; it++){
      if(gROOT->LoadMacro(Form("$ALICE_ROOT/TRD/qaRec/macros/Add%s.C+", TString(fgkTRDtaskClassName[it])(3,20).Data()))) {
        Error("run.C", Form("Error loading %s task.", fgkTRDtaskClassName[it]));
        return;
      } 
  
      switch(it){
      case kCheckESD:
        AddTRDcheckESD(mgr); break;
      case kInfoGen:
        AddTRDinfoGen(mgr, trd, 0x0, ci); break;
      case kCheckDET:
        AddTRDcheckDET(mgr, trd, ci); break;
      case kEfficiency:
        AddTRDefficiency(mgr, trd, ci); break;
      case kResolution:
        AddTRDresolution(mgr, trd, ci); break;
      case kCheckPID:
        AddTRDcheckPID(mgr, trd, ci); break;
      default:
        Warning("run.C", Form("No performance task registered at slot %d.", it)); 
      }
    }
  }

///////////////////////////////////////////////////////////
///////////////         TPC                     ///////////
///////////////////////////////////////////////////////////
  if(gSystem->Load("libPWG1.so")<0) return;
  
  // BUILD STEERING TASK FOR TPC
  if(tpc){
    if(gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddPerformanceTask.C+")) {
      Error("run.C", "Error loading AliPerformanceTask task.");
      return;
    } 
    AddPerformanceTask(mgr, tpc);
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
