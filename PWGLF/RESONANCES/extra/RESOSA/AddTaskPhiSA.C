
#ifdef __ECLIPSE_IDE
//  few includes and external declarations just for the IDE
#include "Riostream.h"
#include "TSystem.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliAnalysisManager.h"
#include "runAnalysis.H"

#endif // ifdef __ECLIPSE_IDE declaration and includes for the ECLIPSE IDE



void AddTaskPhiSA(Int_t RunNo, char *analysislevel, Bool_t useShift, Bool_t bMCtruth, char* runtype, Float_t pairrapidity, char *systematiccut, Int_t dataperiod)

{
  // standard with task
  printf("===================================================================================\n");
  printf("\n                PID: Initialising ADDTASKPHISA                            \n");
  printf("===================================================================================\n");



  gSystem->Load("libCore.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libEventMixing.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGLFresonances.so");
  gSystem->Load("libPWGPPevcharQn.so");
  gSystem->Load("libPWGPPevcharQnInterface.so");







  gROOT->SetStyle("Plain");
  //const char *configpath = ".";
  const char *configpath = "alien:///alice/cern.ch/user/p/prottay/checktxt";

  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/extra/RESOSA/runAnalysis.H");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/extra/RESOSA/loadRunOptions.C");

    
  if (!loadRunOptions(kFALSE, configpath,dataperiod)) {
    cout << "ERROR: configuration options not loaded. ABORTING!!!" << endl;
    return;
  }
  
  // sets name of grid generated macros
  const char *taskname = Form("NewTpcTofPidTaskPhiWt%d",RunNo);
  TString addincpath = "-I$ALICE_PHYSICS/include -I$ALICE_ROOT/include";

  // add aliroot indlude path
  gSystem->AddIncludePath(addincpath.Data());
  
  
  //gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");


  AliAnalysisGrid      *alienHandler   =   NULL;
 
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"); // Needed for LHC2015o                        
  AliMultSelectionTask *taskM = AddTaskMultSelection(kFALSE);            // kFALSE == User mode, kTRUE == Calibration mode                  
  taskM->SetSelectedTriggerClass(AliVEvent::kINT7);
  mgr->AddTask(taskM);






  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/macros/AddTaskFlowQnVectorCorrections.C");
  AliAnalysisDataContainer *corrTask = AddTaskFlowQnVectorCorrections();

  cout<<"*******************************************************************"<<endl;
  
  



  
  // create task
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/extra/RESOSA/AliAnalysisTaskPhiSA.cxx++g");
  
  
  AliAnalysisTaskPhiSA* task = new AliAnalysisTaskPhiSA(taskname,RunNo,useShift);
  task->SetAnalysisLevel(analysislevel);
  task->CheckPileUp(kFALSE);
  task->AcceptTPCVertex(kTRUE);
  task->SetPOIAndRPTrackType("GLOBAL","TPC");
  task->SetPairRapidityCut(pairrapidity);
  task->SetSystemaicCutType(systematiccut);  
  // Open external input files
  //=====================================================================

  //shift correction
  if(useShift) {
    TFile *shiftFile = NULL;
    TList *shiftList = NULL;
    task->SetUseShifting(kTRUE);
    //open the file with the weights:
    shiftFile = TFile::Open(Form("/Users/ranbirsingh/SpinAlignment/EPlane/phi/systematic/data/shiftFiles/shift_%d.root",RunNo),"READ");
    if(shiftFile) {
      //access the list which holds the profile with averages:
      shiftList = (TList*)shiftFile->Get("avShift");
      task->SetShiftList(shiftList);
    }
    else {
      cout<<" WARNING: the file <shift.root> with sin and cosine averages from the previous run was not available."<<endl;
      break;
    }
  }
  
  mgr->AddTask(task);

  if(useShift){
    TString outfilename = Form("AnalysisResults.root");
    cout << "Name of the output file    : " << outfilename.Data() << endl;
  }
  else { TString outfilename = Form("AnalysisResults.root");
    cout << "Name of the output file    : " << outfilename.Data() << endl;
  }

  //Connect input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  
  if(useShift) {    
    AliAnalysisDataContainer *cinputShift = mgr->CreateContainer("avShift",TList::Class(),AliAnalysisManager::kInputContainer); 
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("ShiftFlowOut", TList::Class(), AliAnalysisManager::kOutputContainer, outfilename);
  } else {
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("FlowOut", TList::Class(), AliAnalysisManager::kOutputContainer, outfilename);
  }

  // connect input/output
  mgr->ConnectInput(task, 0, cinput);

  if(useShift) {
    mgr->ConnectInput(task,1,cinputShift);
    cinputShift->SetData(shiftList);
  }
  mgr->ConnectOutput(task, 1, coutput1);
  
  
}
  
