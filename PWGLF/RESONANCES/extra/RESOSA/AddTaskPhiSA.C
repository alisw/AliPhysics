
#ifdef __ECLIPSE_IDE
//  few includes and external declarations just for the IDE
#include "Riostream.h"
#include "TSystem.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliAnalysisManager.h"
AliAnalysisDataContainer* AddTaskFlowQnVectorCorrections();
#endif // ifdef __ECLIPSE_IDE declaration and includes for the ECLIPSE IDE


#ifdef __CLING__
// ROOT6 macro inclusion
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/macros/runAnalysis.H>
#include <PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/macros/loadRunOptions.C>
#include <PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/macros/AddTaskFlowQnVectorCorrections.C>
#include <PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/macros/AddTaskQnVectorAnalysis.C>
#endif


using std::cout;
using std::endl;



void AddTaskPhiSA(char *analysislevel, const Bool_t useShift, Bool_t bMCtruth, char* runtype, Float_t pairrapidity, char *systematiccut, Int_t dataperiod, Int_t sw, Int_t MT)

{
  // standard with task
  printf("===================================================================================\n");
  printf("\n                PID: Initialising ADDTASKPHISA                            \n");
  printf("===================================================================================\n");



  AliAnalysisManager *mgr2 = AliAnalysisManager::GetAnalysisManager();  
  const char *taskname = Form("NewTpcTofPidTaskKStarWt%d",20);  
  std::cout<<"The name of task is:"<<taskname<<"****************"<<std::endl;
  


  if (sw==0)
#ifndef __CLING__
  gROOT->LoadMacro("AliAnalysisTaskPhiSA.cxx++g");                                                                                         
#endif                                                                                        
  else if (sw==1)
#ifndef __CLING__
  gROOT->LoadMacro("AliAnalysisTaskPhiSA.cxx"); 
#endif                                                                                        
  AliAnalysisTaskPhiSA* task = new AliAnalysisTaskPhiSA(taskname,useShift);
  task->SetAnalysisLevel(analysislevel);
  task->CheckPileUp(kFALSE);
  task->AcceptTPCVertex(kTRUE);
  task->SetPOIAndRPTrackType("GLOBAL","TPC");
  task->SetPairRapidityCut(pairrapidity);
  task->SetSystemaicCutType(systematiccut);
  



  gROOT->SetStyle("Plain");                                                                                                                  
  //const char *configpath = ".";
  const char *configpath = "alien:///alice/cern.ch/user/p/prottay/checktxt";                  
  #ifndef __CLING__
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/extra/RESOSA/runAnalysis.H");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/extra/RESOSA/loadRunOptions.C");
  //gROOT->LoadMacro("runAnalysis.H");
  //gROOT->LoadMacro("loadRunOptions.C");
  #endif                                                                                        

  if (!loadRunOptions(kFALSE, configpath,dataperiod)) {
    std::cout << "ERROR: configuration options not loaded. ABORTING!!!" << std::endl;
    return;
  }


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
      std::cout<<" WARNING: the file <shift.root> with sin and cosine averages from the previous run was not available."<<std::endl;
      break;
    }
  }





  if(useShift){
    TString outfilename = Form("AnalysisResults.root");
    std::cout << "Name of the output file    : " << outfilename.Data() << std::endl;
  }
  else { TString outfilename = Form("AnalysisResults.root");
    std::cout << "Name of the output file    : " << outfilename.Data() << std::endl;
  }

  mgr2->AddTask(task);
  
 
  #ifndef __CLING__
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/macros/AddTaskFlowQnVectorCorrections.C");
#endif
  AliAnalysisDataContainer *corrTask = AddTaskFlowQnVectorCorrections();




  if (bRunQnVectorAnalysisTask) {


    printf("===================================================================================\n");
    printf("\n                Hi i am in RunQnVectorAnalysisTask                            \n");
    printf("===================================================================================\n");



#ifndef __CLING__
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/macros/AddTaskQnVectorAnalysis.C");
#endif
    AliAnalysisTaskQnVectorAnalysis* taskQn = AddTaskQnVectorAnalysis(bUseMultiplicity, b2015DataSet);
    taskQn->SetExpectedCorrectionPass(szCorrectionPass.Data());
    taskQn->SetAlternativeCorrectionPass(szAltCorrectionPass.Data());

    mgr2->AddTask(taskQn);


    //create output container                                                                                                                 
    AliAnalysisDataContainer *cOutputQnAnaEventQA =
      mgr2->CreateContainer("QnAnalysisEventQA",
			    TList::Class(),
			    AliAnalysisManager::kOutputContainer,
			    "QnAnalysisEventQA.root");
    AliAnalysisDataContainer *coutput1 = mgr2->CreateContainer("FlowOut", TList::Class(), AliAnalysisManager::kOutputContainer, outfilename);


    mgr2->ConnectInput(taskQn,  0, mgr2->GetCommonInputContainer());
    mgr2->ConnectInput(taskQn,  1, corrTask);
    mgr2->ConnectOutput(taskQn, 1, cOutputQnAnaEventQA );
    if (MT==1)
      {
	std::cout<<"*****************hi i am inside MT loop*************************"<<std::endl;
    mgr2->ConnectInput(task,  0, mgr2->GetCommonInputContainer());
    mgr2->ConnectOutput(task, 1, coutput1 );
      }

  }



  else

  {


 
    printf("===================================================================================\n");
    printf("\n                Hi i am not in RunQnVectorAnalysisTask                            \n");
    printf("===================================================================================\n");


  AliAnalysisDataContainer *cinput = mgr2->GetCommonInputContainer();                                                                        
  if(useShift) {                                                                                             

    AliAnalysisDataContainer *cinputShift = mgr->CreateContainer("avShift",TList::Class(),AliAnalysisManager::kInputContainer);
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("ShiftFlowOut", TList::Class(), AliAnalysisManager::kOutputContainer, outfilename); 

  }

else
  {
  AliAnalysisDataContainer *coutput1 = mgr2->CreateContainer("AnOut", TList::Class(), AliAnalysisManager::kOutputContainer, outfilename);  
  }



  
  // connect input/output                                                                                                                     
  mgr2->ConnectInput(task, 0, cinput);                                                                                                        
  if(useShift) {
    mgr2->ConnectInput(task,1,cinputShift);
    cinputShift->SetData(shiftList);
  }     
  mgr2->ConnectOutput(task, 1, coutput1);                                                                                                    

}
}
