#if defined(__CLING__)
R__ADD_INCLUDE_PATH($ALICE_ROOT)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
R__ADD_INCLUDE_PATH(./)
#include "AliAODInputHandler.h"
#include "ANALYSIS/macros/AddTaskPIDResponse.C"
#include <OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C>
#include <PWGGA/PHOSTasks/PHOS_PbPb/AddAODPHOSTender.C>
#include "AliAnalysisPHOSFluctuations.h"
R__LOAD_LIBRARY(AliAnalysisPHOSFluctuations_cxx.so)
#endif
void Fluct(const char* dataset="collection.xml")
{

  gROOT->LoadMacro("AliAnalysisPHOSFluctuations_cxx.so");

  
  // Create the chain
  TChain* chain = new TChain("aodTree");

  // Connect to alien
  TGrid::Connect("alien://");

//   cout << "Pi0Analysis: processing collection " << dataset << endl;
  cout << " processing collection " << dataset << endl;
  

  TXMLEngine a ;
  XMLDocPointer_t doc = a.ParseFile(dataset);
  XMLNodePointer_t rn= a.DocGetRootElement (doc) ;
  XMLNodePointer_t nn= a.GetChild(rn) ;
  XMLNodePointer_t nn2= a.GetChild(nn) ;
  while(nn2){
    XMLNodePointer_t nn3= a.GetChild(nn2) ;
    const char * rawFile = a.GetAttr(nn3,"turl") ;
    printf("Processing %s\n", rawFile) ;
    if(!rawFile) break ;
    chain->Add(rawFile);
    printf("Chain: %lld entries.\n",chain->GetEntries());
    nn2= a.GetNext(nn2) ;
  }

  // chain->Add("AliAOD.root") ;
    

  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("Pi0Spectrum");
  
  // ESD input handler
  AliAODInputHandler* esdH = new AliAODInputHandler();
  mgr->SetInputEventHandler(esdH);
  
  AddTaskMultSelection() ;

  
  Bool_t isMC=kTRUE; // kTRUE in case of MC
  Bool_t cachePIDResponse                     = kFALSE;
  TString tenderPassData                          = "2";
  Int_t customSplinesPass                         = -1;
  TString PIDresponseSetting                    = "";
  TString splinesYear                                  = "";

//  AddTaskPIDResponse(kTRUE,kFALSE,kTRUE,tenderPassData,kFALSE,PIDresponseSetting, kTRUE, kTRUE, customSplinesPass ) ; 
//   AddTaskPIDResponse(kFALSE,kTRUE,kTRUE,tenderPassData,kFALSE,PIDresponseSetting, kTRUE, kTRUE, customSplinesPass ) ; 
  AddTaskPIDResponse();
  
 //Tender Supplies
  AliPHOSTenderTask *tenderPHOS = AddAODPHOSTender("PHOSTenderTask","PHOStender","Run2NoNcellCut",1,kTRUE,"Run2TuneMCNoNcell",0.020) ;   //for MC
  AliPHOSTenderSupply* supply =  tenderPHOS->GetPHOSTenderSupply() ;   //for MC
  supply->SetPrivateOADBBadMap("$ALICE_PHYSICS/OADB/PHOS/PHOSBadMapsSoft.root") ;   //for MC

 //supply->ForceUsingCalibration("MCCalibration.root") ;   //for MC only (recalibration)

  //AliPHOSTenderTask *tenderPHOS = AddAODPHOSTender("PHOSTenderTask","PHOStender","Run2NoNcellCut",1,kFALSE,"Run2Tune") ;  // for data
  //AliPHOSTenderSupply* supply =  tenderPHOS->GetPHOSTenderSupply() ; //for data
  //supply->SetPrivateOADBBadMap("$ALICE_PHYSICS/OADB/PHOS/PHOSBadMapsSoft.root") ; //for data
  //supply->ForceUsingDummyRunNumber(265305); //!!!!!!!!!! 

 
  // Add my task
  AliAnalysisPHOSFluctuations *task1 = new AliAnalysisPHOSFluctuations("PhiFullEtaFullPt03_10");
  task1->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task1->SetChargedCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task1->SetChargedPtCut(0.3,1.) ;
  task1->SetPhotonCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi);
  task1->SetPhotonPtCut(0.3,1.) ;
  task1->SetPi0Cut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi) ;
  task1->SetPi0PtCut(0.3,1.) ;
  mgr->AddTask(task1);
  
  AliAnalysisPHOSFluctuations *task2 = new AliAnalysisPHOSFluctuations("PhiFullEtaGPHOSPt03_10");
  task2->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task2->SetChargedCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task2->SetChargedPtCut(0.3,1.) ;
  task2->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kFullPhi);
  task2->SetPhotonPtCut(0.3,1.) ;
  task2->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kFullPhi) ;
  task2->SetPi0PtCut(0.3,1.) ;
  mgr->AddTask(task2);

  AliAnalysisPHOSFluctuations *task2a = new AliAnalysisPHOSFluctuations("PhiPHOSEtaGPHOSPt03_10");
  task2a->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task2a->SetChargedCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kPhosPhi); 
  task2a->SetChargedPtCut(0.3,1.) ;
  task2a->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task2a->SetPhotonPtCut(0.3,1.) ;
  task2a->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task2a->SetPi0PtCut(0.3,1.) ;
  mgr->AddTask(task2a);

  AliAnalysisPHOSFluctuations *task3 = new AliAnalysisPHOSFluctuations("PhiPHOSEtaFullPt03_10");
  task3->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task3->SetChargedCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kPhosPhi); 
  task3->SetChargedPtCut(0.3,1.) ;
  task3->SetPhotonCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task3->SetPhotonPtCut(0.3,1.) ;
  task3->SetPi0Cut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task3->SetPi0PtCut(0.3,1.) ;
  mgr->AddTask(task3);


  AliAnalysisPHOSFluctuations *task3a = new AliAnalysisPHOSFluctuations("PhiGPHOSEtaFullPt03_10");
  task3a->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task3a->SetChargedCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task3a->SetChargedPtCut(0.3,1.) ;
  task3a->SetPhotonCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task3a->SetPhotonPtCut(0.3,1.) ;
  task3a->SetPi0Cut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task3a->SetPi0PtCut(0.3,1.) ;
  mgr->AddTask(task3a);

  AliAnalysisPHOSFluctuations *task4 = new AliAnalysisPHOSFluctuations("PhiGPHOSEtaGPHOSPt03_10");
  task4->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task4->SetChargedCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task4->SetChargedPtCut(0.3,1.) ;
  task4->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task4->SetPhotonPtCut(0.3,1.) ;
  task4->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task4->SetPi0PtCut(0.3,1.) ;
  mgr->AddTask(task4);

  //pi+- in PHOS rapidity, full Phi
  AliAnalysisPHOSFluctuations *task5 = new AliAnalysisPHOSFluctuations("PhiGPHOSEtaPHOSPt03_10");
  task5->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task5->SetChargedCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task5->SetChargedPtCut(0.3,1.) ;
  task5->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task5->SetPhotonPtCut(0.3,1.) ;
  task5->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task5->SetPi0PtCut(0.3,1.) ;
  mgr->AddTask(task5);

  //pi+- in PHOS rapidity, PHOS Phi
  AliAnalysisPHOSFluctuations *task6 = new AliAnalysisPHOSFluctuations("PhiPHOSEtaPHOSPt03_10");
  task6->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task6->SetChargedCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task6->SetChargedPtCut(0.3,1.) ;
  task6->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task6->SetPhotonPtCut(0.3,1.) ;
  task6->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task6->SetPi0PtCut(0.3,1.) ;
  mgr->AddTask(task6);

  //pi+- in PHOS rapidity, PHOS Phi
  AliAnalysisPHOSFluctuations *task7 = new AliAnalysisPHOSFluctuations("PhiPHOSEtaPHOSPt04_10");
  task7->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task7->SetChargedCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task7->SetChargedPtCut(0.4,1.) ;
  task7->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task7->SetPhotonPtCut(0.4,1.) ;
  task7->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task7->SetPi0PtCut(0.4,1.) ;
  mgr->AddTask(task7);

  //pi+- in PHOS rapidity, PHOS Phi
  AliAnalysisPHOSFluctuations *task8 = new AliAnalysisPHOSFluctuations("PhiPHOSEtaPHOSPtPi004_10");
  task8->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task8->SetChargedCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task8->SetChargedPtCut(0.4,1.) ;
  task8->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task8->SetPhotonPtCut(0.4,1.) ;
  task8->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task8->SetPi0PtCut(0.42,1.) ;
  mgr->AddTask(task8);


  //pi+- in PHOS rapidity, PHOS Phi
  AliAnalysisPHOSFluctuations *task9 = new AliAnalysisPHOSFluctuations("PhiPHOSEtaPHOSPt02_05");
  task9->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task9->SetChargedCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi); 
  task9->SetChargedPtCut(0.2,0.5) ;
  task9->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task9->SetPhotonPtCut(0.2,0.5) ;
  task9->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task9->SetPi0PtCut(0.2,0.5) ;
  mgr->AddTask(task9);

  //pi+- in PHOS rapidity, PHOS Phi
  AliAnalysisPHOSFluctuations *task10 = new AliAnalysisPHOSFluctuations("PhiFullEtaFullPt02_05");
  task10->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task10->SetChargedCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task10->SetChargedPtCut(0.2,0.5) ;
  task10->SetPhotonCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi);
  task10->SetPhotonPtCut(0.2,0.5) ;
  task10->SetPi0Cut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi) ;
  task10->SetPi0PtCut(0.2,0.5) ;
  mgr->AddTask(task10);

  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer(); 

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("PhiFullEtaFullPt03_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("PhiFullEtaGPHOSPt03_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput2a = mgr->CreateContainer("PhiPHOSEtaGPHOSPt03_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("PhiGPhosEtaFullPt03_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput3a = mgr->CreateContainer("PhiGPHOSEtaFullPt03_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("PhiGPhosEtaGPHOSPt03_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput5 = mgr->CreateContainer("PhiGPhosEtaPHOSPt03_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput6 = mgr->CreateContainer("PhiPhosEtaPHOSPt03_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput7 = mgr->CreateContainer("PhiPhosEtaPHOSPt04_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput8 = mgr->CreateContainer("PhiPhosEtaPHOSPtPi004_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput9 = mgr->CreateContainer("PhiPHOSEtaPHOSPt02_05",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput10 = mgr->CreateContainer("PhiFullEtaFullPt02_05",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  
  // Connect input/output
  mgr->ConnectInput(task1 , 0, cinput);
  mgr->ConnectOutput(task1, 1, coutput1);
  mgr->ConnectInput(task2 , 0, cinput);
  mgr->ConnectOutput(task2, 1, coutput2);
  mgr->ConnectInput(task2a , 0, cinput);
  mgr->ConnectOutput(task2a, 1, coutput2a);
  mgr->ConnectInput(task3 , 0, cinput);
  mgr->ConnectOutput(task3, 1, coutput3);
  mgr->ConnectInput(task3a , 0, cinput);
  mgr->ConnectOutput(task3a, 1, coutput3a);
  mgr->ConnectInput(task4 , 0, cinput);
  mgr->ConnectOutput(task4, 1, coutput4);
  mgr->ConnectInput(task5 , 0, cinput);
  mgr->ConnectOutput(task5, 1, coutput5);
  mgr->ConnectInput(task6 , 0, cinput);
  mgr->ConnectOutput(task6, 1, coutput6);
  mgr->ConnectInput(task7 , 0, cinput);
  mgr->ConnectOutput(task7, 1, coutput7);
  mgr->ConnectInput(task8 , 0, cinput);
  mgr->ConnectOutput(task8, 1, coutput8);
   
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("LOCAL", chain);
    // mgr->StartAnalysis("GRID", chain);
  }
  
}
