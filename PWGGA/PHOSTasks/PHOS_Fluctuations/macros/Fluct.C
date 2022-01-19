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
  
  AliAnalysisPHOSFluctuations *task2 = new AliAnalysisPHOSFluctuations("PhiFullEtaGPhosPt03_10");
  task2->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task2->SetChargedCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task2->SetChargedPtCut(0.3,1.) ;
  task2->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kFullPhi);
  task2->SetPhotonPtCut(0.3,1.) ;
  task2->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kFullPhi) ;
  task2->SetPi0PtCut(0.3,1.) ;
  mgr->AddTask(task2);

  AliAnalysisPHOSFluctuations *task3 = new AliAnalysisPHOSFluctuations("PhiFullEtaPhosPt03_10");
  task3->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task3->SetChargedCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task3->SetChargedPtCut(0.3,1.) ;
  task3->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kFullPhi);
  task3->SetPhotonPtCut(0.3,1.) ;
  task3->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kFullPhi) ;
  task3->SetPi0PtCut(0.3,1.) ;
  mgr->AddTask(task3);

  AliAnalysisPHOSFluctuations *task4 = new AliAnalysisPHOSFluctuations("PhiPhosEtaFullPt03_10");
  task4->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task4->SetChargedCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kPhosPhi); 
  task4->SetChargedPtCut(0.3,1.) ;
  task4->SetPhotonCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task4->SetPhotonPtCut(0.3,1.) ;
  task4->SetPi0Cut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task4->SetPi0PtCut(0.3,1.) ;
  mgr->AddTask(task4);

  AliAnalysisPHOSFluctuations *task5 = new AliAnalysisPHOSFluctuations("PhiPhosEtaGPhosPt03_10");
  task5->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task5->SetChargedCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kPhosPhi); 
  task5->SetChargedPtCut(0.3,1.) ;
  task5->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task5->SetPhotonPtCut(0.3,1.) ;
  task5->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task5->SetPi0PtCut(0.3,1.) ;
  mgr->AddTask(task5);

  AliAnalysisPHOSFluctuations *task6 = new AliAnalysisPHOSFluctuations("PhiPhosEtaPhosPt03_10");
  task6->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task6->SetChargedCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi); 
  task6->SetChargedPtCut(0.3,1.) ;
  task6->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task6->SetPhotonPtCut(0.3,1.) ;
  task6->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task6->SetPi0PtCut(0.3,1.) ;
  mgr->AddTask(task6);

  AliAnalysisPHOSFluctuations *task7 = new AliAnalysisPHOSFluctuations("PhiGPhosEtaFullPt03_10");
  task7->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task7->SetChargedCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task7->SetChargedPtCut(0.3,1.) ;
  task7->SetPhotonCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task7->SetPhotonPtCut(0.3,1.) ;
  task7->SetPi0Cut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task7->SetPi0PtCut(0.3,1.) ;
  mgr->AddTask(task7);

  AliAnalysisPHOSFluctuations *task8 = new AliAnalysisPHOSFluctuations("PhiGPhosEtaGPhosPt03_10");
  task8->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task8->SetChargedCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task8->SetChargedPtCut(0.3,1.) ;
  task8->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task8->SetPhotonPtCut(0.3,1.) ;
  task8->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task8->SetPi0PtCut(0.3,1.) ;
  mgr->AddTask(task8);

  //pi+- in PHOS rapidity, full Phi
  AliAnalysisPHOSFluctuations *task9 = new AliAnalysisPHOSFluctuations("PhiGPhosEtaPhosPt03_10");
  task9->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task9->SetChargedCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task9->SetChargedPtCut(0.3,1.) ;
  task9->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task9->SetPhotonPtCut(0.3,1.) ;
  task9->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task9->SetPi0PtCut(0.3,1.) ;
  mgr->AddTask(task9);



  AliAnalysisPHOSFluctuations *task10 = new AliAnalysisPHOSFluctuations("PhiFullEtaFullPt04_10");
  task10->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task10->SetChargedCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task10->SetChargedPtCut(0.4,1.) ;
  task10->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task10->SetPhotonPtCut(0.4,1.) ;
  task10->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task10->SetPi0PtCut(0.4,1.) ;
  mgr->AddTask(task10);

  //pi+- in PHOS rapidity, PHOS Phi
  AliAnalysisPHOSFluctuations *task11 = new AliAnalysisPHOSFluctuations("PhiPhosEtaPhosPt04_10");
  task11->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task11->SetChargedCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task11->SetChargedPtCut(0.4,1.) ;
  task11->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task11->SetPhotonPtCut(0.4,1.) ;
  task11->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task11->SetPi0PtCut(0.4,1.) ;
  mgr->AddTask(task11);


  //pi+- in PHOS rapidity, PHOS Phi
  AliAnalysisPHOSFluctuations *task12 = new AliAnalysisPHOSFluctuations("PhiGPhosEtaPhosPt04_10");
  task12->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task12->SetChargedCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task12->SetChargedPtCut(0.4,1.) ;
  task12->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task12->SetPhotonPtCut(0.4,1.) ;
  task12->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task12->SetPi0PtCut(0.4,1.) ;
  mgr->AddTask(task12);


  AliAnalysisPHOSFluctuations *task13 = new AliAnalysisPHOSFluctuations("PhiFullEtaFullPt01_03");
  task13->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task13->SetChargedCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task13->SetChargedPtCut(0.1,0.3) ;
  task13->SetPhotonCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi);
  task13->SetPhotonPtCut(0.1,0.3) ;
  task13->SetPi0Cut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi) ;
  task13->SetPi0PtCut(0.1,0.3) ;
  mgr->AddTask(task13);

  AliAnalysisPHOSFluctuations *task14 = new AliAnalysisPHOSFluctuations("PhiPhosEtaPhosPt01_03");
  task14->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task14->SetChargedCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi); 
  task14->SetChargedPtCut(0.1,0.3) ;
  task14->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task14->SetPhotonPtCut(0.1,0.3) ;
  task14->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task14->SetPi0PtCut(0.1,0.3) ;
  mgr->AddTask(task14);


  AliAnalysisPHOSFluctuations *task15 = new AliAnalysisPHOSFluctuations("PhiFullEtaFullPt02_05");
  task15->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task15->SetChargedCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task15->SetChargedPtCut(0.2,0.5) ;
  task15->SetPhotonCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi);
  task15->SetPhotonPtCut(0.2,0.5) ;
  task15->SetPi0Cut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi) ;
  task15->SetPi0PtCut(0.2,0.5) ;
  mgr->AddTask(task15);

  AliAnalysisPHOSFluctuations *task16 = new AliAnalysisPHOSFluctuations("PhiPhosEtaPhosPt02_05");
  task16->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task16->SetChargedCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi); 
  task16->SetChargedPtCut(0.2,0.5) ;
  task16->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task16->SetPhotonPtCut(0.2,0.5) ;
  task16->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task16->SetPi0PtCut(0.2,0.5) ;
  mgr->AddTask(task16);


 AliAnalysisPHOSFluctuations *task17 = new AliAnalysisPHOSFluctuations("PhiFullEtaFullPt05_10");
  task17->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task17->SetChargedCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi); 
  task17->SetChargedPtCut(0.5,1.) ;
  task17->SetPhotonCut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi);
  task17->SetPhotonPtCut(0.5,1.) ;
  task17->SetPi0Cut(AliAnalysisPHOSFluctuations::kTpcEta,AliAnalysisPHOSFluctuations::kFullPhi) ;
  task17->SetPi0PtCut(0.5,1.) ;
  mgr->AddTask(task17);

  AliAnalysisPHOSFluctuations *task18 = new AliAnalysisPHOSFluctuations("PhiPhosEtaPhosPt05_10");
  task18->SetRunType(AliAnalysisPHOSFluctuations::kpp); 
  task18->SetChargedCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi); 
  task18->SetChargedPtCut(0.5,1.) ;
  task18->SetPhotonCut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi);
  task18->SetPhotonPtCut(0.5,1.) ;
  task18->SetPi0Cut(AliAnalysisPHOSFluctuations::kPhosEta,AliAnalysisPHOSFluctuations::kPhosPhi) ;
  task18->SetPi0PtCut(0.5,1.) ;
  mgr->AddTask(task18);


  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer(); 

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("PhiFullEtaFullPt03_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("PhiFullEtaGPhosPt03_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("PhiFullEtaPhosPt03_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("PhiPhosEtaFullPt03_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput5 = mgr->CreateContainer("PhiPhosEtaGPhosPt03_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput6 = mgr->CreateContainer("PhiPhosEtaPhosPt03_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput7 = mgr->CreateContainer("PhiGPhosEtaFullPt03_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput8 = mgr->CreateContainer("PhiGPhosEtaGPhosPt03_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput9 = mgr->CreateContainer("PhiGPhosEtaPhosPt03_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput10 = mgr->CreateContainer("PhiFullEtaFullPt04_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput11 = mgr->CreateContainer("PhiPhosEtaPhosPt04_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput12 = mgr->CreateContainer("PhiGPhosEtaPhosPt04_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput13 = mgr->CreateContainer("PhiFullEtaFullPt01_03",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput14 = mgr->CreateContainer("PhiPhosEtaPhosPt01_03",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput15 = mgr->CreateContainer("PhiFullEtaFullPt02_05",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput16 = mgr->CreateContainer("PhiPhosEtaPhosPt02_05",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput17 = mgr->CreateContainer("PhiFullEtaFullPt05_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  AliAnalysisDataContainer *coutput18 = mgr->CreateContainer("PhiPhosEtaPhosPt05_10",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  
  // Connect input/output
  mgr->ConnectInput(task1 , 0, cinput);
  mgr->ConnectOutput(task1, 1, coutput1);
  mgr->ConnectInput(task2 , 0, cinput);
  mgr->ConnectOutput(task2, 1, coutput2);
  mgr->ConnectInput(task3 , 0, cinput);
  mgr->ConnectOutput(task3, 1, coutput3);
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
  mgr->ConnectInput(task9 , 0, cinput);
  mgr->ConnectOutput(task9, 1, coutput9);
  mgr->ConnectInput(task10 , 0, cinput);
  mgr->ConnectOutput(task10, 1, coutput10);
  mgr->ConnectInput(task11 , 0, cinput);
  mgr->ConnectOutput(task11, 1, coutput11);
  mgr->ConnectInput(task12 , 0, cinput);
  mgr->ConnectOutput(task12, 1, coutput12);
  mgr->ConnectInput(task13 , 0, cinput);
  mgr->ConnectOutput(task13, 1, coutput13);
  mgr->ConnectInput(task14 , 0, cinput);
  mgr->ConnectOutput(task14, 1, coutput14);
  mgr->ConnectInput(task15 , 0, cinput);
  mgr->ConnectOutput(task15, 1, coutput15);
  mgr->ConnectInput(task16 , 0, cinput);
  mgr->ConnectOutput(task16, 1, coutput16);
  mgr->ConnectInput(task17 , 0, cinput);
  mgr->ConnectOutput(task17, 1, coutput17);
  mgr->ConnectInput(task18 , 0, cinput);
  mgr->ConnectOutput(task18, 1, coutput18);
   
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("LOCAL", chain);
    // mgr->StartAnalysis("GRID", chain);
  }
  
}
