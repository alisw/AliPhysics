#if defined(__CLING__)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
R__ADD_INCLUDE_PATH(./)
#include <TAlienCollection.h>
#include <TAlienResult.h>
#include "AliAODInputHandler.h"
#include <OADB/macros/AddTaskPhysicsSelection.C>
#include <OADB/macros/AddTaskCentrality.C>
#include <PWGGA/PHOSTasks/PHOS_PbPb/AddAODPHOSTender.C>
#endif
void Tagging()
{
  gErrorIgnoreLevel=2001 ;  
  
  TChain* chain = new TChain("aodTree");

  chain->Add("AliAOD.root");

  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("Pi0Spectrum");
  
  // ESD input handler
  AliAODInputHandler* esdH = new AliAODInputHandler();
  mgr->SetInputEventHandler(esdH);
  
 
  //Set run number for which bad map and calibrations will be used
  Int_t irun = 265594 ; 
  
  //Tender Supplies
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/PHOSTasks/PHOS_PbPb/AddAODPHOSTender.C");
  AliPHOSTenderTask *tenderPHOS = AddAODPHOSTender("PHOSTenderTask","PHOStender","Run2NoNcellCut",1,kTRUE,"Run2TuneMCNoNcellHighPtFix",0.020) ;  

  AliPHOSTenderSupply* supply =  tenderPHOS->GetPHOSTenderSupply() ;
  supply->SetPrivateOADBBadMap("$ALICE_PHYSICS/OADB/PHOS/PHOSBadMapsSoft.root") ;
  supply->ForceUsingDummyRunNumber(irun);
  
  Float_t timeCut = 100.e-9 ;
  //Centrality estimator: 1: V0A/C, 2: V0M, 3: ZNA/C,  4: CL1
  Int_t centralityEstinator=13 ;
  Int_t binLimits[1]={100};
  TArrayI multBins(1,binLimits) ;
  
  //Primary particles weight to account the difference
  //between PHOSGenLib spectrum and true one
  TArrayD ar(7) ;
  //PHOSGenLib pi0 spectrum 17.03.2020
  ar[0]= 8.14456e-08 ;
  ar[1]= 7.55352e-01 ;
  ar[2]= 6.23998e+01 ;
  ar[3]=-3.51140e+00 ;
  ar[4]= 1.55292e+01 ;
  ar[5]= 5.09556e+00 ;
  ar[6]= 0. ;  

   
  // Add my task
  AliAnalysisTaskTaggedPhotons *task1 = new AliAnalysisTaskTaggedPhotons("Pi0SpectrumkINT7");

  task1->SetTimeCut(30.e-9) ;
  task1->SetMultiplicityBins(&multBins) ;
  task1->SetTrigger(kTRUE) ;
  task1->SetCluCutType(AliAnalysisTaskTaggedPhotons::kLowECut) ;
  task1->SetMC(kTRUE);
  task1->SetFastMC();
  task1->SetCentralityEstimator(centralityEstinator) ; 
  task1->SetRunNumber(irun); 
  task1->SetMCType(AliAnalysisTaskTaggedPhotons::kSinglePi0) ;
  task1->SetPi0WeightParameters(&ar) ;
  task1->SetNonLinearity(1.0245,0.013,0.5) ;
  mgr->AddTask(task1);

  AliAnalysisTaskTaggedPhotons *task2 = new AliAnalysisTaskTaggedPhotons("Pi0SpectrumkPHI7L0");
  task2->SetTimeCut(30.e-9) ;
  task2->SetMultiplicityBins(&multBins) ;
  task2->SetTrigger(kFALSE) ;
  task2->SetPHOSTrigger(AliAnalysisTaskTaggedPhotons::kPHOSL0) ;
  task2->SetCluCutType(AliAnalysisTaskTaggedPhotons::kLowECut) ;
  task2->SetMC(kTRUE);
  task2->SetFastMC();
  task2->SetCentralityEstimator(centralityEstinator) ;
  task2->SetRunNumber(irun); 
  task2->SetMCType(AliAnalysisTaskTaggedPhotons::kSinglePi0) ;
  task2->SetPi0WeightParameters(&ar) ;
  task2->SetNonLinearity(1.0245,0.013,0.5) ;
  mgr->AddTask(task2);  
  
  AliAnalysisTaskTaggedPhotons *task3 = new AliAnalysisTaskTaggedPhotons("Pi0SpectrumkPHI77L1");
  task3->SetTimeCut(30.e-9) ;
  task3->SetMultiplicityBins(&multBins) ;
  task3->SetTrigger(kFALSE) ;
  task3->SetPHOSTrigger(AliAnalysisTaskTaggedPhotons::kPHOSL1low) ;
  task3->SetCluCutType(AliAnalysisTaskTaggedPhotons::kLowECut) ;
  task3->SetMC(kTRUE);
  task3->SetFastMC();
  task3->SetCentralityEstimator(centralityEstinator) ;
  task3->SetRunNumber(irun); 
  task3->SetMCType(AliAnalysisTaskTaggedPhotons::kSinglePi0) ;
  task3->SetPi0WeightParameters(&ar) ;
  task3->SetNonLinearity(1.0245,0.013,0.5) ;
  mgr->AddTask(task3);  
  
  
  
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer(); 

//   Create containers for input/output
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("TaggingkINT7",TList::Class(),AliAnalysisManager::kOutputContainer,"histos16r.root");
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("TaggingkPHI7L0",TList::Class(),AliAnalysisManager::kOutputContainer,"histos16r.root");
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("TaggingkPHI7L1",TList::Class(),AliAnalysisManager::kOutputContainer,"histos16r.root");

  // Connect input/output
  mgr->ConnectInput(task1 , 0, cinput);
  mgr->ConnectOutput(task1, 1, coutput1);
  mgr->ConnectInput(task2 , 0, cinput);
  mgr->ConnectOutput(task2, 1, coutput2);
  mgr->ConnectInput(task3 , 0, cinput);
  mgr->ConnectOutput(task3, 1, coutput3);

  
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
  }
  
}
