/*********************************************
fbellini@cern.ch - created on 20 Nov 2017
Macro to add task for phi analysis in XeXe
*********************************************/

#if !defined (__CINT__) || defined (__CLING__)
#include "ConfigPhiVn.C"
#endif

AliRsnMiniTaskPhiVn * AddTaskPhiVn(Int_t selectTaskConfig = 0, Bool_t isMC = kFALSE, TString multEstimator = "AliMultSelection_V0M");
AliRsnMiniTaskPhiVn * AddTaskPhiVn( Bool_t      isMC = kFALSE,
					 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s,
					 Float_t     nsigmaTPC = 2.0,
					 Float_t     nsigmaTOF = 2.0,
					 Int_t       aodFilterBit = 5,
					 TString     multEstimator = "AliMultSelection_V0M",  
					 Int_t       nmix = 5,
					 Bool_t      enableMonitor = kTRUE,
					 TString     outNameSuffix = "tpc2s_tof3sveto");

AliRsnMiniTaskPhiVn * AddTaskPhiVn(Int_t selectTaskConfig, Bool_t isMC, TString multEstimator)
{
  //Select cuts configuration and returns the corresponding add task
  AliRsnMiniTaskPhiVn * task = 0x0;
  
  switch (selectTaskConfig){
  case 1 :
    task = (AliRsnMiniTaskPhiVn *) AddTaskPhiVn(isMC, AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s, 3.0, 3.0, 5, multEstimator.Data(), 5, kTRUE, "tpc3s_tof3sveto");
    break;
  case 2 :
    task = (AliRsnMiniTaskPhiVn *) AddTaskPhiVn(isMC, AliRsnCutSetDaughterParticle::kTPCpidTOFveto4s, 2.0, 4.0, 5, multEstimator.Data(), 5, kTRUE, "tpc2s_tof4sveto");
    break;
  case 3 :
    task = (AliRsnMiniTaskPhiVn *) AddTaskPhiVn(isMC, AliRsnCutSetDaughterParticle::kFastTPCpidNsigma, 3.0, 10.0,  5, multEstimator.Data(), 5, kTRUE, "tpc3s");
    break;
  case 4 :
    task = (AliRsnMiniTaskPhiVn *) AddTaskPhiVn(isMC, AliRsnCutSetDaughterParticle::kFastTOFpidNsigma, 10.0, 3.0, 5, multEstimator.Data(), 5, kTRUE, "tof3s");
    break;
  case 5 :
    task = (AliRsnMiniTaskPhiVn *) AddTaskPhiVn(isMC, AliRsnCutSetDaughterParticle::kTPCpidphipp2015, 2.0, 10.0, 5, multEstimator.Data(), 5, kTRUE, "tpc2sPtDep");    
    break;
  case 6 :
    task = (AliRsnMiniTaskPhiVn *) AddTaskPhiVn(isMC, AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015, 3.0, 3.0, 5, multEstimator.Data(), 5, kTRUE, "tpc3sPtDep_tof3sveto");    
    break;
  case 7 :
    task = (AliRsnMiniTaskPhiVn *) AddTaskPhiVn(isMC, AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015, 2.0, 3.0, 5, multEstimator.Data(), 5, kTRUE, "tpc2sPtDep_tof3sveto");    
    break;
  case 8 :
    task = (AliRsnMiniTaskPhiVn *) AddTaskPhiVn(isMC, AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015, 2.0, 2.0, 5, multEstimator.Data(), 5, kTRUE, "tpc2sPtDep_tof2sveto");    
    break;
  default:
    task = (AliRsnMiniTaskPhiVn *) AddTaskPhiVn(isMC, AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s, 2.0, 3.0, 5, multEstimator.Data(), 5, kTRUE, "tpc2s_tof3sveto");
  }

  return task;
}


AliRsnMiniTaskPhiVn * AddTaskPhiVn(Bool_t isMC, AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate, Float_t nsigmaTPC, Float_t nsigmaTOF,
					Int_t aodFilterBit, TString multEstimator, Int_t nmix, Bool_t enableMonitor, TString outNameSuffix)
{  
  //-------------------------------------------
  // event cuts
  //-------------------------------------------
  UInt_t   triggerMask  = AliVEvent::kINT7;
  Bool_t   rejectPileUp = kTRUE;
  Double_t vtxZcut = 10.0;

  //-------------------------------------------
  // event mixing settings
  //-------------------------------------------
  Float_t     maxDiffVzMix = 1.0;
  Float_t     maxDiffMultMix = 5.0;
  
  //-------------------------------------------
  // pair cuts
  //-------------------------------------------
  Double_t    minYlab = -0.5;
  Double_t    maxYlab = 0.5;

  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskPhiVn", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   TString taskName = Form("PhiXeXe%s", (isMC ? "MC" : "Data"));
   
   AliRsnMiniTaskPhiVn *task = new AliRsnMiniTaskPhiVn(taskName.Data(), isMC);
   task->SetNharmToProcess(2);//new
   
   task->SetSelCharge('+', '-');
   task->SetSelPid(AliRsnDaughter::kKaon,AliRsnDaughter::kKaon);

   //task->SelectCollisionCandidates(triggerMask);//AOD
   task->UseESDTriggerMask(triggerMask);//ESD
   task->UseMultiplicity("AliMultSelection_V0M");
   
   // set event mixing options
   task->UseContinuousMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   TString message = Form("\nevents to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n", nmix, maxDiffVzMix, maxDiffMultMix);
   Printf("AddTaskPhiVn :::: Event mixing configuration: %s", message.Data());
   mgr->AddTask(task);

   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //   
   AliRsnCutEventUtils* cutEventUtils = new AliRsnCutEventUtils("cutEventUtils", kTRUE, rejectPileUp);
   cutEventUtils->SetRemovePileUppA2013(kFALSE);
   cutEventUtils->SetCheckAcceptedMultSelection();
   Printf("AddTaskPhiVn :::: Centrality estimator: %s", multEstimator.Data());
   
   AliRsnCutSet* eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutEventUtils);
   eventCuts->SetCutScheme(Form("%s", cutEventUtils->GetName()));

   task->SetEventCuts(eventCuts); 

   //
   // -- EVENT-ONLY COMPUTATIONS -------------------------------------------------------------------
   //   
   //vertex
   Int_t vtxID = task->CreateValue(AliRsnMiniValue::kVz, kFALSE);
   //multiplicity or centrality
   Int_t multID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);

   AliRsnMiniOutput *outVtx = task->CreateOutput("eventVtx", "HIST", "EVENT");
   outVtx->AddAxis(vtxID, 240, -12.0, 12.0);
   outVtx->AddAxis(multID, 20, 0.0, 100.0);
   
   AliRsnMiniOutput *outMult = task->CreateOutput("eventMult", "HIST", "EVENT");
   outMult->AddAxis(multID, 100, 0.0, 100.0);
   
   // ------------------------------------------------------
   // PAIR CUTS (common to all resonances) 
   // ------------------------------------------------------

   AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
   cutY->SetRangeD(minYlab, maxYlab);
   
   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());
   
   // ------------------------------------------------------
   // CONFIG ANALYSIS
   // ------------------------------------------------------
#if !defined (__CINT__) || defined (__CLING__)
   ConfigPhiVn(task, isMC, outNameSuffix.Data(), cutsPair, aodFilterBit, cutKaCandidate, nsigmaTPC, nsigmaTOF, 15.0, enableMonitor,kTRUE, kFALSE);
#else
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigPhiVn.C");
   //gROOT->LoadMacro("ConfigPhiVn.C");
   ConfigPhiVn(task, isMC, outNameSuffix.Data(), cutsPair, aodFilterBit, cutKaCandidate, nsigmaTPC, nsigmaTOF, 15.0, enableMonitor, kTRUE, kFALSE);
#endif
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   Printf("AddTaskPhiVn - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}

