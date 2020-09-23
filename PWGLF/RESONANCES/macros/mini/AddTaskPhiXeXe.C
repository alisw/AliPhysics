/*********************************************
fbellini@cern.ch - created on 20 Nov 2017
Macro to add task for phi analysis in XeXe
*********************************************/

#if !defined (__CINT__) || defined (__CLING__)
#include "ConfigPhiXeXe.C"
#endif

AliRsnMiniAnalysisTask * AddTaskPhiXeXe(Int_t selectTaskConfig = 0, Bool_t isMC = kFALSE, Int_t nmix =5,TString multEstimator = "AliMultSelection_V0M");
AliRsnMiniAnalysisTask * AddTaskPhiXeXe( Bool_t isMC = kFALSE,
					 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s,
					 Float_t     nsigmaTPC = 2.0,
					 Float_t     nsigmaTOF = 2.0,
					 Int_t       aodFilterBit = 5,
           Int_t       customQualityCutsID = 100,
					 TString     multEstimator = "AliMultSelection_V0M",  
					 Int_t       nmix = 5,
					 Bool_t      enableMonitor = kTRUE,
					 TString     outNameSuffix = "default");

AliRsnMiniAnalysisTask * AddTaskPhiXeXe(Int_t selectTaskConfig, Bool_t isMC, Int_t nmix, TString multEstimator)
{
  //Select cuts configuration and returns the corresponding add task
  AliRsnMiniAnalysisTask * task = 0x0;
  switch (selectTaskConfig){
  case 1:
    task = (AliRsnMiniAnalysisTask *) AddTaskPhiXeXe(isMC, AliRsnCutSetDaughterParticle::kTPCTOFvetoPhiXeXe, 2.0, 4.0, 5, -1, multEstimator.Data(), 5, kTRUE, "tpc2sPtDep_tof4sveto5smism");    
    break;
  case 2:
    task = (AliRsnMiniAnalysisTask *) AddTaskPhiXeXe(isMC, AliRsnCutSetDaughterParticle::kTPCTOFvetoPhiXeXe, 3.0, 3.0, 5, -1, multEstimator.Data(), 5, kTRUE, "tpc3sPtDep_tof3sveto5smism");    
    break;
  case 11:
    task = (AliRsnMiniAnalysisTask *) AddTaskPhiXeXe(isMC, AliRsnCutSetDaughterParticle::kTPCTOFvetoPhiXeXe, 2.0, 3.0, 5, -1, multEstimator.Data(), 10, kTRUE, "default10");    
    break;
  default:
    task = (AliRsnMiniAnalysisTask *) AddTaskPhiXeXe(isMC, AliRsnCutSetDaughterParticle::kTPCTOFvetoPhiXeXe, 2.0, 3.0, 5, 101, multEstimator.Data(), 5, kTRUE, "default_LowBdca");    
  }

  return task;
}


AliRsnMiniAnalysisTask * AddTaskPhiXeXe(Bool_t isMC, AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate, Float_t nsigmaTPC, Float_t nsigmaTOF,
					                              Int_t aodFilterBit, Int_t customQualityCutsID, TString multEstimator, Int_t nmix, Bool_t enableMonitor, TString outNameSuffix)
{  
  //-------------------------------------------
  // event cuts
  //-------------------------------------------
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
      ::Error("AddTaskPhiXeXe", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   TString taskName = Form("PhiXeXe%s", (isMC ? "MC" : "Data"));
   
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   //task->SelectCollisionCandidates(triggerMask);//AOD
   task->UseESDTriggerMask(AliVEvent::kINT7);//ESD
   task->UseMultiplicity("AliMultSelection_V0M");
   task->SetUseBuiltinEventCuts();

   // set event mixing options
   task->UseContinuousMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   TString message = Form("\nevents to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n", nmix, maxDiffVzMix, maxDiffMultMix);
   Printf("AddTaskPhiXeXe :::: Event mixing configuration: %s", message.Data());
   
   mgr->AddTask(task);
   
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
   ConfigPhiXeXe(task, isMC, outNameSuffix.Data(), cutsPair, aodFilterBit, customQualityCutsID, cutKaCandidate, nsigmaTPC, nsigmaTOF, 15.0, enableMonitor, kTRUE, kFALSE);
#else
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigPhiXeXe.C");
   //gROOT->LoadMacro("./ConfigPhiXeXe.C");
   ConfigPhiXeXe(task, isMC, outNameSuffix.Data(), cutsPair, aodFilterBit, customQualityCutsID, cutKaCandidate, nsigmaTPC, nsigmaTOF, 15.0, enableMonitor, kTRUE, kFALSE);
#endif
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   Printf("AddTaskPhiXeXe - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
//   mgr->ConnectOutput(task, 2, lEvtQA);
   
   return task;
}
