/***************************************************************************
            adash@cern.ch - last modified on 05/09/2016
//
// General macro to configure the RSN analysis task.
// It calls all configs desired by the user, by means
// of the boolean switches defined in the first lines.
// ---
// Inputs:
//  1) flag to know if running on MC or data
//  2) path where all configs are stored
// ---
// Returns:
//  kTRUE  --> initialization successful
//  kFALSE --> initialization failed (some config gave errors)
//
****************************************************************************/

AliRsnMiniAnalysisTask * AddTaskLStarpPbRunII(
						Bool_t      isMC                = kFALSE,
						Bool_t      isPP                = kFALSE,
						Int_t       Strcut              = 2011,
						Int_t       customQualityCutsID = AliRsnCutSetDaughterParticle::kDisableCustom,
					       // AliRsnCutSetDaughterParticle::ERsnDaughterCutSet /cutKaCandidate=AliRsnCutSetDaughterParticle::kTPCTOFpidphikstarpPb2016,
                                                AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate=AliRsnCutSetDaughterParticle::kTPCTOFpidLstar,
                                                AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate   =AliRsnCutSetDaughterParticle::kTPCTOFpidLstar,


						Float_t     nsigmaP            = 2.0,
						Float_t     nsigmaK            = 2.0,
					       	Float_t     nsigmatofP         = 3.0,
			                        Float_t     nsigmatofK         = 3.0,
						Bool_t      enableMonitor       = kTRUE,
						Int_t       nmix                = 10,
						Float_t     maxDiffVzMix        = 1.0,
						Float_t     maxDiffMultMix      = 5.0,
						TString     outNameSuffix       = "pPb"
						)
{  
  Bool_t      rejectPileUp = kTRUE;
  //if(!isPP || isMC) rejectPileUp = kFALSE;

  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskLStarpPbRunTwo", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   TString taskName = Form("LStar%s%s", (isPP? "pp" : "PPb"), (isMC ? "MC" : "Data"));
   
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   task->UseESDTriggerMask(AliVEvent::kINT7);
   //if (isPP) 
   //task->UseMultiplicity("QUALITY");
   //else
   //  task->UseMultiplicity("AliMultSelection_V0M");//Only for RunII
    task->UseMultiplicity("AliMultSelection_V0A");//Only for RunII

   //task->UseMultiplicity("AliMultSelection_ZNA");//Only for RunII

   // task->UseMultiplicity("AliMultSelection_ZNC");//Only for RunII
   // set event mixing options
   task->UseContinuousMix();
   //task->UseBinnedMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   task->UseMC(isMC);
   ::Info("AddTaskLStarpPbRunTwo", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n", nmix, maxDiffVzMix, maxDiffMultMix));
   
   mgr->AddTask(task);
      
   AliRsnCutEventUtils* cutEventUtils=new AliRsnCutEventUtils("cutEventUtils",kTRUE,rejectPileUp);
   cutEventUtils->SetCheckAcceptedMultSelection();
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutEventUtils);
   eventCuts->SetCutScheme(Form("%s",cutEventUtils->GetName()));
   task->SetEventCuts(eventCuts);
   
   //
   // -- EVENT-ONLY COMPUTATIONS -------------------------------------------------------------------
   //   
   //vertex
   Int_t vtxID = task->CreateValue(AliRsnMiniValue::kVz, kFALSE);
   AliRsnMiniOutput *outVtx = task->CreateOutput("eventVtx", "HIST", "EVENT");
   outVtx->AddAxis(vtxID, 400, -20.0, 20.0);
   
   //multiplicity or centrality
   Int_t multID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
   AliRsnMiniOutput *outMult = task->CreateOutput("eventMult", "HIST", "EVENT");
   if (isPP) 
     outMult->AddAxis(multID, 400, 0.0, 400.0);
   else
     outMult->AddAxis(multID, 100, 0.0, 100.0);
   
   //
   // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
   //
   AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
   //      cutY->SetRangeD(-0.9, 0.9);



   //   cutY->SetRangeD(0.5, 0.8);
   //cutY->SetRangeD(0.035, 0.535);
   //cutY->SetRangeD(-0.465, 0.035);// 0 < y_cm < 0.5; y_cm = y_lab + 0.465
   //  cutY->SetRangeD(-0.8, -0.5);// 0 < y_cm < 0.5; y_cm = y_lab + 0.465
    cutY->SetRangeD(-0.465, 0.035);// 0 < y_cm < 0.5; y_cm = y_lab + 0.465
   //  cutY->SetRangeD(-0.765, -0.465);// 0 < y_cm < 0.5; y_cm = y_lab + 0.465

   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());
   
   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------
   //
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigLStarpPbRunII.C");
  // gROOT->LoadMacro("ConfigLStarpPbRunII.C");
   if (!ConfigLStarpPbRunII(task, isMC, isPP, cutsPair,Strcut, customQualityCutsID,cutPrCandidate,cutKaCandidate,nsigmaP,nsigmaK,nsigmatofP,nsigmatofK,enableMonitor)) return 0x0;

   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //outputFileName += ":Rsn";
   Printf("AddTaskLStarpPbRunTwo - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}
