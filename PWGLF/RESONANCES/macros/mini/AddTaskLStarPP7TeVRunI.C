/***************************************************************************
            adash@cern.ch - last modified on 03/04/2013
//  modified by rsingh@cern.ch -- 19-11-2016
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

AliRsnMiniAnalysisTask * AddTaskLStarPP7TeVRunI
	(
	 Bool_t      isMC           = kFALSE,
	 Bool_t      isPP           = kTRUE,
	 Int_t       trCut         = 2010,
	 TString     outNameSuffix  = "tofveto3stpc3s",
	 // Int_t       customQualityCutsID = 1,
	 Int_t       customQualityCutsID = AliRsnCutSetDaughterParticle::kDisableCustom,
	 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate=AliRsnCutSetDaughterParticle::kTPCTOFpidLstar_test1,
	 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate=AliRsnCutSetDaughterParticle::kTPCTOFpidLstar_test1,
	 Float_t     nsigmaP       = 3.0,
	 Float_t     nsigmaK        = 3.0,
	 Float_t     nsigmatofP = 3.0,
	 Float_t     nsigmatofK  = 3.0,
	 Bool_t      enableMonitor  = kTRUE,
	 Int_t       nmix           = 10,
	 Float_t     maxDiffVzMix   = 1.0,
	 Float_t     maxDiffMultMix = 5.0
	 )
   {  

  UInt_t      triggerMask  = AliVEvent::kMB;
  Bool_t      rejectPileUp = kTRUE;
  Double_t    vtxZcut      = 10.0;//cm, default cut on vtx z
  if(!isPP || isMC) rejectPileUp = kFALSE;
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskLStarPP8TeV", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   
   TString taskName = Form("LStar%s%s", (isPP? "pp" : "PbPb"), (isMC ? "MC" : "Data"));
   
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
       task->UseESDTriggerMask(triggerMask);    //ESD
   //   task->SelectCollisionCandidates(triggerMask); //AOD
   
  
   if (isPP) 
         task->UseMultiplicity("MULTV0M");
        else
	  // task->UseMultiplicity("QUALITY");//Only for RunII
      task->UseCentrality("V0M");   
   // set event mixing options
   task->UseContinuousMix();
   //task->UseBinnedMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   task->UseMC(isMC);
   ::Info("AddTaskLStarPP8TeV", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n", nmix, maxDiffVzMix, maxDiffMultMix));
   
   mgr->AddTask(task);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //  
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
   if (isPP && (!isMC)) cutVertex->SetCheckPileUp(rejectPileUp);   // set the check for pileup

   // define and fill cut set for event cut
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutVertex);
   eventCuts->SetCutScheme(Form("%s",cutVertex->GetName()));
   // set cuts in task
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
   cutY->SetRangeD(-0.5, 0.5);
  // cutY->SetRangeD(-0.9, 0.9);
   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());
   
   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------
   //
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigLStarPP7TeVRunI.C");
   //   gROOT->LoadMacro("ConfigLStarPP7TeVRunI.C");
   if (!ConfigLStarPP7TeVRunI(task, isMC, isPP, cutsPair, trCut, customQualityCutsID,cutKaCandidate,cutPrCandidate, nsigmaP,nsigmaK,nsigmatofP,nsigmatofK, enableMonitor)) return 0x0;

   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddTaskLStarPP8TeV - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}
