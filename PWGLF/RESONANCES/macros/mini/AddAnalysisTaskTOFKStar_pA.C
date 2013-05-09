/****************************************************************************************************
	
     Variation of AddAnalysisTaskTOFKStar.C by F. Bellini (fbellini@cern.ch) for pA analysis

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
****************************************************************************************************/

AliRsnMiniAnalysisTask * AddAnalysisTaskTOFKStar_pA
(
   Bool_t      isMC,
   Bool_t      isUseCent,
   Int_t       aodFilterBit = 5,
   AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPiCandidate = AliRsnCutSetDaughterParticle::kTOFpidKstarPbPb2010,
   AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTOFpidKstarPbPb2010,
   Float_t     nsigmaPi = 2.0,
   Float_t     nsigmaKa = 2.0,
   Bool_t      enableMonitor = kTRUE,
   Bool_t      IsMcTrueOnly = kFALSE,
   Int_t       nmix = 0,
   Float_t     maxDiffVzMix = 1.0,
   Float_t     maxDiffMultMix = 10.0,
   Float_t     maxDiffAngleMixDeg = 20.0,
   Bool_t      isAOD49 = kFALSE,
   TString     outNameSuffix = "",
   Bool_t      useMixLS = 0,
   Int_t       signedPdg = 313,
   TString     monitorOpt = ""
)
{  
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddAnalysisTaskTOFKStar", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   TString taskName = Form("TOFKStar%s%s_%i%i", "pA", (isMC ? "MC" : "Data"), (Int_t)cutPiCandidate,(Int_t)cutKaCandidate );
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
    
   if (!isUseCent) 
     task->UseMultiplicity("QUALITY");
   else
   task->UseCentrality("V0A");   
   // set event mixing options
   task->UseContinuousMix();
   //task->UseBinnedMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   ::Info("AddAnalysisTaskTOFKStar", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f", nmix, maxDiffVzMix, maxDiffMultMix));
   
   mgr->AddTask(task);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //  
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
   if (ispA) cutVertex->SetCheckPileUp(kTRUE);   // set the check for pileup
   
   // define and fill cut set for event cut
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutVertex);
   eventCuts->SetCutScheme(cutVertex->GetName());
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
   if (!isUseCent) 
     outMult->AddAxis(multID, 400, 0.0, 400.0);
   else
     outMult->AddAxis(multID, 100, 0.0, 100.0);
   
   }

   TH2F* hvz=new TH2F("hVzVsCent","",100,0.,100., 220,-11.,11.);
   task->SetEventQAHist("vz",hvz);//plugs this histogram into the fHAEventVz data member

   TH2F* hmc=new TH2F("MultiVsCent","",100,0.,100.,4000,0.,4000.);
   hmc->GetYaxis()->SetTitle("QUALITY");
   task->SetEventQAHist("multicent",hmc);//plugs this histogram into the fHAEventMultiCent data member

   TH2F* hep=new TH2F("hEventPlaneVsCent","",100,0.,100., 180,0.,TMath::Pi());
   task->SetEventQAHist("eventplane",hep);//plugs this histogram into the fHAEventPlane data member

   //
   // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
   //
   AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
   cutY->SetRangeD(-0.032, 0.968); //   cutY->SetRangeD(-0.5, 0.5);
   
   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());
   
   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------
   gROOT->LoadMacro("$ALICE_ROOT/PWGLF/RESONANCES/macros/mini/ConfigTOFanalysisKStar.C");
   if (isMC) {
   
     if (((Int_t)cutPiCandidate<4) && ((Int_t)cutKaCandidate<4))
       Printf("========================== MC analysis - no PID used for efficiency estimation");
     else 
       Printf("========================== MC analysis - PID cuts used");
   } else 
     Printf("========================== DATA analysis - PID cuts used");
   if (!ConfigTOFanalysisKStar(task, isMC, !isUseCent, "", cutsPair, aodFilterBit, cutPiCandidate, cutKaCandidate, nsigmaPi, nsigmaKa, enableMonitor, isMC&IsMcTrueOnly, useMixLS, signedPdg, monitorOpt.Data())) return 0x0;
   
   //n.
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddAnalysisTaskTOFKStar - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}
