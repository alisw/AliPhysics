/***************************************************************************
              fbellini@cern.ch - created and last modified on 20/09/2012
//
//Lauches KStar analysis with rsn mini package
//Allows basic configuration of pile-up check
//
****************************************************************************/

AliRsnMiniAnalysisTask * AddTaskKStar_pAPileUp
(
   Bool_t      isMC,
   Bool_t      isPP,
   TString     outNameSuffix = "tof2s_pileup",
   Bool_t      rmFirstEvtChunk = kFALSE,
   Bool_t      rejectPileUp = kTRUE,
   Int_t       MinPlpContribSPD = 1,
   Bool_t      useMVPileUpSelection = kFALSE,
   Int_t       MinPlpContribMV = 5,
   Int_t       aodFilterBit = 5,
   AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPiCandidate = AliRsnCutSetDaughterParticle::kTOFpidKstarPPB2011,
   AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTOFpidKstarPPB2011,
   Float_t     nsigmaPi = 2.0,
   Float_t     nsigmaKa = 2.0,
   Bool_t      enableMonitor = kTRUE,
   Bool_t      IsMcTrueOnly = kFALSE,
   Double_t    minYlab =  -0.465,
   Double_t    maxYlab =  0.035,
   Int_t       nmix = 10,
   Float_t     maxDiffVzMix = 1.0,
   Float_t     maxDiffMultMix = 10.0,
   Int_t       signedPdg = 313,
   TString     monitorOpt = "",
   Bool_t      useMixLS = 0,
   AliRsnMiniValue::EType yaxisvar = AliRsnMiniValue::kPt
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
   TString taskName = Form("TOFKStar%s%s_%i%i", (isPP? "pp" : "PbPb"), (isMC ? "MC" : "Data"), (Int_t)cutPiCandidate,(Int_t)cutKaCandidate );
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   if (isPP) 
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
   //if (isPP) cutVertex->SetCheckPileUp(kTRUE);   // set the check for pileup 
   
   //set check for pileup in 2013
   AliRsnCutEventUtils *cutEventUtils = new AliRsnCutEventUtils("cutEventUtils", rmFirstEvtChunk, rejectPileUp);
   if (useMVPileUpSelection){
     cutEventUtils->SetUseMVPlpSelection(useMVPileUpSelection);
     cutEventUtils->SetMinPlpContribMV(MinPlpContribMV);
     cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
     ::Info("AddAnalysisTaskTOFKStar", Form("Multiple-vtx Pile-up rejection settings: MinPlpContribMV = %i, MinPlpContribSPD = %i", MinPlpContribMV, MinPlpContribSPD));
   } else {
     cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
     ::Info("AddAnalysisTaskTOFKStar", Form("SPD Pile-up rejection settings: MinPlpContribSPD = %i", MinPlpContribSPD));
   }
   ::Info("AddAnalysisTaskTOFKStar", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp?"ON":"OFF")));   
   ::Info("AddAnalysisTaskTOFKStar", Form("::::::::::::: Remove first event in chunk: %s", (rmFirstEvtChunk?"ON":"OFF")));   
   
   // define and fill cut set for event cut
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutVertex);
   eventCuts->AddCut(cutEventUtils);
   eventCuts->SetCutScheme(Form("%s&%s", cutEventUtils->GetName(), cutVertex->GetName()));
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
   
   TH2F* hvz=new TH2F("hVzVsCent","",100,0.,100., 220,-11.,11.);
   task->SetEventQAHist("vz",hvz);//plugs this histogram into the fHAEventVz data member

   TH2F* hmc=new TH2F("MultiVsCent","",100,0.,100.,4000,0.,4000.);
   hmc->GetYaxis()->SetTitle("QUALITY");
   task->SetEventQAHist("multicent",hmc);//plugs this histogram into the fHAEventMultiCent data member

   //
   // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
   //
   AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
   cutY->SetRangeD(minYlab, maxYlab);
   
   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());
   
   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------
   //   
   gROOT->LoadMacro("$ALICE_ROOT/PWGLF/RESONANCES/macros/mini/ConfigTOFanalysisKStar.C");
   if (!ConfigTOFanalysisKStar(task, isMC, isPP, "", cutsPair, aodFilterBit, cutPiCandidate, cutKaCandidate, nsigmaPi, nsigmaKa, enableMonitor, isMC&IsMcTrueOnly, useMixLS, signedPdg, monitorOpt.Data(), yaxisvar)) return 0x0;
   
   
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
