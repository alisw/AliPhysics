/***************************************************************************
  dmallick@lxplus.cern.ch, prottay@lxplus.cern.ch (18/12/2020)
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
AliRsnMiniAnalysisTask * AddTaskKStarPPRun2_HighPt
(
 Bool_t      isMC           = kFALSE,
 Bool_t      isPP           = kTRUE,
 Int_t       Strcut         = 2011,
 Int_t       MultId         = 0,
 TString     outNameSuffix  = "tofveto3stpc2s",
 Int_t       customQualityCutsID = AliRsnCutSetDaughterParticle::kDisableCustom,
 Bool_t                 enableMonitor = kTRUE,
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate=AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s,
 Float_t                ThetaStar =AliRsnMiniValue::kCosThetaStarAbs,
 UInt_t      triggerMask =   AliVEvent::kINT7,
 Float_t                nsigmaPi  = 2.0,
 Float_t                nsigmaK   = 2.0,
 Float_t                nsigmaTOF = 3.0,
 Int_t       nmix                = 5,
 Float_t     maxDiffVzMix        = 1.0,
 Float_t     maxDiffMultMix      = 5.0,
 Int_t                   Multbin=100,
 Int_t                   lMultbin=0,
 Int_t                   hMultbin=100,
 Int_t                   Ptbin=300,
 Float_t                 lPtbin=0.0,
 Float_t                   hPtbin=30.0,
 Int_t                   InvMbin=150,
 Float_t                 lInvMbin=0.6,
 Float_t                 hInvMbin=1.5,
 Int_t                   Cosbin=20,
 Float_t                 lCosbin=-1.0,
 Float_t                  hCosbin=1.0
 )
{  

  //  UInt_t      triggerMask  = AliVEvent::kINT7;
  Bool_t      rejectPileUp = kTRUE;
  Double_t    vtxZcut      = 10.0;//cm, default cut on vtx z
  if(!isPP || isMC) rejectPileUp = kFALSE;
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskKStarMultDepPP5TeV", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   
   TString taskName = Form("KStar%s%s", (isPP? "pp" : "PbPb"), (isMC ? "MC" : "Data"));
   
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   //task->SelectCollisionCandidates(triggerMask);
   task->UseESDTriggerMask(triggerMask);
   if(isPP) {
     if(MultId==1) task->UseMultiplicity("AliMultSelection_V0M");
     else if(MultId==2) task->UseMultiplicity("AliMultSelection_RefMult08");
     else task->UseMultiplicity("QUALITY");
   } else task->UseCentrality("V0M");   
   // set event mixing options
   task->UseContinuousMix();
   //task->UseBinnedMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   task->UseMC(isMC);
   ::Info("AddTaskKStarMultDepPP5TeV", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n", nmix, maxDiffVzMix, maxDiffMultMix));
   
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
   if(MultId==0){
     cutVertex->SetCheckZResolutionSPD();
     cutVertex->SetCheckDispersionSPD();
     cutVertex->SetCheckZDifferenceSPDTrack();
   }
   
   AliRsnCutEventUtils* cutEventUtils=new AliRsnCutEventUtils("cutEventUtils",kTRUE,rejectPileUp);
   if(MultId==0){
     cutEventUtils->SetCheckIncompleteDAQ();
     cutEventUtils->SetCheckSPDClusterVsTrackletBG();
   }
   else if(MultId==1 || MultId==2){
     cutEventUtils->SetRemovePileUppA2013(kFALSE);
     cutEventUtils->SetCheckAcceptedMultSelection();
   }

     
   // define and fill cut set for event cut
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutEventUtils);
   eventCuts->AddCut(cutVertex);
   eventCuts->SetCutScheme(Form("%s&%s",cutEventUtils->GetName(),cutVertex->GetName()));
   // set cuts in task
   task->SetEventCuts(eventCuts);
   
   //
   // -- EVENT-ONLY COMPUTATIONS -------------------------------------------------------------------
   //   
   //vertex
   Int_t vtxID = task->CreateValue(AliRsnMiniValue::kVz, kFALSE);
   AliRsnMiniOutput *outVtx = task->CreateOutput("eventVtx", "HIST", "EVENT");
   outVtx->AddAxis(vtxID, 240, -12.0, 12.0);
   
   //multiplicity or centrality
   Int_t multID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
   AliRsnMiniOutput *outMult = task->CreateOutput("eventMult", "HIST", "EVENT");
   if (isPP) 
     outMult->AddAxis(multID, 100, 0.0, 100.0);
   else
     outMult->AddAxis(multID, 100, 0.0, 100.0);
   
   //
   // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
   //

   AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
   cutY->SetRangeD(-0.5, 0.5);
   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());
   
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigKStarPPRun2_HighPt.C");
      if (!ConfigKStarPPRun2_HighPt(task, isMC, isPP, "", cutsPair, Strcut, customQualityCutsID,kTRUE,cutKaCandidate,ThetaStar,nsigmaPi,nsigmaK, nsigmaTOF, Multbin, lMultbin,hMultbin, Ptbin,lPtbin,hPtbin,InvMbin,lInvMbin,hInvMbin,Cosbin, lCosbin, hCosbin)) return 0x0;

   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddTaskKStarMultDepPP5TeV - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}
