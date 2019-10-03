/***************************************************************************
   sarita.sahoo@cern.ch - last modified on 10/02/2015
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

AliRsnMiniAnalysisTask * AddAnalysisTaskPhiMassStudy
(
   Bool_t      isMC,
   Bool_t      isPP,
   Int_t       aodFilterBit = 5,
   AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,
   Float_t     nsigmaKa = 2.0,
   Bool_t      enableSyst = kFALSE,
   Char_t      DCAxyFormula[100] = "0.0182+0.035/pt^1.01",
   Double_t    dcazmax = 2,
   Double_t    minNcls = 70,
   Double_t    maxX2cls = 4.0,
   Double_t    globalX2cls = 36.0,
   Double_t    minCrossedRows = 70.0,
   Double_t    maxClsCrossedRows = 0.8,
   Bool_t      enableMonitor = kTRUE,
   Bool_t      IsMcTrueOnly = kFALSE,
   UInt_t      triggerMask = AliVEvent::kMB,
   Int_t       signedPdg = 333,
   TString     monitorOpt = "NoSIGN",  //Flag for AddMonitorOutput.C e.g."NoSIGN"
   Bool_t      useCrossedRows = kFALSE,
   const char *yaxisVar = "PtDaughter_PDaughter_cent_hello",  //yaxisVar = "PtDaughter_PDaughter_cent"
   Bool_t      useMixLS = 0,

   Int_t       nmix = 5,
   Float_t     maxDiffVzMix = 1.0,
   Float_t     maxDiffMultMix = 10.0,
   TString     outNameSuffix = ""
)
{  
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddAnalysisTaskPhi", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   TString taskName = Form("Phi%s%s_%s", (isPP? "pp" : "pPb"), (isMC ? "MC" : "Data"), outNameSuffix.Data() );
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);


   //if(is2011PbPb)
   //task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
   //else
   task->SelectCollisionCandidates(triggerMask);


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
   //if (!isPP) task->SetMaxDiffAngle(maxDiffAngleMixDeg*TMath::DegToRad()); //set angle diff in rad
   ::Info("AddAnalysisTasPhi", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n ", nmix, maxDiffVzMix, maxDiffMultMix));
   
   mgr->AddTask(task);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //  
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
   if (isPP) {
     AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
     cutVertex->SetCheckPileUp(kTRUE);   // set the check for pileup
     
     // define and fill cut set for event cut
     AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
     eventCuts->AddCut(cutVertex);
     eventCuts->SetCutScheme(cutVertex->GetName());
   }
   else {
   AliRsnCutEventUtils *cutEventUtils = new AliRsnCutEventUtils("cutEventUtils", kTRUE, kTRUE);
   cutEventUtils->SetUseVertexSelection2013pA(kTRUE);
   cutEventUtils->SetMinPlpContribSPD(5);      
   // define and fill cut set for event cut
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutEventUtils);
   eventCuts->SetCutScheme(cutEventUtils->GetName());
   }

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
     outMult->AddAxis(multID, 300, 0.0, 300.0);
   else
     outMult->AddAxis(multID, 100, 0.0, 100.0);
   
   //
   // -- PAIR CUTS  -------------------------------------------------------------------------------
   //
   AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
   if(isPP) cutY->SetRangeD(-0.5, 0.5);
   else     cutY->SetRangeD(-0.465, 0.035);

   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());
   
   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------

   
   //for systematic checks
     {
       gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigPhiMassStudy.C");

   if (!ConfigPhiMassStudy(task, isMC, isPP, "", cutsPair, aodFilterBit, cutKaCandidate, nsigmaKa, enableSyst, DCAxyFormula, dcazmax, minNcls, maxX2cls, globalX2cls, minCrossedRows, maxClsCrossedRows, enableMonitor, isMC&IsMcTrueOnly, signedPdg,monitorOpt,useCrossedRows,yaxisVar ,useMixLS)) 


return 0x0;  
     }

   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddAnalysisTaskPhi - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}
