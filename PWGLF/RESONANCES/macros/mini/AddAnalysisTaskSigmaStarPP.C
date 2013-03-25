/***************************************************************************
              massimo.venaruzzo@cern.ch - last modified on 11/10/2012
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

AliRsnMiniAnalysisTask * AddAnalysisTaskSigmaStarPP
(
   Bool_t      isMC,
   Bool_t      isPP,
   Float_t     cutV = 10.0,
   Int_t       aodFilterBit = 5,
   Int_t       piPIDCut = 3.0,
   Int_t       pPIDCut = 3.0,
   Float_t     piDCAcut = 0.05,
   Float_t     massTol = 0.01,
   Float_t     lambdaDCA = 0.03,
   Float_t     lambdaCosPoinAn = 0.99,
   Float_t     lambdaDaughDCA = 0.5,
   Int_t       NTPCcluster = 70,
   Int_t       nmix = 5,
   Float_t     maxDiffVzMix = 1.0,
   Float_t     maxDiffMultMix = 10.0,
   Float_t     maxDiffAngleMixDeg = 20.0,
   Int_t       aodN = 0,
   TString     outNameSuffix = "Sigma1385"
)
{  
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddAnalysisTaskSigmaStarPP", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   TString taskName = Form("SigmaStar%s%s_%i%i", (isPP? "pp" : "PbPb"), (isMC ? "MC" : "Data"));
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   if (!isMC && !isPP){
     Printf(Form("========== SETTING USE CENTRALITY PATCH AOD049 : %s", (aodN==49)? "yes" : "no"));
     task->SetUseCentralityPatch(aodN==49);
   }
   if (isPP) 
     task->UseMultiplicity("QUALITY");
   else
     task->UseCentrality("V0M");   
   // set event mixing options
   task->UseContinuousMix();
   //task->UseBinnedMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   if (!isPP) task->SetMaxDiffAngle(maxDiffAngleMixDeg*TMath::DegToRad()); //set angle diff in rad
   ::Info("AddAnalysisTaskSigmaStarPP", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n max diff EP angle = %5.3f deg", nmix, maxDiffVzMix, maxDiffMultMix, (isPP ? 0.0 : maxDiffAngleMixDeg)));
   
   mgr->AddTask(task);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //  
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", cutV, 0, kFALSE);
   if (isPP) cutVertex->SetCheckPileUp(kTRUE);   // set the check for pileup
   
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
   if (isPP) 
     outMult->AddAxis(multID, 400, 0.0, 400.0);
   else
     outMult->AddAxis(multID, 100, 0.0, 100.0);
   
   //event plane (only for PbPb)
   Int_t planeID = task->CreateValue(AliRsnMiniValue::kPlaneAngle, kFALSE);
   AliRsnMiniOutput *outPlane = 0x0; //task->CreateOutput("eventPlane", "HIST", "EVENT");
   if (!isPP){
     outPlane = task->CreateOutput("eventPlane", "HIST", "EVENT");
     outPlane->AddAxis(planeID, 180, 0.0, TMath::Pi());
     }
   
   //
   // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
   //
   AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
   cutY->SetRangeD(-0.5, 0.5);
   
   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());
   
   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------
   gROOT->LoadMacro("$ALICE_ROOT/PWGLF/RESONANCES/macros/mini/ConfigSigmaStar.C");
   if (isMC) {
       Printf("========================== MC analysis - PID cuts not used");
   } else 
     Printf("========================== DATA analysis - PID cuts used");
   if (!ConfigSigmaStar(task, isPP, isMC, piPIDCut, pPIDCut, aodFilterBit, piDCAcut, massTol, lambdaDCA, lambdaCosPoinAn, lambdaDaughDCA, NTPCcluster, "", cutsPair)) return 0x0;
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddAnalysisTaskSigmaStarPP - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}
