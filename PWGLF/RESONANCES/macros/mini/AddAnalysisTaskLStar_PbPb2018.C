/***************************************************************************
// rama.chandra.baral@cern.ch & sarita.sahoo@cern.ch - last modified on 12/06/2014
// adapted for Pb-Pb analysis: Neelima Agrawal, Roberto Preghenella on 11/06/2016
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

AliRsnMiniAnalysisTask *
AddAnalysisTaskLStar_PbPb2018(
			      UInt_t      triggerMask       = AliVEvent::kINT7,
			      Float_t     yCut              = 0.5,
			      Int_t       aodFilterBit      = 5,
			      Bool_t      useTPCCrossedRows = kTRUE,
			      Int_t       qualityCut        = AliRsnCutSetDaughterParticle::kQualityStd2011,
			      Int_t       pidCut            = AliRsnCutSetDaughterParticle::kTPCTOFpidTunedPbPbTOFneed_2018,
			      Float_t     nsPr              = 1.0, // 3.0
			      Float_t     nsKa              = 1.0, // 3.0
			      Int_t       nMix              = 15,
			      Bool_t      isMC              = kFALSE, 
			      const char *suffix            = ""
			      )
{  
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddAnalysisTaskLStar_PbPb2018", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   TString taskName = "Lambda1520_PbPb";
   if (strlen(suffix) > 0) taskName += Form("_%s", suffix);
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   task->SelectCollisionCandidates(triggerMask);
   task->UseMultiplicity("AliMultSelection_V0M");
   // set event mixing options
   task->UseContinuousMix();
   //task->UseBinnedMix();
   task->SetNMix(nMix);
   task->SetMaxDiffVz(1.0);
   task->SetMaxDiffMult(10.);
   task->SetMaxDiffAngle(20.*TMath::DegToRad());
   mgr->AddTask(task);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //  
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
   //   cutVertex->SetCheckPileUp(kTRUE);   // set the check for pileup   
   // define and fill cut set for event cut
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutVertex);
   eventCuts->SetCutScheme(cutVertex->GetName());
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
   AliRsnMiniOutput *outMult = task->CreateOutput("eventCentrality", "HIST", "EVENT");
   outMult->AddAxis(multID, 100, 0.0, 100.0);

   //tracklets Vs centrality                                              
   Int_t trklID = task->CreateValue(AliRsnMiniValue::kTracklets, kFALSE);
   AliRsnMiniOutput *outTrcent = task->CreateOutput("trackletsCentrality", "HIST", "EVENT");
   outTrcent->AddAxis(multID, 100, 0.0, 100.0);
   outTrcent->AddAxis(trklID, 500, 0.0, 5000.);
   
   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------
   
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigLStar_PbPb2018.C");
   ConfigLStar_PbPb2018(task, yCut, aodFilterBit, useTPCCrossedRows, qualityCut, pidCut, nsPr, nsPr, isMC, suffix);
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s", taskName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}
