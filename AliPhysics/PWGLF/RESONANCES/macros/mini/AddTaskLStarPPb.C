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

Bool_t useLStar = 1;


AliRsnMiniAnalysisTask * AddTaskLStarPPb
(
   Bool_t      isMC,
   Bool_t      isPP,
   const char *path,
   Int_t       nmix = 0
)
{  
   //
   // -- INITIALIZATION ----------------------------------------------------------------------------
   //
   
  //-------------------------------------------
  // event cuts
  //-------------------------------------------
  UInt_t      triggerMask = AliVEvent::kINT7;
  Bool_t      rmFirstEvtChunk = kTRUE; //needed for pA 2013
  Bool_t      rejectPileUp = kTRUE; //best if used, for pA 2013
  Int_t       MinPlpContribSPD = 5; //default value if used
  Bool_t      useMVPileUpSelection = kFALSE; //
  Int_t       MinPlpContribMV = 5; //default value if used
  Bool_t      useVtxCut2013pA = kTRUE; //default use recommended 2013 pA vtx selection
  Double_t    vtxZcut = 10.0; //cm, default cut on vtx z



   // retrieve analysis manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

   // create the task and connect with physics selection
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask("RSN", isMC);

   //task->UseESDTriggerMask(triggerMask); //ESD
   task->SelectCollisionCandidates(triggerMask); //AOD	
   
   // settings
   if (isPP) 
      task->UseMultiplicity("QUALITY");
   else
      task->UseCentrality("V0M");
   
   // set mixing
   task->UseContinuousMix();
   //task->UseBinnedMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(1.0);
   task->SetMaxDiffMult(10.0);
   //task->SetMaxDiffAngle(20);
   //   if (!isPP) task->SetMaxDiffAngle(20.0*TMath::DegToRad()); //set angle diff in rad

   mgr->AddTask(task);

   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //
   
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", vtxZcut, 0, kFALSE);
   
   // set the check for pileup
   //if (isPP) cutVertex->SetCheckPileUp(kTRUE);

   //set check for pileup in 2013  #############
 
   AliRsnCutEventUtils *cutEventUtils = new AliRsnCutEventUtils("cutEventUtils", rmFirstEvtChunk, rejectPileUp);
   cutEventUtils->SetUseVertexSelection2013pA(useVtxCut2013pA);
   ::Info("AddTaskLStarPPb", Form(":::::::::::::::::: Vertex cut as pA 2013: %s", (useVtxCut2013pA?"ON":"OFF")));   
   if (useMVPileUpSelection){
     cutEventUtils->SetUseMVPlpSelection(useMVPileUpSelection);
     cutEventUtils->SetMinPlpContribMV(MinPlpContribMV);
     cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
     ::Info("AddTaskLStarPPb", Form("Multiple-vtx Pile-up rejection settings: MinPlpContribMV = %i, MinPlpContribSPD = %i", MinPlpContribMV, MinPlpContribSPD));
   } else {
     cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
     ::Info("AddTaskLStarPPb", Form("SPD Pile-up rejection settings: MinPlpContribSPD = %i", MinPlpContribSPD));
   }
   ::Info("AddTaskLStarPPb", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp?"ON":"OFF")));   
   ::Info("AddTaskLStarPPb", Form("::::::::::::: Remove first event in chunk: %s", (rmFirstEvtChunk?"ON":"OFF")));   
  



   // define and fill cut set
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutEventUtils);
   eventCuts->AddCut(cutVertex);
   eventCuts->SetCutScheme(Form("%s&%s", cutEventUtils->GetName(), cutVertex->GetName()));
   //eventCuts->SetCutScheme(cutVertex->GetName());
   
   // set cuts in task
   task->SetEventCuts(eventCuts);
   
   //
   // -- EVENT-ONLY COMPUTATIONS -------------------------------------------------------------------
   //
   
   // second argument tells if the value must be taken from MC
   // (when this can be done)
   // after creating the value, the task returns its ID

   // VERTEX POSITION    
   Int_t vertexID = task->CreateValue(AliRsnMiniValue::kVz, kFALSE);
   AliRsnMiniOutput *outVertex = task->CreateOutput("eventVertex", "HIST", "EVENT");
   outVertex->AddAxis(vertexID, 400, -20.0, 20.0);


   // initialize value computation for multiplicity/centrality
   Int_t multID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
   AliRsnMiniOutput *outMult = task->CreateOutput("eventMult", "HIST", "EVENT");
   // set axes, by passing value ID and defining the binning
   if (isPP) 
      outMult->AddAxis(multID, 300, 0.0, 300.0);
   else
      outMult->AddAxis(multID, 100, 0.0, 100.0);
   
//   TH2F* hvz=new TH2F("hVzVsCent","", 100, 0., 100., 240, -12.0, 12.0);
//   task->SetEventQAHist("vz",hvz);//plugs this histogram into the fHAEventVz data member

//   TH2F* hmc=new TH2F("MultiVsCent","", 400, 0., 400., 400, 0., 400.);
//   hmc->GetYaxis()->SetTitle("QUALITY");
//   task->SetEventQAHist("multicent",hmc);//plugs this histogram into the fHAEventMultiCent data member


   //
   // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
   //
   
   AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
   //cutY->SetRangeD(-0.5, 0.5);
      cutY->SetRangeD(-0.465, 0.035);// 0 < y_cm < 0.5; y_cm = y_lab + 0.465

   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());
   
   //
   // -- CONFIGS -----------------------------------------------------------------------------------
   //
   

   if (useLStar) {
      gROOT->LoadMacro(Form("%s/ConfigLStarPPb.C", path));
      ConfigLStarPPb(task, isMC, isPP, "LStar", cutsPair);
   }
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //

   const char *file = AliAnalysisManager::GetCommonFileName();
   AliAnalysisDataContainer *output = mgr->CreateContainer("RsnOut", TList::Class(), AliAnalysisManager::kOutputContainer, file);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);

   return task;
}
