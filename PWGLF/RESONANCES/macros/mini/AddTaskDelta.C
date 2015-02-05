//
// General macro to configure the Delta++ analysis with RSN analysis package.
// It calls all configs desired by the user, by means
// of the boolean switches defined in the first lines.
// ---
// Inputs:
//  1) flag to know if running on MC or data
//  2) flag to know if running on pp or AA data
//  3) number of events in the mixing-event pool
//  4) flag to enable systematics
// ---
// Returns:
//  kTRUE  --> initialization successful
//  kFALSE --> initialization failed (some config gave errors)
//

//set to kTRUE if using data AOD049 - needed to enable centrality patch
Bool_t isAOD049 = 0;

AliRsnMiniAnalysisTask * AddTaskDelta
(
 Bool_t      isMC,
 Bool_t      isPP,
 Int_t       nmix,
 Bool_t      enableSyst = 0
 )
{  
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  //
   
  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  // create the task and connect with physics selection
  AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask("RSN", isMC);
  if (isAOD049 && !isMC && !isPP){
    task->SetUseCentralityPatch(kTRUE);
  }

  UInt_t triggerMask = AliVEvent::kMB; // minimum bias selection for pp
  if (isPP) {
    task->UseMultiplicity("QUALITY");
  } else { 
    //is pA
    triggerMask = AliVEvent::kINT7;
    task->UseCentrality("V0A");
  }
  //set trigger selection  
  task->SelectCollisionCandidates(triggerMask); //AOD
  // set mixing
  task->UseContinuousMix();
  //task->UseBinnedMix();
  task->SetNMix(nmix);
  task->SetMaxDiffVz(1.0);
  task->SetMaxDiffMult(5.0);
  task->SetMaxDiffAngle(1E20);
  mgr->AddTask(task);
   
  //
  // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
  //
  Bool_t      rmFirstEvtChunk = kTRUE; //needed for pA 2013
  Bool_t      rejectPileUp = kTRUE; //best if used, for pA 2013
  Int_t       MinPlpContribSPD = 5; //default value if used
  Bool_t      useMVPileUpSelection = kFALSE; //
  Int_t       MinPlpContribMV = 5; //default value if used
  Bool_t      useVtxCut2013pA = kTRUE; //default use recommended 2013 pA vtx selection
  Double_t    vtxZcut = 10.0; //cm, default cut on vtx z
  
  // define the event cut set
  AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent); 
  // cut on primary vertex:
  // - 2nd argument --> |Vz| range
  // - 3rd argument --> minimum required number of contributors
  // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
  AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", vtxZcut, 0, kFALSE);
  AliRsnCutEventUtils *cutEventUtils = 0x0;
  if (isPP) { 
    cutVertex->SetCheckPileUp(kTRUE);
    eventCuts->AddCut(cutVertex);
    eventCuts->SetCutScheme(Form("%s", cutVertex->GetName()));
  } else {
    //cutEventUtils configured and used only for pA
    cutEventUtils = new AliRsnCutEventUtils("cutEventUtils", rmFirstEvtChunk, rejectPileUp);
    cutEventUtils->SetUseVertexSelection2013pA(useVtxCut2013pA);
    ::Info("AddAnalysisTaskTOFDelta", Form(":::::::::::::::::: Vertex cut as pA 2013: %s", (useVtxCut2013pA?"ON":"OFF")));   
    if (useMVPileUpSelection){
      cutEventUtils->SetUseMVPlpSelection(useMVPileUpSelection);
      cutEventUtils->SetMinPlpContribMV(MinPlpContribMV);
      cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
      ::Info("AddAnalysisTaskTOFDelta", Form("Multiple-vtx Pile-up rejection settings: MinPlpContribMV = %i, MinPlpContribSPD = %i", MinPlpContribMV, MinPlpContribSPD));
    } else {
      cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
      ::Info("AddAnalysisTaskTOFDelta", Form("SPD Pile-up rejection settings: MinPlpContribSPD = %i", MinPlpContribSPD));
    }
    ::Info("AddAnalysisTaskTOFDelta", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp?"ON":"OFF")));   
    ::Info("AddAnalysisTaskTOFDelta", Form("::::::::::::: Remove first event in chunk: %s", (rmFirstEvtChunk?"ON":"OFF")));   
    //add event utils cuts
    eventCuts->AddCut(cutVertex);
    eventCuts->AddCut(cutEventUtils);
    //define cut scheme
    eventCuts->SetCutScheme(Form("%s&%s", cutEventUtils->GetName(), cutVertex->GetName()));
  } 
  
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
    outMult->AddAxis(multID, 400, 0.0, 400.0);
  else
    outMult->AddAxis(multID, 100, 0.0, 100.0);
   
  TH2F* hvz=new TH2F("hVzVsCent","", 100, 0., 100., 240, -12.0, 12.0);
  task->SetEventQAHist("vz",hvz);//plugs this histogram into the fHAEventVz data member

  TH2F* hmc=new TH2F("MultiVsCent","", 100, 0., 100., 400, 0., 400.);
  hmc->GetYaxis()->SetTitle("QUALITY");
  task->SetEventQAHist("multicent",hmc);//plugs this histogram into the fHAEventMultiCent data member
      
      
   
  //
  // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
  //
   
  AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  if (isPP) cutY->SetRangeD(-0.5, 0.5);
  else cutY->SetRangeD(-0.465, 0.035);  // 0 < y_cm < 0.5; y_cm = y_lab + 0.465
   
  AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
  cutsPair->AddCut(cutY);
  cutsPair->SetCutScheme(cutY->GetName());
   
  //
  // -- CONFIGS -----------------------------------------------------------------------------------
  //
  
  if (isPP){
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigDeltaPP7TeV.C");
    ConfigDeltaPP7TeV(task, isMC, "", cutsPair);
    if (isMC) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigDeltaPP7TeV_MC.C");
      ConfigDeltaPP7TeV_MC(task, isPP, "", cutsPair);
    } 
  } else {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigDeltaPPb.C");
    ConfigDeltaPPb(task, isMC, "", cutsPair); 
    if (isMC) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigDeltaPP7TeV_MC.C");
      ConfigDeltaPP7TeV_MC(task, isPP, "", cutsPair);
    } 
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
