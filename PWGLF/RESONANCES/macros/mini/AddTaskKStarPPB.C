/***************************************************************************
              fbellini@cern.ch - last modified on 28/11/2013
//
//Lauches KStar analysis with rsn mini package
//Allows basic configuration of pile-up check and event cuts
//
****************************************************************************/
enum pairYCutSet { kPairDefault,
		   kNegative,
		   kCentral,
		   kWide1,
		   kWide2
                 };

enum eventCutSet { kEvtDefault = 0,
		   kNSD,
		   kNSDpA2013,
		   kNSDpA2013spectra,
		   kNSDpA2013DefaultSpectra,
		   kNoPileUpCut,
		   kPileUpMV,
		   kPileUpSPD3,		      
		   kDefaultVtx8,
		   kDefaultVtx5,
		   kMCEvt,
		   kMCEvtpA2013,
		   kMCEvtDefault,
		   kMCEvtNSD,
		   kMCEvtNSDpA2013, 
		   kMCEvtNSDdefault
};

enum eventMixConfig { kDisabled = -1,
		      kMixDefault,//=0 //10 events, Dvz = 1cm, DC = 10
		      k5Evts, //=1 //5 events, Dvz = 1cm, DC = 10
		      k5Cent,  //=2 //10 events, Dvz = 1cm, DC = 5
                    };


AliRsnMiniAnalysisTask * AddTaskKStarPPB
(
 Bool_t      isMC = kFALSE,
 Bool_t      isPP = kFALSE,
 TString     outNameSuffix = "tof2s",
 Int_t       evtCutSetID = 0,
 Int_t       pairCutSetID = 0,
 Int_t       mixingConfigID = 0,
 Int_t       aodFilterBit = 5,
 Int_t       customQualityCutsID = -1,
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPiCandidate = AliRsnCutSetDaughterParticle::kTOFpidKstarPPB2011,
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTOFpidKstarPPB2011,
 Float_t     nsigmaPi = 2.0,
 Float_t     nsigmaKa = 2.0,
 Bool_t      enableMonitor = kTRUE,
 Bool_t      IsMcTrueOnly = kFALSE,
 TString     monitorOpt = "NoSIGN",
 Bool_t      useMixLS = 0,
 Bool_t      checkReflex = 0,
 AliRsnMiniValue::EType yaxisvar = AliRsnMiniValue::kPt
)
{  

  
  //-------------------------------------------
  // event cuts
  //-------------------------------------------
  UInt_t      triggerMask = AliVEvent::kINT7;
  Bool_t      useVtxCut2013pA = kTRUE; //default use recommended 2013 pA vtx selection from AliAnalysisUtils
  Bool_t      useVtxCut2013pAspectra = kFALSE; //set to use 2013 pA vtx selection applied for pi/K/p analysis
  Bool_t      rmFirstEvtChunk = kTRUE; //needed for pA 2013
  Bool_t      rejectPileUp = kTRUE; //best if used, for pA 2013
  Int_t       MinPlpContribSPD = 5; //default value if used
  Bool_t      useMVPileUpSelection = kFALSE; //
  Int_t       MinPlpContribMV = 5; //default value if used
  Double_t    vtxZcut = 10.0; //cm, default cut on vtx z
  Bool_t      selectDPMJETevtNSDpA = kFALSE; //cut to select in DPMJET true NSD events

  if (evtCutSetID==eventCutSet::kEvtDefault) {
    useVtxCut2013pA = kTRUE;
    useVtxCut2013pAspectra = kFALSE;
    rmFirstEvtChunk = kTRUE;
    rejectPileUp = kTRUE;
    vtxZcut = 10.0; //cm
  }

  if (evtCutSetID==eventCutSet::kNSD) {
    useVtxCut2013pA = kFALSE;
    useVtxCut2013pAspectra = kFALSE;
    rejectPileUp = kFALSE;
    vtxZcut = 1.0e3; //cm
  }
  
  if (evtCutSetID==eventCutSet::kNSDpA2013) {
    useVtxCut2013pA = kFALSE;
    useVtxCut2013pAspectra = kFALSE;
    rejectPileUp = kFALSE;
    vtxZcut = 1.0e3; //cm
  }
  
  if (evtCutSetID==eventCutSet::kNSDpA2013spectra) {
    useVtxCut2013pAspectra = kTRUE;
    useVtxCut2013pA = kFALSE;
    vtxZcut = 1.0e3; //cm
    rejectPileUp = kFALSE; 
  }

  if (evtCutSetID==eventCutSet::kNSDpA2013DefaultSpectra) {
    useVtxCut2013pAspectra = kTRUE;
    useVtxCut2013pA = kFALSE;
    vtxZcut = 10.0; //cm
    rejectPileUp = kFALSE; 
  }

  if (evtCutSetID==eventCutSet::kNoPileUpCut) {
    rmFirstEvtChunk = kTRUE;
    rejectPileUp = kFALSE;
  }
  
  if (evtCutSetID==eventCutSet::kPileUpMV) {
    useMVPileUpSelection = kTRUE;
    MinPlpContribSPD = 3;
    //MinPlpContribMV = 5; //already set as default
  }
  
  if (evtCutSetID==eventCutSet::kPileUpSPD3) {
    MinPlpContribSPD = 3;
  }
  
  if (evtCutSetID==eventCutSet::kDefaultVtx8){
    vtxZcut = 8.0; //cm
  } 
  
  if (evtCutSetID==eventCutSet::kDefaultVtx5){
    vtxZcut = 5.0; //cm
  }
  
  if (evtCutSetID==eventCutSet::kMCEvt) {
    rmFirstEvtChunk = kFALSE;
    rejectPileUp = kFALSE;
    useVtxCut2013pA = kFALSE;
    vtxZcut = 1.0e3; //cm
    selectDPMJETevtNSDpA = kFALSE;
  }

  if (evtCutSetID==eventCutSet::kMCEvtpA2013) {
    rmFirstEvtChunk = kFALSE;
    rejectPileUp = kFALSE;
    useVtxCut2013pA = kTRUE;
    vtxZcut = 1.0e3; //cm
    selectDPMJETevtNSDpA = kFALSE;
  }
  
  if (evtCutSetID==eventCutSet::kMCEvtDefault) {
    rmFirstEvtChunk = kFALSE;
    rejectPileUp = kFALSE;
    useVtxCut2013pA = kTRUE;
    vtxZcut = 10.0; //cm
    selectDPMJETevtNSDpA = kFALSE;
  }
  
  if (evtCutSetID>=eventCutSet::kMCEvtNSD) {
    rmFirstEvtChunk = kFALSE;
    rejectPileUp = kFALSE;
    useVtxCut2013pA = kFALSE;
    vtxZcut = 1.0e3; //cm
    selectDPMJETevtNSDpA = kTRUE;
  }
  
  if (evtCutSetID==eventCutSet::kMCEvtNSDpA2013) {
    rmFirstEvtChunk = kFALSE;
    rejectPileUp = kFALSE;
    vtxZcut = 1.0e3; //cm
    selectDPMJETevtNSDpA = kTRUE;
    useVtxCut2013pA = kTRUE;
  }
  
  if (evtCutSetID==eventCutSet::kMCEvtNSDdefault) {
    rmFirstEvtChunk = kFALSE;
    rejectPileUp = kFALSE;
    selectDPMJETevtNSDpA = kTRUE;
    useVtxCut2013pA = kTRUE;
    vtxZcut = 10.0; //cm
  }
  
  //-------------------------------------------
  //pair cuts
  //-------------------------------------------
  Double_t    minYlab =  -0.465;
  Double_t    maxYlab =  0.035;

  if (pairCutSetID==pairYCutSet::kNegative) { //-0.5<y_cm<0.0
    minYlab = -0.965;    maxYlab = -0.465;
  }
  
  if (pairCutSetID==pairYCutSet::kCentral) { //|y_cm|<0.3
    minYlab = -0.765;    maxYlab = -0.165;
  }

  if (pairCutSetID==pairYCutSet::kWide1) { //|y_cm|<1
    minYlab = -1.465;    maxYlab = 0.535;
  }

  if (pairCutSetID==pairYCutSet::kWide2) { //|y_cm|<2
    minYlab = -2.465;    maxYlab = -1.535;
  }
  //-------------------------------------------
  //mixing settings
  //-------------------------------------------
  Int_t       nmix = 0;
  Float_t     maxDiffVzMix = 1.0;
  Float_t     maxDiffMultMix = 10.0;
  
  if (mixingConfigID == eventMixConfig::kMixDefault) {
    nmix = 10;
  }

  if (mixingConfigID == eventMixConfig::k5Evts) {
    nmix = 5;
  }
  
  if (mixingConfigID == eventMixConfig::k5Cent) {
    maxDiffMultMix = 5;
  }

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
   //task->UseESDTriggerMask(triggerMask); //ESD
   task->SelectCollisionCandidates(triggerMask); //AOD
   
   if (isPP) 
     task->UseMultiplicity("QUALITY");
   else
     task->UseCentrality("V0A");   
   // set event mixing options
   task->UseContinuousMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   ::Info("AddAnalysisTaskTOFKStar", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f", nmix, maxDiffVzMix, maxDiffMultMix));
   
   mgr->AddTask(task);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //  
   // cut on primary vertex:
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", vtxZcut, 0, kFALSE);
   //if (isPP) cutVertex->SetCheckPileUp(kTRUE);   // set the check for pileup 
   //set check for pileup in 2013
   AliRsnCutEventUtils *cutEventUtils = new AliRsnCutEventUtils("cutEventUtils", rmFirstEvtChunk, rejectPileUp);
   cutEventUtils->SetUseVertexSelection2013pA(useVtxCut2013pA, vtxZcut);
   ::Info("AddTaskKStarPPB", Form(":::::::::::::::::: Vertex cut as pA 2013 (max Vz = %4.2f cm): %s", vtxZcut, (useVtxCut2013pA?"ON":"OFF")));  
   cutEventUtils->SetUseVertexSelection2013pAIDspectra(useVtxCut2013pAspectra, vtxZcut);
   ::Info("AddTaskKStarPPB", Form(":::::::::::::::::: Vertex cut as pA 2013 (max Vz = %4.2f cm): %s", vtxZcut, (useVtxCut2013pA?"ON":"OFF")));  
   
   if (isMC) {
     cutEventUtils->SetFilterNSDeventsDPMJETpA2013(selectDPMJETevtNSDpA);
     ::Info("AddTaskKStarPPB", Form(":::::::::::::::::: NSD selection in DPMJET pA: %s", (selectDPMJETevtNSDpA?"ON":"OFF")));  
   }
   
   if (useMVPileUpSelection){
     cutEventUtils->SetUseMVPlpSelection(useMVPileUpSelection);
     cutEventUtils->SetMinPlpContribMV(MinPlpContribMV);
     cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
     ::Info("AddTaskKStarPPB", Form("Multiple-vtx Pile-up rejection settings: MinPlpContribMV = %i, MinPlpContribSPD = %i", MinPlpContribMV, MinPlpContribSPD));
   } else {
     cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
     ::Info("AddTaskKStarPPB", Form("SPD Pile-up rejection settings: MinPlpContribSPD = %i", MinPlpContribSPD));
   }
   ::Info("AddTaskKStarPPB", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp?"ON":"OFF")));   
   ::Info("AddTaskKStarPPB", Form("::::::::::::: Remove first event in chunk: %s", (rmFirstEvtChunk?"ON":"OFF")));   
   
   // define and fill cut set for event cut
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutEventUtils);
   if (!useVtxCut2013pAspectra && !useVtxCut2013pA) {
     eventCuts->AddCut(cutVertex);
     eventCuts->SetCutScheme(Form("%s&%s", cutEventUtils->GetName(), cutVertex->GetName()));
   } else {
     eventCuts->AddCut(cutEventUtils);
     eventCuts->SetCutScheme(Form("%s", cutEventUtils->GetName()));
   }
   
     // set cuts in task
   task->SetEventCuts(eventCuts);
   
   //
   // -- EVENT-ONLY COMPUTATIONS -------------------------------------------------------------------
   //   
   //vertex
   Int_t vtxID = task->CreateValue(AliRsnMiniValue::kVz, kFALSE);
   //multiplicity or centrality
   Int_t multID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
   //reference multiplicity (default with global tracks with good quality, if not available uses tracklets)
   Int_t multRefID = task->CreateValue(AliRsnMiniValue::kRefMult, kFALSE);

   AliRsnMiniOutput *outVtx = task->CreateOutput("eventVtx", "HIST", "EVENT");
   outVtx->AddAxis(vtxID, 500, -50.0, 50.0);
   
   AliRsnMiniOutput *outMult = task->CreateOutput("eventMult", "HIST", "EVENT");
   if (isPP) 
     outMult->AddAxis(multID, 400, 0.0, 400.0);
   else
     outMult->AddAxis(multID, 101, 0.0, 101.0);
   
   AliRsnMiniOutput *outRefMult = task->CreateOutput("eventRefMult", "HIST", "EVENT");
   outRefMult->AddAxis(multRefID, 400, 0.0, 400.0);
   
   TH2F* hvz = new TH2F("hVzVsCent",Form("Vertex position vs centrality"), 101, 0., 101., 500, -50.0, 50.0);
   hvz->GetXaxis()->SetTitle("V0A");
   hvz->GetYaxis()->SetTitle("z_{vtx} (cm)");
   task->SetEventQAHist("vz", hvz);//plugs this histogram into the fHAEventVz data member

   TH2F* hRefMultiVsCent = new TH2F("hRefMultiVsCent",Form("Reference multiplicity vs centrality"), 101, 0., 101., 400, 0., 400.);
   hRefMultiVsCent->GetXaxis()->SetTitle("V0A");
   hRefMultiVsCent->GetYaxis()->SetTitle("GLOBAL");
   task->SetEventQAHist("refmulti",hRefMultiVsCent);//plugs this histogram into the fHAEventRefMultiCent data member

   TH2F* hMultiVsCent = new TH2F("hMultiVsCent",Form("Multiplicity vs centrality"), 101, 0., 101., 400, 0., 400.);
   hMultiVsCent->GetXaxis()->SetTitle("V0A");
   hMultiVsCent->GetYaxis()->SetTitle("QUALITY");
   task->SetEventQAHist("multicent",hMultiVsCent);//plugs this histogram into the fHAEventMultiCent data member
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
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigKStarPPb.C");
   if (!ConfigKStarPPb(task, isMC, isPP, "", cutsPair, aodFilterBit, customQualityCutsID, cutPiCandidate, cutKaCandidate, nsigmaPi, nsigmaKa, enableMonitor, isMC&IsMcTrueOnly,  monitorOpt.Data(), useMixLS, isMC&checkReflex, yaxisvar)) return 0x0;
   
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddTaskKStarPPB - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}
