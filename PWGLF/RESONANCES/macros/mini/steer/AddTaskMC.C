/***************************************************************************
              fbellini@cern.ch - last modified on 28/11/2013
//
//Lauches KStar analysis with rsn mini package
//Allows basic configuration of pile-up check and event cuts
//
****************************************************************************/
enum pairYCutSet { kPairDefault,
		   kpADefault,
		   kpANegative,
		   kpACentral,
		   kCentralTight
                 };

enum eventCutSet { kOld = -1, 
		   kEvtDefault,
		   kDefaultVtx8,
		   kDefaultVtx5,
		   kMCEvt,
		   kMCEvtDefault,
		   kpAEvtDefault,		
		   kpANoPileUpCut,
		   kpAPileUpMV,
		   kpAPileUpSPD3,		      
		   kMCEvtpA2013,
		   kMCpAEvtDefault,
		   kMCEvtNSD,
		   kMCEvtNSDpA2013, 
		   kMCEvtNSDdefault
};

AliRsnMiniAnalysisTask * AddTaskMC
(
 Bool_t      isPP = kFALSE,
 TString     outNameSuffix = "q2010",
 Int_t       evtCutSetID = 0,
 Int_t       pairCutSetID = 0,
 Int_t       aodFilterBit = 5,
 TString                partname="kStar",
 Int_t                  pdgCode=313,
 Float_t                mass = 0.89445,
 Float_t                masslow = 0.7,
 Float_t                massup = 1.1,
 Int_t                  nbins = 400,
 Char_t                 charge1 = '+',
 Char_t                 charge2 = '-',
 RSNPID                 d1 = AliRsnDaughter::kKaon,
 RSNPID                 d2 = AliRsnDaughter::kPion,
 Bool_t                 enableMonitor = kTRUE,
 TString                monitorOpt = "NoSIGN")
{  

  //-------------------------------------------
  // event cuts
  //-------------------------------------------
  UInt_t      triggerMask = AliVEvent::kMB;
  Bool_t      rmFirstEvtChunk = kFALSE; //needed for pA 2013
  Bool_t      rejectPileUp = kFALSE; //best if used, for pA 2013
  Int_t       MinPlpContribSPD = 5; //default value if used
  Bool_t      useMVPileUpSelection = kFALSE; //
  Int_t       MinPlpContribMV = 5; //default value if used
  Bool_t      useVtxCut2013pA = kFALSE; //default use recommended 2013 pA vtx selection
  Double_t    vtxZcut = 10.0; //cm, default cut on vtx z
  Bool_t      selectDPMJETevtNSDpA = kFALSE; //cut to select in DPMJET true NSD events
  
  if (evtCutSetID==eventCutSet::kOld) {
    triggerMask = AliVEvent::kAnyINT;
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

 if (evtCutSetID==eventCutSet::kMCEvtDefault) {
   rmFirstEvtChunk = kFALSE;
   rejectPileUp = kFALSE;
   useVtxCut2013pA = kFALSE;
   vtxZcut = 10.0; //cm
   selectDPMJETevtNSDpA = kFALSE;
 }
 
 if (evtCutSetID>=eventCutSet::kpAEvtDefault) {
   triggerMask = AliVEvent::kINT7;
   rmFirstEvtChunk = kTRUE;
   rejectPileUp = kTRUE;
   useVtxCut2013pA = kTRUE;
 }
 
 if (evtCutSetID==eventCutSet::kpANoPileUpCut) {
   rmFirstEvtChunk = kTRUE;
   rejectPileUp = kFALSE;
 }
  
 if (evtCutSetID==eventCutSet::kpAPileUpMV) {
   useMVPileUpSelection = kTRUE;
   MinPlpContribSPD = 3;
   //MinPlpContribMV = 5; //already set as default
 }
  
 if (evtCutSetID==eventCutSet::kpAPileUpSPD3) {
   MinPlpContribSPD = 3;
 }
  
 if (evtCutSetID==eventCutSet::kMCEvtpA2013) {
   rmFirstEvtChunk = kFALSE;
   rejectPileUp = kFALSE;
   useVtxCut2013pA = kTRUE;
   vtxZcut = 1.0e3; //cm
   selectDPMJETevtNSDpA = kFALSE;
 }
  
 if (evtCutSetID==eventCutSet::kMCpAEvtDefault) {
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
  Double_t    minYlab =  -0.5;
  Double_t    maxYlab =  0.5;

  if (pairCutSetID==pairYCutSet::kpADefault) { //-0.5<y_cm<0.0
    minYlab = -0.465;    maxYlab = 0.035;
  }
  
  if (pairCutSetID==pairYCutSet::kpANegative) { //-0.5<y_cm<0.0
    minYlab = -0.965;    maxYlab = -0.465;
  }
  
  if (pairCutSetID==pairYCutSet::kpACentral) { //|y_cm|<0.3
    minYlab = -0.765;    maxYlab = -0.165;
  }
  
  if (pairCutSetID==pairYCutSet::kCentralTight) { //|y_cm|<0.3
    minYlab = -0.3;    maxYlab = 0.3;
  }

  //-------------------------------------------
  //mixing settings
  //-------------------------------------------
  Int_t       nmix = 5;
  Float_t     maxDiffVzMix = 1.0;
  Float_t     maxDiffMultMix = 5.0;

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
   TString taskName = Form("taskRsnMC");
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), kTRUE);
   task->UseESDTriggerMask(triggerMask); //ESD
   // task->SelectCollisionCandidates(triggerMask); //AOD
   
   if (isPP) 
     task->UseMultiplicity("QUALITY");
   else
     task->UseCentrality("V0M");   
   // set event mixing options
   task->UseContinuousMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   ::Info("AddTaskMC", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f", nmix, maxDiffVzMix, maxDiffMultMix));
   mgr->AddTask(task);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", vtxZcut, 0, kFALSE);
   //if (isPP) cutVertex->SetCheckPileUp(kTRUE);   // set the check for pileup 
   
   //set check for pileup in 2013
   AliRsnCutEventUtils *cutEventUtils = new AliRsnCutEventUtils("cutEventUtils", rmFirstEvtChunk, rejectPileUp);
   cutEventUtils->SetUseVertexSelection2013pA(useVtxCut2013pA, vtxZcut);
   ::Info("AddTaskKStarPPB", Form(":::::::::::::::::: Vertex cut as pA 2013 (max Vz = %4.2f cm): %s", vtxZcut, (useVtxCut2013pA?"ON":"OFF")));  
   cutEventUtils->SetFilterNSDeventsDPMJETpA2013(selectDPMJETevtNSDpA);
   ::Info("AddTaskKStarPPB", Form(":::::::::::::::::: NSD selection in DPMJET pA: %s", (selectDPMJETevtNSDpA?"ON":"OFF")));  
 
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
   
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutEventUtils);
   eventCuts->AddCut(cutVertex);
   eventCuts->SetCutScheme(Form("%s&%s", cutEventUtils->GetName(), cutVertex->GetName()));
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
   task->SetEventQAHist("vz", hvz);

   TH2F* hRefMultiVsCent = new TH2F("hRefMultiVsCent",Form("Reference multiplicity vs centrality"), 101, 0., 101., 400, 0., 400.);
   hRefMultiVsCent->GetXaxis()->SetTitle("V0A");
   hRefMultiVsCent->GetYaxis()->SetTitle("GLOBAL");
   task->SetEventQAHist("refmulti",hRefMultiVsCent);

   TH2F* hMultiVsCent = new TH2F("hMultiVsCent",Form("Multiplicity vs centrality"), 101, 0., 101., 400, 0., 400.);
   hMultiVsCent->GetXaxis()->SetTitle("V0A");
   hMultiVsCent->GetYaxis()->SetTitle("QUALITY");
   task->SetEventQAHist("multicent",hMultiVsCent);

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
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/steer/ConfigRsnMC.C");
   if (!ConfigRsnMC(task, kTRUE, isPP, cutsPair, 
		    partname.Data(), pdgCode, mass, masslow, massup, nbins, charge1, charge2, 
		    d1, d2, aodFilterBit, enableMonitor, "NoSIGN") ) return 0x0;
   
   // -- CONTAINERS --------------------------------------------------------------------------------
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   Printf("AddTaskMC - Set OutputFileName : \n %s\n", outputFileName.Data() );
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							  TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}
