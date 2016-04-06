/***************************************************************************
              fbellini@cern.ch - last modified on 28/11/2013
//            inayat.rasool.bhat@cern.ch   last modified on  07/03/2016 
//Lauches Rho analysis with rsn mini package
//Allows basic configuration of pile-up check and event cuts
//
****************************************************************************/
enum pairYCutSet { kPairDefault,
		   kNegative,
		   kCentral
                 };

enum eventCutSet { kOld = -1, 
		   kEvtDefault, //=0
		   kNoPileUpCut, //=1
		   kPileUpMV, //=2
		   kPileUpSPD3, //=3		      
		   kDefaultVtx8, //=4
		   kDefaultVtx5, //=5                    
		   kMCEvtDefault //=6
};

enum eventMixConfig { kDisabled = -1,
		      kMixDefault,//=0 //10 events, Dvz = 1cm, DC = 10
		      k5Evts, //=1 //5 events, Dvz = 1cm, DC = 10
		      k5Cent,  //=2 //10 events, Dvz = 1cm, DC = 5
                    };


AliRsnMiniAnalysisTask * AddTaskRhoPPB
(
 Bool_t      isMC = kTRUE,
 Bool_t      isPP = kFALSE,
 TString     outNameSuffix = "tpc3s",
 Int_t       evtCutSetID = 0,
 Int_t       pairCutSetID = 0,
 Int_t       mixingConfigID = 0,
 Int_t       aodFilterBit = 5,
 Int_t       customQualityCutsID = -1,
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPiCandidate = AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,
 Float_t     nsigmaPi = 3.0,
 Float_t     nsigmaKa = 3.0,
 Bool_t      enableMonitor = kTRUE,
 Bool_t      IsMcTrueOnly = kFALSE,
 TString     monitorOpt = "NoSIGN",
 // Bool_t      useMixLS = 0,
 Bool_t      checkReflex = 1,
 AliRsnMiniValue::EType yaxisvar = AliRsnMiniValue::kPt
)
{  

  
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
  
  if (evtCutSetID==eventCutSet::kOld) {
    triggerMask = AliVEvent::kAnyINT;
    rmFirstEvtChunk = kFALSE;
    rejectPileUp = kFALSE;
    useVtxCut2013pA = kFALSE;
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
  
  if (evtCutSetID==eventCutSet::kMCEvtDefault) {
    rmFirstEvtChunk = kFALSE;
    rejectPileUp = kFALSE;
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
      ::Error("AddAnalysisTaskTPCRho", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   TString taskName = Form("TPCRho%s%s_%i%i", (isPP? "pp" : "PBPB"), (isMC ? "MC" : "Data"), (Int_t)cutPiCandidate,(Int_t)cutKaCandidate );
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   //task->UseESDTriggerMask(triggerMask); //ESD
   task->SelectCollisionCandidates(triggerMask); //AOD
   
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
   ::Info("AddAnalysisTaskTPCRho", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f", nmix, maxDiffVzMix, maxDiffMultMix));
   
   mgr->AddTask(task);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //  
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors to vtx
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", vtxZcut, 0, kFALSE);
   //if (isPP) cutVertex->SetCheckPileUp(kTRUE);   // set the check for pileup 
   //set check for pileup in 2013
   AliRsnCutEventUtils *cutEventUtils = new AliRsnCutEventUtils("cutEventUtils", rmFirstEvtChunk, rejectPileUp);
   cutEventUtils->SetUseVertexSelection2013pA(useVtxCut2013pA);
   ::Info("AddAnalysisTaskTPCRho", Form(":::::::::::::::::: Vertex cut as pA 2013: %s", (useVtxCut2013pA?"ON":"OFF")));   
   if (useMVPileUpSelection){
     cutEventUtils->SetUseMVPlpSelection(useMVPileUpSelection);
     cutEventUtils->SetMinPlpContribMV(MinPlpContribMV);
     cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
     ::Info("AddAnalysisTaskTPCRho", Form("Multiple-vtx Pile-up rejection settings: MinPlpContribMV = %i, MinPlpContribSPD = %i", MinPlpContribMV, MinPlpContribSPD));
   } else {
     cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
     ::Info("AddAnalysisTaskTPCRho", Form("SPD Pile-up rejection settings: MinPlpContribSPD = %i", MinPlpContribSPD));
   }
   ::Info("AddAnalysisTaskTPCRho", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp?"ON":"OFF")));   
   ::Info("AddAnalysisTaskTPCRho", Form("::::::::::::: Remove first event in chunk: %s", (rmFirstEvtChunk?"ON":"OFF")));   
   
   // define and fill cut set for event cut
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutEventUtils);
   eventCuts->AddCut(cutVertex);
   eventCuts->SetCutScheme(Form("%s&%s", cutEventUtils->GetName(), cutVertex->GetName()));
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
   cutY->SetRangeD(minYlab, maxYlab);
   
   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());
   
   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------
   //   
   // gROOT->LoadMacro("ConfigRhoPPb.C");
   gROOT->LoadMacro("${ALICE_PHYSICS}/PWGLF/RESONANCES/macros/mini/ConfigRhoPPb.C");
   if (!ConfigRhoPPb(task, isMC, isPP, "", cutsPair, aodFilterBit, customQualityCutsID, cutPiCandidate, cutKaCandidate, nsigmaPi, nsigmaKa, enableMonitor, isMC&IsMcTrueOnly,  monitorOpt.Data(),/*  useMixLS, */isMC&checkReflex, yaxisvar)) return 0x0;
   
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddAnalysisTaskTPCRho - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}
