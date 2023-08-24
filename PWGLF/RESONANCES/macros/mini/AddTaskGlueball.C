/***************************************************************************
// fbellini@cern.ch - last modified on 17/04/2019
//
//Launches f0(980) analysis with rsn mini package for pp, p-Pb and Pb-Pb analysis
//Allows basic configuration of pile-up check and event cuts
//Rsn output is instead configured in the ConfigGlueball.C macro
//
****************************************************************************/
#ifdef __CLING__
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGLF/RESONANCES/macros/mini/ConfigGlueball.C>
#endif

enum pairYCutSet { kPairDefault,
		   kCentralTight,
		   kpADefault,
		   kpACentral
		 };

enum eventCutSet { kEvtDefault = 0,  //data
		   kEvtDefaultVtx8,
		   kEvtDefaultVtx5,
		   kNoPileUpCut,
		   kMCEvt,         //MC events cuts
		   kMCEvtDefault};

AliRsnMiniAnalysisTask * AddTaskGlueball
(
 TString     outNameSuffix = "glueball",                         //suffix for output container
 Bool_t      isMC          = 0,                            //MC flag
 Bool_t      isAOD          = 1,                            //AOD = 1, ESD = 0
 AliPIDResponse::EBeamType collSys = AliPIDResponse::kPP, //=0, kPPB=1, kPBPB=2 outNameSuffix
 UInt_t      triggerMask   = AliVEvent::kINT7,             //trigger selection
 Bool_t      enaMultSel    = kFALSE,                        //enable multiplicity axis
 Int_t       evtCutSetID   = eventCutSet::kEvtDefault,     //event selection
 Int_t       pairCutSetID  = pairYCutSet::kPairDefault,    //pair cuts
 Int_t       aodFilterBit  = 32,                            //AOD filter bit for AOD only
 Int_t       customQualityCutsID = 1,                      //track quality cuts
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPiPid = AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s, //PID cut
 Float_t     nsigmaTPC        = 2.0,   //PID cut of TPC
 Float_t     nsigmaTOF        = 2.0,   //PID cut of TOF
 Float_t     masslow       = 0.9,   //inv mass range lower boundary
 Float_t     massup        = 2.5,   //inv mass range upper boundary
 Int_t       nbins         = 1600,  //inv mass: N bins
 Float_t     ptlow         = 0.0,   //pT axis low edge
 Float_t     ptup          = 10.0,   //pT axis upper edge
 Int_t       nbinspt       = 20,   //pT axis n bins
 Float_t     ptMClowCut    = 0.001, // lower pT cut on MC
 Float_t     ptMCupCut     = 20.0,  // upper pT cut on MC
 Bool_t      enableTrackQA = kTRUE, //enable single track QA
 Bool_t      enableAdvEvtQA = kFALSE) //enable advanced QA for multiplicity and event properties
{

  //-------------------------------------------
  // event cuts
  //-------------------------------------------
  /*Trigger selection
  Std minimum bias trigger AliVEvent::kMB, corresponding to (V0A | V0C | SPD) to be used with
  - pp 7 TeV (2010 data)
  - PbPb 2.76 TeV (2010 data and 2011 data)
  - pp 2.76 TeV (2011 data)
  Centrality triggers AliVEvent::kCentral, AliVEvent::kSemicentral to be used with
  - PbPb 2.76 TeV (2011 data)
  Std minimum bias trigger AliVEvent::kINT7, corrsesponding to (V0A & V0C) to be used with
  - pA 5.02 TeV (2013 data)
  - pp 13 TeV (2015 data)
  */
  Double_t vtxZcut            = 10.0; //cm, default cut on vtx z
  Bool_t rejectPileUp         = kTRUE; //rejects pileup from SPD
  Bool_t useMVPileUpSelection = kFALSE; //alternative pile-up rejection, default is SPD
  Int_t  MinPlpContribSPD     = 5; //default value if used
  Int_t  MinPlpContribMV      = 5; //default value if used
  // Bool_t selectDPMJETevtNSDpA = kFALSE; //cut to select true NSD events in DPMJET

  if (evtCutSetID==eventCutSet::kEvtDefaultVtx8) vtxZcut = 8.0; //cm
  if (evtCutSetID==eventCutSet::kEvtDefaultVtx8) vtxZcut = 5.0; //cm
  if (evtCutSetID==eventCutSet::kNoPileUpCut) rejectPileUp = kFALSE;

  if (evtCutSetID>=eventCutSet::kMCEvt) {
    vtxZcut = 1.0e3; //cm
  }

  if (evtCutSetID==eventCutSet::kMCEvtDefault) {
    vtxZcut = 10.0; //cm
  }

  //-------------------------------------------
  //pair rapidity cut
  //-------------------------------------------
  Double_t minYlab = -0.5;
  Double_t maxYlab =  0.5;

  if (pairCutSetID==pairYCutSet::kCentralTight) { //|y_cm|<0.3
    minYlab = -0.3;    maxYlab = 0.3;
  }

  if (pairCutSetID==pairYCutSet::kpADefault) { //-0.5<y_cm<0.0
    minYlab = -0.465;    maxYlab = 0.035;
  }

  if (pairCutSetID==pairYCutSet::kpACentral) { //|y_cm|<0.3
    minYlab = -0.765;    maxYlab = -0.165;
  }

  //-------------------------------------------
  //mixing settings
  //-------------------------------------------
  Int_t       nmix = 5;
  Float_t     maxDiffVzMix = 1.0;
  Float_t     maxDiffMultMix = 5.0;
  if (collSys==AliPIDResponse::kPBPB) maxDiffMultMix = 10.0;
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  //
  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskGlueball", "No analysis manager to connect to.");
    return NULL;
  }

  // create the task and configure
  TString taskName = Form("RsnTaskGlueball");
  AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);

  //trigger
  if (!isAOD) task->UseESDTriggerMask(triggerMask); //ESD
  else task->SelectCollisionCandidates(triggerMask); //AOD

  //-----------------------------------------------------------------------------------------------
  // -- MULTIPLICITY/CENTRALITY -------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  if (collSys==AliPIDResponse::kPP) task->UseMultiplicity("AliMultSelection_V0M");
  if (collSys==AliPIDResponse::kPPB) task->UseMultiplicity("AliMultSelection_V0A");
  if (collSys==AliPIDResponse::kPBPB) task->UseMultiplicity("AliMultSelection_V0M");

  //-----------------------------------------------------------------------------------------------
  // -- EVENT MIXING CONFIG -----------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  task->UseContinuousMix();
  task->SetNMix(nmix);
  task->SetMaxDiffVz(maxDiffVzMix);
  task->SetMaxDiffMult(maxDiffMultMix);
  //::Info("AddTaskGlueball", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f", nmix, maxDiffVzMix, maxDiffMultMix));

  //-----------------------------------------------------------------------------------------------
  // -- EVENT SELECTION ---------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
  if ((collSys==AliPIDResponse::kPP) && (!isMC)) cutVertex->SetCheckPileUp(rejectPileUp);   // set the check for pileup
  cutVertex->SetCheckZResolutionSPD();
  Printf("AddTaskGlueball - CheckZResolutionSPD:              ON");
  cutVertex->SetCheckDispersionSPD();
  Printf("AddTaskGlueball - CheckDispersionSPD:               ON");
  cutVertex->SetCheckZDifferenceSPDTrack();
  Printf("AddTaskGlueball - CheckZDifferenceSPDTrack:         ON");

  //set check for pileup in 2013
  AliRsnCutEventUtils *cutEventUtils = new AliRsnCutEventUtils("cutEventUtils", kFALSE, rejectPileUp);
  cutEventUtils->SetCheckIncompleteDAQ(kTRUE);
  cutEventUtils->SetCheckSPDClusterVsTrackletBG();
  Printf("AddTaskGlueball - CheckIncompleteDAQ:                  ON");
  Printf("AddTaskGlueball - SetCheckSPDClusterVsTrackletBG:      ON");

  if (useMVPileUpSelection){
    cutEventUtils->SetUseMVPlpSelection(useMVPileUpSelection);
    cutEventUtils->SetMinPlpContribMV(MinPlpContribMV);
    cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
    //::Info("AddTaskGlueball", Form("Multiple-vtx Pile-up rejection:      ON \nSettings: MinPlpContribMV = %i, MinPlpContribSPD = %i", MinPlpContribMV, MinPlpContribSPD));
  } else {
    cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
    Printf("AddTaskGlueball - SPD Pile-up rejection:     ON \nSettings: MinPlpContribSPD = %i", MinPlpContribSPD);
  }
  Printf("AddTaskGlueball - Pile-up rejection mode:     %i", rejectPileUp);

  AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
  eventCuts->AddCut(cutEventUtils);
  eventCuts->AddCut(cutVertex);
  eventCuts->SetCutScheme(Form("%s&%s", cutEventUtils->GetName(), cutVertex->GetName()));
  task->SetEventCuts(eventCuts);

  //connect task
  mgr->AddTask(task);

  //-----------------------------------------------------------------------------------------------
  // -- EVENT SELECTION MONITOR -------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  //vertex position monitoring
  Int_t vtxID = task->CreateValue(AliRsnMiniValue::kVz, kFALSE);
  AliRsnMiniOutput *outVtx = task->CreateOutput("eventVtx", "HIST", "EVENT");
  outVtx->AddAxis(vtxID, 500, -50.0, 50.0);

  //multiplicity or centrality monitoring
  //multiplicity or centrality with forward estimator from AliMulSelectionTask
  Int_t multID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  //reference multiplicity (default with global tracks with good quality, if not available uses tracklets)
  Int_t multRefID = task->CreateValue(AliRsnMiniValue::kRefMult, kFALSE);

  AliRsnMiniOutput *outMult = task->CreateOutput("eventMult", "HIST", "EVENT");
  outMult->AddAxis(multID, 100, 0.0, 100.0);

  AliRsnMiniOutput *outRefMult = task->CreateOutput("eventRefMult", "HIST", "EVENT");
  outRefMult->AddAxis(multRefID, 400, 0.0, 400.0);

  if (enaMultSel) {
    TH2F* hvz = new TH2F("hVzVsCent","Vertex position vs centrality; multiplicity (%); z_{vtx} (cm); Counts", 101, 0., 101., 240, -12.0, 12.0);
    if (collSys==AliPIDResponse::kPPB)
      hvz->GetXaxis()->SetTitle("V0A");
    else
      hvz->GetXaxis()->SetTitle("V0M");
    task->SetEventQAHist("vz", hvz);

    TH2F* hRefMultiVsCent = new TH2F("hRefMultiVsCent","Reference multiplicity vs centrality; multiplicity (%); GLOBAL; Counts", 101, 0., 101., 400, 0., 400.);
    if (collSys==AliPIDResponse::kPPB)
      hRefMultiVsCent->GetXaxis()->SetTitle("V0A");
    else
      hRefMultiVsCent->GetXaxis()->SetTitle("V0M");
    if (enableAdvEvtQA) task->SetEventQAHist("refmulti",hRefMultiVsCent);

    TH2F* hMultiVsCent = new TH2F("hMultiVsCent","Multiplicity vs centrality; multiplicity (%); QUALITY (%); Counts", 101, 0., 101., 400, 0., 400.);
    if (collSys==AliPIDResponse::kPPB)
      hMultiVsCent->GetXaxis()->SetTitle("V0A");
    else
      hMultiVsCent->GetXaxis()->SetTitle("V0M");
    if (enableAdvEvtQA) task->SetEventQAHist("multicent",hMultiVsCent);
  }

  //-----------------------------------------------------------------------------------------------
  // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  cutY->SetRangeD(minYlab, maxYlab);
  AliRsnCutMiniPair *cutPtMC = new AliRsnCutMiniPair("cutPtMC", AliRsnCutMiniPair::kPtMC);
  cutPtMC->SetRangeD(ptMClowCut,ptMCupCut);

  AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
  cutsPair->AddCut(cutY);
  if (isMC) cutsPair->AddCut(cutPtMC);
  cutsPair->SetCutScheme(Form("%s&%s",cutY->GetName(),cutPtMC->GetName()));

  //-----------------------------------------------------------------------------------------------
  // -- CONFIG ANALYSIS --------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
#ifdef __CINT__
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigGlueball.C");
  //gROOT->LoadMacro("./ConfigGlueball.C");
#endif
  if (!ConfigGlueball(task, isMC, collSys, cutsPair, enaMultSel, masslow, massup, nbins, ptlow, ptup, nbinspt, aodFilterBit,customQualityCutsID, cutPiPid, nsigmaTPC, nsigmaTOF, enableTrackQA) ) return 0x0;

  //-----------------------------------------------------------------------------------------------
  // -- CONTAINERS --------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  Printf("AddTaskGlueball - Set OutputFileName : \n %s\n", outputFileName.Data() );
  AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()),
							  TList::Class(),
							  AliAnalysisManager::kOutputContainer,
							  outputFileName);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);

  return task;
}
