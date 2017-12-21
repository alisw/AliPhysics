/***************************************************************************
              fbellini@cern.ch - last modified on 28/11/2013
//
//Lauches KStar analysis with rsn mini package
//Allows basic configuration of pile-up check and event cuts
//
****************************************************************************/
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

AliRsnMiniAnalysisTask * AddTaskF0
(
 TString     outNameSuffix = "f0",                         //suffix for output container
 Bool_t      isMC          = 0,                            //MC flag
 AliPIDResponse::EBeamType collSys = AliPIDResponse::kPP, //=0, kPPB=1, kPBPB=2 outNameSuffix
 UInt_t      triggerMask   = AliVEvent::kINT7,             //trigger selection
 Bool_t      enaMultSel    = kTRUE,                        //enable multiplicity axis
 Int_t       evtCutSetID   = eventCutSet::kEvtDefault,     //event selection
 Int_t       pairCutSetID  = pairYCutSet::kPairDefault,    //pair cuts
 Int_t       aodFilterBit  = 5,                            //AOD filter bit for AOD only
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPiPid = AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s, //PID cut
 Float_t     nsigma        = 3.0,   //PID cut
 Float_t     masslow       = 0.3,   //inv mass range lower boundary
 Float_t     massup        = 1.3,   //inv mass range upper boundary
 Int_t       nbins         = 1000,  //inv mass: N bins
 Bool_t      enableMonitor = kTRUE) //enable single track QA
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

  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  //
  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskF0", "No analysis manager to connect to.");
    return NULL;
  } 

  // create the task and configure 
  TString taskName = Form("RsnTaskF0");
  AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);

  //trigger 
  task->UseESDTriggerMask(triggerMask); //ESD
  //task->SelectCollisionCandidates(triggerMask); //AOD
  
  //-----------------------------------------------------------------------------------------------
  // -- MULTIPLICITY/CENTRALITY -------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  if (collSys==AliPIDResponse::kPP) task->UseMultiplicity("AliMultSelection_V0M");
  if (collSys==AliPIDResponse::kPPB) task->UseCentrality("AliMultSelection_V0A");
  if (collSys==AliPIDResponse::kPBPB) task->UseCentrality("AliMultSelection_V0M");
  
  //-----------------------------------------------------------------------------------------------
  // -- EVENT MIXING CONFIG -----------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  task->UseContinuousMix();
  task->SetNMix(nmix);
  task->SetMaxDiffVz(maxDiffVzMix);
  task->SetMaxDiffMult(maxDiffMultMix);
  ::Info("AddTaskF0", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f", nmix, maxDiffVzMix, maxDiffMultMix));
   
  //-----------------------------------------------------------------------------------------------
  // -- EVENT SELECTION ---------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
   if ((collSys==AliPIDResponse::kPP) && (!isMC)) cutVertex->SetCheckPileUp(rejectPileUp);   // set the check for pileup
   cutVertex->SetCheckZResolutionSPD(); ::Info("AddTaskF0", Form("CheckZResolutionSPD:              ON"));
  cutVertex->SetCheckDispersionSPD(); ::Info("AddTaskF0", Form("CheckDispersionSPD:               ON"));
  cutVertex->SetCheckZDifferenceSPDTrack(); ::Info("AddTaskF0", Form("CheckZDifferenceSPDTrack:         ON"));
  
  //set check for pileup in 2013
  AliRsnCutEventUtils *cutEventUtils = new AliRsnCutEventUtils("cutEventUtils", kFALSE, rejectPileUp);
  cutEventUtils->SetCheckIncompleteDAQ(kTRUE);
  cutEventUtils->SetCheckSPDClusterVsTrackletBG();
  ::Info("AddTaskF0", Form("CheckIncompleteDAQ:                  ON"));
  ::Info("AddTaskF0", Form("SetCheckSPDClusterVsTrackletBG:      ON"));
  
  if (useMVPileUpSelection){
    cutEventUtils->SetUseMVPlpSelection(useMVPileUpSelection);
    cutEventUtils->SetMinPlpContribMV(MinPlpContribMV);
    cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
    ::Info("AddTaskF0", Form("Multiple-vtx Pile-up rejection:      ON \nSettings: MinPlpContribMV = %i, MinPlpContribSPD = %i", MinPlpContribMV, MinPlpContribSPD));
  } else {
    cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
    ::Info("AddTaskF0", Form("SPD Pile-up rejection:                      ON \nSettings: MinPlpContribSPD = %i", MinPlpContribSPD));
  }
  ::Info("AddTaskF0", Form("Pile-up rejection mode:                %s", (rejectPileUp?"ON":"OFF")));   
  
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
  
  //multiplicity or centrality monitoring -- if enabled
  if (enaMultSel){

    //multiplicity or centrality with forward estimator from AliMulSelectionTask 
    Int_t multID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
    //reference multiplicity (default with global tracks with good quality, if not available uses tracklets)
    Int_t multRefID = task->CreateValue(AliRsnMiniValue::kRefMult, kFALSE);
    
    AliRsnMiniOutput *outMult = task->CreateOutput("eventMult", "HIST", "EVENT");
    outMult->AddAxis(multID, 100, 0.0, 100.0);
    
    AliRsnMiniOutput *outRefMult = task->CreateOutput("eventRefMult", "HIST", "EVENT");
    outRefMult->AddAxis(multRefID, 400, 0.0, 400.0);
    
    TH2F* hvz = new TH2F("hVzVsCent",Form("Vertex position vs centrality"), 101, 0., 101., 500, -50.0, 50.0);
    if (collSys==AliPIDResponse::kPPB) 
      hvz->GetXaxis()->SetTitle("V0A");
    else
      hvz->GetXaxis()->SetTitle("V0M");
    hvz->GetYaxis()->SetTitle("z_{vtx} (cm)");
    task->SetEventQAHist("vz", hvz);
    
    TH2F* hRefMultiVsCent = new TH2F("hRefMultiVsCent",Form("Reference multiplicity vs centrality"), 101, 0., 101., 400, 0., 400.);
    if (collSys==AliPIDResponse::kPPB) 
      hRefMultiVsCent->GetXaxis()->SetTitle("V0A");
    else 
      hRefMultiVsCent->GetXaxis()->SetTitle("V0M");
    hRefMultiVsCent->GetYaxis()->SetTitle("GLOBAL");
    task->SetEventQAHist("refmulti",hRefMultiVsCent);
  
    TH2F* hMultiVsCent = new TH2F("hMultiVsCent",Form("Multiplicity vs centrality"), 101, 0., 101., 400, 0., 400.);
    if (collSys==AliPIDResponse::kPPB) 
      hMultiVsCent->GetXaxis()->SetTitle("V0A");
    else 
      hMultiVsCent->GetXaxis()->SetTitle("V0M");
    hMultiVsCent->GetYaxis()->SetTitle("QUALITY");
    task->SetEventQAHist("multicent",hMultiVsCent);
  }

  //-----------------------------------------------------------------------------------------------
  // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  cutY->SetRangeD(minYlab, maxYlab);
  
  AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
  cutsPair->AddCut(cutY);
  cutsPair->SetCutScheme(cutY->GetName());

  //-----------------------------------------------------------------------------------------------
  // -- CONFIG ANALYSIS --------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigF0.C");
  if (!ConfigF0(task, isMC, collSys, cutsPair, enaMultSel, masslow, massup, nbins, aodFilterBit, cutPiPid, nsigma, enableMonitor) ) return 0x0;
  
  //-----------------------------------------------------------------------------------------------
  // -- CONTAINERS --------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  Printf("AddTaskF0 - Set OutputFileName : \n %s\n", outputFileName.Data() );
  AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							  TList::Class(), 
							  AliAnalysisManager::kOutputContainer, 
							  outputFileName);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);
   
  return task;
}
