/***************************************************************************
              Paraskevi.Ganoti@cern.ch - last modified on //2019
//
//Lauches SigmaStar analysis with rsn mini package
//Allows basic configuration of pile-up check and event cuts
//
****************************************************************************/

AliRsnMiniAnalysisTask * AddTaskSigPM
(
 TString     outNameSuffix = "SigPM",                         //suffix for output container
 Bool_t      isMC          = kFALSE,                            //MC flag
 AliPIDResponse::EBeamType collSys = AliPIDResponse::kPBPB, //=0, kPPB=1, kPBPB=2 outNameSuffix
 UInt_t      triggerMask   = AliVEvent::kCentral,             //trigger selection
 Bool_t      enaMultSel    = kTRUE,                        //enable multiplicity axis
 Float_t     nsigma        = 3.0,   //PID cut
 Float_t     masslow       = 1.2,   //inv mass range lower boundary
 Float_t     massup        = 3.2,   //inv mass range upper boundary
 Int_t       nbins         = 2000,  //inv mass: N bins
 Bool_t      enableMonitor = kTRUE, //enable single track QA
 Float_t     pi_Ls_PIDCut=4.,        //nsigma V0 daughters 
 Float_t     LsDCA = 0.3,             //V0 vtx to PV
 Float_t     LsCosPoinAn = 0.98,      // cos of Pointing Angle
 Float_t     LsDaughDCA=0.8,          // dca V0 daughters
 Float_t     massTol = 0.006,         // mass tolerance 6 MeV
 Float_t     massTolVeto = 0.0043,    // mass tol veto
 Bool_t      Switch = kFALSE,         // switch
 Float_t     pLife = 25,              // life
 Float_t     v0rapidity= 0.5,         // V0 rapidity
 Float_t     radiuslow=5.,            // radius low 
 Bool_t      doCustomDCAcuts=kTRUE,  //custom dca cuts for V0 daughters
 Double_t    dcaProton=0.1,           // proton dca
 Double_t    dcaPion=0.1)             //pion dca

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
  // // Bool_t selectDPMJETevtNSDpA = kFALSE; //cut to select true NSD events in DPMJET

  //-------------------------------------------
  //pair rapidity cut
  //-------------------------------------------
  Double_t minYlab = -0.5;
  Double_t maxYlab =  0.5;
      
  // if (pairCutSetID==pairYCutSet::kCentralTight) { //|y_cm|<0.3
  //   minYlab = -0.3;    maxYlab = 0.3;
  // }

  // if (pairCutSetID==pairYCutSet::kpADefault) { //-0.5<y_cm<0.0
  //   minYlab = -0.465;    maxYlab = 0.035;
  // }
  
  // if (pairCutSetID==pairYCutSet::kpACentral) { //|y_cm|<0.3
  //   minYlab = -0.765;    maxYlab = -0.165;
  // }
  
  //-------------------------------------------
  //mixing settings
  //-------------------------------------------
  Int_t       nmix = 10;
  Float_t     maxDiffVzMix = 1.0;
  Float_t     maxDiffMultMix = 5.0;

  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  //
  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskSigPM", "No analysis manager to connect to.");
    return NULL;
  } 

  // create the task and configure 
  TString taskName = Form("RsnTaskSigPM");
  AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);

  //trigger 
  if (!isMC) task->UseESDTriggerMask(triggerMask); //ESD
  if(isMC) task->UseESDTriggerMask(AliVEvent::kAnyINT);
//  task->SelectCollisionCandidates(triggerMask); //AOD
  
  //-----------------------------------------------------------------------------------------------
  // -- MULTIPLICITY/CENTRALITY -------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  if (collSys==AliPIDResponse::kPP) task->UseMultiplicity("AliMultSelection_V0M");
  if (collSys==AliPIDResponse::kPPB) task->UseCentrality("AliMultSelection_V0A");
  if (collSys==AliPIDResponse::kPBPB)  task->UseMultiplicity("AliMultSelection_V0M"); //task->UseCentrality("AliMultSelection_V0M");
  
  //-----------------------------------------------------------------------------------------------
  // -- EVENT MIXING CONFIG -----------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  task->UseContinuousMix();
  task->SetNMix(nmix);
  task->SetMaxDiffVz(maxDiffVzMix);
  task->SetMaxDiffMult(maxDiffMultMix);
  ::Info("AddTaskSstar", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f", nmix, maxDiffVzMix, maxDiffMultMix));
   
  //-----------------------------------------------------------------------------------------------
  // -- EVENT SELECTION ---------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);

  AliRsnEventCuts * rsnEventCuts = new AliRsnEventCuts("rsnEventCuts");
  rsnEventCuts->ForceSetupPbPb2018();
  rsnEventCuts->SetUseMultSelectionEvtSel();
  eventCuts->AddCut(rsnEventCuts);
  eventCuts->SetCutScheme(Form("%s", rsnEventCuts->GetName()));

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

  //AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0", AliRsnCutMiniPair::kContainsV0Daughter);
  
  AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
  cutsPair->AddCut(cutY);
  //cutsPair->AddCut(cutV0);
  //cutsPair->SetCutScheme(TString::Format("%s&(!%s)",cutY->GetName(),cutV0->GetName()).Data());
  cutsPair->SetCutScheme(cutY->GetName());
  
  // AliRsnCutSet* PairCutsMix=new AliRsnCutSet("PairCutsMix",AliRsnTarget::kMother);
  // PairCutsMix->AddCut(cutY);
  // PairCutsMix->SetCutScheme(cutY->GetName());

  //-----------------------------------------------------------------------------------------------
  // -- CONFIG ANALYSIS --------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
 gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigSigPM.C");
//  gROOT->LoadMacro("ConfigSigPM.C");
 if (!ConfigSigPM(task, isMC, collSys, cutsPair, enaMultSel, masslow, massup, nbins, nsigma, 
enableMonitor, pi_Ls_PIDCut, LsDCA, LsCosPoinAn, LsDaughDCA, massTol, massTolVeto, Switch, pLife, v0rapidity, radiuslow, doCustomDCAcuts, dcaProton, dcaPion)) 
return 0x0;
  
  //-----------------------------------------------------------------------------------------------
  // -- CONTAINERS --------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  Printf("AddTaskSstar - Set OutputFileName : \n %s\n", outputFileName.Data() );
  AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							  TList::Class(), 
							  AliAnalysisManager::kOutputContainer, 
							  outputFileName);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);
   
  return task;
}
