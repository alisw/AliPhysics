/***************************************************************************
// Created on 22.03.2017 by
// Francesca Bellini f(bellini@cern.ch) - Sourav Kundu (sourav.kundu@cern.ch)
//
//Lauches Phi analysis with rsn mini package for QA
//Allows basic configuration of pile-up check and event cuts
//
****************************************************************************/
enum pairYCutSet { kYcentral = 0,
		   kYpPb5TeV,
		   // kYpPb8TeV,
		   // kYPbp5TeV,
		   // kYPbp8TeV,
		   kYcentralTight
                 };

enum eventCutSet { kEvtDefault = 0,
		   kMCEvtDefault,
		   kPileUpCut
                 };

AliRsnMiniAnalysisTask * AddTaskRsnQA(
 Bool_t      isMC = kFALSE,
 Bool_t      useGeoCutsPbPb2015 = kFALSE,
 TString     multEstimator = "AliMultSelection_V0M",
 UInt_t      triggerMask = AliVEvent::kINT7,
 TString     outNameSuffix = "phi",
 Int_t       evtCutSetID = 0,  //0 for data, 1 for MC, 2 for data with pile-up rejection
 Int_t       pairCutSetID = 0, //selects on pair rapidity: 0 for symmetric system, 1 for p-Pb 5 TeV, 2 for p-Pb 8 TeV
 Int_t       aodFilterBit = 5, //filter bit 5 corresponds to StdITSTPCtrackCuts2011 with TPC crossed rows
 Bool_t      enableMonitor = kTRUE,
 TString     monitorOpt = "NoSIGN")
{
  //-------------------------------------------
  // event cuts
  //-------------------------------------------
  if (evtCutSetID == eventCutSet::kPileUpCut)
    rejectPileUp = kTRUE;
  else
    rejectPileUp = kFALSE;

  //-------------------------------------------
  //pair cuts
  //-------------------------------------------
  Double_t    minYlab =  -0.5;
  Double_t    maxYlab =  0.5;

  if (pairCutSetID==pairYCutSet::kYpPb5TeV) {
    //-0.5 < y_cm < 0.0
    minYlab = -0.465;    maxYlab = 0.035;
  }

  // if (pairCutSetID==pairYCutSet::kYpPb8TeV) {
  //   //to be calculated
  //   minYlab = -0.9;    maxYlab = 0.9;
  // }

  // if (pairCutSetID==pairYCutSet::kYPbp5TeV) {
  //   //to be calculated
  //   minYlab = -0.9;    maxYlab = 0.9;
  // }

  // if (pairCutSetID==pairYCutSet::kYPbp5TeV) {
  //   //to be calculated
  //   minYlab = -0.765;    maxYlab = -0.165;
  // }

  if (pairCutSetID==pairYCutSet::kYcentralTight) {
    //|y_cm| < 0.3
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
      ::Error("AddTaskRsnQA", "No analysis manager to connect to.");
      return NULL;
   }

   // create the task and configure
   TString taskName = Form("taskRsnQA");
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), kTRUE);

   //Set trigger selection
   task->UseESDTriggerMask(triggerMask);
   //task->SelectCollisionCandidate(triggerMask);

   /*
   //Set multiplicity/centrality estimator
   //if (isPP)
   task->UseMultiplicity("AliMultSelection_V0M");
     //else
     //task->UseCentrality("V0M");
     */

   if (multEstimator.IsNull()){
     ::Error("AddTaskRsnQA", "No multiplicity selection estimator set.");
     return NULL;
   }

   if (multEstimator.Contains("AliMultSelection"))
     	task->UseMultiplicity(multEstimator.Data());
   else
     	task->UseCentrality(multEstimator.Data());


   //Set event mixing options
   task->UseContinuousMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   ::Info("AddTaskRsnQA", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f", nmix, maxDiffVzMix, maxDiffMultMix));

   //Add task
   mgr->AddTask(task);

   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //
   //  AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE); //cut on z_vtx < 10 cm //
   AliRsnCutEventUtils* cutEventUtils=new AliRsnCutEventUtils("cutEventUtils", kTRUE, rejectPileUp);
   cutEventUtils->SetCheckAcceptedMultSelection();
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutEventUtils);
   //  eventCuts->AddCut(cutVertex);
   eventCuts->SetCutScheme(Form("%s", cutEventUtils->GetName()));
   task->SetEventCuts(eventCuts);

   //
   // -- EVENT-ONLY COMPUTATIONS -------------------------------------------------------------------
   //
   //vertex
   Int_t vtxID = task->CreateValue(AliRsnMiniValue::kVz, kFALSE);
   //multiplicity or centrality
   Int_t multID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);

   AliRsnMiniOutput *outVtx = task->CreateOutput("eventVtx", "HIST", "EVENT");
   outVtx->AddAxis(vtxID, 220, -11.0, 11.0);

   AliRsnMiniOutput *outMult = task->CreateOutput("eventMult", "HIST", "EVENT");
   outMult->AddAxis(multID, 101, 0.0, 101.0); //also in pp, p-Pb the percentile is returned

   TH2F* hvz = new TH2F("hVzVsCent",Form("Vertex position vs centrality"), 101, 0., 101., 220, -11.0, 11.0);
   hvz->GetXaxis()->SetTitle("multiplicity %");
   hvz->GetYaxis()->SetTitle("z_{vtx} (cm)");
   task->SetEventQAHist("vz", hvz);

   // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
   AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
   cutY->SetRangeD(minYlab, maxYlab);

   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());

   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------
   //
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/ConfigRsnQA.C");
   if (!ConfigRsnQA(task, isMC, cutsPair, aodFilterBit, enableMonitor, "NoSIGN", useGeoCutsPbPb2015) ) return 0x0;

   // -- CONTAINERS --------------------------------------------------------------------------------
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   Printf("AddTaskRsnQA - Set OutputFileName : \n %s\n", outputFileName.Data() );
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnQA_%s", outNameSuffix.Data()),
							   TList::Class(),
							   AliAnalysisManager::kOutputContainer,
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);

   return task;
}
