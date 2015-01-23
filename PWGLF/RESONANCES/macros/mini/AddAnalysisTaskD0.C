/***************************************************************************
              fbellini@cern.ch - last modified on 06/08/2012
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

AliRsnMiniAnalysisTask * AddAnalysisTaskD0
(
   Bool_t      isMC,
   Bool_t      isPP,
   Bool_t      ispPb,
   Bool_t      monitor = kTRUE,
   Bool_t      centortracklets = kTRUE,
   Bool_t      sanityhistos = kTRUE,
   TString     centrality = "V0M",
   Int_t       aodFilterBit = 5,  
   Float_t     nsigmaTPCPi = 3.0,
   Float_t     nsigmaTPCKa = 3.0,
   Float_t     nsigmaTOFPi = 2.0,
   Float_t     nsigmaTOFKa = 2.0,
   Float_t     trackDCAcutMax = 7.0,   
   Float_t     trackDCAZcutMax = 2.0,
   Int_t       NTPCcluster = 70,
   Double_t    NTPCcrratio = 0.8,
   Int_t       minSPDclt = 0,
   Double_t    minpt = 0.15,
   TString     triggerMask = AliVEvent::kMB,
   Bool_t      useNTPCclt = kTRUE,   
   Bool_t      maxDCAcutFixed = kFALSE,
   Bool_t      ptdepPIDcut = kFALSE,
   Bool_t      fixedYcut = kTRUE,
   Bool_t      checkpileup = kFALSE,
   Bool_t      SPDpileup = kTRUE,
   Bool_t      doCalculationInMC = kTRUE,
   UShort_t    originDselection = 0,
   Int_t       nmix = 5,
   Double_t    minYlab = -0.5,
   Double_t    maxYlab = 0.5,
   Float_t     mineta = -0.8,
   Float_t     maxeta = 0.8,
   Float_t     min_inv_mass = 0.6,
   Float_t     max_inv_mass = 2.2,
   Int_t       bins = 320,
   Float_t     maxDiffVzMix = 1.0,
   Float_t     maxDiffMultMix = 10.0,
   Float_t     maxDiffAngleMixDeg = 20.0,
   TString     outNameSuffix = "D0"  
)
{  
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //
  Float_t     trackDCAcutMin = 0.0;
  Bool_t      minDCAcutFixed = kTRUE;
  Double_t    dcaProduct = 100.0;
  Float_t     cutV = 10.0;
  Short_t     maxSisters = 2;
  Bool_t      checkP = kTRUE;
  Bool_t      checkFeedDown = kTRUE;
  Bool_t      checkQuark = kTRUE;
  Int_t       aodN = 0;
 
  
  UInt_t trigger = 0;
  
  
  
  TString eventType = "";
   if(triggerMask=="AliVEvent::kMB | AliVEvent::kCentral") {trigger = AliVEvent::kMB | AliVEvent::kCentral; eventType+="Central";}
   else if(triggerMask=="AliVEvent::kMB | AliVEvent::kSemiCentral") {trigger = AliVEvent::kMB | AliVEvent::kSemiCentral; eventType+="SemiCentral";}
   else if(triggerMask=="AliVEvent::kINT7") {trigger = AliVEvent::kINT7; eventType+="kINT7";}
   else if(triggerMask=="AliVEvent::kMB") {trigger = AliVEvent::kMB; eventType+="MinimumBias";}

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddAnalysisTaskD0", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   TString taskName = Form("D0%s%s_%.1f_%d_%.2f_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.1f_%.2f_%d_%s", (isPP? "pp" : ispPb? "pPB": "PbPb"), (isMC ? "MC" : "Data"), cutV, NTPCcluster, NTPCcrratio, minSPDclt, nsigmaTPCPi, nsigmaTPCKa, nsigmaTOFPi, nsigmaTOFKa, trackDCAcutMax, trackDCAZcutMax, minpt, originDselection, eventType.Data());
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   if (!isMC && !isPP){
     Printf(Form("========== SETTING USE CENTRALITY PATCH AOD049 : %s", (aodN==49)? "yes" : "no"));
     task->SetUseCentralityPatch(aodN==49);
   }
   
   ::Info("AddAnalysisTaskD0", Form("TriggerMask: %i",trigger));

   task->SelectCollisionCandidates(trigger);
   task->SetMaxNDaughters(maxSisters);
   task->SetCheckMomentumConservation(checkP);
   task->SetCheckFeedDown(checkFeedDown);
   task->SetRejectCandidateIfNotFromQuark(checkQuark);
   task->SetDselection(originDselection);
   task->KeepMotherInAcceptance(kTRUE);
   task->SetMotherAcceptanceCutMinPt(minpt);
   task->SetMotherAcceptanceCutMaxEta(maxeta);
   
      
   ::Info("AddAnalysisTaskD0", Form("Maximum numbers of daughters allowed (-1 means cut not applied): %i",maxSisters));
   ::Info("AddAnalysisTaskD0", Form("Are we checking the momentum conservation? %s", checkP? "yes" : "no"));
   ::Info("AddAnalysisTaskD0", Form("Are we checking the feeddown? %s", checkFeedDown? "yes" : "no"));
   ::Info("AddAnalysisTaskD0", Form("Are we rejecting the Hijing generated? %s", checkQuark? "yes" : "no"));
   ::Info("AddAnalysisTaskD0", Form("Which D0 are we keeping? %s", (originDselection==0? "only from c quark" : originDselection==1? "only from b quark" : "both from c and b quark") ));
   ::Info("AddAnalysisTaskD0", Form("Selecting Mother in Acceptance: Min pT %.1f, Eta Range %.1f - %.1f", minpt, mineta, maxeta));


   if (isPP) 
     task->UseMultiplicity("QUALITY");
   else
     task->UseCentrality(centrality);   
   // set event mixing options
   task->UseContinuousMix();
   //task->UseBinnedMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   if (!isPP && !ispPb) task->SetMaxDiffAngle(maxDiffAngleMixDeg*TMath::DegToRad()); //set angle diff in rad
   ::Info("AddAnalysisTaskD0", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n max diff EP angle = %5.3f deg", nmix, maxDiffVzMix, maxDiffMultMix, (isPP ? 0.0 : ispPb ? 0.0 : maxDiffAngleMixDeg)));
   
   mgr->AddTask(task);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //  
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", cutV, 1, kFALSE, kFALSE);
   AliRsnCutEventUtils *cutEventUtils = new AliRsnCutEventUtils("cutEventUtils", kFALSE, kTRUE);
   
   if(checkpileup == kTRUE){
   	
	if(SPDpileup == kTRUE) {cutEventUtils->SetRemovePileUppA2013(kTRUE); 
				//cutEventUtils->SetRemoveFirstEvtInChunk(kTRUE);
				//cutEventUtils->SetUseVertexSelection2013pA(kTRUE);
				cutEventUtils->SetUseMVPlpSelection(kFALSE); 
				cutEventUtils->SetMinPlpContribSPD(5);
									  }
   	if(SPDpileup == kFALSE){cutEventUtils->SetRemovePileUppA2013(kTRUE); 
				//cutEventUtils->SetRemoveFirstEvtInChunk(kTRUE);
				//cutEventUtils->SetUseVertexSelection2013pA(kTRUE);
				cutEventUtils->SetUseMVPlpSelection(kTRUE); 
				cutEventUtils->SetMinPlpContribMV(5);
									  }
   } 
		 
   ::Info("AddAnalysisTaskD0", Form("Checking Pile up? %s", (checkpileup ? "yes" : "no") ));
   ::Info("AddAnalysisTaskD0", Form("Which Method? %s", (checkpileup? (SPDpileup ? "SPD" : "Multi Vertex") : "none" )));
   
   // define and fill cut set for event cut
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutVertex);
   //if(SPDpileup == kFALSE) eventCuts->AddCut(cutEventUtils);
   eventCuts->AddCut(cutEventUtils);
   eventCuts->SetCutScheme(Form("%s&%s", cutEventUtils->GetName(), cutVertex->GetName() ));
   eventCuts->ShowCuts();
   eventCuts->PrintSetInfo();
   // set cuts in task
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
   AliRsnMiniOutput *outMult = task->CreateOutput("eventMult", "HIST", "EVENT");
   if (isPP) 
     outMult->AddAxis(multID, 400, 0.0, 400.0);
   else
     outMult->AddAxis(multID, 100, 0.0, 100.0);
     
   //tracklets
   Int_t trackletID = task->CreateValue(AliRsnMiniValue::kTracklets, kFALSE);
   AliRsnMiniOutput *outTracklets = task->CreateOutput("eventTracklets", "HIST", "EVENT");
   outTracklets->AddAxis(trackletID, 250, -0.5, 249.5);
   
   
   //event plane (only for PbPb)
   Int_t planeID = task->CreateValue(AliRsnMiniValue::kPlaneAngle, kFALSE);
   AliRsnMiniOutput *outPlane = 0x0;
   if (!isPP){
     outPlane = task->CreateOutput("eventPlane", "HIST", "EVENT");
     outPlane->AddAxis(planeID, 180, 0.0, TMath::Pi());
     }
     
   TH2F* hvz=new TH2F("hVzVsCent","",100,0.,100., 220,-11.,11.);
   task->SetEventQAHist("vz",hvz);//plugs this histogram into the fHAEventVz data member

   TH2F* hmc=new TH2F("MultiVsCent","",100,0.,100.,1000,0.,1000.);
   hmc->GetYaxis()->SetTitle("QUALITY");
   task->SetEventQAHist("multicent",hmc);//plugs this histogram into the fHAEventMultiCent data member

   TH2F* hep=new TH2F("hEventPlaneVsCent","",100,0.,100., 180,0.,TMath::Pi());
   task->SetEventQAHist("eventplane",hep);//plugs this histogram into the fHAEventPlane data member
   
   //
   // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
   //
   AliRsnCutMiniPair *cutY = 0x0;
   
   if(fixedYcut == kTRUE){
   cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
   cutY->SetRangeD(minYlab, maxYlab);
   }
   else {
   cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityFiducialRegion);
   cutY->SetRangeD(-0.8, 0.8);
   cutY->SetPtDepCut(kTRUE);
   cutY->SetMinPt(0.0);
   cutY->SetMaxPt(5.0);
   cutY->SetPtDepCutMaxFormula("-0.2/15*pt*pt+1.9/15*pt+0.5");
   cutY->SetPtDepCutMinFormula("0.2/15*pt*pt-1.9/15*pt-0.5");
   }
   
   ::Info("AddAnalysisTaskD0", Form("Fixed Y cut? %s", (fixedYcut ? "yes" : "no") ));
   
   AliRsnCutMiniPair *cutDCAproduct = new AliRsnCutMiniPair("cutDCAproduct", AliRsnCutMiniPair::kDCAproduct);
   cutDCAproduct->SetRangeD(-1E20, dcaProduct);
   
   AliRsnCutSet *cutsPairY = new AliRsnCutSet("pairCutsY", AliRsnTarget::kMother);
   cutsPairY->AddCut(cutY);
   cutsPairY->UseMonitor(kTRUE);
   cutsPairY->SetCutScheme("setPairD0_Y");
   cutsPairY->ShowCuts();
   cutsPairY->PrintSetInfo();
   
   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->AddCut(cutDCAproduct);
   cutsPair->UseMonitor(kTRUE);
   cutsPair->SetCutScheme(Form("%s&%s", cutY->GetName(), cutDCAproduct->GetName()));
   cutsPair->ShowCuts();
   cutsPair->PrintSetInfo();
   
   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigD0.C");

   if (isMC) {
       Printf("========================== MC analysis - PID cuts used");
   } else 
     Printf("========================== DATA analysis - PID cuts used");
   if (!ConfigD0(task, isPP, isMC, monitor, centortracklets, sanityhistos, nsigmaTPCPi, nsigmaTPCKa, nsigmaTOFPi, nsigmaTOFKa, aodFilterBit, trackDCAcutMax, trackDCAcutMin, trackDCAZcutMax, NTPCcluster, NTPCcrratio, minSPDclt, minpt, maxSisters, checkP, useNTPCclt, minDCAcutFixed, maxDCAcutFixed, ptdepPIDcut, checkFeedDown, checkQuark, doCalculationInMC, originDselection, mineta, maxeta, min_inv_mass, max_inv_mass, bins, "", cutsPairY, cutsPair)) return 0x0;
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   Printf("AddAnalysisTaskD0 - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("%s_%.1f_%d_%.2f_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.1f_%.2f_%d_%s",outNameSuffix.Data(),cutV,NTPCcluster,NTPCcrratio,minSPDclt,nsigmaTPCPi,nsigmaTPCKa,nsigmaTOFPi,nsigmaTOFKa,trackDCAcutMax,trackDCAZcutMax,minpt,originDselection,eventType.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}
