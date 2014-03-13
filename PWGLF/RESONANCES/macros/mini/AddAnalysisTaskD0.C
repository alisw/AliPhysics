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
   Float_t     cutV = 10.0,
   Int_t       aodFilterBit = 5,  
   Float_t     nsigmaTPCPi = 3.0,
   Float_t     nsigmaTPCKa = 3.0,
   Float_t     nsigmaTOFPi = 2.0,
   Float_t     nsigmaTOFKa = 2.0,
   Float_t     trackDCAcutMax = 7.0,
   Float_t     trackDCAcutMin = 0.0,
   Int_t       NTPCcluster = 70,
   Double_t    minpt = 0.15,
   TString     triggerMask = AliVEvent::kMB,
   Short_t     maxSisters = 2,
   Bool_t      checkP = kTRUE,
   Bool_t      minDCAcutFixed = kFALSE,
   Bool_t      maxDCAcutFixed = kFALSE,
   Bool_t      ptdepPIDcut = kFALSE,
   Int_t       nmix = 5,
   Double_t    minYlab =  -0.5,
   Double_t    maxYlab =  0.5,
   Double_t    dcaProduct = -1E-4,
   Float_t     maxDiffVzMix = 1.0,
   Float_t     maxDiffMultMix = 10.0,
   Float_t     maxDiffAngleMixDeg = 20.0,
   Int_t       aodN = 0,
   TString     outNameSuffix = "D0"
)
{  
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //
  
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
   TString taskName = Form("D0%s%s_%.1f_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.4f_%.5f_%.2f_%s", (isPP? "pp" : "PbPb"), (isMC ? "MC" : "Data"), cutV, NTPCcluster, nsigmaTPCPi, nsigmaTPCKa, nsigmaTOFPi, nsigmaTOFKa, trackDCAcutMax, trackDCAcutMin, dcaProduct, minpt, eventType.Data());
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   if (!isMC && !isPP){
     Printf(Form("========== SETTING USE CENTRALITY PATCH AOD049 : %s", (aodN==49)? "yes" : "no"));
     task->SetUseCentralityPatch(aodN==49);
   }
   
   ::Info("AddAnalysisTaskD0", Form("TriggerMask: %i",trigger));

   task->SelectCollisionCandidates(trigger);
   task->SetMaxNDaughters(maxSisters);
   task->SetCheckMomentumConservation(checkP);
   
   ::Info("AddAnalysisTaskD0", Form("Maximum numbers of daughters allowed (-1 means cut not applied): %i",maxSisters));
   ::Info("AddAnalysisTaskD0", Form("Are we checking the momentum conservation?: %s", checkP? "yes" : "no"));


   if (isPP) 
     task->UseMultiplicity("QUALITY");
   else
     task->UseCentrality("V0M");   
   // set event mixing options
   task->UseContinuousMix();
   //task->UseBinnedMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   if (!isPP) task->SetMaxDiffAngle(maxDiffAngleMixDeg*TMath::DegToRad()); //set angle diff in rad
   ::Info("AddAnalysisTaskD0", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n max diff EP angle = %5.3f deg", nmix, maxDiffVzMix, maxDiffMultMix, (isPP ? 0.0 : maxDiffAngleMixDeg)));
   
   mgr->AddTask(task);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //  
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", cutV, 0, kFALSE);
   if (isPP) cutVertex->SetCheckPileUp(kTRUE);   // set the check for pileup
   
   // define and fill cut set for event cut
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutVertex);
   eventCuts->SetCutScheme(cutVertex->GetName());
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
   
   //event plane (only for PbPb)
   Int_t planeID = task->CreateValue(AliRsnMiniValue::kPlaneAngle, kFALSE);
   AliRsnMiniOutput *outPlane = 0x0;
   if (!isPP){
     outPlane = task->CreateOutput("eventPlane", "HIST", "EVENT");
     outPlane->AddAxis(planeID, 180, 0.0, TMath::Pi());
     }
     
   TH2F* hvz=new TH2F("hVzVsCent","",100,0.,100., 220,-11.,11.);
   task->SetEventQAHist("vz",hvz);//plugs this histogram into the fHAEventVz data member

   TH2F* hmc=new TH2F("MultiVsCent","",100,0.,100.,4000,0.,4000.);
   hmc->GetYaxis()->SetTitle("QUALITY");
   task->SetEventQAHist("multicent",hmc);//plugs this histogram into the fHAEventMultiCent data member

   TH2F* hep=new TH2F("hEventPlaneVsCent","",100,0.,100., 180,0.,TMath::Pi());
   task->SetEventQAHist("eventplane",hep);//plugs this histogram into the fHAEventPlane data member
   
   //
   // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
   //
   AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
   cutY->SetRangeD(minYlab, maxYlab);
   
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
   gROOT->LoadMacro("$ALICE_ROOT/PWGLF/RESONANCES/macros/mini/ConfigD0.C");

   if (isMC) {
       Printf("========================== MC analysis - PID cuts used");
   } else 
     Printf("========================== DATA analysis - PID cuts used");
   if (!ConfigD0(task, isPP, isMC, nsigmaTPCPi, nsigmaTPCKa, nsigmaTOFPi, nsigmaTOFKa, aodFilterBit, trackDCAcutMax, trackDCAcutMin, NTPCcluster, minpt, maxSisters, checkP,  minDCAcutFixed, maxDCAcutFixed, ptdepPIDcut,"", cutsPairY, cutsPair)) return 0x0;
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   Printf("AddAnalysisTaskD0 - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("%s_%.1f_%d_%.1f_%.1f_%.1f_%.1f_%.1f_%.4f_%.5f_%.2f_%s",outNameSuffix.Data(),cutV,NTPCcluster,nsigmaTPCPi,nsigmaTPCKa,nsigmaTOFPi,nsigmaTOFKa,trackDCAcutMax,trackDCAcutMin,dcaProduct,minpt,eventType.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}
