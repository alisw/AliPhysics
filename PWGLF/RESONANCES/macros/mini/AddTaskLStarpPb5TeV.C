/***************************************************************************
   sarita.sahoo@cern.ch - last modified on 25/06/2015
   rbaral@cern.ch - last modified on 28/07/2016
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
// RAMA check the outnamesuffix = "Default_LStar" here

AliRsnMiniAnalysisTask * AddAnalysisTask
(
 Bool_t      isMC ="",
 Bool_t      isPP ="",
 Int_t       aodFilterBit = 5,
 Int_t       customQualityCutsID=1,
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate = AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,
 Float_t     nsigmaPr = 3.0,
 Float_t     nsigmaKa = 3.0,
 Bool_t      enableMonitor = kTRUE, //kFALSE,
 Bool_t      IsMcTrueOnly = kFALSE, //kTRUE,
 UInt_t      triggerMask = AliVEvent::kINT7,
 Int_t       signedPdg = 3124,
 TString                monitorOpt = "NoSIGN",  //Flag for AddMonitorOutput.C e.g."NoSIGN"
 Bool_t                 useCrossedRows = kFALSE,
 const char *yaxisVar = "PtDaughter",  //yaxisVar = "PtDaughter_PDaughter_cent"
 Bool_t      useMixLS = 0,
 Int_t       nmix = 5,
 Float_t     maxDiffVzMix = 1.0,
 Float_t     maxDiffMultMix = 10.0,
 TString     outNameSuffix = "", 
 Int_t       MultBins
)
{

   AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar;
   AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar;
   
   Float_t     nsigmaPr = 3.0;
   Float_t     nsigmaKa = 3.0;
   Bool_t      enableSyst = kTRUE ;//kFALSE  || kTRUE
   Char_t      DCAxyFormula[100] = "0.0105+0.035/pt^1.1";
   Char_t      DCAxyFormula4[100] = "0.0060+0.020/pt^1.1";
   Char_t      DCAxyFormula5[100] = "0.0075+0.025/pt^1.1";
   Double_t    dcazmax = 2;
   Double_t    minNcls = 50;
   Double_t    maxX2cls = 4.0;
   Float_t     TPCGlobalchi2 = 36;
   Double_t    ITSchi2 = 36;

   Double_t    minCrossedRows = 70.0;
   Double_t    maxClsCrossedRows = 0.8;

   UInt_t      triggerMask = AliVEvent::kINT7;
   TString                monitorOpt = "NoSIGN";  //Flag for AddMonitorOutput.C e.g."NoSIGN"
   Bool_t                 useCrossedRows = kTRUE;
   //   AliRsnMiniValue::EType yaxisVar = AliRsnMiniValue::kPt; //
   Bool_t                 useMixLS = 0;
   Float_t     maxDiffVzMix = 1.0;
   Float_t     maxDiffMultMix = 5.0;


 

  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //


   TString outputFileName = Form("%s", AliAnalysisManager::GetCommonFileName());
  // Objects name


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddAnalysisTaskTPCLStarSyst", "No analysis manager to connect to.");
      return NULL;
   } 


   // create the task and configure 
   TString taskName = outNameSuffix.Data();
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);



   task->SelectCollisionCandidates(triggerMask);


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
   //if (!isPP) task->SetMaxDiffAngle(maxDiffAngleMixDeg*TMath::DegToRad()); //set angle diff in rad
   ::Info("AddAnalysisTasLStar", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n ", nmix, maxDiffVzMix, maxDiffMultMix));
   
   mgr->AddTask(task);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //  
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);

   AliRsnCutEventUtils *cutEventUtils = new AliRsnCutEventUtils("cutEventUtils", kFALSE, kTRUE);

   cutEventUtils->SetCheckAcceptedMultSelection();
       
   // define and fill cut set for event cut
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutEventUtils);
   eventCuts->AddCut(cutVertex);  
   //   eventCuts->SetCutScheme(cutEventUtils->GetName()); // RAMA
   eventCuts->SetCutScheme(Form("%s&%s", cutEventUtils->GetName(), cutVertex->GetName()));



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
     outMult->AddAxis(multID, 300, 0.0, 300.0);
   else
     outMult->AddAxis(multID, 101, 0.0, 101.0);
   
   // RAMA added lines below
   // =============>>
   Int_t multRefID = task->CreateValue(AliRsnMiniValue::kRefMult, kFALSE);
   AliRsnMiniOutput *outRefMult = task->CreateOutput("eventRefMult", "HIST", "EVENT");
   outRefMult->AddAxis(multRefID, 300, 0.0, 300.0);
   /*
   TH2F* hvz = new TH2F("hVzVsCent",Form("Vertex position vs centrality"), 101, 0., 101., 500, -50.0, 50.0);
   hvz->GetXaxis()->SetTitle("V0A");
   hvz->GetYaxis()->SetTitle("z_{vtx} (cm)");
   task->SetEventQAHist("vz", hvz);//plugs this histogram into the fHAEventVz data member

   TH2F* hRefMultiVsCent = new TH2F("hRefMultiVsCent",Form("Reference multiplicity vs centrality"), 101, 0., 101., 400, 0., 400.);
   hRefMultiVsCent->GetXaxis()->SetTitle("V0A");
   hRefMultiVsCent->GetYaxis()->SetTitle("GLOBAL");
   task->SetEventQAHist("refmulti",hRefMultiVsCent);//plugs this histogram into the fHAEventRefMultiCent data member
   */
   TH2F* hMultiVsCent = new TH2F("hMultiVsCent",Form("Multiplicity vs centrality"), 101, 0., 101., 400, 0., 400.);
   hMultiVsCent->GetXaxis()->SetTitle("V0A");
   hMultiVsCent->GetYaxis()->SetTitle("QUALITY");
   task->SetEventQAHist("multicent",hMultiVsCent);//plugs this histogram into the fHAEventMultiCent data member
   //  <<==============

   //
   // -- PAIR CUTS  -------------------------------------------------------------------------------
   //
   AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
   if(isPP) cutY->SetRangeD(-0.5, 0.5);
   else     cutY->SetRangeD(-0.465, 0.035);

   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());
   
   
   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------
   
   
   //for systematic checks

   gROOT->LoadMacro("ConfigLStarpPb5TeV.C");
   if (!ConfigLStar(task, isMC, isPP, "", cutsPair, aodFilterBit, customQualityCutsID, cutPrCandidate, cutKaCandidate, nsigmaPr, nsigmaKa,  enableMonitor, isMC&IsMcTrueOnly, signedPdg,monitorOpt, useCrossedRows, yaxisVar ,useMixLS)) 
     return 0x0;  


   /* gROOT->LoadMacro("ConfigLStar.C");
   
   if (!ConfigLStar(task, isMC, isPP, "", cutsPair, aodFilterBit,cutPrCandidate,cutKaCandidate, nsigmaPr,nsigmaKa, enableSyst, DCAxyFormula, dcazmax, minNcls, maxX2cls, TPCGlobalchi2,ITSchi2, minCrossedRows, maxClsCrossedRows, enableMonitor, isMC&IsMcTrueOnly, signedPdg,monitorOpt,useCrossedRows,yaxisVar ,useMixLS))
     return 0x0;
   */

   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()),  TList::Class(), AliAnalysisManager::kOutputContainer,outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);

   cout<<" taskname  =  "<<taskName.Data()<<endl;
   return task;
}
