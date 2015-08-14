/***************************************************************************
   sarita.sahoo@cern.ch - last modified on 25/06/2015
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

AliRsnMiniAnalysisTask * AddAnalysisTaskLStar
(
 Bool_t      isMC ="",
 Bool_t      isPP ="",
 Int_t       aodFilterBit = 10,
 Bool_t      IsMcTrueOnly = kTRUE,
 Bool_t      enableMonitor = kFALSE,
 const char* systuncert = "Default_PID25_Cls70_AntiLStar"
)
{
   TString opt(systuncert);
   opt.ToUpper();
   
   Bool_t isDefault      = opt.Contains("DEFAULT");
   Bool_t isPID          = opt.Contains("PID");
   Bool_t isCls          = opt.Contains("CLS");
   Bool_t isDCAxy        = opt.Contains("DCAXY");
   Bool_t isDCAz         = opt.Contains("DCAZ");
   Bool_t isChi2TPC      = opt.Contains("CHI2TPC");
   Bool_t isChi2TPCGlobal= opt.Contains("CHI2TPCGLOBAL");
   Bool_t isChi2ITS      = opt.Contains("CHI2ITS");
   Bool_t isAntiLstar    = opt.Contains("ANTILSTAR");
 

   AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar;
   AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar;
   
   Float_t     nsigmaPr = 2.0;
   Float_t     nsigmaKa = 2.0;
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
   Int_t                  signedPdg = 3124;
   TString                monitorOpt = "NoSIGN";  //Flag for AddMonitorOutput.C e.g."NoSIGN"
   Bool_t                 useCrossedRows = kFALSE;
   AliRsnMiniValue::EType yaxisVar = AliRsnMiniValue::kPt; //
   Bool_t                 useMixLS = 0;
   int nmix = 5;



   if(isDefault)
     AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,dcazmax,minNcls,maxX2cls,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,signedPdg,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"Default_LStar");
   
   if(isDefault && isAntiLstar && isMC)
     {
       AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,dcazmax,minNcls,maxX2cls,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,-3124,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"Default_ALStar");
     }
   if(isPID)
     {
       AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,1.5,1.5,enableSyst,DCAxyFormula,dcazmax,minNcls,maxX2cls,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,signedPdg,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"PID15Sig_LStar");
       AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,2.5,2.5,enableSyst,DCAxyFormula,dcazmax,minNcls,maxX2cls,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,signedPdg,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"PID25Sig_LStar");
     }
   if(isPID && isAntiLstar && isMC)
     {
       AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,1.5,1.5,enableSyst,DCAxyFormula,dcazmax,minNcls,maxX2cls,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,-3124,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"PID15Sig_ALStar");
       AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,2.5,2.5,enableSyst,DCAxyFormula,dcazmax,minNcls,maxX2cls,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,-3124,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"PID25Sig_ALStar");
     }
   if(isCls)
     {
       AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,dcazmax,60,maxX2cls,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,signedPdg,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"Cls60_LStar");
       AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,dcazmax,70,maxX2cls,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,signedPdg,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"Cls70_LStar");
     }
   if(isCls && isAntiLstar && isMC)
     {
       AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,dcazmax,60,maxX2cls,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,-3124,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"Cls60_ALStar");
       AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,dcazmax,70,maxX2cls,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,-3124,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"Cls70_ALStar");
     }
   
   if(isDCAxy)
     {
     AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula5,dcazmax,minNcls,maxX2cls,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,signedPdg,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"DCAxy5Sig_LStar");
     AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula4,dcazmax,minNcls,maxX2cls,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,signedPdg,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"DCAxy4Sig_LStar");
     }
   if(isDCAxy && isAntiLstar && isMC)
     {
     AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula5,dcazmax,minNcls,maxX2cls,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,-3124,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"DCAxy5Sig_ALStar");
     AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula4,dcazmax,minNcls,maxX2cls,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,-3124,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"DCAxy4Sig_ALStar");
     }

   if(isDCAz)
     {
       AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,1.5,minNcls,maxX2cls,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,signedPdg,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"DCAz15cm_LStar");
       AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,1.0,minNcls,maxX2cls,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,signedPdg,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"DCAz10cm_LStar");
     }
   if(isDCAz && isAntiLstar && isMC)
     {
     AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,1.5,minNcls,maxX2cls,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,-3124,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"DCAz15cm_ALStar");
     AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,1.0,minNcls,maxX2cls,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,-3124,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"DCAz10cm_ALStar");
     }


   if(isChi2TPC)
     {
     AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,dcazmax,minNcls,2.5,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,signedPdg,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"Chi2TPC25_LStar");
     AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,dcazmax,minNcls,3.0,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,signedPdg,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"Chi2TPC30_LStar");

     }

   if(isChi2TPC && isAntiLstar && isMC)
     {
     AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,dcazmax,minNcls,2.5,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,-3124,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"Chi2TPC25_ALStar");
     AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,dcazmax,minNcls,3.0,TPCGlobalchi2,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,-3124,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"Chi2TPC30_ALStar");
     }

   if(isChi2TPCGlobal)
     {
     AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,dcazmax,minNcls,maxX2cls,25,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,signedPdg,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"Chi2TPCGlobal25_LStar");
     AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,dcazmax,minNcls,maxX2cls,30,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,signedPdg,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"Chi2TPCGlobal30_LStar");

     }

   if(isChi2TPCGlobal && isAntiLstar && isMC)
     {
     AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,dcazmax,minNcls,maxX2cls,25,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,-3124,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"Chi2TPCGlobal25_ALStar");
     AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,dcazmax,minNcls,maxX2cls,30,ITSchi2,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,-3124,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"Chi2TPCGlobal30_ALStar");
     }

   if(isChi2ITS)
     {
     AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,dcazmax,minNcls,maxX2cls,TPCGlobalchi2,25,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,signedPdg,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"Chi2ITS25_LStar");
     AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,dcazmax,minNcls,maxX2cls,TPCGlobalchi2,30,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,signedPdg,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"Chi2ITS30_LStar");

     }

   if(isChi2ITS && isAntiLstar && isMC)
     {
     AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,dcazmax,minNcls,maxX2cls,TPCGlobalchi2,25,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,-3124,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"Chi2ITS25_ALStar");
     AddAnalysisTask(isMC,isPP,aodFilterBit,cutPrCandidate,cutKaCandidate,nsigmaPr,nsigmaKa,enableSyst,DCAxyFormula,dcazmax,minNcls,maxX2cls,TPCGlobalchi2,30,minCrossedRows,maxClsCrossedRows,enableMonitor,IsMcTrueOnly,triggerMask,-3124,monitorOpt,useCrossedRows,yaxisVar,useMixLS,nmix,1.0,10.0,"Chi2ITS30_ALStar");
     }
}


AliRsnMiniAnalysisTask  * AddAnalysisTask
(
   Bool_t      isMC,
   Bool_t      isPP,
   Int_t       aodFilterBit = 5,
   AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate = AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,
   AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,
   Float_t     nsigmaPr = 2.0,
   Float_t     nsigmaKa = 2.0,
   Bool_t      enableSyst = kFALSE,
   Char_t      DCAxyFormula[100] = "0.0182+0.035/pt^1.01",
   Double_t    dcazmax = 2,
   Double_t    minNcls = 70,
   Double_t    maxX2cls = 4.0,
   Float_t     TPCGlobalchi2 = 36,
   Double_t    ITSchi2 = 36,
   Double_t    minCrossedRows = 70.0,
   Double_t    maxClsCrossedRows = 0.8,
   Bool_t      enableMonitor = kTRUE,
   Bool_t      IsMcTrueOnly = kFALSE,
   UInt_t      triggerMask = AliVEvent::kMB,
   Int_t                  signedPdg = 3124,
   TString                monitorOpt = "NoSIGN",  //Flag for AddMonitorOutput.C e.g."NoSIGN"
   Bool_t                 useCrossedRows = kFALSE,
   AliRsnMiniValue::EType yaxisVar = AliRsnMiniValue::kPt,
   Bool_t                 useMixLS = 0,
   Int_t       nmix = 5,
   Float_t     maxDiffVzMix = 1.0,
   Float_t     maxDiffMultMix = 10.0,
   TString     outNameSuffix = ""
 )

{  
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //


   TString outputFileName = Form("%s", AliAnalysisManager::GetCommonFileName());
  // Objects name


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddAnalysisTaskTPCKStarSyst", "No analysis manager to connect to.");
      return NULL;
   } 


   // create the task and configure 
   //TString taskName = Form("LStar%s%s_%s", (isPP? "pp" : "pPb"), (isMC ? "MC" : "Data"), outNameSuffix.Data() );
   //TString taskName = Form("LStar%s%s_%s", (isPP? "pp" : "pPb"), (isMC ? "MC" : "Data"), outNameSuffix.Data() );
   TString taskName = outNameSuffix.Data();
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);


   //if(is2011PbPb)
   //task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
   //else
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
   if (isPP) {
     AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
     cutVertex->SetCheckPileUp(kTRUE);   // set the check for pileup
     
     // define and fill cut set for event cut
     AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
     eventCuts->AddCut(cutVertex);
     eventCuts->SetCutScheme(cutVertex->GetName());
   }
   else {
   AliRsnCutEventUtils *cutEventUtils = new AliRsnCutEventUtils("cutEventUtils", kTRUE, kTRUE);
   cutEventUtils->SetUseVertexSelection2013pA(kTRUE);
   cutEventUtils->SetMinPlpContribSPD(5);      
   // define and fill cut set for event cut
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutEventUtils);
   eventCuts->SetCutScheme(cutEventUtils->GetName());
   }

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
     outMult->AddAxis(multID, 100, 0.0, 100.0);
   
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
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigLStar.C");
   
   if (!ConfigLStar(task, isMC, isPP, "", cutsPair, aodFilterBit,cutPrCandidate,cutKaCandidate, nsigmaPr,nsigmaKa, enableSyst, DCAxyFormula, dcazmax, minNcls, maxX2cls, TPCGlobalchi2,ITSchi2, minCrossedRows, maxClsCrossedRows, enableMonitor, isMC&IsMcTrueOnly, signedPdg,monitorOpt,useCrossedRows,yaxisVar ,useMixLS))
     return 0x0;
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()),  TList::Class(), AliAnalysisManager::kOutputContainer,outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);

   cout<<" taskname  =  "<<taskName.Data()<<endl;
   return task;
}
