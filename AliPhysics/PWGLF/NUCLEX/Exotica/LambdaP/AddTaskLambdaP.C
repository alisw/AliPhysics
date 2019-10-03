/***************************************************************************
//            Modified by me - 
//
// Macro to configure the LambdaP analysis task 
// It calls all configs desired by the user, by means
// of the boolean switches defined in the first lines.
// ---
// Inputs:
//  1) flag to know if running on MC or data
//  2) collision system, whether pp, pPb or PbPb
// --
// Returns:
//  kTRUE  --> initialization successful
//  kFALSE --> initialization failed (some config gave errors)
//
****************************************************************************/
 
enum ERsnCollType_t { kPP=0,
		      kPPb,
		      kPbPb};

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
		   kDefaultVtx5 //=5                    
                 };

enum eventMixConfig { kDisabled = -1,
		      kMixDefault,//=0 //10 events, Dvz = 1cm, DC = 10
		      k5Evts, //=1 //5 events, Dvz = 1cm, DC = 10
		      k5Cent,  //=2 //10 events, Dvz = 1cm, DC = 5
                    };

AliRsnMiniAnalysisTask *AddTaskLambdaP
(
 Bool_t      isMC,
 Int_t       collSyst,
 Float_t     cutV = 10.0,
 Int_t       evtCutSetID = 0,
 Int_t       pairCutSetID = 0,
 Int_t       mixingConfigID = 0,
 Int_t       aodFilterBit = 5,
 Float_t     piPIDCut = 3.0,
 Float_t     pPIDCut = 3.0,
 Float_t     trackDCAcut = 7.0,
 Float_t     massTol = 0.01,
 Float_t     lambdaDCA = 0.3,
 Float_t     lambdaCosPoinAn = 0.97,
 Float_t     lambdaDaughDCA = 0.5,
 Int_t       NTPCcluster = 70,
 Float_t     maxDiffVzMix = 1.0,
 Float_t     maxDiffMultMix = 10.0,
 Float_t     maxDiffAngleMixDeg = 20.0,
 Float_t     minDCAToVertexXYlambdaDaugh = 0.15,
 Int_t       aodN = 68,
 TString     outNameSuffix = "LambdaP",
 Int_t       centr = 0
 )
{  
  
  //-------------------------------------------
  // event cuts
  // Note that some default values refer to pPb data 2013
  // settings from AddTaskKStarPPB.C by Francesca Bellini
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

  //-------------------------------------------
  //pair cuts
  //-------------------------------------------
  Double_t    minYlab =  -0.5;
  Double_t    maxYlab =  0.5;

  if( collSyst==kPPb ) {
    if (pairCutSetID==pairYCutSet::kPairDefault) { //0<y_cm<0.5
      minYlab =  -0.465;  maxYlab =  0.035;
    }
    
    if (pairCutSetID==pairYCutSet::kNegative) { //-0.5<y_cm<0.0
      minYlab = -0.965;    maxYlab = -0.465;
    }
    
    if (pairCutSetID==pairYCutSet::kCentral) { //|y_cm|<0.3
      minYlab = -0.765;    maxYlab = -0.165;
    }
  } 
  ::Info("AddTaskLambdaP", Form("Min rapidity = %6.3f, Max rapidity = %6.3f",  minYlab, maxYlab) );

  //-------------------------------------------
  //mixing settings
  //-------------------------------------------

  Int_t       nmix = 0;
  //  Float_t     maxDiffVzMix = 1.0;
  //  Float_t     maxDiffMultMix = 10.0;
  
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
      ::Error("AddTaskLambdaP", "No analysis manager to connect to.");
      return NULL;
   } 
   
   // create the task and configure 
   TString collSystName;
   if(collSyst==kPP) collSystName="pp";
   else if(collSyst==kPPb) collSystName="pPb";
   else collSystName="PbPb";

   TString taskName = Form("LambdaP%s%s_%.1f_%d_%.1f_%.1f_%.2f_%.4f_%.2f_%.2f_%.1f_%.2f", 
			   collSystName.Data(), (isMC ? "MC" : "Data"),cutV,NTPCcluster,piPIDCut,pPIDCut,trackDCAcut,massTol,lambdaDCA,lambdaCosPoinAn,lambdaDaughDCA,minDCAToVertexXYlambdaDaugh);

   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   if (!isMC && (collSyst==kPP) ){
     Printf(Form("========== SETTING USE CENTRALITY PATCH AOD049 : %s", (aodN==49)? "yes" : "no"));
     task->SetUseCentralityPatch(aodN==49);
   }
   
   if(collSyst==kPPb)  task->UseESDTriggerMask(triggerMask);
   else if(collSyst==kPbPb) {
     if (centr == 1) { task->UseESDTriggerMask(AliVEvent::kCentral); }
     else  if (centr == 2) {  task->UseESDTriggerMask(AliVEvent::kSemiCentral);}
     else  if (centr == 3) {  task->UseESDTriggerMask(AliVEvent::kMB); }
     else { task->UseESDTriggerMask(AliVEvent::kMB  | AliVEvent::kCentral | AliVEvent::kSemiCentral); }
   }
   
   if(collSyst==kPPb) 
     task->SelectCollisionCandidates(triggerMask); //
   else if ( collSyst == kPP ) 
     task->SelectCollisionCandidates(AliVEvent::kMB); //
   else {
     if (centr == 1) { 
       task->SelectCollisionCandidates(AliVEvent::kCentral); }
     if (centr == 2) { 
       task->SelectCollisionCandidates(AliVEvent::kSemiCentral); }
     if (centr == 3) { 
       task->SelectCollisionCandidates(AliVEvent::kMB); }
     else { 
       task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral); }
   }
 
   
   if ( collSyst == kPP ) 
     task->UseMultiplicity("QUALITY");
   else if(collSyst==kPPb) 
    task->UseCentrality("V0A");   
   else
     task->UseCentrality("V0M");   

   // set event mixing options
   task->UseContinuousMix();
   //task->UseBinnedMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);

   if (! (collSyst==kPP) ) task->SetMaxDiffAngle(maxDiffAngleMixDeg*TMath::DegToRad()); //set angle diff in rad
  ::Info("AddTaskLambdaP", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n max diff EP angle = %5.3f deg", nmix, maxDiffVzMix, maxDiffMultMix, ( !collSyst ? 0.0 : maxDiffAngleMixDeg)));
   
   mgr->AddTask(task);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //  
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted

   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", cutV, 0, kFALSE);
   if ( collSyst == kPP ) cutVertex->SetCheckPileUp(kTRUE);   // set the check for pileup
   
   if (collSyst==kPPb) { 
     //set check for pileup in 2013
     AliRsnCutEventUtils *cutEventUtils = new AliRsnCutEventUtils("cutEventUtils", rmFirstEvtChunk, rejectPileUp);
     cutEventUtils->SetUseVertexSelection2013pA(useVtxCut2013pA);
     ::Info("AddTaskLambdaP", Form(":::::::::::::::::: Vertex cut as pA 2013: %s", (useVtxCut2013pA?"ON":"OFF")));   
     if (useMVPileUpSelection){
       cutEventUtils->SetUseMVPlpSelection(useMVPileUpSelection);
       cutEventUtils->SetMinPlpContribMV(MinPlpContribMV);
       cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
       ::Info("AddTaskLambdaP", Form("Multiple-vtx Pile-up rejection settings: MinPlpContribMV = %i, MinPlpContribSPD = %i", 
				       MinPlpContribMV, MinPlpContribSPD));
     } else {
       cutEventUtils->SetMinPlpContribSPD(MinPlpContribSPD);
       ::Info("AddTaskLambdaP", Form("SPD Pile-up rejection settings: MinPlpContribSPD = %i", MinPlpContribSPD));
     }
     ::Info("AddTaskLambdaP", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp?"ON":"OFF")));   
     ::Info("AddTaskLambdaP", Form("::::::::::::: Remove first event in chunk: %s", (rmFirstEvtChunk?"ON":"OFF")));   
   }
   
   // define and fill cut set for event cut
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   if (collSyst==kPPb) eventCuts->AddCut(cutEventUtils);
   else eventCuts->AddCut(cutVertex);
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
   if (collSyst==kPP) 
     outMult->AddAxis(multID, 400, 0.0, 400.0);
   else
     outMult->AddAxis(multID, 100, 0.0, 100.0);
   
   //event plane (only for PbPb)
   Int_t planeID = task->CreateValue(AliRsnMiniValue::kPlaneAngle, kFALSE);
   AliRsnMiniOutput *outPlane = 0x0; //task->CreateOutput("eventPlane", "HIST", "EVENT");
   if ( collSyst==kPbPb ){
     outPlane = task->CreateOutput("eventPlane", "HIST", "EVENT");
     outPlane->AddAxis(planeID, 180, 0.0, TMath::Pi());
   }
   

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
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/STRANGENESS/Hypernuclei/ConfigLambdaP.C");
   if (isMC) {
     Printf("========================== MC analysis - PID cuts not used");
     
   } else 
     Printf("========================== DATA analysis - PID cuts used");
   
   if (!ConfigLambdaP(task, collSyst, isMC, piPIDCut, pPIDCut, aodFilterBit, trackDCAcut, massTol, lambdaDCA, lambdaCosPoinAn, lambdaDaughDCA, NTPCcluster, minDCAToVertexXYlambdaDaugh,"", cutsPair)) return 0x0;
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //

   
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   outputFileName += ":Rsn";
   Printf("AddTaskLambdaP - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s_%.1f_%d_%.1f_%.1f_%.2f_%.4f_%.2f_%.2f_%.1f_%.2f",
								outNameSuffix.Data(),cutV,NTPCcluster,piPIDCut,pPIDCut,
								trackDCAcut,massTol,lambdaDCA,lambdaCosPoinAn,lambdaDaughDCA,minDCAToVertexXYlambdaDaugh), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
  return task;
}
