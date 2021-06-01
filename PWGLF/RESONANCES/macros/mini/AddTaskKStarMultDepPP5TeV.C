/***************************************************************************
            adash@cern.ch - last modified on 03/04/2013
//  modified by rsingh@cern.ch -- 31-10-2020
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

enum pairYCutSet { kPairDefault=0,
		   kCentral //=1
                  };

enum eventCutSet { kEvtDefault=0,
		   kNoPileUpCut, //=1
		   kDefaultVtx12,//=2
		   kDefaultVtx8, //=3
		   kDefaultVtx5, //=4
		   kMCEvtDefault, //=5
		   kTriggered, //=6
		   kNoVzCut, //=7
		   kNoEvtSel, //=8
		   kINEL10, //=9
		   kIGZ10, //=10
           kIGZ //=11
                 };

AliRsnMiniAnalysisTask * AddTaskKStarMultDepPP5TeV
(
 Bool_t      isMC           = kFALSE,
 Bool_t      isPP           = kTRUE,
 Int_t       Strcut         = 2011,
 Int_t       MultId         = 0,
 Int_t       evtCutSetID    = 0,
 Int_t       pairCutSetID   = 0,
 TString     outNameSuffix  = "tofveto3stpc2s",
 Int_t       customQualityCutsID = AliRsnCutSetDaughterParticle::kDisableCustom,
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate=AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s,
 Float_t     nsigmaPi       = 2.0,
 Float_t     nsigmaK        = 2.0,
 Float_t     nsigmaTOF       = 3.0,
 Bool_t      enableMonitor  = kTRUE,
 Int_t       nmix           = 5,
 Float_t     maxDiffVzMix   = 1.0,
 Float_t     maxDiffMultMix = 5.0
 )
{  

  UInt_t      triggerMask  = AliVEvent::kINT7;
  Bool_t      rejectPileUp = kTRUE;
  Double_t    vtxZcut      = 10.0;//cm, default cut on vtx z

  if(evtCutSetID==eventCutSet::kDefaultVtx12) vtxZcut=12.0; //cm
  if(evtCutSetID==eventCutSet::kDefaultVtx8) vtxZcut=8.0; //cm
  if(evtCutSetID==eventCutSet::kDefaultVtx5) vtxZcut=5.0; //cm
  if(evtCutSetID==eventCutSet::kNoPileUpCut) rejectPileUp=kFALSE;
  if(evtCutSetID==eventCutSet::kNoVzCut) vtxZcut=1.e6;//off
  
  
  if(!isPP || isMC || MultId) rejectPileUp = kFALSE;


  //-------------------------------------------
  //pair cuts
  //-------------------------------------------
  Double_t    minYlab = -0.5;
  Double_t    maxYlab = 0.5;
  
  if(pairCutSetID==pairYCutSet::kCentral){//|y_cm|<0.3
    minYlab=-0.3; maxYlab=0.3;
  }
  
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskKStarMultDepPP5TeV", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   
   TString taskName = Form("KStar%s%s", (isPP? "pp" : "PbPb"), (isMC ? "MC" : "Data"));
   
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   if(evtCutSetID!=eventCutSet::kNoEvtSel && evtCutSetID!=eventCutSet::kINEL10 && evtCutSetID!=eventCutSet::kIGZ10 && evtCutSetID!=eventCutSet::kIGZ){
     task->UseESDTriggerMask(triggerMask); //ESD
     // task->SelectCollisionCandidates(triggerMask); //AOD
   }
   
   if(isPP) {
     if(MultId==1) task->UseMultiplicity("AliMultSelection_V0M");
     else if(MultId==2) task->UseMultiplicity("AliMultSelection_RefMult08");
     else if(MultId==3) task->UseMultiplicity("AliMultSelection_RefMult05");
     else task->UseMultiplicity("QUALITY");
   } else task->UseCentrality("V0M");   
   // set event mixing options
   task->UseContinuousMix();
   //task->UseBinnedMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   task->UseMC(isMC);
   ::Info("AddTaskKStarMultDepPP5TeV", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n", nmix, maxDiffVzMix, maxDiffMultMix));
   
   mgr->AddTask(task);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //  
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted

   AliRsnCutPrimaryVertex* cutVertex=0;
   if(evtCutSetID!=eventCutSet::kTriggered && evtCutSetID!=eventCutSet::kNoEvtSel && evtCutSetID!=eventCutSet::kIGZ){
     if(evtCutSetID==eventCutSet::kINEL10 || evtCutSetID==eventCutSet::kIGZ10){
       cutVertex=new AliRsnCutPrimaryVertex("cutVertex",vtxZcut,0,kFALSE);
       cutVertex->SetCheckGeneratedVertexZ();
       
     }else if(!MultId || fabs(vtxZcut-10.)>1.e-10){
       cutVertex=new AliRsnCutPrimaryVertex("cutVertex",vtxZcut,0,kFALSE);
       if(!MultId){
	 cutVertex->SetCheckZResolutionSPD();
	 cutVertex->SetCheckDispersionSPD();
	 cutVertex->SetCheckZDifferenceSPDTrack();
       }
     }
   }

  AliRsnCutEventUtils* cutEventUtils=0;
  if(evtCutSetID!=eventCutSet::kNoEvtSel && evtCutSetID!=eventCutSet::kINEL10){
    cutEventUtils=new AliRsnCutEventUtils("cutEventUtils",kTRUE,rejectPileUp);
    if(evtCutSetID==eventCutSet::kIGZ10 || evtCutSetID==eventCutSet::kIGZ) cutEventUtils->SetCheckInelGt0MC();
    else if(!MultId){
      cutEventUtils->SetCheckIncompleteDAQ();
      cutEventUtils->SetCheckSPDClusterVsTrackletBG();
    }else{
      //cutEventUtils->SetCheckInelGt0SPDtracklets();
      cutEventUtils->SetRemovePileUppA2013(kFALSE);
      if(evtCutSetID!=eventCutSet::kTriggered) cutEventUtils->SetCheckAcceptedMultSelection();
    }
  }

  if(isPP && (!isMC) && cutVertex){
    cutVertex->SetCheckPileUp(rejectPileUp);// set the check for pileup
    ::Info("AddTaskKStarMultDepPP5TeV", "%s", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp)?"ON":"OFF"));
  }
  
  
  
  // define and fill cut set for event cut
  AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
  if(cutEventUtils && cutVertex){
    eventCuts->AddCut(cutEventUtils);
    eventCuts->AddCut(cutVertex);
    eventCuts->SetCutScheme(Form("%s&%s",cutEventUtils->GetName(),cutVertex->GetName()));
  }else if(cutEventUtils && !cutVertex){
    eventCuts->AddCut(cutEventUtils);
    eventCuts->SetCutScheme(Form("%s",cutEventUtils->GetName()));
  }else if(!cutEventUtils && cutVertex){
    eventCuts->AddCut(cutVertex);
    eventCuts->SetCutScheme(Form("%s",cutVertex->GetName()));
  }
  
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
    outMult->AddAxis(multID, 110, 0.0, 110.0);
  
  //
  // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
  //
  
  AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  cutY->SetRangeD(minYlab,maxYlab);
  AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
  cutsPair->AddCut(cutY);
  cutsPair->SetCutScheme(cutY->GetName());
   
   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------
   //
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigKStarPP5TeV.C");
   if (!ConfigKStarPP5TeV(task, isMC, isPP, "", cutsPair, Strcut, customQualityCutsID,cutKaCandidate,nsigmaPi,nsigmaK, nsigmaTOF, enableMonitor)) return 0x0;

   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddTaskKStarMultDepPP5TeV - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}
