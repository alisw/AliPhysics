/***************************************************************************
  priyanka.sett@cern.ch - last modified on 27/10/2016
  for L* in pp 13 TeV analysis


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

enum eventCutSet { kDefaultVtx = 0,
		   kDefaultVtx12,//=1
		   kDefaultVtx8 //=2
                 };



AliRsnMiniAnalysisTask * AddTaskLstar13TeVpp
(
   Bool_t      isMC = kFALSE,
   Bool_t      isPP = kTRUE,
   Int_t       aodFilterBit = 5,
   Int_t       evtCutSetID = 0, 
   Int_t       MultBins = 1,// default for V0_M and MultBins = 2 for RefMult0_8 
   Int_t       customQualityCutsID=1, // for default
   AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV,
   AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV,
   Float_t     nsigmaPr = 2.0,
   Float_t     nsigmaKa = 2.0,
   Bool_t      enableMonitor = kTRUE,
   Bool_t      IsMcTrueOnly = kFALSE,
   UInt_t      triggerMask = AliVEvent::kINT7,
   Int_t       signedPdg = 3124,
   TString     monitorOpt = "NoSIGN",  //Flag for AddMonitorOutput.C e.g."NoSIGN"
   Bool_t      useCrossedRows = kTRUE,
   const char *yaxisVar = "",  //yaxisVar = "PtDaughter_PDaughter_cent"
   Bool_t      useMixLS = 0,
   Int_t       nmix = 5,
   Float_t     maxDiffVzMix = 1.0,
   Float_t     maxDiffMultMix = 10.0,
   TString     outNameSuffix = "Default"
)


{  
  //-------------------------------------------
  // event cuts
  //-------------------------------------------
 
  Double_t  vtxZcut = 10.0;//default cut on vtx z
  
  //  if(evtCutSetID==eventCutSet::kDefaultVtx) vtxZcut=10.0; //cm
  if(evtCutSetID==eventCutSet::kDefaultVtx12) vtxZcut=12.0; //cm
  if(evtCutSetID==eventCutSet::kDefaultVtx8) vtxZcut=8.0; //cm
  
  //vtxZcut = 10.0;//default cut on vtx z

  
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskLstar13TeVpp", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   TString taskName = Form("LStar%s%s_%i%i_%s", (isPP? "pp" : "pPb"), (isMC ? "MC" : "Data"), (Int_t)cutPrCandidate,(Int_t)cutKaCandidate, outNameSuffix.Data() );
   
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);


   //if(is2011PbPb)
   //task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
   //else
   task->SelectCollisionCandidates(triggerMask);

   /*
   AliMultSelectionTask *taskm = AddTaskMultSelection();

   if(isMC == 1)
     {
       taskm->SetAlternateOADBforEstimators( "LHC15f" );     
     }
   */

   if (isPP) {
     //     task->UseMultiplicity("QUALITY");
     if (MultBins == 1) task->UseMultiplicity("AliMultSelection_V0M"); // for multiplicity percentile
     else if(MultBins == 2) task->UseMultiplicity("AliMultSelection_RefMult05");
     else task->UseMultiplicity("QUALITY");
   }
   else
     task->UseCentrality("V0M");   





   // set event mixing options
   task->UseContinuousMix();
   //task->UseBinnedMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   //if (!isPP) task->SetMaxDiffAngle(maxDiffAngleMixDeg*TMath::DegToRad()); //set angle diff in rad
   ::Info("AddTaskLstar13TeVpp", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n ", nmix, maxDiffVzMix, maxDiffMultMix));
   
   mgr->AddTask(task);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //  
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
   //   if (isPP)
   //{

   Bool_t  rejectPileUp = kTRUE;

   //   if(!isPP || isMC || MultBins) rejectPileUp=kFALSE;
   if(!isPP || isMC) rejectPileUp=kFALSE;

   
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", vtxZcut, 0, kFALSE);
   cutVertex->SetCheckZResolutionSPD();
   cutVertex->SetCheckDispersionSPD();
   cutVertex->SetCheckZDifferenceSPDTrack();
   
   
   AliRsnCutEventUtils *cutEventUtils = new AliRsnCutEventUtils("cutEventUtils", kTRUE, rejectPileUp);
   //   cutEventUtils->SetCheckIncompleteDAQ();
   //   cutEventUtils->SetCheckSPDClusterVsTrackletBG();   
   if(!MultBins){
     cutEventUtils->SetCheckIncompleteDAQ();
     cutEventUtils->SetCheckSPDClusterVsTrackletBG();
     cutEventUtils->SetCheckInelGt0SPDtracklets();
   }
   else{
     cutEventUtils->SetRemovePileUppA2013(kFALSE);
     cutEventUtils->SetCheckAcceptedMultSelection();
   }


   
   
   // if(isPP && (!isMC)){ 
   if(isPP && (!isMC) && cutVertex){ // modified on 21 July. Ref Anders code on git
     cutVertex->SetCheckPileUp(rejectPileUp);// set the check for pileup  
     ::Info("AddTaskLstar13TeVpp", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp)?"ON":"OFF"));   
   }


   // define and fill cut set for event cut
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutEventUtils);
   eventCuts->AddCut(cutVertex);
   eventCuts->SetCutScheme(Form("%s&%s",cutEventUtils->GetName(),cutVertex->GetName()));
   cout<< "--------------"<< cutVertex <<"cm ---------------- "<<endl;
   
   
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
   if (isPP && !MultBins)   outMult->AddAxis(multID, 300, 0.0, 300.0);
   else  outMult->AddAxis(multID, 110, 0.0, 110.0);
   
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
     {
       //gROOT->LoadMacro("$ALICE_ROOT/PWGLF/RESONANCES/macros/mini/ConfigLStar.C");
       //   gROOT->LoadMacro("ConfigureLstar13TeVpp.C");
       gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigureLstar13TeVpp.C");
       if (!ConfigureLstar13TeVpp(task, isMC, isPP, "", cutsPair, aodFilterBit, customQualityCutsID, cutPrCandidate, cutKaCandidate, nsigmaPr, nsigmaKa,  enableMonitor, isMC&IsMcTrueOnly, signedPdg, monitorOpt, useCrossedRows, yaxisVar ,useMixLS)) 
return 0x0;  
     }

   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddAnalysisTaskTPCKStarSyst - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}
