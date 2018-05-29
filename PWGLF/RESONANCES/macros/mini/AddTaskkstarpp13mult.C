/***************************************************************************
            skundu@cern.ch - last modified on 01/01/2018
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

AliRsnMiniAnalysisTask * AddTaskkstarpp13mult
(
 Bool_t      isMC = kFALSE,
 Bool_t      isPP = kTRUE,
 Bool_t      rejectPileUp = kTRUE,
 Int_t       aodFilterBit=5,
 Int_t       nmix = 0,
 Float_t     nsigmaPi = 2.0,
 Float_t     nsigmaK = 2.0,
 Bool_t      enableMonitor = kTRUE,
 Float_t     maxDiffVzMix = 1.0,
 Float_t     maxDiffMultMix = 10.0,
 TString     outNameSuffix = "",
 TString     optSy="Default",
 Int_t       evtCutSetID=0
 )
{  
  
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskPhiPPb_TPC", "No analysis manager to connect to.");
      return NULL;
   } 




  UInt_t      triggerMask=AliVEvent::kINT7;
  if(evtCutSetID>=100){
    triggerMask=AliVEvent::kHighMultV0;
    // evtCutSetID=evtCutSetID%100;
  }


   
   // create the task and configure 
   
   TString taskName = Form("TPCPhiMeson%s%s", (isPP? "pp" : "pPb"), (isMC ? "MC" : "Data"));
   
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   //  task->SelectCollisionCandidates(triggerMask);
      task->UseESDTriggerMask(triggerMask);
      
   // task->UseESDTriggerMask(AliVEvent::kINT7); i have did this
   if (isPP) 
     task->UseMultiplicity("AliMultSelection_V0M");
   else
     task->UseCentrality("V0M");   
   // set event mixing options
   task->UseContinuousMix();
   //task->UseBinnedMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   task->UseMC(isMC);
   ::Info("AddTaskPhiPPb_TPC", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n", nmix, maxDiffVzMix, maxDiffMultMix));
   
   mgr->AddTask(task);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //  
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
   /*  AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
     if (isPP) cutVertex->SetCheckPileUp(kTRUE);   // set the check for pileup
      cutVertex->SetCheckZResolutionSPD();
     cutVertex->SetCheckDispersionSPD();
   cutVertex->SetCheckZDifferenceSPDTrack();
   AliRsnCutEventUtils* cutEventUtils=new AliRsnCutEventUtils("cutEventUtils",kTRUE,rejectPileUp);
   //cutEventUtils->SetCheckIncompleteDAQ();
   //if(aodFilterBit<200) cutEventUtils->SetCheckIncompleteDAQ();
   //cutEventUtils->SetCheckSPDClusterVsTrackletBG();
   cutEventUtils->SetRemovePileUppA2013(kFALSE);
   cutEventUtils->SetCheckAcceptedMultSelection();*/
   // define and fill cut set for event cut
   /*   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutVertex);
 eventCuts->AddCut(cutEventUtils);
 eventCuts->SetCutScheme(Form("%s&%s",cutEventUtils->GetName(),cutVertex->GetName()));
   // set cuts in task
   task->SetEventCuts(eventCuts);*/
 
   AliRsnCutEventUtils* cutEventUtils=new AliRsnCutEventUtils("cutEventUtils",kTRUE,rejectPileUp);
   cutEventUtils->SetRemovePileUppA2013(kFALSE);
   cutEventUtils->SetCheckAcceptedMultSelection();
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   
   eventCuts->AddCut(cutEventUtils);
   eventCuts->SetCutScheme(Form("%s",cutEventUtils->GetName()));
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
     outMult->AddAxis(multID, 110, 0.0, 110.0);

   /*
   Double_t multbins[200];
   int j,nmult=0;
   for(j=0;j<10;j++){multbins[nmult]=0.0001*j; nmult++;}
   for(j=1;j<10;j++){multbins[nmult]=0.001*j; nmult++;}
   for(j=1;j<10;j++){multbins[nmult]=0.01*j; nmult++;}
   for(j=1;j<10;j++){multbins[nmult]=0.1*j; nmult++;}
   for(j=1;j<=100;j++){multbins[nmult]=j; nmult++;}
   nmult--;
   */
   /*
   Double_t multbins[200];
   int j,nmult=0;
   for(j=0;j<10;j++){multbins[nmult]=0.001*j; nmult++;}
   for(j=1;j<10;j++){multbins[nmult]=0.01*j; nmult++;}
   for(j=1;j<10;j++){multbins[nmult]=0.1*j; nmult++;}
   for(j=1;j<=100;j++){multbins[nmult]=j; nmult++;}
   nmult--;
   */
   
   Double_t multbins[200];
   int j,nmult=0;
   if(triggerMask==AliVEvent::kHighMultV0){
     for(j=0;j<10;j++){multbins[nmult]=0.001*j; nmult++;}
     for(j=1;j<10;j++){multbins[nmult]=0.01*j; nmult++;}
     for(j=1;j<=10;j++){multbins[nmult]=0.1*j; nmult++;}
   }else{
     for(j=0;j<10;j++){multbins[nmult]=0.1*j; nmult++;}
     for(j=1;j<=100;j++){multbins[nmult]=j; nmult++;}
  }
   nmult--;

   
   TH1F* hEventsVsMulti=new TH1F("hAEventsVsMulti","",nmult,multbins);
   task->SetEventQAHist("EventsVsMulti",hEventsVsMulti);//custom binning for fHAEventsVsMulti
   
   double ybins[500];
   for(j=0;j<=401;j++) ybins[j]=j-0.5;
   TH2F* hmc=new TH2F("MultiVsCent","", nmult,multbins, 401,ybins);
   hmc->GetYaxis()->SetTitle("QUALITY");
   task->SetEventQAHist("multicent",hmc);//plugs this histogram into the fHAEventMultiCent data member

   
   
   //
   // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
   //

   AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
   cutY->SetRangeD(-0.5, 0.5);
   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());
   
   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------
   //
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/Configkstarpp13mult.C");
   //gROOT->LoadMacro("/home/sourav/alice/ali-master/AliPhysics/PWGLF/RESONANCES/macros/mini/Configkstarpp13mult.C");
   // gROOT->LoadMacro("Configkstarpp13mult.C");
   if (!Configkstarpp13mult(task, isMC, isPP, "", cutsPair, nsigmaPi,nsigmaK, enableMonitor,optSy,triggerMask)) return 0x0;

   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddTaskPhiPPb_TPC - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}
