/***************************************************************************
              fbellini@cern.ch - last modified on 28/11/2013
//
//Lauches KStar analysis with rsn mini package
//Allows basic configuration of pile-up check and event cuts
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
		   kSpecial1, //=6                   
		   kSpecial2, //=7
		   kNoEvtSel, //=8 
		   kSpecial3//=9
                 };

enum eventMixConfig { kDisabled = -1,
		      kMixDefault,//=0 //10 events, Dvz = 1cm, DC = 10
		      k5Evts, //=1 //5 events, Dvz = 1cm, DC = 10
		      k5Cent,  //=2 //10 events, Dvz = 1cm, DC = 5
		      k5Evts5Cent
                    };


AliRsnMiniAnalysisTask * AddTaskKStarPP8TeV_PID
(
 Bool_t      isMC = kFALSE,
 Bool_t      isPP = kTRUE,
 TString     outNameSuffix = "tpc2stof3sveto",
 Int_t       evtCutSetID = 0,
 Int_t       pairCutSetID = 0,
 Int_t       mixingConfigID = 0,
 Int_t       aodFilterBit = 5,
 Int_t       customQualityCutsID = -1,
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPiCandidate = AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s,
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s,
 Float_t     nsigmaPi = 2.0,
 Float_t     nsigmaKa = 2.0,
 Bool_t      enableMonitor = kTRUE,
 Bool_t      IsMcTrueOnly = kFALSE,
 TString     monitorOpt = "NoSIGN",
 Bool_t      useMixLS = 0,
 Bool_t      checkReflex = 0,
 AliRsnMiniValue::EType yaxisvar = AliRsnMiniValue::kPt
)
{  

  
  //-------------------------------------------
  // event cuts
  //-------------------------------------------
  UInt_t      triggerMask = AliVEvent::kINT7;//A Khuntia
  // if(isMC && (evtCutSetID==eventCutSet::kNoEvtSel || evtCutSetID==eventCutSet::kSpecial3)) triggerMask=AliVEvent::kAny;
  Bool_t      rejectPileUp = kTRUE; //
  Double_t    vtxZcut = 10.0; //cm, default cut on vtx z
  
  if (evtCutSetID==eventCutSet::kDefaultVtx12){vtxZcut = 12.0;} //cm

  if (evtCutSetID==eventCutSet::kDefaultVtx8){vtxZcut = 8.0;} //cm
  
  if (evtCutSetID==eventCutSet::kDefaultVtx5){vtxZcut = 5.0;}//cm
    
  if (evtCutSetID==eventCutSet::kNoPileUpCut){rejectPileUp=kFALSE;}//cm
  
  if(evtCutSetID==eventCutSet::kSpecial2) vtxZcut=1.e6;//off

  //-------------------------------------------
  //pair cuts
  //-------------------------------------------
  Double_t    minYlab =  -0.5;
  Double_t    maxYlab =  0.5;
  
  if (pairCutSetID==pairYCutSet::kCentral) { //|y_cm|<0.3
    minYlab = -0.3;    maxYlab = 0.3;
  }

  //-------------------------------------------
  //mixing settings
  //-------------------------------------------
  Int_t       nmix = 0;
  Float_t     maxDiffVzMix = 1.0;
  Float_t     maxDiffMultMix = 10.0;
  
  if (mixingConfigID == eventMixConfig::kMixDefault) { nmix = 10;}

  if (mixingConfigID == eventMixConfig::k5Evts) {nmix = 5;}
  
  if (mixingConfigID == eventMixConfig::k5Cent) {maxDiffMultMix = 5;}
  
  if(mixingConfigID==eventMixConfig::k5Evts5Cent){nmix=5; maxDiffMultMix=5;}

  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddAnalysisTaskTOFKStar", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   TString taskName = Form("TOFKStar%s%s_%i%i", (isPP? "pp" : "PbPb"), (isMC ? "MC" : "Data"), (Int_t)cutPiCandidate,(Int_t)cutKaCandidate );
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   //task->UseESDTriggerMask(triggerMask); //ESD
   //task->SelectCollisionCandidates(triggerMask); //AOD
   if(evtCutSetID!=eventCutSet::kNoEvtSel && evtCutSetID!=eventCutSet::kSpecial3) task->SelectCollisionCandidates(triggerMask); //AOD

   
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
   ::Info("AddAnalysisTaskTOFKStar", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f", nmix, maxDiffVzMix, maxDiffMultMix));
   
   mgr->AddTask(task);
   
   //
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors to vtx
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
   AliRsnCutPrimaryVertex *cutVertex=0;
    if (evtCutSetID!=eventCutSet::kSpecial1 && evtCutSetID!=eventCutSet::kNoEvtSel){
    cutVertex = new AliRsnCutPrimaryVertex("cutVertex", vtxZcut, 0, kFALSE);
    if (evtCutSetID==eventCutSet::kSpecial3) cutVertex->SetCheckGeneratedVertexZ();

     }//vertex loop


   
   if (isPP && (!isMC)&&cutVertex) { 
     cutVertex->SetCheckPileUp(rejectPileUp);   // set the check for pileup  
     ::Info("AddAnalysisTaskTOFKStar", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp)?"ON":"OFF"));
     //cutVertex->SetCheckZResolutionSPD();//A Khuntia
     //cutVertex->SetCheckDispersionSPD();//A Khuntia
     //cutVertex->SetCheckZDifferenceSPDTrack();//A Khuntia
   }
  

  
   ///////----------AKhuntia----------//////
   /*AliRsnCutEventUtils* cutEventUtils=0;
   cutEventUtils=new AliRsnCutEventUtils("cutEventUtils",kTRUE,rejectPileUp);
   cutEventUtils->SetCheckIncompleteDAQ();
   cutEventUtils->SetCheckSPDClusterVsTrackletBG();*/
   //------------------------------------
   // define and fill cut set for event cut

   AliRsnCutSet* eventCuts=0;
   if(cutVertex){
   eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutVertex);
   eventCuts->SetCutScheme(Form("%s", cutVertex->GetName()));
   task->SetEventCuts(eventCuts);
   }
   //
   // -- EVENT-ONLY COMPUTATIONS -------------------------------------------------------------------
   //   
   //vertex
   Int_t vtxID = task->CreateValue(AliRsnMiniValue::kVz, kFALSE);
   AliRsnMiniOutput *outVtx = task->CreateOutput("eventVtx", "HIST", "EVENT");
   outVtx->AddAxis(vtxID, 240, -12.0, 12.0);
   
   //multiplicity or centrality
   Int_t multID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
   AliRsnMiniOutput *outMult = task->CreateOutput("eventMult", "HIST", "EVENT");
   if (isPP) 
     outMult->AddAxis(multID, 400, 0.0, 400.0);
   else
     outMult->AddAxis(multID, 100, 0.0, 100.0);
   
   TH2F* hvz=new TH2F("hVzVsCent","", 100, 0., 100., 240, -12.0, 12.0);
   task->SetEventQAHist("vz",hvz);//plugs this histogram into the fHAEventVz data member

   TH2F* hmc=new TH2F("MultiVsCent","", 100, 0., 100., 400, 0., 400.);
   hmc->GetYaxis()->SetTitle("QUALITY");
   task->SetEventQAHist("multicent",hmc);//plugs this histogram into the fHAEventMultiCent data member

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
   //   
     gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigKStarPP8TeV_PID.C");
   //gROOT->LoadMacro("ConfigKStarPP8TeV_PID.C");
   if (!ConfigKStarPP8TeV_PID(task, isMC, isPP, "", cutsPair, aodFilterBit, customQualityCutsID, cutPiCandidate, cutKaCandidate, nsigmaPi, nsigmaKa, enableMonitor, isMC&IsMcTrueOnly,  monitorOpt.Data(), useMixLS, isMC&checkReflex, yaxisvar)) return 0x0;
   
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddAnalysisTaskTOFKStar - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}
