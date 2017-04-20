/***************************************************************************
//            Modified by Pragati Sahoo - 20/4/2017
//            Modified by Enrico Fragiacomo - 15/01/2014
//            Based on AddAnalysisTaskRsnMini
//
// Macro to configure the KStarPlusMinus analysis task 
// It calls all configs desired by the user, by means
// of the boolean switches defined in the first lines.
// ---
// Inputs:
//  1) flag to know if running on MC or data
//  2) collision system pp Run2 Analysis 
// --
// Returns:
//  kTRUE  --> initialization successful
//  kFALSE --> initialization failed (some config gave errors)
//
****************************************************************************/

enum eventMixConfig { kDisabled = -1,
		      kMixDefault,     //=0 //10 events, Dvz = 1cm, DC = 10
		      k5Evts,          //=1 //5 events, Dvz = 1cm, DC = 10
		      k5Cent,          //=2 //10 events, Dvz = 1cm, DC = 5
};


AliRsnMiniAnalysisTask *AddTaskKStarPlusMinus8TeVpp
(
 Bool_t      isMC,
 Bool_t      isPP,
 UInt_t      triggerMask=AliVEvent::kINT7,
 Float_t     cutV = 10.0,
 Int_t       mixingConfigID = 0,
 Int_t       aodFilterBit = 5,
 Bool_t      enableSys = kFALSE,
 Int_t       Sys = 0,
 Bool_t      enableMonitor=kTRUE,
 TString     monitorOpt="NoSIGN",
 Float_t     piPIDCut = 3.0,
 Float_t     pi_k0s_PIDCut = 5.0,
 Float_t     massTol = 0.03,
 Float_t     massTolVeto = 0.006,//here
 Float_t     pLife = 20,  
 Float_t     radiuslow = 0.5,
 Float_t     radiushigh = 200,    
 Float_t     MinDCAXY = 0.15,
 Bool_t      Switch = kFALSE,//here
 Float_t     k0sDCA = 0.3,
 Float_t     k0sCosPoinAn = 0.97,
 Float_t     k0sDaughDCA = 1.0,
 Int_t       NTPCcluster = 70,
 Float_t     maxDiffVzMix = 1.0,
 Float_t     maxDiffMultMix = 10.0,
 Float_t     maxDiffAngleMixDeg = 20.0,
 Int_t       aodN = 68,
 TString     outNameSuffix = "KStarPlusMinus",
 Int_t       centr = 0
 )
{
  //
  //-------------------------------------------                                                                                          
  // event cuts                                                                                                                             
  //-------------------------------------------                                                                                               
  //  UInt_t      triggerMask=AliVEvent::kINT7;
  Bool_t      rejectPileUp=kTRUE;
  Double_t    vtxZcut=10.0;//cm, default cut on vtx z                                                   
  if(isMC) rejectPileUp=kFALSE;
  
  //-------------------------------------------
  //mixing settings
  //-------------------------------------------

  Int_t       nmix = 10;
  if (mixingConfigID == eventMixConfig::kMixDefault) nmix = 10;  
  if (mixingConfigID == eventMixConfig::k5Evts)      nmix = 5;  
  if (mixingConfigID == eventMixConfig::k5Cent)      maxDiffMultMix = 5;
  
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskKStarPlusMinus8TeVpp", "No analysis manager to connect to.");
      return NULL;
   } 
   
   // create the task and configure 
   TString taskName = Form("KStarPlusMinus%s%s_%.1f_%d_%.1f_%.1f_%.2f_%.4f_%.2f_%.2f_%.1f", (isPP? "pp" : "PbPb"), (isMC ? "MC" : "Data"),cutV,NTPCcluster,piPIDCut,pi_k0s_PIDCut,massTol,k0sDCA,k0sCosPoinAn,k0sDaughDCA);
   //TString taskName=Form("TOFKstar%s%s_%i%i",(isPP? "pp" : "PbPb"),(isMC ? "MC" : "Data"),(Int_t)cutKaCandidate);
   AliRsnMiniAnalysisTask* task = new AliRsnMiniAnalysisTask(taskName.Data(),isMC);                   
   task->UseESDTriggerMask(triggerMask);

   //if(isPP) 
   task->UseMultiplicity("QUALITY");
   //else task->UseCentrality("V0M");

   // set event mixing options                                                                                                               
   task->UseContinuousMix();
   //task->UseBinnedMix();                                                                                                                   
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   ::Info("AddAnalysisTaskTOFKStar", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %\5.3f", nmix, maxDiffVzMix, maxDiffMultMix));
   
   mgr->AddTask(task);
     
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //  
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted

   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", cutV, 0, kFALSE, kTRUE);
   
   if(isPP && (!isMC)){ //assume pp data
     cutVertex->SetCheckPileUp(rejectPileUp);// set the check for pileup                                                                  
     ::Info("AddAnalysisTaskTOFKStar", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp)?"ON":"OFF"));
   }
   
   
   // define and fill cut set for event cut                                                                                          
   AliRsnCutSet* eventCuts=new AliRsnCutSet("eventCuts",AliRsnTarget::kEvent);
   //eventCuts->AddCut(cutEventUtils);
   eventCuts->AddCut(cutVertex);
   eventCuts->SetCutScheme(Form("%s",cutVertex->GetName()));
   task->SetEventCuts(eventCuts);

   // -- EVENT-ONLY COMPUTATIONS -------------------------------------------------------------------                                       
   //vertex                                                                                                                                
   Int_t vtxID=task->CreateValue(AliRsnMiniValue::kVz,kFALSE);
   AliRsnMiniOutput* outVtx=task->CreateOutput("eventVtx","HIST","EVENT");
   outVtx->AddAxis(vtxID,240,-12.0,12.0);

   //multiplicity
   Int_t multID=task->CreateValue(AliRsnMiniValue::kMult,kFALSE);
   AliRsnMiniOutput* outMult=task->CreateOutput("eventMult","HIST","EVENT");
   //if(isPP) 
   outMult->AddAxis(multID,400,0.5,400.5);
   //else outMult->AddAxis(multID,100,0.,100.);
   
   TH2F* hvz=new TH2F("hVzVsCent","",100,0.,100., 240,-12.0,12.0);
   task->SetEventQAHist("vz",hvz);//plugs this histogram into the fHAEventVz data member                                                    

   TH2F* hmc=new TH2F("MultiVsCent","", 100,0.,100., 400,0.5,400.5);
   hmc->GetYaxis()->SetTitle("QUALITY");
   task->SetEventQAHist("multicent",hmc);//plugs this histogram into the fHAEventMultiCent data member   
   
   //
   // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
   Double_t    minYlab =  -0.5;
   Double_t    maxYlab =  0.5;
   
   AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
   cutY->SetRangeD(minYlab, maxYlab);
   
   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());
   
   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigKStarPlusMinus8TeVpp.C");
   //gROOT->LoadMacro("ConfigKStarPlusMinus8TeVpp.C");
   if (isMC) {
     Printf("========================== MC analysis - PID cuts not used"); 
   } else 
     Printf("========================== DATA analysis - PID cuts used");
   
   if (!ConfigKStarPlusMinus8TeVpp(task, isPP, isMC, piPIDCut, pi_k0s_PIDCut, aodFilterBit,enableSys,Sys,enableMonitor,monitorOpt.Data(),massTol, massTolVeto, pLife, radiuslow, radiushigh,MinDCAXY, Switch, k0sDCA, k0sCosPoinAn, k0sDaughDCA, NTPCcluster, "", cutsPair)) return 0x0;
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddTaskKStarPlusMinus8TeVpp - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s", outNameSuffix.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}

