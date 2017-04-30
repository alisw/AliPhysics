/***************************************************************************
//            Modified by Kishora Nayak - 14/06/2016
//            Modified by Enrico Fragiacomo - 15/01/2014
//            Modified by Kunal Garg - 04/02/2017  
//            Based on AddAnalysisTaskRsnMini
//            pPb specific settings from AddTaskKStarPPB.C
//
// Macro to configure the KStarPlusMinus analysis task 
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

//enum ERsnCollType_t { kPP=0,      kPPb,      kPbPb};

enum pairYCutSet { kPairDefault,    // USED ONLY FOR pA
		   kNegative,       // USED ONLY FOR pA
		   kCentral         // USED ONLY FOR pA
                 };

/*enum eventCutSet { kOld = -1, 
		   kEvtDefault, //=0
		   kNoPileUpCut, //=1
		   kPileUpMV, //=2
		   kPileUpSPD3, //=3		      
		   kDefaultVtx8, //=4
		   kDefaultVtx5 //=5                    
};*/

enum eventCutSet { kEvtDefault=0,
                   kNoPileUpCut, //=1                                                                                                         
                   kDefaultVtx12,//=2                                                                                                         
                   kDefaultVtx8, //=3                                                                                                         
                   kDefaultVtx5, //=4                                                                                                         
                   kMCEvtDefault //=5                                                                                                         
};


enum eventMixConfig { kDisabled = -1,
		      kMixDefault,     //=0 //10 events, Dvz = 1cm, DC = 10
		      k5Evts,          //=1 //5 events, Dvz = 1cm, DC = 10
		      k5Cent,          //=2 //10 events, Dvz = 1cm, DC = 5
};


Ali:RsnMiniAnalysisTask *AddTaskKStarPlusMinusRun2
(
 Bool_t      isMC,
 Bool_t      isPP,
 // Int_t     collSyst,
 Float_t     cutV = 10.0,
 Int_t       evtCutSetID = 0,
 Int_t       pairCutSetID = 0,
 Int_t       mixingConfigID = 0,
 Int_t       aodFilterBit = 5,
 Bool_t      enableMonitor=kTRUE,
 TString     monitorOpt="pp", 
 Float_t     piPIDCut = 3.0,
 AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPiCandidate = AliRsnCutSetDaughterParticle::kTPCpidphipp2015,    
 Float_t     pi_k0s_PIDCut = 5.0,
 Float_t     massTol = 0.03,
 Float_t     massTolVeto = 0.004,
 Float_t     pLife = 20,  
 Float_t     radiuslow = 0.5,
 Float_t     radiushigh = 200,    
 Bool_t      Switch = kFALSE,
 Float_t     k0sDCA = 0.3,
 Float_t     k0sCosPoinAn = 0.97,
 Float_t     k0sDaughDCA = 1.0,
 Int_t       NTPCcluster = 70,
 Float_t     maxDiffVzMix = 1.0,
 Float_t     maxDiffMultMix = 10.0,
 Float_t     maxDiffAngleMixDeg = 20.0,
 Int_t       aodN = 68,
 TString     outNameSuffix = "KStarPlusMinus_TestPID",
 Int_t       centr = 0,
 Bool_t      ptDep = kTRUE,
 Float_t     DCAxy = 0.06,
 Bool_t      enableSys = kFALSE,
 Int_t       Sys= 0
 )
{  
  //-------------------------------------------                                                                                          
  // event cuts                                                                                                                             
  //-------------------------------------------                                                                                               
  UInt_t      triggerMask=AliVEvent::kINT7;
  Bool_t      rejectPileUp=kTRUE;
  Double_t    vtxZcut=10.0;//cm, default cut on vtx z                                                                                       
  //   cout<<"EVENTCUTID is    "<<evtCutSetID<<endl;                                                                                        
  if(evtCutSetID==eventCutSet::kDefaultVtx12) vtxZcut=12.0; //cm                                                                              
  if(evtCutSetID==eventCutSet::kDefaultVtx8) vtxZcut=8.0; //cm                                                                                
  if(evtCutSetID==eventCutSet::kDefaultVtx5) vtxZcut=5.0; //cm                                                                                
  if(evtCutSetID==eventCutSet::kNoPileUpCut) rejectPileUp=kFALSE;


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
      ::Error("AddTaskKStarPlusMinus", "No analysis manager to connect to.");
      return NULL;
   } 
   
   // create the task and configure 
    TString taskName = Form("KStarPlusMinus%s%s", (isPP? "pp" : "PbPb"), (isMC ? "MC" : "Data"));
 
   AliRsnMiniAnalysisTask* task = new AliRsnMiniAnalysisTask(taskName.Data(),isMC);
   
   //task->UseESDTriggerMask(AliVEvent::kINT7); //ESD ****** check this *****                    
   task->SelectCollisionCandidates(triggerMask); //AOD                                                                                        

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

   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", cutV, 0, kFALSE);
   cutVertex->SetCheckZResolutionSPD();
   cutVertex->SetCheckDispersionSPD(); 
   cutVertex->SetCheckZDifferenceSPDTrack();
   
   AliRsnCutEventUtils* cutEventUtils=new AliRsnCutEventUtils("cutEventUtils",kTRUE,rejectPileUp);
   cutEventUtils->SetCheckIncompleteDAQ();
   cutEventUtils->SetCheckSPDClusterVsTrackletBG();
   
   if(!isMC){ //assume pp data
     cutVertex->SetCheckPileUp(rejectPileUp);// set the check for pileup                                                                  
     ::Info("AddAnalysisTaskTOFKStar", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp)?"ON":"OFF"));
   }
   
   
   // define and fill cut set for event cut                                                                                          
   AliRsnCutSet* eventCuts=new AliRsnCutSet("eventCuts",AliRsnTarget::kEvent);
   eventCuts->AddCut(cutEventUtils);
   eventCuts->AddCut(cutVertex);
   eventCuts->SetCutScheme(Form("%s&%s",cutEventUtils->GetName(),cutVertex->GetName()));
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
   if (ptDep) {
     cutsPair->SetCutScheme(cutY->GetName()); 
   } else {
     AliRsnCutMiniPair *cutV0 = new AliRsnCutMiniPair("cutV0", AliRsnCutMiniPair::kContainsV0Daughter);
     cutsPair->AddCut(cutV0);
     cutsPair->SetCutScheme(TString::Format("%s&!%s",cutY->GetName(),cutV0->GetName()).Data());
   }
   
   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigKStarPlusMinusRun2.C");
   //gROOT->LoadMacro("ConfigKStarPlusMinusRun2.C");
   if (isMC) {
     Printf("========================== MC analysis - PID cuts not used"); 
   } else 
     Printf("========================== DATA analysis - PID cuts used");
   
   if (!ConfigKStarPlusMinusRun2(task, isPP, isMC, piPIDCut, cutPiCandidate, pi_k0s_PIDCut, aodFilterBit, enableMonitor, monitorOpt.Data(), massTol, massTolVeto, pLife, radiuslow, radiushigh, Switch, k0sDCA, k0sCosPoinAn, k0sDaughDCA, NTPCcluster, "", cutsPair, ptDep, DCAxy, enableSys, Sys)) return 0x0;
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddTaskKStarPlusMinus - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s_%.1f_%.1f_%.2f_%.3f_%.f_%.f_%.f_%.1f_%.2f_%.1f_%.3f_%.1f", outNameSuffix.Data(),piPIDCut,pi_k0s_PIDCut,massTol,massTolVeto,pLife,radiuslow,radiushigh,k0sDCA,k0sCosPoinAn,k0sDaughDCA, DCAxy, Sys), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}

