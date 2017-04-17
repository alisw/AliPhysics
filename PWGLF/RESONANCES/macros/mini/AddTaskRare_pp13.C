/***************************************************************************
              Anders Knospe: anders.knospe@cern.ch
                  last modified on 24/2/2017
  Macro to configure the resonance package for searches for rare resonances.

****************************************************************************/

AliRsnMiniAnalysisTask* AddTaskRare_pp13(TString lname, Bool_t isMC, Int_t system, AliRsnDaughter::ESpecies d1, AliRsnDaughter::ESpecies d2, Int_t EventCuts=0, Int_t TrackCuts1=0, Int_t TrackCuts2=0){
  if(d1==AliRsnDaughter::kPion && d2==AliRsnDaughter::kKaon){
    return AddTask_pikx(lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kPion && d1==AliRsnDaughter::kKaon){
    return AddTask_pikx(lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kPion && d2==AliRsnDaughter::kKaon0){
    return AddTask_pik0(lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kPion && d1==AliRsnDaughter::kKaon0){
    return AddTask_pik0(lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kKaon && d2==AliRsnDaughter::kKaon0){
    return AddTask_kxk0(lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kKaon && d1==AliRsnDaughter::kKaon0){
    return AddTask_kxk0(lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kProton && d2==AliRsnDaughter::kKaon){
    return AddTask_pkx(lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kProton && d1==AliRsnDaughter::kKaon){
    return AddTask_pkx(lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kProton && d2==AliRsnDaughter::kKaon0){
    return AddTask_pk0(lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kProton && d1==AliRsnDaughter::kKaon0){
    return AddTask_pk0(lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kLambda && d2==AliRsnDaughter::kPion){
    return AddTask_Lambdapi(lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kLambda && d1==AliRsnDaughter::kPion){
    return AddTask_Lambdapi(lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kLambda && d2==AliRsnDaughter::kKaon){
    return AddTask_Lambdakx(lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kLambda && d1==AliRsnDaughter::kKaon){
    return AddTask_Lambdakx(lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kLambda && d2==AliRsnDaughter::kKaon0){
    return AddTask_Lambdak0(lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kLambda && d1==AliRsnDaughter::kKaon0){
    return AddTask_Lambdak0(lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kLambda && d2==AliRsnDaughter::kProton){
    return AddTask_Lambdap(lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kLambda && d1==AliRsnDaughter::kProton){
    return AddTask_Lambdap(lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
  }else return 0;
}


//=============================


AliRsnMiniAnalysisTask * AddTask_pikx
(
 TString     lname = "pikx",
 Bool_t      isMC = kFALSE,
 Int_t       system = 0,
 Int_t       EventCuts = 0,
 Int_t       TrackCutsPi = 0,
 Int_t       TrackCutsK = 0
 )
{
  Bool_t      isPP=(system==0)?true:false;
  Bool_t      rejectPileUp = kTRUE;
  Int_t       nmix = 5;
  Float_t     nsigmaPi = 3002.0;
  Float_t     nsigmaK = 3002.0;
  Float_t     vertexz = 10.0;
  Bool_t      enableMonitor = kTRUE;
  Float_t     maxDiffVzMix = 1.0;
  Float_t     maxDiffMultMix = 5.0;
  TString     outNameSuffix = lname;
  
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskRare_pp13::AddTask_pikx", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   
   TString taskName = Form("%s_%s%s", lname.Data(), (isPP? "pp" : "pPb"), (isMC ? "MC" : "Data"));
   
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   //task->SelectCollisionCandidates(AliVEvent::kMB);
   int trigger=EventCuts%10;
   if(!trigger) task->UseESDTriggerMask(AliVEvent::kINT7);
   else if(trigger==1) task->UseESDTriggerMask(AliVEvent::kHighMult);
   else if(trigger==2) task->UseESDTriggerMask(AliVEvent::kHighMultV0);
   else task->UseESDTriggerMask(AliVEvent::kAny);

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
   task->UseMC(isMC);
   ::Info("AddTaskRare_pp13::AddTask_pikx", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n", nmix, maxDiffVzMix, maxDiffMultMix));
   
   mgr->AddTask(task);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //  
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
   
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", vertexz, 0, kFALSE);
   //    if (isPP) cutVertex->SetCheckPileUp(kTRUE); // set the check for pileup// check in AliRsnCutEventUtils
   cutVertex->SetCheckZResolutionSPD();
   cutVertex->SetCheckDispersionSPD();
   cutVertex->SetCheckZDifferenceSPDTrack();
   
   AliRsnCutEventUtils* cutEventUtils=new AliRsnCutEventUtils("cutEventUtils",kTRUE,rejectPileUp);
   cutEventUtils->SetCheckIncompleteDAQ();
   cutEventUtils->SetCheckSPDClusterVsTrackletBG();
   
   
   // define and fill cut set for event cut
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutVertex);
   eventCuts->SetCutScheme(Form("%s",cutVertex->GetName()));
      
 
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
   if (isPP) outMult->AddAxis(multID, 400, 0.0, 400.0);
   else outMult->AddAxis(multID, 100, 0.0, 100.0);
   
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

   if (!Config_pikx(task,isMC,isPP,"",cutsPair,TrackCutsPi,TrackCutsK,nsigmaPi,nsigmaK,enableMonitor)) return 0x0;

   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddTaskRare_pp13::AddTask_pikx - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}

Bool_t Config_pikx(  
			AliRsnMiniAnalysisTask *task, 
			Bool_t                 isMC, 
			Bool_t                 isPP,
			const char             *suffix,
			AliRsnCutSet           *cutsPair,
			Int_t                  TrackCutsPi = 0,
			Int_t                  TrackCutsK = 0,
			Float_t                nsigmaPi = 3.0,
			Float_t                nsigmaK  = 3.0,
			Bool_t                 enableMonitor = kTRUE
		     )
{ 
  //These are the Default values for 2011 ESD track cuts
  
  AliPID::EParticleType  type1   = AliPID::kPion;
  AliPID::EParticleType  type2   = AliPID::kKaon;

  // manage suffix
  if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
  
  // retrieve mass from PDG database
  Int_t         pdg  = 313;
  TDatabasePDG *db   = TDatabasePDG::Instance();
  TParticlePDG *part = db->GetParticle(pdg);
  Double_t mass      = part->Mass();
  
  // set daughter cuts

  Float_t nsigmaPiTPC=fmod(nsigmaPi,1000.);
  Float_t nsigmaPiTOF=(nsigmaPi-fmod(nsigmaPi,1000.))/1000.;
  Float_t nsigmaKTPC=fmod(nsigmaK,1000.);
  Float_t nsigmaKTOF=(nsigmaK-fmod(nsigmaK,1000.))/1000.;

  AliRsnCutTrackQuality* trkQualityCut= new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  AliRsnCutSetDaughterParticle* cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaPi),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kPion,nsigmaPiTPC,nsigmaPiTOF);
  AliRsnCutSetDaughterParticle* cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015, nsigmaK),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kKaon,nsigmaKTPC,nsigmaKTOF);
  
  Int_t iCutPi = task->AddTrackCuts(cutSetPi);
  Int_t iCutK  = task->AddTrackCuts(cutSetK);
  
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetPi->GetMonitorOutput());
    AddMonitorOutput(isMC, cutSetK->GetMonitorOutput());
  }  
  // -- Values ------------------------------------------------------------------------------------
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass,    kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);
  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges
  // [0] = unlike
  // [1] = mixing
  // [2] = like ++
  // [3] = like --

  Bool_t  use     [12] = {1         ,1         ,1         ,1         ,1       ,1       ,isMC     ,isMC     ,isMC     ,isMC     ,isMC    ,isMC    };
  Bool_t  useIM   [12] = {1         ,1         ,1         ,1         ,1       ,1       ,1        ,1        ,1        ,1        ,0       ,0       };
  TString name    [12] = {"UnlikePM","UnlikeMP","MixingPM","MixingMP","LikePP","LikeMM","MCGenPM","MCGenMP","TruesPM","TruesMP","ResPM" ,"ResMP" };
  TString comp    [12] = {"PAIR"    ,"PAIR"    ,"MIX"     ,"MIX"     ,"PAIR"  ,"PAIR"  ,"MOTHER" ,"MOTHER" ,"TRUE"   ,"TRUE"   ,"TRUE"  ,"TRUE"  };
  TString output  [12] = {"SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE"  ,"SPARSE","SPARSE","SPARSE" ,"SPARSE" ,"SPARSE" ,"SPARSE" ,"SPARSE","SPARSE"};
  Char_t  charge1 [12] = {'+'       ,'-'       ,'+'       ,'-'       ,'+'     ,'-'     ,'+'      ,'-'      ,'+'      ,'-'      ,'+'     ,'-'     };
  Char_t  charge2 [12] = {'-'       ,'+'       ,'-'       ,'+'       ,'+'     ,'-'     ,'-'      ,'+'      ,'_'      ,'+'      ,'-'     ,'+'     };
  Int_t   cutIDPi [12] = {iCutPi    ,iCutPi    ,iCutPi    ,iCutPi    ,iCutPi  ,iCutPi  ,iCutPi   ,iCutPi   ,iCutPi   ,iCutPi   ,iCutPi  ,iCutPi  };
  Int_t   cutIDK  [12] = {iCutK     ,iCutK     ,iCutK     ,iCutK     ,iCutK   ,iCutK   ,iCutK    ,iCutK    ,iCutK    ,iCutK    ,iCutK   ,iCutK   };
  Int_t   PDGCode [12] = {313       ,313       ,313       ,313       ,313     ,313     ,313      ,-313     ,313      ,313      ,313     ,-313    };

  for (Int_t i = 0; i < 12; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("pikx_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    out->SetDaughter(0, AliRsnDaughter::kKaon);
    out->SetDaughter(1, AliRsnDaughter::kPion);
    out->SetCutID(0, cutIDK[i]);
    out->SetCutID(1, cutIDPi[i]);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(PDGCode[i]);
    out->SetMotherMass(mass);
    out->SetPairCuts(cutsPair);

    // axis X: invmass (or resolution)
    if (useIM[i]) 
      out->AddAxis(imID, 1370, 0.63, 2.);
    else
      out->AddAxis(resID, 200, -0.02, 0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID, 200, 0.0, 20.0);
    
    // axis Z: centrality-multiplicity
    /*
    if (!isPP)
      out->AddAxis(centID, 100, 0.0, 100.0);
    else 
      out->AddAxis(centID, 400, 0.0, 400.0);
    */
      
    // axis W: pseudorapidity
    // out->AddAxis(etaID, 20, -1.0, 1.0);
    // axis J: rapidity
    //out->AddAxis(yID, 90, -4.5, 4.5);
    
  }
  return kTRUE;
}


//=============================


enum pairYCutSet { kPairDefault,    // USED ONLY FOR pA
		   kNegative,       // USED ONLY FOR pA
		   kCentral         // USED ONLY FOR pA
                 };

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

AliRsnMiniAnalysisTask *AddTask_pik0
(
 TString     lname,
 Bool_t      isMC,
 Int_t       system,
 Int_t       EventCuts,
 Int_t       TrackCutsPi,
 Int_t       TrackCutsK
 )
{
  Float_t     cutV = 10.0;
  Int_t       evtCutSetID = 0;
  Int_t       pairCutSetID = 0;
  Int_t       mixingConfigID = 1;
  Int_t       aodFilterBit = 5;
  Bool_t      enableMonitor=kTRUE;
  TString     monitorOpt="pp";
  Float_t     piPIDCut = 3.0;
  Float_t     pi_k0s_PIDCut = 5.0;
  //Float_t     trackDCAcut = 7.0;
  Float_t     massTol = 0.03;
  Float_t     massTolVeto = 0.004;
  Float_t     pLife = 20;
  Float_t     radiuslow = 0.5;
  Float_t     radiushigh = 200;   
  Bool_t      Switch = kFALSE;
  Float_t     k0sDCA = 0.3;
  Float_t     k0sCosPoinAn = 0.97;
  Float_t     k0sDaughDCA = 1.0;
  Int_t       NTPCcluster = 70;
  Float_t     maxDiffVzMix = 1.0;
  Float_t     maxDiffMultMix = 10.0;
  Float_t     maxDiffAngleMixDeg = 20.0;
  Int_t       aodN = 68;
  TString     outNameSuffix = lname;
  Int_t       centr = 0;
  Bool_t      ptDep = kTRUE;
   
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
      ::Error("AddTaskRare_pp13::AddTask_pik0", "No analysis manager to connect to.");
      return NULL;
   } 
   
   // create the task and configure 
   TString taskName = Form("%s_%s%s", lname.Data(), ((system==0)? "pp" : "PbPb"), (isMC ? "MC" : "Data"));
 
   AliRsnMiniAnalysisTask* task = new AliRsnMiniAnalysisTask(taskName.Data(),isMC);
   
   //task->SelectCollisionCandidates(triggerMask); //AOD
   int trigger=EventCuts%10;
   if(!trigger) task->UseESDTriggerMask(AliVEvent::kINT7);
   else if(trigger==1) task->UseESDTriggerMask(AliVEvent::kHighMult);
   else if(trigger==2) task->UseESDTriggerMask(AliVEvent::kHighMultV0);
   else task->UseESDTriggerMask(AliVEvent::kAny);

   //if(isPP) 
   task->UseMultiplicity("QUALITY");
   //else task->UseCentrality("V0M");

   // set event mixing options                                                                                                               
   task->UseContinuousMix();
   //task->UseBinnedMix();                                                                                                                   
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   ::Info("AddTaskRare_pp13::AddTask_pik0", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %\5.3f", nmix, maxDiffVzMix, maxDiffMultMix));
   
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
     ::Info("AddTaskRare_pp13::AddTask_pik0", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp)?"ON":"OFF"));
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
   cutsPair->SetCutScheme(cutY->GetName());
   
   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------

   Bool_t isPP=kTRUE;
   if (isMC) {
     Printf("========================== MC analysis - PID cuts not used"); 
   } else 
     Printf("========================== DATA analysis - PID cuts used");
   
   if (!Config_pik0(task, isPP, isMC, TrackCutsPi, TrackCutsK, piPIDCut, pi_k0s_PIDCut, aodFilterBit, enableMonitor, monitorOpt.Data(), massTol, massTolVeto, pLife, radiuslow, radiushigh, Switch, k0sDCA, k0sCosPoinAn, k0sDaughDCA, NTPCcluster, "", cutsPair, ptDep)) return 0x0;
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddTaskRare_pp13::AddTask_pik0 - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s", outNameSuffix.Data(),), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}


//
// *** Configuration script for K*+-->K0Short-Pi analysis ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
//
Bool_t Config_pik0
(  
   AliRsnMiniAnalysisTask *task,
   Bool_t                  isPP,
   //Int_t		           collSyst, 
   Bool_t                  isMC,
   Int_t                   TrackCutsPi,
   Int_t                   TrackCutsK,
   Float_t                 piPIDCut,
   Float_t                 pi_k0s_PIDCut,
   Int_t                   aodFilterBit,
   //Float_t               trackDCAcut,
   Bool_t                  enableMonitor=kTRUE  ,
   TString                 monitorOpt="",
   Float_t                 massTol,
   Float_t                 massTolVeto, 
   Float_t                 pLife, 
   Float_t                 radiuslow,
   Float_t                 radiushigh,    
   Bool_t                  Switch,     
   Float_t                 k0sDCA,
   Float_t                 k0sCosPoinAn,
   Float_t                 k0sDaughDCA,
   Int_t                   NTPCcluster,
   const char             *suffix,
   AliRsnCutSet           *cutsPair,
   Bool_t                  ptDep    
)
{
  Float_t nsigmaPiTOF=3.;
   // manage suffix
   if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
 
   /////////////////////////////////////////////////////
   // selections for the pion from the decay of KStarPlusMinus*
   /////////////////////////////////////////////////////
   //
   AliRsnCutTrackQuality* trkQualityCut= new AliRsnCutTrackQuality("myQualityCut");
   trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
   AliRsnCutSetDaughterParticle* cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,piPIDCut),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kPion,piPIDCut,nsigmaPiTOF);
   Int_t iCutPi = task->AddTrackCuts(cutSetPi);

   /*
   AliRsnCutDaughterSigmaStar2010PP *cutPi = new AliRsnCutDaughterSigmaStar2010PP("cutPionFork0pi", AliPID::kPion);
   cutPi->SetPIDCut(piPIDCut);    // fPIDCut used in IsSelected() after the call to cutQuality
   AliRsnCutTrackQuality *cutQuality = (AliRsnCutTrackQuality*) cutPi->CutQuality();
   cutQuality->SetDefaults2011();
   //cutQuality->SetDefaults2010(0,1);  // 1st par. not default (0 -> use TPC clusters). 2nd par. default (-> standard Pt and eta range)
     
   AliRsnCutSet *cutSetPi = new AliRsnCutSet("setPionFork0pi", AliRsnTarget::kDaughter);
   cutSetPi->AddCut(cutPi);
   cutSetPi->SetCutScheme(cutPi->GetName());
   Int_t iCutPi = task->AddTrackCuts(cutSetPi);
   */
   //
   /////////////////////////////////////////////////////////////
   // selections for K0s and for the daughters of K0s
   /////////////////////////////////////////////////////////////
   // 
   // selections for pion daugthers of K0s
   AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("qualityDaughterK0s");   
   esdTrackCuts->SetPtRange(0.15,1.E10);
   esdTrackCuts->SetEtaRange(-0.8,0.8);
   esdTrackCuts->SetRequireTPCRefit();
   esdTrackCuts->SetAcceptKinkDaughters(0); //
   esdTrackCuts->SetMinNClustersTPC(NTPCcluster);
   esdTrackCuts->SetMaxChi2PerClusterTPC(4);
    
    if(ptDep){
   esdTrackCuts->SetMinDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    }else
   esdTrackCuts->SetMinDCAToVertexXY("0.06"); //Use one of the two - pt dependent or fixed value cut.
  
   //
   /////////////////////////////////////////////////
   // selections for K0s
   AliRsnCutV0 *cutK0s = new AliRsnCutV0("cutK0s", kK0Short, AliPID::kPion, AliPID::kPion);
   //cutK0s->SetPIDCutProton(pPIDCut);       // PID for the proton daughter of K0s
   cutK0s->SetPIDCutPion(pi_k0s_PIDCut);        // PID for the pion daughter of K0s
   cutK0s->SetESDtrackCuts(esdTrackCuts);  // all the other selections (defined above) for proton and pion daughters of K0s
   cutK0s->SetMaxDaughtersDCA(k0sDaughDCA);
   cutK0s->SetMaxDCAVertex(k0sDCA);
   cutK0s->SetMinCosPointingAngle(k0sCosPoinAn);
   cutK0s->SetTolerance(massTol);
   cutK0s->SetToleranceVeto(massTolVeto);   //Rejection range for Competing V0 Rejection
   cutK0s->SetSwitch(Switch);    
   cutK0s->SetfLife(pLife); 
   cutK0s->SetfLowRadius(radiuslow); 
   cutK0s->SetfHighRadius(radiushigh); 
   cutK0s->SetMaxRapidity(2.0);
   //
   AliRsnCutSet *cutSetK0s = new AliRsnCutSet("setK0s", AliRsnTarget::kDaughter);
   cutSetK0s->AddCut(cutK0s);
   cutSetK0s->SetCutScheme(cutK0s->GetName());
   Int_t iCutK0s = task->AddTrackCuts(cutSetK0s);
   //
   
   //
    
  if(enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetPi->GetMonitorOutput(), monitorOpt.Data());
  } 
    
   //
   // -- Values ------------------------------------------------------------------------------------
   //
   
   /// -- Values ------------------------------------------------------------------------------------                                        
    /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
   /* transv. momentum */ Int_t ptID    = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
   /* centrality       */ Int_t centID  = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
   /* pseudorapidity   */ Int_t etaID   = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
   /* rapidity         */ Int_t yID     = task->CreateValue(AliRsnMiniValue::kY, kFALSE);

   
   //
   // -- Create all needed outputs -----------------------------------------------------------------
   //
   
   // use an array for more compact writing, which are different on mixing and charges
   // [0] = unlike
   // [1] = mixing
   // [2] = like ++
   // [3] = like --
   Bool_t  use     [6] = {1               ,1                ,1                  ,1                   ,1                ,1                 };
   Bool_t  useIM   [6] = {1               ,1                ,1                  ,1                   ,1                ,1                 };
   TString name    [6] = {"K0Pip"         ,"K0Pim"          ,"K0PipMix"         ,"K0PimMix"          ,"KStarPlusMinust","AKStarPlusMinust"};
   TString comp    [6] = {"PAIR"          ,"PAIR"           ,"MIX"              ,"MIX"               ,"TRUE"           ,"TRUE"            };
   TString output  [6] = {"SPARSE"        ,"SPARSE"         ,"SPARSE"           ,"SPARSE"            ,"SPARSE"         ,"SPARSE"          };
   Char_t  charge1 [6] = {'0'             ,'0'              ,'0'                ,'0'                 ,'0'              ,'0'               };
   Char_t  charge2 [6] = {'+'             ,'-'              ,'+'                ,'-'                 ,'+'              ,'-'               };
   Int_t   cutID1  [6] = { iCutK0s        ,iCutK0s           ,iCutK0s            ,iCutK0s            ,iCutK0s          ,iCutK0s           };
   Int_t   cutID2  [6] = { iCutPi         ,iCutPi           ,iCutPi             ,iCutPi              ,iCutPi           ,iCutPi            };
   Int_t   ipdg    [6] = {323             ,-323             ,323                ,-323                ,323              ,-323              };
   Double_t mass   [6] = { 0.89166        ,0.89166          ,0.89166            ,0.89166             ,0.89166          ,0.89166           };
   
   for (Int_t i = 0; i < 6; i++) {
     if (!use[i]) continue;
     //if (collSyst) output[i] = "SPARSE";
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("pik0_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
     // selection settings
     out->SetCutID(0, cutID1[i]);
     out->SetCutID(1, cutID2[i]);
     out->SetDaughter(0, AliRsnDaughter::kKaon0);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetCharge(0, charge1[i]);
     out->SetCharge(1, charge2[i]);
     out->SetMotherPDG(ipdg[i]);
     out->SetMotherMass(mass[i]);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // axis X: invmass
     if (useIM[i]) 
       out->AddAxis(imID, 1370, 0.63, 2.);
     // axis Y: transverse momentum
     out->AddAxis(ptID, 300, 0.0, 30.0);
   //  out->AddAxis(k0sDCA, 10, 0.0, 1.0);
     
     //if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     //if(isPP) out->AddAxis(centID, 400, 0.5, 400.5);
     //else out->AddAxis(centID, 100, 0.0, 100.);
   } 

   TString pname="pi";
   AddMonitorOutput_Pt(pname,cutSetPi->GetMonitorOutput());
   AddMonitorOutput_Eta(pname,cutSetPi->GetMonitorOutput());
   AddMonitorOutput_DCAxy(pname,cutSetPi->GetMonitorOutput());
   AddMonitorOutput_DCAz(pname,cutSetPi->GetMonitorOutput());
   AddMonitorOutput_TPCpi(pname,cutSetPi->GetMonitorOutput());
   AddMonitorOutput_NclTPC(pname,cutSetPi->GetMonitorOutput());
   AddMonitorOutput_chi2TPC(pname,cutSetPi->GetMonitorOutput());

   pname.Form("k0");
   AddMonitorOutput_P(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_Pt(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0NPt(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0PPt(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0Mass(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0DCA(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0Radius(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0Lifetime(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0DaughterDCA(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0CPA(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0DCA2TPV(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0TPCpim(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0TPCpip(pname,cutSetK0s->GetMonitorOutput());
   
   if (isMC) {
     
     TString mode = "SPARSE";
     //TString mode = "HIST";
     //if (collSyst) mode = "SPARSE";
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("KStarPlusMinus_MotherMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kKaon0);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetMotherPDG(323);
     out->SetMotherMass(0.89166);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 90, 0.6, 1.5);
     out->AddAxis(ptID, 300, 0.0, 30.0);
     //out->AddAxis(k0sDCA, 10, 0.0, 1.0);
     
     if(isPP) out->AddAxis(centID, 400, 0.5, 400.5);
     else out->AddAxis(centID, 100, 0.0, 100.);
     //if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("AKStarPlusMinus_MotherMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kKaon0);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetMotherPDG(-323);
     out->SetMotherMass(0.89166);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 90, 0.6, 1.5);
     out->AddAxis(ptID, 300, 0.0, 30.0);
     //out->AddAxis(k0sDCA, 10, 0.0, 1.0);
     //if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     if(isPP) out->AddAxis(centID, 400, 0.5, 400.5);
     else out->AddAxis(centID, 100, 0.0, 100.);
   }
   
   return kTRUE;
}


//=============================


AliRsnMiniAnalysisTask *AddTask_kxk0
(
 TString     lname,
 Bool_t      isMC,
 Int_t       system,
 Int_t       EventCuts,
 Int_t       TrackCutsKx,
 Int_t       TrackCutsK0
 )
{
  Float_t     cutV = 10.0;
  Int_t       evtCutSetID = 0;
  Int_t       pairCutSetID = 0;
  Int_t       mixingConfigID = 1;
  Int_t       aodFilterBit = 5;
  Bool_t      enableMonitor=kTRUE;
  TString     monitorOpt="pp";
  Float_t     kxPIDCut = 3.0;
  Float_t     pi_k0s_PIDCut = 5.0;
  //Float_t     trackDCAcut = 7.0;
  Float_t     massTol = 0.03;
  Float_t     massTolVeto = 0.004;
  Float_t     pLife = 20;
  Float_t     radiuslow = 0.5;
  Float_t     radiushigh = 200;   
  Bool_t      Switch = kFALSE;
  Float_t     k0sDCA = 0.3;
  Float_t     k0sCosPoinAn = 0.97;
  Float_t     k0sDaughDCA = 1.0;
  Int_t       NTPCcluster = 70;
  Float_t     maxDiffVzMix = 1.0;
  Float_t     maxDiffMultMix = 10.0;
  Float_t     maxDiffAngleMixDeg = 20.0;
  Int_t       aodN = 68;
  TString     outNameSuffix = lname;
  Int_t       centr = 0;
  Bool_t      ptDep = kTRUE;
   
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
      ::Error("AddTaskRare_pp13::AddTask_kxk0", "No analysis manager to connect to.");
      return NULL;
   } 
   
   // create the task and configure 
   TString taskName = Form("%s_%s%s", lname.Data(), ((system==0)? "pp" : "PbPb"), (isMC ? "MC" : "Data"));
 
   AliRsnMiniAnalysisTask* task = new AliRsnMiniAnalysisTask(taskName.Data(),isMC);

   //task->SelectCollisionCandidates(triggerMask); //AOD
   int trigger=EventCuts%10;
   if(!trigger) task->UseESDTriggerMask(AliVEvent::kINT7);
   else if(trigger==1) task->UseESDTriggerMask(AliVEvent::kHighMult);
   else if(trigger==2) task->UseESDTriggerMask(AliVEvent::kHighMultV0);
   else task->UseESDTriggerMask(AliVEvent::kAny);

   //if(isPP) 
   task->UseMultiplicity("QUALITY");
   //else task->UseCentrality("V0M");

   // set event mixing options                                                                                                               
   task->UseContinuousMix();
   //task->UseBinnedMix();                                                                                                                   
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   ::Info("AddTaskRare_pp13::AddTask_kxk0", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %\5.3f", nmix, maxDiffVzMix, maxDiffMultMix));
   
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
     ::Info("AddTaskRare_pp13::AddTask_kxk0", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp)?"ON":"OFF"));
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
   cutsPair->SetCutScheme(cutY->GetName());
   
   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------

   Bool_t isPP=kTRUE;
   if (isMC) {
     Printf("========================== MC analysis - PID cuts not used"); 
   } else 
     Printf("========================== DATA analysis - PID cuts used");
   
   if (!Config_kxk0(task, isPP, isMC, TrackCutsKx, TrackCutsK0, kxPIDCut, pi_k0s_PIDCut, aodFilterBit, enableMonitor, monitorOpt.Data(), massTol, massTolVeto, pLife, radiuslow, radiushigh, Switch, k0sDCA, k0sCosPoinAn, k0sDaughDCA, NTPCcluster, "", cutsPair, ptDep)) return 0x0;
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddTaskRare_pp13::AddTask_kxk0 - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s", outNameSuffix.Data(),), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}


Bool_t Config_kxk0
(  
   AliRsnMiniAnalysisTask *task,
   Bool_t                  isPP,
   //Int_t		           collSyst, 
   Bool_t                  isMC,
   Int_t                   TrackCutsKx,
   Int_t                   TrackCutsK0,
   Float_t                 kxPIDCut,
   Float_t                 pi_k0s_PIDCut,
   Int_t                   aodFilterBit,
   //Float_t               trackDCAcut,
   Bool_t                  enableMonitor=kTRUE  ,
   TString                 monitorOpt="",
   Float_t                 massTol,
   Float_t                 massTolVeto, 
   Float_t                 pLife, 
   Float_t                 radiuslow,
   Float_t                 radiushigh,    
   Bool_t                  Switch,     
   Float_t                 k0sDCA,
   Float_t                 k0sCosPoinAn,
   Float_t                 k0sDaughDCA,
   Int_t                   NTPCcluster,
   const char             *suffix,
   AliRsnCutSet           *cutsPair,
   Bool_t                  ptDep    
)
{
   Float_t nsigmaKxTOF=3.;

   // manage suffix
   if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
 
   /////////////////////////////////////////////////////
   // selections for the pion from the decay of KStarPlusMinus*
   /////////////////////////////////////////////////////
   //

   AliRsnCutTrackQuality* trkQualityCut= new AliRsnCutTrackQuality("myQualityCut");
   trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
   AliRsnCutSetDaughterParticle* cutSetKx=new AliRsnCutSetDaughterParticle(Form("cutKx%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,kxPIDCut),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kKaon,kxPIDCut,nsigmaKxTOF);
   Int_t iCutKx = task->AddTrackCuts(cutSetKx);

   /*
   AliRsnCutDaughterSigmaStar2010PP *cutKx = new AliRsnCutDaughterSigmaStar2010PP("cutKaonForkxk0", AliPID::kKaon);
   cutKx->SetPIDCut(kxPIDCut);    // fPIDCut used in IsSelected() after the call to cutQuality
   AliRsnCutTrackQuality *cutQuality = (AliRsnCutTrackQuality*) cutKx->CutQuality();
   cutQuality->SetDefaults2011();
   //cutQuality->SetDefaults2010(0,1);  // 1st par. not default (0 -> use TPC clusters). 2nd par. default (-> standard Pt and eta range)
     
   AliRsnCutSet *cutSetKx = new AliRsnCutSet("setKaonForkxk0", AliRsnTarget::kDaughter);
   cutSetKx->AddCut(cutKx);
   cutSetKx->SetCutScheme(cutKx->GetName());
   Int_t iCutKx = task->AddTrackCuts(cutSetKx);
   */
   //
   /////////////////////////////////////////////////////////////
   // selections for K0s and for the daughters of K0s
   /////////////////////////////////////////////////////////////
   // 
   // selections for pion daugthers of K0s
   AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("qualityDaughterK0s");   
   esdTrackCuts->SetPtRange(0.15,1.E10);
   esdTrackCuts->SetEtaRange(-0.8,0.8);
   esdTrackCuts->SetRequireTPCRefit();
   esdTrackCuts->SetAcceptKinkDaughters(0); //
   esdTrackCuts->SetMinNClustersTPC(NTPCcluster);
   esdTrackCuts->SetMaxChi2PerClusterTPC(4);
    
    if(ptDep){
   esdTrackCuts->SetMinDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    }else
   esdTrackCuts->SetMinDCAToVertexXY("0.06"); //Use one of the two - pt dependent or fixed value cut.
  
   //
   /////////////////////////////////////////////////
   // selections for K0s
   AliRsnCutV0 *cutK0s = new AliRsnCutV0("cutK0s", kK0Short, AliPID::kPion, AliPID::kPion);
   //cutK0s->SetPIDCutProton(pPIDCut);       // PID for the proton daughter of K0s
   cutK0s->SetPIDCutPion(pi_k0s_PIDCut);        // PID for the pion daughter of K0s
   cutK0s->SetESDtrackCuts(esdTrackCuts);  // all the other selections (defined above) for proton and pion daughters of K0s
   cutK0s->SetMaxDaughtersDCA(k0sDaughDCA);
   cutK0s->SetMaxDCAVertex(k0sDCA);
   cutK0s->SetMinCosPointingAngle(k0sCosPoinAn);
   cutK0s->SetTolerance(massTol);
   cutK0s->SetToleranceVeto(massTolVeto);   //Rejection range for Competing V0 Rejection
   cutK0s->SetSwitch(Switch);    
   cutK0s->SetfLife(pLife); 
   cutK0s->SetfLowRadius(radiuslow); 
   cutK0s->SetfHighRadius(radiushigh); 
   cutK0s->SetMaxRapidity(2.0);
   //
   AliRsnCutSet *cutSetK0s = new AliRsnCutSet("setK0s", AliRsnTarget::kDaughter);
   cutSetK0s->AddCut(cutK0s);
   cutSetK0s->SetCutScheme(cutK0s->GetName());
   Int_t iCutK0s = task->AddTrackCuts(cutSetK0s);
   //
   
   //
    
  if(enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetKx->GetMonitorOutput(), monitorOpt.Data());
  } 
    
   //
   // -- Values ------------------------------------------------------------------------------------
   //
   
   /// -- Values ------------------------------------------------------------------------------------                                        
    /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
   /* transv. momentum */ Int_t ptID    = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
   /* centrality       */ Int_t centID  = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
   /* pseudorapidity   */ Int_t etaID   = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
   /* rapidity         */ Int_t yID     = task->CreateValue(AliRsnMiniValue::kY, kFALSE);

   
   //
   // -- Create all needed outputs -----------------------------------------------------------------
   //
   
   // use an array for more compact writing, which are different on mixing and charges
   // [0] = unlike
   // [1] = mixing
   // [2] = like ++
   // [3] = like --
   Bool_t  use     [6] = {1               ,1                ,1                  ,1                   ,1                ,1                 };
   Bool_t  useIM   [6] = {1               ,1                ,1                  ,1                   ,1                ,1                 };
   TString name    [6] = {"K0Kp"          ,"K0Km"           ,"K0KpMix"          ,"K0KmMix"           ,"KStarPlusMinust","AKStarPlusMinust"};
   TString comp    [6] = {"PAIR"          ,"PAIR"           ,"MIX"              ,"MIX"               ,"TRUE"           ,"TRUE"            };
   TString output  [6] = {"SPARSE"        ,"SPARSE"         ,"SPARSE"           ,"SPARSE"            ,"SPARSE"         ,"SPARSE"          };
   Char_t  charge1 [6] = {'0'             ,'0'              ,'0'                ,'0'                 ,'0'              ,'0'               };
   Char_t  charge2 [6] = {'+'             ,'-'              ,'+'                ,'-'                 ,'+'              ,'-'               };
   Int_t   cutID1  [6] = { iCutK0s        ,iCutK0s           ,iCutK0s            ,iCutK0s            ,iCutK0s          ,iCutK0s           };
   Int_t   cutID2  [6] = { iCutKx         ,iCutKx           ,iCutKx             ,iCutKx              ,iCutKx           ,iCutKx            };
   Int_t   ipdg    [6] = {323             ,-323             ,323                ,-323                ,323              ,-323              };
   Double_t mass   [6] = { 0.89166        ,0.89166          ,0.89166            ,0.89166             ,0.89166          ,0.89166           };
   
   for (Int_t i = 0; i < 6; i++) {
     if (!use[i]) continue;
     //if (collSyst) output[i] = "SPARSE";
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("kxk0_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
     // selection settings
     out->SetCutID(0, cutID1[i]);
     out->SetCutID(1, cutID2[i]);
     out->SetDaughter(0, AliRsnDaughter::kKaon0);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetCharge(0, charge1[i]);
     out->SetCharge(1, charge2[i]);
     out->SetMotherPDG(ipdg[i]);
     out->SetMotherMass(mass[i]);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // axis X: invmass
     if (useIM[i]) 
       out->AddAxis(imID, 1005, 0.99, 3.);
     // axis Y: transverse momentum
     out->AddAxis(ptID, 300, 0.0, 30.0);
   //  out->AddAxis(k0sDCA, 10, 0.0, 1.0);
     
     //if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     //if(isPP) out->AddAxis(centID, 400, 0.5, 400.5);
     //else out->AddAxis(centID, 100, 0.0, 100.);
   } 

   TString pname="kx";
   AddMonitorOutput_Pt(pname,cutSetKx->GetMonitorOutput());
   AddMonitorOutput_Eta(pname,cutSetKx->GetMonitorOutput());
   AddMonitorOutput_DCAxy(pname,cutSetKx->GetMonitorOutput());
   AddMonitorOutput_DCAz(pname,cutSetKx->GetMonitorOutput());
   AddMonitorOutput_TPCpi(pname,cutSetKx->GetMonitorOutput());
   AddMonitorOutput_NclTPC(pname,cutSetKx->GetMonitorOutput());
   AddMonitorOutput_chi2TPC(pname,cutSetKx->GetMonitorOutput());

   pname.Form("k0");
   AddMonitorOutput_P(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_Pt(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0NPt(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0PPt(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0Mass(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0DCA(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0Radius(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0Lifetime(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0DaughterDCA(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0CPA(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0DCA2TPV(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0TPCpim(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0TPCpip(pname,cutSetK0s->GetMonitorOutput());
   
   if (isMC) {
     
     TString mode = "SPARSE";
     //TString mode = "HIST";
     //if (collSyst) mode = "SPARSE";
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("KStarPlusMinus_MotherMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kKaon0);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetMotherPDG(323);
     out->SetMotherMass(0.89166);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 90, 0.6, 1.5);
     out->AddAxis(ptID, 300, 0.0, 30.0);
     //out->AddAxis(k0sDCA, 10, 0.0, 1.0);
     
     if(isPP) out->AddAxis(centID, 400, 0.5, 400.5);
     else out->AddAxis(centID, 100, 0.0, 100.);
     //if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("AKStarPlusMinus_MotherMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kKaon0);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetMotherPDG(-323);
     out->SetMotherMass(0.89166);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 90, 0.6, 1.5);
     out->AddAxis(ptID, 300, 0.0, 30.0);
     //out->AddAxis(k0sDCA, 10, 0.0, 1.0);
     //if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     if(isPP) out->AddAxis(centID, 400, 0.5, 400.5);
     else out->AddAxis(centID, 100, 0.0, 100.);
   }
   
   return kTRUE;
}


//=============================


AliRsnMiniAnalysisTask * AddTask_pkx
(
 TString     lname = "pkx",
 Bool_t      isMC = kFALSE,
 Int_t       system = 0,
 Int_t       EventCuts = 0,
 Int_t       TrackCutsP = 0,
 Int_t       TrackCutsK = 0
)
{
  Bool_t      isPP = (system==0)?true:false;
  Int_t       aodFilterBit = 5;
  Int_t       evtCutSetID = 0;
  Int_t       MultBins = 0;// default for V0_M and MultBins = 2 for RefMult0_8 
  Int_t       customQualityCutsID=1; // for default
  AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV;
  AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV;
  Float_t     nsigmaPr = 2.0;
  Float_t     nsigmaKa = 2.0;
  Bool_t      enableMonitor = kTRUE;
  Bool_t      IsMcTrueOnly = kFALSE;
  UInt_t      triggerMask = AliVEvent::kINT7;
  Int_t       signedPdg = 3124;
  TString     monitorOpt = "NoSIGN";  //Flag for AddMonitorOutput.C e.g."NoSIGN"
  Bool_t      useCrossedRows = kTRUE;
  const char *yaxisVar = "";  //yaxisVar = "PtDaughter_PDaughter_cent"
  Bool_t      useMixLS = 0;
  Int_t       nmix = 5;
  Float_t     maxDiffVzMix = 1.0;
  Float_t     maxDiffMultMix = 5.0;
  TString     outNameSuffix = lname;
   
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
      ::Error("AddTaskRare_pp13::AddTask_pkx", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   TString taskName = Form("%s_%s%s", lname.Data(), (isPP? "pp" : "pPb"), (isMC ? "MC" : "Data"));
   
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);


   //if(is2011PbPb)
   //task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
   //else
   //task->SelectCollisionCandidates(triggerMask);
   int trigger=EventCuts%10;
   if(!trigger) task->UseESDTriggerMask(AliVEvent::kINT7);
   else if(trigger==1) task->UseESDTriggerMask(AliVEvent::kHighMult);
   else if(trigger==2) task->UseESDTriggerMask(AliVEvent::kHighMultV0);
   else task->UseESDTriggerMask(AliVEvent::kAny);

   /*
   AliMultSelectionTask *taskm = AddTaskMultSelection();

   if(isMC == 1)
     {
       taskm->SetAlternateOADBforEstimators( "LHC15f" );     
     }
   */

   //   task->UseMultiplicity("QUALITY");

   
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
   ::Info("AddTaskRare_pp13::AddTask_pkx", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n ", nmix, maxDiffVzMix, maxDiffMultMix));
   
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
     ::Info("AddTaskRare_pp13::AddTask_pkx", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp)?"ON":"OFF"));   
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
   

   //   AliRsnValuePair * openangle = new AliRsnValuePair("openangle", AliRsnValuePair::kDipAngle);
   //   openangle->DipAngle(1);


   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------

   
   //for systematic checks
     {
       //       gROOT->LoadMacro("$ALICE_ROOT/PWGLF/RESONANCES/macros/mini/ConfigLStar.C");
       //      gROOT->LoadMacro("ConfigureLstar13TeVpp.C");
       //gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigureLstar13TeVpp.C");
       if (!Config_pkx(task, isMC, isPP, "", cutsPair, TrackCutsP, TrackCutsK, aodFilterBit, customQualityCutsID, cutPrCandidate, cutKaCandidate, nsigmaPr, nsigmaKa,  enableMonitor, isMC&IsMcTrueOnly, signedPdg, monitorOpt, useCrossedRows, yaxisVar ,useMixLS)) 
return 0x0;  
     }

   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddTaskRare_pp13::AddTask_pkx - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",outNameSuffix.Data()), 
							   TList::Class(), 
							   AliAnalysisManager::kOutputContainer, 
							   outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}


Bool_t Config_pkx
(  
    AliRsnMiniAnalysisTask *task, 
    Bool_t                 isMC, 
    Bool_t                 isPP,
    const char             *suffix,
    AliRsnCutSet           *cutsPair,
    Int_t                  TrackCutsP = 0,
    Int_t                  TrackCutsK = 0,
    Int_t                  aodFilterBit = 5,
    Int_t                  customQualityCutsID=AliRsnCutSetDaughterParticle::kDisableCustom,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV,
    AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutKaCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV,
    Float_t                nsigmaPr = 2.0,
    Float_t                nsigmaKa = 2.0,
    Bool_t                 enableMonitor = kTRUE,
    Bool_t                 IsMcTrueOnly = kFALSE,
    Int_t                  signedPdg = 3124,
    TString                monitorOpt = "",  //Flag for AddMonitorOutput.C e.g."NoSIGN"
    Bool_t                 useCrossedRows = kTRUE,
    const char*            yaxisVar = "",  //yaxisVar = "PtDaughter_PDaughter_cent"
    Bool_t                 useMixLS = 0
)
{
  TString opt(yaxisVar);
  opt.ToUpper();

  Bool_t isPtDaughter  = opt.Contains("PTDAUGHTER") ;
  Bool_t isPDaughter   = opt.Contains("PDAUGHTER") ;
  Bool_t iscent        = opt.Contains("CENT") ;
  Bool_t iseta         = opt.Contains("ETA") ;
  Bool_t israpidity    = opt.Contains("RAPIDITY") ;


  // manage suffix
  if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
  
  // set daughter cuts
  AliRsnCutSetDaughterParticle * cutSetQ;
  AliRsnCutSetDaughterParticle * cutSetP;
  AliRsnCutSetDaughterParticle * cutSetK;


  AliRsnCutTrackQuality * trkQualityCut =  new AliRsnCutTrackQuality("myQualityCut");
  if(SetCustomQualityCut(trkQualityCut,customQualityCutsID,aodFilterBit)){
    //Set custom quality cuts for systematic checks
    
    cutSetQ=new AliRsnCutSetDaughterParticle(Form("cutQuality_bit%i",aodFilterBit),trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2011,AliPID::kPion,-1.);
    cutSetK=new AliRsnCutSetDaughterParticle(Form("cutKaon_%i_%2.1fsigma",cutKaCandidate, nsigmaKa),trkQualityCut,cutKaCandidate,AliPID::kKaon,nsigmaKa);
    cutSetP=new AliRsnCutSetDaughterParticle(Form("cutProton_%i_%2.1fsigma",cutPrCandidate, nsigmaPr),trkQualityCut,cutPrCandidate,AliPID::kProton,nsigmaPr);
    
    /*

    if(isPP) {
      Bool_t useCrossedRows = 1;
      cutSetQ  = new AliRsnCutSetDaughterParticle("cutQuality", AliRsnCutSetDaughterParticle::kQualityStd2010, AliPID::kPion, -1.0, aodFilterBit, useCrossedRows);
      cutSetP  = new AliRsnCutSetDaughterParticle(Form("cutProton_%2.1fsigma",nsigmaPr), cutPrCandidate, AliPID::kProton, nsigmaPr, aodFilterBit, useCrossedRows);
      cutSetK  = new AliRsnCutSetDaughterParticle(Form("cutKaon_%2.1f2sigma",nsigmaKa), cutKaCandidate, AliPID::kKaon, nsigmaKa, aodFilterBit, useCrossedRows);
    }
    else {
      Bool_t useCrossedRows = 1;
      cutSetQ  = new AliRsnCutSetDaughterParticle("cutQuality", AliRsnCutSetDaughterParticle::kQualityStd2011, AliPID::kPion, -1.0, aodFilterBit, useCrossedRows);
      cutSetQ->SetUse2011StdQualityCuts(kTRUE);
      cutSetP = new AliRsnCutSetDaughterParticle(Form("cutProton2011_%2.1fsigma",nsigmaPr), cutPrCandidate, AliPID::kProton, nsigmaPr, aodFilterBit, useCrossedRows);
      cutSetP->SetUse2011StdQualityCuts(kTRUE);
      cutSetK  = new AliRsnCutSetDaughterParticle(Form("cutKaon2011_%2.1f2sigma",nsigmaKa), cutKaCandidate, AliPID::kKaon, nsigmaKa, aodFilterBit, useCrossedRows);
      cutSetK->SetUse2011StdQualityCuts(kTRUE);
    }
    */
  }
    
  else{
    
    Bool_t useCrossedRows = 1;
    cutSetQ=new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit),AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.,aodFilterBit,useCrossedRows);
    cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",cutKaCandidate,nsigmaKa),cutKaCandidate,AliPID::kKaon,nsigmaKa,aodFilterBit,useCrossedRows);
    cutSetP=new AliRsnCutSetDaughterParticle(Form("cutPr%i_%2.1fsigma",cutPrCandidate,nsigmaPr),cutPrCandidate,AliPID::kProton,nsigmaPr,aodFilterBit,useCrossedRows);
  }

    

  
  Int_t iCutQ = task->AddTrackCuts(cutSetQ);
  Int_t iCutP = task->AddTrackCuts(cutSetP);
  Int_t iCutK = task->AddTrackCuts(cutSetK);
  
  if (enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetQ->GetMonitorOutput(), monitorOpt.Data());
    AddMonitorOutput(isMC, cutSetP->GetMonitorOutput(), monitorOpt.Data());
    AddMonitorOutput(isMC, cutSetK->GetMonitorOutput(), monitorOpt.Data());
  }  
  
  // -- Values ------------------------------------------------------------------------------------
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);
  /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt, kFALSE);
  /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt, kFALSE);
  /* 1st daughter p   */ Int_t fdp    = task->CreateValue(AliRsnMiniValue::kFirstDaughterP, kFALSE);
  /* 2nd daughter p   */ Int_t sdp    = task->CreateValue(AliRsnMiniValue::kSecondDaughterP, kFALSE);

  /* Dip Angle        */ Int_t OpAn   = task->CreateValue(AliRsnMiniValue::kDipAngle, kFALSE);
  /* kPtRatio         */ Int_t PtRat  = task->CreateValue(AliRsnMiniValue::kPtRatio, kFALSE);


  /* kptresolution    */ Int_t respT  = task->CreateValue(AliRsnMiniValue::kPairPtRes, kFALSE);
  /* kmassresolution  */ Int_t resmass  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kFALSE);



  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges
  // [0] = unlike
  // [1] = mixing
  // [2] = like ++
  // [3] = like --

  Bool_t  use     [12] = { !IsMcTrueOnly, !IsMcTrueOnly, !IsMcTrueOnly,!IsMcTrueOnly,!IsMcTrueOnly,!IsMcTrueOnly, isMC    ,  isMC   , 0      ,   0    , useMixLS , useMixLS};
  Bool_t  useIM   [12] = { 1            ,  1           ,  1           ,  1          ,  1          , 1           ,  1      ,   1     , 0      ,   0    , 1        , 1       };
  TString name    [12] = {"UnlikePM"    , "UnlikeMP"   , "MixingPM"   , "MixingMP"  , "LikePP"    , "LikeMM"    ,"TruesPM","TruesMP","ResPM" ,"ResMP" ,"MixingPP","MixingMM"};
  TString comp    [12] = {"PAIR"        , "PAIR"       , "MIX"        , "MIX"       , "PAIR"      , "PAIR"      , "TRUE"  , "TRUE"  , "TRUE" , "TRUE" , "MIX"    ,"MIX"   };
  TString output  [12] = {"SPARSE"      , "SPARSE"     , "SPARSE"     , "SPARSE"    , "SPARSE"    , "SPARSE"    , "SPARSE","SPARSE" ,"SPARSE","SPARSE","SPARSE"  ,"SPARSE"};
  Char_t  charge1 [12] = {'+'           , '-'          , '+'          , '-'         , '+'         , '-'         , '+'     ,  '-'    , '+'    ,  '-'   , '+'      , '-'    };
  Char_t  charge2 [12] = {'-'           , '+'          , '-'          , '+'         , '+'         , '-'         , '-'     ,  '+'    , '-'    ,  '+'   ,'+'       , '-'    };
  Int_t   cutID1  [12] = { iCutP        ,  iCutP       ,  iCutP       ,  iCutP      ,  iCutP      ,  iCutP      ,  iCutP  ,  iCutP  ,  iCutP ,  iCutP , iCutP    , iCutP };
  Int_t   cutID2  [12] = { iCutK        ,  iCutK       ,  iCutK       ,  iCutK      ,  iCutK      ,  iCutK      ,  iCutK  ,  iCutK  ,  iCutK ,  iCutK , iCutK    , iCutK };
  
  //TString output  [10] = {"HIST"   , "HIST"   , "HIST"   , "HIST"   , "HIST"  , "HIST"  , "HIST"  ,  "HIST"  , "HIST"  ,  "HIST"  };


  for (Int_t i = 0; i < 12; i++) {
    if (!use[i]) continue;
    AliRsnMiniOutput *out = task->CreateOutput(Form("pkx_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
    out->SetCutID(0, cutID1[i]);
    out->SetCutID(1, cutID2[i]);
    out->SetDaughter(0, AliRsnDaughter::kProton);
    out->SetDaughter(1, AliRsnDaughter::kKaon);
    out->SetCharge(0, charge1[i]);
    out->SetCharge(1, charge2[i]);
    out->SetMotherPDG(signedPdg);
    out->SetMotherMass(1.51953);
    out->SetPairCuts(cutsPair);

    // axis X: invmass (or resolution)
    if (useIM[i]) 
      //      out->AddAxis(imID, 300, 1.4, 1.7);
      //     out->AddAxis(imID, 300, 1.4, 2.2);
      out->AddAxis(imID, 800, 1.4, 3.);
    else
      out->AddAxis(resID, 200, -0.02, 0.02);
    

    // axis Y: transverse momentum of pair as default - else chosen value



    out->AddAxis(ptID, 100, 0.0, 10.0); //default use mother pt
    //out->AddAxis(OpAn,  50, 0.0, 25.0);
    //out->AddAxis(PtRat, 20, 0.0, 1.0);


    //out->AddAxis(centID, 110, 0.0, 110.0); // adding multiplicity axis


    if(isPtDaughter) {
      out->AddAxis(fdpt, 100, 0.0, 10.0);
      out->AddAxis(sdpt, 100, 0.0, 10.0);     }
    
    if(isPDaughter) {
	out->AddAxis(fdp, 100, 0.0, 10.0);
	out->AddAxis(sdp, 100, 0.0, 10.0);    }
    
    // axis Z: centrality-multiplicity
    if(iscent) {
      if (!isPP)	out->AddAxis(centID, 100, 0.0, 100.0);
      else       	out->AddAxis(centID, 400, 0.0, 400.0);      }
    // axis W: pseudorapidity

    if(iseta)       out->AddAxis(etaID, 20, -1.0, 1.0);
    // axis J: rapidity
    if(israpidity)  out->AddAxis(yID, 12, -0.6, 0.6);
    
  }
  
  if (isMC){   
    // create output
    AliRsnMiniOutput *outm = task->CreateOutput(Form("Lstar_Mother%s", suffix), "SPARSE", "MOTHER");
    outm->SetDaughter(0, AliRsnDaughter::kProton);
    outm->SetDaughter(1, AliRsnDaughter::kKaon);
    outm->SetMotherPDG(signedPdg);
    outm->SetMotherMass(1.51953);
    // pair cuts
    outm->SetPairCuts(cutsPair);
    // binnings
    outm->AddAxis(imID, 800, 1.4, 2.2);

    outm->AddAxis(ptID, 100, 0.0, 10.0); //default use mother pt
    //outm->AddAxis(OpAn, 20, -10.0, 10.0);
    //outm->AddAxis(PtRat, 10, -1.0, 1.0);

    // axis Y: transverse momentum of pair as default - else chosen value
    outm->AddAxis(centID, 110, 0.0, 110.0); // adding multiplicity axis


    if(isPtDaughter) {
      outm->AddAxis(fdpt, 100, 0.0, 10.0);
      outm->AddAxis(sdpt, 100, 0.0, 10.0);     }
    
    if(isPDaughter) {
	outm->AddAxis(fdp, 100, 0.0, 10.0);
	outm->AddAxis(sdp, 100, 0.0, 10.0);    }
    
    // axis Z: centrality-multiplicity
    if(iscent) {
      if (!isPP)	outm->AddAxis(centID, 100, 0.0, 100.0);
      else       	outm->AddAxis(centID, 400, 0.0, 400.0);      }
    // axis W: pseudorapidity

    if(iseta)       outm->AddAxis(etaID, 20, -1.0, 1.0);
    // axis J: rapidity
    if(israpidity)  outm->AddAxis(yID, 12, -0.6, 0.6);



     AliRsnMiniOutput *outam = task->CreateOutput(Form("AntiLstar_Mother%s", suffix), "SPARSE", "MOTHER");
    outam->SetDaughter(0, AliRsnDaughter::kProton);
    outam->SetDaughter(1, AliRsnDaughter::kKaon);
    outam->SetMotherPDG(-3124);
    outam->SetMotherMass(1.51953);
    // pair cuts
    outam->SetPairCuts(cutsPair);
    // binnings
    outam->AddAxis(imID, 800, 1.4, 2.2);

    outam->AddAxis(ptID, 100, 0.0, 10.0); //default use mother pt
    //outm->AddAxis(OpAn, 20, -10.0, 10.0);
    //outm->AddAxis(PtRat, 10, -1.0, 1.0);

    // axis Y: transverse momentum of pair as default - else chosen value
    outam->AddAxis(centID, 110, 0.0, 110.0); // adding multiplicity axis



    AliRsnMiniOutput *outnew1 = task->CreateOutput(Form("Lstar_resolution1%s", suffix), "SPARSE", "TRUE");

    outnew1->SetCutID(0, cutID1[0]);
    outnew1->SetCutID(1, cutID2[0]);
    outnew1->SetDaughter(0, AliRsnDaughter::kProton);
    outnew1->SetDaughter(1, AliRsnDaughter::kKaon);
    outnew1->SetCharge(0, charge1[0]);
    outnew1->SetCharge(1, charge2[0]);
    outnew1->SetMotherPDG(signedPdg);
    outnew1->SetMotherMass(1.51953);
    outnew1->SetPairCuts(cutsPair);



    //    outnew1->SetDaughter(0, AliRsnDaughter::kProton);
    //    outnew1->SetDaughter(1, AliRsnDaughter::kKaon);
    //    outnew1->SetMotherPDG(signedPdg);
    //    outnew1->SetMotherMass(1.51953);
    //    outnew1->SetPairCuts(cutsPair);

    outnew1->AddAxis(imID, 800, 1.4, 2.2); // mass axis
    outnew1->AddAxis(ptID, 100, 0.0, 10.0); //pt axis
    outnew1->AddAxis(resmass, 400, -0.2, 0.2);// mass resolution
    
    
    
    AliRsnMiniOutput *outnew2 = task->CreateOutput(Form("Lstar_resolution2%s", suffix), "SPARSE", "TRUE");
  
    outnew2->SetCutID(0, cutID1[0]);
    outnew2->SetCutID(1, cutID2[0]);
    outnew2->SetDaughter(0, AliRsnDaughter::kProton);
    outnew2->SetDaughter(1, AliRsnDaughter::kKaon);
    outnew2->SetCharge(0, charge1[0]);
    outnew2->SetCharge(1, charge2[0]);
    outnew2->SetMotherPDG(signedPdg);
    outnew2->SetMotherMass(1.51953);
    outnew2->SetPairCuts(cutsPair);

    outnew2->AddAxis(imID, 800, 1.4, 2.2); // mass axis
    outnew2->AddAxis(ptID, 100, 0.0, 10.0); //pt axis
    outnew2->AddAxis(resmass, 400, -0.2, 0.2);// mass resolution
  }  // if(isMC)

  return kTRUE;
}


Bool_t SetCustomQualityCut(AliRsnCutTrackQuality * trkQualityCut, Int_t customQualityCutsID = 0, Int_t customFilterBit = 0)
{
  //Sets configuration for track quality object different from std quality cuts.
  //Returns kTRUE if track quality cut object is successfully defined,
  //returns kFALSE if an invalid set of cuts (customQualityCutsID) is chosen or if the
  //object to be configured does not exist.

  if ((!trkQualityCut)){
    Printf("::::: SetCustomQualityCut:: use default quality cuts specified in task configuration.");
    return kFALSE;
  }

  
  //for pA 2013
  //trkQualityCut->SetDefaults2011();//with filter bit=10
  //reset filter bit to very loose cuts 
  trkQualityCut->SetAODTestFilterBit(customFilterBit); 
  //apply all other cuts "by hand"
  trkQualityCut->SetCheckOnlyFilterBit(kFALSE);
  trkQualityCut->SetMinNCrossedRowsTPC(70, kTRUE);
  trkQualityCut->SetTPCminNClusters(70);
  trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(0.8, kTRUE);
  trkQualityCut->SetMaxChi2TPCConstrainedGlobal(36);//used for ESD only - for AOD does not correspond to any cut
  trkQualityCut->SetTPCmaxChi2(4.0); //already in filter bit 0
  trkQualityCut->SetRejectKinkDaughters(kTRUE); //already in filter bit 0
  trkQualityCut->SetSPDminNClusters(AliESDtrackCuts::kAny);
  trkQualityCut->SetITSmaxChi2(36);
  trkQualityCut->SetDCARPtFormula("0.0105+0.0350/pt^1.1");
  trkQualityCut->SetDCAZmax(2.0); 
  trkQualityCut->AddStatusFlag(AliESDtrack::kTPCin   , kTRUE);//already in defaults 2011
  trkQualityCut->AddStatusFlag(AliESDtrack::kTPCrefit, kTRUE);//already in defaults 2011
  trkQualityCut->AddStatusFlag(AliESDtrack::kITSrefit, kTRUE);//already in defaults 2011

   if (customQualityCutsID==AliRsnCutSetDaughterParticle::kFilterBitCustom) {
    trkQualityCut->SetCheckOnlyFilterBit(kTRUE);
  } 
  

  if (customQualityCutsID==2){
    trkQualityCut->SetDCARPtFormula("0.006+0.0200/pt^1.1");
  }

  if (customQualityCutsID==3){
    trkQualityCut->SetDCARPtFormula("0.0150+0.0500/pt^1.1");
  }


  if (customQualityCutsID==4){
    trkQualityCut->SetDCAZmax(3.2);
  } 

  if (customQualityCutsID==5) {
    trkQualityCut->SetDCAZmax(5.0); 
  }

  
  if (customQualityCutsID==6){
    trkQualityCut->SetMinNCrossedRowsTPC(90, kTRUE);
  }

  if (customQualityCutsID==7){
    trkQualityCut->SetMinNCrossedRowsTPC(80, kTRUE);
  }

  


  if (customQualityCutsID==8){
    trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(0.9, kTRUE);
  }
  
  if (customQualityCutsID==9){
    trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(0.85, kTRUE);
  }



  if (customQualityCutsID==10){
    trkQualityCut->SetTPCmaxChi2(5.0);
  }
  
  if (customQualityCutsID==11){
    trkQualityCut->SetTPCmaxChi2(6.0);
  }
  



  if (customQualityCutsID==12){
    trkQualityCut->SetITSmaxChi2(38.0);
  }
  
  if (customQualityCutsID==13){
    trkQualityCut->SetITSmaxChi2(40.0);
  }



  /*
  if (customQualityCutsID==3){
    trkQualityCut->SetDCAZmax(3.2);
  } else {
    trkQualityCut->SetDCAZmax(2.0); 
  }
  
  if (customQualityCutsID==4){
    trkQualityCut->SetMinNCrossedRowsTPC(60, kTRUE);
  }
  
  if (customQualityCutsID==5){
    trkQualityCut->SetMinNCrossedRowsTPC(80, kTRUE);
  }
  
  if (customQualityCutsID==6){
    trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(0.75, kTRUE);
  }
  
  if (customQualityCutsID==7){
    trkQualityCut->SetMinNCrossedRowsOverFindableClsTPC(0.85, kTRUE);
  }
  
  if (customQualityCutsID==8){
    trkQualityCut->SetAODTestFilterBit(10);
    trkQualityCut->SetTPCminNClusters(70);
  }
  
  if (customQualityCutsID==9){
    trkQualityCut->SetTPCmaxChi2(3.5);
  }
  */

  trkQualityCut->SetPtRange(0.15, 100.0);
  trkQualityCut->SetEtaRange(-0.8, 0.8);
  
  Printf(Form("::::: SetCustomQualityCut:: using custom track quality cuts #%i",customQualityCutsID));

  trkQualityCut->Print();
  return kTRUE;
}


//=============================


AliRsnMiniAnalysisTask *AddTask_pk0
(
 TString     lname,
 Bool_t      isMC,
 Int_t       system,
 Int_t       EventCuts,
 Int_t       TrackCutsP,
 Int_t       TrackCutsK
 )
{
  Float_t     cutV = 10.0;
  Int_t       evtCutSetID = 0;
  Int_t       pairCutSetID = 0;
  Int_t       mixingConfigID = 1;
  Int_t       aodFilterBit = 5;
  Bool_t      enableMonitor=kTRUE;
  TString     monitorOpt="pp";
  Float_t     pPIDCut = 2.0;
  Float_t     pi_k0s_PIDCut = 5.0;
  //Float_t     trackDCAcut = 7.0;
  Float_t     massTol = 0.03;
  Float_t     massTolVeto = 0.004;
  Float_t     pLife = 20;
  Float_t     radiuslow = 0.5;
  Float_t     radiushigh = 200;   
  Bool_t      Switch = kFALSE;
  Float_t     k0sDCA = 0.3;
  Float_t     k0sCosPoinAn = 0.97;
  Float_t     k0sDaughDCA = 1.0;
  Int_t       NTPCcluster = 70;
  Float_t     maxDiffVzMix = 1.0;
  Float_t     maxDiffMultMix = 10.0;
  Float_t     maxDiffAngleMixDeg = 20.0;
  Int_t       aodN = 68;
  TString     outNameSuffix = lname;
  Int_t       centr = 0;
  Bool_t      ptDep = kTRUE;
   
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
      ::Error("AddTaskRare_pp13::AddTask_pk0", "No analysis manager to connect to.");
      return NULL;
   } 
   
   // create the task and configure 
   TString taskName = Form("%s_%s%s", lname.Data(), ((system==0)? "pp" : "PbPb"), (isMC ? "MC" : "Data"));
 
   AliRsnMiniAnalysisTask* task = new AliRsnMiniAnalysisTask(taskName.Data(),isMC);

   //task->SelectCollisionCandidates(triggerMask); //AOD
   int trigger=EventCuts%10;
   if(!trigger) task->UseESDTriggerMask(AliVEvent::kINT7);
   else if(trigger==1) task->UseESDTriggerMask(AliVEvent::kHighMult);
   else if(trigger==2) task->UseESDTriggerMask(AliVEvent::kHighMultV0);
   else task->UseESDTriggerMask(AliVEvent::kAny);

   //if(isPP) 
   task->UseMultiplicity("QUALITY");
   //else task->UseCentrality("V0M");

   // set event mixing options                                                                                                               
   task->UseContinuousMix();
   //task->UseBinnedMix();                                                                                                                   
   task->SetNMix(nmix);
   task->SetMaxDiffVz(maxDiffVzMix);
   task->SetMaxDiffMult(maxDiffMultMix);
   ::Info("AddTaskRare_pp13::AddTask_pk0", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %\5.3f", nmix, maxDiffVzMix, maxDiffMultMix));
   
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
     ::Info("AddTaskRare_pp13::AddTask_pk0", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp)?"ON":"OFF"));
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
   cutsPair->SetCutScheme(cutY->GetName());
   
   //
   // -- CONFIG ANALYSIS --------------------------------------------------------------------------

   Bool_t isPP=kTRUE;
   if (isMC) {
     Printf("========================== MC analysis - PID cuts not used"); 
   } else 
     Printf("========================== DATA analysis - PID cuts used");
   
   if (!Config_pk0(task, isPP, isMC, TrackCutsP, TrackCutsK, pPIDCut, pi_k0s_PIDCut, aodFilterBit, enableMonitor, monitorOpt.Data(), massTol, massTolVeto, pLife, radiuslow, radiushigh, Switch, k0sDCA, k0sCosPoinAn, k0sDaughDCA, NTPCcluster, "", cutsPair, ptDep)) return 0x0;
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddTaskRare_pp13::AddTask_pk0 - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",lname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}


Bool_t Config_pk0
(  
   AliRsnMiniAnalysisTask *task,
   Bool_t                  isPP,
   //Int_t		           collSyst, 
   Bool_t                  isMC,
   Int_t                   TrackCutsP,
   Int_t                   TrackCutsK,
   Float_t                 nsigmaPr,
   Float_t                 pi_k0s_PIDCut,
   Int_t                   aodFilterBit,
   //Float_t               trackDCAcut,
   Bool_t                  enableMonitor=kTRUE  ,
   TString                 monitorOpt="",
   Float_t                 massTol,
   Float_t                 massTolVeto, 
   Float_t                 pLife, 
   Float_t                 radiuslow,
   Float_t                 radiushigh,    
   Bool_t                  Switch,     
   Float_t                 k0sDCA,
   Float_t                 k0sCosPoinAn,
   Float_t                 k0sDaughDCA,
   Int_t                   NTPCcluster,
   const char             *suffix,
   AliRsnCutSet           *cutsPair,
   Bool_t                  ptDep    
)
{
  Int_t customQualityCutsID=0;
  AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate=AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV;

   // manage suffix
   if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
 
   /////////////////////////////////////////////////////
   // selections for the proton
   /////////////////////////////////////////////////////

  AliRsnCutSetDaughterParticle* cutSetQ;
  AliRsnCutSetDaughterParticle* cutSetP;

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  if(SetCustomQualityCut(trkQualityCut,customQualityCutsID,aodFilterBit)){
    //Set custom quality cuts for systematic checks
    cutSetQ=new AliRsnCutSetDaughterParticle(Form("cutQuality_bit%i",aodFilterBit),trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2011,AliPID::kPion,-1.);
    cutSetP=new AliRsnCutSetDaughterParticle(Form("cutProton_%i_%2.1fsigma",cutPrCandidate, nsigmaPr),trkQualityCut,cutPrCandidate,AliPID::kProton,nsigmaPr);
  }else{
    Bool_t useCrossedRows = 1;
    cutSetQ=new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit),AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.,aodFilterBit,useCrossedRows);
    cutSetP=new AliRsnCutSetDaughterParticle(Form("cutPr%i_%2.1fsigma",cutPrCandidate,nsigmaPr),cutPrCandidate,AliPID::kProton,nsigmaPr,aodFilterBit,useCrossedRows);
  }

  Int_t iCutQ = task->AddTrackCuts(cutSetQ);
  Int_t iCutP = task->AddTrackCuts(cutSetP);
  
  if (enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetQ->GetMonitorOutput(), monitorOpt.Data());
    AddMonitorOutput(isMC, cutSetP->GetMonitorOutput(), monitorOpt.Data());
  }

   //
   /////////////////////////////////////////////////////////////
   // selections for K0s and for the daughters of K0s
   /////////////////////////////////////////////////////////////
   // 
   // selections for pion daugthers of K0s
   AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("qualityDaughterK0s");   
   esdTrackCuts->SetPtRange(0.15,1.E10);
   esdTrackCuts->SetEtaRange(-0.8,0.8);
   esdTrackCuts->SetRequireTPCRefit();
   esdTrackCuts->SetAcceptKinkDaughters(0); //
   esdTrackCuts->SetMinNClustersTPC(NTPCcluster);
   esdTrackCuts->SetMaxChi2PerClusterTPC(4);
    
    if(ptDep){
   esdTrackCuts->SetMinDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    }else
   esdTrackCuts->SetMinDCAToVertexXY("0.06"); //Use one of the two - pt dependent or fixed value cut.
  
   //
   /////////////////////////////////////////////////
   // selections for K0s
   AliRsnCutV0 *cutK0s = new AliRsnCutV0("cutK0s", kK0Short, AliPID::kPion, AliPID::kPion);
   //cutK0s->SetPIDCutProton(pPIDCut);       // PID for the proton daughter of K0s
   cutK0s->SetPIDCutPion(pi_k0s_PIDCut);        // PID for the pion daughter of K0s
   cutK0s->SetESDtrackCuts(esdTrackCuts);  // all the other selections (defined above) for proton and pion daughters of K0s
   cutK0s->SetMaxDaughtersDCA(k0sDaughDCA);
   cutK0s->SetMaxDCAVertex(k0sDCA);
   cutK0s->SetMinCosPointingAngle(k0sCosPoinAn);
   cutK0s->SetTolerance(massTol);
   cutK0s->SetToleranceVeto(massTolVeto);   //Rejection range for Competing V0 Rejection
   cutK0s->SetSwitch(Switch);    
   cutK0s->SetfLife(pLife); 
   cutK0s->SetfLowRadius(radiuslow); 
   cutK0s->SetfHighRadius(radiushigh); 
   cutK0s->SetMaxRapidity(2.0);
   //
   AliRsnCutSet *cutSetK0s = new AliRsnCutSet("setK0s", AliRsnTarget::kDaughter);
   cutSetK0s->AddCut(cutK0s);
   cutSetK0s->SetCutScheme(cutK0s->GetName());
   Int_t iCutK0s = task->AddTrackCuts(cutSetK0s);

   //
   // -- Values ------------------------------------------------------------------------------------
   //
   
   /// -- Values ------------------------------------------------------------------------------------                                        
    /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
   /* transv. momentum */ Int_t ptID    = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
   /* centrality       */ Int_t centID  = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
   /* pseudorapidity   */ Int_t etaID   = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
   /* rapidity         */ Int_t yID     = task->CreateValue(AliRsnMiniValue::kY, kFALSE);

   
   //
   // -- Create all needed outputs -----------------------------------------------------------------
   //
   
   // use an array for more compact writing, which are different on mixing and charges
   // [0] = unlike
   // [1] = mixing
   // [2] = like ++
   // [3] = like --
   Bool_t  use     [6] = {1               ,1                ,1                  ,1                   ,1                ,1                 };
   Bool_t  useIM   [6] = {1               ,1                ,1                  ,1                   ,1                ,1                 };
   TString name    [6] = {"K0Pp"          ,"K0Pm"           ,"K0PpMix"          ,"K0PmMix"           ,"KStarPlusMinust","AKStarPlusMinust"};
   TString comp    [6] = {"PAIR"          ,"PAIR"           ,"MIX"              ,"MIX"               ,"TRUE"           ,"TRUE"            };
   TString output  [6] = {"SPARSE"        ,"SPARSE"         ,"SPARSE"           ,"SPARSE"            ,"SPARSE"         ,"SPARSE"          };
   Char_t  charge1 [6] = {'0'             ,'0'              ,'0'                ,'0'                 ,'0'              ,'0'               };
   Char_t  charge2 [6] = {'+'             ,'-'              ,'+'                ,'-'                 ,'+'              ,'-'               };
   Int_t   cutID1  [6] = { iCutK0s        ,iCutK0s          ,iCutK0s            ,iCutK0s             ,iCutK0s          ,iCutK0s           };
   Int_t   cutID2  [6] = { iCutP          ,iCutP            ,iCutP              ,iCutP               ,iCutP            ,iCutP             };
   Int_t   ipdg    [6] = {323             ,-323             ,323                ,-323                ,323              ,-323              };
   Double_t mass   [6] = { 0.89166        ,0.89166          ,0.89166            ,0.89166             ,0.89166          ,0.89166           };
   
   for (Int_t i = 0; i < 6; i++) {
     if (!use[i]) continue;
     //if (collSyst) output[i] = "SPARSE";
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("pk0_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
     // selection settings
     out->SetCutID(0, cutID1[i]);
     out->SetCutID(1, cutID2[i]);
     out->SetDaughter(0, AliRsnDaughter::kKaon0);
     out->SetDaughter(1, AliRsnDaughter::kProton);
     out->SetCharge(0, charge1[i]);
     out->SetCharge(1, charge2[i]);
     out->SetMotherPDG(ipdg[i]);
     out->SetMotherMass(mass[i]);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // axis X: invmass
     if (useIM[i]) 
       out->AddAxis(imID, 785, 1.43, 3.);
     // axis Y: transverse momentum
     out->AddAxis(ptID, 300, 0.0, 30.0);
   //  out->AddAxis(k0sDCA, 10, 0.0, 1.0);
     
     //if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     //if(isPP) out->AddAxis(centID, 400, 0.5, 400.5);
     //else out->AddAxis(centID, 100, 0.0, 100.);
   } 

   TString pname="p";
   AddMonitorOutput_Pt(pname,cutSetP->GetMonitorOutput());
   AddMonitorOutput_Eta(pname,cutSetP->GetMonitorOutput());
   AddMonitorOutput_DCAxy(pname,cutSetP->GetMonitorOutput());
   AddMonitorOutput_DCAz(pname,cutSetP->GetMonitorOutput());
   AddMonitorOutput_TPCpi(pname,cutSetP->GetMonitorOutput());
   AddMonitorOutput_NclTPC(pname,cutSetP->GetMonitorOutput());
   AddMonitorOutput_chi2TPC(pname,cutSetP->GetMonitorOutput());

   pname.Form("k0");
   AddMonitorOutput_P(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_Pt(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0NPt(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0PPt(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0Mass(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0DCA(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0Radius(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0Lifetime(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0DaughterDCA(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0CPA(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0DCA2TPV(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0TPCpim(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0TPCpip(pname,cutSetK0s->GetMonitorOutput());
   
   if (isMC) {
     TString mode = "SPARSE";
     //TString mode = "HIST";
     //if (collSyst) mode = "SPARSE";
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("KStarPlusMinus_MotherMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kKaon0);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetMotherPDG(323);
     out->SetMotherMass(0.89166);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 90, 0.6, 1.5);
     out->AddAxis(ptID, 300, 0.0, 30.0);
     //out->AddAxis(k0sDCA, 10, 0.0, 1.0);
     
     if(isPP) out->AddAxis(centID, 400, 0.5, 400.5);
     else out->AddAxis(centID, 100, 0.0, 100.);
     //if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("AKStarPlusMinus_MotherMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kKaon0);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetMotherPDG(-323);
     out->SetMotherMass(0.89166);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 90, 0.6, 1.5);
     out->AddAxis(ptID, 300, 0.0, 30.0);
     //out->AddAxis(k0sDCA, 10, 0.0, 1.0);
     //if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     if(isPP) out->AddAxis(centID, 400, 0.5, 400.5);
     else out->AddAxis(centID, 100, 0.0, 100.);
   }
   
   return kTRUE;
}


//=============================


AliRsnMiniAnalysisTask* AddTask_Lambdapi(TString lname,Bool_t isMC,Int_t system,Int_t EventCuts,Int_t TrackCutsLambda,Int_t TrackCutsPi){
  Bool_t      isPP = true;
  if(system) isPP = false;
  Float_t     cutV = 10.0;
  Int_t       aodFilterBit = 5;
  Float_t     piPIDCut = 3.0;
  Float_t     pPIDCut = 3.0;
  Float_t     trackDCAcut = 7.0;
  Float_t     massTol = 0.006;
  Float_t     lambdaDCA = 0.3;
  Float_t     lambdaCosPoinAn = 0.99;
  Float_t     lambdaDaughDCA = 0.5;
  Int_t       NTPCcluster = 70;
  Int_t       nmix = 5;
  Float_t     maxDiffVzMix = 1.0;
  Float_t     maxDiffMultMix = 2.0;
  Float_t     maxDiffAngleMixDeg = 20.0;
  Int_t       aodN = 0;
  TString     outNameSuffix = lname;
     
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskRare_pp13::AddTask_Lambdapi", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   TString taskName = Form("%s_%s%s", lname.Data(), (isPP? "pp" : "PbPb"), (isMC ? "MC" : "Data"));
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   if (!isMC && !isPP){
     Printf(Form("========== SETTING USE CENTRALITY PATCH AOD049 : %s", (aodN==49)? "yes" : "no"));
     task->SetUseCentralityPatch(aodN==49);
   }
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
   ::Info("AddTaskRare_pp13::AddTask_Lambdapi", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n max diff EP angle = %5.3f deg", nmix, maxDiffVzMix, maxDiffMultMix, (isPP ? 0.0 : maxDiffAngleMixDeg)));

   int trigger=EventCuts%10;
   if(!trigger) task->UseESDTriggerMask(AliVEvent::kINT7);
   else if(trigger==1) task->UseESDTriggerMask(AliVEvent::kHighMult);
   else if(trigger==2) task->UseESDTriggerMask(AliVEvent::kHighMultV0);
   else task->UseESDTriggerMask(AliVEvent::kAny);
   
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
   AliRsnMiniOutput *outPlane = 0x0; //task->CreateOutput("eventPlane", "HIST", "EVENT");
   if (!isPP){
     outPlane = task->CreateOutput("eventPlane", "HIST", "EVENT");
     outPlane->AddAxis(planeID, 180, 0.0, TMath::Pi());
     }
   
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
   //gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigSigmaStar.C");
   if (isMC) {
       Printf("========================== MC analysis - PID cuts not used");
   } else 
     Printf("========================== DATA analysis - PID cuts used");
   if (!Config_Lambdapi(task, system, isMC, TrackCutsLambda, TrackCutsPi, piPIDCut, pPIDCut, aodFilterBit, trackDCAcut, massTol, lambdaDCA, lambdaCosPoinAn, lambdaDaughDCA, NTPCcluster, "", cutsPair)) return 0x0;
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddTaskRare_pp13::AddTask_Lambdapi - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",lname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}


Bool_t Config_Lambdapi
(  
   AliRsnMiniAnalysisTask *task,
   Int_t		   collSyst, 
   Bool_t                  isMC,
   Int_t                   TrackCutsLambda,
   Int_t                   TrackCutsPi,
   Float_t                 piPIDCut,
   Float_t                 pPIDCut,
   Int_t                   aodFilterBit,
   Float_t                 trackDCAcut,
   Float_t                 massTol,
   Float_t                 lambdaDCA,
   Float_t                 lambdaCosPoinAn,
   Float_t                 lambdaDaughDCA,
   Int_t                   NTPCcluster,
   const char             *suffix,
   AliRsnCutSet           *cutsPair
)
{
  Float_t     pLife = 20;
  Float_t     radiuslow = 0.5;
  Float_t     radiushigh = 200;
  Bool_t      enableMonitor=kTRUE;
  TString     monitorOpt="pp";

   // manage suffix
   if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
 
   /////////////////////////////////////////////////////
   // selections for the pion from the decay of Sigma*
   /////////////////////////////////////////////////////
   //
   AliRsnCutDaughterSigmaStar2010PP *cutPi = new AliRsnCutDaughterSigmaStar2010PP("cutPionForLambdapi", AliPID::kPion);
   cutPi->SetPIDCut(piPIDCut);    // fPIDCut used in IsSelected() after the call to cutQuality
   //cutPi->SetMinTPCcluster(NTPCcluster);

   AliRsnCutTrackQuality *cutQuality = (AliRsnCutTrackQuality*) cutPi->CutQuality();
   cutQuality->SetDefaults2011();
   //cutQuality->SetDCARmax(trackDCAcut);
   //cutQuality->SetDefaults2010(0,1);  // 1st par. not default (0 -> use TPC clusters). 2nd par. default (-> standard Pt and eta range)
   // SetDefaults2010 contains the following selections:
   //     SetPtRange(0.15, 1E+20);
   //     SetEtaRange(-0.8, 0.8);
   //     and from aliroot/master/src/ANALYSIS/ANALYSISalice/AliESDtrackCuts.cxx
   //     AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(1,0)
   //         esdTrackCuts->SetMinNClustersTPC(70);
   //         esdTrackCuts->SetMaxChi2PerClusterTPC(4);
   //         esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
   //         esdTrackCuts->SetRequireTPCRefit(kTRUE);
   //         esdTrackCuts->SetRequireITSRefit(kTRUE);
   //         esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
   //         esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");    // NB. With pt_min=0.15 (see above) -> DCAxy_max = 0.2560
   //         esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
   //         esdTrackCuts->SetMaxDCAToVertexZ(2);
   //         esdTrackCuts->SetDCAToVertex2D(kFALSE);
   //         esdTrackCuts->SetRequireSigmaToVertex(kFALSE);  
   //         esdTrackCuts->SetMaxChi2PerClusterITS(36);
   //  
   AliRsnCutSet *cutSetPi = new AliRsnCutSet("setPionForLambdapi", AliRsnTarget::kDaughter);
   cutSetPi->AddCut(cutPi);
   cutSetPi->SetCutScheme(cutPi->GetName());
   Int_t iCutPi = task->AddTrackCuts(cutSetPi);
   //
   /////////////////////////////////////////////////////////////
   // selections for Lambda and for the daughters of Lambda 
   /////////////////////////////////////////////////////////////
   // 
   // selections for the proton and pion daugthers of Lambda and AntiLambda
   AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("qualityDaughterLambda");   
   esdTrackCuts->SetPtRange(0.15,1.E10);
   esdTrackCuts->SetEtaRange(-0.8,0.8);
   esdTrackCuts->SetRequireTPCRefit();
   esdTrackCuts->SetAcceptKinkDaughters(0);
   esdTrackCuts->SetMinNClustersTPC(NTPCcluster);
   esdTrackCuts->SetMaxChi2PerClusterTPC(4);
   esdTrackCuts->SetMinDCAToVertexXY(0.15);   
   //
   /////////////////////////////////////////////////
   // selections for Lambda
   AliRsnCutV0 *cutLambda = new AliRsnCutV0("cutLambda", kLambda0, AliPID::kProton, AliPID::kPion);
   cutLambda->SetPIDCutProton(pPIDCut);       // PID for the proton daughter of Lambda
   cutLambda->SetPIDCutPion(piPIDCut);        // PID for the pion daughter of Lambda 
   cutLambda->SetESDtrackCuts(esdTrackCuts);  // all the other selections (defined above) for proton and pion daughters of Lambda
   cutLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
   cutLambda->SetMaxDCAVertex(lambdaDCA);
   cutLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
   cutLambda->SetTolerance(massTol);
   cutLambda->SetSwitch(false);
   cutLambda->SetfLife(pLife); 
   cutLambda->SetfLowRadius(radiuslow); 
   cutLambda->SetfHighRadius(radiushigh); 
   cutLambda->SetMaxRapidity(2.);
   cutLambda->SetMinTPCcluster(NTPCcluster);
   //
   AliRsnCutSet *cutSetLambda = new AliRsnCutSet("setLambda", AliRsnTarget::kDaughter);
   cutSetLambda->AddCut(cutLambda);
   cutSetLambda->SetCutScheme(cutLambda->GetName());
   Int_t iCutLambda = task->AddTrackCuts(cutSetLambda);
   //
   /////////////////////////////////////////////////
   // selections for AntiLambda
   AliRsnCutV0 *cutAntiLambda = new AliRsnCutV0("cutAntiLambda", kLambda0Bar, AliPID::kProton, AliPID::kPion);
   cutAntiLambda->SetPIDCutProton(pPIDCut);
   cutAntiLambda->SetPIDCutPion(piPIDCut);
   cutAntiLambda->SetESDtrackCuts(esdTrackCuts);
   cutAntiLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
   cutAntiLambda->SetMaxDCAVertex(lambdaDCA);
   cutAntiLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
   cutAntiLambda->SetTolerance(massTol);
   cutAntiLambda->SetSwitch(false);
   cutAntiLambda->SetfLife(pLife); 
   cutAntiLambda->SetfLowRadius(radiuslow); 
   cutAntiLambda->SetfHighRadius(radiushigh); 
   cutAntiLambda->SetMaxRapidity(2.);
   cutAntiLambda->SetMinTPCcluster(NTPCcluster);
   // 
   AliRsnCutSet *cutSetAntiLambda = new AliRsnCutSet("setAntiLambda", AliRsnTarget::kDaughter);
   cutSetAntiLambda->AddCut(cutAntiLambda);
   cutSetAntiLambda->SetCutScheme(cutAntiLambda->GetName());
   Int_t iCutAntiLambda = task->AddTrackCuts(cutSetAntiLambda); 
   //
   /////////////////////////////////////////////////
    
  if(enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetPi->GetMonitorOutput(), monitorOpt.Data());
  } 
   
   
   //######################################################################################################  
    
   //
   // -- Values ------------------------------------------------------------------------------------
   //
   
   /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
   /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
   /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
   
   //
   // -- Create all needed outputs -----------------------------------------------------------------
   //
   
   // use an array for more compact writing, which are different on mixing and charges
   // [0] = unlike
   // [1] = mixing
   // [2] = like ++
   // [3] = like --
   Bool_t   use     [18] = { 1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1          ,  1              ,  1              ,  1              ,  1              ,  1              };
   Bool_t   useIM   [18] = { 1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1          ,  1              ,  1              ,  1              ,  1              ,  1              };
   TString  name    [18] = {"LambdapPip"   , "LambdapPim"   , "LambdaaPim"      , "LambdaaPip"      , "LambdapPipMix", "LambdapPimMix", "LambdaaPimMix"   , "LambdaaPipMix"   , "SigmaPt"  , "SigmaMt"  , "ASigmaPt"     , "ASigmaMt"     , "XiM"       , "XiP"           , "Lambda1520P"   , "Lambda1520M"   , "Lambda1520PBar", "Lambda1520MBar"};
   TString  comp    [18] = {"PAIR"     , "PAIR"     , "PAIR"         , "PAIR"         , "MIX"      , "MIX"      , "MIX"          , "MIX"          , "TRUE"     , "TRUE"     , "TRUE"         , "TRUE"         , "TRUE"      , "TRUE"          , "TRUE"          , "TRUE"          , "TRUE"          , "TRUE"          };
   TString  output  [18] = {"SPARSE"     , "SPARSE"     , "SPARSE"         , "SPARSE"         , "SPARSE"     , "SPARSE"     , "SPARSE"         , "SPARSE"         , "SPARSE"     , "SPARSE"     , "SPARSE"         , "SPARSE"         , "SPARSE"      , "SPARSE"          , "SPARSE"          , "SPARSE"          , "SPARSE"          , "SPARSE"          };
   Char_t   charge1 [18] = {'0'        , '0'        , '0'            , '0'            , '0'        , '0'        , '0'            , '0'            , '0'        , '0'        , '0'            , '0'            , '0'         , '0'             , '0'             , '0'             , '0'             , '0'             };
   Char_t   charge2 [18] = {'+'        , '-'        , '-'            , '+'            , '+'        , '-'        , '-'            , '+'            , '+'        , '-'        , '-'            , '+'            , '-'         , '+'             , '+'             , '-'             , '-'             , '+'             };
   Int_t    cutID1  [18] = { iCutLambda,  iCutLambda,  iCutAntiLambda,  iCutAntiLambda,  iCutLambda,  iCutLambda,  iCutAntiLambda,  iCutAntiLambda,  iCutLambda,  iCutLambda,  iCutAntiLambda,  iCutAntiLambda,  iCutLambda ,  iCutAntiLambda ,  iCutLambda     ,  iCutLambda     ,  iCutAntiLambda ,  iCutAntiLambda };
   Int_t    cutID2  [18] = { iCutPi    ,  iCutPi    ,  iCutPi        ,  iCutPi        ,  iCutPi    ,  iCutPi    ,  iCutPi        ,  iCutPi        ,  iCutPi    ,  iCutPi    ,  iCutPi        ,  iCutPi        ,  iCutPi     ,  iCutPi         ,  iCutPi         ,  iCutPi         ,  iCutPi         ,  iCutPi         };
   Int_t    ipdg    [18] = { 3224      ,  3114      , -3224          , -3114          ,  3224      ,  3114      , -3224          , -3114          ,  3224      ,  3114      , -3224          , -3114          ,  3312       , -3312           ,  3124           ,  3124           , -3124           , -3124           };
   Double_t mass    [18] = { 1.3828    ,  1.3872    ,  1.3828        ,  1.3872        ,  1.3828    ,  1.3872    ,  1.3828        ,  1.3872        ,  1.3828    ,  1.3872    ,  1.3828        ,  1.3872        ,  1.32171    ,  1.32171        ,  1.5195         ,  1.5195         ,  1.5195         ,  1.5195         };
   
   for (Int_t i = 0; i < 18; i++) {
      if (!use[i]) continue;
      if (collSyst) output[i] = "SPARSE";
      // create output
      AliRsnMiniOutput *out = task->CreateOutput(Form("Lambdapi_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
      // selection settings
      out->SetCutID(0, cutID1[i]);
      out->SetCutID(1, cutID2[i]);
      out->SetDaughter(0, AliRsnDaughter::kLambda);
      out->SetDaughter(1, AliRsnDaughter::kPion);
      out->SetCharge(0, charge1[i]);
      out->SetCharge(1, charge2[i]);
      out->SetMotherPDG(ipdg[i]);
      out->SetMotherMass(mass[i]);
      // pair cuts
      out->SetPairCuts(cutsPair);
      // axis X: invmass
      if (useIM[i]) 
             out->AddAxis(imID, 875, 1.25, 3.);
      // axis Y: transverse momentum
	  out->AddAxis(ptID, 100, 0.0, 10.0);
	 //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
	 
	  //if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
      
    }

   TString pname="pi";
   AddMonitorOutput_Pt(pname,cutSetPi->GetMonitorOutput());
   AddMonitorOutput_Eta(pname,cutSetPi->GetMonitorOutput());
   AddMonitorOutput_DCAxy(pname,cutSetPi->GetMonitorOutput());
   AddMonitorOutput_DCAz(pname,cutSetPi->GetMonitorOutput());
   AddMonitorOutput_TPCpi(pname,cutSetPi->GetMonitorOutput());
   AddMonitorOutput_NclTPC(pname,cutSetPi->GetMonitorOutput());
   AddMonitorOutput_chi2TPC(pname,cutSetPi->GetMonitorOutput());

   pname.Form("lambdap");
   AddMonitorOutput_P(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_Pt(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0NPt(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0PPt(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0Mass(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0DCA(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0Radius(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0Lifetime(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0DaughterDCA(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0CPA(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0DCA2TPV(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0TPCpim(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaProtonPID(cutSetLambda->GetMonitorOutput());

   pname.Form("lambdaa");
   AddMonitorOutput_P(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_Pt(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0NPt(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0PPt(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0Mass(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0DCA(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0Radius(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0Lifetime(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0DaughterDCA(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0CPA(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0DCA2TPV(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0TPCpip(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaAntiProtonPID(cutSetAntiLambda->GetMonitorOutput());

   if (isMC) {
     TString mode = "HIST";
     if (collSyst) mode = "SPARSE";
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarP_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetMotherPDG(3224);
     out->SetMotherMass(1.3828);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarM_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetMotherPDG(3114);
     out->SetMotherMass(1.3872);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarPBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetMotherPDG(-3224);
     out->SetMotherMass(1.3828);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarMBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetMotherPDG(-3114);
     out->SetMotherMass(1.3872);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("XiP_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetMotherPDG(-3312);
     out->SetMotherMass(1.32171);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("XiM_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetMotherPDG(3312);
     out->SetMotherMass(1.32171);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520P_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetCharge(0, 0);
     out->SetCharge(1, 1);
     out->SetMotherPDG(3124);
     out->SetMotherMass(1.5195);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520M_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetCharge(0, 0);
     out->SetCharge(1, -1);
     out->SetMotherPDG(3124);
     out->SetMotherMass(1.5195);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520PBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetCharge(0, 0);
     out->SetCharge(1, 1);
     out->SetMotherPDG(-3124);
     out->SetMotherMass(1.5195);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520MBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetCharge(0, 0);
     out->SetCharge(1, -1);
     out->SetMotherPDG(-3124);
     out->SetMotherMass(1.5195);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     
   }
   
   return kTRUE;
}


//=============================


AliRsnMiniAnalysisTask* AddTask_Lambdakx(TString lname,Bool_t isMC,Int_t system,Int_t EventCuts,Int_t TrackCutsLambda,Int_t TrackCutsK){
  Bool_t      isPP = true;
  if(system) isPP = false;
  Float_t     cutV = 10.0;
  Int_t       aodFilterBit = 5;
  Float_t     piPIDCut = 3.0;
  Float_t     pPIDCut = 3.0;
  Float_t     trackDCAcut = 7.0;
  Float_t     massTol = 0.006;
  Float_t     lambdaDCA = 0.3;
  Float_t     lambdaCosPoinAn = 0.99;
  Float_t     lambdaDaughDCA = 0.5;
  Int_t       NTPCcluster = 70;
  Int_t       nmix = 5;
  Float_t     maxDiffVzMix = 1.0;
  Float_t     maxDiffMultMix = 2.0;
  Float_t     maxDiffAngleMixDeg = 20.0;
  Int_t       aodN = 0;
  TString     outNameSuffix = lname;
     
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskRare_pp13::AddTask_Lambdapi", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   TString taskName = Form("%s_%s%s", lname.Data(), (isPP? "pp" : "PbPb"), (isMC ? "MC" : "Data"));
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   if (!isMC && !isPP){
     Printf(Form("========== SETTING USE CENTRALITY PATCH AOD049 : %s", (aodN==49)? "yes" : "no"));
     task->SetUseCentralityPatch(aodN==49);
   }
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
   ::Info("AddTaskRare_pp13::AddTask_Lambdakx", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n max diff EP angle = %5.3f deg", nmix, maxDiffVzMix, maxDiffMultMix, (isPP ? 0.0 : maxDiffAngleMixDeg)));

   int trigger=EventCuts%10;
   if(!trigger) task->UseESDTriggerMask(AliVEvent::kINT7);
   else if(trigger==1) task->UseESDTriggerMask(AliVEvent::kHighMult);
   else if(trigger==2) task->UseESDTriggerMask(AliVEvent::kHighMultV0);
   else task->UseESDTriggerMask(AliVEvent::kAny);
   
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
   AliRsnMiniOutput *outPlane = 0x0; //task->CreateOutput("eventPlane", "HIST", "EVENT");
   if (!isPP){
     outPlane = task->CreateOutput("eventPlane", "HIST", "EVENT");
     outPlane->AddAxis(planeID, 180, 0.0, TMath::Pi());
     }
   
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
   //gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigSigmaStar.C");
   if (isMC) {
       Printf("========================== MC analysis - PID cuts not used");
   } else 
     Printf("========================== DATA analysis - PID cuts used");
   if (!Config_Lambdakx(task, system, isMC, TrackCutsLambda, TrackCutsK, piPIDCut, pPIDCut, aodFilterBit, trackDCAcut, massTol, lambdaDCA, lambdaCosPoinAn, lambdaDaughDCA, NTPCcluster, "", cutsPair)) return 0x0;
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddTaskRare_pp13::AddTask_Lambdakx - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",lname.Data(),), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}


Bool_t Config_Lambdakx
(  
   AliRsnMiniAnalysisTask *task,
   Int_t		   collSyst, 
   Bool_t                  isMC,
   Int_t                   TrackCutsLambda,
   Int_t                   TrackCutsK,
   Float_t                 piPIDCut,
   Float_t                 pPIDCut,
   Int_t                   aodFilterBit,
   Float_t                 trackDCAcut,
   Float_t                 massTol,
   Float_t                 lambdaDCA,
   Float_t                 lambdaCosPoinAn,
   Float_t                 lambdaDaughDCA,
   Int_t                   NTPCcluster,
   const char             *suffix,
   AliRsnCutSet           *cutsPair
)
{
  Float_t     KPIDCut = 2.;
  Float_t     nsigmaKxTOF=3.;
  Float_t     pLife = 20;
  Float_t     radiuslow = 0.5;
  Float_t     radiushigh = 200;
  Bool_t      enableMonitor=kTRUE;
  TString     monitorOpt="pp";

   // manage suffix
   if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
 
   /////////////////////////////////////////////////////
   // selections for the pion from the decay of Sigma*
   /////////////////////////////////////////////////////
   //

   AliRsnCutTrackQuality* trkQualityCut= new AliRsnCutTrackQuality("myQualityCut");
   trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
   AliRsnCutSetDaughterParticle* cutSetKaon=new AliRsnCutSetDaughterParticle(Form("cutKaon%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,KPIDCut),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kKaon,KPIDCut,nsigmaKxTOF);
   Int_t iCutKaon = task->AddTrackCuts(cutSetKaon);

   /*
   AliRsnCutDaughterSigmaStar2010PP *cutKaon = new AliRsnCutDaughterSigmaStar2010PP("cutKaonForLambdakx", AliPID::kKaon);
   cutKaon->SetPIDCut(KPIDCut);    // fPIDCut used in IsSelected() after the call to cutQuality

   AliRsnCutTrackQuality *cutQuality = (AliRsnCutTrackQuality*) cutKaon->CutQuality();
   cutQuality->SetDefaults2011();
   //cutQuality->SetDCARmax(trackDCAcut);
   //cutQuality->SetDefaults2010(0,1);  // 1st par. not default (0 -> use TPC clusters). 2nd par. default (-> standard Pt and eta range)
   // SetDefaults2010 contains the following selections:
   //     SetPtRange(0.15, 1E+20);
   //     SetEtaRange(-0.8, 0.8);
   //     and from aliroot/master/src/ANALYSIS/ANALYSISalice/AliESDtrackCuts.cxx
   //     AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(1,0)
   //         esdTrackCuts->SetMinNClustersTPC(70);
   //         esdTrackCuts->SetMaxChi2PerClusterTPC(4);
   //         esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
   //         esdTrackCuts->SetRequireTPCRefit(kTRUE);
   //         esdTrackCuts->SetRequireITSRefit(kTRUE);
   //         esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
   //         esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");    // NB. With pt_min=0.15 (see above) -> DCAxy_max = 0.2560
   //         esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
   //         esdTrackCuts->SetMaxDCAToVertexZ(2);
   //         esdTrackCuts->SetDCAToVertex2D(kFALSE);
   //         esdTrackCuts->SetRequireSigmaToVertex(kFALSE);  
   //         esdTrackCuts->SetMaxChi2PerClusterITS(36);
   //  
   AliRsnCutSet *cutSetKaon = new AliRsnCutSet("setKaonForLambdakx", AliRsnTarget::kDaughter);
   cutSetKaon->AddCut(cutKaon);
   cutSetKaon->SetCutScheme(cutKaon->GetName());
   Int_t iCutKaon = task->AddTrackCuts(cutSetKaon);
   */
   //
   /////////////////////////////////////////////////////////////
   // selections for Lambda and for the daughters of Lambda 
   /////////////////////////////////////////////////////////////
   // 
   // selections for the proton and pion daugthers of Lambda and AntiLambda
   AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("qualityDaughterLambda");   
   esdTrackCuts->SetPtRange(0.15,1.E10);
   esdTrackCuts->SetEtaRange(-0.8,0.8);
   esdTrackCuts->SetRequireTPCRefit();
   esdTrackCuts->SetAcceptKinkDaughters(0);
   esdTrackCuts->SetMinNClustersTPC(NTPCcluster);
   esdTrackCuts->SetMaxChi2PerClusterTPC(4);
   esdTrackCuts->SetMinDCAToVertexXY(0.15);   
   //
   /////////////////////////////////////////////////
   // selections for Lambda
   AliRsnCutV0 *cutLambda = new AliRsnCutV0("cutLambda", kLambda0, AliPID::kProton, AliPID::kPion);
   cutLambda->SetPIDCutProton(pPIDCut);       // PID for the proton daughter of Lambda
   cutLambda->SetPIDCutPion(piPIDCut);        // PID for the pion daughter of Lambda 
   cutLambda->SetESDtrackCuts(esdTrackCuts);  // all the other selections (defined above) for proton and pion daughters of Lambda
   cutLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
   cutLambda->SetMaxDCAVertex(lambdaDCA);
   cutLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
   cutLambda->SetTolerance(massTol);
   cutLambda->SetSwitch(false);
   cutLambda->SetfLife(pLife); 
   cutLambda->SetfLowRadius(radiuslow); 
   cutLambda->SetfHighRadius(radiushigh); 
   cutLambda->SetMaxRapidity(2.);
   cutLambda->SetMinTPCcluster(NTPCcluster);
   //
   AliRsnCutSet *cutSetLambda = new AliRsnCutSet("setLambda", AliRsnTarget::kDaughter);
   cutSetLambda->AddCut(cutLambda);
   cutSetLambda->SetCutScheme(cutLambda->GetName());
   Int_t iCutLambda = task->AddTrackCuts(cutSetLambda);
   //
   /////////////////////////////////////////////////
   // selections for AntiLambda
   AliRsnCutV0 *cutAntiLambda = new AliRsnCutV0("cutAntiLambda", kLambda0Bar, AliPID::kProton, AliPID::kPion);
   cutAntiLambda->SetPIDCutProton(pPIDCut);
   cutAntiLambda->SetPIDCutPion(piPIDCut);
   cutAntiLambda->SetESDtrackCuts(esdTrackCuts);
   cutAntiLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
   cutAntiLambda->SetMaxDCAVertex(lambdaDCA);
   cutAntiLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
   cutAntiLambda->SetTolerance(massTol);
   cutAntiLambda->SetSwitch(false);
   cutAntiLambda->SetfLife(pLife); 
   cutAntiLambda->SetfLowRadius(radiuslow); 
   cutAntiLambda->SetfHighRadius(radiushigh); 
   cutAntiLambda->SetMaxRapidity(2.);
   cutAntiLambda->SetMinTPCcluster(NTPCcluster);
   // 
   AliRsnCutSet *cutSetAntiLambda = new AliRsnCutSet("setAntiLambda", AliRsnTarget::kDaughter);
   cutSetAntiLambda->AddCut(cutAntiLambda);
   cutSetAntiLambda->SetCutScheme(cutAntiLambda->GetName());
   Int_t iCutAntiLambda = task->AddTrackCuts(cutSetAntiLambda); 
   //
   /////////////////////////////////////////////////
    
  if(enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetKaon->GetMonitorOutput(), monitorOpt.Data());
  } 
   
   
   //######################################################################################################  
    
   //
   // -- Values ------------------------------------------------------------------------------------
   //
   
   /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
   /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
   /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
   
   //
   // -- Create all needed outputs -----------------------------------------------------------------
   //
   
   // use an array for more compact writing, which are different on mixing and charges
   // [0] = unlike
   // [1] = mixing
   // [2] = like ++
   // [3] = like --
   Bool_t   use     [18] = { 1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1          ,  1              ,  1              ,  1              ,  1              ,  1              };
   Bool_t   useIM   [18] = { 1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1          ,  1              ,  1              ,  1              ,  1              ,  1              };
   TString  name    [18] = {"LambdapKp"   , "LambdapKm"   , "LambdaaKm"      , "LambdaaKp"      , "LambdapKpMix", "LambdapKmMix", "LambdaaKmMix"   , "LambdaaKpMix"   , "SigmaPt"  , "SigmaMt"  , "ASigmaPt"     , "ASigmaMt"     , "XiM"       , "XiP"           , "Lambda1520P"   , "Lambda1520M"   , "Lambda1520PBar", "Lambda1520MBar"};
   TString  comp    [18] = {"PAIR"     , "PAIR"     , "PAIR"         , "PAIR"         , "MIX"      , "MIX"      , "MIX"          , "MIX"          , "TRUE"     , "TRUE"     , "TRUE"         , "TRUE"         , "TRUE"      , "TRUE"          , "TRUE"          , "TRUE"          , "TRUE"          , "TRUE"          };
   TString  output  [18] = {"SPARSE"     , "SPARSE"     , "SPARSE"         , "SPARSE"         , "SPARSE"     , "SPARSE"     , "SPARSE"         , "SPARSE"         , "SPARSE"     , "SPARSE"     , "SPARSE"         , "SPARSE"         , "SPARSE"      , "SPARSE"          , "SPARSE"          , "SPARSE"          , "SPARSE"          , "SPARSE"          };
   Char_t   charge1 [18] = {'0'        , '0'        , '0'            , '0'            , '0'        , '0'        , '0'            , '0'            , '0'        , '0'        , '0'            , '0'            , '0'         , '0'             , '0'             , '0'             , '0'             , '0'             };
   Char_t   charge2 [18] = {'+'        , '-'        , '-'            , '+'            , '+'        , '-'        , '-'            , '+'            , '+'        , '-'        , '-'            , '+'            , '-'         , '+'             , '+'             , '-'             , '-'             , '+'             };
   Int_t    cutID1  [18] = { iCutLambda,  iCutLambda,  iCutAntiLambda,  iCutAntiLambda,  iCutLambda,  iCutLambda,  iCutAntiLambda,  iCutAntiLambda,  iCutLambda,  iCutLambda,  iCutAntiLambda,  iCutAntiLambda,  iCutLambda ,  iCutAntiLambda ,  iCutLambda     ,  iCutLambda     ,  iCutAntiLambda ,  iCutAntiLambda };
   Int_t    cutID2  [18] = { iCutKaon  ,  iCutKaon  ,  iCutKaon      ,  iCutKaon      ,  iCutKaon  ,  iCutKaon  ,  iCutKaon      ,  iCutKaon      ,  iCutKaon  ,  iCutKaon  ,  iCutKaon      ,  iCutKaon      ,  iCutKaon   ,  iCutKaon       ,  iCutKaon       ,  iCutKaon       ,  iCutKaon       ,  iCutKaon       };
   Int_t    ipdg    [18] = { 3224      ,  3114      , -3224          , -3114          ,  3224      ,  3114      , -3224          , -3114          ,  3224      ,  3114      , -3224          , -3114          ,  3312       , -3312           ,  3124           ,  3124           , -3124           , -3124           };
   Double_t mass    [18] = { 1.3828    ,  1.3872    ,  1.3828        ,  1.3872        ,  1.3828    ,  1.3872    ,  1.3828        ,  1.3872        ,  1.3828    ,  1.3872    ,  1.3828        ,  1.3872        ,  1.32171    ,  1.32171        ,  1.5195         ,  1.5195         ,  1.5195         ,  1.5195         };
   
   for (Int_t i = 0; i < 18; i++) {
      if (!use[i]) continue;
      if (collSyst) output[i] = "SPARSE";
      // create output
      AliRsnMiniOutput *out = task->CreateOutput(Form("Lambdakx_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
      // selection settings
      out->SetCutID(0, cutID1[i]);
      out->SetCutID(1, cutID2[i]);
      out->SetDaughter(0, AliRsnDaughter::kLambda);
      out->SetDaughter(1, AliRsnDaughter::kKaon);
      out->SetCharge(0, charge1[i]);
      out->SetCharge(1, charge2[i]);
      out->SetMotherPDG(ipdg[i]);
      out->SetMotherMass(mass[i]);
      // pair cuts
      out->SetPairCuts(cutsPair);
      // axis X: invmass
      if (useIM[i]) 
             out->AddAxis(imID, 700, 1.6, 3.);
      // axis Y: transverse momentum
	  out->AddAxis(ptID, 100, 0.0, 10.0);
	 //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
	 
	  //if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
      
    }

   TString pname="kx";
   AddMonitorOutput_Pt(pname,cutSetKaon->GetMonitorOutput());
   AddMonitorOutput_Eta(pname,cutSetKaon->GetMonitorOutput());
   AddMonitorOutput_DCAxy(pname,cutSetKaon->GetMonitorOutput());
   AddMonitorOutput_DCAz(pname,cutSetKaon->GetMonitorOutput());
   AddMonitorOutput_TPCpi(pname,cutSetKaon->GetMonitorOutput());
   AddMonitorOutput_NclTPC(pname,cutSetKaon->GetMonitorOutput());
   AddMonitorOutput_chi2TPC(pname,cutSetKaon->GetMonitorOutput());

   pname.Form("lambdap");
   AddMonitorOutput_P(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_Pt(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0NPt(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0PPt(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0Mass(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0DCA(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0Radius(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0Lifetime(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0DaughterDCA(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0CPA(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0DCA2TPV(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0TPCpim(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaProtonPID(cutSetLambda->GetMonitorOutput());

   pname.Form("lambdaa");
   AddMonitorOutput_P(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_Pt(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0NPt(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0PPt(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0Mass(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0DCA(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0Radius(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0Lifetime(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0DaughterDCA(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0CPA(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0DCA2TPV(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0TPCpip(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaAntiProtonPID(cutSetAntiLambda->GetMonitorOutput());

   if (isMC) {
     TString mode = "HIST";
     if (collSyst) mode = "SPARSE";
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarP_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetMotherPDG(3224);
     out->SetMotherMass(1.3828);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarM_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetMotherPDG(3114);
     out->SetMotherMass(1.3872);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarPBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetMotherPDG(-3224);
     out->SetMotherMass(1.3828);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarMBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetMotherPDG(-3114);
     out->SetMotherMass(1.3872);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("XiP_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetMotherPDG(-3312);
     out->SetMotherMass(1.32171);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("XiM_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetMotherPDG(3312);
     out->SetMotherMass(1.32171);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520P_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetCharge(0, 0);
     out->SetCharge(1, 1);
     out->SetMotherPDG(3124);
     out->SetMotherMass(1.5195);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520M_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetCharge(0, 0);
     out->SetCharge(1, -1);
     out->SetMotherPDG(3124);
     out->SetMotherMass(1.5195);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520PBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetCharge(0, 0);
     out->SetCharge(1, 1);
     out->SetMotherPDG(-3124);
     out->SetMotherMass(1.5195);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520MBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetCharge(0, 0);     out->SetCharge(1, -1);
     out->SetMotherPDG(-3124);
     out->SetMotherMass(1.5195);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
   }

   return kTRUE;
}


//=============================


AliRsnMiniAnalysisTask* AddTask_Lambdak0(TString lname,Bool_t isMC,Int_t system,Int_t EventCuts,Int_t TrackCutsLambda,Int_t TrackCutsK){
  Bool_t      isPP = true;
  if(system) isPP = false;
  Float_t     cutV = 10.0;
  Int_t       aodFilterBit = 5;
  Float_t     piPIDCut = 3.0;
  Float_t     pPIDCut = 3.0;
  Float_t     trackDCAcut = 7.0;
  Float_t     massTol = 0.006;
  Float_t     lambdaDCA = 0.3;
  Float_t     lambdaCosPoinAn = 0.99;
  Float_t     lambdaDaughDCA = 0.5;
  Int_t       NTPCcluster = 70;
  Int_t       nmix = 5;
  Float_t     maxDiffVzMix = 1.0;
  Float_t     maxDiffMultMix = 2.0;
  Float_t     maxDiffAngleMixDeg = 20.0;
  Int_t       aodN = 0;
  TString     outNameSuffix = lname;
     
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskRare_pp13::AddTask_Lambdak0", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   TString taskName = Form("%s_%s%s", lname.Data(), (isPP? "pp" : "PbPb"), (isMC ? "MC" : "Data"));
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   if (!isMC && !isPP){
     Printf(Form("========== SETTING USE CENTRALITY PATCH AOD049 : %s", (aodN==49)? "yes" : "no"));
     task->SetUseCentralityPatch(aodN==49);
   }
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
   ::Info("AddTaskRare_pp13::AddTask_Lambdak0", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n max diff EP angle = %5.3f deg", nmix, maxDiffVzMix, maxDiffMultMix, (isPP ? 0.0 : maxDiffAngleMixDeg)));

   int trigger=EventCuts%10;
   if(!trigger) task->UseESDTriggerMask(AliVEvent::kINT7);
   else if(trigger==1) task->UseESDTriggerMask(AliVEvent::kHighMult);
   else if(trigger==2) task->UseESDTriggerMask(AliVEvent::kHighMultV0);
   else task->UseESDTriggerMask(AliVEvent::kAny);
   
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
   AliRsnMiniOutput *outPlane = 0x0; //task->CreateOutput("eventPlane", "HIST", "EVENT");
   if (!isPP){
     outPlane = task->CreateOutput("eventPlane", "HIST", "EVENT");
     outPlane->AddAxis(planeID, 180, 0.0, TMath::Pi());
     }
   
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
   //gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigSigmaStar.C");
   if (isMC) {
       Printf("========================== MC analysis - PID cuts not used");
   } else 
     Printf("========================== DATA analysis - PID cuts used");
   if (!Config_Lambdak0(task, system, isMC, TrackCutsLambda, TrackCutsK, piPIDCut, pPIDCut, aodFilterBit, trackDCAcut, massTol, lambdaDCA, lambdaCosPoinAn, lambdaDaughDCA, NTPCcluster, "", cutsPair)) return 0x0;
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddTaskRare_pp13::AddTask_Lambdak0 - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",lname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}


Bool_t Config_Lambdak0
(  
   AliRsnMiniAnalysisTask *task,
   Int_t		   collSyst, 
   Bool_t                  isMC,
   Int_t                   TrackCutsLambda,
   Int_t                   TrackCutsK,
   Float_t                 piPIDCut,
   Float_t                 pPIDCut,
   Int_t                   aodFilterBit,
   Float_t                 trackDCAcut,
   Float_t                 massTol,
   Float_t                 lambdaDCA,
   Float_t                 lambdaCosPoinAn,
   Float_t                 lambdaDaughDCA,
   Int_t                   NTPCcluster,
   const char             *suffix,
   AliRsnCutSet           *cutsPair
)
{
  Float_t     massTolVeto = 0.004;
  Float_t     pLife = 20;
  Float_t     radiuslow = 0.5;
  Float_t     radiushigh = 200;
  Bool_t      enableMonitor=kTRUE;
  TString     monitorOpt="pp";
  Bool_t      ptDep=true;
  Float_t     pi_k0s_PIDCut=3.;

  Float_t     k0s_massTol = 0.03;
  Float_t     k0s_massTolVeto = 0.004;
  Float_t     pLife = 20;
  Float_t     radiuslow = 0.5;
  Float_t     radiushigh = 200;   
  Bool_t      Switch = kFALSE;
  Float_t     k0sDCA = 0.3;
  Float_t     k0sCosPoinAn = 0.97;
  Float_t     k0sDaughDCA = 1.0;

   // manage suffix
   if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
 
   /////////////////////////////////////////////////////
   // selections for the K0S
   /////////////////////////////////////////////////////

   // selections for pion daugthers of K0s
   AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("qualityDaughterK0s");   
   esdTrackCuts->SetPtRange(0.15,1.E10);
   esdTrackCuts->SetEtaRange(-0.8,0.8);
   esdTrackCuts->SetRequireTPCRefit();
   esdTrackCuts->SetAcceptKinkDaughters(0); //
   esdTrackCuts->SetMinNClustersTPC(NTPCcluster);
   esdTrackCuts->SetMaxChi2PerClusterTPC(4);
    
    if(ptDep){
   esdTrackCuts->SetMinDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    }else
   esdTrackCuts->SetMinDCAToVertexXY("0.06"); //Use one of the two - pt dependent or fixed value cut.
  
   //
   /////////////////////////////////////////////////
   // selections for K0s
   AliRsnCutV0 *cutK0s = new AliRsnCutV0("cutK0s", kK0Short, AliPID::kPion, AliPID::kPion);
   cutK0s->SetPIDCutPion(pi_k0s_PIDCut);        // PID for the pion daughter of K0s
   cutK0s->SetESDtrackCuts(esdTrackCuts);  // all the other selections (defined above) for proton and pion daughters of K0s
   cutK0s->SetMaxDaughtersDCA(k0sDaughDCA);
   cutK0s->SetMaxDCAVertex(k0sDCA);
   cutK0s->SetMinCosPointingAngle(k0sCosPoinAn);
   cutK0s->SetTolerance(k0s_massTol);
   cutK0s->SetToleranceVeto(k0s_massTolVeto);   //Rejection range for Competing V0 Rejection
   cutK0s->SetSwitch(Switch);    
   cutK0s->SetfLife(pLife); 
   cutK0s->SetfLowRadius(radiuslow); 
   cutK0s->SetfHighRadius(radiushigh); 
   cutK0s->SetMaxRapidity(2.0);
   //
   AliRsnCutSet *cutSetK0s = new AliRsnCutSet("setK0s", AliRsnTarget::kDaughter);
   cutSetK0s->AddCut(cutK0s);
   cutSetK0s->SetCutScheme(cutK0s->GetName());
   Int_t iCutK0s = task->AddTrackCuts(cutSetK0s);

   /////////////////////////////////////////////////////////////
   // selections for Lambda and for the daughters of Lambda 
   /////////////////////////////////////////////////////////////

   /////////////////////////////////////////////////
   // selections for Lambda
   AliRsnCutV0 *cutLambda = new AliRsnCutV0("cutLambda", kLambda0, AliPID::kProton, AliPID::kPion);
   cutLambda->SetPIDCutProton(pPIDCut);       // PID for the proton daughter of Lambda
   cutLambda->SetPIDCutPion(piPIDCut);        // PID for the pion daughter of Lambda 
   cutLambda->SetESDtrackCuts(esdTrackCuts);  // all the other selections (defined above) for proton and pion daughters of Lambda
   cutLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
   cutLambda->SetMaxDCAVertex(lambdaDCA);
   cutLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
   cutLambda->SetTolerance(massTol);
   cutLambda->SetToleranceVeto(massTolVeto);   //Rejection range for Competing V0 Rejection
   cutLambda->SetSwitch(Switch);
   cutLambda->SetfLife(pLife); 
   cutLambda->SetfLowRadius(radiuslow); 
   cutLambda->SetfHighRadius(radiushigh); 
   cutLambda->SetMaxRapidity(2.);
   cutLambda->SetMinTPCcluster(NTPCcluster);
   //
   AliRsnCutSet *cutSetLambda = new AliRsnCutSet("setLambda", AliRsnTarget::kDaughter);
   cutSetLambda->AddCut(cutLambda);
   cutSetLambda->SetCutScheme(cutLambda->GetName());
   Int_t iCutLambda = task->AddTrackCuts(cutSetLambda);
   //
   /////////////////////////////////////////////////
   // selections for AntiLambda
   AliRsnCutV0 *cutAntiLambda = new AliRsnCutV0("cutAntiLambda", kLambda0Bar, AliPID::kProton, AliPID::kPion);
   cutAntiLambda->SetPIDCutProton(pPIDCut);
   cutAntiLambda->SetPIDCutPion(piPIDCut);
   cutAntiLambda->SetESDtrackCuts(esdTrackCuts);
   cutAntiLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
   cutAntiLambda->SetMaxDCAVertex(lambdaDCA);
   cutAntiLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
   cutAntiLambda->SetTolerance(massTol);
   cutAntiLambda->SetToleranceVeto(massTolVeto);   //Rejection range for Competing V0 Rejection
   cutAntiLambda->SetSwitch(Switch);
   cutAntiLambda->SetfLife(pLife); 
   cutAntiLambda->SetfLowRadius(radiuslow); 
   cutAntiLambda->SetfHighRadius(radiushigh); 
   cutAntiLambda->SetMaxRapidity(2.);
   cutAntiLambda->SetMinTPCcluster(NTPCcluster);
   // 
   AliRsnCutSet *cutSetAntiLambda = new AliRsnCutSet("setAntiLambda", AliRsnTarget::kDaughter);
   cutSetAntiLambda->AddCut(cutAntiLambda);
   cutSetAntiLambda->SetCutScheme(cutAntiLambda->GetName());
   Int_t iCutAntiLambda = task->AddTrackCuts(cutSetAntiLambda); 
   //
   /////////////////////////////////////////////////
   
   
   //######################################################################################################  
    
   //
   // -- Values ------------------------------------------------------------------------------------
   //
   
   /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
   /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
   /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
   
   //
   // -- Create all needed outputs -----------------------------------------------------------------
   //
   
   // use an array for more compact writing, which are different on mixing and charges
   // [0] = unlike
   // [1] = mixing
   // [2] = like ++
   // [3] = like --
   Bool_t   use     [4] = { 1         ,  1         ,  1             ,  1};
   Bool_t   useIM   [4] = { 1         ,  1         ,  1             ,  1};
   TString  name    [4] = {"LambdapK0"   , "LambdaaK0"   , "LambdapK0Mix"      , "LambdaaK0Mix"};
   TString  comp    [4] = {"PAIR"     , "PAIR"     ,  "MIX"      , "MIX"};
   TString  output  [4] = {"SPARSE"     , "SPARSE"     , "SPARSE"         , "SPARSE"};
   Char_t   charge1 [4] = {'0'        , '0'        , '0'            , '0'};
   Char_t   charge2 [4] = {'0'        , '0'        , '0'            , '0'};
   Int_t    cutID1  [4] = { iCutLambda,  iCutAntiLambda,  iCutLambda, iCutAntiLambda};
   Int_t    cutID2  [4] = { iCutK0s  ,  iCutK0s  ,  iCutK0s      ,  iCutK0s};
   Int_t    ipdg    [4] = { 3224      ,  3114      , -3224          , -3114};
   Double_t mass    [4] = { 1.3828    ,  1.3872    ,  1.3828        ,  1.3872};
   
   for (Int_t i = 0; i < 4; i++) {
      if (!use[i]) continue;
      if (collSyst) output[i] = "SPARSE";
      // create output
      AliRsnMiniOutput *out = task->CreateOutput(Form("Lambdak0_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
      // selection settings
      out->SetCutID(0, cutID1[i]);
      out->SetCutID(1, cutID2[i]);
      out->SetDaughter(0, AliRsnDaughter::kLambda);
      out->SetDaughter(1, AliRsnDaughter::kKaon0);
      out->SetCharge(0, charge1[i]);
      out->SetCharge(1, charge2[i]);
      out->SetMotherPDG(ipdg[i]);
      out->SetMotherMass(mass[i]);
      // pair cuts
      out->SetPairCuts(cutsPair);
      // axis X: invmass
      if (useIM[i]) 
             out->AddAxis(imID, 700, 1.6, 3.);
      // axis Y: transverse momentum
	  out->AddAxis(ptID, 100, 0.0, 10.0);
	 //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
	 
	  //if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
      
    }

   TString pname="k0";
   AddMonitorOutput_P(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_Pt(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0NPt(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0PPt(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0Mass(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0DCA(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0Radius(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0Lifetime(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0DaughterDCA(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0CPA(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0DCA2TPV(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0TPCpim(pname,cutSetK0s->GetMonitorOutput());
   AddMonitorOutput_V0TPCpip(pname,cutSetK0s->GetMonitorOutput());

   pname.Form("lambdap");
   AddMonitorOutput_P(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_Pt(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0NPt(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0PPt(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0Mass(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0DCA(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0Radius(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0Lifetime(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0DaughterDCA(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0CPA(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0DCA2TPV(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0TPCpim(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaProtonPID(cutSetLambda->GetMonitorOutput());

   pname.Form("lambdaa");
   AddMonitorOutput_P(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_Pt(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0NPt(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0PPt(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0Mass(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0DCA(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0Radius(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0Lifetime(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0DaughterDCA(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0CPA(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0DCA2TPV(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0TPCpip(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaAntiProtonPID(cutSetAntiLambda->GetMonitorOutput());

   if (isMC) {
     TString mode = "HIST";
     if (collSyst) mode = "SPARSE";
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarP_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetMotherPDG(3224);
     out->SetMotherMass(1.3828);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarM_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetMotherPDG(3114);
     out->SetMotherMass(1.3872);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarPBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetMotherPDG(-3224);
     out->SetMotherMass(1.3828);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarMBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetMotherPDG(-3114);
     out->SetMotherMass(1.3872);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("XiP_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetMotherPDG(-3312);
     out->SetMotherMass(1.32171);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("XiM_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetMotherPDG(3312);
     out->SetMotherMass(1.32171);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520P_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetCharge(0, 0);
     out->SetCharge(1, 1);
     out->SetMotherPDG(3124);
     out->SetMotherMass(1.5195);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520M_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetCharge(0, 0);
     out->SetCharge(1, -1);
     out->SetMotherPDG(3124);
     out->SetMotherMass(1.5195);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520PBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetCharge(0, 0);
     out->SetCharge(1, 1);
     out->SetMotherPDG(-3124);
     out->SetMotherMass(1.5195);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520MBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kKaon);
     out->SetCharge(0, 0);     out->SetCharge(1, -1);
     out->SetMotherPDG(-3124);
     out->SetMotherMass(1.5195);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
   }

   return kTRUE;
}


//=============================


AliRsnMiniAnalysisTask* AddTask_Lambdap(TString lname,Bool_t isMC,Int_t system,Int_t EventCuts,Int_t TrackCutsLambda,Int_t TrackCutsP){
  Bool_t      isPP = true;
  if(system) isPP = false;
  Float_t     cutV = 10.0;
  Int_t       aodFilterBit = 5;
  Float_t     piPIDCut = 3.0;
  Float_t     pPIDCut = 3.0;
  Float_t     trackDCAcut = 7.0;
  Float_t     massTol = 0.01;
  Float_t     lambdaDCA = 0.3;
  Float_t     lambdaCosPoinAn = 0.99;
  Float_t     lambdaDaughDCA = 0.5;
  Int_t       NTPCcluster = 70;
  Int_t       nmix = 5;
  Float_t     maxDiffVzMix = 1.0;
  Float_t     maxDiffMultMix = 10.0;
  Float_t     maxDiffAngleMixDeg = 20.0;
  Int_t       aodN = 0;
  TString     outNameSuffix = lname;
     
  //
  // -- INITIALIZATION ----------------------------------------------------------------------------
  // retrieve analysis manager
  //

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskRare_pp13::AddTask_Lambdap", "No analysis manager to connect to.");
      return NULL;
   } 

   // create the task and configure 
   TString taskName = Form("%s_%s%s", lname.Data(), (isPP? "pp" : "PbPb"), (isMC ? "MC" : "Data"));
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName.Data(), isMC);
   if (!isMC && !isPP){
     Printf(Form("========== SETTING USE CENTRALITY PATCH AOD049 : %s", (aodN==49)? "yes" : "no"));
     task->SetUseCentralityPatch(aodN==49);
   }
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
   ::Info("AddTaskRare_pp13::AddTask_Lambdap", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f \n max diff EP angle = %5.3f deg", nmix, maxDiffVzMix, maxDiffMultMix, (isPP ? 0.0 : maxDiffAngleMixDeg)));

   int trigger=EventCuts%10;
   if(!trigger) task->UseESDTriggerMask(AliVEvent::kINT7);
   else if(trigger==1) task->UseESDTriggerMask(AliVEvent::kHighMult);
   else if(trigger==2) task->UseESDTriggerMask(AliVEvent::kHighMultV0);
   else task->UseESDTriggerMask(AliVEvent::kAny);
   
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
   AliRsnMiniOutput *outPlane = 0x0; //task->CreateOutput("eventPlane", "HIST", "EVENT");
   if (!isPP){
     outPlane = task->CreateOutput("eventPlane", "HIST", "EVENT");
     outPlane->AddAxis(planeID, 180, 0.0, TMath::Pi());
     }
   
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
   //gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/ConfigSigmaStar.C");
   if (isMC) {
       Printf("========================== MC analysis - PID cuts not used");
   } else 
     Printf("========================== DATA analysis - PID cuts used");
   if (!Config_Lambdap(task, system, isMC, TrackCutsLambda, TrackCutsP, piPIDCut, pPIDCut, aodFilterBit, trackDCAcut, massTol, lambdaDCA, lambdaCosPoinAn, lambdaDaughDCA, NTPCcluster, "", cutsPair)) return 0x0;
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   //  outputFileName += ":Rsn";
   Printf("AddTaskRare_pp13::AddTask_Lambdap - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",lname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);
   
   return task;
}


Bool_t Config_Lambdap
(  
   AliRsnMiniAnalysisTask *task,
   Int_t		   collSyst, 
   Bool_t                  isMC,
   Int_t                   TrackCutsLambda,
   Int_t                   TrackCutsP,
   Float_t                 piPIDCut,
   Float_t                 pPIDCut,
   Int_t                   aodFilterBit,
   Float_t                 trackDCAcut,
   Float_t                 massTol,
   Float_t                 lambdaDCA,
   Float_t                 lambdaCosPoinAn,
   Float_t                 lambdaDaughDCA,
   Int_t                   NTPCcluster,
   const char             *suffix,
   AliRsnCutSet           *cutsPair
)
{
  Float_t     pLife = 20;
  Float_t     radiuslow = 0.5;
  Float_t     radiushigh = 200;
  Bool_t      enableMonitor=kTRUE;
  TString     monitorOpt="pp";
  Float_t     p1PIDCut=2.;
  Int_t       customQualityCutsID=0;
  AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutPrCandidate = AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV;

   // manage suffix
   if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
 
   /////////////////////////////////////////////////////
   // selections for the proton
   /////////////////////////////////////////////////////

  AliRsnCutSetDaughterParticle* cutSetQ;
  AliRsnCutSetDaughterParticle* cutSetP;

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  if(SetCustomQualityCut(trkQualityCut,customQualityCutsID,aodFilterBit)){
    //Set custom quality cuts for systematic checks
    cutSetQ=new AliRsnCutSetDaughterParticle(Form("cutQuality_bit%i",aodFilterBit),trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2011,AliPID::kPion,-1.);
    cutSetP=new AliRsnCutSetDaughterParticle(Form("cutProton_%i_%2.1fsigma",cutPrCandidate, p1PIDCut),trkQualityCut,cutPrCandidate,AliPID::kProton,p1PIDCut);
  }else{
    Bool_t useCrossedRows = 1;
    cutSetQ=new AliRsnCutSetDaughterParticle(Form("cutQ_bit%i",aodFilterBit),AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.,aodFilterBit,useCrossedRows);
    cutSetP=new AliRsnCutSetDaughterParticle(Form("cutPr%i_%2.1fsigma",cutPrCandidate,p1PIDCut),cutPrCandidate,AliPID::kProton,p1PIDCut,aodFilterBit,useCrossedRows);
  }

  Int_t iCutQ = task->AddTrackCuts(cutSetQ);
  Int_t iCutP = task->AddTrackCuts(cutSetP);
  
  if (enableMonitor){
    Printf("======== Cut monitoring enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC, cutSetQ->GetMonitorOutput(), monitorOpt.Data());
    AddMonitorOutput(isMC, cutSetP->GetMonitorOutput(), monitorOpt.Data());
  }

   /////////////////////////////////////////////////////////////
   // selections for Lambda and for the daughters of Lambda 
   /////////////////////////////////////////////////////////////
   // 
   // selections for the proton and pion daugthers of Lambda and AntiLambda
   AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("qualityDaughterLambda");   
   esdTrackCuts->SetPtRange(0.15,1.E10);
   esdTrackCuts->SetEtaRange(-0.8,0.8);
   esdTrackCuts->SetRequireTPCRefit();
   esdTrackCuts->SetAcceptKinkDaughters(0);
   esdTrackCuts->SetMinNClustersTPC(NTPCcluster);
   esdTrackCuts->SetMaxChi2PerClusterTPC(4);
   esdTrackCuts->SetMinDCAToVertexXY(0.15);   
   //
   /////////////////////////////////////////////////
   // selections for Lambda
   AliRsnCutV0 *cutLambda = new AliRsnCutV0("cutLambda", kLambda0, AliPID::kProton, AliPID::kPion);
   cutLambda->SetPIDCutProton(pPIDCut);       // PID for the proton daughter of Lambda
   cutLambda->SetPIDCutPion(piPIDCut);        // PID for the pion daughter of Lambda 
   cutLambda->SetESDtrackCuts(esdTrackCuts);  // all the other selections (defined above) for proton and pion daughters of Lambda
   cutLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
   cutLambda->SetMaxDCAVertex(lambdaDCA);
   cutLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
   cutLambda->SetTolerance(massTol);
   cutLambda->SetSwitch(false);
   cutLambda->SetfLife(pLife); 
   cutLambda->SetfLowRadius(radiuslow); 
   cutLambda->SetfHighRadius(radiushigh); 
   cutLambda->SetMaxRapidity(2.);
   cutLambda->SetMinTPCcluster(NTPCcluster);
   //
   AliRsnCutSet *cutSetLambda = new AliRsnCutSet("setLambda", AliRsnTarget::kDaughter);
   cutSetLambda->AddCut(cutLambda);
   cutSetLambda->SetCutScheme(cutLambda->GetName());
   Int_t iCutLambda = task->AddTrackCuts(cutSetLambda);
   //
   /////////////////////////////////////////////////
   // selections for AntiLambda
   AliRsnCutV0 *cutAntiLambda = new AliRsnCutV0("cutAntiLambda", kLambda0Bar, AliPID::kProton, AliPID::kPion);
   cutAntiLambda->SetPIDCutProton(pPIDCut);
   cutAntiLambda->SetPIDCutPion(piPIDCut);
   cutAntiLambda->SetESDtrackCuts(esdTrackCuts);
   cutAntiLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
   cutAntiLambda->SetMaxDCAVertex(lambdaDCA);
   cutAntiLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
   cutAntiLambda->SetTolerance(massTol);
   cutAntiLambda->SetSwitch(false);
   cutAntiLambda->SetfLife(pLife); 
   cutAntiLambda->SetfLowRadius(radiuslow); 
   cutAntiLambda->SetfHighRadius(radiushigh); 
   cutAntiLambda->SetMaxRapidity(2.);
   cutAntiLambda->SetMinTPCcluster(NTPCcluster);
   // 
   AliRsnCutSet *cutSetAntiLambda = new AliRsnCutSet("setAntiLambda", AliRsnTarget::kDaughter);
   cutSetAntiLambda->AddCut(cutAntiLambda);
   cutSetAntiLambda->SetCutScheme(cutAntiLambda->GetName());
   Int_t iCutAntiLambda = task->AddTrackCuts(cutSetAntiLambda);

   //######################################################################################################  
    
   //
   // -- Values ------------------------------------------------------------------------------------
   //
   
   /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
   /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
   /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
   
   //
   // -- Create all needed outputs -----------------------------------------------------------------
   //
   
   // use an array for more compact writing, which are different on mixing and charges
   // [0] = unlike
   // [1] = mixing
   // [2] = like ++
   // [3] = like --
   Bool_t   use     [18] = { 1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1          ,  1              ,  1              ,  1              ,  1              ,  1              };
   Bool_t   useIM   [18] = { 1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1          ,  1              ,  1              ,  1              ,  1              ,  1              };
   TString  name    [18] = {"LambdapPp"   , "LambdapPm"   , "LambdaaPm"      , "LambdaaPp"      , "LambdapPpMix", "LambdapPmMix", "LambdaaPmMix"   , "LambdaaPpMix"   , "SigmaPt"  , "SigmaMt"  , "ASigmaPt"     , "ASigmaMt"     , "XiM"       , "XiP"           , "Lambda1520P"   , "Lambda1520M"   , "Lambda1520PBar", "Lambda1520MBar"};
   TString  comp    [18] = {"PAIR"     , "PAIR"     , "PAIR"         , "PAIR"         , "MIX"      , "MIX"      , "MIX"          , "MIX"          , "TRUE"     , "TRUE"     , "TRUE"         , "TRUE"         , "TRUE"      , "TRUE"          , "TRUE"          , "TRUE"          , "TRUE"          , "TRUE"          };
   TString  output  [18] = {"SPARSE"     , "SPARSE"     , "SPARSE"         , "SPARSE"         , "SPARSE"     , "SPARSE"     , "SPARSE"         , "SPARSE"         , "SPARSE"     , "SPARSE"     , "SPARSE"         , "SPARSE"         , "SPARSE"      , "SPARSE"          , "SPARSE"          , "SPARSE"          , "SPARSE"          , "SPARSE"          };
   Char_t   charge1 [18] = {'0'        , '0'        , '0'            , '0'            , '0'        , '0'        , '0'            , '0'            , '0'        , '0'        , '0'            , '0'            , '0'         , '0'             , '0'             , '0'             , '0'             , '0'             };
   Char_t   charge2 [18] = {'+'        , '-'        , '-'            , '+'            , '+'        , '-'        , '-'            , '+'            , '+'        , '-'        , '-'            , '+'            , '-'         , '+'             , '+'             , '-'             , '-'             , '+'             };
   Int_t    cutID1  [18] = { iCutLambda,  iCutLambda,  iCutAntiLambda,  iCutAntiLambda,  iCutLambda,  iCutLambda,  iCutAntiLambda,  iCutAntiLambda,  iCutLambda,  iCutLambda,  iCutAntiLambda,  iCutAntiLambda,  iCutLambda ,  iCutAntiLambda ,  iCutLambda     ,  iCutLambda     ,  iCutAntiLambda ,  iCutAntiLambda };
   Int_t    cutID2  [18] = { iCutP     ,  iCutP     ,  iCutP         ,  iCutP         ,  iCutP     ,  iCutP     ,  iCutP         ,  iCutP         ,  iCutP     ,  iCutP     ,  iCutP         ,  iCutP         ,  iCutP      ,  iCutP          ,  iCutP          ,  iCutP          ,  iCutP          ,  iCutP          };
   Int_t    ipdg    [18] = { 3224      ,  3114      , -3224          , -3114          ,  3224      ,  3114      , -3224          , -3114          ,  3224      ,  3114      , -3224          , -3114          ,  3312       , -3312           ,  3124           ,  3124           , -3124           , -3124           };
   Double_t mass    [18] = { 1.3828    ,  1.3872    ,  1.3828        ,  1.3872        ,  1.3828    ,  1.3872    ,  1.3828        ,  1.3872        ,  1.3828    ,  1.3872    ,  1.3828        ,  1.3872        ,  1.32171    ,  1.32171        ,  1.5195         ,  1.5195         ,  1.5195         ,  1.5195         };
   
   for (Int_t i = 0; i < 18; i++) {
      if (!use[i]) continue;
      if (collSyst) output[i] = "SPARSE";
      // create output
      AliRsnMiniOutput *out = task->CreateOutput(Form("Lambdap_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
      // selection settings
      out->SetCutID(0, cutID1[i]);
      out->SetCutID(1, cutID2[i]);
      out->SetDaughter(0, AliRsnDaughter::kLambda);
      out->SetDaughter(1, AliRsnDaughter::kProton);
      out->SetCharge(0, charge1[i]);
      out->SetCharge(1, charge2[i]);
      out->SetMotherPDG(ipdg[i]);
      out->SetMotherMass(mass[i]);
      // pair cuts
      out->SetPairCuts(cutsPair);
      // axis X: invmass
      if (useIM[i]) 
             out->AddAxis(imID, 975, 2.05, 4.);
      // axis Y: transverse momentum
	  out->AddAxis(ptID, 100, 0.0, 10.0);
	 //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
	 
	  //if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
      
    }

   TString pname="p";
   AddMonitorOutput_Pt(pname,cutSetP->GetMonitorOutput());
   AddMonitorOutput_Eta(pname,cutSetP->GetMonitorOutput());
   AddMonitorOutput_DCAxy(pname,cutSetP->GetMonitorOutput());
   AddMonitorOutput_DCAz(pname,cutSetP->GetMonitorOutput());
   AddMonitorOutput_TPCp(pname,cutSetP->GetMonitorOutput());
   AddMonitorOutput_NclTPC(pname,cutSetP->GetMonitorOutput());
   AddMonitorOutput_chi2TPC(pname,cutSetP->GetMonitorOutput());

   pname.Form("lambdap");
   AddMonitorOutput_P(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_Pt(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0NPt(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0PPt(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0Mass(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0DCA(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0Radius(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0Lifetime(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0DaughterDCA(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0CPA(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0DCA2TPV(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_V0TPCpim(pname,cutSetLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaProtonPID(cutSetLambda->GetMonitorOutput());

   pname.Form("lambdaa");
   AddMonitorOutput_P(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_Pt(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0NPt(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0PPt(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0Mass(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0DCA(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0Radius(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0Lifetime(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0DaughterDCA(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0CPA(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0DCA2TPV(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_V0TPCpip(pname,cutSetAntiLambda->GetMonitorOutput());
   AddMonitorOutput_LambdaAntiProtonPID(cutSetAntiLambda->GetMonitorOutput());

   if (isMC) {
     TString mode = "HIST";
     if (collSyst) mode = "SPARSE";
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarP_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetMotherPDG(3224);
     out->SetMotherMass(1.3828);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     //if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarM_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetMotherPDG(3114);
     out->SetMotherMass(1.3872);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarPBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetMotherPDG(-3224);
     out->SetMotherMass(1.3828);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("sigmastarMBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetMotherPDG(-3114);
     out->SetMotherMass(1.3872);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     //out->AddAxis(lambdaDCA, 10, 0.0, 1.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("XiP_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetMotherPDG(-3312);
     out->SetMotherMass(1.32171);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("XiM_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetMotherPDG(3312);
     out->SetMotherMass(1.32171);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520P_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetCharge(0, 0);
     out->SetCharge(1, 1);
     out->SetMotherPDG(3124);
     out->SetMotherMass(1.5195);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520M_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetCharge(0, 0);
     out->SetCharge(1, -1);
     out->SetMotherPDG(3124);
     out->SetMotherMass(1.5195);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     // create output
     AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520PBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetCharge(0, 0);
     out->SetCharge(1, 1);
     out->SetMotherPDG(-3124);
     out->SetMotherMass(1.5195);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     AliRsnMiniOutput *out = task->CreateOutput(Form("Lambda1520MBar_TrueMC%s", suffix), mode.Data(), "MOTHER");
     // selection settings
     out->SetDaughter(0, AliRsnDaughter::kLambda);
     out->SetDaughter(1, AliRsnDaughter::kPion);
     out->SetCharge(0, 0);
     out->SetCharge(1, -1);
     out->SetMotherPDG(-3124);
     out->SetMotherMass(1.5195);
     // pair cuts
     out->SetPairCuts(cutsPair);
     // binnings
     out->AddAxis(imID, 800, 1.2, 2.0);
     out->AddAxis(ptID, 100, 0.0, 10.0);
     
     if (collSyst) out->AddAxis(centID, 10, 0.0, 100.0);
     
     
     
   }
   
   return kTRUE;
}


//=============================


void AddMonitorOutput_P(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_mom",name.Data()),AliRsnValueDaughter::kP);
  a->SetBins(0.,10.0,0.05);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}


void AddMonitorOutput_Pt(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_pt",name.Data()),AliRsnValueDaughter::kPt);
  a->SetBins(0.,10.0,0.05);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_Eta(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_eta",name.Data()),AliRsnValueDaughter::kEta);
  a->SetBins(-2.,2.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_DCAxy(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_dcaxy",name.Data()),AliRsnValueDaughter::kDCAXY);
  a->SetBins(-0.5,0.5,0.001);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_DCAz(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_dcaz",name.Data()),AliRsnValueDaughter::kDCAZ);
  a->SetBins(-2.5,2.5,0.005);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_TPCpi(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_TPCpi",name.Data()),AliRsnValueDaughter::kTPCnsigmaPi);
  a->SetBins(-10.,10.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_TPCK(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_TPCK",name.Data()),AliRsnValueDaughter::kTPCnsigmaK);
  a->SetBins(-10.,10.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_TPCp(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_TPCp",name.Data()),AliRsnValueDaughter::kTPCnsigmaP);
  a->SetBins(-10.,10.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_NclTPC(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_NclTPC",name.Data()),AliRsnValueDaughter::kNTPCclusters);
  a->SetBins(-0.5,199.5,1);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_chi2TPC(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_chi2TPC",name.Data()),AliRsnValueDaughter::kTPCchi2);
  a->SetBins(0.0,6,.1);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0NPt(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0npt",name.Data()),AliRsnValueDaughter::kV0NPt);
  a->SetBins(0.,10.0,0.05);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0PPt(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0ppt",name.Data()),AliRsnValueDaughter::kV0PPt);
  a->SetBins(0.,10.0,0.05);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0Mass(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0mass",name.Data()),AliRsnValueDaughter::kV0Mass);
  name.ToLower();
  if(name.Contains("k0")) a->SetBins(0.4,0.6,0.001);
  else if(name.Contains("lambda")) a->SetBins(1.08,1.16,0.001);
  else a->SetBins(0.,3.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0DCA(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0dca",name.Data()),AliRsnValueDaughter::kV0DCA);
  a->SetBins(0.0,0.4,0.001);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0Radius(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0radius",name.Data()),AliRsnValueDaughter::kV0Radius);
  a->SetBins(0.0,200,0.2);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0Lifetime(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0lifetime",name.Data()),AliRsnValueDaughter::kV0Lifetime);
  a->SetBins(0.0,200,0.1);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0DaughterDCA(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0ddca",name.Data()),AliRsnValueDaughter::kDaughterDCA);
  a->SetBins(0.0,2,0.001);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0DCA2TPV(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{//DCA of secondary tracks to primary vertex
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0dca2tpv",name.Data()),AliRsnValueDaughter::kV0DCAXY);
  a->SetBins(-10.,10.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0CPA(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0cpa",name.Data()),AliRsnValueDaughter::kCosPointAng);
  a->SetBins(0.96,1.,0.001);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0TPCpim(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0TPCpim",name.Data()),AliRsnValueDaughter::kLambdaPionPIDCut);
  a->SetBins(0.,5.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0TPCpip(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0)
{
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0TPCpip",name.Data()),AliRsnValueDaughter::kAntiLambdaAntiPionPIDCut);
  a->SetBins(-0.,5.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_LambdaProtonPID(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lpPID=0)
{

   // Lambda Cosine of the Pointing Angle
   AliRsnValueDaughter *axisLambdaProtonPID = new AliRsnValueDaughter("lambda_protonPID", AliRsnValueDaughter::kLambdaProtonPIDCut);
   axisLambdaProtonPID->SetBins(0.0,5,0.01);

   // output: 2D histogram
   AliRsnListOutput *outMonitorLambdaProtonPID = new AliRsnListOutput("Lambda_ProtonPID", AliRsnListOutput::kHistoDefault);
   outMonitorLambdaProtonPID->AddValue(axisLambdaProtonPID);

   // add outputs to loop
   if (mon) mon->Add(outMonitorLambdaProtonPID);
   if (lpPID) lpPID->AddOutput(outMonitorLambdaProtonPID);
  
}

void AddMonitorOutput_LambdaAntiProtonPID(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lapPID=0)
{

   // Lambda Cosine of the Pointing Angle
   AliRsnValueDaughter *axisLambdaAntiProtonPID = new AliRsnValueDaughter("lambda_antiprotonPID", AliRsnValueDaughter::kAntiLambdaAntiProtonPIDCut);
   axisLambdaAntiProtonPID->SetBins(0.0,5,0.01);

   // output: 2D histogram
   AliRsnListOutput *outMonitorLambdaAntiProtonPID = new AliRsnListOutput("Lambda_AntiProtonPID", AliRsnListOutput::kHistoDefault);
   outMonitorLambdaAntiProtonPID->AddValue(axisLambdaAntiProtonPID);

   // add outputs to loop
   if (mon) mon->Add(outMonitorLambdaAntiProtonPID);
   if (lapPID) lapPID->AddOutput(outMonitorLambdaAntiProtonPID);
  
}
