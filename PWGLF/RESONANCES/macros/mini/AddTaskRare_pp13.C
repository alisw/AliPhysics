/***************************************************************************
              Anders Knospe: anders.knospe@cern.ch
  Macro to configure the resonance package for searches for rare resonances.

****************************************************************************/

AliRsnMiniAnalysisTask* AddTaskRare_pp13(
  TString lname,
  Bool_t isMC,
  Int_t system,
  AliRsnDaughter::ESpecies d1,
  AliRsnDaughter::ESpecies d2,
  Int_t EventCuts=0,
  Int_t TrackCuts1=0,
  Int_t TrackCuts2=0
){
  // ----- INITIALIZATION -----

  // retrieve analysis manager
  AliAnalysisManager* mgr=AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    ::Error("AddTaskRare_pp13", "No analysis manager to connect to.");
    return NULL;
  }

  // create the task and configure
  AliRsnMiniAnalysisTask* task=new AliRsnMiniAnalysisTask(lname,isMC);

  // trigger
  int trigger=EventCuts%10;
  if(!trigger) task->UseESDTriggerMask(AliVEvent::kINT7);
  else if(trigger==1) task->UseESDTriggerMask(AliVEvent::kHighMultV0);

  // multiplicity
  bool isPP=false;
  if(!system) isPP=true;
  int MultBins=(EventCuts/10)%10;
  if(system==1 || system==2) MultBins=1;

  if(isPP){
    if(MultBins==1) task->UseMultiplicity("AliMultSelection_V0M");
    else if(MultBins==2) task->UseMultiplicity("AliMultSelection_RefMult08");
    else task->UseMultiplicity("QUALITY");
  }else if(system==1) task->UseMultiplicity("AliMultSelection_V0A");
  else if(system==2) task->UseMultiplicity("AliMultSelection_V0M");
  else task->UseCentrality("V0M");

  // set event mixing options
  int nmix=5;
  if((EventCuts%10000)/1000==1) nmix=0;
  float maxDiffVzMix=1;
  float maxDiffMultMix=5;
  task->UseContinuousMix();
  task->SetNMix(nmix);
  task->SetMaxDiffVz(maxDiffVzMix);
  task->SetMaxDiffMult(maxDiffMultMix);
  ::Info("AddTaskRare_pp13", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f", nmix, maxDiffVzMix, maxDiffMultMix));

  // vertex cuts
  float vtxZcut=10;
  Bool_t rejectPileUp=kTRUE;
  AliRsnCutPrimaryVertex* cutVertex=0;
  if(!MultBins || fabs(vtxZcut-10.)>1.e-10){
    cutVertex=new AliRsnCutPrimaryVertex("cutVertex",vtxZcut,0,kFALSE);
    if(!MultBins){
      cutVertex->SetCheckZResolutionSPD();
      cutVertex->SetCheckDispersionSPD();
      cutVertex->SetCheckZDifferenceSPDTrack();
    }
    if(0) cutVertex->SetCheckGeneratedVertexZ();
  }

  // other event selection cuts
  AliRsnCutEventUtils* cutEventUtils=0;
  if(1){
    cutEventUtils=new AliRsnCutEventUtils("cutEventUtils",kTRUE,rejectPileUp);
    if(!MultBins){
      cutEventUtils->SetCheckIncompleteDAQ();
      cutEventUtils->SetCheckSPDClusterVsTrackletBG();
    }else{
      cutEventUtils->SetRemovePileUppA2013(kFALSE);
      cutEventUtils->SetCheckAcceptedMultSelection();
    }
  }

  // set the check for pileup 
  if(isPP && (!isMC) && cutVertex){ 
    cutVertex->SetCheckPileUp(rejectPileUp); 
    ::Info("AddTaskRare_pp13", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp)?"ON":"OFF"));
  }

  // define and fill cut set for event cuts
  AliRsnCutSet* eventCuts=0;
  if(cutEventUtils || cutVertex){
    eventCuts=new AliRsnCutSet("eventCuts",AliRsnTarget::kEvent);

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
  }

  // ----- EVENT-ONLY COMPUTATIONS -----
    
  Double_t multbins[1000];
  int j,nmult=0;
  if(!MultBins){
    for(j=0;j<=401;j++){multbins[nmult]=j-0.5; nmult++;}
  }else if(!trigger){
    for(j=0;j<=100;j++){multbins[nmult]=j; nmult++;}
  }else{
    for(j=0;j<10;j++){multbins[nmult]=0.0001*j; nmult++;}
    for(j=1;j<10;j++){multbins[nmult]=0.001*j; nmult++;}
    for(j=1;j<10;j++){multbins[nmult]=0.01*j; nmult++;}
    for(j=1;j<10;j++){multbins[nmult]=0.1*j; nmult++;}
    for(j=1;j<=100;j++){multbins[nmult]=j; nmult++;}
  }
  nmult--;

  //vertex
  Int_t vtxID=task->CreateValue(AliRsnMiniValue::kVz,kFALSE);
  AliRsnMiniOutput* outVtx=task->CreateOutput("eventVtx","HIST","EVENT");
  outVtx->AddAxis(vtxID,240,-12.0,12.0);

  //multiplicity or centrality
  Int_t multID=task->CreateValue(AliRsnMiniValue::kMult,kFALSE);
  AliRsnMiniOutput* outMult=task->CreateOutput("eventMult","HIST","EVENT");
  outMult->AddAxis(multID,nmult+1,multbins);

  TH1F* hEventsVsMulti=new TH1F("hAEventsVsMulti","",nmult,multbins);
  task->SetEventQAHist("EventsVsMulti",hEventsVsMulti);//custom binning for fHAEventsVsMulti
    
  double ybins[1000];
  for(j=0;j<=240;j++) ybins[j]=-12+0.1*j;

  TH2F* hvz=new TH2F("hVzVsCent","",nmult,multbins, 240,ybins);
  task->SetEventQAHist("vz",hvz);//plugs this histogram into the fHAEventVz data member

  for(j=0;j<=401;j++) ybins[j]=j-0.5;

  TH2F* hmc=new TH2F("MultiVsCent","", nmult,multbins, 401,ybins);
  hmc->GetYaxis()->SetTitle("QUALITY");
  task->SetEventQAHist("multicent",hmc);//plugs this histogram into the fHAEventMultiCent data member

  // ----- CONFIGURE -----

  cerr<<"configuring"<<endl;
  if(d1==AliRsnDaughter::kPion && d2==AliRsnDaughter::kPion){
    Config_pipi(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
        
  }else if(d1==AliRsnDaughter::kPion && d2==AliRsnDaughter::kKaon){
    Config_pikx(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kPion && d1==AliRsnDaughter::kKaon){
    Config_pikx(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kPion && d2==AliRsnDaughter::kKaon0){
    Config_pik0(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kPion && d1==AliRsnDaughter::kKaon0){
    Config_pik0(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kKaon && d2==AliRsnDaughter::kKaon){
    Config_kxkx(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);

  }else if(d1==AliRsnDaughter::kKaon && d2==AliRsnDaughter::kKaon0){
    Config_kxk0(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kKaon && d1==AliRsnDaughter::kKaon0){
    Config_kxk0(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kProton && d2==AliRsnDaughter::kKaon){
    Config_pkx(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kProton && d1==AliRsnDaughter::kKaon){
    Config_pkx(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kProton && d2==AliRsnDaughter::kKaon0){
    Config_pk0(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kProton && d1==AliRsnDaughter::kKaon0){
    Config_pk0(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kLambda && d2==AliRsnDaughter::kPion){
    Config_Lambdapi(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kLambda && d1==AliRsnDaughter::kPion){
    Config_Lambdapi(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kLambda && d2==AliRsnDaughter::kKaon){
    Config_Lambdakx(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kLambda && d1==AliRsnDaughter::kKaon){
    Config_Lambdakx(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kLambda && d2==AliRsnDaughter::kKaon0){
    Config_Lambdak0(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kLambda && d1==AliRsnDaughter::kKaon0){
    Config_Lambdak0(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kLambda && d2==AliRsnDaughter::kProton){
    Config_Lambdap(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kLambda && d1==AliRsnDaughter::kProton){
    Config_Lambdap(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
  }
  cerr<<"done configuring"<<endl;
  
  // ----- CONTAINERS -----

  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  Printf("AddTaskRare_pp13 - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
  AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",lname.Data()), 
							  TList::Class(), 
							  AliAnalysisManager::kOutputContainer, 
							  outputFileName);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);
   
  return task;
}


//=============================


Bool_t Config_pipi(
  AliRsnMiniAnalysisTask *task,
  TString     lname="pipi",
  Bool_t      isMC=kFALSE,
  Int_t       system=0,
  Int_t       EventCuts=0,
  Int_t       TrackCutsPi=0,
  Int_t       TrackCuts2=0
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;
  if(system==1 || system==2) MultBins=1;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;
  
  // retrieve mass from PDG database
  Int_t pdg=TrackCuts2;
  TDatabasePDG* db=TDatabasePDG::Instance();
  TParticlePDG* part=db->GetParticle(pdg);
  Double_t mass=part->Mass();
  
  // set daughter cuts
  if(!(TrackCutsPi%10000)) TrackCutsPi+=3020;//default settings
  Float_t nsigmaPiTPC=0.1*(TrackCutsPi%100);
  Float_t nsigmaPiTOF=0.1*((TrackCutsPi/100)%100);
  Int_t CutTypePi=(TrackCutsPi/10000)%100;//0=TPC+TOFveto (default), 1=TPC only, 2=TOF only

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);

  AliRsnCutSetDaughterParticle* cutSetPi=0;
  if(!CutTypePi) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaPiTPC),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kPion,nsigmaPiTPC,nsigmaPiTOF);
  else if(CutTypePi==1) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,nsigmaPiTPC),trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kPion,nsigmaPiTPC,-1.);
  else if(CutTypePi==2) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,nsigmaPiTOF),trkQualityCut,AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,AliPID::kPion,-1.,nsigmaPiTOF);
  if(!cutSetPi){cerr<<"Error in AddTaskRare_pp13::Config_pipi(): missing cutSetPi"<<endl; return kFALSE;}

  Int_t iCutQ=task->AddTrackCuts(cutSetQ);
  Int_t iCutPi=task->AddTrackCuts(cutSetPi);

  // monitoring
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
    AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());
  }

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  if(system!=1) cutY->SetRangeD(-0.5,0.5);
  else cutY->SetRangeD(-0.465,0.035);
  AliRsnCutSet* cutsPair=new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
  cutsPair->AddCut(cutY);
  cutsPair->SetCutScheme(cutY->GetName());
    
  // multiplicity binning
  Double_t multbins[200];
  int j,nmult=0;
  if(!MultBins){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.e6; nmult++;
  }else if(!trigger){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.; nmult++;
    multbins[nmult]=5.; nmult++;
    multbins[nmult]=10.; nmult++;
    multbins[nmult]=15.; nmult++;
    for(j=2;j<=10;j++){multbins[nmult]=j*10; nmult++;}
  }else{
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=0.001; nmult++;
    multbins[nmult]=0.005; nmult++;
    multbins[nmult]=0.01; nmult++;
    multbins[nmult]=0.05; nmult++;
    multbins[nmult]=0.1; nmult++;
    multbins[nmult]=1.; nmult++;
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
    
    
  Bool_t  use    [9]={ 1      ,  1     , 1      ,  1     , isMC   , isMC , isMC ,  0       ,  0       };
  Int_t   useIM  [9]={ 1      ,  1     , 1      ,  1     ,  1     ,  1   ,  0   ,  1       ,  1       };
  TString name   [9]={"Unlike","Mixing","LikePP","LikeMM","gen"   ,"true","res" ,"MixingPP","MixingMM"};
  TString comp   [9]={"PAIR"  , "MIX"  ,"PAIR"  ,"PAIR"  ,"MOTHER","TRUE","TRUE","MIX"     ,"MIX"     };
  TString output [9]={"HIST"  ,"HIST"  ,"HIST"  ,"HIST"  ,"HIST"  ,"HIST","HIST","HIST"    ,"HIST"    };
  Char_t  charge1[9]={'+'     , '+'    ,'+'     ,'-'     , '+'    , '+'  ,'+'   ,'+'       ,'-'       };
  Char_t  charge2[9]={'-'     , '-'    ,'+'     ,'-'     , '-'    , '-'  ,'-'   ,'+'       ,'-'       };

  for(Int_t i=0;i<9;i++){
    if(!use[i]) continue;
    AliRsnMiniOutput *out=task->CreateOutput(Form("pipi_%s%s",name[i].Data(),suffix),output[i].Data(),comp[i].Data());
    out->SetDaughter(0,AliRsnDaughter::kPion);
    out->SetDaughter(1,AliRsnDaughter::kPion);
    out->SetCutID(0,iCutPi);
    out->SetCutID(1,iCutPi);
    out->SetCharge(0,charge1[i]);
    out->SetCharge(1,charge2[i]);
    out->SetMotherPDG(pdg);
    out->SetMotherMass(mass);
    ////out->SetPairCuts(cutsPair);

    // axis X: invmass or resolution
    if(useIM[i]) out->AddAxis(imID,173,0.27,2.);
    else out->AddAxis(resID,200,-0.02,0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID,200,0.0,20.0);
    
    // axis Z: centrality-multiplicity
    out->AddAxis(centID,nmult,multbins);
  }
  return kTRUE;
}


//=============================


Bool_t Config_pikx(
  AliRsnMiniAnalysisTask *task,
  TString     lname="pikx",
  Bool_t      isMC=kFALSE,
  Int_t       system=0,
  Int_t       EventCuts=0,
  Int_t       TrackCutsPi=0,
  Int_t       TrackCutsK=0
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;
  if(system==1 || system==2) MultBins=1;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;
  
  // retrieve mass from PDG database
  Int_t pdg=313;
  TDatabasePDG* db=TDatabasePDG::Instance();
  TParticlePDG* part=db->GetParticle(pdg);
  Double_t mass=part->Mass();
  
  // set daughter cuts
  if(!(TrackCutsPi%10000)) TrackCutsPi+=3020;//default settings
  Float_t nsigmaPiTPC=0.1*(TrackCutsPi%100);
  Float_t nsigmaPiTOF=0.1*((TrackCutsPi/100)%100);
  Int_t CutTypePi=(TrackCutsPi/10000)%100;//0=TPC+TOF (default), 1=TPC only, 2=TOF only

  if(!(TrackCutsK%10000)) TrackCutsK+=3020;//default settings
  Float_t nsigmaKTPC=0.1*(TrackCutsK%100);
  Float_t nsigmaKTOF=0.1*((TrackCutsK/100)%100);
  Int_t CutTypeK=(TrackCutsK/10000)%100;//0=TPC+TOFveto (default), 1=TPC only, 2=TOF only

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);

  AliRsnCutSetDaughterParticle* cutSetPi=0;
  if(!CutTypePi) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaPiTPC),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kPion,nsigmaPiTPC,nsigmaPiTOF);
  else if(CutTypePi==1) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,nsigmaPiTPC),trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kPion,nsigmaPiTPC,-1.);
  else if(CutTypePi==2) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,nsigmaPiTOF),trkQualityCut,AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,AliPID::kPion,-1.,nsigmaPiTOF);
  if(!cutSetPi){cerr<<"Error in AddTaskRare_pp13::Config_pikx(): missing cutSetPi"<<endl; return kFALSE;}

  AliRsnCutSetDaughterParticle* cutSetK=0;
  if(!CutTypeK) cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaKTPC),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kKaon,nsigmaKTPC,nsigmaKTOF);
  else if(CutTypeK==1) cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,nsigmaKTPC),trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kKaon,nsigmaKTPC,-1.);
  else if(CutTypeK==2) cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,nsigmaKTOF),trkQualityCut,AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,AliPID::kKaon,-1.,nsigmaKTOF);
  if(!cutSetK){cerr<<"Error in AddTaskRare_pp13::Config_pikx(): missing cutSetK"<<endl; return kFALSE;}

  Int_t iCutQ=task->AddTrackCuts(cutSetQ);
  Int_t iCutPi=task->AddTrackCuts(cutSetPi);
  Int_t iCutK=task->AddTrackCuts(cutSetK);

  // monitoring
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
    AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());
    AddMonitorOutput(isMC,cutSetK->GetMonitorOutput());
  }

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  if(system!=1) cutY->SetRangeD(-0.5,0.5);
  else cutY->SetRangeD(-0.465,0.035);
  AliRsnCutSet* cutsPair=new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
  cutsPair->AddCut(cutY);
  cutsPair->SetCutScheme(cutY->GetName());
    
  // multiplicity binning
  Double_t multbins[200];
  int j,nmult=0;
  if(!MultBins){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.e6; nmult++;
  }else if(!trigger){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.; nmult++;
    multbins[nmult]=5.; nmult++;
    multbins[nmult]=10.; nmult++;
    multbins[nmult]=15.; nmult++;
    for(j=2;j<=10;j++){multbins[nmult]=j*10; nmult++;}
  }else{
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=0.001; nmult++;
    multbins[nmult]=0.005; nmult++;
    multbins[nmult]=0.01; nmult++;
    multbins[nmult]=0.05; nmult++;
    multbins[nmult]=0.1; nmult++;
    multbins[nmult]=1.; nmult++;
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

  Bool_t  use     [12] = {1         ,1         ,1         ,1         ,1       ,1       ,isMC     ,isMC     ,isMC     ,isMC     ,isMC    ,isMC    };
  Bool_t  useIM   [12] = {1         ,1         ,1         ,1         ,1       ,1       ,1        ,1        ,1        ,1        ,0       ,0       };
  TString name    [12] = {"UnlikePM","UnlikeMP","MixingPM","MixingMP","LikePP","LikeMM","MCGenPM","MCGenMP","TruesPM","TruesMP","ResPM" ,"ResMP" };
  TString comp    [12] = {"PAIR"    ,"PAIR"    ,"MIX"     ,"MIX"     ,"PAIR"  ,"PAIR"  ,"MOTHER" ,"MOTHER" ,"TRUE"   ,"TRUE"   ,"TRUE"  ,"TRUE"  };
  Char_t  charge1 [12] = {'+'       ,'-'       ,'+'       ,'-'       ,'+'     ,'-'     ,'+'      ,'-'      ,'+'      ,'-'      ,'+'     ,'-'     };//K
  Char_t  charge2 [12] = {'-'       ,'+'       ,'-'       ,'+'       ,'+'     ,'-'     ,'-'      ,'+'      ,'_'      ,'+'      ,'-'     ,'+'     };//pi
  Int_t   PDGCode [12] = {313       ,313       ,313       ,313       ,313     ,313     ,313      ,-313     ,313      ,-313     ,313     ,-313    };

  for(Int_t i=0;i<12;i++){
    if(!use[i]) continue;
    AliRsnMiniOutput *out=task->CreateOutput(Form("pikx_%s%s",name[i].Data(),suffix),"HIST",comp[i].Data());
    out->SetDaughter(0,AliRsnDaughter::kKaon);
    out->SetDaughter(1,AliRsnDaughter::kPion);
    out->SetCutID(0,iCutK);
    out->SetCutID(1,iCutPi);
    out->SetCharge(0,charge1[i]);
    out->SetCharge(1,charge2[i]);
    out->SetMotherPDG(PDGCode[i]);
    out->SetMotherMass(mass);
    out->SetPairCuts(cutsPair);

    // axis X: invmass or resolution
    if(useIM[i]) out->AddAxis(imID,137,0.63,2.);
    else out->AddAxis(resID,200,-0.02,0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID,200,0.0,20.0);
    
    // axis Z: centrality-multiplicity
    out->AddAxis(centID,nmult,multbins);
  }
  return kTRUE;
}


//=============================


Bool_t Config_pik0(
  AliRsnMiniAnalysisTask *task,
  TString     lname="pik0",
  Bool_t      isMC=kFALSE,
  Int_t       system=0,
  Int_t       EventCuts=0,
  Int_t       TrackCutsPi=0,
  Int_t       TrackCutsK=0
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;
  if(system==1 || system==2) MultBins=1;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;
  
  // set cuts for primary pion
  if(!(TrackCutsPi%10000)) TrackCutsPi+=3020;//default settings
  Float_t nsigmaPiTPC=0.1*(TrackCutsPi%100);
  Float_t nsigmaPiTOF=0.1*((TrackCutsPi/100)%100);

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);
  AliRsnCutSetDaughterParticle* cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaPiTPC),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kPion,nsigmaPiTPC,nsigmaPiTOF);
  
  Int_t iCutQ=task->AddTrackCuts(cutSetQ);
  Int_t iCutPi=task->AddTrackCuts(cutSetPi);
    
  // selections for V0 daughters
  Int_t v0d_xrows=70;
  Float_t v0d_rtpc=0.8;
  Float_t v0d_dcaxy=0.06;
    
  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts("qualityDaughterK0s");
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetRequireTPCRefit();
  esdTrackCuts->SetAcceptKinkDaughters(0);
  esdTrackCuts->SetMinNCrossedRowsTPC(v0d_xrows);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(v0d_rtpc);
  esdTrackCuts->SetMinDCAToVertexXY(v0d_dcaxy);
    
  // selections for K0S
  Float_t k0s_piPIDCut=5.;
  Float_t k0sDaughDCA=1.;
  Float_t k0sDCA=1.e10;//0.3;
  Float_t k0s_pLife=20.;
  Float_t k0s_radiuslow=0.5;
  Float_t k0s_radiushigh=200.;
  Float_t k0s_massTol=0.03;
  Float_t k0s_massTolVeto=0.004;
  Bool_t  k0sSwitch=kFALSE;
  Float_t k0sCosPoinAn=0.97;
    
  AliRsnCutV0* cutK0s=new AliRsnCutV0("cutK0s",kK0Short,AliPID::kPion,AliPID::kPion);
  cutK0s->SetPIDCutPion(k0s_piPIDCut);// PID for the pion daughters of K0S
  cutK0s->SetESDtrackCuts(esdTrackCuts);// all the other selections (defined above) for pion daughters of K0S
  cutK0s->SetMaxDaughtersDCA(k0sDaughDCA);
  cutK0s->SetMaxDCAVertex(k0sDCA);
  cutK0s->SetfLife(k0s_pLife);
  cutK0s->SetfLowRadius(k0s_radiuslow);
  cutK0s->SetfHighRadius(k0s_radiushigh);
  cutK0s->SetTolerance(k0s_massTol);
  cutK0s->SetToleranceVeto(k0s_massTolVeto);//Rejection range for Competing V0 Rejection
  cutK0s->SetSwitch(k0sSwitch);
  cutK0s->SetMinCosPointingAngle(k0sCosPoinAn);
  cutK0s->SetMaxRapidity(2.);
  cutK0s->SetMinTPCcluster(-1);
    
  AliRsnCutSet* cutSetK0s=new AliRsnCutSet("setK0s",AliRsnTarget::kDaughter);
  cutSetK0s->AddCut(cutK0s);
  cutSetK0s->SetCutScheme(cutK0s->GetName());
  Int_t iCutK0s=task->AddTrackCuts(cutSetK0s);

  // monitoring
  TString pname="k0";
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
    AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());

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
  }

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  if(system!=1) cutY->SetRangeD(-0.5,0.5);
  else cutY->SetRangeD(-0.465,0.035);

  AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
   
  AliRsnCutSet* cutsPairSame=new AliRsnCutSet("pairCutsSame",AliRsnTarget::kMother);
  cutsPairSame->AddCut(cutY);
  cutsPairSame->AddCut(cutV0);
  cutsPairSame->SetCutScheme(TString::Format("%s&(!%s)",cutY->GetName(),cutV0->GetName()).Data());
  
  AliRsnCutSet* cutsPairMix=new AliRsnCutSet("pairCutsMix", AliRsnTarget::kMother);
  cutsPairMix->AddCut(cutY);
  cutsPairMix->SetCutScheme(cutY->GetName());
    
  // multiplicity binning
  Double_t multbins[200];
  int j,nmult=0;
  if(!MultBins){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.e6; nmult++;
  }else if(!trigger){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.; nmult++;
    multbins[nmult]=5.; nmult++;
    multbins[nmult]=10.; nmult++;
    multbins[nmult]=15.; nmult++;
    for(j=2;j<=10;j++){multbins[nmult]=j*10; nmult++;}
  }else{
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=0.001; nmult++;
    multbins[nmult]=0.005; nmult++;
    multbins[nmult]=0.01; nmult++;
    multbins[nmult]=0.05; nmult++;
    multbins[nmult]=0.1; nmult++;
    multbins[nmult]=1.; nmult++;
  }

  // -- Values ------------------------------------------------------------------------------------                                    
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);

  // -- Create all needed outputs ----------------------------------------------------------------- 
  // use an array for more compact writing, which are different on mixing and charges

  Bool_t  use     [10] = {1      ,1      ,1         ,1         ,isMC        ,isMC        ,isMC         ,isMC         ,isMC        ,isMC        };
  Bool_t  useIM   [10] = {1      ,1      ,1         ,1         ,1           ,1           ,1            ,1            ,0           ,0           };
  TString name    [10] = {"K0Pip","K0Pim","K0PipMix","K0PimMix","Kstarp_gen","Kstarm_gen","Kstarp_true","Kstarm_true","Kstarp_res","Kstarm_res"};
  TString comp    [10] = {"PAIR" ,"PAIR" ,"MIX"     ,"MIX"     ,"MOTHER"    ,"MOTHER"    ,"TRUE"       ,"TRUE"       ,"TRUE"      ,"TRUE"      };
  Char_t  charge1 ='0';
  Char_t  charge2 [10] = {'+'    ,'-'    ,'+'       ,'-'       ,'+'         ,'-'         ,'+'          ,'-'          ,'+'         ,'-'         };
  Int_t   ipdg    [10] = {323    ,-323   ,323       ,-323      ,323         ,-323        ,323          ,-323         ,323         ,-323        };
  Double_t mass = 0.89166;
  Int_t   pairID  [10] = { 0     ,  0    ,  1       ,  1       ,  1         ,  1         ,  1          ,  1          ,  1         ,  1         };
   
  for(Int_t i=0;i<10;i++){
    if (!use[i]) continue;
    // create output
    AliRsnMiniOutput* out=task->CreateOutput(Form("pik0_%s%s",name[i].Data(),suffix),"HIST",comp[i].Data());
    // selection settings
    out->SetCutID(0,iCutK0s);
    out->SetCutID(1,iCutPi);
    out->SetDaughter(0,AliRsnDaughter::kKaon0);
    out->SetDaughter(1,AliRsnDaughter::kPion);
    out->SetCharge(0,charge1);
    out->SetCharge(1,charge2[i]);
    out->SetMotherPDG(ipdg[i]);
    out->SetMotherMass(mass);
    // pair cuts
    if(TrackCutsK & 1024){
      out->SetPairCuts(cutsPairMix);
    }else if(TrackCutsK & 2048){
      out->SetPairCuts(cutsPairSame);
    }else{
      if(!pairID[i]) out->SetPairCuts(cutsPairSame);
      else out->SetPairCuts(cutsPairMix);
    }

    // axis X: invmass or resolution
    if(useIM[i]) out->AddAxis(imID,1370,0.63,2.);
    else out->AddAxis(resID,200,-0.02,0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID,200,0.,20.);
    
    // axis Z: centrality-multiplicity
    out->AddAxis(centID,nmult,multbins);
   }

  return kTRUE;
}


//=============================


Bool_t Config_kxkx(
  AliRsnMiniAnalysisTask *task,
  TString     lname="kxkx",
  Bool_t      isMC=kFALSE,
  Int_t       system=0,
  Int_t       EventCuts=0,
  Int_t       TrackCutsK=0,
  Int_t       TrackCutsDummy=0
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;
  if(system==1 || system==2) MultBins=1;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;
  
  // retrieve mass from PDG database
  Int_t pdg=333;
  TDatabasePDG* db=TDatabasePDG::Instance();
  TParticlePDG* part=db->GetParticle(pdg);
  Double_t mass=part->Mass();
  
  // set daughter cuts
  if(!(TrackCutsK%10000)) TrackCutsK+=3020;//default settings
  Float_t nsigmaKTPC=0.1*(TrackCutsK%100);
  Float_t nsigmaKTOF=0.1*((TrackCutsK/100)%100);

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);
  AliRsnCutSetDaughterParticle* cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaKTPC),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kKaon,nsigmaKTPC,nsigmaKTOF);

  Int_t iCutQ=task->AddTrackCuts(cutSetQ);
  Int_t iCutK=task->AddTrackCuts(cutSetK);

  // monitoring
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
    AddMonitorOutput(isMC,cutSetK->GetMonitorOutput());
  }

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  if(system!=1) cutY->SetRangeD(-0.5,0.5);
  else cutY->SetRangeD(-0.465,0.035);
  AliRsnCutSet* cutsPair=new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
  cutsPair->AddCut(cutY);
  cutsPair->SetCutScheme(cutY->GetName());
    
  // multiplicity binning
  Double_t multbins[200];
  int j,nmult=0;
  if(!MultBins){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.e6; nmult++;
  }else if(!trigger){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.; nmult++;
    multbins[nmult]=5.; nmult++;
    multbins[nmult]=10.; nmult++;
    multbins[nmult]=15.; nmult++;
    for(j=2;j<=10;j++){multbins[nmult]=j*10; nmult++;}
  }else{
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=0.001; nmult++;
    multbins[nmult]=0.005; nmult++;
    multbins[nmult]=0.01; nmult++;
    multbins[nmult]=0.05; nmult++;
    multbins[nmult]=0.1; nmult++;
    multbins[nmult]=1.; nmult++;
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
    
  Bool_t  use    [9]={ 1      ,  1     , 1      ,  1     , isMC   , isMC , isMC ,  0       ,  0       };
  Int_t   useIM  [9]={ 1      ,  1     , 1      ,  1     ,  1     ,  1   ,  0   ,  1       ,  1       };
  TString name   [9]={"Unlike","Mixing","LikePP","LikeMM","gen"   ,"true","res" ,"MixingPP","MixingMM"};
  TString comp   [9]={"PAIR"  , "MIX"  ,"PAIR"  ,"PAIR"  ,"MOTHER","TRUE","TRUE","MIX"     ,"MIX"     };
  Char_t  charge1[9]={'+'     , '+'    ,'+'     ,'-'     , '+'    , '+'  ,'+'   ,'+'       ,'-'       };
  Char_t  charge2[9]={'-'     , '-'    ,'+'     ,'-'     , '-'    , '-'  ,'-'   ,'+'       ,'-'       };

  for(Int_t i=0;i<9;i++){
    if(!use[i]) continue;
    AliRsnMiniOutput* out=task->CreateOutput(Form("kxkx_%s%s",name[i].Data(),suffix),"HIST",comp[i].Data());
    out->SetCutID(0,iCutK);
    out->SetCutID(1,iCutK);
    out->SetDaughter(0,AliRsnDaughter::kKaon);
    out->SetDaughter(1,AliRsnDaughter::kKaon);
    out->SetCharge(0,charge1[i]);
    out->SetCharge(1,charge2[i]);
    out->SetMotherPDG(pdg);
    out->SetMotherMass(mass);
    out->SetPairCuts(cutsPair);

    // axis X: invmass or resolution
    if(useIM[i]) out->AddAxis(imID,1015,0.985,2.);
    else out->AddAxis(resID,200,-0.02,0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID,200,0.0,20.0);
    
    // axis Z: centrality-multiplicity
    out->AddAxis(centID,nmult,multbins);
  }
  return kTRUE;
}


//=============================


Bool_t Config_kxk0(
  AliRsnMiniAnalysisTask *task,
  TString     lname="kxk0",
  Bool_t      isMC=kFALSE,
  Int_t       system=0,
  Int_t       EventCuts=0,
  Int_t       TrackCutsKx=0,
  Int_t       TrackCutsK0=0
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;
  if(system==1 || system==2) MultBins=1;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;
  
  // set cuts for primary K+/-
  if(!(TrackCutsKx%10000)) TrackCutsKx+=3020;//default settings
  Float_t nsigmaKxTPC=0.1*(TrackCutsKx%100);
  Float_t nsigmaKxTOF=0.1*((TrackCutsKx/100)%100);

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);
  AliRsnCutSetDaughterParticle* cutSetKx=new AliRsnCutSetDaughterParticle(Form("cutKx%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaKxTPC),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kKaon,nsigmaKxTPC,nsigmaKxTOF);
  
  Int_t iCutQ=task->AddTrackCuts(cutSetQ);
  Int_t iCutKx=task->AddTrackCuts(cutSetKx);
    
  // selections for V0 daughters
  Int_t v0d_xrows=70;
  Float_t v0d_rtpc=0.8;
  Float_t v0d_dcaxy=0.06;
    
  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts("qualityDaughterK0s");
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetRequireTPCRefit();
  esdTrackCuts->SetAcceptKinkDaughters(0);
  esdTrackCuts->SetMinNCrossedRowsTPC(v0d_xrows);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(v0d_rtpc);
  esdTrackCuts->SetMinDCAToVertexXY(v0d_dcaxy);
    
  // selections for K0S
  Float_t k0s_piPIDCut=5.;
  Float_t k0sDaughDCA=1.;
  Float_t k0sDCA=1.e10;//0.3;
  Float_t k0s_pLife=20.;
  Float_t k0s_radiuslow=0.5;
  Float_t k0s_radiushigh=200.;
  Float_t k0s_massTol=0.03;
  Float_t k0s_massTolVeto=0.004;
  Bool_t  k0sSwitch=kFALSE;
  Float_t k0sCosPoinAn=0.97;
    
  AliRsnCutV0* cutK0s=new AliRsnCutV0("cutK0s",kK0Short,AliPID::kPion,AliPID::kPion);
  cutK0s->SetPIDCutPion(k0s_piPIDCut);// PID for the pion daughters of K0S
  cutK0s->SetESDtrackCuts(esdTrackCuts);// all the other selections (defined above) for pion daughters of K0S
  cutK0s->SetMaxDaughtersDCA(k0sDaughDCA);
  cutK0s->SetMaxDCAVertex(k0sDCA);
  cutK0s->SetfLife(k0s_pLife);
  cutK0s->SetfLowRadius(k0s_radiuslow);
  cutK0s->SetfHighRadius(k0s_radiushigh);
  cutK0s->SetTolerance(k0s_massTol);
  cutK0s->SetToleranceVeto(k0s_massTolVeto);//Rejection range for Competing V0 Rejection
  cutK0s->SetSwitch(k0sSwitch);
  cutK0s->SetMinCosPointingAngle(k0sCosPoinAn);
  cutK0s->SetMaxRapidity(2.);
  cutK0s->SetMinTPCcluster(-1);
    
  AliRsnCutSet* cutSetK0s=new AliRsnCutSet("setK0s",AliRsnTarget::kDaughter);
  cutSetK0s->AddCut(cutK0s);
  cutSetK0s->SetCutScheme(cutK0s->GetName());
  Int_t iCutK0s=task->AddTrackCuts(cutSetK0s);

  // monitoring
  TString pname="k0";
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
    AddMonitorOutput(isMC,cutSetKx->GetMonitorOutput());

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
  }

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  if(system!=1) cutY->SetRangeD(-0.5,0.5);
  else cutY->SetRangeD(-0.465,0.035);

  AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
   
  AliRsnCutSet* cutsPairSame=new AliRsnCutSet("pairCutsSame",AliRsnTarget::kMother);
  cutsPairSame->AddCut(cutY);
  cutsPairSame->AddCut(cutV0);
  cutsPairSame->SetCutScheme(TString::Format("%s&(!%s)",cutY->GetName(),cutV0->GetName()).Data());
  
  AliRsnCutSet* cutsPairMix=new AliRsnCutSet("pairCutsMix", AliRsnTarget::kMother);
  cutsPairMix->AddCut(cutY);
  cutsPairMix->SetCutScheme(cutY->GetName());
    
  // multiplicity binning
  Double_t multbins[200];
  int j,nmult=0;
  if(!MultBins){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.e6; nmult++;
  }else if(!trigger){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.; nmult++;
    multbins[nmult]=5.; nmult++;
    multbins[nmult]=10.; nmult++;
    multbins[nmult]=15.; nmult++;
    for(j=2;j<=10;j++){multbins[nmult]=j*10; nmult++;}
  }else{
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=0.001; nmult++;
    multbins[nmult]=0.005; nmult++;
    multbins[nmult]=0.01; nmult++;
    multbins[nmult]=0.05; nmult++;
    multbins[nmult]=0.1; nmult++;
    multbins[nmult]=1.; nmult++;
  }

  // -- Values ------------------------------------------------------------------------------------                                    
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* transv. momentum */ Int_t ptID    = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality       */ Int_t centID  = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity   */ Int_t etaID   = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity         */ Int_t yID     = task->CreateValue(AliRsnMiniValue::kY, kFALSE);

  // -- Create all needed outputs ----------------------------------------------------------------- 
  // use an array for more compact writing, which are different on mixing and charges

  Bool_t  use     [6] = {1               ,1                ,1                  ,1                   ,0                ,0                 };
  Bool_t  useIM   [6] = {1               ,1                ,1                  ,1                   ,1                ,1                 };
  TString name    [6] = {"K0Kp"         ,"K0Km"          ,"K0KpMix"         ,"K0KmMix"          ,"KStarPlusMinust","AKStarPlusMinust"};
  TString comp    [6] = {"PAIR"          ,"PAIR"           ,"MIX"              ,"MIX"               ,"TRUE"           ,"TRUE"            };
  TString output  [6] = {"HIST"        ,"HIST"         ,"HIST"           ,"HIST"            ,"HIST"         ,"HIST"          };
  Char_t  charge1 [6] = {'0'             ,'0'              ,'0'                ,'0'                 ,'0'              ,'0'               };
  Char_t  charge2 [6] = {'+'             ,'-'              ,'+'                ,'-'                 ,'+'              ,'-'               };
  Int_t   cutID1  [6] = { iCutK0s        ,iCutK0s           ,iCutK0s            ,iCutK0s            ,iCutK0s          ,iCutK0s           };
  Int_t   cutID2  [6] = { iCutKx         ,iCutKx           ,iCutKx             ,iCutKx              ,iCutKx           ,iCutKx            };
  Int_t   ipdg    [6] = {323             ,-323             ,323                ,-323                ,323              ,-323              };
  Double_t mass   [6] = { 0.89166        ,0.89166          ,0.89166            ,0.89166             ,0.89166          ,0.89166           };
   
  for(Int_t i=0;i<6;i++){
    if (!use[i]) continue;
    // create output
    AliRsnMiniOutput* out=task->CreateOutput(Form("kxk0_%s%s",name[i].Data(),suffix),output[i].Data(),comp[i].Data());
    // selection settings
    out->SetCutID(0,cutID1[i]);
    out->SetCutID(1,cutID2[i]);
    out->SetDaughter(0,AliRsnDaughter::kKaon0);
    out->SetDaughter(1,AliRsnDaughter::kKaon);
    out->SetCharge(0,charge1[i]);
    out->SetCharge(1,charge2[i]);
    out->SetMotherPDG(ipdg[i]);
    out->SetMotherMass(mass[i]);
    // pair cuts
    if(TrackCutsK0 & 1024){
      out->SetPairCuts(cutsPairMix);
    }else if(TrackCutsK0 & 2048){
      out->SetPairCuts(cutsPairSame);
    }else{
      if(i==0 || i==1 || i==4 || i==5) out->SetPairCuts(cutsPairSame);
      else out->SetPairCuts(cutsPairMix);
    }

    // axis X: invmass or resolution
    
    if(useIM[i]) out->AddAxis(imID, 1005, 0.99, 3.);
    else out->AddAxis(resID,200,-0.02,0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID,200,0.,20.);
    
    // axis Z: centrality-multiplicity
    out->AddAxis(centID,nmult,multbins);
   }

  return kTRUE;
}


//=============================


Bool_t Config_pkx(
  AliRsnMiniAnalysisTask *task,
  TString     lname="pkx",
  Bool_t      isMC=kFALSE,
  Int_t       system=0,
  Int_t       EventCuts=0,
  Int_t       TrackCutsP=0,
  Int_t       TrackCutsK=0
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;
  if(system==1 || system==2) MultBins=1;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;

  Double_t mass=1.51953;
  
  // set daughter cuts
  if(!(TrackCutsP%10000)) TrackCutsP+=3020;//default settings
  Float_t nsigmaPTPC=0.1*(TrackCutsP%100);
  Float_t nsigmaPTOF=0.1*((TrackCutsP/100)%100);
  Int_t CutTypeP=(TrackCutsP/10000)%100;//0=TPC+TOF (default), 1=TPC only, 2=TOF only

  if(!(TrackCutsK%10000)) TrackCutsK+=3020;//default settings
  Float_t nsigmaKTPC=0.1*(TrackCutsK%100);
  Float_t nsigmaKTOF=0.1*((TrackCutsK/100)%100);
  Int_t CutTypeK=(TrackCutsK/10000)%100;//0=TPC+TOF (default), 1=TPC only, 2=TOF only

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);

  AliRsnCutSetDaughterParticle* cutSetP=0;
  if(!CutTypeP) cutSetP=new AliRsnCutSetDaughterParticle(Form("cutProton_%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV,nsigmaPTPC),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV,AliPID::kProton,nsigmaPTPC);
  else if(CutTypeP==1) cutSetP=new AliRsnCutSetDaughterParticle(Form("cutProton%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,nsigmaPTPC),trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kProton,nsigmaPTPC,-1.);
  else if(CutTypeP==2) cutSetP=new AliRsnCutSetDaughterParticle(Form("cutProton%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,nsigmaPTOF),trkQualityCut,AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,AliPID::kProton,-1.,nsigmaPTOF);
  if(!cutSetP){cerr<<"Error in AddTaskRare_pp13::Config_pkx(): missing cutSetP"<<endl; return kFALSE;}

  AliRsnCutSetDaughterParticle* cutSetK=0;
  if(!CutTypeK) cutSetK=new AliRsnCutSetDaughterParticle(Form("cutKaon_%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV,nsigmaKTPC),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV,AliPID::kKaon,nsigmaKTPC);
  else if(CutTypeK==1) cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,nsigmaKTPC),trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kKaon,nsigmaKTPC,-1.);
  else if(CutTypeK==2) cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,nsigmaKTOF),trkQualityCut,AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,AliPID::kKaon,-1.,nsigmaKTOF);
  if(!cutSetK){cerr<<"Error in AddTaskRare_pp13::Config_pkx(): missing cutSetK"<<endl; return kFALSE;}

  Int_t iCutQ=task->AddTrackCuts(cutSetQ);
  Int_t iCutP=task->AddTrackCuts(cutSetP);
  Int_t iCutK=task->AddTrackCuts(cutSetK);

  // monitoring
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
    AddMonitorOutput(isMC,cutSetP->GetMonitorOutput());
    AddMonitorOutput(isMC,cutSetK->GetMonitorOutput());
  }

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  if(system!=1) cutY->SetRangeD(-0.5,0.5);
  else cutY->SetRangeD(-0.465,0.035);
  AliRsnCutSet* cutsPair=new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
  cutsPair->AddCut(cutY);
  cutsPair->SetCutScheme(cutY->GetName());
    
  // multiplicity binning
  Double_t multbins[200];
  int j,nmult=0;
  if(!MultBins){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.e6; nmult++;
  }else if(!trigger){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.; nmult++;
    multbins[nmult]=5.; nmult++;
    multbins[nmult]=10.; nmult++;
    multbins[nmult]=15.; nmult++;
    for(j=2;j<=10;j++){multbins[nmult]=j*10; nmult++;}
  }else{
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=0.001; nmult++;
    multbins[nmult]=0.005; nmult++;
    multbins[nmult]=0.01; nmult++;
    multbins[nmult]=0.05; nmult++;
    multbins[nmult]=0.1; nmult++;
    multbins[nmult]=1.; nmult++;
  }

  // -- Values ------------------------------------------------------------------------------------
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass,    kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);

  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges

  Bool_t  use     [12] = {1         ,1         ,1         ,1         ,1       ,1       ,isMC     ,isMC     ,isMC     ,isMC     ,isMC    ,isMC   };
  Bool_t  useIM   [12] = {1         ,1         ,1         ,1         ,1       ,1       ,1        ,1        ,1        ,1        ,0       ,0      };
  TString name    [12] = {"UnlikePM","UnlikeMP","MixingPM","MixingMP","LikePP","LikeMM","MCGenPM","MCGenMP","TruesPM","TruesMP","ResPM" ,"ResMP"};
  TString comp    [12] = {"PAIR"    ,"PAIR"    ,"MIX"     ,"MIX"     ,"PAIR"  ,"PAIR"  ,"MOTHER" ,"MOTHER" ,"TRUE"   ,"TRUE"   ,"TRUE"  ,"TRUE" };
  Char_t  charge1 [12] = {'+'       ,'-'       ,'+'       ,'-'       ,'+'     ,'-'     ,'+'      ,'-'      ,'+'      ,'-'      ,'+'     ,'-'    };
  Char_t  charge2 [12] = {'-'       ,'+'       ,'-'       ,'+'       ,'+'     ,'-'     ,'-'      ,'+'      ,'_'      ,'+'      ,'-'     ,'+'    };
  Int_t   PDGCode [12] = {3124      ,3124      ,3124      ,3124      ,3124    ,3124    ,-3124    ,3124     ,-3124    ,3124     ,-3124   ,3124   };

  for(Int_t i=0;i<12;i++){
    if(!use[i]) continue;
    AliRsnMiniOutput *out=task->CreateOutput(Form("pkx_%s%s",name[i].Data(),suffix),"HIST",comp[i].Data());
    out->SetDaughter(0,AliRsnDaughter::kKaon);
    out->SetDaughter(1,AliRsnDaughter::kProton);
    out->SetCutID(0,iCutK);
    out->SetCutID(1,iCutP);
    out->SetCharge(0,charge1[i]);
    out->SetCharge(1,charge2[i]);
    out->SetMotherPDG(PDGCode[i]);
    out->SetMotherMass(mass);
    out->SetPairCuts(cutsPair);

    // axis X: invmass or resolution
    if(useIM[i]) out->AddAxis(imID,160,1.4,3.);
    else out->AddAxis(diffID,200,-0.02,0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID,200,0.0,20.0);
    
    // axis Z: centrality-multiplicity
    out->AddAxis(centID,nmult,multbins);
  }
  return kTRUE;
}


//=============================


Bool_t Config_pk0(
  AliRsnMiniAnalysisTask *task,
  TString     lname="pk0",
  Bool_t      isMC=kFALSE,
  Int_t       system=0,
  Int_t       EventCuts=0,
  Int_t       TrackCutsP=0,
  Int_t       TrackCutsK=0
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;
  if(system==1 || system==2) MultBins=1;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;

  Double_t mass=1.51953;
  
  // set cuts for primary proton
  if(!(TrackCutsP%10000)) TrackCutsP+=3020;//default settings
  Float_t nsigmaPTPC=0.1*(TrackCutsP%100);
  Float_t nsigmaPTOF=0.1*((TrackCutsP/100)%100);

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);
  AliRsnCutSetDaughterParticle* cutSetP=new AliRsnCutSetDaughterParticle(Form("cutProton_%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV,nsigmaPTPC),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV,AliPID::kProton,nsigmaPTPC);
  
  Int_t iCutQ=task->AddTrackCuts(cutSetQ);
  Int_t iCutP=task->AddTrackCuts(cutSetP);
    
  // selections for V0 daughters
  Int_t v0d_xrows=70;
  Float_t v0d_rtpc=0.8;
  Float_t v0d_dcaxy=0.06;
    
  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts("qualityDaughterK0s");
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetRequireTPCRefit();
  esdTrackCuts->SetAcceptKinkDaughters(0);
  esdTrackCuts->SetMinNCrossedRowsTPC(v0d_xrows);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(v0d_rtpc);
  esdTrackCuts->SetMinDCAToVertexXY(v0d_dcaxy);
    
  // selections for K0S
  Float_t k0s_piPIDCut=5.;
  Float_t k0sDaughDCA=1.;
  Float_t k0sDCA=1.e10;//0.3;
  Float_t k0s_pLife=20.;
  Float_t k0s_radiuslow=0.5;
  Float_t k0s_radiushigh=200.;
  Float_t k0s_massTol=0.03;
  Float_t k0s_massTolVeto=0.004;
  Bool_t  k0sSwitch=kFALSE;
  Float_t k0sCosPoinAn=0.97;
    
  AliRsnCutV0* cutK0s=new AliRsnCutV0("cutK0s",kK0Short,AliPID::kPion,AliPID::kPion);
  cutK0s->SetPIDCutPion(k0s_piPIDCut);// PID for the pion daughters of K0S
  cutK0s->SetESDtrackCuts(esdTrackCuts);// all the other selections (defined above) for pion daughters of K0S
  cutK0s->SetMaxDaughtersDCA(k0sDaughDCA);
  cutK0s->SetMaxDCAVertex(k0sDCA);
  cutK0s->SetfLife(k0s_pLife);
  cutK0s->SetfLowRadius(k0s_radiuslow);
  cutK0s->SetfHighRadius(k0s_radiushigh);
  cutK0s->SetTolerance(k0s_massTol);
  cutK0s->SetToleranceVeto(k0s_massTolVeto);//Rejection range for Competing V0 Rejection
  cutK0s->SetSwitch(k0sSwitch);
  cutK0s->SetMinCosPointingAngle(k0sCosPoinAn);
  cutK0s->SetMaxRapidity(2.);
  cutK0s->SetMinTPCcluster(-1);
    
  AliRsnCutSet* cutSetK0s=new AliRsnCutSet("setK0s",AliRsnTarget::kDaughter);
  cutSetK0s->AddCut(cutK0s);
  cutSetK0s->SetCutScheme(cutK0s->GetName());
  Int_t iCutK0s=task->AddTrackCuts(cutSetK0s);

  // monitoring
  TString pname="k0";
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
    AddMonitorOutput(isMC,cutSetP->GetMonitorOutput());

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
  }

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  if(system!=1) cutY->SetRangeD(-0.5,0.5);
  else cutY->SetRangeD(-0.465,0.035);

  AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
   
  AliRsnCutSet* cutsPairSame=new AliRsnCutSet("pairCutsSame",AliRsnTarget::kMother);
  cutsPairSame->AddCut(cutY);
  cutsPairSame->AddCut(cutV0);
  cutsPairSame->SetCutScheme(TString::Format("%s&(!%s)",cutY->GetName(),cutV0->GetName()).Data());
  
  AliRsnCutSet* cutsPairMix=new AliRsnCutSet("pairCutsMix", AliRsnTarget::kMother);
  cutsPairMix->AddCut(cutY);
  cutsPairMix->SetCutScheme(cutY->GetName());
    
  // multiplicity binning
  Double_t multbins[200];
  int j,nmult=0;
  if(!MultBins){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.e6; nmult++;
  }else if(!trigger){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.; nmult++;
    multbins[nmult]=5.; nmult++;
    multbins[nmult]=10.; nmult++;
    multbins[nmult]=15.; nmult++;
    for(j=2;j<=10;j++){multbins[nmult]=j*10; nmult++;}
  }else{
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=0.001; nmult++;
    multbins[nmult]=0.005; nmult++;
    multbins[nmult]=0.01; nmult++;
    multbins[nmult]=0.05; nmult++;
    multbins[nmult]=0.1; nmult++;
    multbins[nmult]=1.; nmult++;
  }

  // -- Values ------------------------------------------------------------------------------------                                    
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* transv. momentum */ Int_t ptID    = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality       */ Int_t centID  = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity   */ Int_t etaID   = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity         */ Int_t yID     = task->CreateValue(AliRsnMiniValue::kY, kFALSE);

  // -- Create all needed outputs ----------------------------------------------------------------- 
  // use an array for more compact writing, which are different on mixing and charges

  Bool_t  use     [6] = {1               ,1                ,1                  ,1                   ,0                ,0                 };
  Bool_t  useIM   [6] = {1               ,1                ,1                  ,1                   ,1                ,1                 };
  TString name    [6] = {"K0Pp"         ,"K0Pm"          ,"K0PpMix"         ,"K0PmMix"          ,"KStarPlusMinust","AKStarPlusMinust"};
  TString comp    [6] = {"PAIR"          ,"PAIR"           ,"MIX"              ,"MIX"               ,"TRUE"           ,"TRUE"            };
  TString output  [6] = {"HIST"        ,"HIST"         ,"HIST"           ,"HIST"            ,"HIST"         ,"HIST"          };
  Char_t  charge1 [6] = {'0'             ,'0'              ,'0'                ,'0'                 ,'0'              ,'0'               };
  Char_t  charge2 [6] = {'+'             ,'-'              ,'+'                ,'-'                 ,'+'              ,'-'               };
  Int_t   cutID1  [6] = { iCutK0s        ,iCutK0s           ,iCutK0s            ,iCutK0s            ,iCutK0s          ,iCutK0s           };
  Int_t   cutID2  [6] = { iCutP         ,iCutP           ,iCutP             ,iCutP              ,iCutP           ,iCutP            };
  Int_t   ipdg    [6] = {3124             ,-3124             ,3124                ,-3124                ,3124              ,-3124              };
   
  for(Int_t i=0;i<6;i++){
    if (!use[i]) continue;
    // create output
    AliRsnMiniOutput* out=task->CreateOutput(Form("pk0_%s%s",name[i].Data(),suffix),output[i].Data(),comp[i].Data());
    // selection settings
    out->SetCutID(0,cutID1[i]);
    out->SetCutID(1,cutID2[i]);
    out->SetDaughter(0,AliRsnDaughter::kKaon0);
    out->SetDaughter(1,AliRsnDaughter::kProton);
    out->SetCharge(0,charge1[i]);
    out->SetCharge(1,charge2[i]);
    out->SetMotherPDG(ipdg[i]);
    out->SetMotherMass(mass);
    // pair cuts
    if(TrackCutsK & 1024){
      out->SetPairCuts(cutsPairMix);
    }else if(TrackCutsK & 2048){
      out->SetPairCuts(cutsPairSame);
    }else{
      if(i==0 || i==1 || i==4 || i==5) out->SetPairCuts(cutsPairSame);
      else out->SetPairCuts(cutsPairMix);
    }

    // axis X: invmass or resolution
    if(useIM[i]) out->AddAxis(imID, 785, 1.43, 3.);
    else out->AddAxis(resID,200,-0.02,0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID,200,0.,20.);
    
    // axis Z: centrality-multiplicity
    out->AddAxis(centID,nmult,multbins);
   }

  return kTRUE;
}


//=============================


Bool_t Config_Lambdapi(
  AliRsnMiniAnalysisTask *task,
  TString     lname="Lambdapi",
  Bool_t      isMC=kFALSE,
  Int_t       system=0,
  Int_t       EventCuts=0,
  Int_t       TrackCutsLambda=0,
  Int_t       TrackCutsPi=0
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;
  if(system==1 || system==2) MultBins=1;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;
  
  // set cuts for primary pion
  if(!(TrackCutsPi%10000)) TrackCutsPi+=3020;//default settings
  Float_t nsigmaPiTPC=0.1*(TrackCutsPi%100);
  Float_t nsigmaPiTOF=0.1*((TrackCutsPi/100)%100);
  Int_t CutTypePi=(TrackCutsPi/10000)%100;//0=TPC+TOFveto (default), 1=TPC only, 2=TOF only, 3 TPC+TOFcut
  Int_t MisidentifiedAsKaon=(TrackCutsPi/1000000)%10;//0=pion assigned pion mass, 1=pion assigned kaon mass (for Xi(1820)- analysis)

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);

  AliRsnCutSetDaughterParticle* cutSetPi=0;
  if(!CutTypePi) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaPiTPC),
                                                           trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kPion,nsigmaPiTPC,nsigmaPiTOF);
  else if(CutTypePi==1) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,nsigmaPiTPC),
                                                                  trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kPion,nsigmaPiTPC,-1.);
  else if(CutTypePi==2) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,nsigmaPiTOF),
                                                                  trkQualityCut,AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,AliPID::kPion,-1.,nsigmaPiTOF);
  else if(CutTypePi==3) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFtightPidKStarPPB2011,nsigmaPiTPC),
                                                                  trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFtightPidKStarPPB2011,AliPID::kPion,nsigmaPiTPC,-1.);
  if(!cutSetPi){cerr<<"Error in AddTaskRare_pp13::Config_Lambdapi(): missing cutSetPi"<<endl; return kFALSE;}
  
  Int_t iCutQ=task->AddTrackCuts(cutSetQ);
  Int_t iCutPi=task->AddTrackCuts(cutSetPi);
    
  // selections for V0 daughters
  Int_t v0d_xrows=70;
  Float_t v0d_rtpc=0.8;
  Float_t v0d_dcaxy=0.06;
    
  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts("qualityDaughterK0s");
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetRequireTPCRefit();
  esdTrackCuts->SetAcceptKinkDaughters(0);
  esdTrackCuts->SetMinNCrossedRowsTPC(v0d_xrows);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(v0d_rtpc);
  esdTrackCuts->SetMinDCAToVertexXY(v0d_dcaxy);
    
  // selections for Lambda
  Float_t lambda_piPIDCut=5.;
  Float_t lambda_pPIDCut=5.;
  Float_t lambdaDaughDCA=1.0;//0.5;
  Float_t lambdaDCA=1.e10;//0.3;
  Float_t lambda_pLife=30.;
  Float_t lambda_radiuslow=0.5;
  Float_t lambda_radiushigh=200.;
  Float_t lambda_massTol=0.006;
  Float_t lambda_massTolVeto=0.004;
  Bool_t  lambdaSwitch=kFALSE;
  Float_t lambdaCosPoinAn=0.99;//0.995 for Lambda analysis
    
  // selections for the proton and pion daugthers of Lambda
    
  AliRsnCutV0* cutLambda=new AliRsnCutV0("cutLambda",kLambda0,AliPID::kProton,AliPID::kPion);
  cutLambda->SetPIDCutProton(lambda_pPIDCut); // PID for the proton daughter of Lambda
  cutLambda->SetPIDCutPion(lambda_piPIDCut);  // PID for the pion daughter of Lambda
  cutLambda->SetESDtrackCuts(esdTrackCuts);  // all the other selections (defined above) for proton and pion daughters of Lambda
  cutLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
  cutLambda->SetMaxDCAVertex(lambdaDCA);
  cutLambda->SetfLife(lambda_pLife);
  cutLambda->SetfLowRadius(lambda_radiuslow);
  cutLambda->SetfHighRadius(lambda_radiushigh);
  cutLambda->SetTolerance(lambda_massTol);
  cutLambda->SetToleranceVeto(lambda_massTolVeto);//Rejection range for Competing V0 Rejection
  cutLambda->SetSwitch(lambdaSwitch);
  cutLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
  cutLambda->SetMaxRapidity(2.);
  cutLambda->SetMinTPCcluster(-1);
    
  AliRsnCutSet* cutSetLambda=new AliRsnCutSet("setLambda",AliRsnTarget::kDaughter);
  cutSetLambda->AddCut(cutLambda);
  cutSetLambda->SetCutScheme(cutLambda->GetName());
  Int_t iCutLambda=task->AddTrackCuts(cutSetLambda);
    
  // selections for AntiLambda
  AliRsnCutV0* cutAntiLambda=new AliRsnCutV0("cutAntiLambda",kLambda0Bar,AliPID::kProton,AliPID::kPion);
  cutAntiLambda->SetPIDCutProton(lambda_pPIDCut);
  cutAntiLambda->SetPIDCutPion(lambda_piPIDCut);
  cutAntiLambda->SetESDtrackCuts(esdTrackCuts);
  cutAntiLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
  cutAntiLambda->SetMaxDCAVertex(lambdaDCA);
  cutAntiLambda->SetfLife(lambda_pLife);
  cutAntiLambda->SetfLowRadius(lambda_radiuslow);
  cutAntiLambda->SetfHighRadius(lambda_radiushigh);
  cutAntiLambda->SetTolerance(lambda_massTol);
  cutAntiLambda->SetToleranceVeto(lambda_massTolVeto);//Rejection range for Competing V0 Rejection
  cutAntiLambda->SetSwitch(lambdaSwitch);
  cutAntiLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
  cutAntiLambda->SetMaxRapidity(2.);
  cutAntiLambda->SetMinTPCcluster(-1);
    
  AliRsnCutSet* cutSetAntiLambda=new AliRsnCutSet("setAntiLambda",AliRsnTarget::kDaughter);
  cutSetAntiLambda->AddCut(cutAntiLambda);
  cutSetAntiLambda->SetCutScheme(cutAntiLambda->GetName());
  Int_t iCutAntiLambda=task->AddTrackCuts(cutSetAntiLambda);

  // monitoring
  TString pname="lambdap";
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
    AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());

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
  }

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  if(system!=1) cutY->SetRangeD(-0.5,0.5);
  else cutY->SetRangeD(-0.465,0.035);

  AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
   
  AliRsnCutSet* cutsPairSame=new AliRsnCutSet("pairCutsSame",AliRsnTarget::kMother);
  cutsPairSame->AddCut(cutY);
  cutsPairSame->AddCut(cutV0);
  cutsPairSame->SetCutScheme(TString::Format("%s&(!%s)",cutY->GetName(),cutV0->GetName()).Data());
  
  AliRsnCutSet* cutsPairMix=new AliRsnCutSet("pairCutsMix", AliRsnTarget::kMother);
  cutsPairMix->AddCut(cutY);
  cutsPairMix->SetCutScheme(cutY->GetName());
    
  // multiplicity binning
  Double_t multbins[200];
  int j,nmult=0;
  if(!MultBins){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.e6; nmult++;
  }else if(!trigger){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.; nmult++;
    multbins[nmult]=5.; nmult++;
    multbins[nmult]=10.; nmult++;
    multbins[nmult]=15.; nmult++;
    for(j=2;j<=10;j++){multbins[nmult]=j*10; nmult++;}
  }else{
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=0.001; nmult++;
    multbins[nmult]=0.005; nmult++;
    multbins[nmult]=0.01; nmult++;
    multbins[nmult]=0.05; nmult++;
    multbins[nmult]=0.1; nmult++;
    multbins[nmult]=1.; nmult++;
  }

  // -- Values ------------------------------------------------------------------------------------                                    
  /* invariant mass   */ Int_t imID    = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* transv. momentum */ Int_t ptID    = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* IM difference    */ Int_t diffID  = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
  /* centrality       */ Int_t centID  = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity   */ Int_t etaID   = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity         */ Int_t yID     = task->CreateValue(AliRsnMiniValue::kY, kFALSE);

  // -- Create all needed outputs ----------------------------------------------------------------- 
  // use an array for more compact writing, which are different on mixing and charges

  Bool_t   use     [22] = { 1          ,  1         ,  1           ,  1           ,  1            ,  1            ,  1            ,  1            ,
                               isMC         ,  isMC           ,  isMC           ,  isMC           ,  isMC            ,  isMC            ,  isMC            ,  isMC            ,
                            isMC     ,  isMC       ,  isMC            ,  isMC          ,  isMC           ,  isMC           };
  Bool_t   useIM   [22] = { 1          ,  1         ,  1           ,  1           ,  1            ,  1            ,  1            ,  1            ,
                               1            ,  1              ,  1             ,  1               ,  1               ,  1               ,  1               ,  1               ,
                            1        ,  1          ,  1               ,  1             ,                 ,                 };
  TString  name    [22] = {"LambdapPip","LambdapPim", "LambdaaPim" , "LambdaaPip" ,"LambdapPipMix","LambdapPimMix","LambdaaPimMix","LambdaaPipMix",
                           "Sigmastarpp_gen","Sigmastarmp_gen","Sigmastarma_gen","Sigmastarpa_gen","Sigmastarpp_true","Sigmastarmp_true","Sigmastarma_true","Sigmastarpa_true",
                           "Xim"     , "Xip"        ,"Lambda1520p_pip","Lambda1520_pim","Lambda1520a_pim","Lambda1520a_pip"};
  TString  comp    [22] = {"PAIR"      , "PAIR"     , "PAIR"       , "PAIR"       , "MIX"         , "MIX"         , "MIX"         , "MIX"         ,
                              "MOTHER"      , "MOTHER"        , "MOTHER"        , "MOTHER"        , "TRUE"           , "TRUE"           , "TRUE"           , "TRUE"           ,
                           "TRUE"    , "TRUE"       ,"TRUE"           , "TRUE"         ,"TRUE"           , "TRUE"          };
  Char_t   charge1 ='0';
  Char_t   charge2 [22] = {'+'         , '-'        , '-'          , '+'          , '+'           , '-'           , '-'           , '+'           ,
                              '+'           , '-'             , '-'             , '+'             , '+'              ,  '-'             ,  '-'             , '+'              ,
                           '-'       , '+'          , '+'             , '-'            , '-'             , '+'             };
  Int_t    cutID1  [22] = {  iCutLambda,  iCutLambda,iCutAntiLambda,iCutAntiLambda,  iCutLambda   ,  iCutLambda   , iCutAntiLambda, iCutAntiLambda,
                              iCutLambda    ,  iCutLambda     ,  iCutAntiLambda ,  iCutAntiLambda ,  iCutLambda      ,  iCutLambda      ,  iCutAntiLambda  ,  iCutAntiLambda  ,
                           iCutLambda,iCutAntiLambda,  iCutLambda     ,iCutLambda      ,  iCutAntiLambda ,  iCutAntiLambda };
  Int_t    cutID2 = iCutPi;
  Int_t    ipdg    [22] = { 3224      ,  3114       , -3224        , -3114        ,  3224         ,  3114         , -3224         , -3114         ,
                               3224         ,  3114           , -3224           , -3114           ,  3224            ,  3114            , -3224            , -3114            ,
                            3312     , -3312        ,  3124           ,  3124          , -3124           , -3124           };
  Double_t mass    [22] = { 1.3828    ,  1.3872     ,  1.3828      ,  1.3872      ,  1.3828       ,  1.3872       ,  1.3828       ,  1.3872       ,
                               1.3828       ,  1.3872         ,  1.3828         ,  1.3872         ,  1.3828          ,  1.3872          ,  1.3828          ,  1.3872          ,
                            1.32171  ,  1.32171     ,  1.5195         ,  1.5195        ,  1.5195         ,  1.5195         };
  Int_t    pairID  [22] = { 0         , 0           , 0            , 0            , 1             , 1             , 1             , 1             ,
                               1            , 1               , 1               , 1               , 1                , 1                , 1                , 1                ,
                            1         , 1           , 1               , 1              , 1               , 1               };
   
  for(Int_t i=0;i<22;i++){
    if(!use[i]) continue;
    // create output
    AliRsnMiniOutput *out = task->CreateOutput(Form("Lambdapi_%s%s",name[i].Data(),suffix),"HIST",comp[i].Data());
    // selection settings
    out->SetCutID(0,cutID1[i]);
    out->SetCutID(1,cutID2);
    out->SetDaughter(0,AliRsnDaughter::kLambda);
    out->SetDaughter(1,AliRsnDaughter::kPion);
    if(MisidentifiedAsKaon){
      out->SetDaughter(1,AliRsnDaughter::kKaon);
      out->SetDaughterTrue(1,AliRsnDaughter::kPion);
    }
    out->SetCharge(0,charge1);
    out->SetCharge(1,charge2[i]);
    out->SetMotherPDG(ipdg[i]);
    out->SetMotherMass(mass[i]);
    // pair cuts
    if(TrackCutsLambda & 1024){
      out->SetPairCuts(cutsPairMix);
    }else if(TrackCutsLambda & 2048){
      out->SetPairCuts(cutsPairSame);
    }else{
      if(!pairID[i]) out->SetPairCuts(cutsPairSame);
      else out->SetPairCuts(cutsPairMix);
    }

    // axis X: invmass or resolution
    if(useIM[i]) out->AddAxis(imID, 875, 1.25, 3.);
    else out->AddAxis(diffID,200,-0.02,0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID,200,0.,20.);
    
    // axis Z: centrality-multiplicity
    out->AddAxis(centID,nmult,multbins);
  }

  return kTRUE;
}


//=============================


Bool_t Config_Lambdakx(
  AliRsnMiniAnalysisTask *task,
  TString     lname="Lambdakx",
  Bool_t      isMC=kFALSE,
  Int_t       system=0,
  Int_t       EventCuts=0,
  Int_t       TrackCutsLambda=0,
  Int_t       TrackCutsK=0
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;
  if(system==1 || system==2) MultBins=1;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;
  
  // set cuts for primary kaon
  if(!(TrackCutsK%10000)) TrackCutsK+=3020;//default settings
  Float_t nsigmaKTPC=0.1*(TrackCutsK%100);
  Float_t nsigmaKTOF=0.1*((TrackCutsK/100)%100);

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);
  AliRsnCutSetDaughterParticle* cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaKTPC),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kKaon,nsigmaKTPC,nsigmaKTOF);
  
  Int_t iCutQ=task->AddTrackCuts(cutSetQ);
  Int_t iCutK=task->AddTrackCuts(cutSetK);
    
  V0Cuts=TrackCutsLambda%1024;

  // selections for V0 daughters
  Int_t v0d_xrows=70;
  Float_t v0d_rtpc=0.8;
  Float_t v0d_dcaxy=0.06;

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts("qualityDaughterK0s");
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetRequireTPCRefit();
  esdTrackCuts->SetAcceptKinkDaughters(0);
  esdTrackCuts->SetMinNCrossedRowsTPC(v0d_xrows);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(v0d_rtpc);
  esdTrackCuts->SetMinDCAToVertexXY(v0d_dcaxy);

  // selections for Lambda
  Float_t lambda_piPIDCut=5.;
  Float_t lambda_pPIDCut=5.;
  Float_t lambdaDaughDCA=1.0;//0.5;
  Float_t lambdaDCA=0.4;//1.e10 0.3
  Float_t lambda_pLife=30.;
  Float_t lambda_radiuslow=0.5;
  Float_t lambda_radiushigh=200.;
  Float_t lambda_massTol=0.006;
  Float_t lambda_massTolVeto=0.004;
  Bool_t  lambdaSwitch=kFALSE;
  Float_t lambdaCosPoinAn=0.99;//0.995 for Lambda analysis

  if(V0Cuts==1) lambdaDCA=1.e10;
  else if(V0Cuts==2) lambdaDaughDCA=0.5;

  // selections for the proton and pion daugthers of Lambda

  AliRsnCutV0* cutLambda=new AliRsnCutV0("cutLambda",kLambda0,AliPID::kProton,AliPID::kPion);
  cutLambda->SetPIDCutProton(lambda_pPIDCut); // PID for the proton daughter of Lambda
  cutLambda->SetPIDCutPion(lambda_piPIDCut);  // PID for the pion daughter of Lambda 
  cutLambda->SetESDtrackCuts(esdTrackCuts);  // all the other selections (defined above) for proton and pion daughters of Lambda
  cutLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
  cutLambda->SetMaxDCAVertex(lambdaDCA);
  cutLambda->SetfLife(lambda_pLife); 
  cutLambda->SetfLowRadius(lambda_radiuslow); 
  cutLambda->SetfHighRadius(lambda_radiushigh); 
  cutLambda->SetTolerance(lambda_massTol);
  cutLambda->SetToleranceVeto(lambda_massTolVeto);//Rejection range for Competing V0 Rejection
  cutLambda->SetSwitch(lambdaSwitch);
  cutLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
  cutLambda->SetMaxRapidity(2.);
  cutLambda->SetMinTPCcluster(-1);

  AliRsnCutSet* cutSetLambda=new AliRsnCutSet("setLambda",AliRsnTarget::kDaughter);
  cutSetLambda->AddCut(cutLambda);
  cutSetLambda->SetCutScheme(cutLambda->GetName());
  Int_t iCutLambda=task->AddTrackCuts(cutSetLambda);

  // selections for AntiLambda
  AliRsnCutV0* cutAntiLambda=new AliRsnCutV0("cutAntiLambda",kLambda0Bar,AliPID::kProton,AliPID::kPion);
  cutAntiLambda->SetPIDCutProton(lambda_pPIDCut);
  cutAntiLambda->SetPIDCutPion(lambda_piPIDCut);
  cutAntiLambda->SetESDtrackCuts(esdTrackCuts);
  cutAntiLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
  cutAntiLambda->SetMaxDCAVertex(lambdaDCA);
  cutAntiLambda->SetfLife(lambda_pLife); 
  cutAntiLambda->SetfLowRadius(lambda_radiuslow); 
  cutAntiLambda->SetfHighRadius(lambda_radiushigh);
  cutAntiLambda->SetTolerance(lambda_massTol);
  cutAntiLambda->SetToleranceVeto(lambda_massTolVeto);//Rejection range for Competing V0 Rejection
  cutAntiLambda->SetSwitch(lambdaSwitch);
  cutAntiLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
  cutAntiLambda->SetMaxRapidity(2.);
  cutAntiLambda->SetMinTPCcluster(-1);

  AliRsnCutSet* cutSetAntiLambda=new AliRsnCutSet("setAntiLambda",AliRsnTarget::kDaughter);
  cutSetAntiLambda->AddCut(cutAntiLambda);
  cutSetAntiLambda->SetCutScheme(cutAntiLambda->GetName());
  Int_t iCutAntiLambda=task->AddTrackCuts(cutSetAntiLambda);

  // monitoring
  TString pname="lambdap";
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
    AddMonitorOutput(isMC,cutSetK->GetMonitorOutput());

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
  }

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  if(system!=1) cutY->SetRangeD(-0.5,0.5);
  else cutY->SetRangeD(-0.465,0.035);

  AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
   
  AliRsnCutSet* cutsPairSame=new AliRsnCutSet("pairCutsSame",AliRsnTarget::kMother);
  cutsPairSame->AddCut(cutY);
  cutsPairSame->AddCut(cutV0);
  cutsPairSame->SetCutScheme(TString::Format("%s&(!%s)",cutY->GetName(),cutV0->GetName()).Data());
  
  AliRsnCutSet* cutsPairMix=new AliRsnCutSet("pairCutsMix", AliRsnTarget::kMother);
  cutsPairMix->AddCut(cutY);
  cutsPairMix->SetCutScheme(cutY->GetName());

  // multiplicity binning
  Double_t multbins[200];
  int j,nmult=0;
  if(!MultBins){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.e6; nmult++;
  }else if(!trigger){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.; nmult++;
    multbins[nmult]=5.; nmult++;
    multbins[nmult]=10.; nmult++;
    multbins[nmult]=15.; nmult++;
    for(j=2;j<=10;j++){multbins[nmult]=j*10; nmult++;}
  }else{
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=0.001; nmult++;
    multbins[nmult]=0.005; nmult++;
    multbins[nmult]=0.01; nmult++;
    multbins[nmult]=0.05; nmult++;
    multbins[nmult]=0.1; nmult++;
    multbins[nmult]=1.; nmult++;
  }

  // -- Values ------------------------------------------------------------------------------------                                    
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);
  /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt,kFALSE);
  /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt,kFALSE);

  // -- Create all needed outputs ----------------------------------------------------------------- 
  // use an array for more compact writing, which are different on mixing and charges
  
  Bool_t   use     [14] = { 1         ,  1        ,  1           ,  1           ,  1           ,  1           ,  1           ,  1           ,
                              isMC         , isMC         , isMC          , isMC         , isMC          , isMC         };
  Bool_t   useIM   [14] = { 1         ,  1        ,  1           ,  1           ,  1           ,  1           ,  1           ,  1           ,
                              1            , 1            , 1             , 0            , 1             , 0            };
  TString  name    [14] = {"LambdapKp","LambdapKm","LambdaaKm"   , "LambdaaKp"  ,"LambdapKpMix","LambdapKmMix","LambdaaKmMix","LambdaaKpMix",
                             "Xi1820_m_gen","Xi1820_p_gen","Xi1820_m_true","Xi1820_m_res","Xi1820_p_true","Xi1820_p_res"};
  TString  comp    [14] = {"PAIR"     , "PAIR"    , "PAIR"       , "PAIR"       , "MIX"        , "MIX"        , "MIX"        , "MIX"        ,
                             "MOTHER"      , "MOTHER"     , "TRUE"        , "TRUE"       , "TRUE"        , "TRUE"       };
  Char_t   charge1 = '0';
  Char_t   charge2 [14] = {'+'        , '-'       , '-'          , '+'          , '+'          , '-'          , '-'          , '+'          ,
                             '-'           , '+'          , '-'           , '-'          , '+'           , '+'          };
  Int_t    cutID1  [14] = { iCutLambda, iCutLambda,iCutAntiLambda,iCutAntiLambda,  iCutLambda  ,  iCutLambda  ,iCutAntiLambda,iCutAntiLambda,
                             iCutLambda    ,iCutAntiLambda,  iCutLambda   ,  iCutLambda  , iCutAntiLambda,iCutAntiLambda};
  Int_t    cutID2  = iCutK;
  Int_t    ipdg    [14] = { 123314    ,  123314   , -123314      , -123314      ,  123314      ,  123314      , -123314      , -123314      ,
                             123314        , -123314      ,  123314       , 123314       ,  -123314      ,  -123314     };
  Double_t mass    = 1.8234;
  Int_t    pairID  [14] = { 0         ,  0        ,  0           ,  0           ,  1           ,  1           ,  1           ,  1           ,
                              1            ,  1           ,  1            ,  1           ,  1            ,  1           };

  for(Int_t i=0;i<14;i++){
    if(!use[i]) continue;
    // create output
    AliRsnMiniOutput* out=task->CreateOutput(Form("Lambdakx_%s%s",name[i].Data(),suffix),"HIST",comp[i].Data());
    // selection settings
    out->SetCutID(0,cutID1[i]);
    out->SetCutID(1,cutID2);
    out->SetDaughter(0,AliRsnDaughter::kLambda);
    out->SetDaughter(1,AliRsnDaughter::kKaon);
    out->SetCharge(0,charge1);
    out->SetCharge(1,charge2[i]);
      out->SetMotherPDG(ipdg[i]);
    out->SetMotherMass(mass);
    // pair cuts
    if(TrackCutsLambda & 1024){
      out->SetPairCuts(cutsPairMix);
    }else if(TrackCutsLambda & 2048){
      out->SetPairCuts(cutsPairSame);
    }else{
      if(!pairID[i]) out->SetPairCuts(cutsPairSame);
      else out->SetPairCuts(cutsPairMix);
    }
        
    // axis X: invmass or resolution
    if(useIM[i]) out->AddAxis(imID,400,1.6,2.4);
    else out->AddAxis(diffID,200,-0.02,0.02);
        
    // axis Y: transverse momentum
    out->AddAxis(ptID,200,0.,20.);
        
    // axis Z: centrality-multiplicity
    out->AddAxis(centID,nmult,multbins);
  }
    
  if(isMC){
    //phase-space histograms
    //Xi(1820)-
    AliRsnMiniOutput* out=task->CreateOutput(Form("Xi1820_m_mother_ps%s", suffix),"HIST","TRUE");
    out->SetCutID(0,iCutLambda);
    out->SetCutID(1,cutID2);
    out->SetDaughter(0,AliRsnDaughter::kLambda);
    out->SetDaughter(1,AliRsnDaughter::kKaon);
    out->SetCharge(0,charge1);
    out->SetCharge(1,'-');
    out->SetMotherPDG(123314);
    out->SetMotherMass(mass);
    out->SetPairCuts(cutsPairMix);//just rapidity, no autocorrelation check
    out->AddAxis(fdpt,100,0.,10.);
    out->AddAxis(sdpt,100,0.,10.);
    out->AddAxis(ptID,40,0.,20.);
      
    //anti-Xi(1820)+
    AliRsnMiniOutput* out=task->CreateOutput(Form("Xi1820_p_mother_ps%s", suffix),"HIST","TRUE");
    out->SetCutID(0,iCutAntiLambda);
    out->SetCutID(1,cutID2);
    out->SetDaughter(0,AliRsnDaughter::kLambda);
    out->SetDaughter(1,AliRsnDaughter::kKaon);
    out->SetCharge(0,charge1);
    out->SetCharge(1,'+');
    out->SetMotherPDG(-123314);
    out->SetMotherMass(mass);
    out->SetPairCuts(cutsPairMix);
    out->AddAxis(fdpt,100,0.,10.);
    out->AddAxis(sdpt,100,0.,10.);
    out->AddAxis(ptID,40,0.,20.);
  }

  return kTRUE;
}


//=============================


Bool_t Config_Lambdak0(
  AliRsnMiniAnalysisTask *task,
  TString     lname="Lambdak0",
  Bool_t      isMC=kFALSE,
  Int_t       system=0,
  Int_t       EventCuts=0,
  Int_t       TrackCutsLambda=0,
  Int_t       TrackCutsK=0
){
    
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;
  if(system==1 || system==2) MultBins=1;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;
    
  V0Cuts=TrackCutsLambda%1024;

  // selections for V0 daughters
  Int_t v0d_xrows=70;
  Float_t v0d_rtpc=0.8;
  Float_t v0d_dcaxy=0.06;

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts("qualityDaughterK0s");
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetRequireTPCRefit();
  esdTrackCuts->SetAcceptKinkDaughters(0);
  esdTrackCuts->SetMinNCrossedRowsTPC(v0d_xrows);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(v0d_rtpc);
  esdTrackCuts->SetMinDCAToVertexXY(v0d_dcaxy);

  // selections for K0S 
  Float_t k0s_piPIDCut=5.;
  Float_t k0sDaughDCA=1.;
  Float_t k0sDCA=0.3;//1.e10 0.3
  Float_t k0s_pLife=20.;
  Float_t k0s_radiuslow=0.5;
  Float_t k0s_radiushigh=200.;
  Float_t k0s_massTol=0.02;//0.03
  Float_t k0s_massTolVeto=0.004;
  Bool_t  k0sSwitch=kFALSE;
  Float_t k0sCosPoinAn=0.97;
    
  if(V0Cuts==1){
    k0sDCA=1.e10;
    k0s_massTol=0.03;
  }else if(V0Cuts==2) k0sDaughDCA=0.5;
  else if(V0Cuts==3) k0s_massTol=0.03;
  else if(V0Cuts==4) task->SetCheckDecay(false);

  AliRsnCutV0* cutK0s=new AliRsnCutV0("cutK0s",kK0Short,AliPID::kPion,AliPID::kPion);
  cutK0s->SetPIDCutPion(k0s_piPIDCut);// PID for the pion daughters of K0S
  cutK0s->SetESDtrackCuts(esdTrackCuts);// all the other selections (defined above) for pion daughters of K0S
  cutK0s->SetMaxDaughtersDCA(k0sDaughDCA);
  cutK0s->SetMaxDCAVertex(k0sDCA);
  cutK0s->SetfLife(k0s_pLife);
  cutK0s->SetfLowRadius(k0s_radiuslow);
  cutK0s->SetfHighRadius(k0s_radiushigh);
  cutK0s->SetTolerance(k0s_massTol);
  cutK0s->SetToleranceVeto(k0s_massTolVeto);//Rejection range for Competing V0 Rejection
  cutK0s->SetSwitch(k0sSwitch);
  cutK0s->SetMinCosPointingAngle(k0sCosPoinAn);
  cutK0s->SetMaxRapidity(2.);
  cutK0s->SetMinTPCcluster(-1);

  AliRsnCutSet* cutSetK0s=new AliRsnCutSet("setK0s",AliRsnTarget::kDaughter);
  cutSetK0s->AddCut(cutK0s);
  cutSetK0s->SetCutScheme(cutK0s->GetName());
  Int_t iCutK0s=task->AddTrackCuts(cutSetK0s);

  // selections for Lambda
  Float_t lambda_piPIDCut=5.;
  Float_t lambda_pPIDCut=5.;
  Float_t lambdaDaughDCA=1.;//0.5
  Float_t lambdaDCA=0.4;//1.e10 0.3
  Float_t lambda_pLife=30.;
  Float_t lambda_radiuslow=0.5;
  Float_t lambda_radiushigh=200.;
  Float_t lambda_massTol=0.006;
  Float_t lambda_massTolVeto=0.004;
  Bool_t  lambdaSwitch=kFALSE;
  Float_t lambdaCosPoinAn=0.99;//0.995 for Lambda analysis
    
  if(V0Cuts==1) lambdaDCA=1.e10;
  else if(V0Cuts==2) lambdaDaughDCA=0.5;

  // selections for the proton and pion daugthers of Lambda

  AliRsnCutV0* cutLambda=new AliRsnCutV0("cutLambda",kLambda0,AliPID::kProton,AliPID::kPion);
  cutLambda->SetPIDCutProton(lambda_pPIDCut); // PID for the proton daughter of Lambda
  cutLambda->SetPIDCutPion(lambda_piPIDCut);  // PID for the pion daughter of Lambda 
  cutLambda->SetESDtrackCuts(esdTrackCuts);  // all the other selections (defined above) for proton and pion daughters of Lambda
  cutLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
  cutLambda->SetMaxDCAVertex(lambdaDCA);
  cutLambda->SetfLife(lambda_pLife); 
  cutLambda->SetfLowRadius(lambda_radiuslow); 
  cutLambda->SetfHighRadius(lambda_radiushigh); 
  cutLambda->SetTolerance(lambda_massTol);
  cutLambda->SetToleranceVeto(lambda_massTolVeto);//Rejection range for Competing V0 Rejection
  cutLambda->SetSwitch(lambdaSwitch);
  cutLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
  cutLambda->SetMaxRapidity(2.);
  cutLambda->SetMinTPCcluster(-1);

  AliRsnCutSet* cutSetLambda=new AliRsnCutSet("setLambda",AliRsnTarget::kDaughter);
  cutSetLambda->AddCut(cutLambda);
  cutSetLambda->SetCutScheme(cutLambda->GetName());
  Int_t iCutLambda=task->AddTrackCuts(cutSetLambda);

  // selections for AntiLambda
  AliRsnCutV0* cutAntiLambda=new AliRsnCutV0("cutAntiLambda",kLambda0Bar,AliPID::kProton,AliPID::kPion);
  cutAntiLambda->SetPIDCutProton(lambda_pPIDCut);
  cutAntiLambda->SetPIDCutPion(lambda_piPIDCut);
  cutAntiLambda->SetESDtrackCuts(esdTrackCuts);
  cutAntiLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
  cutAntiLambda->SetMaxDCAVertex(lambdaDCA);
  cutAntiLambda->SetfLife(lambda_pLife); 
  cutAntiLambda->SetfLowRadius(lambda_radiuslow); 
  cutAntiLambda->SetfHighRadius(lambda_radiushigh);
  cutAntiLambda->SetTolerance(lambda_massTol);
  cutAntiLambda->SetToleranceVeto(lambda_massTolVeto);//Rejection range for Competing V0 Rejection
  cutAntiLambda->SetSwitch(lambdaSwitch);
  cutAntiLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
  cutAntiLambda->SetMaxRapidity(2.);
  cutAntiLambda->SetMinTPCcluster(-1);

  AliRsnCutSet* cutSetAntiLambda=new AliRsnCutSet("setAntiLambda",AliRsnTarget::kDaughter);
  cutSetAntiLambda->AddCut(cutAntiLambda);
  cutSetAntiLambda->SetCutScheme(cutAntiLambda->GetName());
  Int_t iCutAntiLambda=task->AddTrackCuts(cutSetAntiLambda);

  // monitoring
  TString pname="lambdap";
  if(enableMonitor){
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
  }

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  if(system!=1) cutY->SetRangeD(-0.5,0.5);
  else cutY->SetRangeD(-0.465,0.035);

  AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
   
  AliRsnCutSet* cutsPairSame=new AliRsnCutSet("pairCutsSame",AliRsnTarget::kMother);
  cutsPairSame->AddCut(cutY);
  cutsPairSame->AddCut(cutV0);
  cutsPairSame->SetCutScheme(TString::Format("%s&(!%s)",cutY->GetName(),cutV0->GetName()).Data());
  
  AliRsnCutSet* cutsPairMix=new AliRsnCutSet("pairCutsMix", AliRsnTarget::kMother);
  cutsPairMix->AddCut(cutY);
  cutsPairMix->SetCutScheme(cutY->GetName());
    
  // multiplicity binning
  Double_t multbins[200];
  int j,nmult=0;
  if(!MultBins){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.e6; nmult++;
  }else if(!trigger){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.; nmult++;
    multbins[nmult]=5.; nmult++;
    multbins[nmult]=10.; nmult++;
    multbins[nmult]=15.; nmult++;
    for(j=2;j<=10;j++){multbins[nmult]=j*10; nmult++;}
  }else{
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=0.001; nmult++;
    multbins[nmult]=0.005; nmult++;
    multbins[nmult]=0.01; nmult++;
    multbins[nmult]=0.05; nmult++;
    multbins[nmult]=0.1; nmult++;
    multbins[nmult]=1.; nmult++;
  }

  // -- Values ------------------------------------------------------------------------------------                                    
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);
  /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt,kFALSE);
  /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt,kFALSE);

  // -- Create all needed outputs ----------------------------------------------------------------- 
  // use an array for more compact writing, which are different on mixing and charges
   
  Bool_t   use     [10] = { 1         ,  1           ,  1           ,  1           ,  isMC         ,  isMC         ,  isMC          ,  isMC         ,  isMC          ,  isMC         };
  Bool_t   useIM   [10] = { 1         ,  1           ,  1           ,  1           ,  1            ,  1            ,  1             ,  0            ,  1             ,  0            };
  TString  name    [10] = {"LambdapK0", "LambdaaK0"  ,"LambdapK0Mix","LambdaaK0Mix","Xi1820_0p_gen","Xi1820_0a_gen","Xi1820_0p_true","Xi1820_0p_res","Xi1820_0a_true","Xi1820_0a_res"};
  TString  comp    [10] = {"PAIR"     , "PAIR"       ,  "MIX"       , "MIX"        , "MOTHER"      , "MOTHER"      , "TRUE"         , "TRUE"        , "TRUE"         , "TRUE"        };
  Char_t   charge1 = '0';
  Char_t   charge2 = '0';
  Int_t    cutID1  [10] = {iCutLambda ,iCutAntiLambda, iCutLambda   ,iCutAntiLambda,  iCutLambda   , iCutAntiLambda,  iCutLambda    ,  iCutLambda   , iCutAntiLambda , iCutAntiLambda};
  Int_t    cutID2  = iCutK0s;
  Int_t    ipdg    [10] = { 123324    ,  -123324     , 123324       , -123324      , 123324        , -123324       , 123324         , 123324        , -123324        , -123324       };
  Double_t mass    = 1.8234;
  Int_t    pairID  [10] = { 0         ,  0           ,  1           ,  1           ,  1            ,  1            ,  1             ,  1            ,  1             ,  1            };
    
  for(Int_t i=0;i<10;i++){
    if(!use[i]) continue;
    // create output
    AliRsnMiniOutput* out=task->CreateOutput(Form("Lambdak0_%s%s",name[i].Data(),suffix),"HIST",comp[i].Data());
    // selection settings
    out->SetCutID(0,cutID1[i]);
    out->SetCutID(1,cutID2);
    out->SetDaughter(0,AliRsnDaughter::kLambda);
    out->SetDaughter(1,AliRsnDaughter::kKaon0);
    out->SetCharge(0,charge1);
    out->SetCharge(1,charge2);
    out->SetMotherPDG(ipdg[i]);
    out->SetMotherMass(mass);
    // pair cuts
    if(TrackCutsLambda & 1024){
      out->SetPairCuts(cutsPairMix);
    }else if(TrackCutsLambda & 2048){
      out->SetPairCuts(cutsPairSame);
    }else{
      if(!pairID[i]) out->SetPairCuts(cutsPairSame);
      else out->SetPairCuts(cutsPairMix);
    }
        
    // axis X: invmass or resolution
    if(useIM[i]) out->AddAxis(imID,400,1.6,2.4);
    else out->AddAxis(diffID,200,-0.02,0.02);
        
    // axis Y: transverse momentum
    out->AddAxis(ptID,200,0.,20.);
        
    // axis Z: centrality-multiplicity
    out->AddAxis(centID,nmult,multbins);
  }
    
  if(isMC){
    //phase-space histograms
    //Xi(1820)0
    AliRsnMiniOutput* out=task->CreateOutput(Form("Xi1820_0p_mother_ps%s", suffix),"HIST","TRUE");
    out->SetCutID(0,iCutLambda);
    out->SetCutID(1,cutID2);
    out->SetDaughter(0,AliRsnDaughter::kLambda);
    out->SetDaughter(1,AliRsnDaughter::kKaon0);
    out->SetCharge(0,charge1);
    out->SetCharge(1,charge2);
    out->SetMotherPDG(123324);
    out->SetMotherMass(mass);
    out->SetPairCuts(cutsPairMix);//just rapidity, no autocorrelation check
    out->AddAxis(fdpt,100,0.,10.);
    out->AddAxis(sdpt,100,0.,10.);
    out->AddAxis(ptID,40,0.,20.);
      
    //anti-Xi(1820)0
    AliRsnMiniOutput* out=task->CreateOutput(Form("Xi1820_0a_mother_ps%s", suffix),"HIST","TRUE");
    out->SetCutID(0,iCutAntiLambda);
    out->SetCutID(1,cutID2);
    out->SetDaughter(0,AliRsnDaughter::kLambda);
    out->SetDaughter(1,AliRsnDaughter::kKaon0);
    out->SetCharge(0,charge1);
    out->SetCharge(1,charge2);
    out->SetMotherPDG(-123324);
    out->SetMotherMass(mass);
    out->SetPairCuts(cutsPairMix);
    out->AddAxis(fdpt,100,0.,10.);
    out->AddAxis(sdpt,100,0.,10.);
    out->AddAxis(ptID,40,0.,20.);
  }

  return kTRUE;
}


//=============================


Bool_t Config_Lambdap(
  AliRsnMiniAnalysisTask *task,
  TString     lname="Lambdap",
  Bool_t      isMC=kFALSE,
  Int_t       system=0,
  Int_t       EventCuts=0,
  Int_t       TrackCutsLambda=0,
  Int_t       TrackCutsP=0
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;
  if(system==1 || system==2) MultBins=1;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;
  
  // set cuts for primary proton
  if(!(TrackCutsP%10000)) TrackCutsP+=3020;//default settings
  Float_t nsigmaPTPC=0.1*(TrackCutsP%100);
  Float_t nsigmaPTOF=0.1*((TrackCutsP/100)%100);

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);
  AliRsnCutSetDaughterParticle* cutSetP=new AliRsnCutSetDaughterParticle(Form("cutProton_%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV,nsigmaPTPC),trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV,AliPID::kProton,nsigmaPTPC);
  
  Int_t iCutQ=task->AddTrackCuts(cutSetQ);
  Int_t iCutP=task->AddTrackCuts(cutSetP);
    
  // selections for V0 daughters
  Int_t v0d_xrows=70;
  Float_t v0d_rtpc=0.8;
  Float_t v0d_dcaxy=0.06;
    
  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts("qualityDaughterK0s");
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetRequireTPCRefit();
  esdTrackCuts->SetAcceptKinkDaughters(0);
  esdTrackCuts->SetMinNCrossedRowsTPC(v0d_xrows);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(v0d_rtpc);
  esdTrackCuts->SetMinDCAToVertexXY(v0d_dcaxy);
    
  // selections for Lambda
  Float_t lambda_piPIDCut=5.;
  Float_t lambda_pPIDCut=5.;
  Float_t lambdaDaughDCA=1.0;//0.5;
  Float_t lambdaDCA=1.e10;//0.3;
  Float_t lambda_pLife=30.;
  Float_t lambda_radiuslow=0.5;
  Float_t lambda_radiushigh=200.;
  Float_t lambda_massTol=0.006;
  Float_t lambda_massTolVeto=0.004;
  Bool_t  lambdaSwitch=kFALSE;
  Float_t lambdaCosPoinAn=0.99;//0.995 for Lambda analysis
    
  // selections for the proton and pion daugthers of Lambda
    
  AliRsnCutV0* cutLambda=new AliRsnCutV0("cutLambda",kLambda0,AliPID::kProton,AliPID::kPion);
  cutLambda->SetPIDCutProton(lambda_pPIDCut); // PID for the proton daughter of Lambda
  cutLambda->SetPIDCutPion(lambda_piPIDCut);  // PID for the pion daughter of Lambda
  cutLambda->SetESDtrackCuts(esdTrackCuts);  // all the other selections (defined above) for proton and pion daughters of Lambda
  cutLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
  cutLambda->SetMaxDCAVertex(lambdaDCA);
  cutLambda->SetfLife(lambda_pLife);
  cutLambda->SetfLowRadius(lambda_radiuslow);
  cutLambda->SetfHighRadius(lambda_radiushigh);
  cutLambda->SetTolerance(lambda_massTol);
  cutLambda->SetToleranceVeto(lambda_massTolVeto);//Rejection range for Competing V0 Rejection
  cutLambda->SetSwitch(lambdaSwitch);
  cutLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
  cutLambda->SetMaxRapidity(2.);
  cutLambda->SetMinTPCcluster(-1);
    
  AliRsnCutSet* cutSetLambda=new AliRsnCutSet("setLambda",AliRsnTarget::kDaughter);
  cutSetLambda->AddCut(cutLambda);
  cutSetLambda->SetCutScheme(cutLambda->GetName());
  Int_t iCutLambda=task->AddTrackCuts(cutSetLambda);
    
  // selections for AntiLambda
  AliRsnCutV0* cutAntiLambda=new AliRsnCutV0("cutAntiLambda",kLambda0Bar,AliPID::kProton,AliPID::kPion);
  cutAntiLambda->SetPIDCutProton(lambda_pPIDCut);
  cutAntiLambda->SetPIDCutPion(lambda_piPIDCut);
  cutAntiLambda->SetESDtrackCuts(esdTrackCuts);
  cutAntiLambda->SetMaxDaughtersDCA(lambdaDaughDCA);
  cutAntiLambda->SetMaxDCAVertex(lambdaDCA);
  cutAntiLambda->SetfLife(lambda_pLife);
  cutAntiLambda->SetfLowRadius(lambda_radiuslow);
  cutAntiLambda->SetfHighRadius(lambda_radiushigh);
  cutAntiLambda->SetTolerance(lambda_massTol);
  cutAntiLambda->SetToleranceVeto(lambda_massTolVeto);//Rejection range for Competing V0 Rejection
  cutAntiLambda->SetSwitch(lambdaSwitch);
  cutAntiLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
  cutAntiLambda->SetMaxRapidity(2.);
  cutAntiLambda->SetMinTPCcluster(-1);
    
  AliRsnCutSet* cutSetAntiLambda=new AliRsnCutSet("setAntiLambda",AliRsnTarget::kDaughter);
  cutSetAntiLambda->AddCut(cutAntiLambda);
  cutSetAntiLambda->SetCutScheme(cutAntiLambda->GetName());
  Int_t iCutAntiLambda=task->AddTrackCuts(cutSetAntiLambda);

  // monitoring
  TString pname="lambdap";
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
    AddMonitorOutput(isMC,cutSetP->GetMonitorOutput());

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
  }

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  if(system!=1) cutY->SetRangeD(-0.5,0.5);
  else cutY->SetRangeD(-0.465,0.035);

  AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
   
  AliRsnCutSet* cutsPairSame=new AliRsnCutSet("pairCutsSame",AliRsnTarget::kMother);
  cutsPairSame->AddCut(cutY);
  cutsPairSame->AddCut(cutV0);
  cutsPairSame->SetCutScheme(TString::Format("%s&(!%s)",cutY->GetName(),cutV0->GetName()).Data());
  
  AliRsnCutSet* cutsPairMix=new AliRsnCutSet("pairCutsMix", AliRsnTarget::kMother);
  cutsPairMix->AddCut(cutY);
  cutsPairMix->SetCutScheme(cutY->GetName());
    
  // multiplicity binning
  Double_t multbins[200];
  int j,nmult=0;
  if(!MultBins){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.e6; nmult++;
  }else if(!trigger){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.; nmult++;
    multbins[nmult]=5.; nmult++;
    multbins[nmult]=10.; nmult++;
    multbins[nmult]=15.; nmult++;
    for(j=2;j<=10;j++){multbins[nmult]=j*10; nmult++;}
  }else{
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=0.001; nmult++;
    multbins[nmult]=0.005; nmult++;
    multbins[nmult]=0.01; nmult++;
    multbins[nmult]=0.05; nmult++;
    multbins[nmult]=0.1; nmult++;
    multbins[nmult]=1.; nmult++;
  }

  // -- Values ------------------------------------------------------------------------------------                                    
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* transv. momentum */ Int_t ptID    = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality       */ Int_t centID  = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity   */ Int_t etaID   = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity         */ Int_t yID     = task->CreateValue(AliRsnMiniValue::kY, kFALSE);

  // -- Create all needed outputs ----------------------------------------------------------------- 
  // use an array for more compact writing, which are different on mixing and charges

  Bool_t   use     [18] = { 1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  0         ,  0         ,  0             ,  0             ,  0          ,  0              ,  0              ,  0              ,  0              ,  0              };
  Bool_t   useIM   [18] = { 1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1         ,  1         ,  1             ,  1             ,  1          ,  1              ,  1              ,  1              ,  1              ,  1              };
  TString  name    [18] = {"LambdapPp"   , "LambdapPm"   , "LambdaaPm"      , "LambdaaPp"      , "LambdapPpMix", "LambdapPmMix", "LambdaaPmMix"   , "LambdaaPpMix"   , "SigmaPt"  , "SigmaMt"  , "ASigmaPt"     , "ASigmaMt"     , "XiM"       , "XiP"           , "Lambda1520P"   , "Lambda1520M"   , "Lambda1520PBar", "Lambda1520MBar"};
  TString  comp    [18] = {"PAIR"     , "PAIR"     , "PAIR"         , "PAIR"         , "MIX"      , "MIX"      , "MIX"          , "MIX"          , "TRUE"     , "TRUE"     , "TRUE"         , "TRUE"         , "TRUE"      , "TRUE"          , "TRUE"          , "TRUE"          , "TRUE"          , "TRUE"          };
  TString  output  [18] = {"HIST"     , "HIST"     , "HIST"         , "HIST"         , "HIST"     , "HIST"     , "HIST"         , "HIST"         , "HIST"     , "HIST"     , "HIST"         , "HIST"         , "HIST"      , "HIST"          , "HIST"          , "HIST"          , "HIST"          , "HIST"          };
  Char_t   charge1 [18] = {'0'        , '0'        , '0'            , '0'            , '0'        , '0'        , '0'            , '0'            , '0'        , '0'        , '0'            , '0'            , '0'         , '0'             , '0'             , '0'             , '0'             , '0'             };
  Char_t   charge2 [18] = {'+'        , '-'        , '-'            , '+'            , '+'        , '-'        , '-'            , '+'            , '+'        , '-'        , '-'            , '+'            , '-'         , '+'             , '+'             , '-'             , '-'             , '+'             };
  Int_t    cutID1  [18] = { iCutLambda,  iCutLambda,  iCutAntiLambda,  iCutAntiLambda,  iCutLambda,  iCutLambda,  iCutAntiLambda,  iCutAntiLambda,  iCutLambda,  iCutLambda,  iCutAntiLambda,  iCutAntiLambda,  iCutLambda ,  iCutAntiLambda ,  iCutLambda     ,  iCutLambda     ,  iCutAntiLambda ,  iCutAntiLambda };
  Int_t    cutID2  [18] = { iCutP    ,  iCutP    ,  iCutP        ,  iCutP        ,  iCutP    ,  iCutP    ,  iCutP        ,  iCutP        ,  iCutP    ,  iCutP    ,  iCutP        ,  iCutP        ,  iCutP     ,  iCutP         ,  iCutP         ,  iCutP         ,  iCutP         ,  iCutP         };
  Int_t    ipdg    [18] = { 3224      ,  3114      , -3224          , -3114          ,  3224      ,  3114      , -3224          , -3114          ,  3224      ,  3114      , -3224          , -3114          ,  3312       , -3312           ,  3124           ,  3124           , -3124           , -3124           };
  Double_t mass    [18] = { 1.3828    ,  1.3872    ,  1.3828        ,  1.3872        ,  1.3828    ,  1.3872    ,  1.3828        ,  1.3872        ,  1.3828    ,  1.3872    ,  1.3828        ,  1.3872        ,  1.32171    ,  1.32171        ,  1.5195         ,  1.5195         ,  1.5195         ,  1.5195         };
   
  for(Int_t i=0;i<18;i++){
    if(!use[i]) continue;
    // create output
    AliRsnMiniOutput* out=task->CreateOutput(Form("Lambdap_%s%s",name[i].Data(),suffix),output[i].Data(),comp[i].Data());
    // selection settings
    out->SetCutID(0,cutID1[i]);
    out->SetCutID(1,cutID2[i]);
    out->SetDaughter(0,AliRsnDaughter::kLambda);
    out->SetDaughter(1,AliRsnDaughter::kProton);
    out->SetCharge(0,charge1[i]);
    out->SetCharge(1,charge2[i]);
    out->SetMotherPDG(ipdg[i]);
    out->SetMotherMass(mass[i]);
    // pair cuts
    if(TrackCutsLambda & 1024){
      out->SetPairCuts(cutsPairMix);
    }else if(TrackCutsLambda & 2048){
      out->SetPairCuts(cutsPairSame);
    }else{
      if(!(i>=4 && i<=7)) out->SetPairCuts(cutsPairSame);
      else out->SetPairCuts(cutsPairMix);
    }

    // axis X: invmass or resolution
    if(useIM[i]) out->AddAxis(imID,975,2.05,4.);
    else out->AddAxis(resID,200,-0.02,0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID,200,0.,20.);
    
    // axis Z: centrality-multiplicity
    out->AddAxis(centID,nmult,multbins);
  }

  return kTRUE;
}


//=============================


void AddMonitorOutput_P(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_mom",name.Data()),AliRsnValueDaughter::kP);
  a->SetBins(0.,10.0,0.05);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}


void AddMonitorOutput_Pt(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_pt",name.Data()),AliRsnValueDaughter::kPt);
  a->SetBins(0.,10.0,0.05);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_Eta(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_eta",name.Data()),AliRsnValueDaughter::kEta);
  a->SetBins(-2.,2.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_DCAxy(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_dcaxy",name.Data()),AliRsnValueDaughter::kDCAXY);
  a->SetBins(-0.5,0.5,0.001);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_DCAz(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_dcaz",name.Data()),AliRsnValueDaughter::kDCAZ);
  a->SetBins(-2.5,2.5,0.005);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_TPCpi(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_TPCpi",name.Data()),AliRsnValueDaughter::kTPCnsigmaPi);
  a->SetBins(-10.,10.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_TPCK(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_TPCK",name.Data()),AliRsnValueDaughter::kTPCnsigmaK);
  a->SetBins(-10.,10.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_TPCp(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_TPCp",name.Data()),AliRsnValueDaughter::kTPCnsigmaP);
  a->SetBins(-10.,10.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_NclTPC(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_NclTPC",name.Data()),AliRsnValueDaughter::kNTPCclusters);
  a->SetBins(-0.5,199.5,1);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_chi2TPC(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_chi2TPC",name.Data()),AliRsnValueDaughter::kTPCchi2);
  a->SetBins(0.0,6,.1);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0NPt(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0npt",name.Data()),AliRsnValueDaughter::kV0NPt);
  a->SetBins(0.,10.0,0.05);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0PPt(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0ppt",name.Data()),AliRsnValueDaughter::kV0PPt);
  a->SetBins(0.,10.0,0.05);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0Mass(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
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

void AddMonitorOutput_V0DCA(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0dca",name.Data()),AliRsnValueDaughter::kV0DCA);
  a->SetBins(0.0,0.4,0.001);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0Radius(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0radius",name.Data()),AliRsnValueDaughter::kV0Radius);
  a->SetBins(0.0,200,0.2);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0Lifetime(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0lifetime",name.Data()),AliRsnValueDaughter::kV0Lifetime);
  a->SetBins(0.0,200,0.1);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0DaughterDCA(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0ddca",name.Data()),AliRsnValueDaughter::kDaughterDCA);
  a->SetBins(0.0,2,0.001);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0DCA2TPV(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  //DCA of secondary tracks to primary vertex
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0dca2tpv",name.Data()),AliRsnValueDaughter::kV0DCAXY);
  a->SetBins(-10.,10.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0CPA(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0cpa",name.Data()),AliRsnValueDaughter::kCosPointAng);
  a->SetBins(0.96,1.,0.001);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0TPCpim(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0TPCpim",name.Data()),AliRsnValueDaughter::kLambdaPionPIDCut);
  a->SetBins(0.,5.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0TPCpip(TString name="",TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *loop=0){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0TPCpip",name.Data()),AliRsnValueDaughter::kAntiLambdaAntiPionPIDCut);
  a->SetBins(-0.,5.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_LambdaProtonPID(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lpPID=0){
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

void AddMonitorOutput_LambdaAntiProtonPID(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *lapPID=0){
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
