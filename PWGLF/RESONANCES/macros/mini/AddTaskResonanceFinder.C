/***************************************************************************
              Anders Knospe: anders.knospe@cern.ch
  Macro to configure the resonance package for analyses using
  AliRsnMiniResonanceFinder (where one decay product is itself a resonance).

****************************************************************************/


Bool_t Config_piphi(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
Bool_t Config_kxphi(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
Bool_t Config_k0phi(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
Bool_t Config_pphi(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
Bool_t Config_phiphi(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
Bool_t Config_Lambdaphi(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);


void AddMonitorOutput_P(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_Pt(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_Eta(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_DCAxy(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_DCAz(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_TPCpi(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_TPCK(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_TPCp(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_NclTPC(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_chi2TPC(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_V0NPt(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_V0PPt(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_V0Mass(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_V0DCA(TString n="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_V0Radius(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_V0Lifetime(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_V0DaughterDCA(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_V0DCA2TPV(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_V0CPA(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_V0TPCpim(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_V0TPCpip(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_LambdaProtonPID(TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_LambdaAntiProtonPID(TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);


AliRsnMiniAnalysisTask* AddTaskResonanceFinder(
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
    ::Error("AddTaskResonanceFinder", "No analysis manager to connect to.");
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
  if(isPP){
    if(MultBins==1) task->UseMultiplicity("AliMultSelection_V0M");
    else if(MultBins==2) task->UseMultiplicity("AliMultSelection_RefMult08");
    else task->UseMultiplicity("QUALITY");
  }else task->UseCentrality("V0M");

  // set event mixing options
  int nmix=5;
  float maxDiffVzMix=1;
  float maxDiffMultMix=5;
  task->UseContinuousMix();
  task->SetNMix(nmix);
  task->SetMaxDiffVz(maxDiffVzMix);
  task->SetMaxDiffMult(maxDiffMultMix);
  ::Info("AddTaskResonanceFinder", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f", nmix, maxDiffVzMix, maxDiffMultMix));

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
    ::Info("AddTaskResonanceFinder", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp)?"ON":"OFF"));
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

  //vertex
  Int_t vtxID=task->CreateValue(AliRsnMiniValue::kVz,kFALSE);
  AliRsnMiniOutput* outVtx=task->CreateOutput("eventVtx","HIST","EVENT");
  outVtx->AddAxis(vtxID,240,-12.0,12.0);

  //multiplicity or centrality
  Int_t multID=task->CreateValue(AliRsnMiniValue::kMult,kFALSE);
  AliRsnMiniOutput* outMult=task->CreateOutput("eventMult","HIST","EVENT");
  if(isPP && !MultBins) outMult->AddAxis(multID,400,0.5,400.5);
  else outMult->AddAxis(multID,110,0.,110.);

  Double_t multbins[200];
  int j,nmult=0;
  for(j=0;j<10;j++){multbins[nmult]=0.0001*j; nmult++;}
  for(j=1;j<10;j++){multbins[nmult]=0.001*j; nmult++;}
  for(j=1;j<10;j++){multbins[nmult]=0.01*j; nmult++;}
  for(j=1;j<10;j++){multbins[nmult]=0.1*j; nmult++;}
  for(j=1;j<=100;j++){multbins[nmult]=j; nmult++;}
  nmult--;
  TH1F* hEventsVsMulti=new TH1F("hAEventsVsMulti","",nmult,multbins);
  task->SetEventQAHist("EventsVsMulti",hEventsVsMulti);//custom binning for fHAEventsVsMulti

  TH2F* hvz=new TH2F("hVzVsCent","",110,0.,110., 240,-12.0,12.0);
  task->SetEventQAHist("vz",hvz);//plugs this histogram into the fHAEventVz data member

  double ybins[500];
  for(j=0;j<=401;j++) ybins[j]=j-0.5;

  TH2F* hmc=new TH2F("MultiVsCent","", nmult,multbins, 401,ybins);
  hmc->GetYaxis()->SetTitle("QUALITY");
  task->SetEventQAHist("multicent",hmc);//plugs this histogram into the fHAEventMultiCent data member

  // ----- CONFIGURE -----

  cerr<<"configuring"<<endl;
  if(d1==AliRsnDaughter::kPion && d2==AliRsnDaughter::kPhi){
    Config_piphi(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kPhi && d1==AliRsnDaughter::kPion){
    Config_piphi(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
        
  }else if(d1==AliRsnDaughter::kKaon && d2==AliRsnDaughter::kPhi){
    Config_kxphi(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kPhi && d1==AliRsnDaughter::kKaon){
    Config_kxphi(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
      
  }else if(d1==AliRsnDaughter::kKaon0 && d2==AliRsnDaughter::kPhi){
    Config_k0phi(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kPhi && d1==AliRsnDaughter::kKaon0){
    Config_k0phi(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kProton && d2==AliRsnDaughter::kPhi){
    Config_pphi(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kPhi && d1==AliRsnDaughter::kProton){
    Config_pphi(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
      
  }else if(d1==AliRsnDaughter::kPhi && d2==AliRsnDaughter::kPhi){
    Config_phiphi(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);

  }else if(d1==AliRsnDaughter::kLambda && d2==AliRsnDaughter::kPhi){
    Config_Lambdaphi(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kPhi && d1==AliRsnDaughter::kLambda){
    Config_Lambdaphi(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }
  cerr<<"done configuring"<<endl;
  
  // ----- CONTAINERS -----

  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  Printf("AddTaskResonanceFinder - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
  AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",lname.Data()), 
							  TList::Class(), 
							  AliAnalysisManager::kOutputContainer, 
							  outputFileName);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);
   
  return task;
}


//=============================


Bool_t Config_piphi(
  AliRsnMiniAnalysisTask *task,
  TString     lname,
  Bool_t      isMC,
  Int_t       system,
  Int_t       EventCuts,
  Int_t       TrackCutsPi,
  Int_t       TrackCutsPhi
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;

  Double_t mass=0.139571+1.019460;
  
  // set daughter cuts
  if(!(TrackCutsPi%10000)) TrackCutsPi+=3020;//default settings
  Float_t nsigmaPiTPC=0.1*(TrackCutsPi%100);
  Float_t nsigmaPiTOF=0.1*((TrackCutsPi/100)%100);
  Int_t CutTypePi=(TrackCutsPi/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only

  if(!(TrackCutsPhi%10000)) TrackCutsPhi+=3020;//default settings
  Float_t nsigmaKTPC=0.1*(TrackCutsPhi%100);
  Float_t nsigmaKTOF=0.1*((TrackCutsPhi/100)%100);
  Int_t CutTypeK=(TrackCutsPhi/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only
  Int_t SideBand=(TrackCutsPhi/100000)%10;

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,
    AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);

  AliRsnCutSetDaughterParticle* cutSetPi=0;
  if(!CutTypePi) cutSetPi=new AliRsnCutSetDaughterParticle(
    Form("cutPion_%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaPiTPC),
    trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kPion,nsigmaPiTPC,nsigmaPiTOF);
  else if(CutTypePi==1) cutSetPi=new AliRsnCutSetDaughterParticle(
    Form("cutPion%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,nsigmaPiTPC),
    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kPion,nsigmaPiTPC,-1.);
  else if(CutTypePi==2) cutSetPi=new AliRsnCutSetDaughterParticle(
    Form("cutPion%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,nsigmaPiTOF),
    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,AliPID::kPion,-1.,nsigmaPiTOF);
  if(!cutSetPi){cerr<<"Error in AddTaskResonanceFinder::Config_piphi(): missing cutSetPi"<<endl; return kFALSE;}

  AliRsnCutSetDaughterParticle* cutSetK=0;
  if(!CutTypeK) cutSetK=new AliRsnCutSetDaughterParticle(
    Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaKTPC),
    trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kKaon,nsigmaKTPC,nsigmaKTOF);
  else if(CutTypeK==1) cutSetK=new AliRsnCutSetDaughterParticle(
    Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,nsigmaKTPC),
    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kKaon,nsigmaKTPC,-1.);
  else if(CutTypeK==2) cutSetK=new AliRsnCutSetDaughterParticle(
    Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,nsigmaKTOF),
    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,AliPID::kKaon,-1.,nsigmaKTOF);
  if(!cutSetK){cerr<<"Error in AddTaskResonanceFinder::Config_piphi(): missing cutSetK"<<endl; return kFALSE;}

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
    
  // AliRsnMiniResonanceFinder
  AliRsnCutMiniPair* cutMassPhi=new AliRsnCutMiniPair("cutMassPhi",AliRsnCutMiniPair::kMassRange);
  if(!SideBand) cutMassPhi->SetRangeD(1.01,1.03);
  else cutMassPhi->SetRangeD(1.04,1.06);
  AliRsnCutMiniPair* cutYPhi=new AliRsnCutMiniPair("cutRapidityPhi",AliRsnCutMiniPair::kRapidityRange);
  cutYPhi->SetRangeD(-0.6,0.6);
  AliRsnCutSet* cutsPhi=new AliRsnCutSet("pairCutsPhi",AliRsnTarget::kMother);
  cutsPhi->AddCut(cutMassPhi);
  cutsPhi->AddCut(cutYPhi);
  cutsPhi->SetCutScheme(TString::Format("%s&%s",cutMassPhi->GetName(),cutYPhi->GetName()).Data());

  AliRsnMiniResonanceFinder* rsnfinder=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder",task->GetName()));
  rsnfinder->SetCutID(0,iCutK);
  rsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
  rsnfinder->SetCharge(0,'-');
  rsnfinder->SetCutID(1,iCutK);
  rsnfinder->SetDaughter(1,AliRsnDaughter::kKaon);
  rsnfinder->SetCharge(1,'+');
  rsnfinder->SetResonanceMass(1.01946);
  rsnfinder->SetResonancePDG(333);
  rsnfinder->SetPairCuts(cutsPhi);
  Int_t iCutPhi=task->SetResonanceFinder(rsnfinder);

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  cutY->SetRangeD(-0.5,0.5);
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
  if(!trigger){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.; nmult++;
    multbins[nmult]=5.; nmult++;
    for(j=1;j<=10;j++){multbins[nmult]=j*10; nmult++;}
  }else{
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=0.001; nmult++;
    multbins[nmult]=0.01; nmult++;
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

  Bool_t  use     [4] = {1       , 1       , 1          , 1          };
  Bool_t  useIM   [4] = {1       , 1       , 1          , 1          };
  TString name    [4] = {"PipPhi", "PimPhi", "PipPhiMix", "PimPhiMix"};
  TString comp    [4] = {"PAIR"  , "PAIR"  , "MIX"      , "MIX"      };
  Char_t  charge1 [4] = {'+'     , '-'     , '+'        , '-'        };
  Int_t   ipdg    [4] = {3124    , -3124   , 3124       , -3124      };

  Int_t i;
  AliRsnMiniOutput* out;
  for(i=0;i<4;i++){
    if(!use[i]) continue;
    out=task->CreateOutput(Form("piphi_%s%s",name[i].Data(),suffix),"HIST",comp[i].Data());
    out->SetDaughter(0,AliRsnDaughter::kPion);
    out->SetCutID(0,iCutPi);
    out->SetCharge(0,charge1[i]);

    out->SetDaughter(1,AliRsnDaughter::kPhi);
    out->SetCutID(1,iCutPhi);
    out->SetCharge(1,'0');
    if(!SideBand) out->SetUseStoredMass(1);
    out->SetMotherPDG(ipdg[i]);
    out->SetMotherMass(mass);

    if(i==0 || i==1) out->SetPairCuts(cutsPairSame);
    else out->SetPairCuts(cutsPairMix);

    // axis X: invmass or resolution
    if(useIM[i]) out->AddAxis(imID,170,1.15,2.);
    else out->AddAxis(resID,200,-0.02,0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID,200,0.0,20.0);
    
    // axis Z: centrality-multiplicity
    out->AddAxis(centID,nmult,multbins);
  }

  // fill monitoring histogram for the resonance (phi)

  for(i=0;i<3;i++){
    if(!i){
      name[0].Form("piphi_phimass");
      comp[0].Form("PAIR");
    }else if(i==1){
      if(!isMC) continue;
      name[0].Form("piphi_phimass_gen");
      comp[0].Form("MOTHER");
    }else if(i==2){
      if(!isMC) continue;
      name[0].Form("piphi_phimass_rec");
      comp[0].Form("TRUE");
    }

    out=task->CreateOutput(name[0].Data(),"HIST",comp[0].Data());
    for(j=0;j<2;j++){
      out->SetDaughter(j,rsnfinder->GetDaughter(j));
      out->SetCutID(j,rsnfinder->GetCutID(j));
      out->SetCharge(j,rsnfinder->GetCharge(j));
    }

    out->SetMotherPDG(333);
    out->SetMotherMass(rsnfinder->GetResonanceMass());
    out->SetPairCuts(cutsPhi);
      
    out->AddAxis(imID,70,1.,1.07);
    out->AddAxis(ptID,200,0.0,20.0);
    out->AddAxis(centID,nmult,multbins);
  }

  return kTRUE;
}


//=============================


Bool_t Config_kxphi(
  AliRsnMiniAnalysisTask *task,
  TString     lname,
  Bool_t      isMC,
  Int_t       system,
  Int_t       EventCuts,
  Int_t       TrackCutsK,
  Int_t       TrackCutsPhi
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;

  Double_t mass=0.493677+1.019460;
  
  // set daughter cuts
  if(!(TrackCutsK%10000)) TrackCutsK+=3020;//default settings
  Float_t nsigmaK1TPC=0.1*(TrackCutsK%100);
  Float_t nsigmaK1TOF=0.1*((TrackCutsK/100)%100);
  Int_t CutTypeK1=(TrackCutsK/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only

  if(!(TrackCutsPhi%10000)) TrackCutsPhi+=3020;//default settings
  Float_t nsigmaK2TPC=0.1*(TrackCutsPhi%100);
  Float_t nsigmaK2TOF=0.1*((TrackCutsPhi/100)%100);
  Int_t CutTypeK2=(TrackCutsPhi/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only
  Int_t SideBand=(TrackCutsPhi/100000)%10;

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,
    AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);

  AliRsnCutSetDaughterParticle* cutSetK1=0;
  if(!CutTypeK1) cutSetK1=new AliRsnCutSetDaughterParticle(
    Form("cutK1%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaK1TPC),
    trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kKaon,nsigmaK1TPC,nsigmaK1TOF);
  else if(CutTypeK1==1) cutSetK1=new AliRsnCutSetDaughterParticle(
    Form("cutK1%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,nsigmaK1TPC),
    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kKaon,nsigmaK1TPC,-1.);
  else if(CutTypeK1==2) cutSetK1=new AliRsnCutSetDaughterParticle(
    Form("cutK1%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,nsigmaK1TOF),
    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,AliPID::kKaon,-1.,nsigmaK1TOF);
  if(!cutSetK1){cerr<<"Error in AddTaskResonanceFinder::Config_piphi(): missing cutSetK1"<<endl; return kFALSE;}

  AliRsnCutSetDaughterParticle* cutSetK2=0;
  if(!CutTypeK2) cutSetK2=new AliRsnCutSetDaughterParticle(
    Form("cutK2%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaK2TPC),
    trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kKaon,nsigmaK2TPC,nsigmaK2TOF);
  else if(CutTypeK2==1) cutSetK2=new AliRsnCutSetDaughterParticle(
    Form("cutK2%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,nsigmaK2TPC),
    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kKaon,nsigmaK2TPC,-1.);
  else if(CutTypeK2==2) cutSetK2=new AliRsnCutSetDaughterParticle(
    Form("cutK2%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,nsigmaK2TOF),
    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,AliPID::kKaon,-1.,nsigmaK2TOF);
  if(!cutSetK2){cerr<<"Error in AddTaskResonanceFinder::Config_kxphi(): missing cutSetK2"<<endl; return kFALSE;}

  Int_t iCutQ=task->AddTrackCuts(cutSetQ);
  Int_t iCutK1=task->AddTrackCuts(cutSetK1);
  Int_t iCutK2=task->AddTrackCuts(cutSetK2);

  // monitoring
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
    AddMonitorOutput(isMC,cutSetK1->GetMonitorOutput());
    AddMonitorOutput(isMC,cutSetK2->GetMonitorOutput());
  }
    
  // AliRsnMiniResonanceFinder
  AliRsnCutMiniPair* cutMassPhi=new AliRsnCutMiniPair("cutMassPhi",AliRsnCutMiniPair::kMassRange);
  if(!SideBand) cutMassPhi->SetRangeD(1.01,1.03);
  else cutMassPhi->SetRangeD(1.04,1.06);
  AliRsnCutMiniPair* cutYPhi=new AliRsnCutMiniPair("cutRapidityPhi",AliRsnCutMiniPair::kRapidityRange);
  cutYPhi->SetRangeD(-0.6,0.6);
  AliRsnCutSet* cutsPhi=new AliRsnCutSet("pairCutsPhi",AliRsnTarget::kMother);
  cutsPhi->AddCut(cutMassPhi);
  cutsPhi->AddCut(cutYPhi);
  cutsPhi->SetCutScheme(TString::Format("%s&%s",cutMassPhi->GetName(),cutYPhi->GetName()).Data());

  AliRsnMiniResonanceFinder* rsnfinder=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder",task->GetName()));
  rsnfinder->SetCutID(0,iCutK2);
  rsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
  rsnfinder->SetCharge(0,'-');
  rsnfinder->SetCutID(1,iCutK2);
  rsnfinder->SetDaughter(1,AliRsnDaughter::kKaon);
  rsnfinder->SetCharge(1,'+');
  rsnfinder->SetResonanceMass(1.01946);
  rsnfinder->SetResonancePDG(333);
  rsnfinder->SetPairCuts(cutsPhi);
  Int_t iCutPhi=task->SetResonanceFinder(rsnfinder);

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  cutY->SetRangeD(-0.5,0.5);
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
  if(!trigger){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.; nmult++;
    multbins[nmult]=5.; nmult++;
    for(j=1;j<=10;j++){multbins[nmult]=j*10; nmult++;}
  }else{
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=0.001; nmult++;
    multbins[nmult]=0.01; nmult++;
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

  Bool_t  use     [4] = {1       , 1       , 1          , 1          };
  Bool_t  useIM   [4] = {1       , 1       , 1          , 1          };
  TString name    [4] = {"KpPhi" , "KmPhi" , "KpPhiMix" , "KmPhiMix" };
  TString comp    [4] = {"PAIR"  , "PAIR"  , "MIX"      , "MIX"      };
  Char_t  charge1 [4] = {'+'     , '-'     , '+'        , '-'        };
  Int_t   ipdg    [4] = {3124    , -3124   , 3124       , -3124      };

  Int_t i;
  AliRsnMiniOutput* out;
  for(i=0;i<4;i++){
    if(!use[i]) continue;
    out=task->CreateOutput(Form("kxphi_%s%s",name[i].Data(),suffix),"HIST",comp[i].Data());
    out->SetDaughter(0,AliRsnDaughter::kKaon);
    out->SetCutID(0,iCutK1);
    out->SetCharge(0,charge1[i]);

    out->SetDaughter(1,AliRsnDaughter::kPhi);
    out->SetCutID(1,iCutPhi);
    out->SetCharge(1,'0');
    if(!SideBand) out->SetUseStoredMass(1);
    out->SetMotherPDG(ipdg[i]);
    out->SetMotherMass(mass);

    if(i==0 || i==1) out->SetPairCuts(cutsPairSame);
    else out->SetPairCuts(cutsPairMix);

    // axis X: invmass or resolution
    if(useIM[i]) out->AddAxis(imID,200,1.5,2.5);
    else out->AddAxis(resID,200,-0.02,0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID,200,0.0,20.0);
    
    // axis Z: centrality-multiplicity
    out->AddAxis(centID,nmult,multbins);
  }

  // fill monitoring histogram for the resonance (phi)

  for(i=0;i<3;i++){
    if(!i){
      name[0].Form("kxphi_phimass");
      comp[0].Form("PAIR");
    }else if(i==1){
      if(!isMC) continue;
      name[0].Form("kxphi_phimass_gen");
      comp[0].Form("MOTHER");
    }else if(i==2){
      if(!isMC) continue;
      name[0].Form("kxphi_phimass_rec");
      comp[0].Form("TRUE");
    }

    out=task->CreateOutput(name[0].Data(),"HIST",comp[0].Data());
    for(j=0;j<2;j++){
      out->SetDaughter(j,rsnfinder->GetDaughter(j));
      out->SetCutID(j,rsnfinder->GetCutID(j));
      out->SetCharge(j,rsnfinder->GetCharge(j));
    }

    out->SetMotherPDG(333);
    out->SetMotherMass(rsnfinder->GetResonanceMass());
    out->SetPairCuts(cutsPhi);
      
    out->AddAxis(imID,70,1.,1.07);
    out->AddAxis(ptID,200,0.0,20.0);
    out->AddAxis(centID,nmult,multbins);
  }

  return kTRUE;
}


//=============================


Bool_t Config_k0phi(
  AliRsnMiniAnalysisTask *task,
  TString     lname,
  Bool_t      isMC,
  Int_t       system,
  Int_t       EventCuts,
  Int_t       TrackCutsK0,
  Int_t       TrackCutsPhi
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;

  Double_t mass=0.497611+1.019460;

  if(!(TrackCutsPhi%10000)) TrackCutsPhi+=3020;//default settings
  Float_t nsigmaKTPC=0.1*(TrackCutsPhi%100);
  Float_t nsigmaKTOF=0.1*((TrackCutsPhi/100)%100);
  Int_t CutTypeKx=(TrackCutsPhi/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only
  Int_t SideBand=(TrackCutsPhi/100000)%10;

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,
    AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);
  Int_t iCutQ=task->AddTrackCuts(cutSetQ);
    
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
  Float_t k0sDCA=0.3;
  Float_t k0s_pLife=20.;
  Float_t k0s_radiuslow=0.5;
  Float_t k0s_radiushigh=200.;
  Float_t k0s_massTol=0.02;
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
    
  // selection for phi daughters

  AliRsnCutSetDaughterParticle* cutSetKx=0;
  if(!CutTypeKx) cutSetKx=new AliRsnCutSetDaughterParticle(
    Form("cutKx%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaKTPC),
    trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kKaon,nsigmaKTPC,nsigmaKTOF);
  else if(CutTypeKx==1) cutSetKx=new AliRsnCutSetDaughterParticle(
    Form("cutKx%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,nsigmaKTPC),
    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kKaon,nsigmaKTPC,-1.);
  else if(CutTypeKx==2) cutSetKx=new AliRsnCutSetDaughterParticle(
    Form("cutKx%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,nsigmaKTOF),
    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,AliPID::kKaon,-1.,nsigmaKTOF);
  if(!cutSetKx){cerr<<"Error in AddTaskResonanceFinder::Config_k0phi(): missing cutSetKx"<<endl; return kFALSE;}
  Int_t iCutKx=task->AddTrackCuts(cutSetKx);

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
    
  // AliRsnMiniResonanceFinder
  AliRsnCutMiniPair* cutMassPhi=new AliRsnCutMiniPair("cutMassPhi",AliRsnCutMiniPair::kMassRange);
  if(!SideBand) cutMassPhi->SetRangeD(1.01,1.03);
  else cutMassPhi->SetRangeD(1.04,1.06);
  AliRsnCutMiniPair* cutYPhi=new AliRsnCutMiniPair("cutRapidityPhi",AliRsnCutMiniPair::kRapidityRange);
  cutYPhi->SetRangeD(-0.6,0.6);
  AliRsnCutSet* cutsPhi=new AliRsnCutSet("pairCutsPhi",AliRsnTarget::kMother);
  cutsPhi->AddCut(cutMassPhi);
  cutsPhi->AddCut(cutYPhi);
  cutsPhi->SetCutScheme(TString::Format("%s&%s",cutMassPhi->GetName(),cutYPhi->GetName()).Data());

  AliRsnMiniResonanceFinder* rsnfinder=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder",task->GetName()));
  rsnfinder->SetCutID(0,iCutKx);
  rsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
  rsnfinder->SetCharge(0,'-');
  rsnfinder->SetCutID(1,iCutKx);
  rsnfinder->SetDaughter(1,AliRsnDaughter::kKaon);
  rsnfinder->SetCharge(1,'+');
  rsnfinder->SetResonanceMass(1.01946);
  rsnfinder->SetResonancePDG(333);
  rsnfinder->SetPairCuts(cutsPhi);
  Int_t iCutPhi=task->SetResonanceFinder(rsnfinder);

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  cutY->SetRangeD(-0.5,0.5);
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
  if(!trigger){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.; nmult++;
    multbins[nmult]=5.; nmult++;
    for(j=1;j<=10;j++){multbins[nmult]=j*10; nmult++;}
  }else{
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=0.001; nmult++;
    multbins[nmult]=0.01; nmult++;
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

  Bool_t  use     [2] = {1       , 1        };
  Bool_t  useIM   [2] = {1       , 1        };
  TString name    [2] = {"K0Phi" ,"K0PhiMix"};
  TString comp    [2] = {"PAIR"  ,"MIX"     };
  Int_t   ipdg    [2] = {3124    , -3124    };

  Int_t i;
  AliRsnMiniOutput* out;
  for(i=0;i<2;i++){
    if(!use[i]) continue;
    out=task->CreateOutput(Form("k0phi_%s%s",name[i].Data(),suffix),"HIST",comp[i].Data());
    out->SetDaughter(0,AliRsnDaughter::kKaon0);
    out->SetCutID(0,iCutK0s);
    out->SetCharge(0,'0');

    out->SetDaughter(1,AliRsnDaughter::kPhi);
    out->SetCutID(1,iCutPhi);
    out->SetCharge(1,'0');
    if(!SideBand) out->SetUseStoredMass(1);
    out->SetMotherPDG(ipdg[i]);
    out->SetMotherMass(mass);

    if(i==0 || i==1) out->SetPairCuts(cutsPairSame);
    else out->SetPairCuts(cutsPairMix);

    // axis X: invmass or resolution
    if(useIM[i]) out->AddAxis(imID,200,1.5,2.5);
    else out->AddAxis(resID,200,-0.02,0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID,200,0.0,20.0);
    
    // axis Z: centrality-multiplicity
    out->AddAxis(centID,nmult,multbins);
  }

  // fill monitoring histogram for the resonance (phi)

  for(i=0;i<3;i++){
    if(!i){
      name[0].Form("k0phi_phimass");
      comp[0].Form("PAIR");
    }else if(i==1){
      if(!isMC) continue;
      name[0].Form("k0phi_phimass_gen");
      comp[0].Form("MOTHER");
    }else if(i==2){
      if(!isMC) continue;
      name[0].Form("k0phi_phimass_rec");
      comp[0].Form("TRUE");
    }

    out=task->CreateOutput(name[0].Data(),"HIST",comp[0].Data());
    for(j=0;j<2;j++){
      out->SetDaughter(j,rsnfinder->GetDaughter(j));
      out->SetCutID(j,rsnfinder->GetCutID(j));
      out->SetCharge(j,rsnfinder->GetCharge(j));
    }

    out->SetMotherPDG(333);
    out->SetMotherMass(rsnfinder->GetResonanceMass());
    out->SetPairCuts(cutsPhi);
      
    out->AddAxis(imID,70,1.,1.07);
    out->AddAxis(ptID,200,0.0,20.0);
    out->AddAxis(centID,nmult,multbins);
  }

  return kTRUE;
}


//=============================


Bool_t Config_pphi(
  AliRsnMiniAnalysisTask *task,
  TString     lname,
  Bool_t      isMC,
  Int_t       system,
  Int_t       EventCuts,
  Int_t       TrackCutsP,
  Int_t       TrackCutsPhi
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;

  Double_t mass=0.938272+1.019460;
  
  // set daughter cuts
  if(!(TrackCutsP%10000)) TrackCutsP+=3020;//default settings
  Float_t nsigmaPTPC=0.1*(TrackCutsP%100);
  Float_t nsigmaPTOF=0.1*((TrackCutsP/100)%100);
  Int_t CutTypeP=(TrackCutsP/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only

  if(!(TrackCutsPhi%10000)) TrackCutsPhi+=3020;//default settings
  Float_t nsigmaKTPC=0.1*(TrackCutsPhi%100);
  Float_t nsigmaKTOF=0.1*((TrackCutsPhi/100)%100);
  Int_t CutTypeK=(TrackCutsPhi/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only
  Int_t SideBand=(TrackCutsPhi/100000)%10;

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,
    AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);

  AliRsnCutSetDaughterParticle* cutSetP=0;
  if(!CutTypeP) cutSetP=new AliRsnCutSetDaughterParticle(
    Form("cutProton_%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV,nsigmaPTPC),
    trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV,AliPID::kProton,nsigmaPTPC);
  else if(CutTypeP==1) cutSetP=new AliRsnCutSetDaughterParticle(
    Form("cutProton%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,nsigmaPTPC),
    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kProton,nsigmaPTPC,-1.);
  else if(CutTypeP==2) cutSetP=new AliRsnCutSetDaughterParticle(
    Form("cutProton%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,nsigmaPTOF),
    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,AliPID::kProton,-1.,nsigmaPTOF);
  if(!cutSetP){cerr<<"Error in AddTaskResonanceFinder::Config_pphi(): missing cutSetP"<<endl; return kFALSE;}

  AliRsnCutSetDaughterParticle* cutSetK=0;
  if(!CutTypeK) cutSetK=new AliRsnCutSetDaughterParticle(
    Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaKTPC),
    trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kKaon,nsigmaKTPC,nsigmaKTOF);
  else if(CutTypeK==1) cutSetK=new AliRsnCutSetDaughterParticle(
    Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,nsigmaKTPC),
    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kKaon,nsigmaKTPC,-1.);
  else if(CutTypeK==2) cutSetK=new AliRsnCutSetDaughterParticle(
    Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,nsigmaKTOF),
    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,AliPID::kKaon,-1.,nsigmaKTOF);
  if(!cutSetK){cerr<<"Error in AddTaskResonanceFinder::Config_pphi(): missing cutSetK"<<endl; return kFALSE;}

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
    
  // AliRsnMiniResonanceFinder
  AliRsnCutMiniPair* cutMassPhi=new AliRsnCutMiniPair("cutMassPhi",AliRsnCutMiniPair::kMassRange);
  if(!SideBand) cutMassPhi->SetRangeD(1.01,1.03);
  else cutMassPhi->SetRangeD(1.04,1.06);
  AliRsnCutMiniPair* cutYPhi=new AliRsnCutMiniPair("cutRapidityPhi",AliRsnCutMiniPair::kRapidityRange);
  cutYPhi->SetRangeD(-0.6,0.6);
  AliRsnCutSet* cutsPhi=new AliRsnCutSet("pairCutsPhi",AliRsnTarget::kMother);
  cutsPhi->AddCut(cutMassPhi);
  cutsPhi->AddCut(cutYPhi);
  cutsPhi->SetCutScheme(TString::Format("%s&%s",cutMassPhi->GetName(),cutYPhi->GetName()).Data());

  AliRsnMiniResonanceFinder* rsnfinder=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder",task->GetName()));
  rsnfinder->SetCutID(0,iCutK);
  rsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
  rsnfinder->SetCharge(0,'-');
  rsnfinder->SetCutID(1,iCutK);
  rsnfinder->SetDaughter(1,AliRsnDaughter::kKaon);
  rsnfinder->SetCharge(1,'+');
  rsnfinder->SetResonanceMass(1.01946);
  rsnfinder->SetResonancePDG(333);
  rsnfinder->SetPairCuts(cutsPhi);
  Int_t iCutPhi=task->SetResonanceFinder(rsnfinder);

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  cutY->SetRangeD(-0.5,0.5);
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
  if(!trigger){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.; nmult++;
    multbins[nmult]=5.; nmult++;
    for(j=1;j<=10;j++){multbins[nmult]=j*10; nmult++;}
  }else{
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=0.001; nmult++;
    multbins[nmult]=0.01; nmult++;
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

  Bool_t  use     [4] = {1       ,1       ,1          ,1         };
  Bool_t  useIM   [4] = {1       ,1       ,1          ,1         };
  TString name    [4] = {"PpPhi" ,"PmPhi" ,"PpPhiMix" ,"PmPhiMix"};
  TString comp    [4] = {"PAIR"  ,"PAIR"  ,"MIX"      ,"MIX"     };
  Char_t  charge1 [4] = {'+'     ,'-'     ,'+'        ,'-'       };
  Int_t   ipdg    [4] = {3124    ,-3124   ,3124       ,-3124     };

  Int_t i;
  AliRsnMiniOutput* out;
  for(i=0;i<4;i++){
    if(!use[i]) continue;
    out=task->CreateOutput(Form("pphi_%s%s",name[i].Data(),suffix),"HIST",comp[i].Data());
    out->SetDaughter(0,AliRsnDaughter::kProton);
    out->SetCutID(0,iCutP);
    out->SetCharge(0,charge1[i]);

    out->SetDaughter(1,AliRsnDaughter::kPhi);
    out->SetCutID(1,iCutPhi);
    out->SetCharge(1,'0');
    if(!SideBand) out->SetUseStoredMass(1);
    out->SetMotherPDG(ipdg[i]);
    out->SetMotherMass(mass);

    if(i==0 || i==1) out->SetPairCuts(cutsPairSame);
    else out->SetPairCuts(cutsPairMix);

    // axis X: invmass or resolution
    if(useIM[i]) out->AddAxis(imID,210,1.95,3.);
    else out->AddAxis(resID,200,-0.02,0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID,200,0.0,20.0);
    
    // axis Z: centrality-multiplicity
    out->AddAxis(centID,nmult,multbins);
  }

  // fill monitoring histogram for the resonance (phi)

  for(i=0;i<3;i++){
    if(!i){
      name[0].Form("pphi_phimass");
      comp[0].Form("PAIR");
    }else if(i==1){
      if(!isMC) continue;
      name[0].Form("pphi_phimass_gen");
      comp[0].Form("MOTHER");
    }else if(i==2){
      if(!isMC) continue;
      name[0].Form("pphi_phimass_rec");
      comp[0].Form("TRUE");
    }

    out=task->CreateOutput(name[0].Data(),"HIST",comp[0].Data());
    for(j=0;j<2;j++){
      out->SetDaughter(j,rsnfinder->GetDaughter(j));
      out->SetCutID(j,rsnfinder->GetCutID(j));
      out->SetCharge(j,rsnfinder->GetCharge(j));
    }

    out->SetMotherPDG(333);
    out->SetMotherMass(rsnfinder->GetResonanceMass());
    out->SetPairCuts(cutsPhi);
        
    out->AddAxis(imID,70,1.,1.07);
    out->AddAxis(ptID,200,0.0,20.0);
    out->AddAxis(centID,nmult,multbins);
  }

  return kTRUE;
}


//=============================


Bool_t Config_phiphi(
  AliRsnMiniAnalysisTask *task,
  TString     lname,
  Bool_t      isMC,
  Int_t       system,
  Int_t       EventCuts,
  Int_t       TrackCutsPhi,
  Int_t       TrackCutsDummy
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;

  Double_t mass=2*1.019460;
  
  // set daughter cuts

  if(!(TrackCutsPhi%10000)) TrackCutsPhi+=3020;//default settings
  Float_t nsigmaKTPC=0.1*(TrackCutsPhi%100);
  Float_t nsigmaKTOF=0.1*((TrackCutsPhi/100)%100);
  Int_t CutTypeK=(TrackCutsPhi/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,
    AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);

  AliRsnCutSetDaughterParticle* cutSetK=0;
  if(!CutTypeK) cutSetK=new AliRsnCutSetDaughterParticle(
    Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaKTPC),
    trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kKaon,nsigmaKTPC,nsigmaKTOF);
  else if(CutTypeK==1) cutSetK=new AliRsnCutSetDaughterParticle(
    Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,nsigmaKTPC),
    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kKaon,nsigmaKTPC,-1.);
  else if(CutTypeK==2) cutSetK=new AliRsnCutSetDaughterParticle(
    Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,nsigmaKTOF),
    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,AliPID::kKaon,-1.,nsigmaKTOF);
  if(!cutSetK){cerr<<"Error in AddTaskResonanceFinder::Config_pphi(): missing cutSetK"<<endl; return kFALSE;}

  Int_t iCutQ=task->AddTrackCuts(cutSetQ);
  Int_t iCutK=task->AddTrackCuts(cutSetK);

  // monitoring
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
    AddMonitorOutput(isMC,cutSetK->GetMonitorOutput());
  }
    
  // AliRsnMiniResonanceFinder
  AliRsnCutMiniPair* cutMassPhi=new AliRsnCutMiniPair("cutMassPhi",AliRsnCutMiniPair::kMassRange);
  cutMassPhi->SetRangeD(1.01,1.03);
  AliRsnCutMiniPair* cutYPhi=new AliRsnCutMiniPair("cutRapidityPhi",AliRsnCutMiniPair::kRapidityRange);
  cutYPhi->SetRangeD(-0.6,0.6);
  AliRsnCutSet* cutsPhi=new AliRsnCutSet("pairCutsPhi",AliRsnTarget::kMother);
  cutsPhi->AddCut(cutMassPhi);
  cutsPhi->AddCut(cutYPhi);
  cutsPhi->SetCutScheme(TString::Format("%s&%s",cutMassPhi->GetName(),cutYPhi->GetName()).Data());

  AliRsnMiniResonanceFinder* rsnfinder=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder",task->GetName()));
  rsnfinder->SetCutID(0,iCutK);
  rsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
  rsnfinder->SetCharge(0,'-');
  rsnfinder->SetCutID(1,iCutK);
  rsnfinder->SetDaughter(1,AliRsnDaughter::kKaon);
  rsnfinder->SetCharge(1,'+');
  rsnfinder->SetResonanceMass(1.01946);
  rsnfinder->SetResonancePDG(333);
  rsnfinder->SetPairCuts(cutsPhi);
  Int_t iCutPhi=task->SetResonanceFinder(rsnfinder);

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  cutY->SetRangeD(-0.5,0.5);
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
  if(!trigger){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.; nmult++;
    multbins[nmult]=5.; nmult++;
    for(j=1;j<=10;j++){multbins[nmult]=j*10; nmult++;}
  }else{
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=0.001; nmult++;
    multbins[nmult]=0.01; nmult++;
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

  Bool_t  use     [2] = {1       , 1         };
  Bool_t  useIM   [2] = {1       , 1         };
  TString name    [2] = {"PhiPhi","PhiPhiMix"};
  TString comp    [2] = {"PAIR"  ,"MIX"      };
  Int_t   ipdg    [2] = {3124    , -3124     };

  Int_t i;
  AliRsnMiniOutput* out;
  for(i=0;i<2;i++){
    if(!use[i]) continue;
    out=task->CreateOutput(Form("phihphi_%s%s",name[i].Data(),suffix),"HIST",comp[i].Data());
    out->SetDaughter(0,AliRsnDaughter::kPhi);
    out->SetCutID(0,iCutPhi);
    out->SetCharge(0,'0');
    out->SetUseStoredMass(0);

    out->SetDaughter(1,AliRsnDaughter::kPhi);
    out->SetCutID(1,iCutPhi);
    out->SetCharge(1,'0');
    out->SetUseStoredMass(1);
    out->SetMotherPDG(ipdg[i]);
    out->SetMotherMass(mass);

    if(i==0 || i==1) out->SetPairCuts(cutsPairSame);
    else out->SetPairCuts(cutsPairMix);

    // axis X: invmass or resolution
    if(useIM[i]) out->AddAxis(imID,200,2.,3.);
    else out->AddAxis(resID,200,-0.02,0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID,200,0.0,20.0);
    
    // axis Z: centrality-multiplicity
    out->AddAxis(centID,nmult,multbins);
  }

  // fill monitoring histogram for the resonance (phi)

  for(i=0;i<3;i++){
    if(!i){
      name[0].Form("phiphi_phimass");
      comp[0].Form("PAIR");
    }else if(i==1){
      if(!isMC) continue;
      name[0].Form("phiphi_phimass_gen");
      comp[0].Form("MOTHER");
    }else if(i==2){
      if(!isMC) continue;
      name[0].Form("phiphi_phimass_rec");
      comp[0].Form("TRUE");
    }

    out=task->CreateOutput(name[0].Data(),"HIST",comp[0].Data());
    for(j=0;j<2;j++){
      out->SetDaughter(j,rsnfinder->GetDaughter(j));
      out->SetCutID(j,rsnfinder->GetCutID(j));
      out->SetCharge(j,rsnfinder->GetCharge(j));
    }

    out->SetMotherPDG(333);
    out->SetMotherMass(rsnfinder->GetResonanceMass());
    out->SetPairCuts(cutsPhi);
      
    out->AddAxis(imID,70,1.,1.07);
    out->AddAxis(ptID,200,0.0,20.0);
    out->AddAxis(centID,nmult,multbins);
  }

  return kTRUE;
}


//=============================


Bool_t Config_Lambdaphi(
  AliRsnMiniAnalysisTask *task,
  TString     lname,
  Bool_t      isMC,
  Int_t       system,
  Int_t       EventCuts,
  Int_t       TrackCutsLambda,
  Int_t       TrackCutsPhi
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;

  Double_t mass=1.115683+1.019460;

  if(!(TrackCutsPhi%10000)) TrackCutsPhi+=3020;//default settings
  Float_t nsigmaKTPC=0.1*(TrackCutsPhi%100);
  Float_t nsigmaKTOF=0.1*((TrackCutsPhi/100)%100);
  Int_t CutTypeKx=(TrackCutsPhi/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only
  Int_t SideBand=(TrackCutsPhi/100000)%10;

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,
    AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);
  Int_t iCutQ=task->AddTrackCuts(cutSetQ);
    
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
  Float_t lambdaDaughDCA=1.;//0.5
  Float_t lambdaDCA=0.4;//1.e10 0.3
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
    
  // selection for phi daughters

  AliRsnCutSetDaughterParticle* cutSetKx=0;
  if(!CutTypeKx) cutSetKx=new AliRsnCutSetDaughterParticle(
    Form("cutKx%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaKTPC),
    trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kKaon,nsigmaKTPC,nsigmaKTOF);
  else if(CutTypeKx==1) cutSetKx=new AliRsnCutSetDaughterParticle(
    Form("cutKx%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,nsigmaKTPC),
    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kKaon,nsigmaKTPC,-1.);
  else if(CutTypeKx==2) cutSetKx=new AliRsnCutSetDaughterParticle(
    Form("cutKx%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,nsigmaKTOF),
    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,AliPID::kKaon,-1.,nsigmaKTOF);
  if(!cutSetKx){cerr<<"Error in AddTaskResonanceFinder::Config_k0phi(): missing cutSetKx"<<endl; return kFALSE;}
  Int_t iCutKx=task->AddTrackCuts(cutSetKx);

  // monitoring
  TString pname="lambdap";
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
    AddMonitorOutput(isMC,cutSetKx->GetMonitorOutput());

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
    
  // AliRsnMiniResonanceFinder
  AliRsnCutMiniPair* cutMassPhi=new AliRsnCutMiniPair("cutMassPhi",AliRsnCutMiniPair::kMassRange);
  if(!SideBand) cutMassPhi->SetRangeD(1.01,1.03);
  else cutMassPhi->SetRangeD(1.04,1.06);
  AliRsnCutMiniPair* cutYPhi=new AliRsnCutMiniPair("cutRapidityPhi",AliRsnCutMiniPair::kRapidityRange);
  cutYPhi->SetRangeD(-0.6,0.6);
  AliRsnCutSet* cutsPhi=new AliRsnCutSet("pairCutsPhi",AliRsnTarget::kMother);
  cutsPhi->AddCut(cutMassPhi);
  cutsPhi->AddCut(cutYPhi);
  cutsPhi->SetCutScheme(TString::Format("%s&%s",cutMassPhi->GetName(),cutYPhi->GetName()).Data());

  AliRsnMiniResonanceFinder* rsnfinder=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder",task->GetName()));
  rsnfinder->SetCutID(0,iCutKx);
  rsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
  rsnfinder->SetCharge(0,'-');
  rsnfinder->SetCutID(1,iCutKx);
  rsnfinder->SetDaughter(1,AliRsnDaughter::kKaon);
  rsnfinder->SetCharge(1,'+');
  rsnfinder->SetResonanceMass(1.01946);
  rsnfinder->SetResonancePDG(333);
  rsnfinder->SetPairCuts(cutsPhi);
  Int_t iCutPhi=task->SetResonanceFinder(rsnfinder);

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  cutY->SetRangeD(-0.5,0.5);
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
  if(!trigger){
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=1.; nmult++;
    multbins[nmult]=5.; nmult++;
    for(j=1;j<=10;j++){multbins[nmult]=j*10; nmult++;}
  }else{
    multbins[nmult]=0.; nmult++;
    multbins[nmult]=0.001; nmult++;
    multbins[nmult]=0.01; nmult++;
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

  Bool_t  use     [4] = {1           ,1             ,1              ,1              };
  Bool_t  useIM   [4] = {1           ,1             ,1              ,1              };
  TString name    [4] = {"LambdapPhi","LambdaaPhi"  ,"LambdapPhiMix","LambdaaPhiMix"};
  TString comp    [4] = {"PAIR"      ,"PAIR"        ,"MIX"          ,"MIX"          };
  Int_t cutID1    [4] = {iCutLambda  ,iCutAntiLambda,iCutLambda     ,iCutAntiLambda };
  Int_t   ipdg    [4] = {3124        ,-3124         ,3124           ,-3124          };

  Int_t i;
  AliRsnMiniOutput* out;
  for(i=0;i<4;i++){
    if(!use[i]) continue;
    out=task->CreateOutput(Form("Lambdaphi_%s%s",name[i].Data(),suffix),"HIST",comp[i].Data());
    out->SetDaughter(0,AliRsnDaughter::kLambda);
    out->SetCutID(0,cutID1[i]);
    out->SetCharge(0,'0');

    out->SetDaughter(1,AliRsnDaughter::kPhi);
    out->SetCutID(1,iCutPhi);
    out->SetCharge(1,'0');
    if(!SideBand) out->SetUseStoredMass(1);
    out->SetMotherPDG(ipdg[i]);
    out->SetMotherMass(mass);

    if(i==0 || i==1) out->SetPairCuts(cutsPairSame);
    else out->SetPairCuts(cutsPairMix);

    // axis X: invmass or resolution
    if(useIM[i]) out->AddAxis(imID,280,2.1,3.5);
    else out->AddAxis(resID,200,-0.02,0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID,200,0.0,20.0);
    
    // axis Z: centrality-multiplicity
    out->AddAxis(centID,nmult,multbins);
  }

  // fill monitoring histogram for the resonance (phi)

  for(i=0;i<3;i++){
    if(!i){
      name[0].Form("Lambdaphi_phimass");
      comp[0].Form("PAIR");
    }else if(i==1){
      if(!isMC) continue;
      name[0].Form("Lambdaphi_phimass_gen");
      comp[0].Form("MOTHER");
    }else if(i==2){
      if(!isMC) continue;
      name[0].Form("Lambdaphi_phimass_rec");
      comp[0].Form("TRUE");
    }

    out=task->CreateOutput(name[0].Data(),"HIST",comp[0].Data());
    for(j=0;j<2;j++){
      out->SetDaughter(j,rsnfinder->GetDaughter(j));
      out->SetCutID(j,rsnfinder->GetCutID(j));
      out->SetCharge(j,rsnfinder->GetCharge(j));
    }

    out->SetMotherPDG(333);
    out->SetMotherMass(rsnfinder->GetResonanceMass());
    out->SetPairCuts(cutsPhi);
      
    out->AddAxis(imID,70,1.,1.07);
    out->AddAxis(ptID,200,0.0,20.0);
    out->AddAxis(centID,nmult,multbins);
  }

  return kTRUE;
}


//=============================


void AddMonitorOutput_P(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_mom",name.Data()),AliRsnValueDaughter::kP);
  a->SetBins(0.,10.0,0.05);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}


void AddMonitorOutput_Pt(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_pt",name.Data()),AliRsnValueDaughter::kPt);
  a->SetBins(0.,10.0,0.05);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_Eta(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_eta",name.Data()),AliRsnValueDaughter::kEta);
  a->SetBins(-2.,2.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_DCAxy(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_dcaxy",name.Data()),AliRsnValueDaughter::kDCAXY);
  a->SetBins(-0.5,0.5,0.001);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_DCAz(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_dcaz",name.Data()),AliRsnValueDaughter::kDCAZ);
  a->SetBins(-2.5,2.5,0.005);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_TPCpi(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_TPCpi",name.Data()),AliRsnValueDaughter::kTPCnsigmaPi);
  a->SetBins(-10.,10.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_TPCK(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_TPCK",name.Data()),AliRsnValueDaughter::kTPCnsigmaK);
  a->SetBins(-10.,10.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_TPCp(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_TPCp",name.Data()),AliRsnValueDaughter::kTPCnsigmaP);
  a->SetBins(-10.,10.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_NclTPC(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_NclTPC",name.Data()),AliRsnValueDaughter::kNTPCclusters);
  a->SetBins(-0.5,199.5,1);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_chi2TPC(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_chi2TPC",name.Data()),AliRsnValueDaughter::kTPCchi2);
  a->SetBins(0.0,6,.1);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0NPt(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0npt",name.Data()),AliRsnValueDaughter::kV0NPt);
  a->SetBins(0.,10.0,0.05);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0PPt(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0ppt",name.Data()),AliRsnValueDaughter::kV0PPt);
  a->SetBins(0.,10.0,0.05);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0Mass(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
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

void AddMonitorOutput_V0DCA(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0dca",name.Data()),AliRsnValueDaughter::kV0DCA);
  a->SetBins(0.0,0.4,0.001);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0Radius(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0radius",name.Data()),AliRsnValueDaughter::kV0Radius);
  a->SetBins(0.0,200,0.2);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0Lifetime(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0lifetime",name.Data()),AliRsnValueDaughter::kV0Lifetime);
  a->SetBins(0.0,200,0.1);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0DaughterDCA(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0ddca",name.Data()),AliRsnValueDaughter::kDaughterDCA);
  a->SetBins(0.0,2,0.001);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0DCA2TPV(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  //DCA of secondary tracks to primary vertex
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0dca2tpv",name.Data()),AliRsnValueDaughter::kV0DCAXY);
  a->SetBins(-10.,10.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0CPA(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0cpa",name.Data()),AliRsnValueDaughter::kCosPointAng);
  a->SetBins(0.96,1.,0.001);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0TPCpim(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0TPCpim",name.Data()),AliRsnValueDaughter::kLambdaPionPIDCut);
  a->SetBins(0.,5.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_V0TPCpip(TString name,TObjArray *mon,TString opt,AliRsnLoopDaughter *loop){
  AliRsnValueDaughter* a=new AliRsnValueDaughter(Form("%s_v0TPCpip",name.Data()),AliRsnValueDaughter::kAntiLambdaAntiPionPIDCut);
  a->SetBins(-0.,5.,0.01);
  AliRsnListOutput* o=new AliRsnListOutput(Form("out_%s",a->GetName()),AliRsnListOutput::kHistoDefault);
  o->AddValue(a);
  if (mon) mon->Add(o);
  if (loop) loop->AddOutput(o);
}

void AddMonitorOutput_LambdaProtonPID(TObjArray *mon,TString opt,AliRsnLoopDaughter *lpPID){
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

void AddMonitorOutput_LambdaAntiProtonPID(TObjArray *mon,TString opt,AliRsnLoopDaughter *lapPID){
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
