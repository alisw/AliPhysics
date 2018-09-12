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

Bool_t Config_kxSigmastar(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
Bool_t Config_k0Sigmastar(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
Bool_t Config_kstar0Sigmastar(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
Bool_t Config_kstarxSigmastar(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);


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
  gSystem->Load("libPWGLFresonances.so");
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
  if(d1==AliRsnDaughter::kPion && d2==AliRsnDaughter::kPhi){
    Config_piphi(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kPion && d1==AliRsnDaughter::kPhi){
    Config_piphi(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
        
  }else if(d1==AliRsnDaughter::kKaon && d2==AliRsnDaughter::kPhi){
    Config_kxphi(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kKaon && d1==AliRsnDaughter::kPhi){
    Config_kxphi(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
      
  }else if(d1==AliRsnDaughter::kKaon0 && d2==AliRsnDaughter::kPhi){
    Config_k0phi(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kKaon0 && d1==AliRsnDaughter::kPhi){
    Config_k0phi(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kProton && d2==AliRsnDaughter::kPhi){
    Config_pphi(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kProton && d1==AliRsnDaughter::kPhi){
    Config_pphi(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
      
  }else if(d1==AliRsnDaughter::kPhi && d2==AliRsnDaughter::kPhi){
    Config_phiphi(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);

  }else if(d1==AliRsnDaughter::kLambda && d2==AliRsnDaughter::kPhi){
    Config_Lambdaphi(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kLambda && d1==AliRsnDaughter::kPhi){
    Config_Lambdaphi(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kKaon && (d2==AliRsnDaughter::kSigmastarp || d2==AliRsnDaughter::kSigmastarm)){
    Config_kxSigmastar(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kKaon && (d1==AliRsnDaughter::kSigmastarp || d1==AliRsnDaughter::kSigmastarm)){
    Config_kxSigmastar(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kKaon0 && (d2==AliRsnDaughter::kSigmastarp || d2==AliRsnDaughter::kSigmastarm)){
    Config_k0Sigmastar(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kKaon0 && (d1==AliRsnDaughter::kSigmastarp || d1==AliRsnDaughter::kSigmastarm)){
    Config_k0Sigmastar(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kKstar0 && (d2==AliRsnDaughter::kSigmastarp || d2==AliRsnDaughter::kSigmastarm)){
    Config_kstar0Sigmastar(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kKstar0 && (d1==AliRsnDaughter::kSigmastarp || d1==AliRsnDaughter::kSigmastarm)){
    Config_kstar0Sigmastar(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);

  }else if(d1==AliRsnDaughter::kKstarpm && (d2==AliRsnDaughter::kSigmastarp || d2==AliRsnDaughter::kSigmastarm)){
    Config_kstarxSigmastar(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  }else if(d2==AliRsnDaughter::kKstarpm && (d1==AliRsnDaughter::kSigmastarp || d1==AliRsnDaughter::kSigmastarm)){
    Config_kstarxSigmastar(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
      
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
  if(system==1 || system==2) MultBins=1;

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
  cutMassPhi->SetRangeD(1.01,1.03);
  AliRsnCutMiniPair* cutYPhi=new AliRsnCutMiniPair("cutRapidityPhi",AliRsnCutMiniPair::kRapidityRange);
  cutYPhi->SetRangeD(-0.6,0.6);
  AliRsnCutSet* cutsPhi=new AliRsnCutSet("pairCutsPhi",AliRsnTarget::kMother);
  cutsPhi->AddCut(cutMassPhi);
  cutsPhi->AddCut(cutYPhi);
  cutsPhi->SetCutScheme(TString::Format("%s&%s",cutMassPhi->GetName(),cutYPhi->GetName()).Data());
    
  AliRsnMiniResonanceFinder* finder[4];
  int i,iCutPhi[4];
    
  i=0;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderPhi",task->GetName()));
  finder[i]->SetCutID(0,iCutK);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'-');
  finder[i]->SetCutID(1,iCutK);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsPhi);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);
    
  i=1;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderKpKp",task->GetName()));
  finder[i]->SetCutID(0,iCutK);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'+');
  finder[i]->SetCutID(1,iCutK);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsPhi);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);
    
  i=2;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderKmKm",task->GetName()));
  finder[i]->SetCutID(0,iCutK);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'-');
  finder[i]->SetCutID(1,iCutK);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'-');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsPhi);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);

  AliRsnCutMiniPair* cutMassSB=new AliRsnCutMiniPair("cutMassSB",AliRsnCutMiniPair::kMassRange);
  cutMassSB->SetRangeD(1.04,1.06);
  AliRsnCutSet* cutsSB=new AliRsnCutSet("pairCutsSB",AliRsnTarget::kMother);
  cutsSB->AddCut(cutMassSB);
  cutsSB->AddCut(cutYPhi);
  cutsSB->SetCutScheme(TString::Format("%s&%s",cutMassSB->GetName(),cutYPhi->GetName()).Data());

  i=3;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderSB",task->GetName()));
  finder[i]->SetCutID(0,iCutK);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'-');
  finder[i]->SetCutID(1,iCutK);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsSB);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);

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

  Int_t xID,cut2,pairID,ipdg;
  TString name,comp;
  Char_t charge1;
  AliRsnMiniOutput* out;
    
  for(i=0;i<10;i++){
    if(!i){
      xID=imID;
      name.Form("PipPhi");
      comp.Form("PAIR");
      charge1='+';
      cut2=0;
      pairID=0;
      ipdg=3124;
    }else if(i==1){
      xID=imID;
      name.Form("PimPhi");
      comp.Form("PAIR");
      charge1='-';
      cut2=0;
      pairID=0;
      ipdg=-3124;
    }else if(i==2){
      xID=imID;
      name.Form("PipKpKp");
      comp.Form("PAIR");
      charge1='+';
      cut2=1;
      pairID=0;
      ipdg=3124;
    }else if(i==3){
      xID=imID;
      name.Form("PimKpKp");
      comp.Form("PAIR");
      charge1='-';
      cut2=1;
      pairID=0;
      ipdg=3124;
    }else if(i==4){
      xID=imID;
      name.Form("PipKmKm");
      comp.Form("PAIR");
      charge1='+';
      cut2=2;
      pairID=0;
      ipdg=3124;
    }else if(i==5){
      xID=imID;
      name.Form("PimKmKm");
      comp.Form("PAIR");
      charge1='-';
      cut2=2;
      pairID=0;
      ipdg=3124;
    }else if(i==6){
      xID=imID;
      name.Form("PipSB");
      comp.Form("PAIR");
      charge1='+';
      cut2=3;
      pairID=0;
      ipdg=3124;
    }else if(i==7){
      xID=imID;
      name.Form("PimSB");
      comp.Form("PAIR");
      charge1='-';
      cut2=3;
      pairID=0;
      ipdg=-3124;
    }else if(i==8){
      xID=imID;
      name.Form("PipPhiMix");
      comp.Form("MIX");
      charge1='+';
      cut2=0;
      pairID=1;
      ipdg=3124;
    }else if(i==9){
      xID=imID;
      name.Form("PimPhiMix");
      comp.Form("MIX");
      charge1='-';
      cut2=0;
      pairID=1;
      ipdg=-3124;
    }
        
    out=task->CreateOutput(Form("piphi_%s%s",name.Data(),suffix),"HIST",comp.Data());
    out->SetDaughter(0,AliRsnDaughter::kPion);
    out->SetCutID(0,iCutPi);
    out->SetCharge(0,charge1);

    out->SetDaughter(1,AliRsnDaughter::kPhi);
    out->SetCutID(1,iCutPhi[cut2]);
    out->SetCharge(1,'0');
    if(cut2!=3) out->SetUseStoredMass(1);

    if(!pairID) out->SetPairCuts(cutsPairSame);
    else out->SetPairCuts(cutsPairMix);
    out->SetMotherPDG(ipdg);
    out->SetMotherMass(mass);

    if(xID==imID) out->AddAxis(imID,170,1.15,2.);// axis X: invmass or resolution
    else out->AddAxis(resID,200,-0.02,0.02);
    out->AddAxis(ptID,200,0.0,20.0);// axis Y: transverse momentum
    out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
  }

  // fill monitoring histogram for the resonance (phi)
  for(i=0;i<6;i++){
    if(!i){
      name.Form("phimass");
      comp.Form("PAIR");
      cut2=0;
    }else if(i==1){
      name.Form("KpKpmass");
      comp.Form("PAIR");
      cut2=1;
    }else if(i==2){
      name.Form("KmKmmass");
      comp.Form("PAIR");
      cut2=2;
    }else if(i==3){
      name.Form("SBmass");
      comp.Form("PAIR");
      cut2=3;
    }else if(i==4){
      if(!isMC) continue;
      name.Form("phimass_gen");
      comp.Form("MOTHER");
      cut2=0;
    }else if(i==5){
      if(!isMC) continue;
      name.Form("phimass_rec");
      comp.Form("TRUE");
      cut2=0;
    }

    out=task->CreateOutput(Form("piphi_%s",name.Data()),"HIST",comp.Data());
    out->SetMotherPDG(333);
      
    for(j=0;j<2;j++){
      out->SetDaughter(j,finder[cut2]->GetDaughter(j));
      out->SetCutID(j,finder[cut2]->GetCutID(j));
      out->SetCharge(j,finder[cut2]->GetCharge(j));
    }
    out->SetMotherMass(finder[cut2]->GetResonanceMass());
    if(cut2!=3) out->SetPairCuts(cutsPhi);
    else out->SetPairCuts(cutsSB);
      
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
  if(system==1 || system==2) MultBins=1;

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
  cutMassPhi->SetRangeD(1.01,1.03);
  AliRsnCutMiniPair* cutYPhi=new AliRsnCutMiniPair("cutRapidityPhi",AliRsnCutMiniPair::kRapidityRange);
  cutYPhi->SetRangeD(-0.6,0.6);
  AliRsnCutSet* cutsPhi=new AliRsnCutSet("pairCutsPhi",AliRsnTarget::kMother);
  cutsPhi->AddCut(cutMassPhi);
  cutsPhi->AddCut(cutYPhi);
  cutsPhi->SetCutScheme(TString::Format("%s&%s",cutMassPhi->GetName(),cutYPhi->GetName()).Data());
    
  AliRsnMiniResonanceFinder* finder[4];
  int i,iCutPhi[4];
    
  i=0;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderPhi",task->GetName()));
  finder[i]->SetCutID(0,iCutK2);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'-');
  finder[i]->SetCutID(1,iCutK2);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsPhi);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);
    
  i=1;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderKpKp",task->GetName()));
  finder[i]->SetCutID(0,iCutK2);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'+');
  finder[i]->SetCutID(1,iCutK2);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsPhi);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);
    
  i=2;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderKmKm",task->GetName()));
  finder[i]->SetCutID(0,iCutK2);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'-');
  finder[i]->SetCutID(1,iCutK2);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'-');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsPhi);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);

  AliRsnCutMiniPair* cutMassSB=new AliRsnCutMiniPair("cutMassSB",AliRsnCutMiniPair::kMassRange);
  cutMassSB->SetRangeD(1.04,1.06);
  AliRsnCutSet* cutsSB=new AliRsnCutSet("pairCutsSB",AliRsnTarget::kMother);
  cutsSB->AddCut(cutMassSB);
  cutsSB->AddCut(cutYPhi);
  cutsSB->SetCutScheme(TString::Format("%s&%s",cutMassSB->GetName(),cutYPhi->GetName()).Data());

  i=3;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderSB",task->GetName()));
  finder[i]->SetCutID(0,iCutK2);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'-');
  finder[i]->SetCutID(1,iCutK2);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsSB);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);

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
    
  Int_t xID,cut2,pairID,ipdg;
  TString name,comp;
  Char_t charge1;
  AliRsnMiniOutput* out;
    
  for(i=0;i<6;i++){
    if(!i){
      xID=imID;
      name.Form("KpPhi");
      comp.Form("PAIR");
      charge1='+';
      cut2=0;
      pairID=0;
      ipdg=3124;
    }else if(i==1){
      xID=imID;
      name.Form("KmPhi");
      comp.Form("PAIR");
      charge1='-';
      cut2=0;
      pairID=0;
      ipdg=-3124;
    }else if(i==2){
      xID=imID;
      name.Form("KpKpKp");
      comp.Form("PAIR");
      charge1='+';
      cut2=1;
      pairID=0;
      ipdg=3124;
    }else if(i==3){
      xID=imID;
      name.Form("KmKpKp");
      comp.Form("PAIR");
      charge1='-';
      cut2=1;
      pairID=0;
      ipdg=3124;
    }else if(i==4){
      xID=imID;
      name.Form("KpKmKm");
      comp.Form("PAIR");
      charge1='+';
      cut2=2;
      pairID=0;
      ipdg=3124;
    }else if(i==5){
      xID=imID;
      name.Form("KmKmKm");
      comp.Form("PAIR");
      charge1='-';
      cut2=2;
      pairID=0;
      ipdg=3124;
    }else if(i==6){
      xID=imID;
      name.Form("KpSB");
      comp.Form("PAIR");
      charge1='+';
      cut2=3;
      pairID=0;
      ipdg=3124;
    }else if(i==7){
      xID=imID;
      name.Form("KmSB");
      comp.Form("PAIR");
      charge1='-';
      cut2=3;
      pairID=0;
      ipdg=-3124;
    }else if(i==8){
      xID=imID;
      name.Form("KpPhiMix");
      comp.Form("MIX");
      charge1='+';
      cut2=0;
      pairID=1;
      ipdg=3124;
    }else if(i==9){
      xID=imID;
      name.Form("KmPhiMix");
      comp.Form("MIX");
      charge1='-';
      cut2=0;
      pairID=1;
      ipdg=-3124;
    }
      
    out=task->CreateOutput(Form("kxphi_%s%s",name.Data(),suffix),"HIST",comp.Data());
    out->SetDaughter(0,AliRsnDaughter::kKaon);
    out->SetCutID(0,iCutK1);
    out->SetCharge(0,charge1);

    out->SetDaughter(1,AliRsnDaughter::kPhi);
    out->SetCutID(1,iCutPhi[cut2]);
    out->SetCharge(1,'0');
    if(cut2!=3) out->SetUseStoredMass(1);

    if(!pairID) out->SetPairCuts(cutsPairSame);
    else out->SetPairCuts(cutsPairMix);
    out->SetMotherPDG(ipdg);
    out->SetMotherMass(mass);

    if(xID==imID) out->AddAxis(imID,200,1.5,2.5);// axis X: invmass or resolution
    else out->AddAxis(resID,200,-0.02,0.02);
    out->AddAxis(ptID,200,0.0,20.0);// axis Y: transverse momentum
    out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
  }

  // fill monitoring histogram for the resonance (phi)
  for(i=0;i<6;i++){
    if(!i){
      name.Form("phimass");
      comp.Form("PAIR");
      cut2=0;
    }else if(i==1){
      name.Form("KpKpmass");
      comp.Form("PAIR");
      cut2=1;
    }else if(i==2){
      name.Form("KmKmmass");
      comp.Form("PAIR");
      cut2=2;
    }else if(i==3){
      name.Form("SBmass");
      comp.Form("PAIR");
      cut2=3;
    }else if(i==4){
      if(!isMC) continue;
      name.Form("phimass_gen");
      comp.Form("MOTHER");
      cut2=0;
    }else if(i==5){
      if(!isMC) continue;
      name.Form("phimass_rec");
      comp.Form("TRUE");
      cut2=0;
    }

    out=task->CreateOutput(Form("kxphi_%s",name.Data()),"HIST",comp.Data());
    out->SetMotherPDG(333);
      
    for(j=0;j<2;j++){
      out->SetDaughter(j,finder[cut2]->GetDaughter(j));
      out->SetCutID(j,finder[cut2]->GetCutID(j));
      out->SetCharge(j,finder[cut2]->GetCharge(j));
    }
    out->SetMotherMass(finder[cut2]->GetResonanceMass());
    if(cut2!=3) out->SetPairCuts(cutsPhi);
    else out->SetPairCuts(cutsSB);
      
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
  if(system==1 || system==2) MultBins=1;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;

  Double_t mass=0.497611+1.019460;

  if(!(TrackCutsPhi%10000)) TrackCutsPhi+=3020;//default settings
  Float_t nsigmaKTPC=0.1*(TrackCutsPhi%100);
  Float_t nsigmaKTOF=0.1*((TrackCutsPhi/100)%100);
  Int_t CutTypeKx=(TrackCutsPhi/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only

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
  }  // AliRsnMiniResonanceFinder
  AliRsnCutMiniPair* cutMassPhi=new AliRsnCutMiniPair("cutMassPhi",AliRsnCutMiniPair::kMassRange);
  cutMassPhi->SetRangeD(1.01,1.03);
  AliRsnCutMiniPair* cutYPhi=new AliRsnCutMiniPair("cutRapidityPhi",AliRsnCutMiniPair::kRapidityRange);
  cutYPhi->SetRangeD(-0.6,0.6);
  AliRsnCutSet* cutsPhi=new AliRsnCutSet("pairCutsPhi",AliRsnTarget::kMother);
  cutsPhi->AddCut(cutMassPhi);
  cutsPhi->AddCut(cutYPhi);
  cutsPhi->SetCutScheme(TString::Format("%s&%s",cutMassPhi->GetName(),cutYPhi->GetName()).Data());
    
  AliRsnMiniResonanceFinder* finder[4];
  int i,iCutPhi[4];
    
  i=0;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderPhi",task->GetName()));
  finder[i]->SetCutID(0,iCutKx);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'-');
  finder[i]->SetCutID(1,iCutKx);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsPhi);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);
    
  i=1;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderKpKp",task->GetName()));
  finder[i]->SetCutID(0,iCutKx);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'+');
  finder[i]->SetCutID(1,iCutKx);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsPhi);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);
    
  i=2;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderKmKm",task->GetName()));
  finder[i]->SetCutID(0,iCutKx);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'-');
  finder[i]->SetCutID(1,iCutKx);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'-');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsPhi);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);

  AliRsnCutMiniPair* cutMassSB=new AliRsnCutMiniPair("cutMassSB",AliRsnCutMiniPair::kMassRange);
  cutMassSB->SetRangeD(1.04,1.06);
  AliRsnCutSet* cutsSB=new AliRsnCutSet("pairCutsSB",AliRsnTarget::kMother);
  cutsSB->AddCut(cutMassSB);
  cutsSB->AddCut(cutYPhi);
  cutsSB->SetCutScheme(TString::Format("%s&%s",cutMassSB->GetName(),cutYPhi->GetName()).Data());

  i=3;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderSB",task->GetName()));
  finder[i]->SetCutID(0,iCutKx);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'-');
  finder[i]->SetCutID(1,iCutKx);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsSB);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);

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

  Int_t xID,cut2,pairID,ipdg;
  TString name,comp;
  AliRsnMiniOutput* out;
    
  for(i=0;i<5;i++){
    if(!i){
      xID=imID;
      name.Form("K0Phi");
      comp.Form("PAIR");
      cut2=0;
      pairID=0;
      ipdg=3124;
    }else if(i==1){
      xID=imID;
      name.Form("K0KpKp");
      comp.Form("PAIR");
      cut2=1;
      pairID=0;
      ipdg=3124;
    }else if(i==2){
      xID=imID;
      name.Form("K0KmKm");
      comp.Form("PAIR");
      cut2=2;
      pairID=0;
      ipdg=3124;
    }else if(i==3){
      xID=imID;
      name.Form("K0SB");
      comp.Form("PAIR");
      cut2=3;
      pairID=0;
      ipdg=3124;
    }else if(i==4){
      xID=imID;
      name.Form("K0PhiMix");
      comp.Form("MIX");
      cut2=0;
      pairID=1;
      ipdg=3124;
    }
      
    out=task->CreateOutput(Form("k0phi_%s%s",name.Data(),suffix),"HIST",comp.Data());
    out->SetDaughter(0,AliRsnDaughter::kKaon0);
    out->SetCutID(0,iCutK0s);
    out->SetCharge(0,'0');

    out->SetDaughter(1,AliRsnDaughter::kPhi);
    out->SetCutID(1,iCutPhi[cut2]);
    out->SetCharge(1,'0');
    if(cut2!=3) out->SetUseStoredMass(1);

    if(!pairID) out->SetPairCuts(cutsPairSame);
    else out->SetPairCuts(cutsPairMix);
    out->SetMotherPDG(ipdg);
    out->SetMotherMass(mass);

    if(xID==imID) out->AddAxis(imID,200,1.5,2.5);// axis X: invmass or resolution
    else out->AddAxis(resID,200,-0.02,0.02);
    out->AddAxis(ptID,200,0.0,20.0);// axis Y: transverse momentum
    out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
  }

  // fill monitoring histogram for the resonance (phi)
  
  for(i=0;i<6;i++){
    if(!i){
      name.Form("phimass");
      comp.Form("PAIR");
      cut2=0;
    }else if(i==1){
      name.Form("KpKpmass");
      comp.Form("PAIR");
      cut2=1;
    }else if(i==2){
      name.Form("KmKmmass");
      comp.Form("PAIR");
      cut2=2;
    }else if(i==3){
      name.Form("SBmass");
      comp.Form("PAIR");
      cut2=3;
    }else if(i==4){
      if(!isMC) continue;
      name.Form("phimass_gen");
      comp.Form("MOTHER");
      cut2=0;
    }else if(i==5){
      if(!isMC) continue;
      name.Form("phimass_rec");
      comp.Form("TRUE");
      cut2=0;
    }

    out=task->CreateOutput(Form("k0phi_%s",name.Data()),"HIST",comp.Data());
    out->SetMotherPDG(333);
      
    for(j=0;j<2;j++){
      out->SetDaughter(j,finder[cut2]->GetDaughter(j));
      out->SetCutID(j,finder[cut2]->GetCutID(j));
      out->SetCharge(j,finder[cut2]->GetCharge(j));
    }
    out->SetMotherMass(finder[cut2]->GetResonanceMass());
    if(cut2!=3) out->SetPairCuts(cutsPhi);
    else out->SetPairCuts(cutsSB);
      
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
  if(system==1 || system==2) MultBins=1;

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
  cutMassPhi->SetRangeD(1.01,1.03);
  AliRsnCutMiniPair* cutYPhi=new AliRsnCutMiniPair("cutRapidityPhi",AliRsnCutMiniPair::kRapidityRange);
  cutYPhi->SetRangeD(-0.6,0.6);
  AliRsnCutSet* cutsPhi=new AliRsnCutSet("pairCutsPhi",AliRsnTarget::kMother);
  cutsPhi->AddCut(cutMassPhi);
  cutsPhi->AddCut(cutYPhi);
  cutsPhi->SetCutScheme(TString::Format("%s&%s",cutMassPhi->GetName(),cutYPhi->GetName()).Data());
    
  AliRsnMiniResonanceFinder* finder[4];
  int i,iCutPhi[4];
    
  i=0;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderPhi",task->GetName()));
  finder[i]->SetCutID(0,iCutK);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'-');
  finder[i]->SetCutID(1,iCutK);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsPhi);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);
    
  i=1;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderKpKp",task->GetName()));
  finder[i]->SetCutID(0,iCutK);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'+');
  finder[i]->SetCutID(1,iCutK);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsPhi);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);
    
  i=2;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderKmKm",task->GetName()));
  finder[i]->SetCutID(0,iCutK);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'-');
  finder[i]->SetCutID(1,iCutK);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'-');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsPhi);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);

  AliRsnCutMiniPair* cutMassSB=new AliRsnCutMiniPair("cutMassSB",AliRsnCutMiniPair::kMassRange);
  cutMassSB->SetRangeD(1.04,1.06);
  AliRsnCutSet* cutsSB=new AliRsnCutSet("pairCutsSB",AliRsnTarget::kMother);
  cutsSB->AddCut(cutMassSB);
  cutsSB->AddCut(cutYPhi);
  cutsSB->SetCutScheme(TString::Format("%s&%s",cutMassSB->GetName(),cutYPhi->GetName()).Data());

  i=3;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderSB",task->GetName()));
  finder[i]->SetCutID(0,iCutK);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'-');
  finder[i]->SetCutID(1,iCutK);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsSB);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);

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

  Int_t xID,cut2,pairID,ipdg;
  TString name,comp;
  Char_t charge1;
  AliRsnMiniOutput* out;
    
  for(i=0;i<10;i++){
    if(!i){
      xID=imID;
      name.Form("PpPhi");
      comp.Form("PAIR");
      charge1='+';
      cut2=0;
      pairID=0;
      ipdg=3124;
    }else if(i==1){
      xID=imID;
      name.Form("PmPhi");
      comp.Form("PAIR");
      charge1='-';
      cut2=0;
      pairID=0;
      ipdg=-3124;
    }else if(i==2){
      xID=imID;
      name.Form("PpKpKp");
      comp.Form("PAIR");
      charge1='+';
      cut2=1;
      pairID=0;
      ipdg=3124;
    }else if(i==3){
      xID=imID;
      name.Form("PmKpKp");
      comp.Form("PAIR");
      charge1='-';
      cut2=1;
      pairID=0;
      ipdg=3124;
    }else if(i==4){
      xID=imID;
      name.Form("PpKmKm");
      comp.Form("PAIR");
      charge1='+';
      cut2=2;
      pairID=0;
      ipdg=3124;
    }else if(i==5){
      xID=imID;
      name.Form("PmKmKm");
      comp.Form("PAIR");
      charge1='-';
      cut2=2;
      pairID=0;
      ipdg=3124;
    }else if(i==6){
      xID=imID;
      name.Form("PpSB");
      comp.Form("PAIR");
      charge1='+';
      cut2=3;
      pairID=0;
      ipdg=3124;
    }else if(i==7){
      xID=imID;
      name.Form("PmSB");
      comp.Form("PAIR");
      charge1='-';
      cut2=3;
      pairID=0;
      ipdg=-3124;
    }else if(i==8){
      xID=imID;
      name.Form("PpPhiMix");
      comp.Form("MIX");
      charge1='+';
      cut2=0;
      pairID=1;
      ipdg=3124;
    }else if(i==9){
      xID=imID;
      name.Form("PmPhiMix");
      comp.Form("MIX");
      charge1='-';
      cut2=0;
      pairID=1;
      ipdg=-3124;
    }
      
    out=task->CreateOutput(Form("pphi_%s%s",name.Data(),suffix),"HIST",comp.Data());
    out->SetDaughter(0,AliRsnDaughter::kProton);
    out->SetCutID(0,iCutP);
    out->SetCharge(0,charge1);

    out->SetDaughter(1,AliRsnDaughter::kPhi);
    out->SetCutID(1,iCutPhi[cut2]);
    out->SetCharge(1,'0');
    if(cut2!=3) out->SetUseStoredMass(1);

    if(!pairID) out->SetPairCuts(cutsPairSame);
    else out->SetPairCuts(cutsPairMix);
    out->SetMotherPDG(ipdg);
    out->SetMotherMass(mass);

    if(xID==imID) out->AddAxis(imID,210,1.95,3.);// axis X: invmass or resolution
    else out->AddAxis(resID,200,-0.02,0.02);
    out->AddAxis(ptID,200,0.0,20.0);// axis Y: transverse momentum
    out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
  }

  // fill monitoring histogram for the resonance (phi)
  for(i=0;i<6;i++){
    if(!i){
      name.Form("phimass");
      comp.Form("PAIR");
      cut2=0;
    }else if(i==1){
      name.Form("KpKpmass");
      comp.Form("PAIR");
      cut2=1;
    }else if(i==2){
      name.Form("KmKmmass");
      comp.Form("PAIR");
      cut2=2;
    }else if(i==3){
      name.Form("SBmass");
      comp.Form("PAIR");
      cut2=3;
    }else if(i==4){
      if(!isMC) continue;
      name.Form("phimass_gen");
      comp.Form("MOTHER");
      cut2=0;
    }else if(i==5){
      if(!isMC) continue;
      name.Form("phimass_rec");
      comp.Form("TRUE");
      cut2=0;
    }

    out=task->CreateOutput(Form("pphi_%s",name.Data()),"HIST",comp.Data());
    out->SetMotherPDG(333);
      
    for(j=0;j<2;j++){
      out->SetDaughter(j,finder[cut2]->GetDaughter(j));
      out->SetCutID(j,finder[cut2]->GetCutID(j));
      out->SetCharge(j,finder[cut2]->GetCharge(j));
    }
    out->SetMotherMass(finder[cut2]->GetResonanceMass());
    if(cut2!=3) out->SetPairCuts(cutsPhi);
    else out->SetPairCuts(cutsSB);
      
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
  if(system==1 || system==2) MultBins=1;

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
    
  AliRsnMiniResonanceFinder* finder[4];
  int i,iCutPhi[4];
    
  i=0;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderPhi",task->GetName()));
  finder[i]->SetCutID(0,iCutK);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'-');
  finder[i]->SetCutID(1,iCutK);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsPhi);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);
    
  i=1;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderKpKp",task->GetName()));
  finder[i]->SetCutID(0,iCutK);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'+');
  finder[i]->SetCutID(1,iCutK);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsPhi);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);
    
  i=2;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderKmKm",task->GetName()));
  finder[i]->SetCutID(0,iCutK);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'-');
  finder[i]->SetCutID(1,iCutK);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'-');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsPhi);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);

  AliRsnCutMiniPair* cutMassSB=new AliRsnCutMiniPair("cutMassSB",AliRsnCutMiniPair::kMassRange);
  cutMassSB->SetRangeD(1.04,1.06);
  AliRsnCutSet* cutsSB=new AliRsnCutSet("pairCutsSB",AliRsnTarget::kMother);
  cutsSB->AddCut(cutMassSB);
  cutsSB->AddCut(cutYPhi);
  cutsSB->SetCutScheme(TString::Format("%s&%s",cutMassSB->GetName(),cutYPhi->GetName()).Data());

  i=3;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderSB",task->GetName()));
  finder[i]->SetCutID(0,iCutK);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'-');
  finder[i]->SetCutID(1,iCutK);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsSB);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);

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

  Int_t xID,cut2,pairID,ipdg;
  TString name,comp;
  AliRsnMiniOutput* out;
    
  for(i=0;i<5;i++){
    if(!i){
      xID=imID;
      name.Form("PhiPhi");
      comp.Form("PAIR");
      cut2=0;
      pairID=0;
      ipdg=3124;
    }else if(i==1){
      xID=imID;
      name.Form("PhiKpKp");
      comp.Form("PAIR");
      cut2=1;
      pairID=0;
      ipdg=3124;
    }else if(i==2){
      xID=imID;
      name.Form("PhiKmKm");
      comp.Form("PAIR");
      cut2=2;
      pairID=0;
      ipdg=3124;
    }else if(i==3){
      xID=imID;
      name.Form("PhiSB");
      comp.Form("PAIR");
      cut2=3;
      pairID=0;
      ipdg=3124;
    }else if(i==4){
      xID=imID;
      name.Form("PhiPhiMix");
      comp.Form("MIX");
      cut2=0;
      pairID=1;
      ipdg=3124;
    }
      
    out=task->CreateOutput(Form("phiphi_%s%s",name.Data(),suffix),"HIST",comp.Data());
    out->SetDaughter(0,AliRsnDaughter::kPhi);
    out->SetCutID(0,iCutPhi[0]);
    out->SetCharge(0,'0');
    out->SetUseStoredMass(0);

    out->SetDaughter(1,AliRsnDaughter::kPhi);
    out->SetCutID(1,iCutPhi[cut2]);
    out->SetCharge(1,'0');
    if(cut2!=3) out->SetUseStoredMass(1);

    if(!pairID) out->SetPairCuts(cutsPairSame);
    else out->SetPairCuts(cutsPairMix);
    out->SetMotherPDG(ipdg);
    out->SetMotherMass(mass);

    if(xID==imID) out->AddAxis(imID,200,2.,3.);// axis X: invmass or resolution
    else out->AddAxis(resID,200,-0.02,0.02);
    out->AddAxis(ptID,200,0.0,20.0);// axis Y: transverse momentum
    out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
  }

  // fill monitoring histogram for the resonance (phi)
  for(i=0;i<6;i++){
    if(!i){
      name.Form("phimass");
      comp.Form("PAIR");
      cut2=0;
    }else if(i==1){
      name.Form("KpKpmass");
      comp.Form("PAIR");
      cut2=1;
    }else if(i==2){
      name.Form("KmKmmass");
      comp.Form("PAIR");
      cut2=2;
    }else if(i==3){
      name.Form("SBmass");
      comp.Form("PAIR");
      cut2=3;
    }else if(i==4){
      if(!isMC) continue;
      name.Form("phimass_gen");
      comp.Form("MOTHER");
      cut2=0;
    }else if(i==5){
      if(!isMC) continue;
      name.Form("phimass_rec");
      comp.Form("TRUE");
      cut2=0;
    }

    out=task->CreateOutput(Form("phiphi_%s",name.Data()),"HIST",comp.Data());
    out->SetMotherPDG(333);
      
    for(j=0;j<2;j++){
      out->SetDaughter(j,finder[cut2]->GetDaughter(j));
      out->SetCutID(j,finder[cut2]->GetCutID(j));
      out->SetCharge(j,finder[cut2]->GetCharge(j));
    }
    out->SetMotherMass(finder[cut2]->GetResonanceMass());
    if(cut2!=3) out->SetPairCuts(cutsPhi);
    else out->SetPairCuts(cutsSB);
      
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
  if(system==1 || system==2) MultBins=1;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;

  Double_t mass=1.115683+1.019460;

  if(!(TrackCutsPhi%10000)) TrackCutsPhi+=3020;//default settings
  Float_t nsigmaKTPC=0.1*(TrackCutsPhi%100);
  Float_t nsigmaKTOF=0.1*((TrackCutsPhi/100)%100);
  Int_t CutTypeKx=(TrackCutsPhi/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only

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
  if(!cutSetKx){cerr<<"Error in AddTaskResonanceFinder::Config_Lambdaphi(): missing cutSetKx"<<endl; return kFALSE;}
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
  cutMassPhi->SetRangeD(1.01,1.03);
  AliRsnCutMiniPair* cutYPhi=new AliRsnCutMiniPair("cutRapidityPhi",AliRsnCutMiniPair::kRapidityRange);
  cutYPhi->SetRangeD(-0.6,0.6);
  AliRsnCutSet* cutsPhi=new AliRsnCutSet("pairCutsPhi",AliRsnTarget::kMother);
  cutsPhi->AddCut(cutMassPhi);
  cutsPhi->AddCut(cutYPhi);
  cutsPhi->SetCutScheme(TString::Format("%s&%s",cutMassPhi->GetName(),cutYPhi->GetName()).Data());
    
  AliRsnMiniResonanceFinder* finder[4];
  int i,iCutPhi[4];
    
  i=0;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderPhi",task->GetName()));
  finder[i]->SetCutID(0,iCutKx);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'-');
  finder[i]->SetCutID(1,iCutKx);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsPhi);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);
    
  i=1;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderKpKp",task->GetName()));
  finder[i]->SetCutID(0,iCutKx);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'+');
  finder[i]->SetCutID(1,iCutKx);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsPhi);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);
    
  i=2;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderKmKm",task->GetName()));
  finder[i]->SetCutID(0,iCutKx);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'-');
  finder[i]->SetCutID(1,iCutKx);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'-');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsPhi);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);

  AliRsnCutMiniPair* cutMassSB=new AliRsnCutMiniPair("cutMassSB",AliRsnCutMiniPair::kMassRange);
  cutMassSB->SetRangeD(1.04,1.06);
  AliRsnCutSet* cutsSB=new AliRsnCutSet("pairCutsSB",AliRsnTarget::kMother);
  cutsSB->AddCut(cutMassSB);
  cutsSB->AddCut(cutYPhi);
  cutsSB->SetCutScheme(TString::Format("%s&%s",cutMassSB->GetName(),cutYPhi->GetName()).Data());

  i=3;
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderSB",task->GetName()));
  finder[i]->SetCutID(0,iCutKx);
  finder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(0,'-');
  finder[i]->SetCutID(1,iCutKx);
  finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.01946);
  finder[i]->SetResonancePDG(333);
  finder[i]->SetPairCuts(cutsSB);
  iCutPhi[i]=task->AddResonanceFinder(finder[i]);

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

  Int_t xID,cut1,cut2,pairID,ipdg;
  TString name,comp;
  AliRsnMiniOutput* out;
    
  for(i=0;i<10;i++){
    if(!i){
      xID=imID;
      name.Form("LambdapPhi");
      comp.Form("PAIR");
      cut1=iCutLambda;
      cut2=0;
      pairID=0;
      ipdg=3124;
    }else if(i==1){
      xID=imID;
      name.Form("LambdaaPhi");
      comp.Form("PAIR");
      cut1=iCutAntiLambda;
      cut2=0;
      pairID=0;
      ipdg=-3124;
    }else if(i==2){
      xID=imID;
      name.Form("LambdapKpKp");
      comp.Form("PAIR");
      cut1=iCutLambda;
      cut2=1;
      pairID=0;
      ipdg=3124;
    }else if(i==3){
      xID=imID;
      name.Form("LambdaaKpKp");
      comp.Form("PAIR");
      cut1=iCutAntiLambda;
      cut2=1;
      pairID=0;
      ipdg=-3124;
    }else if(i==4){
      xID=imID;
      name.Form("LambdapKmKm");
      comp.Form("PAIR");
      cut1=iCutLambda;
      cut2=2;
      pairID=0;
      ipdg=3124;
    }else if(i==5){
      xID=imID;
      name.Form("LambdaaKmKm");
      comp.Form("PAIR");
      cut1=iCutAntiLambda;
      cut2=2;
      pairID=0;
      ipdg=-3124;
    }else if(i==6){
      xID=imID;
      name.Form("LambdapSB");
      comp.Form("PAIR");
      cut1=iCutLambda;
      cut2=3;
      pairID=0;
      ipdg=3124;
    }else if(i==7){
      xID=imID;
      name.Form("LambdaaSB");
      comp.Form("PAIR");
      cut1=iCutAntiLambda;
      cut2=3;
      pairID=0;
      ipdg=-3124;
    }else if(i==8){
      xID=imID;
      name.Form("LambdapPhiMix");
      comp.Form("MIX");
      cut1=iCutLambda;
      cut2=0;
      pairID=1;
      ipdg=3124;
    }else if(i==9){
      xID=imID;
      name.Form("LambdaaPhiMix");
      comp.Form("MIX");
      cut1=iCutAntiLambda;
      cut2=0;
      pairID=1;
      ipdg=-3124;
    }
      
    out=task->CreateOutput(Form("Lambdaphi_%s%s",name.Data(),suffix),"HIST",comp.Data());
    out->SetDaughter(0,AliRsnDaughter::kLambda);
    out->SetCutID(0,cut1);
    out->SetCharge(0,'0');

    out->SetDaughter(1,AliRsnDaughter::kPhi);
    out->SetCutID(1,iCutPhi[cut2]);
    out->SetCharge(1,'0');
    if(cut2!=3) out->SetUseStoredMass(1);

    if(!pairID) out->SetPairCuts(cutsPairSame);
    else out->SetPairCuts(cutsPairMix);
    out->SetMotherPDG(ipdg);
    out->SetMotherMass(mass);

    if(xID==imID) out->AddAxis(imID,280,2.1,3.5);// axis X: invmass or resolution
    else out->AddAxis(resID,200,-0.02,0.02);
    out->AddAxis(ptID,200,0.0,20.0);// axis Y: transverse momentum
    out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
  }

  // fill monitoring histogram for the resonance (phi)
  for(i=0;i<6;i++){
    if(!i){
      name.Form("phimass");
      comp.Form("PAIR");
      cut2=0;
    }else if(i==1){
      name.Form("KpKpmass");
      comp.Form("PAIR");
      cut2=1;
    }else if(i==2){
      name.Form("KmKmmass");
      comp.Form("PAIR");
      cut2=2;
    }else if(i==3){
      name.Form("SBmass");
      comp.Form("PAIR");
      cut2=3;
    }else if(i==4){
      if(!isMC) continue;
      name.Form("phimass_gen");
      comp.Form("MOTHER");
      cut2=0;
    }else if(i==5){
      if(!isMC) continue;
      name.Form("phimass_rec");
      comp.Form("TRUE");
      cut2=0;
    }

    out=task->CreateOutput(Form("Lambdaphi_%s",name.Data()),"HIST",comp.Data());
    out->SetMotherPDG(333);
      
    for(j=0;j<2;j++){
      out->SetDaughter(j,finder[cut2]->GetDaughter(j));
      out->SetCutID(j,finder[cut2]->GetCutID(j));
      out->SetCharge(j,finder[cut2]->GetCharge(j));
    }
    out->SetMotherMass(finder[cut2]->GetResonanceMass());
    if(cut2!=3) out->SetPairCuts(cutsPhi);
    else out->SetPairCuts(cutsSB);
      
    out->AddAxis(imID,70,1.,1.07);
    out->AddAxis(ptID,200,0.0,20.0);
    out->AddAxis(centID,nmult,multbins);
  }

  return kTRUE;
}


//=============================


Bool_t Config_kxSigmastar(
  AliRsnMiniAnalysisTask *task,
  TString     lname,
  Bool_t      isMC,
  Int_t       system,
  Int_t       EventCuts,
  Int_t       TrackCutsK,
  Int_t       TrackCutsS
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;
  if(system==1 || system==2) MultBins=1;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;

  Double_t mass=0.493677+1.385;

  // set cuts for primary bachelor kaon
  if(!(TrackCutsK%10000)) TrackCutsK+=3020;//default settings
  Float_t nsigmaKTPC=0.1*(TrackCutsK%100);
  Float_t nsigmaKTOF=0.1*((TrackCutsK/100)%100);

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);
  AliRsnCutSetDaughterParticle* cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaKTPC),
                                                                         trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kKaon,nsigmaKTPC,nsigmaKTOF);

  Int_t iCutQ=task->AddTrackCuts(cutSetQ);
  Int_t iCutK=task->AddTrackCuts(cutSetK);

  // set cuts for primary pion from Sigma*
  Int_t TrackCutsPi=TrackCutsS/10000;
  if(!(TrackCutsPi%10000)) TrackCutsPi+=3020;//default settings
  Float_t nsigmaPiTPC=0.1*(TrackCutsPi%100);
  Float_t nsigmaPiTOF=0.1*((TrackCutsPi/100)%100);
  Int_t CutTypePi=(TrackCutsPi/10000)%100;//0=TPC+TOFveto (default), 1=TPC only, 2=TOF only, 3 TPC+TOFcut
  Int_t MisidentifiedAsKaon=(TrackCutsPi/1000000)%10;//0=pion assigned pion mass, 1=pion assigned kaon mass

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

  Int_t iCutPi=task->AddTrackCuts(cutSetPi);

  // set cuts for Lambda from Sigma*

  Int_t V0Cuts=TrackCutsS%10000;

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

  // AliRsnMiniResonanceFinder

  AliRsnMiniResonanceFinder* finder[8];
  Int_t i,iCutSig[8];

  AliRsnCutMiniPair* cutMassS=new AliRsnCutMiniPair("cutMassSigmastar",AliRsnCutMiniPair::kMassRange);
  cutMassS->SetRangeD(1.308,1.466);
  AliRsnCutMiniPair* cutYS=new AliRsnCutMiniPair("cutRapiditySigmastar",AliRsnCutMiniPair::kRapidityRange);
  cutYS->SetRangeD(-0.6,0.6);
  AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
  AliRsnCutSet* cutsS=new AliRsnCutSet("pairCutsSigmastar",AliRsnTarget::kMother);
  cutsS->AddCut(cutMassS);
  cutsS->AddCut(cutYS);
  cutsS->AddCut(cutV0);
  cutsS->SetCutScheme(TString::Format("%s&%s&(!%s)",cutMassS->GetName(),cutYS->GetName(),cutV0->GetName()).Data());

  i=0; // Sigma*+
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_Sigmastarp",task->GetName()));
  finder[i]->SetCutID(0,iCutLambda);
  finder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  finder[i]->SetCharge(0,'0');
  finder[i]->SetCutID(1,iCutPi);
  finder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.3828);
  finder[i]->SetResonancePDG(3224);
  finder[i]->SetPairCuts(cutsS);
  iCutSig[i]=task->AddResonanceFinder(finder[i]);

  i=1; // anti-Sigma*-
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_AntiSigmastarm",task->GetName()));
  finder[i]->SetCutID(0,iCutAntiLambda);
  finder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  finder[i]->SetCharge(0,'0');
  finder[i]->SetCutID(1,iCutPi);
  finder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  finder[i]->SetCharge(1,'-');
  finder[i]->SetResonanceMass(1.3828);
  finder[i]->SetResonancePDG(-3224);
  finder[i]->SetPairCuts(cutsS);
  iCutSig[i]=task->AddResonanceFinder(finder[i]);

  i=2; // Sigma*-
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_Sigmastarm",task->GetName()));
  finder[i]->SetCutID(0,iCutLambda);
  finder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  finder[i]->SetCharge(0,'0');
  finder[i]->SetCutID(1,iCutPi);
  finder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  finder[i]->SetCharge(1,'-');
  finder[i]->SetResonanceMass(1.3872);
  finder[i]->SetResonancePDG(3114);
  finder[i]->SetPairCuts(cutsS);
  iCutSig[i]=task->AddResonanceFinder(finder[i]);

  i=3; // anti-Sigma*+
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_AntiSigmastarp",task->GetName()));
  finder[i]->SetCutID(0,iCutAntiLambda);
  finder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  finder[i]->SetCharge(0,'0');
  finder[i]->SetCutID(1,iCutPi);
  finder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.3872);
  finder[i]->SetResonancePDG(-3114);
  finder[i]->SetPairCuts(cutsS);
  iCutSig[i]=task->AddResonanceFinder(finder[i]);
    
  // sidebands
  AliRsnCutMiniPair* cutMassSB=new AliRsnCutMiniPair("cutMassSB",AliRsnCutMiniPair::kMassRange);
  cutMassSB->SetRangeD(1.505,1.585);
  AliRsnCutSet* cutsSB=new AliRsnCutSet("pairCutsSB",AliRsnTarget::kMother);
  cutsSB->AddCut(cutMassSB);
  cutsSB->AddCut(cutYS);
  cutsSB->AddCut(cutV0);
  cutsSB->SetCutScheme(TString::Format("%s&%s&(!%s)",cutMassSB->GetName(),cutYS->GetName(),cutV0->GetName()).Data());

  i=4; // Sigma*+ sideband
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBpp",task->GetName()));
  finder[i]->SetCutID(0,iCutLambda);
  finder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  finder[i]->SetCharge(0,'0');
  finder[i]->SetCutID(1,iCutPi);
  finder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.3828);
  finder[i]->SetResonancePDG(3224);
  finder[i]->SetPairCuts(cutsSB);
  iCutSig[i]=task->AddResonanceFinder(finder[i]);

  i=5; // anti-Sigma*- sideband
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBma",task->GetName()));
  finder[i]->SetCutID(0,iCutAntiLambda);
  finder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  finder[i]->SetCharge(0,'0');
  finder[i]->SetCutID(1,iCutPi);
  finder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  finder[i]->SetCharge(1,'-');
  finder[i]->SetResonanceMass(1.3828);
  finder[i]->SetResonancePDG(-3224);
  finder[i]->SetPairCuts(cutsSB);
  iCutSig[i]=task->AddResonanceFinder(finder[i]);

  i=6; // Sigma*- sideband
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBmp",task->GetName()));
  finder[i]->SetCutID(0,iCutLambda);
  finder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  finder[i]->SetCharge(0,'0');
  finder[i]->SetCutID(1,iCutPi);
  finder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  finder[i]->SetCharge(1,'-');
  finder[i]->SetResonanceMass(1.3872);
  finder[i]->SetResonancePDG(3114);
  finder[i]->SetPairCuts(cutsSB);
  iCutSig[i]=task->AddResonanceFinder(finder[i]);

  i=7; // anti-Sigma*+ sideband
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBpa",task->GetName()));
  finder[i]->SetCutID(0,iCutAntiLambda);
  finder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  finder[i]->SetCharge(0,'0');
  finder[i]->SetCutID(1,iCutPi);
  finder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.3872);
  finder[i]->SetResonancePDG(-3114);
  finder[i]->SetPairCuts(cutsSB);
  iCutSig[i]=task->AddResonanceFinder(finder[i]);

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  cutY->SetRangeD(-0.5,0.5);
   
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
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass,    kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);

  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges

  Int_t k,xID,cut2,pairID,ipdg;
  AliRsnDaughter::ESpecies d2=AliRsnDaughter::kUnknown;
  TString name,comp;
  Char_t charge1,charge2;
  AliRsnMiniOutput* out;

  for(i=0;i<2;i++) for(j=0;j<4;j++) for(k=0;k<3;k++){
    if(!i){
      name.Form("Kp");
      charge1='+';
    }else{
      name.Form("Km");
      charge1='-';
    }

    if(!j){
      name.Append("Sigmastarpp");
      d2=AliRsnDaughter::kSigmastarp;
      charge2='+';
      cut2=0;
    }else if(j==1){
      name.Append("Sigmastarma");
      d2=AliRsnDaughter::kSigmastarp;
      charge2='-';
      cut2=1;
    }else if(j==2){
      name.Append("Sigmastarmp");
      d2=AliRsnDaughter::kSigmastarm;
      charge2='-';
      cut2=2;
    }else{
      name.Append("Sigmastarpa");
      d2=AliRsnDaughter::kSigmastarm;
      charge2='+';
      cut2=3;
    }

    if(!k){
      comp.Form("PAIR");
    }else if(k==1){
      name.Append("SB");
      comp.Form("PAIR");
      cut2+=4;
    }else if(k==2){
      name.Append("Mix");
      comp.Form("MIX");
    }

    out=task->CreateOutput(Form("kxSigmastar_%s%s",name.Data(),suffix),"HIST",comp.Data());
    out->SetDaughter(0,AliRsnDaughter::kKaon);
    out->SetCutID(0,iCutK);
    out->SetCharge(0,charge1);

    out->SetDaughter(1,d2);
    out->SetCutID(1,iCutSig[cut2]);
    out->SetCharge(1,charge2);
    if(k!=1) out->SetUseStoredMass(1);

    if(k<=1) out->SetPairCuts(cutsPairSame);
    else out->SetPairCuts(cutsPairMix);
    out->SetMotherPDG(3124);
    out->SetMotherMass(mass);

    out->AddAxis(imID,260,1.7,3);// axis X: invmass or resolution
    //out->AddAxis(resID,200,-0.02,0.02);
    out->AddAxis(ptID,50,0.0,20.0);// axis Y: transverse momentum
    out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
  }

  // fill monitoring histogram for the resonance (Sigma*)
  for(i=0;i<4;i++) for(j=0;j<4;j++){
    if(!i) name.Form("Sigmastarpp");
    else if(i==1) name.Form("Sigmastarma");
    else if(i==2) name.Form("Sigmastarmp");
    else if(i==3) name.Form("Sigmastarpa");

    if(!j){
      name.Append("_mass");
      comp.Form("PAIR");
    }else if(j==1){
      name.Append("_SBmass");
      comp.Form("PAIR");
    }else if(j==2){
      if(!isMC) continue;
      name.Append("_genmass");
      comp.Form("MOTHER");
    }else if(j==3){
      if(!isMC) continue;
      name.Append("_recmass");
      comp.Form("TRUE");
    }

    k=i;
    if(j==1) k+=4;

    out=task->CreateOutput(Form("kxSigmastar_%s",name.Data()),"HIST",comp.Data());
    out->SetMotherPDG(finder[k]->GetResonancePDG());

    out->SetDaughter(0,finder[k]->GetDaughter(0));
    out->SetCutID(0,finder[k]->GetCutID(0));
    out->SetCharge(0,finder[k]->GetCharge(0));

    out->SetDaughter(1,finder[k]->GetDaughter(1));
    out->SetCutID(1,finder[k]->GetCutID(1));
    out->SetCharge(1,finder[k]->GetCharge(1));

    out->SetMotherMass(finder[k]->GetResonanceMass());
    if(j!=1) out->SetPairCuts(cutsS);
    else out->SetPairCuts(cutsSB);

    out->AddAxis(imID,150,1.3,1.6);
    out->AddAxis(ptID,50,0.0,20.0);
    out->AddAxis(centID,nmult,multbins);
  }

  return kTRUE;
}


//=============================


Bool_t Config_k0Sigmastar(
  AliRsnMiniAnalysisTask *task,
  TString     lname,
  Bool_t      isMC,
  Int_t       system,
  Int_t       EventCuts,
  Int_t       TrackCutsK,
  Int_t       TrackCutsS
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;
  if(system==1 || system==2) MultBins=1;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;

  Double_t mass=0.497611+1.385;

  Int_t V0Cuts=TrackCutsK;

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

  // set cuts for primary pion from Sigma*
  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
  AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);
  Int_t iCutQ=task->AddTrackCuts(cutSetQ);

  Int_t TrackCutsPi=TrackCutsS/10000;
  if(!(TrackCutsPi%10000)) TrackCutsPi+=3020;//default settings
  Float_t nsigmaPiTPC=0.1*(TrackCutsPi%100);
  Float_t nsigmaPiTOF=0.1*((TrackCutsPi/100)%100);
  Int_t CutTypePi=(TrackCutsPi/10000)%100;//0=TPC+TOFveto (default), 1=TPC only, 2=TOF only, 3 TPC+TOFcut
  Int_t MisidentifiedAsKaon=(TrackCutsPi/1000000)%10;//0=pion assigned pion mass, 1=pion assigned kaon mass

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

  Int_t iCutPi=task->AddTrackCuts(cutSetPi);

  // set cuts for Lambda from Sigma*

  Int_t LambdaCuts=TrackCutsS%10000;

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

  if(LambdaCuts==1) lambdaDCA=1.e10;
  else if(LambdaCuts==2) lambdaDaughDCA=0.5;

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
  }

  // AliRsnMiniResonanceFinder

  AliRsnMiniResonanceFinder* finder[8];
  Int_t i,iCutSig[8];

  AliRsnCutMiniPair* cutMassS=new AliRsnCutMiniPair("cutMassSigmastar",AliRsnCutMiniPair::kMassRange);
  cutMassS->SetRangeD(1.308,1.466);
  AliRsnCutMiniPair* cutYS=new AliRsnCutMiniPair("cutRapiditySigmastar",AliRsnCutMiniPair::kRapidityRange);
  cutYS->SetRangeD(-0.6,0.6);
  AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
  AliRsnCutSet* cutsS=new AliRsnCutSet("pairCutsSigmastar",AliRsnTarget::kMother);
  cutsS->AddCut(cutMassS);
  cutsS->AddCut(cutYS);
  cutsS->AddCut(cutV0);
  cutsS->SetCutScheme(TString::Format("%s&%s&(!%s)",cutMassS->GetName(),cutYS->GetName(),cutV0->GetName()).Data());

  i=0; // Sigma*+
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_Sigmastarp",task->GetName()));
  finder[i]->SetCutID(0,iCutLambda);
  finder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  finder[i]->SetCharge(0,'0');
  finder[i]->SetCutID(1,iCutPi);
  finder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.3828);
  finder[i]->SetResonancePDG(3224);
  finder[i]->SetPairCuts(cutsS);
  iCutSig[i]=task->AddResonanceFinder(finder[i]);

  i=1; // anti-Sigma*-
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_AntiSigmastarm",task->GetName()));
  finder[i]->SetCutID(0,iCutAntiLambda);
  finder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  finder[i]->SetCharge(0,'0');
  finder[i]->SetCutID(1,iCutPi);
  finder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  finder[i]->SetCharge(1,'-');
  finder[i]->SetResonanceMass(1.3828);
  finder[i]->SetResonancePDG(-3224);
  finder[i]->SetPairCuts(cutsS);
  iCutSig[i]=task->AddResonanceFinder(finder[i]);

  i=2; // Sigma*-
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_Sigmastarm",task->GetName()));
  finder[i]->SetCutID(0,iCutLambda);
  finder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  finder[i]->SetCharge(0,'0');
  finder[i]->SetCutID(1,iCutPi);
  finder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  finder[i]->SetCharge(1,'-');
  finder[i]->SetResonanceMass(1.3872);
  finder[i]->SetResonancePDG(3114);
  finder[i]->SetPairCuts(cutsS);
  iCutSig[i]=task->AddResonanceFinder(finder[i]);

  i=3; // anti-Sigma*+
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_AntiSigmastarp",task->GetName()));
  finder[i]->SetCutID(0,iCutAntiLambda);
  finder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  finder[i]->SetCharge(0,'0');
  finder[i]->SetCutID(1,iCutPi);
  finder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.3872);
  finder[i]->SetResonancePDG(-3114);
  finder[i]->SetPairCuts(cutsS);
  iCutSig[i]=task->AddResonanceFinder(finder[i]);
    
  // sidebands
  AliRsnCutMiniPair* cutMassSB=new AliRsnCutMiniPair("cutMassSB",AliRsnCutMiniPair::kMassRange);
  cutMassSB->SetRangeD(1.505,1.585);
  AliRsnCutSet* cutsSB=new AliRsnCutSet("pairCutsSB",AliRsnTarget::kMother);
  cutsSB->AddCut(cutMassSB);
  cutsSB->AddCut(cutYS);
  cutsSB->AddCut(cutV0);
  cutsSB->SetCutScheme(TString::Format("%s&%s&(!%s)",cutMassSB->GetName(),cutYS->GetName(),cutV0->GetName()).Data());

  i=4; // Sigma*+ sideband
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBpp",task->GetName()));
  finder[i]->SetCutID(0,iCutLambda);
  finder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  finder[i]->SetCharge(0,'0');
  finder[i]->SetCutID(1,iCutPi);
  finder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.3828);
  finder[i]->SetResonancePDG(3224);
  finder[i]->SetPairCuts(cutsSB);
  iCutSig[i]=task->AddResonanceFinder(finder[i]);

  i=5; // anti-Sigma*- sideband
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBma",task->GetName()));
  finder[i]->SetCutID(0,iCutAntiLambda);
  finder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  finder[i]->SetCharge(0,'0');
  finder[i]->SetCutID(1,iCutPi);
  finder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  finder[i]->SetCharge(1,'-');
  finder[i]->SetResonanceMass(1.3828);
  finder[i]->SetResonancePDG(-3224);
  finder[i]->SetPairCuts(cutsSB);
  iCutSig[i]=task->AddResonanceFinder(finder[i]);

  i=6; // Sigma*- sideband
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBmp",task->GetName()));
  finder[i]->SetCutID(0,iCutLambda);
  finder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  finder[i]->SetCharge(0,'0');
  finder[i]->SetCutID(1,iCutPi);
  finder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  finder[i]->SetCharge(1,'-');
  finder[i]->SetResonanceMass(1.3872);
  finder[i]->SetResonancePDG(3114);
  finder[i]->SetPairCuts(cutsSB);
  iCutSig[i]=task->AddResonanceFinder(finder[i]);

  i=7; // anti-Sigma*+ sideband
  finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBpa",task->GetName()));
  finder[i]->SetCutID(0,iCutAntiLambda);
  finder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  finder[i]->SetCharge(0,'0');
  finder[i]->SetCutID(1,iCutPi);
  finder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  finder[i]->SetCharge(1,'+');
  finder[i]->SetResonanceMass(1.3872);
  finder[i]->SetResonancePDG(-3114);
  finder[i]->SetPairCuts(cutsSB);
  iCutSig[i]=task->AddResonanceFinder(finder[i]);

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  cutY->SetRangeD(-0.5,0.5);
   
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
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass,    kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);

  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges

  Int_t k,xID,cut2,pairID,ipdg;
  AliRsnDaughter::ESpecies d2=AliRsnDaughter::kUnknown;
  TString name,comp;
  Char_t charge1='0',charge2;
  AliRsnMiniOutput* out;

  for(j=0;j<4;j++) for(k=0;k<3;k++){
    if(!j){
      name.Form("K0Sigmastarpp");
      d2=AliRsnDaughter::kSigmastarp;
      charge2='+';
      cut2=0;
    }else if(j==1){
      name.Form("K0Sigmastarma");
      d2=AliRsnDaughter::kSigmastarp;
      charge2='-';
      cut2=1;
    }else if(j==2){
      name.Form("K0Sigmastarmp");
      d2=AliRsnDaughter::kSigmastarm;
      charge2='-';
      cut2=2;
    }else{
      name.Form("K0Sigmastarpa");
      d2=AliRsnDaughter::kSigmastarm;
      charge2='+';
      cut2=3;
    }

    if(!k){
      comp.Form("PAIR");
    }else if(k==1){
      name.Append("SB");
      comp.Form("PAIR");
      cut2+=4;
    }else if(k==2){
      name.Append("Mix");
      comp.Form("MIX");
    }

    out=task->CreateOutput(Form("k0Sigmastar_%s%s",name.Data(),suffix),"HIST",comp.Data());
    out->SetDaughter(0,AliRsnDaughter::kKaon0);
    out->SetCutID(0,iCutK0s);
    out->SetCharge(0,charge1);

    out->SetDaughter(1,d2);
    out->SetCutID(1,iCutSig[cut2]);
    out->SetCharge(1,charge2);
    if(k!=1) out->SetUseStoredMass(1);

    if(k<=1) out->SetPairCuts(cutsPairSame);
    else out->SetPairCuts(cutsPairMix);
    out->SetMotherPDG(3124);
    out->SetMotherMass(mass);

    out->AddAxis(imID,260,1.7,3);// axis X: invmass or resolution
    //out->AddAxis(resID,200,-0.02,0.02);
    out->AddAxis(ptID,50,0.0,20.0);// axis Y: transverse momentum
    out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
  }

  // fill monitoring histogram for the resonance (Sigma*)
  for(i=0;i<4;i++) for(j=0;j<4;j++){
    if(!i) name.Form("Sigmastarpp");
    else if(i==1) name.Form("Sigmastarma");
    else if(i==2) name.Form("Sigmastarmp");
    else if(i==3) name.Form("Sigmastarpa");

    if(!j){
      name.Append("_mass");
      comp.Form("PAIR");
    }else if(j==1){
      name.Append("_SBmass");
      comp.Form("PAIR");
    }else if(j==2){
      if(!isMC) continue;
      name.Append("_genmass");
      comp.Form("MOTHER");
    }else if(j==3){
      if(!isMC) continue;
      name.Append("_recmass");
      comp.Form("TRUE");
    }

    k=i;
    if(j==1) k+=4;

    out=task->CreateOutput(Form("k0Sigmastar_%s",name.Data()),"HIST",comp.Data());
    out->SetMotherPDG(finder[k]->GetResonancePDG());

    out->SetDaughter(0,finder[k]->GetDaughter(0));
    out->SetCutID(0,finder[k]->GetCutID(0));
    out->SetCharge(0,finder[k]->GetCharge(0));

    out->SetDaughter(1,finder[k]->GetDaughter(1));
    out->SetCutID(1,finder[k]->GetCutID(1));
    out->SetCharge(1,finder[k]->GetCharge(1));

    out->SetMotherMass(finder[k]->GetResonanceMass());
    if(j!=1) out->SetPairCuts(cutsS);
    else out->SetPairCuts(cutsSB);

    out->AddAxis(imID,150,1.3,1.6);
    out->AddAxis(ptID,50,0.0,20.0);
    out->AddAxis(centID,nmult,multbins);
  }

  return kTRUE;
}


//=============================


Bool_t Config_kstar0Sigmastar(
  AliRsnMiniAnalysisTask *task,
  TString     lname,
  Bool_t      isMC,
  Int_t       system,
  Int_t       EventCuts,
  Int_t       TrackCutsK,
  Int_t       TrackCutsS
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;
  if(system==1 || system==2) MultBins=1;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;

  Double_t mass=0.896+1.385;

  // set cuts for pions and kaons
  Int_t TrackCutsPi=TrackCutsS/10000;
  if(!(TrackCutsPi%10000)) TrackCutsPi+=3020;//default settings
  Float_t nsigmaPiTPC=0.1*(TrackCutsPi%100);
  Float_t nsigmaPiTOF=0.1*((TrackCutsPi/100)%100);
  Int_t CutTypePi=(TrackCutsPi/10000)%100;//0=TPC+TOFveto (default), 1=TPC only, 2=TOF only, 3 TPC+TOFcut

  Int_t SidebandKstar=(TrackCutsK/10000)%10;
  Int_t TrackCutsKx=TrackCutsK%10000;
  if(!(TrackCutsKx)) TrackCutsKx+=3020;//default settings
  Float_t nsigmaKTPC=0.1*(TrackCutsKx%100);
  Float_t nsigmaKTOF=0.1*((TrackCutsKx/100)%100);

  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  //AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);

  AliRsnCutSetDaughterParticle* cutSetPi=0;
  if(!CutTypePi) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaPiTPC),
                                                           trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kPion,nsigmaPiTPC,nsigmaPiTOF);
  else if(CutTypePi==1) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,nsigmaPiTPC),
                                                                  trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kPion,nsigmaPiTPC,-1.);
  else if(CutTypePi==2) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,nsigmaPiTOF),
                                                                  trkQualityCut,AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,AliPID::kPion,-1.,nsigmaPiTOF);
  else if(CutTypePi==3) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFtightPidKStarPPB2011,nsigmaPiTPC),
                                                                  trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFtightPidKStarPPB2011,AliPID::kPion,nsigmaPiTPC,-1.);
  if(!cutSetPi){cerr<<"Error in AddTaskRare_pp13::Config_Kstar0Sigmastar(): missing cutSetPi"<<endl; return kFALSE;}
    
  AliRsnCutSetDaughterParticle* cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaKTPC),
                                                                         trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kKaon,nsigmaKTPC,nsigmaKTOF);

  //Int_t iCutQ=task->AddTrackCuts(cutSetQ);
  Int_t iCutPi=task->AddTrackCuts(cutSetPi);
  Int_t iCutK=task->AddTrackCuts(cutSetK);

  // set cuts for Lambda from Sigma*

  Int_t V0Cuts=TrackCutsS%10000;

  // selections for V0 daughters
  Int_t v0d_xrows=70;
  Float_t v0d_rtpc=0.8;
  Float_t v0d_dcaxy=0.06;

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts("qualityDaughterLambda");
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
    //AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
    AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());
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

  // AliRsnMiniResonanceFinder - K*0
  AliRsnMiniResonanceFinder* Kfinder[4];
  Int_t i,iCutKstar[4];

  AliRsnCutMiniPair* cutYRes=new AliRsnCutMiniPair("cutRapidityResonance",AliRsnCutMiniPair::kRapidityRange);
  cutYRes->SetRangeD(-0.6,0.6);

  AliRsnCutMiniPair* cutMassKstar=new AliRsnCutMiniPair("cutMassKstar",AliRsnCutMiniPair::kMassRange);
  if(!SidebandKstar) cutMassKstar->SetRangeD(0.8,0.99);
  else cutMassKstar->SetRangeD(1.037,1.085);
  AliRsnCutSet* cutsKstar=new AliRsnCutSet("pairCutsKstar",AliRsnTarget::kMother);
  cutsKstar->AddCut(cutMassKstar);
  cutsKstar->AddCut(cutYRes);
  cutsKstar->SetCutScheme(TString::Format("%s&%s",cutMassKstar->GetName(),cutYRes->GetName()).Data());
  // NOTE: An AliRsnMiniParticle only has space for 16 cut bits. For this reason, it was not possible to work
  // with all cuts at once: pi, K, Lambda, K*0 (with sidebands), and Sigma* (with sidebands). Instead, one
  // must run the K*0 sidebands separately.

  i=0; // K*0
  Kfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_Kstar0p",task->GetName()));
  Kfinder[i]->SetCutID(0,iCutK);
  Kfinder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  Kfinder[i]->SetCharge(0,'+');
  Kfinder[i]->SetCutID(1,iCutPi);
  Kfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Kfinder[i]->SetCharge(1,'-');
  Kfinder[i]->SetResonanceMass(0.89555);
  Kfinder[i]->SetResonancePDG(313);
  Kfinder[i]->SetPairCuts(cutsKstar);
  iCutKstar[i]=task->AddResonanceFinder(Kfinder[i]);

  i=1; // anti-K*0
  Kfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_Kstar0a",task->GetName()));
  Kfinder[i]->SetCutID(0,iCutK);
  Kfinder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  Kfinder[i]->SetCharge(0,'-');
  Kfinder[i]->SetCutID(1,iCutPi);
  Kfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Kfinder[i]->SetCharge(1,'+');
  Kfinder[i]->SetResonanceMass(0.89555);
  Kfinder[i]->SetResonancePDG(-313);
  Kfinder[i]->SetPairCuts(cutsKstar);
  iCutKstar[i]=task->AddResonanceFinder(Kfinder[i]);

  i=2; // K+ pi+
  Kfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_KpPip",task->GetName()));
  Kfinder[i]->SetCutID(0,iCutK);
  Kfinder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  Kfinder[i]->SetCharge(0,'+');
  Kfinder[i]->SetCutID(1,iCutPi);
  Kfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Kfinder[i]->SetCharge(1,'+');
  Kfinder[i]->SetResonanceMass(0.89555);
  Kfinder[i]->SetResonancePDG(313);
  Kfinder[i]->SetPairCuts(cutsKstar);
  iCutKstar[i]=task->AddResonanceFinder(Kfinder[i]);

  i=3; // K- pi-
  Kfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_KmPim",task->GetName()));
  Kfinder[i]->SetCutID(0,iCutK);
  Kfinder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
  Kfinder[i]->SetCharge(0,'-');
  Kfinder[i]->SetCutID(1,iCutPi);
  Kfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Kfinder[i]->SetCharge(1,'-');
  Kfinder[i]->SetResonanceMass(0.89555);
  Kfinder[i]->SetResonancePDG(-313);
  Kfinder[i]->SetPairCuts(cutsKstar);
  iCutKstar[i]=task->AddResonanceFinder(Kfinder[i]);

  // AliRsnMiniResonanceFinder - Sigma*
  AliRsnMiniResonanceFinder* Sfinder[8];
  Int_t iCutSig[8];

  AliRsnCutMiniPair* cutMassS=new AliRsnCutMiniPair("cutMassSigmastar",AliRsnCutMiniPair::kMassRange);
  cutMassS->SetRangeD(1.308,1.466);
  AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
  AliRsnCutSet* cutsS=new AliRsnCutSet("pairCutsSigmastar",AliRsnTarget::kMother);
  cutsS->AddCut(cutMassS);
  cutsS->AddCut(cutYRes);
  cutsS->AddCut(cutV0);
  cutsS->SetCutScheme(TString::Format("%s&%s&(!%s)",cutMassS->GetName(),cutYRes->GetName(),cutV0->GetName()).Data());

  i=0; // Sigma*+
  Sfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_Sigmastarp",task->GetName()));
  Sfinder[i]->SetCutID(0,iCutLambda);
  Sfinder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  Sfinder[i]->SetCharge(0,'0');
  Sfinder[i]->SetCutID(1,iCutPi);
  Sfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Sfinder[i]->SetCharge(1,'+');
  Sfinder[i]->SetResonanceMass(1.3828);
  Sfinder[i]->SetResonancePDG(3224);
  Sfinder[i]->SetPairCuts(cutsS);
  iCutSig[i]=task->AddResonanceFinder(Sfinder[i]);

  i=1; // anti-Sigma*-
  Sfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_AntiSigmastarm",task->GetName()));
  Sfinder[i]->SetCutID(0,iCutAntiLambda);
  Sfinder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  Sfinder[i]->SetCharge(0,'0');
  Sfinder[i]->SetCutID(1,iCutPi);
  Sfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Sfinder[i]->SetCharge(1,'-');
  Sfinder[i]->SetResonanceMass(1.3828);
  Sfinder[i]->SetResonancePDG(-3224);
  Sfinder[i]->SetPairCuts(cutsS);
  iCutSig[i]=task->AddResonanceFinder(Sfinder[i]);

  i=2; // Sigma*-
  Sfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_Sigmastarm",task->GetName()));
  Sfinder[i]->SetCutID(0,iCutLambda);
  Sfinder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  Sfinder[i]->SetCharge(0,'0');
  Sfinder[i]->SetCutID(1,iCutPi);
  Sfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Sfinder[i]->SetCharge(1,'-');
  Sfinder[i]->SetResonanceMass(1.3872);
  Sfinder[i]->SetResonancePDG(3114);
  Sfinder[i]->SetPairCuts(cutsS);
  iCutSig[i]=task->AddResonanceFinder(Sfinder[i]);

  i=3; // anti-Sigma*+
  Sfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_AntiSigmastarp",task->GetName()));
  Sfinder[i]->SetCutID(0,iCutAntiLambda);
  Sfinder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  Sfinder[i]->SetCharge(0,'0');
  Sfinder[i]->SetCutID(1,iCutPi);
  Sfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Sfinder[i]->SetCharge(1,'+');
  Sfinder[i]->SetResonanceMass(1.3872);
  Sfinder[i]->SetResonancePDG(-3114);
  Sfinder[i]->SetPairCuts(cutsS);
  iCutSig[i]=task->AddResonanceFinder(Sfinder[i]);
    
  // sidebands
  AliRsnCutMiniPair* cutMassSigmastarSB=new AliRsnCutMiniPair("cutMassSigmastarSB",AliRsnCutMiniPair::kMassRange);
  cutMassSigmastarSB->SetRangeD(1.505,1.585);
  AliRsnCutSet* cutsSigmastarSB=new AliRsnCutSet("pairCutsSB",AliRsnTarget::kMother);
  cutsSigmastarSB->AddCut(cutMassSigmastarSB);
  cutsSigmastarSB->AddCut(cutYRes);
  cutsSigmastarSB->AddCut(cutV0);
  cutsSigmastarSB->SetCutScheme(TString::Format("%s&%s&(!%s)",cutMassSigmastarSB->GetName(),cutYRes->GetName(),cutV0->GetName()).Data());

  i=4; // Sigma*+ sideband
  Sfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBpp",task->GetName()));
  Sfinder[i]->SetCutID(0,iCutLambda);
  Sfinder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  Sfinder[i]->SetCharge(0,'0');
  Sfinder[i]->SetCutID(1,iCutPi);
  Sfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Sfinder[i]->SetCharge(1,'+');
  Sfinder[i]->SetResonanceMass(1.3828);
  Sfinder[i]->SetResonancePDG(3224);
  Sfinder[i]->SetPairCuts(cutsSigmastarSB);
  iCutSig[i]=task->AddResonanceFinder(Sfinder[i]);

  i=5; // anti-Sigma*- sideband
  Sfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBma",task->GetName()));
  Sfinder[i]->SetCutID(0,iCutAntiLambda);
  Sfinder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  Sfinder[i]->SetCharge(0,'0');
  Sfinder[i]->SetCutID(1,iCutPi);
  Sfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Sfinder[i]->SetCharge(1,'-');
  Sfinder[i]->SetResonanceMass(1.3828);
  Sfinder[i]->SetResonancePDG(-3224);
  Sfinder[i]->SetPairCuts(cutsSigmastarSB);
  iCutSig[i]=task->AddResonanceFinder(Sfinder[i]);

  i=6; // Sigma*- sideband
  Sfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBmp",task->GetName()));
  Sfinder[i]->SetCutID(0,iCutLambda);
  Sfinder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  Sfinder[i]->SetCharge(0,'0');
  Sfinder[i]->SetCutID(1,iCutPi);
  Sfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Sfinder[i]->SetCharge(1,'-');
  Sfinder[i]->SetResonanceMass(1.3872);
  Sfinder[i]->SetResonancePDG(3114);
  Sfinder[i]->SetPairCuts(cutsSigmastarSB);
  iCutSig[i]=task->AddResonanceFinder(Sfinder[i]);

  i=7; // anti-Sigma*+ sideband
  Sfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBpa",task->GetName()));
  Sfinder[i]->SetCutID(0,iCutAntiLambda);
  Sfinder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  Sfinder[i]->SetCharge(0,'0');
  Sfinder[i]->SetCutID(1,iCutPi);
  Sfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Sfinder[i]->SetCharge(1,'+');
  Sfinder[i]->SetResonanceMass(1.3872);
  Sfinder[i]->SetResonancePDG(-3114);
  Sfinder[i]->SetPairCuts(cutsSigmastarSB);
  iCutSig[i]=task->AddResonanceFinder(Sfinder[i]);

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  cutY->SetRangeD(-0.5,0.5);
   
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
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass,    kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);

  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges

  Int_t k,xID,cut1,cut2,pairID,ipdg;
  AliRsnDaughter::ESpecies d2=AliRsnDaughter::kUnknown;
  TString name,comp;
  Char_t charge2;
  AliRsnMiniOutput* out;

  for(i=0;i<4;i++) for(j=0;j<4;j++) for(k=0;k<3;k++){
    if(!i){
      name.Form("Kstar0p");
      cut1=0;
    }else if(i==1){
      name.Form("Kstar0a");
      cut1=1;
    }else if(i==2){
      name.Form("KpPip");
      cut1=2;
    }else if(i==3){
      name.Form("KmPim");
      cut1=3;
    }

    if(!j){
      name.Append("Sigmastarpp");
      d2=AliRsnDaughter::kSigmastarp;
      charge2='+';
      cut2=0;
    }else if(j==1){
      name.Append("Sigmastarma");
      d2=AliRsnDaughter::kSigmastarp;
      charge2='-';
      cut2=1;
    }else if(j==2){
      name.Append("Sigmastarmp");
      d2=AliRsnDaughter::kSigmastarm;
      charge2='-';
      cut2=2;
    }else{
      name.Append("Sigmastarpa");
      d2=AliRsnDaughter::kSigmastarm;
      charge2='+';
      cut2=3;
    }

    if(!k){
      comp.Form("PAIR");
    }else if(k==1){
      name.Append("SB");
      comp.Form("PAIR");
      cut2+=4;
    }else if(k==2){
      name.Append("Mix");
      comp.Form("MIX");
    }
      
    if((i>=2 || SidebandKstar) && k) continue;

    out=task->CreateOutput(Form("Kstar0Sigmastar_%s%s",name.Data(),suffix),"HIST",comp.Data());
    out->SetDaughter(0,AliRsnDaughter::kKstar0);
    out->SetCutID(0,iCutKstar[cut1]);
    out->SetCharge(0,'0');
    if(!SidebandKstar) out->SetUseStoredMass(0);

    out->SetDaughter(1,d2);
    out->SetCutID(1,iCutSig[cut2]);
    out->SetCharge(1,charge2);
    if(k!=1) out->SetUseStoredMass(1);

    if(k<=1) out->SetPairCuts(cutsPairSame);
    else out->SetPairCuts(cutsPairMix);
    out->SetMotherPDG(3124);
    out->SetMotherMass(mass);

    out->AddAxis(imID,570,2.36,3.5);// axis X: invmass or resolution
    //out->AddAxis(resID,200,-0.02,0.02);
    out->AddAxis(ptID,50,0.0,20.0);// axis Y: transverse momentum
    out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
  }

  // fill monitoring histogram for the K*0
  for(i=0;i<4;i++) for(j=0;j<3;j++){
    if(!i) name.Form("Kstar0p");
    else if(i==1) name.Form("Kstar0a");
    else if(i==2) name.Form("KpPip");
    else if(i==3) name.Form("KmPim");

    if(!j){
      name.Append("_mass");
      comp.Form("PAIR");
    }else if(j==1){
      if(!isMC) continue;
      name.Append("_genmass");
      comp.Form("MOTHER");
    }else if(j==2){
      if(!isMC) continue;
      name.Append("_recmass");
      comp.Form("TRUE");
    }

    if(i>=2 && (j || SidebandKstar)) continue;

    out=task->CreateOutput(Form("Kstar0Sigmastar_%s",name.Data()),"HIST",comp.Data());
    out->SetMotherPDG(Kfinder[i]->GetResonancePDG());

    out->SetDaughter(0,Kfinder[i]->GetDaughter(0));
    out->SetCutID(0,Kfinder[i]->GetCutID(0));
    out->SetCharge(0,Kfinder[i]->GetCharge(0));

    out->SetDaughter(1,Kfinder[i]->GetDaughter(1));
    out->SetCutID(1,Kfinder[i]->GetCutID(1));
    out->SetCharge(1,Kfinder[i]->GetCharge(1));

    out->SetMotherMass(Kfinder[i]->GetResonanceMass());
    out->SetPairCuts(cutsKstar);

    out->AddAxis(imID,170,0.75,1.09);
    out->AddAxis(ptID,50,0.0,20.0);
    out->AddAxis(centID,nmult,multbins);
  }

  // fill monitoring histogram for the Sigma*
  for(i=0;i<4;i++) for(j=0;j<4;j++){
    if(!i) name.Form("Sigmastarpp");
    else if(i==1) name.Form("Sigmastarma");
    else if(i==2) name.Form("Sigmastarmp");
    else if(i==3) name.Form("Sigmastarpa");

    if(!j){
      name.Append("_mass");
      comp.Form("PAIR");
    }else if(j==1){
      name.Append("_SBmass");
      comp.Form("PAIR");
    }else if(j==2){
      if(!isMC) continue;
      name.Append("_genmass");
      comp.Form("MOTHER");
    }else if(j==3){
      if(!isMC) continue;
      name.Append("_recmass");
      comp.Form("TRUE");
    }

    k=i;
    if(j==1) k+=4;

    out=task->CreateOutput(Form("Kstar0Sigmastar_%s",name.Data()),"HIST",comp.Data());
    out->SetMotherPDG(Sfinder[k]->GetResonancePDG());

    out->SetDaughter(0,Sfinder[k]->GetDaughter(0));
    out->SetCutID(0,Sfinder[k]->GetCutID(0));
    out->SetCharge(0,Sfinder[k]->GetCharge(0));

    out->SetDaughter(1,Sfinder[k]->GetDaughter(1));
    out->SetCutID(1,Sfinder[k]->GetCutID(1));
    out->SetCharge(1,Sfinder[k]->GetCharge(1));

    out->SetMotherMass(Sfinder[k]->GetResonanceMass());
    if(j!=1) out->SetPairCuts(cutsS);
    else out->SetPairCuts(cutsSigmastarSB);

    out->AddAxis(imID,150,1.3,1.6);
    out->AddAxis(ptID,50,0.0,20.0);
    out->AddAxis(centID,nmult,multbins);
  }

  return kTRUE;
}


//=============================


Bool_t Config_kstarxSigmastar(
  AliRsnMiniAnalysisTask *task,
  TString     lname,
  Bool_t      isMC,
  Int_t       system,
  Int_t       EventCuts,
  Int_t       TrackCutsK,
  Int_t       TrackCutsS
){
  bool isPP=false;
  if(!system) isPP=true;
  int trigger=EventCuts%10;
  int MultBins=(EventCuts/10)%10;
  if(system==1 || system==2) MultBins=1;

  char suffix[1000];
  sprintf(suffix,"_%s",lname.Data());
  Bool_t enableMonitor=kTRUE;

  Double_t mass=0.892+1.385;

  // set cuts for pions
  Int_t TrackCutsPi=TrackCutsS/10000;
  if(!(TrackCutsPi%10000)) TrackCutsPi+=3020;//default settings
  Float_t nsigmaPiTPC=0.1*(TrackCutsPi%100);
  Float_t nsigmaPiTOF=0.1*((TrackCutsPi/100)%100);
  Int_t CutTypePi=(TrackCutsPi/10000)%100;//0=TPC+TOFveto (default), 1=TPC only, 2=TOF only, 3 TPC+TOFcut
    
  AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
  trkQualityCut->SetDefaults2011(kTRUE,kTRUE);

  //AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);

  AliRsnCutSetDaughterParticle* cutSetPi=0;
  if(!CutTypePi) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaPiTPC),
                                                           trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kPion,nsigmaPiTPC,nsigmaPiTOF);
  else if(CutTypePi==1) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,nsigmaPiTPC),
                                                                  trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kPion,nsigmaPiTPC,-1.);
  else if(CutTypePi==2) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,nsigmaPiTOF),
                                                                  trkQualityCut,AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,AliPID::kPion,-1.,nsigmaPiTOF);
  else if(CutTypePi==3) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFtightPidKStarPPB2011,nsigmaPiTPC),
                                                                  trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFtightPidKStarPPB2011,AliPID::kPion,nsigmaPiTPC,-1.);
  if(!cutSetPi){cerr<<"Error in AddTaskRare_pp13::Config_Kstar0Sigmastar(): missing cutSetPi"<<endl; return kFALSE;}

  //Int_t iCutQ=task->AddTrackCuts(cutSetQ);
  Int_t iCutPi=task->AddTrackCuts(cutSetPi);

  // set cuts for K0S from K*

  Int_t SidebandKstar=(TrackCutsK/10000)%10;
  Int_t V0Cuts=TrackCutsK%10000;

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

  // set cuts for Lambda from Sigma*

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
  TString pname="k0s";
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
    //AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
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
  }

  // AliRsnMiniResonanceFinder - K*
  AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
  AliRsnCutMiniPair* cutYRes=new AliRsnCutMiniPair("cutRapidityResonance",AliRsnCutMiniPair::kRapidityRange);
  cutYRes->SetRangeD(-0.6,0.6);

  AliRsnMiniResonanceFinder* Kfinder[4];
  Int_t i,iCutKstar[4];

  AliRsnCutMiniPair* cutMassKstar=new AliRsnCutMiniPair("cutMassKstar",AliRsnCutMiniPair::kMassRange);
  cutMassKstar->SetRangeD(0.8,0.99);
  AliRsnCutSet* cutsKstar=new AliRsnCutSet("pairCutsKstar",AliRsnTarget::kMother);
  cutsKstar->AddCut(cutMassKstar);
  cutsKstar->AddCut(cutYRes);
  cutsKstar->AddCut(cutV0);
  cutsKstar->SetCutScheme(TString::Format("%s&%s&(!%s)",cutMassKstar->GetName(),cutYRes->GetName(),cutV0->GetName()).Data());

  i=0; // K*+
  Kfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_Kstarp",task->GetName()));
  Kfinder[i]->SetCutID(0,iCutK0s);
  Kfinder[i]->SetDaughter(0,AliRsnDaughter::kKaon0);
  Kfinder[i]->SetCharge(0,'0');
  Kfinder[i]->SetCutID(1,iCutPi);
  Kfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Kfinder[i]->SetCharge(1,'+');
  Kfinder[i]->SetResonanceMass(0.89176);
  Kfinder[i]->SetResonancePDG(323);
  Kfinder[i]->SetPairCuts(cutsKstar);
  iCutKstar[i]=task->AddResonanceFinder(Kfinder[i]);

  i=1; // K*-
  Kfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_Kstarm",task->GetName()));
  Kfinder[i]->SetCutID(0,iCutK0s);
  Kfinder[i]->SetDaughter(0,AliRsnDaughter::kKaon0);
  Kfinder[i]->SetCharge(0,'0');
  Kfinder[i]->SetCutID(1,iCutPi);
  Kfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Kfinder[i]->SetCharge(1,'-');
  Kfinder[i]->SetResonanceMass(0.89176);
  Kfinder[i]->SetResonancePDG(-323);
  Kfinder[i]->SetPairCuts(cutsKstar);
  iCutKstar[i]=task->AddResonanceFinder(Kfinder[i]);

  AliRsnCutMiniPair* cutMassKstarSB=new AliRsnCutMiniPair("cutMassKstarSB",AliRsnCutMiniPair::kMassRange);
  cutMassKstarSB->SetRangeD(1.037,1.085);
  AliRsnCutSet* cutsKstarSB=new AliRsnCutSet("pairCutsKstarSB",AliRsnTarget::kMother);
  cutsKstarSB->AddCut(cutMassKstarSB);
  cutsKstarSB->AddCut(cutYRes);
  cutsKstarSB->AddCut(cutV0);
  cutsKstarSB->SetCutScheme(TString::Format("%s&%s&(!%s)",cutMassKstarSB->GetName(),cutYRes->GetName(),cutV0->GetName()).Data());

  i=2; // K*+ sideband
  Kfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_KstarpSB",task->GetName()));
  Kfinder[i]->SetCutID(0,iCutK0s);
  Kfinder[i]->SetDaughter(0,AliRsnDaughter::kKaon0);
  Kfinder[i]->SetCharge(0,'0');
  Kfinder[i]->SetCutID(1,iCutPi);
  Kfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Kfinder[i]->SetCharge(1,'+');
  Kfinder[i]->SetResonanceMass(0.89176);
  Kfinder[i]->SetResonancePDG(323);
  Kfinder[i]->SetPairCuts(cutsKstarSB);
  iCutKstar[i]=task->AddResonanceFinder(Kfinder[i]);

  i=3; // K*- sideband
  Kfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_KstarmSB",task->GetName()));
  Kfinder[i]->SetCutID(0,iCutK0s);
  Kfinder[i]->SetDaughter(0,AliRsnDaughter::kKaon0);
  Kfinder[i]->SetCharge(0,'0');
  Kfinder[i]->SetCutID(1,iCutPi);
  Kfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Kfinder[i]->SetCharge(1,'-');
  Kfinder[i]->SetResonanceMass(0.89176);
  Kfinder[i]->SetResonancePDG(-323);
  Kfinder[i]->SetPairCuts(cutsKstarSB);
  iCutKstar[i]=task->AddResonanceFinder(Kfinder[i]);

  // AliRsnMiniResonanceFinder - Sigma*
  AliRsnMiniResonanceFinder* Sfinder[8];
  Int_t iCutSig[8];

  AliRsnCutMiniPair* cutMassS=new AliRsnCutMiniPair("cutMassSigmastar",AliRsnCutMiniPair::kMassRange);
  cutMassS->SetRangeD(1.308,1.466);
  AliRsnCutSet* cutsS=new AliRsnCutSet("pairCutsSigmastar",AliRsnTarget::kMother);
  cutsS->AddCut(cutMassS);
  cutsS->AddCut(cutYRes);
  cutsS->AddCut(cutV0);
  cutsS->SetCutScheme(TString::Format("%s&%s&(!%s)",cutMassS->GetName(),cutYRes->GetName(),cutV0->GetName()).Data());

  i=0; // Sigma*+
  Sfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_Sigmastarp",task->GetName()));
  Sfinder[i]->SetCutID(0,iCutLambda);
  Sfinder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  Sfinder[i]->SetCharge(0,'0');
  Sfinder[i]->SetCutID(1,iCutPi);
  Sfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Sfinder[i]->SetCharge(1,'+');
  Sfinder[i]->SetResonanceMass(1.3828);
  Sfinder[i]->SetResonancePDG(3224);
  Sfinder[i]->SetPairCuts(cutsS);
  iCutSig[i]=task->AddResonanceFinder(Sfinder[i]);

  i=1; // anti-Sigma*-
  Sfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_AntiSigmastarm",task->GetName()));
  Sfinder[i]->SetCutID(0,iCutAntiLambda);
  Sfinder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  Sfinder[i]->SetCharge(0,'0');
  Sfinder[i]->SetCutID(1,iCutPi);
  Sfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Sfinder[i]->SetCharge(1,'-');
  Sfinder[i]->SetResonanceMass(1.3828);
  Sfinder[i]->SetResonancePDG(-3224);
  Sfinder[i]->SetPairCuts(cutsS);
  iCutSig[i]=task->AddResonanceFinder(Sfinder[i]);

  i=2; // Sigma*-
  Sfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_Sigmastarm",task->GetName()));
  Sfinder[i]->SetCutID(0,iCutLambda);
  Sfinder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  Sfinder[i]->SetCharge(0,'0');
  Sfinder[i]->SetCutID(1,iCutPi);
  Sfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Sfinder[i]->SetCharge(1,'-');
  Sfinder[i]->SetResonanceMass(1.3872);
  Sfinder[i]->SetResonancePDG(3114);
  Sfinder[i]->SetPairCuts(cutsS);
  iCutSig[i]=task->AddResonanceFinder(Sfinder[i]);

  i=3; // anti-Sigma*+
  Sfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_AntiSigmastarp",task->GetName()));
  Sfinder[i]->SetCutID(0,iCutAntiLambda);
  Sfinder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  Sfinder[i]->SetCharge(0,'0');
  Sfinder[i]->SetCutID(1,iCutPi);
  Sfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Sfinder[i]->SetCharge(1,'+');
  Sfinder[i]->SetResonanceMass(1.3872);
  Sfinder[i]->SetResonancePDG(-3114);
  Sfinder[i]->SetPairCuts(cutsS);
  iCutSig[i]=task->AddResonanceFinder(Sfinder[i]);
    
  // sidebands
  AliRsnCutMiniPair* cutMassSigmastarSB=new AliRsnCutMiniPair("cutMassSigmastarSB",AliRsnCutMiniPair::kMassRange);
  cutMassSigmastarSB->SetRangeD(1.505,1.585);
  AliRsnCutSet* cutsSigmastarSB=new AliRsnCutSet("pairCutsSB",AliRsnTarget::kMother);
  cutsSigmastarSB->AddCut(cutMassSigmastarSB);
  cutsSigmastarSB->AddCut(cutYRes);
  cutsSigmastarSB->AddCut(cutV0);
  cutsSigmastarSB->SetCutScheme(TString::Format("%s&%s&(!%s)",cutMassSigmastarSB->GetName(),cutYRes->GetName(),cutV0->GetName()).Data());

  i=4; // Sigma*+ sideband
  Sfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBpp",task->GetName()));
  Sfinder[i]->SetCutID(0,iCutLambda);
  Sfinder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  Sfinder[i]->SetCharge(0,'0');
  Sfinder[i]->SetCutID(1,iCutPi);
  Sfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Sfinder[i]->SetCharge(1,'+');
  Sfinder[i]->SetResonanceMass(1.3828);
  Sfinder[i]->SetResonancePDG(3224);
  Sfinder[i]->SetPairCuts(cutsSigmastarSB);
  iCutSig[i]=task->AddResonanceFinder(Sfinder[i]);

  i=5; // anti-Sigma*- sideband
  Sfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBma",task->GetName()));
  Sfinder[i]->SetCutID(0,iCutAntiLambda);
  Sfinder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  Sfinder[i]->SetCharge(0,'0');
  Sfinder[i]->SetCutID(1,iCutPi);
  Sfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Sfinder[i]->SetCharge(1,'-');
  Sfinder[i]->SetResonanceMass(1.3828);
  Sfinder[i]->SetResonancePDG(-3224);
  Sfinder[i]->SetPairCuts(cutsSigmastarSB);
  iCutSig[i]=task->AddResonanceFinder(Sfinder[i]);

  i=6; // Sigma*- sideband
  Sfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBmp",task->GetName()));
  Sfinder[i]->SetCutID(0,iCutLambda);
  Sfinder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  Sfinder[i]->SetCharge(0,'0');
  Sfinder[i]->SetCutID(1,iCutPi);
  Sfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Sfinder[i]->SetCharge(1,'-');
  Sfinder[i]->SetResonanceMass(1.3872);
  Sfinder[i]->SetResonancePDG(3114);
  Sfinder[i]->SetPairCuts(cutsSigmastarSB);
  iCutSig[i]=task->AddResonanceFinder(Sfinder[i]);

  i=7; // anti-Sigma*+ sideband
  Sfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBpa",task->GetName()));
  Sfinder[i]->SetCutID(0,iCutAntiLambda);
  Sfinder[i]->SetDaughter(0,AliRsnDaughter::kLambda);
  Sfinder[i]->SetCharge(0,'0');
  Sfinder[i]->SetCutID(1,iCutPi);
  Sfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
  Sfinder[i]->SetCharge(1,'+');
  Sfinder[i]->SetResonanceMass(1.3872);
  Sfinder[i]->SetResonancePDG(-3114);
  Sfinder[i]->SetPairCuts(cutsSigmastarSB);
  iCutSig[i]=task->AddResonanceFinder(Sfinder[i]);

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  cutY->SetRangeD(-0.5,0.5);
   
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
  /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass,    kFALSE);
  /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
  /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
  /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
  /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);

  // -- Create all needed outputs -----------------------------------------------------------------
  // use an array for more compact writing, which are different on mixing and charges

  Int_t k,xID,cut1,cut2,pairID,ipdg;
  AliRsnDaughter::ESpecies d2=AliRsnDaughter::kUnknown;
  TString name,comp;
  Char_t charge1,charge2;
  AliRsnMiniOutput* out;

  for(i=0;i<4;i++) for(j=0;j<4;j++) for(k=0;k<3;k++){
    if(!i){
      name.Form("Kstarp");
      cut1=0;
      charge1='+';
    }else if(i==1){
      name.Form("Kstarm");
      cut1=1;
      charge1='-';
    }else if(i==2){
      name.Form("KstarpSB");
      cut1=2;
      charge1='+';
    }else if(i==3){
      name.Form("KstarmSB");
      cut1=3;
      charge1='-';
    }

    if(!j){
      name.Append("Sigmastarpp");
      d2=AliRsnDaughter::kSigmastarp;
      charge2='+';
      cut2=0;
    }else if(j==1){
      name.Append("Sigmastarma");
      d2=AliRsnDaughter::kSigmastarp;
      charge2='-';
      cut2=1;
    }else if(j==2){
      name.Append("Sigmastarmp");
      d2=AliRsnDaughter::kSigmastarm;
      charge2='-';
      cut2=2;
    }else{
      name.Append("Sigmastarpa");
      d2=AliRsnDaughter::kSigmastarm;
      charge2='+';
      cut2=3;
    }

    if(!k){
      comp.Form("PAIR");
    }else if(k==1){
      name.Append("SB");
      comp.Form("PAIR");
      cut2+=4;
    }else if(k==2){
      name.Append("Mix");
      comp.Form("MIX");
    }
      
    if(i>=2 && k) continue;

    out=task->CreateOutput(Form("KstarxSigmastar_%s%s",name.Data(),suffix),"HIST",comp.Data());
    out->SetDaughter(0,AliRsnDaughter::kKstarpm);
    out->SetCutID(0,iCutKstar[cut1]);
    out->SetCharge(0,charge1);
    if(i<=1) out->SetUseStoredMass(0);

    out->SetDaughter(1,d2);
    out->SetCutID(1,iCutSig[cut2]);
    out->SetCharge(1,charge2);
    if(k!=1) out->SetUseStoredMass(1);

    if(k<=1) out->SetPairCuts(cutsPairSame);
    else out->SetPairCuts(cutsPairMix);
    out->SetMotherPDG(3124);
    out->SetMotherMass(mass);

    out->AddAxis(imID,570,2.36,3.5);// axis X: invmass or resolution
    //out->AddAxis(resID,200,-0.02,0.02);
    out->AddAxis(ptID,50,0.0,20.0);// axis Y: transverse momentum
    out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
  }

  // fill monitoring histogram for the K*
  for(i=0;i<4;i++) for(j=0;j<3;j++){
    if(!i) name.Form("Kstarp");
    else if(i==1) name.Form("Kstarm");
    else if(i==2) name.Form("KstarpSB");
    else if(i==3) name.Form("KstarmSB");

    if(!j){
      name.Append("_mass");
      comp.Form("PAIR");
    }else if(j==1){
      if(!isMC) continue;
      name.Append("_genmass");
      comp.Form("MOTHER");
    }else if(j==2){
      if(!isMC) continue;
      name.Append("_recmass");
      comp.Form("TRUE");
    }

    if(i>=2 && j) continue;

    out=task->CreateOutput(Form("KstarxSigmastar_%s",name.Data()),"HIST",comp.Data());
    out->SetMotherPDG(Kfinder[i]->GetResonancePDG());

    out->SetDaughter(0,Kfinder[i]->GetDaughter(0));
    out->SetCutID(0,Kfinder[i]->GetCutID(0));
    out->SetCharge(0,Kfinder[i]->GetCharge(0));

    out->SetDaughter(1,Kfinder[i]->GetDaughter(1));
    out->SetCutID(1,Kfinder[i]->GetCutID(1));
    out->SetCharge(1,Kfinder[i]->GetCharge(1));

    out->SetMotherMass(Kfinder[i]->GetResonanceMass());
    if(i<=1) out->SetPairCuts(cutsKstar);
    else out->SetPairCuts(cutsKstarSB);

    out->AddAxis(imID,170,0.75,1.09);
    out->AddAxis(ptID,50,0.0,20.0);
    out->AddAxis(centID,nmult,multbins);
  }

  // fill monitoring histogram for the Sigma*
  for(i=0;i<4;i++) for(j=0;j<4;j++){
    if(!i) name.Form("Sigmastarpp");
    else if(i==1) name.Form("Sigmastarma");
    else if(i==2) name.Form("Sigmastarmp");
    else if(i==3) name.Form("Sigmastarpa");

    if(!j){
      name.Append("_mass");
      comp.Form("PAIR");
    }else if(j==1){
      name.Append("_SBmass");
      comp.Form("PAIR");
    }else if(j==2){
      if(!isMC) continue;
      name.Append("_genmass");
      comp.Form("MOTHER");
    }else if(j==3){
      if(!isMC) continue;
      name.Append("_recmass");
      comp.Form("TRUE");
    }

    k=i;
    if(j==1) k+=4;

    out=task->CreateOutput(Form("KstarxSigmastar_%s",name.Data()),"HIST",comp.Data());
    out->SetMotherPDG(Sfinder[k]->GetResonancePDG());

    out->SetDaughter(0,Sfinder[k]->GetDaughter(0));
    out->SetCutID(0,Sfinder[k]->GetCutID(0));
    out->SetCharge(0,Sfinder[k]->GetCharge(0));

    out->SetDaughter(1,Sfinder[k]->GetDaughter(1));
    out->SetCutID(1,Sfinder[k]->GetCutID(1));
    out->SetCharge(1,Sfinder[k]->GetCharge(1));

    out->SetMotherMass(Sfinder[k]->GetResonanceMass());
    if(j!=1) out->SetPairCuts(cutsS);
    else out->SetPairCuts(cutsSigmastarSB);

    out->AddAxis(imID,150,1.3,1.6);
    out->AddAxis(ptID,50,0.0,20.0);
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
