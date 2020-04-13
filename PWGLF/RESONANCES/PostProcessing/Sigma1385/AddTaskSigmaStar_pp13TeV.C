#ifdef __CLING__
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C>
#endif
/***************************************************************************
//            Modified by Bong-Hwi Lim - 08/05/2019
//            Based on AddTaskRare_pp13
//
// Macro to configure the SigmaStar analysis task 
//
// Input parameters:
// isMC: kTRUE if MC
// system: 0 for pp, 1 for pPb, 2 for PbPb
// EventCuts:      1*tigger(0 for INT7, 1 for HighMultV0) 
//            +   10*Multibins(10 for V0M, 20 for RefMult08)
//            + 1000*EventMixing(1000 for no mxing, n*1000 for n mixing, default:5)
// eg1. 10: INT7, V0M Mutliplicity, 5 mixing
// eg2. 20011: HighMultV0, V0M Multiplicity, 20 mixing
****************************************************************************/
enum ERsnCollType_t { kPP=0,
		      kPPb,
		      kPbPb};

void AddMonitorOutput_P(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_Pt(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_Eta(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_DCAxy(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
void AddMonitorOutput_DCAz(TString s="",TObjArray* m=0,TString o="",AliRsnLoopDaughter* l=0);
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

Bool_t Config_Lambdapi( // From Anders's master macro.
  AliRsnMiniAnalysisTask *task,
  TString     lname="Lambdapi",
  Bool_t      isMC=kFALSE,
  Int_t       system=0,
  Int_t       EventCuts=0,
  Int_t       TrackCutsLambda=0,
  Int_t       TrackCutsPi=0);
AliRsnMiniAnalysisTask* AddTaskSigmaStar_pp13TeV(
  TString lname,
  Bool_t isMC,
  Int_t system,
  Int_t EventCuts=0,
  Int_t TrackCuts1=0,
  Int_t TrackCuts2=0){
  // ----- INITIALIZATION -----

  // retrieve analysis manager
  AliAnalysisManager* mgr=AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    ::Error("AddTaskSigmaStar_pp13TeV", "No analysis manager to connect to.");
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
  ::Info("AddTaskSigmaStar_pp13TeV", "%s", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f", nmix, maxDiffVzMix, maxDiffMultMix));

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
    ::Info("AddTaskSigmaStar_pp13TeV", "%s", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp)?"ON":"OFF"));
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
  Config_Lambdapi(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
  cerr<<"done configuring"<<endl;
  
  // ----- CONTAINERS -----

  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  Printf("AddTaskSigmaStar_pp13TeV - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
  AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",lname.Data()),
							  TList::Class(),
							  AliAnalysisManager::kOutputContainer,
							  outputFileName);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);
   
  return task;
}

Bool_t Config_Lambdapi(
  AliRsnMiniAnalysisTask *task,
  TString     lname,
  Bool_t      isMC,
  Int_t       system,
  Int_t       EventCuts,
  Int_t       TrackCutsLambda,
  Int_t       TrackCutsPi){
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
  if(!cutSetPi){cerr<<"Error in AddTaskSigmaStar_pp13TeV::Config_Lambdapi(): missing cutSetPi"<<endl; return kFALSE;}
  
  Int_t iCutQ=task->AddTrackCuts(cutSetQ);
  Int_t iCutPi=task->AddTrackCuts(cutSetPi);


  // selections for V0 daughters
    
  Int_t V0Cuts=TrackCutsLambda%1000;
  Int_t checkAC=TrackCutsLambda/1000;
    
  Int_t v0d_xrows=70;
  Float_t v0d_rtpc=0.8;
  Float_t v0d_dcaxy=0.06;

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts("qualityDaughterLambda");
  esdTrackCuts->SetPtRange(0.15,1.E10);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetRequireTPCRefit();
  esdTrackCuts->SetAcceptKinkDaughters(0);
  esdTrackCuts->SetMinNCrossedRowsTPC(v0d_xrows);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(v0d_rtpc);
  esdTrackCuts->SetMinDCAToVertexXY(v0d_dcaxy);
    
  Float_t lambda_piPIDCut=3.;
  Float_t lambda_pPIDCut=3.;
  Float_t lambdaDaughDCA=1.6;//0.5;
  Float_t lambdaDCA=0.3;
  Float_t lambda_pLife=30.;
  Float_t lambda_radiuslow=1.4;
  Float_t lambda_radiushigh=100.;//100
  Float_t lambda_massTol=0.010;
  Float_t lambda_massTolVeto=0.004;
  Bool_t  lambdaSwitch=kFALSE;
  Float_t lambdaCosPoinAn=0.99;//0.99

  if(V0Cuts==1) lambdaDCA=1.e10;
  else if(V0Cuts==2) lambdaDaughDCA=0.5;
  else if(V0Cuts==3) lambdaCosPoinAn=0.99;
    
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
  // cutLambda->SetToleranceVeto(lambda_massTolVeto);//Rejection range for Competing V0 Rejection
  cutLambda->SetSwitch(lambdaSwitch);
  cutLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
  cutLambda->SetMaxRapidity(0.8);
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
  // cutAntiLambda->SetToleranceVeto(lambda_massTolVeto);//Rejection range for Competing V0 Rejection
  cutAntiLambda->SetSwitch(lambdaSwitch);
  cutAntiLambda->SetMinCosPointingAngle(lambdaCosPoinAn);
  cutAntiLambda->SetMaxRapidity(0.8);
  cutAntiLambda->SetMinTPCcluster(-1);
    
  AliRsnCutSet* cutSetAntiLambda=new AliRsnCutSet("setAntiLambda",AliRsnTarget::kDaughter);
  cutSetAntiLambda->AddCut(cutAntiLambda);
  cutSetAntiLambda->SetCutScheme(cutAntiLambda->GetName());
  Int_t iCutAntiLambda=task->AddTrackCuts(cutSetAntiLambda);

  // monitoring
  TString pname="lambdap";
  if(enableMonitor){
    Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
#ifdef __CINT__
    gROOT->LoadMacro(
        "$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
#endif
    AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
    AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());
    // AddMonitorOutput_Eta("pi",cutSetPi->GetMonitorOutput());
    // AddMonitorOutput_DCAxy("pi",cutSetPi->GetMonitorOutput());
    // AddMonitorOutput_DCAz("pi",cutSetPi->GetMonitorOutput());
    // AddMonitorOutput_Eta("pi",cutSetPi->GetMonitorOutput());
    // AddMonitorOutput_DCAxy("pi",cutSetPi->GetMonitorOutput());
    // AddMonitorOutput_DCAz("pi",cutSetPi->GetMonitorOutput());
    

    // AddMonitorOutput_P(pname,cutSetLambda->GetMonitorOutput());
    // AddMonitorOutput_Pt(pname,cutSetLambda->GetMonitorOutput());
    // AddMonitorOutput_V0NPt(pname,cutSetLambda->GetMonitorOutput());
    // AddMonitorOutput_V0PPt(pname,cutSetLambda->GetMonitorOutput());
    AddMonitorOutput_V0Mass(pname,cutSetLambda->GetMonitorOutput());
    AddMonitorOutput_V0DCA(pname,cutSetLambda->GetMonitorOutput());
    AddMonitorOutput_V0Radius(pname,cutSetLambda->GetMonitorOutput());
    // AddMonitorOutput_V0Lifetime(pname,cutSetLambda->GetMonitorOutput());
    AddMonitorOutput_V0DaughterDCA(pname,cutSetLambda->GetMonitorOutput());
    AddMonitorOutput_V0CPA(pname,cutSetLambda->GetMonitorOutput());
    AddMonitorOutput_V0DCA2TPV(pname,cutSetLambda->GetMonitorOutput());
    AddMonitorOutput_V0TPCpim(pname,cutSetLambda->GetMonitorOutput());
    AddMonitorOutput_LambdaProtonPID(cutSetLambda->GetMonitorOutput());

    pname.Form("lambdaa");
    // AddMonitorOutput_P(pname,cutSetAntiLambda->GetMonitorOutput());
    // AddMonitorOutput_Pt(pname,cutSetAntiLambda->GetMonitorOutput());
    // AddMonitorOutput_V0NPt(pname,cutSetAntiLambda->GetMonitorOutput());
    // AddMonitorOutput_V0PPt(pname,cutSetAntiLambda->GetMonitorOutput());
    AddMonitorOutput_V0Mass(pname,cutSetAntiLambda->GetMonitorOutput());
    AddMonitorOutput_V0DCA(pname,cutSetAntiLambda->GetMonitorOutput());
    AddMonitorOutput_V0Radius(pname,cutSetAntiLambda->GetMonitorOutput());
    // AddMonitorOutput_V0Lifetime(pname,cutSetAntiLambda->GetMonitorOutput());
    AddMonitorOutput_V0DaughterDCA(pname,cutSetAntiLambda->GetMonitorOutput());
    AddMonitorOutput_V0CPA(pname,cutSetAntiLambda->GetMonitorOutput());
    AddMonitorOutput_V0DCA2TPV(pname,cutSetAntiLambda->GetMonitorOutput());
    AddMonitorOutput_V0TPCpip(pname,cutSetAntiLambda->GetMonitorOutput());
    AddMonitorOutput_LambdaAntiProtonPID(cutSetAntiLambda->GetMonitorOutput());
  }

  // pair cuts
  AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
  if(system!=1) cutY->SetRangeD(-0.5,0.5);
  else cutY->SetRangeD(-0.465,0.035); //pPb

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
                            1        ,  1          ,  1               ,  1             ,              1  ,              1  };
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
    AliRsnMiniOutput* out = task->CreateOutput(
        Form("Lambdapi_%s%s", name[i].Data(), suffix), "HIST", comp[i].Data());
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
    if(checkAC==1){
      out->SetPairCuts(cutsPairMix);
    }else if(checkAC==2){
      out->SetPairCuts(cutsPairSame);
    }else{
      if(!pairID[i]) out->SetPairCuts(cutsPairSame);
      else out->SetPairCuts(cutsPairMix);
    }

    // axis X: invmass or resolution
    if(useIM[i]) out->AddAxis(imID, 1800, 1.2, 3.0);
    else out->AddAxis(diffID,200,-0.02,0.02);
    
    // axis Y: transverse momentum
    out->AddAxis(ptID,200,0.,20.);
    
    // axis Z: centrality-multiplicity
    out->AddAxis(centID,nmult,multbins);
  }

  return kTRUE;
}

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
