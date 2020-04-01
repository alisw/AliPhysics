#ifdef __CLING__
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C>
#endif
/***************************************************************************
 Anders Knospe: anders.knospe@cern.ch
 Macro to configure the resonance package for searches for rare resonances.
 Modified by Corey Myers: corey.james.myers@cern.ch
 Committed by Jihye Song: Jihye.Song@cern.ch (This version is used for test and after testing, it will be removed)
 ****************************************************************************/

/*
 Bool_t Config_pipi(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_pikx(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_pik0(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_kxkx(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_kxk0(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_pkx(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_pk0(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_Lambdapi(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);*/
Bool_t Config_Lambdakx(
                       AliRsnMiniAnalysisTask *task,
                       TString     lname="Lambdakx",
                       Bool_t      isMC=kFALSE,
                       Int_t       system=0,
                       Int_t       EventCuts=0,
                       Int_t       TrackCutsLambda=0,
                       Int_t       TrackCutsK=0
                       );
Bool_t Config_Lambdak0(
                       AliRsnMiniAnalysisTask *task,
                       TString     lname="Lambdak0",
                       Bool_t      isMC=kFALSE,
                       Int_t       system=0,
                       Int_t       EventCuts=0,
                       Int_t       TrackCutsLambda=0,
                       Int_t       TrackCutsK=0
                       );/*
 Bool_t Config_Lambdap(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 */

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

AliRsnMiniAnalysisTask* AddtaskXi1820_test(
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
        ::Error("AddtaskXi1820_test", "No analysis manager to connect to.");
        return NULL;
    }
    
    // create the task and configure
    AliRsnMiniAnalysisTask* task=new AliRsnMiniAnalysisTask(lname,isMC);
    
    // trigger
    int trigger=EventCuts%10;
    if(!trigger) task->UseESDTriggerMask(AliVEvent::kINT7);
    else if(trigger==1) task->UseESDTriggerMask(AliVEvent::kHighMultV0);
    else if(trigger==2) task->UseESDTriggerMask(AliVEvent::kCentral);//PbPb
    else if(trigger==3) task->UseESDTriggerMask(AliVEvent::kSemiCentral);//PbPb
    
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
    float maxDiffMultMix=10;
    task->UseContinuousMix();
    task->SetNMix(nmix);
    task->SetMaxDiffVz(maxDiffVzMix);
    task->SetMaxDiffMult(maxDiffMultMix);
    ::Info("AddtaskXi1820_test", "%s", Form("Event mixing configuration: \n events to mix = %i \n max diff. vtxZ = cm %5.3f \n max diff multi = %5.3f", nmix, maxDiffVzMix, maxDiffMultMix));
    
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
        ::Info("AddtaskXi1820_test", "%s", Form(":::::::::::::::::: Pile-up rejection mode: %s", (rejectPileUp)?"ON":"OFF"));
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
    }else if((!trigger) || (trigger==2) || (trigger==3)){//!trigger
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
    if(d1==AliRsnDaughter::kLambda && d2==AliRsnDaughter::kKaon){
        Config_Lambdakx(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
    }else if(d2==AliRsnDaughter::kLambda && d1==AliRsnDaughter::kKaon){
        Config_Lambdakx(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
        
    }else if(d1==AliRsnDaughter::kLambda && d2==AliRsnDaughter::kKaon0){
        Config_Lambdak0(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
    }else if(d2==AliRsnDaughter::kLambda && d1==AliRsnDaughter::kKaon0){
        Config_Lambdak0(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
    }
        
    cerr<<"done configuring"<<endl;
    
    // ----- CONTAINERS -----
    
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    Printf("AddtaskXi1820_test - Set OutputFileName : \n %s\n", outputFileName.Data() );
    
    AliAnalysisDataContainer *output = mgr->CreateContainer(Form("RsnOut_%s",lname.Data()),
                                                            TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            outputFileName);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, output);
    
    return task;
}

//=============================


Bool_t Config_Lambdakx(AliRsnMiniAnalysisTask* task,
                       TString lname,
                       Bool_t isMC,
                       Int_t system,
                       Int_t EventCuts,
                       Int_t TrackCutsLambda,
                       Int_t TrackCutsK) {
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
    
    // selections for V0 daughters
    //new
    
    //Int_t V0Cuts=TrackCutsLambda%1000;
    Int_t checkAC=0;
    
    Int_t v0d_xrows=70;//70
    Float_t v0d_rtpc=0.8;//0.8
    Float_t v0d_dcaxy=0.06;//0.06
    
    if (((TrackCutsLambda/1000)%10) == 1) {
        v0d_dcaxy=0.06;
    }
    else if(((TrackCutsLambda/1000)%10) == 2){
        v0d_dcaxy=0.07;//tight
    }
    
    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts("qualityDaughterK0s");
    esdTrackCuts->SetEtaRange(-0.8,0.8);
    esdTrackCuts->SetRequireTPCRefit();
    esdTrackCuts->SetAcceptKinkDaughters(0);
    esdTrackCuts->SetMinNCrossedRowsTPC(v0d_xrows);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(v0d_rtpc);
    esdTrackCuts->SetMinDCAToVertexXY(v0d_dcaxy);
    
    // selections for Lambda
    Float_t lambda_piPIDCut=5.;//5.
    Float_t lambda_pPIDCut=5.;//5.
    Float_t lambdaDaughDCA=1.0;//0.5;
    Float_t lambdaDCA=0.4;//1.e10 0.3//0.4
    Float_t lambda_pLife=30.;
    Float_t lambda_radiuslow=0.5;
    Float_t lambda_radiushigh=200.;
    Float_t lambda_massTol=0.006;
    Float_t lambda_massTolVeto=0.004;
    Bool_t  lambdaSwitch=kFALSE;
    Float_t lambdaCosPoinAn=0.99;//0.995 for Lambda analysis
    
    //if(V0Cuts==1) lambdaDCA=1.e10;
    //else if(V0Cuts==2) lambdaDaughDCA=0.5;
    
    if (((TrackCutsLambda/10)%10) == 1) {
        lambdaDCA=0.4;
    }
    else if(((TrackCutsLambda/10)%10) == 2){
        lambdaDCA=0.2;//tight
    }
    
    if ((TrackCutsLambda%10) == 1) {
        lambdaDaughDCA=1.0;
    }
    else if((TrackCutsLambda%10) == 2){
        lambdaDaughDCA=0.3;//tight
    }
    
    if (((TrackCutsLambda/100)%10) == 1) {
        lambdaCosPoinAn=0.99;
    }
    else if(((TrackCutsLambda/100)%10) == 2){
        lambdaCosPoinAn=0.995;//tight
    }
    
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
#ifdef __CINT__
        gROOT->LoadMacro(
                         "$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
#endif
        AddMonitorOutput(isMC, cutSetQ->GetMonitorOutput());
        AddMonitorOutput(isMC, cutSetK->GetMonitorOutput());
        
        
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
    }else if((!trigger) || (trigger==2) || (trigger==3)){
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
    /* mother mass      */ Int_t mmID   = task->CreateValue(AliRsnMiniValue::kInvMassMother,kFALSE);
    /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
    /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
    /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
    /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
    /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);
    /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt,kFALSE);
    /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt,kFALSE);
    
    // -- Create all needed outputs -----------------------------------------------------------------
    
    Int_t xID,cut1,pairID,ipdg;
    TString name,comp;
    Char_t charge2;
    Double_t mass=1.8234;
    AliRsnMiniOutput* out;
    
    for(Int_t i=0;i<20;i++){//0,16
        if(!i){
            name.Form("LambdapKp");
            comp.Form("PAIR");
            charge2='+';
            cut1=iCutLambda;
            ipdg=123314;
            xID=imID;
            pairID=0;
        }else if(i==1){
            name.Form("LambdapKm");
            comp.Form("PAIR");
            charge2='-';
            cut1=iCutLambda;
            ipdg=123314;
            xID=imID;
            pairID=0;
        }else if(i==2){
            name.Form("LambdaaKm");
            comp.Form("PAIR");
            charge2='-';
            cut1=iCutAntiLambda;
            ipdg=-123314;
            xID=imID;
            pairID=0;
        }else if(i==3){
            name.Form("LambdaaKp");
            comp.Form("PAIR");
            charge2='+';
            cut1=iCutAntiLambda;
            ipdg=-123314;
            xID=imID;
            pairID=0;
        }else if(i==4){
            name.Form("LambdapKpMix");
            comp.Form("MIX");
            charge2='+';
            cut1=iCutLambda;
            ipdg=123314;
            xID=imID;
            pairID=1;
        }else if(i==5){
            name.Form("LambdapKmMix");
            comp.Form("MIX");
            charge2='-';
            cut1=iCutLambda;
            ipdg=123314;
            xID=imID;
            pairID=1;
        }else if(i==6){
            name.Form("LambdaaKmMix");
            comp.Form("MIX");
            charge2='-';
            cut1=iCutAntiLambda;
            ipdg=-123314;
            xID=imID;
            pairID=1;
        }else if(i==7){
            name.Form("LambdaaKpMix");
            comp.Form("MIX");
            charge2='+';
            cut1=iCutAntiLambda;
            ipdg=-123314;
            xID=imID;
            pairID=1;
        }else if(i==8){//rotate
            name.Form("LambdapKpROTATE");
            comp.Form("ROTATE1");
            charge2='+';
            cut1=iCutLambda;
            ipdg=123314;
            xID=imID;
            pairID=0;
        }else if(i==9){//rotate
            name.Form("LambdapKmROTATE");
            comp.Form("ROTATE1");
            charge2='-';
            cut1=iCutLambda;
            ipdg=123314;
            xID=imID;
            pairID=0;
        }else if(i==10){//rotate
            name.Form("LambdaaKmROTATE");
            comp.Form("ROTATE1");
            charge2='-';
            cut1=iCutAntiLambda;
            ipdg=-123314;
            xID=imID;
            pairID=0;
        }else if(i==11){//rotate
            name.Form("LambdaaKpROTATE");
            comp.Form("ROTATE1");
            charge2='+';
            cut1=iCutAntiLambda;
            ipdg=-123314;
            xID=imID;
            pairID=0;
        }else if(i==12){//8
            name.Form("Xi1820_m_gen");
            comp.Form("MOTHER");
            charge2='-';
            cut1=iCutLambda;
            ipdg=123314;
            xID=imID;
            pairID=1;
        }else if(i==13){
            name.Form("Xi1820_p_gen");
            comp.Form("MOTHER");
            charge2='+';
            cut1=iCutAntiLambda;
            ipdg=-123314;
            xID=imID;
            pairID=1;
        }else if(i==14){
            name.Form("Xi1820_m_true");
            comp.Form("TRUE");
            charge2='-';
            cut1=iCutLambda;
            ipdg=123314;
            xID=imID;
            pairID=1;
        }else if(i==15){
            name.Form("Xi1820_m_trueMM");
            comp.Form("TRUE");
            charge2='-';
            cut1=iCutLambda;
            ipdg=123314;
            xID=mmID;
            pairID=1;
        }else if(i==16){
            name.Form("Xi1820_m_res");
            comp.Form("TRUE");
            charge2='-';
            cut1=iCutLambda;
            ipdg=123314;
            xID=diffID;
            pairID=1;
        }else if(i==17){
            name.Form("Xi1820_p_true");
            comp.Form("TRUE");
            charge2='+';
            cut1=iCutAntiLambda;
            ipdg=-123314;
            xID=imID;
            pairID=1;
        }else if(i==18){
            name.Form("Xi1820_p_trueMM");
            comp.Form("TRUE");
            charge2='+';
            cut1=iCutAntiLambda;
            ipdg=-123314;
            xID=mmID;
            pairID=1;
        }else if(i==19){//15
            name.Form("Xi1820_p_res");
            comp.Form("TRUE");
            charge2='+';
            cut1=iCutAntiLambda;
            ipdg=-123314;
            xID=diffID;
            pairID=1;
        }
        if(!isMC && i>=12) continue;
        
        out=task->CreateOutput(Form("Lambdakx_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kLambda);
        out->SetCutID(0,cut1);
        out->SetCharge(0,'0');
        
        out->SetDaughter(1,AliRsnDaughter::kKaon);
        out->SetCutID(1,iCutK);
        out->SetCharge(1,charge2);
        
        out->SetMotherPDG(ipdg);
        out->SetMotherMass(mass);
        if(checkAC==1){
            out->SetPairCuts(cutsPairMix);
        }else if(checkAC==2){
            out->SetPairCuts(cutsPairSame);
        }else{
            if(!pairID) out->SetPairCuts(cutsPairSame);
            else out->SetPairCuts(cutsPairMix);
        }
        
        // axis X: invmass or resolution
        if(xID==imID || xID==mmID) out->AddAxis(xID,400,1.6,2.4);
        else out->AddAxis(diffID,200,-0.02,0.02);
        
        // axis Y: transverse momentum
        out->AddAxis(ptID,200,0.,20.);
        
        // axis Z: centrality-multiplicity
        out->AddAxis(centID,nmult,multbins);
    }
    
    if(isMC){//phase-space histograms
        //Xi(1820)-
        out=task->CreateOutput(Form("Xi1820_m_mother_ps%s", suffix),"HIST","TRUE");
        out->SetDaughter(0,AliRsnDaughter::kLambda);
        out->SetCutID(0,iCutLambda);
        out->SetCharge(0,'0');
        
        out->SetDaughter(1,AliRsnDaughter::kKaon);
        out->SetCutID(1,iCutK);
        out->SetCharge(1,'-');
        
        out->SetMotherPDG(123314);
        out->SetMotherMass(mass);
        out->SetPairCuts(cutsPairMix);//just rapidity, no autocorrelation check
        out->AddAxis(fdpt,100,0.,10.);
        out->AddAxis(sdpt,100,0.,10.);
        out->AddAxis(ptID,40,0.,20.);
        
        //anti-Xi(1820)+
        out=task->CreateOutput(Form("Xi1820_p_mother_ps%s", suffix),"HIST","TRUE");
        out->SetDaughter(0,AliRsnDaughter::kLambda);
        out->SetCutID(0,iCutAntiLambda);
        out->SetCharge(0,'0');
        
        out->SetDaughter(1,AliRsnDaughter::kKaon);
        out->SetCutID(1,iCutK);
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


Bool_t Config_Lambdak0(AliRsnMiniAnalysisTask* task,
                       TString lname,
                       Bool_t isMC,
                       Int_t system,
                       Int_t EventCuts,
                       Int_t TrackCutsLambda,
                       Int_t TrackCutsK) {
    bool isPP=false;
    if(!system) isPP=true;
    int trigger=EventCuts%10;
    int MultBins=(EventCuts/10)%10;
    if(system==1 || system==2) MultBins=1;
    
    char suffix[1000];
    sprintf(suffix,"_%s",lname.Data());
    Bool_t enableMonitor=kTRUE;
    
    // selections for V0 daughters
    
    Int_t K0sCuts=((TrackCutsK/1000)%10);
    
    Int_t v0d_xrows=70;//70
    Float_t v0d_rtpc=0.8;//0.8
    Float_t v0d_dcaxy=0.06;//0.06
    
    if (((TrackCutsLambda/1000)%10) == 1) {
        v0d_dcaxy=0.06;
    }
    else if(((TrackCutsLambda/1000)%10) == 2){
        v0d_dcaxy=0.07;//tight
    }
    
    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts("qualityDaughterK0s");
    esdTrackCuts->SetEtaRange(-0.8,0.8);
    esdTrackCuts->SetRequireTPCRefit();
    esdTrackCuts->SetAcceptKinkDaughters(0);
    esdTrackCuts->SetMinNCrossedRowsTPC(v0d_xrows);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(v0d_rtpc);
    esdTrackCuts->SetMinDCAToVertexXY(v0d_dcaxy);
    
    // selections for K0S
    Float_t k0s_piPIDCut=5.0;//5.0
    Float_t k0sDaughDCA=1.;
    Float_t k0sDCA=0.3;//0.3
    Float_t k0s_pLife=20.;
    Float_t k0s_radiuslow=0.5;
    Float_t k0s_radiushigh=200.;
    Float_t k0s_massTolSigma=4;
    Int_t   k0s_massTolID=0;
    Float_t k0s_massTol=0.03;
    Float_t k0s_massTolVeto=0.004;
    Bool_t  k0sSwitch=kFALSE;
    Float_t k0sCosPoinAn=0.97;
    
    if (((TrackCutsK/10)%10) == 1) {
        k0sDCA=0.3;
    }
    else if(((TrackCutsK/10)%10) == 2){
        k0sDCA=0.15;//tight
    }
    
    if ((TrackCutsK%10) == 1) {
        k0sDaughDCA=1.0;
    }
    else if((TrackCutsK%10) == 2){
        k0sDaughDCA=0.3;//tight
    }
    
    if (((TrackCutsK/100)%10) == 1) {
        k0sCosPoinAn=0.97;
    }
    else if(((TrackCutsK/100)%10) == 2){
        k0sCosPoinAn=0.995;//tight
    }
    
    if(K0sCuts==1) k0s_massTolID=1;//use pT-dependent mass tolerance cut
    
    AliRsnCutV0* cutK0s=new AliRsnCutV0("cutK0s",kK0Short,AliPID::kPion,AliPID::kPion);
    cutK0s->SetPIDCutPion(k0s_piPIDCut);// PID for the pion daughters of K0S
    cutK0s->SetESDtrackCuts(esdTrackCuts);// all the other selections (defined above) for pion daughters of K0S
    cutK0s->SetMaxDaughtersDCA(k0sDaughDCA);
    cutK0s->SetMaxDCAVertex(k0sDCA);
    cutK0s->SetfLife(k0s_pLife);
    cutK0s->SetfLowRadius(k0s_radiuslow);
    cutK0s->SetfHighRadius(k0s_radiushigh);
    cutK0s->SetpT_Tolerance(k0s_massTolID);
    cutK0s->SetMassTolSigma(k0s_massTolSigma);
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
    
    //Int_t LambdaCuts=TrackCutsLambda%1000;
    Int_t checkAC=0;
    if((TrackCutsLambda/10000)%10) task->SetCheckDecay(false);
    
    Float_t lambda_piPIDCut=5.0;//5.0
    Float_t lambda_pPIDCut=5.0;//5.0
    Float_t lambdaDaughDCA=1.;//0.5
    Float_t lambdaDCA=0.4;//1.e10 0.3//0.4
    Float_t lambda_pLife=30.;
    Float_t lambda_radiuslow=0.5;
    Float_t lambda_radiushigh=200.;
    Float_t lambda_massTol=0.006;
    Float_t lambda_massTolVeto=0.004;
    Bool_t  lambdaSwitch=kFALSE;
    Float_t lambdaCosPoinAn=0.99;//0.995 for Lambda analysis
    
    if (((TrackCutsLambda/10)%10) == 1) {
        lambdaDCA=0.4;
    }
    else if(((TrackCutsLambda/10)%10) == 2){
        lambdaDCA=0.2;//tight
    }
    
    if ((TrackCutsLambda%10) == 1) {
        lambdaDaughDCA=1.0;
    }
    else if((TrackCutsLambda%10) == 2){
        lambdaDaughDCA=0.3;//tight
    }
    
    if (((TrackCutsLambda/100)%10) == 1) {
        lambdaCosPoinAn=0.99;
    }
    else if(((TrackCutsLambda/100)%10) == 2){
        lambdaCosPoinAn=0.995;//tight
    }
    
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
    }else if((!trigger) || (trigger==2) || (trigger==3)){
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
    /* mother mass      */ Int_t mmID   = task->CreateValue(AliRsnMiniValue::kInvMassMother,kFALSE);
    /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
    /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
    /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
    /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
    /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);
    /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt,kFALSE);
    /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt,kFALSE);
    
    // -- Create all needed outputs -----------------------------------------------------------------
    
    Int_t xID,cut1,pairID,ipdg;
    TString name,comp;
    Double_t mass=1.8234;
    AliRsnMiniOutput* out;
    
    for(Int_t i=0;i<14;i++){
        if(!i){
            name.Form("LambdapK0");
            comp.Form("PAIR");
            cut1=iCutLambda;
            ipdg=123324;
            xID=imID;
            pairID=0;
        }else if(i==1){
            name.Form("LambdaaK0");
            comp.Form("PAIR");
            cut1=iCutAntiLambda;
            ipdg=-123324;
            xID=imID;
            pairID=0;
        }else if(i==2){
            name.Form("LambdapK0Mix");
            comp.Form("MIX");
            cut1=iCutLambda;
            ipdg=123324;
            xID=imID;
            pairID=1;
        }else if(i==3){
            name.Form("LambdaaK0Mix");
            comp.Form("MIX");
            cut1=iCutAntiLambda;
            ipdg=-123324;
            xID=imID;
            pairID=1;
        }else if(i==4){//rotate
            name.Form("LambdapK0ROTATE");
            comp.Form("ROTATE1");
            cut1=iCutLambda;
            ipdg=123324;
            xID=imID;
            pairID=0;
        }else if(i==5){//rotate
            name.Form("LambdaaK0ROTATE");
            comp.Form("ROTATE1");
            cut1=iCutAntiLambda;
            ipdg=-123324;
            xID=imID;
            pairID=0;
        }else if(i==6){//4
            name.Form("Xi1820_0p_gen");
            comp.Form("MOTHER");
            cut1=iCutLambda;
            ipdg=123324;
            xID=imID;
            pairID=1;
        }else if(i==7){
            name.Form("Xi1820_0a_gen");
            comp.Form("MOTHER");
            cut1=iCutAntiLambda;
            ipdg=-123324;
            xID=imID;
            pairID=1;
        }else if(i==8){
            name.Form("Xi1820_0p_true");
            comp.Form("TRUE");
            cut1=iCutLambda;
            ipdg=123324;
            xID=imID;
            pairID=1;
        }else if(i==9){
            name.Form("Xi1820_0p_trueMM");
            comp.Form("TRUE");
            cut1=iCutLambda;
            ipdg=123324;
            xID=mmID;
            pairID=1;
        }else if(i==10){
            name.Form("Xi1820_0p_res");
            comp.Form("TRUE");
            cut1=iCutLambda;
            ipdg=123324;
            xID=diffID;
            pairID=1;
        }else if(i==11){
            name.Form("Xi1820_0a_true");
            comp.Form("TRUE");
            cut1=iCutAntiLambda;
            ipdg=-123324;
            xID=imID;
            pairID=1;
        }else if(i==12){
            name.Form("Xi1820_0a_trueMM");
            comp.Form("TRUE");
            cut1=iCutAntiLambda;
            ipdg=-123324;
            xID=mmID;
            pairID=1;
        }else if(i==13){//11
            name.Form("Xi1820_0a_res");
            comp.Form("TRUE");
            cut1=iCutAntiLambda;
            ipdg=-123324;
            xID=diffID;
            pairID=1;
        }
        if(!isMC && i>=6) continue;
        
        out=task->CreateOutput(Form("Lambdak0_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kLambda);
        out->SetCutID(0,cut1);
        out->SetCharge(0,'0');
        
        out->SetDaughter(1,AliRsnDaughter::kKaon0);
        out->SetCutID(1,iCutK0s);
        out->SetCharge(1,'0');
        
        out->SetMotherPDG(ipdg);
        out->SetMotherMass(mass);
        if(checkAC==1){
            out->SetPairCuts(cutsPairMix);
        }else if(checkAC==2){
            out->SetPairCuts(cutsPairSame);
        }else{
            if(!pairID) out->SetPairCuts(cutsPairSame);
            else out->SetPairCuts(cutsPairMix);
        }
        
        // axis X: invmass or resolution
        if(xID==imID || xID==mmID) out->AddAxis(xID,400,1.6,2.4);
        else out->AddAxis(diffID,200,-0.02,0.02);
        
        // axis Y: transverse momentum
        out->AddAxis(ptID,200,0.,20.);
        
        // axis Z: centrality-multiplicity
        out->AddAxis(centID,nmult,multbins);
    }
    
    if(isMC){ //phase-space histograms
        //Xi(1820)0
        out=task->CreateOutput(Form("Xi1820_0p_mother_ps%s", suffix),"HIST","TRUE");
        out->SetDaughter(0,AliRsnDaughter::kLambda);
        out->SetCutID(0,iCutLambda);
        out->SetCharge(0,'0');
        
        out->SetDaughter(1,AliRsnDaughter::kKaon0);
        out->SetCutID(1,iCutK0s);
        out->SetCharge(1,'0');
        
        out->SetMotherPDG(123324);
        out->SetMotherMass(mass);
        out->SetPairCuts(cutsPairMix);//just rapidity, no autocorrelation check
        out->AddAxis(fdpt,100,0.,10.);
        out->AddAxis(sdpt,100,0.,10.);
        out->AddAxis(ptID,40,0.,20.);
        
        //anti-Xi(1820)0
        out=task->CreateOutput(Form("Xi1820_0a_mother_ps%s", suffix),"HIST","TRUE");
        out->SetDaughter(0,AliRsnDaughter::kLambda);
        out->SetCutID(0,iCutAntiLambda);
        out->SetCharge(0,'0');
        
        out->SetDaughter(1,AliRsnDaughter::kKaon0);
        out->SetCutID(1,iCutK0s);
        out->SetCharge(1,'0');
        
        out->SetMotherPDG(-123324);
        out->SetMotherMass(mass);
        out->SetPairCuts(cutsPairMix);
        out->AddAxis(fdpt,100,0.,10.);
        out->AddAxis(sdpt,100,0.,10.);
        out->AddAxis(ptID,40,0.,20.);
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
