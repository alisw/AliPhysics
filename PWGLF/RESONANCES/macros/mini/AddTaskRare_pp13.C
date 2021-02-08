/***************************************************************************
 Anders Knospe: anders.knospe@cern.ch
 Macro to configure the resonance package for searches for rare resonances.
 ****************************************************************************/

#ifdef __CLING__
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C>
#include <PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C>
#include <PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputCascade.C>

 Bool_t Config_pipi(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_pikx(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_pik0(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_kxkx(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_kxk0(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_pkx(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_pk0(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_Lambdapi(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_Lambdakx(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_Lambdak0(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_Lambdap(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 
 Bool_t Config_LambdaLambda(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_XiPi(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_Xikx(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_Xik0(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_XiP(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
 Bool_t Config_XiLambda(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);

#endif

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
    
    // AODs
    bool isAOD=(system>=100);
    system=system%100;
    
    // trigger
    int trigger=EventCuts%10;
    unsigned int triggerMask=AliVEvent::kINT7;
    if(trigger==1) triggerMask=AliVEvent::kHighMultV0;
    else if(trigger==2) triggerMask=AliVEvent::kMB;
    if(!isAOD) task->UseESDTriggerMask(triggerMask);
    else task->SelectCollisionCandidates(triggerMask);
    
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
    int nmix=10;
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
    
    bool CheckAcceptedMultSelection=true;
    if((EventCuts%100000)/10000==1) CheckAcceptedMultSelection=false;
    
    // other event selection cuts
    AliRsnCutEventUtils* cutEventUtils=0;
    if(1){
        cutEventUtils=new AliRsnCutEventUtils("cutEventUtils",kTRUE,rejectPileUp);
        if(!MultBins){
            cutEventUtils->SetCheckIncompleteDAQ();
            cutEventUtils->SetCheckSPDClusterVsTrackletBG();
        }else{
            cutEventUtils->SetRemovePileUppA2013(kFALSE);
            if(CheckAcceptedMultSelection) cutEventUtils->SetCheckAcceptedMultSelection();
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
    
    //TH2F* hmc=new TH2F("MultiVsCent","", nmult,multbins, 401,ybins);
    //hmc->GetYaxis()->SetTitle("QUALITY");
    TH2F* hmc=new TH2F("TrackletsVsCent","", nmult,multbins, 401,ybins);
    hmc->GetYaxis()->SetTitle("TRACKLETS");
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
        
    }else if(d1==AliRsnDaughter::kLambda && d2==AliRsnDaughter::kLambda){
        Config_LambdaLambda(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
        
    }else if(d1==AliRsnDaughter::kXi && d2==AliRsnDaughter::kPion){
        Config_XiPi(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
    }else if(d2==AliRsnDaughter::kXi && d1==AliRsnDaughter::kPion){
        Config_XiPi(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
        
    }else if(d1==AliRsnDaughter::kXi && d2==AliRsnDaughter::kKaon){
        Config_Xikx(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
    }else if(d2==AliRsnDaughter::kXi && d1==AliRsnDaughter::kKaon){
        Config_Xikx(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
        
    }else if(d1==AliRsnDaughter::kXi && d2==AliRsnDaughter::kKaon0){
        Config_Xik0(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
    }else if(d2==AliRsnDaughter::kXi && d1==AliRsnDaughter::kKaon0){
        Config_Xik0(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
        
    }else if(d1==AliRsnDaughter::kXi && d2==AliRsnDaughter::kProton){
        Config_XiP(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
    }else if(d2==AliRsnDaughter::kXi && d1==AliRsnDaughter::kProton){
        Config_XiP(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
        
    }else if(d1==AliRsnDaughter::kXi && d2==AliRsnDaughter::kLambda){
        Config_XiLambda(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
    }else if(d2==AliRsnDaughter::kXi && d1==AliRsnDaughter::kLambda){
        Config_XiLambda(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
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
    /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
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
        out->SetPairCuts(cutsPair);
        
        // axis X: invmass or resolution
        if(useIM[i]) out->AddAxis(imID,173,0.27,2.);
        else out->AddAxis(diffID,200,-0.02,0.02);
        
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
    /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
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
        else out->AddAxis(diffID,200,-0.02,0.02);
        
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
    
    Int_t V0Cuts=TrackCutsK%1000;
    Int_t checkAC=TrackCutsK/1000;
    
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
    Float_t k0s_massTolSigma=4;
    Int_t   k0s_massTolID=0;
    Float_t k0s_massTol=0.03;
    Float_t k0s_massTolVeto=0.004;
    Bool_t  k0sSwitch=kFALSE;
    Float_t k0sCosPoinAn=0.97;
    
    if(V0Cuts==1) k0s_massTolID=1;//use pT-dependent mass tolerance cut
    
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
    
    // monitoring
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetK0s->GetMonitorOutput(), "K0S");
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
        if(checkAC==1){
            out->SetPairCuts(cutsPairMix);
        }else if(checkAC==2){
            out->SetPairCuts(cutsPairSame);
        }else{
            if(!pairID[i]) out->SetPairCuts(cutsPairSame);
            else out->SetPairCuts(cutsPairMix);
        }
        
        // axis X: invmass or resolution
        if(useIM[i]) out->AddAxis(imID,1370,0.63,2.);
        else out->AddAxis(diffID,200,-0.02,0.02);
        
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
    /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
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
        else out->AddAxis(diffID,200,-0.02,0.02);
        
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
    
    Int_t V0Cuts=TrackCutsK0%1000;
    Int_t checkAC=TrackCutsK0/1000;
    
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
    Float_t k0s_massTolSigma=4;
    Int_t   k0s_massTolID=0;
    Float_t k0s_massTol=0.03;
    Float_t k0s_massTolVeto=0.004;
    Bool_t  k0sSwitch=kFALSE;
    Float_t k0sCosPoinAn=0.97;
    
    if(V0Cuts==1) k0s_massTolID=1;//use pT-dependent mass tolerance cut
    
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
    
    // monitoring
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetKx->GetMonitorOutput());
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetK0s->GetMonitorOutput(), "K0S");
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
        if(checkAC==1){
            out->SetPairCuts(cutsPairMix);
        }else if(checkAC==2){
            out->SetPairCuts(cutsPairSame);
        }else{
            if(i==0 || i==1 || i==4 || i==5) out->SetPairCuts(cutsPairSame);
            else out->SetPairCuts(cutsPairMix);
        }
        
        // axis X: invmass or resolution
        
        if(useIM[i]) out->AddAxis(imID, 1005, 0.99, 3.);
        else out->AddAxis(diffID,200,-0.02,0.02);
        
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
    
    Int_t V0Cuts=TrackCutsK%1000;
    Int_t checkAC=TrackCutsK/1000;
    
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
    Float_t k0s_massTolSigma=4;
    Int_t   k0s_massTolID=0;
    Float_t k0s_massTol=0.03;
    Float_t k0s_massTolVeto=0.004;
    Bool_t  k0sSwitch=kFALSE;
    Float_t k0sCosPoinAn=0.97;
    
    if(V0Cuts==1) k0s_massTolID=1;//use pT-dependent mass tolerance cut
    
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
    
    // monitoring
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetP->GetMonitorOutput());
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetK0s->GetMonitorOutput(), "K0S");
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
    /* IM difference    */ Int_t diffID  = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
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
        if(checkAC==1){
            out->SetPairCuts(cutsPairMix);
        }else if(checkAC==2){
            out->SetPairCuts(cutsPairSame);
        }else{
            if(i==0 || i==1 || i==4 || i==5) out->SetPairCuts(cutsPairSame);
            else out->SetPairCuts(cutsPairMix);
        }
        
        // axis X: invmass or resolution
        if(useIM[i]) out->AddAxis(imID, 785, 1.43, 3.);
        else out->AddAxis(diffID,200,-0.02,0.02);
        
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
    
    Int_t V0Cuts=TrackCutsLambda%1000;
    Int_t checkAC=TrackCutsLambda/1000;
    
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
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetLambda->GetMonitorOutput(), "lambda");
        AddMonitorOutputV0(isMC,cutSetAntiLambda->GetMonitorOutput(), "antilambda");
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
        if(checkAC==1){
            out->SetPairCuts(cutsPairMix);
        }else if(checkAC==2){
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
    
    // selections for V0 daughters
    
    Int_t V0Cuts=TrackCutsLambda%1000;
    Int_t checkAC=TrackCutsLambda/1000;
    
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
    cutLambda->SetMaxPseudorapidity(0.8);
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
    cutAntiLambda->SetMaxPseudorapidity(0.8);
    cutAntiLambda->SetMinTPCcluster(-1);
    
    AliRsnCutSet* cutSetAntiLambda=new AliRsnCutSet("setAntiLambda",AliRsnTarget::kDaughter);
    cutSetAntiLambda->AddCut(cutAntiLambda);
    cutSetAntiLambda->SetCutScheme(cutAntiLambda->GetName());
    Int_t iCutAntiLambda=task->AddTrackCuts(cutSetAntiLambda);
    
    // monitoring
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetK->GetMonitorOutput());
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetLambda->GetMonitorOutput(), "lambda");
        AddMonitorOutputV0(isMC,cutSetAntiLambda->GetMonitorOutput(), "antilambda");
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
    /* mother mass      */ Int_t mmID   = task->CreateValue(AliRsnMiniValue::kInvMassMother,kFALSE);
    /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
    /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
    /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
    /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
    /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);
    /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt,kFALSE);
    /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt,kFALSE);
    
    // -- Create all needed outputs -----------------------------------------------------------------
    
    Int_t i,k,xID,cut1,pairID,ipdg;
    TString name,comp;
    Char_t charge2;
    Double_t mass=1.8234;
    AliRsnMiniOutput* out;
    
    task->SetMotherAcceptanceCutMinPt(0.15);
    task->SetMotherAcceptanceCutMaxEta(0.8);
    task->KeepMotherInAcceptance(true);
    
    for(i=0;i<2;i++) for(j=0;j<2;j++) for(k=0;k<9;k++){
        ipdg=0;
        if(!i){
            name.Form("Lambdap");
            cut1=iCutLambda;
            ipdg=123314;
        }else{
            name.Form("Lambdaa");
            cut1=iCutAntiLambda;
            ipdg=-123314;
        }
        
        if(!j){
            name.Append("Kp");
            charge2='+';
        }else{
            name.Append("Km");
            charge2='-';
        }
        
        if(!isMC && k>=3) continue;
        xID=imID;
        pairID=1;
        if(!k){
            comp.Form("PAIR");
            pairID=0;
        }else if(k==1){
            name.Append("Mix");
            comp.Form("MIX");
        }else if(k==2){
            name.Append("Rotated");
            comp.Form("ROTATE1");
            pairID=0;
        }else if(k==3){
            name.Append("_gen");
            comp.Form("MOTHER");
        }else if(k==4){
            name.Append("_rec");
            comp.Form("TRUE");
        }else if(k==5){
            name.Append("_recMM");
            comp.Form("TRUE");
            xID=mmID;
        }else if(k==6){
            name.Append("_res");
            comp.Form("TRUE");
            xID=diffID;
        }else if(k==7){
            name.Append("_genPS");
            comp.Form("MOTHER_IN_ACC");
        }else if(k==8){
            name.Append("_recPS");
            comp.Form("TRUE");
        }
        
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
        
        if(k<=6){
            if(xID==imID || xID==mmID) out->AddAxis(xID,400,1.6,2.4);// axis X: invmass or resolution
            else out->AddAxis(diffID,200,-0.02,0.02);
            out->AddAxis(ptID,200,0.,20.);// axis Y: transverse momentum
            out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
        }else{//Phase-space histograms
            out->AddAxis(fdpt,100,0.,10.);
            out->AddAxis(sdpt,100,0.,10.);
            out->AddAxis(ptID,40,0.,20.);
        }
    }
    
    if(isMC){
        //for efficiency of K+/-
        AliRsnCutMiniPair* cutetaK=new AliRsnCutMiniPair("cutPseudorapidityK", AliRsnCutMiniPair::kPseudorapidityRangeMC);
        cutetaK->SetRangeD(-0.8,0.8);
        
        AliRsnCutSet* cutsPairK=new AliRsnCutSet("pairCutsK", AliRsnTarget::kMother);
        cutsPairK->AddCut(cutetaK);
        cutsPairK->SetCutScheme(cutetaK->GetName());
        
        out = task->CreateOutput("Lambdakx_Kp_mother", "HIST", "SINGLE");
        out->SetDaughter(0, AliRsnDaughter::kUnknown);
        out->SetDaughter(1, AliRsnDaughter::kUnknown);
        out->SetMotherPDG(321);
        out->SetMotherMass(0.493677);
        out->SetPairCuts(cutsPairK);
        out->AddAxis(ptID,200,0.0,20.0);
        
        out = task->CreateOutput("Lambdakx_Km_mother", "HIST", "SINGLE");
        out->SetDaughter(0, AliRsnDaughter::kUnknown);
        out->SetDaughter(1, AliRsnDaughter::kUnknown);
        out->SetMotherPDG(-321);
        out->SetMotherMass(0.493677);
        out->SetPairCuts(cutsPairK);
        out->AddAxis(ptID,200,0.0,20.0);
        
        //for efficiency of (anti-)Lambda
        AliRsnCutMiniPair* cutetaV0=new AliRsnCutMiniPair("cutPseudorapidityV0", AliRsnCutMiniPair::kPseudorapidityRangeMC);
        cutetaV0->SetRangeD(-0.8,0.8);
        
        AliRsnCutSet* cutsPairV0=new AliRsnCutSet("pairCutsV0", AliRsnTarget::kMother);
        cutsPairV0->AddCut(cutetaV0);
        cutsPairV0->SetCutScheme(cutetaV0->GetName());
        
        out = task->CreateOutput("Lambdakx_Lambdap_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kProton);
        out->SetMotherPDG(3122);
        out->SetMotherMass(1.115683);
        out->SetPairCuts(cutsPairV0);
        out->AddAxis(ptID,200,0.0,20.0);
        
        out = task->CreateOutput("Lambdakx_Lambdaa_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kProton);
        out->SetMotherPDG(-3122);
        out->SetMotherMass(1.115683);
        out->SetPairCuts(cutsPairV0);
        out->AddAxis(ptID,200,0.0,20.0);
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
    
    // selections for V0 daughters

    Int_t K0sCuts=TrackCutsK%1000;
    Bool_t CheckOOBP = true;
    Bool_t CheckOOBPv0 = false;
    if(TrackCutsK>1000){
        CheckOOBP=false;
        CheckOOBPv0=true;
    }
    
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
    Float_t k0s_massTolSigma=4;
    Int_t   k0s_massTolID=0;
    Float_t k0s_massTol=0.03;
    Float_t k0s_massTolVeto=0.004;
    Bool_t  k0sSwitch=kFALSE;
    Float_t k0sCosPoinAn=0.97;
    
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
    cutK0s->SetMaxPseudorapidity(0.8);
    cutK0s->SetMinTPCcluster(-1);
    cutK0s->SetCheckOOBPileup(CheckOOBPv0);
    
    AliRsnCutSet* cutSetK0s=new AliRsnCutSet("setK0s",AliRsnTarget::kDaughter);
    cutSetK0s->AddCut(cutK0s);
    cutSetK0s->SetCutScheme(cutK0s->GetName());
    Int_t iCutK0s=task->AddTrackCuts(cutSetK0s);
    
    // selections for Lambda
    
    Int_t LambdaCuts=TrackCutsLambda%1000;
    Int_t checkAC=(TrackCutsLambda/1000)%10;
    if((TrackCutsLambda/10000)%10) task->SetCheckDecay(false);
    
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
    cutLambda->SetMaxPseudorapidity(0.8);
    cutLambda->SetMinTPCcluster(-1);
    cutLambda->SetCheckOOBPileup(CheckOOBPv0);
    
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
    cutAntiLambda->SetMaxPseudorapidity(0.8);
    cutAntiLambda->SetMinTPCcluster(-1);
    cutAntiLambda->SetCheckOOBPileup(CheckOOBPv0);
    
    AliRsnCutSet* cutSetAntiLambda=new AliRsnCutSet("setAntiLambda",AliRsnTarget::kDaughter);
    cutSetAntiLambda->AddCut(cutAntiLambda);
    cutSetAntiLambda->SetCutScheme(cutAntiLambda->GetName());
    Int_t iCutAntiLambda=task->AddTrackCuts(cutSetAntiLambda);
    
    // monitoring
    if(enableMonitor){
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetK0s->GetMonitorOutput(), "K0S");
        AddMonitorOutputV0(isMC,cutSetLambda->GetMonitorOutput(), "lambda");
        AddMonitorOutputV0(isMC,cutSetAntiLambda->GetMonitorOutput(), "antilambda");
    }
    
    // pair cuts
    AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
    if(system!=1) cutY->SetRangeD(-0.5,0.5);
    else cutY->SetRangeD(-0.465,0.035);
    
    AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
    AliRsnCutMiniPair* cutOOBP=new AliRsnCutMiniPair("cutOOBP",AliRsnCutMiniPair::kPassesOOBPileupCut);
    
    AliRsnCutSet* cutsPairSame=new AliRsnCutSet("pairCutsSame",AliRsnTarget::kMother);
    cutsPairSame->AddCut(cutY);
    cutsPairSame->AddCut(cutV0);
    
    AliRsnCutSet* cutsPairMix=new AliRsnCutSet("pairCutsMix", AliRsnTarget::kMother);
    cutsPairMix->AddCut(cutY);
    
    if(CheckOOBP){
        cutsPairSame->AddCut(cutOOBP);
        cutsPairSame->SetCutScheme(TString::Format("%s&(!%s)&%s",cutY->GetName(),cutV0->GetName(),cutOOBP->GetName()).Data());
        cutsPairMix->AddCut(cutOOBP);
        cutsPairMix->SetCutScheme(TString::Format("%s&%s",cutY->GetName(),cutOOBP->GetName()).Data());
    }else{
        cutsPairSame->SetCutScheme(TString::Format("%s&(!%s)",cutY->GetName(),cutV0->GetName()).Data());
        cutsPairMix->SetCutScheme(cutY->GetName());
    }
    
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
    /* mother mass      */ Int_t mmID   = task->CreateValue(AliRsnMiniValue::kInvMassMother,kFALSE);
    /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
    /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
    /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
    /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
    /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);
    /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt,kFALSE);
    /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt,kFALSE);
    
    // -- Create all needed outputs -----------------------------------------------------------------
    
    Int_t i,xID,cut1,pairID,ipdg;
    TString name,comp;
    Double_t mass=1.8234;
    AliRsnMiniOutput* out;
    
    task->SetMotherAcceptanceCutMinPt(0.15);
    task->SetMotherAcceptanceCutMaxEta(0.8);
    task->KeepMotherInAcceptance(true);
    
    for(i=0;i<2;i++) for(j=0;j<9;j++){
        if(!i){
            name.Form("LambdapK0");
            cut1=iCutLambda;
            ipdg=123324;
        }else{
            name.Form("LambdaaK0");
            cut1=iCutAntiLambda;
            ipdg=-123324;
        }
        
        if(!isMC && j>=3) continue;
        xID=imID;
        pairID=1;
        if(!j){
            comp.Form("PAIR");
            pairID=0;
        }else if(j==1){
            name.Append("Mix");
            comp.Form("MIX");
        }else if(j==2){
            name.Append("Rotated");
            comp.Form("ROTATE1");
            pairID=0;
        }else if(j==3){
            name.Append("_gen");
            comp.Form("MOTHER");
        }else if(j==4){
            name.Append("_rec");
            comp.Form("TRUE");
        }else if(j==5){
            name.Append("_recMM");
            comp.Form("TRUE");
            xID=mmID;
        }else if(j==6){
            name.Append("_res");
            comp.Form("TRUE");
            xID=diffID;
        }else if(j==7){
            name.Append("_genPS");
            comp.Form("MOTHER_IN_ACC");
        }else if(j==8){
            name.Append("_recPS");
            comp.Form("TRUE");
        }
        
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
        
        if(j<=6){
            if(xID==imID || xID==mmID) out->AddAxis(xID,400,1.6,2.4);// axis X: invmass or resolution
            else out->AddAxis(diffID,200,-0.02,0.02);
            out->AddAxis(ptID,200,0.,20.);// axis Y: transverse momentum
            out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
        }else{//Phase-space histograms
            out->AddAxis(fdpt,100,0.,10.);
            out->AddAxis(sdpt,100,0.,10.);
            out->AddAxis(ptID,40,0.,20.);
        }
    }
    
    if(isMC){
        //for efficiency of (anti-)Lambda
        AliRsnCutMiniPair* cutYV0=new AliRsnCutMiniPair("cutRapidityV0", AliRsnCutMiniPair::kPseudorapidityRangeMC);
        cutYV0->SetRangeD(-0.8,0.8);
        
        AliRsnCutSet* cutsPairV0=new AliRsnCutSet("pairCutsV0", AliRsnTarget::kMother);
        cutsPairV0->AddCut(cutYV0);
        cutsPairV0->SetCutScheme(cutYV0->GetName());
        
        out = task->CreateOutput("Lambdak0_Lambdap_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kProton);
        out->SetMotherPDG(3122);
        out->SetMotherMass(1.115683);
        out->SetPairCuts(cutsPairV0);
        out->AddAxis(ptID,200,0.0,20.0);
        
        out = task->CreateOutput("Lambdak0_Lambdaa_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kProton);
        out->SetMotherPDG(-3122);
        out->SetMotherMass(1.115683);
        out->SetPairCuts(cutsPairV0);
        out->AddAxis(ptID,200,0.0,20.0);
        
        //for efficiency of K0S
        out = task->CreateOutput("Lambdak0_K0S_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kPion);
        out->SetMotherPDG(310);
        out->SetMotherMass(0.497611);
        out->SetPairCuts(cutsPairV0);
        out->AddAxis(ptID,200,0.0,20.0);
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
    
    Int_t V0Cuts=TrackCutsLambda%1000;
    Int_t checkAC=TrackCutsLambda/1000;
    
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
    cutLambda->SetMaxPseudorapidity(0.8);
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
    cutAntiLambda->SetMaxPseudorapidity(0.8);
    cutAntiLambda->SetMinTPCcluster(-1);
    
    AliRsnCutSet* cutSetAntiLambda=new AliRsnCutSet("setAntiLambda",AliRsnTarget::kDaughter);
    cutSetAntiLambda->AddCut(cutAntiLambda);
    cutSetAntiLambda->SetCutScheme(cutAntiLambda->GetName());
    Int_t iCutAntiLambda=task->AddTrackCuts(cutSetAntiLambda);
    
    // monitoring
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetP->GetMonitorOutput());
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetLambda->GetMonitorOutput(), "lambda");
        AddMonitorOutputV0(isMC,cutSetAntiLambda->GetMonitorOutput(), "antilambda");
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
    /* mother mass      */ Int_t mmID   = task->CreateValue(AliRsnMiniValue::kInvMassMother,kFALSE);
    /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
    /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
    /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
    /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
    /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY, kFALSE);
    /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt,kFALSE);
    /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt,kFALSE);
    
    // -- Create all needed outputs -----------------------------------------------------------------
    // use an array for more compact writing, which are different on mixing and charges
    
    Int_t i,k,xID,cut1,pairID,ipdg;
    char charge2;
    TString name,comp;
    AliRsnMiniOutput* out;
    Double_t mass=1.115683+0.938272;
    
    task->SetMotherAcceptanceCutMinPt(0.15);
    task->SetMotherAcceptanceCutMaxEta(0.8);
    task->KeepMotherInAcceptance(true);
    
    for(i=0;i<2;i++) for(j=0;j<2;j++) for(k=0;k<9;k++){
        if(!i){
            name.Form("Lambdap");
            cut1=iCutLambda;
            ipdg=1;
        }else{
            name.Form("Lambdaa");
            cut1=iCutAntiLambda;
            ipdg=-1;
        }
        
        if(!j){
            name.Append("Pp");
            charge2='+';
        }else{
            name.Append("Pm");
            charge2='-';
        }
        
        if(!isMC && k>=3) continue;
        xID=imID;
        pairID=1;
        if(!k){
            comp.Form("PAIR");
            pairID=0;
        }else if(k==1){
            name.Append("Mix");
            comp.Form("MIX");
        }else if(k==2){
            name.Append("Rotated");
            comp.Form("ROTATE1");
            pairID=0;
        }else if(k==3){
            name.Append("_gen");
            comp.Form("MOTHER");
        }else if(k==4){
            name.Append("_rec");
            comp.Form("TRUE");
        }else if(k==5){
            name.Append("_recMM");
            comp.Form("TRUE");
            xID=mmID;
        }else if(k==6){
            name.Append("_res");
            comp.Form("TRUE");
            xID=diffID;
        }else if(k==7){
            name.Append("_genPS");
            comp.Form("MOTHER_IN_ACC");
        }else if(k==8){
            name.Append("_recPS");
            comp.Form("TRUE");
        }
        
        out=task->CreateOutput(Form("Lambdap_%s%s",name.Data(),suffix),"HIST",comp.Data());
        // selection settings
        out->SetDaughter(0,AliRsnDaughter::kLambda);
        out->SetCharge(0,'0');
        out->SetCutID(0,cut1);
        
        out->SetDaughter(1,AliRsnDaughter::kProton);
        out->SetCharge(1,charge2);
        out->SetCutID(1,iCutP);
        
        out->SetMotherPDG(ipdg);
        out->SetMotherMass(mass);
        // pair cuts
        if(!pairID) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        
        if(k<=6){
            if(xID==imID || xID==mmID) out->AddAxis(xID,975,2.05,4.);
            else out->AddAxis(diffID,200,-0.02,0.02);
            out->AddAxis(ptID,200,0.,20.);
            out->AddAxis(centID,nmult,multbins);
        }else{
            out->AddAxis(fdpt,100,0.,10.);
            out->AddAxis(sdpt,100,0.,10.);
            out->AddAxis(ptID,40,0.,20.);
        }
    }
    
    if(isMC){
        //for efficiency of (anti)proton
        AliRsnCutMiniPair* cutetaP=new AliRsnCutMiniPair("cutPseudorapidityP", AliRsnCutMiniPair::kPseudorapidityRangeMC);
        cutetaP->SetRangeD(-0.8,0.8);
        
        AliRsnCutSet* cutsPairP=new AliRsnCutSet("pairCutsP", AliRsnTarget::kMother);
        cutsPairP->AddCut(cutetaP);
        cutsPairP->SetCutScheme(cutetaP->GetName());
        
        out = task->CreateOutput("Lambdap_Pp_mother", "HIST", "SINGLE");
        out->SetDaughter(0, AliRsnDaughter::kUnknown);
        out->SetDaughter(1, AliRsnDaughter::kUnknown);
        out->SetMotherPDG(2212);
        out->SetMotherMass(0.938272);
        out->SetPairCuts(cutsPairP);
        out->AddAxis(ptID,200,0.0,20.0);
        
        out = task->CreateOutput("Lambdap_Pm_mother", "HIST", "SINGLE");
        out->SetDaughter(0, AliRsnDaughter::kUnknown);
        out->SetDaughter(1, AliRsnDaughter::kUnknown);
        out->SetMotherPDG(-2212);
        out->SetMotherMass(0.938272);
        out->SetPairCuts(cutsPairP);
        out->AddAxis(ptID,200,0.0,20.0);
        
        //for efficiency of (anti-)Lambda
        AliRsnCutMiniPair* cutetaV0=new AliRsnCutMiniPair("cutPseudorapidityV0", AliRsnCutMiniPair::kPseudorapidityRangeMC);
        cutetaV0->SetRangeD(-0.8,0.8);
        
        AliRsnCutSet* cutsPairV0=new AliRsnCutSet("pairCutsV0", AliRsnTarget::kMother);
        cutsPairV0->AddCut(cutetaV0);
        cutsPairV0->SetCutScheme(cutetaV0->GetName());
        
        out = task->CreateOutput("Lambdap_Lambdap_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kProton);
        out->SetMotherPDG(3122);
        out->SetMotherMass(1.115683);
        out->SetPairCuts(cutsPairV0);
        out->AddAxis(ptID,200,0.0,20.0);
        
        out = task->CreateOutput("Lambdap_Lambdaa_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kProton);
        out->SetMotherPDG(-3122);
        out->SetMotherMass(1.115683);
        out->SetPairCuts(cutsPairV0);
        out->AddAxis(ptID,200,0.0,20.0);
        
    }
    
    return kTRUE;
}


//=============================


Bool_t Config_LambdaLambda(
                           AliRsnMiniAnalysisTask *task,
                           TString     lname,
                           Bool_t      isMC,
                           Int_t       system,
                           Int_t       EventCuts,
                           Int_t       TrackCutsLambda,
                           Int_t       TrackCuts2
                           ){
    bool isPP=false;
    if(!system) isPP=true;
    int trigger=EventCuts%10;
    int MultBins=(EventCuts/10)%10;
    if(system==1 || system==2) MultBins=1;
    
    char suffix[1000];
    sprintf(suffix,"_%s",lname.Data());
    Bool_t enableMonitor=kTRUE;
    
    Double_t mass= 1.115683+1.115683;
    
    //Not sure about CutQ, but better to have it in case
    AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
    trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
    
    AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,
                                                                           AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);
    Int_t iCutQ=task->AddTrackCuts(cutSetQ);
    
    Bool_t CheckOOBP=true;
    if(TrackCuts2==1) CheckOOBP=false;
    
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
    cutLambda->SetMaxPseudorapidity(0.8);
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
    cutAntiLambda->SetMaxPseudorapidity(0.8);
    cutAntiLambda->SetMinTPCcluster(-1);
    
    AliRsnCutSet* cutSetAntiLambda=new AliRsnCutSet("setAntiLambda",AliRsnTarget::kDaughter);
    cutSetAntiLambda->AddCut(cutAntiLambda);
    cutSetAntiLambda->SetCutScheme(cutAntiLambda->GetName());
    Int_t iCutAntiLambda=task->AddTrackCuts(cutSetAntiLambda);
    
    // monitoring
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetLambda->GetMonitorOutput(), "lambda");
        AddMonitorOutputV0(isMC,cutSetAntiLambda->GetMonitorOutput(), "antilambda");
    }
    
    // pair cuts
    AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
    if(system!=1) cutY->SetRangeD(-0.5,0.5);
    else cutY->SetRangeD(-0.465,0.035);
    
    AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
    AliRsnCutMiniPair* cutOOBP=new AliRsnCutMiniPair("cutOOBP",AliRsnCutMiniPair::kPassesOOBPileupCut);
    
    AliRsnCutSet* cutsPairSame=new AliRsnCutSet("pairCutsSame",AliRsnTarget::kMother);
    cutsPairSame->AddCut(cutY);
    cutsPairSame->AddCut(cutV0);
    
    AliRsnCutSet* cutsPairMix=new AliRsnCutSet("pairCutsMix", AliRsnTarget::kMother);
    cutsPairMix->AddCut(cutY);
    
    if(CheckOOBP){
        cutsPairSame->AddCut(cutOOBP);
        cutsPairSame->SetCutScheme(TString::Format("%s&(!%s)&%s",cutY->GetName(),cutV0->GetName(),cutOOBP->GetName()).Data());
        cutsPairMix->AddCut(cutOOBP);
        cutsPairMix->SetCutScheme(TString::Format("%s&%s",cutY->GetName(),cutOOBP->GetName()).Data());
    }else{
        cutsPairSame->SetCutScheme(TString::Format("%s&(!%s)",cutY->GetName(),cutV0->GetName()).Data());
        cutsPairMix->SetCutScheme(cutY->GetName());
    }
    
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
    /* mother mass      */ Int_t mmID   = task->CreateValue(AliRsnMiniValue::kInvMassMother,kFALSE);
    /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
    /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
    /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
    /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
    /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);
    /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt,kFALSE);
    /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt,kFALSE);
    
    // -- Create all needed outputs -----------------------------------------------------------------
    // use an array for more compact writing, which are different on mixing and charges
    
    Int_t i,xID,cut1,cut2,pairID,ipdg;
    TString name,comp;
    AliRsnMiniOutput* out;
    
    task->SetMotherAcceptanceCutMinPt(0.15);
    task->SetMotherAcceptanceCutMaxEta(0.8);
    task->KeepMotherInAcceptance(true);
    
    for(i=0;i<3;i++) for(j=0;j<9;j++){
        if(!i){
            name.Form("LambdapLambdap");
            comp.Form("PAIR");
            cut1=cut2=iCutLambda;
        }else if(i==1){
            name.Form("LambdapLambdaa");
            comp.Form("PAIR");
            cut1=iCutLambda;
            cut2=iCutAntiLambda;
        }else{
            name.Form("LambdaaLambdaa");
            comp.Form("PAIR");
            cut1=cut2=iCutAntiLambda;
        }
        
        if(!isMC && j>=3) continue;
        xID=imID;
        pairID=1;
        if(!j){
            comp.Form("PAIR");
            pairID=0;
        }else if(j==1){
            name.Append("Mix");
            comp.Form("MIX");
        }else if(j==2){
            name.Append("Rotated");
            comp.Form("ROTATE1");
            pairID=0;
        }else if(j==3){
            name.Append("_gen");
            comp.Form("MOTHER");
        }else if(j==4){
            name.Append("_rec");
            comp.Form("TRUE");
        }else if(j==5){
            name.Append("_recMM");
            comp.Form("TRUE");
            xID=mmID;
        }else if(j==6){
            name.Append("_res");
            comp.Form("TRUE");
            xID=diffID;
        }else if(j==7){
            name.Append("_genPS");
            comp.Form("MOTHER_IN_ACC");
        }else if(j==8){
            name.Append("_recPS");
            comp.Form("TRUE");
        }
        
        out=task->CreateOutput(Form("LambdaLambda_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kLambda);
        out->SetCutID(0,cut1);
        out->SetCharge(0,'0');
        
        out->SetDaughter(1,AliRsnDaughter::kLambda);
        out->SetCutID(1,cut2);
        out->SetCharge(1,'0');
        
        if(!pairID) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(ipdg);
        out->SetMotherMass(mass);
        if(name.Contains("LambdapLambdaa")) out->SetCheckSameCutID(true);
        
        if(j<=6){
            if(xID==imID || xID==mmID) out->AddAxis(xID,200,2.2,3.2);// axis X: invmass or resolution
            else out->AddAxis(diffID,200,-0.02,0.02);
            out->AddAxis(ptID,200,0.0,20.0);// axis Y: transverse momentum
            out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
        }else{//Phase-space histograms
            out->AddAxis(fdpt,100,0.,10.);
            out->AddAxis(sdpt,100,0.,10.);
            out->AddAxis(ptID,40,0.,20.);
        }
    }
    
    if(isMC){
        //for efficiency of (anti-)Lambda
        AliRsnCutMiniPair* cutetaV0=new AliRsnCutMiniPair("cutPseudorapidityV0", AliRsnCutMiniPair::kPseudorapidityRangeMC);
        cutetaV0->SetRangeD(-0.8,0.8);
        
        AliRsnCutSet* cutsPairV0=new AliRsnCutSet("pairCutsV0", AliRsnTarget::kMother);
        cutsPairV0->AddCut(cutetaV0);
        cutsPairV0->SetCutScheme(cutetaV0->GetName());
        
        out = task->CreateOutput("Lambdap_Lambdap_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kProton);
        out->SetMotherPDG(3122);
        out->SetMotherMass(1.115683);
        out->SetPairCuts(cutsPairV0);
        out->AddAxis(ptID,200,0.0,20.0);
        
        out = task->CreateOutput("Lambdap_Lambdaa_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kProton);
        out->SetMotherPDG(-3122);
        out->SetMotherMass(1.115683);
        out->SetPairCuts(cutsPairV0);
        out->AddAxis(ptID,200,0.0,20.0);
    }
    
    return kTRUE;
}


//=============================


Bool_t Config_XiPi(
                   AliRsnMiniAnalysisTask *task,
                   TString     lname,
                   Bool_t      isMC,
                   Int_t       system,
                   Int_t       EventCuts,
                   Int_t       TrackCutsXi,
                   Int_t       TrackCutsPi
                   ){
    bool isPP=false;
    if(!system) isPP=true;
    int trigger=EventCuts%10;
    int MultBins=(EventCuts/10)%10;
    if(system==1 || system==2) MultBins=1;
    
    char suffix[1000];
    sprintf(suffix,"_%s",lname.Data());
    Bool_t enableMonitor=kTRUE;
    
    Double_t mass=0.139571+1.67243;//1.32171;
    
    // set cuts for pions
    if(!(TrackCutsPi%10000)) TrackCutsPi+=3020;//default settings
    Float_t nsigmaPiTPC=0.1*(TrackCutsPi%100);
    Float_t nsigmaPiTOF=0.1*((TrackCutsPi/100)%100);
    Int_t CutTypePi=(TrackCutsPi/10000)%100;//0=TPC+TOFveto (default), 1=TPC only, 2=TOF only, 3 TPC+TOFcut
    
    AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
    trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
    
    AliRsnCutSetDaughterParticle* cutSetPi=0;
    if(!CutTypePi) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaPiTPC),
                                                             trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kPion,nsigmaPiTPC,nsigmaPiTOF);
    else if(CutTypePi==1) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,nsigmaPiTPC),
                                                                    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kPion,nsigmaPiTPC,-1.);
    else if(CutTypePi==2) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,nsigmaPiTOF),
                                                                    trkQualityCut,AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,AliPID::kPion,-1.,nsigmaPiTOF);
    else if(CutTypePi==3) cutSetPi=new AliRsnCutSetDaughterParticle(Form("cutPi%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFtightPidKStarPPB2011,nsigmaPiTPC),
                                                                    trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFtightPidKStarPPB2011,AliPID::kPion,nsigmaPiTPC,-1.);
    if(!cutSetPi){cerr<<"Error in AddTaskResonanceFinder::Config_XiPi(): missing cutSetPi"<<endl; return kFALSE;}
    
    Int_t iCutPi=task->AddTrackCuts(cutSetPi);
    
    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts("qualityTrackcut");
    esdTrackCuts->SetEtaRange(-0.8,0.8);
    esdTrackCuts->SetRequireTPCRefit();
    esdTrackCuts->SetAcceptKinkDaughters(0);
    esdTrackCuts->SetMinNCrossedRowsTPC(70);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetMinDCAToVertexXY(0.06);
    
    // selections for Xi
    Float_t XiPIDcut=3.;
    Float_t V0dDCA=1.6;
    Float_t XidDCA=1.6;
    Float_t XiMinDCA=0.07;
    Float_t Xi_massTol=0.007;
    Float_t Xi_massTolVeto=0.007;
    Float_t V0CosPoinAn=0.97;
    Float_t XiCosPoinAn=0.97;
    
    AliRsnCutCascade* cutXi=new AliRsnCutCascade("cutXi",kOmegaMinus);//kXiMinus);
    cutXi->SetPIDCutV0Proton(XiPIDcut);
    cutXi->SetPIDCutV0Pion(XiPIDcut);
    cutXi->SetPIDCutBachelor(XiPIDcut);
    cutXi->SetESDtrackCuts(esdTrackCuts);
    cutXi->SetV0MaxDaughtersDCA(V0dDCA);
    cutXi->SetCascadeMaxDaughtersDCA(XidDCA);
    cutXi->SetV0MaxDCAVertex(1e5); // not using
    cutXi->SetV0MinDCAVertex(XiMinDCA);
    cutXi->SetCascadeMaxDCAVertex(1e5); // not using
    cutXi->SetCascadeMinDCAVertex(-1e5); // not using
    cutXi->SetV0LowRadius(0); // not using
    cutXi->SetV0HighRadius(1e5); // not using
    cutXi->SetCascadeLowRadius(0); // not using
    cutXi->SetCascadeHighRadius(1e5); // not using
    cutXi->SetMassTolerance(Xi_massTol);
    cutXi->SetMassToleranceVeto(Xi_massTolVeto);//Rejection range for Competing Xi Rejection
    cutXi->SetSwitch(kFALSE); // not using
    cutXi->SetV0MinCosPointingAngle(V0CosPoinAn);
    cutXi->SetCascadeMinCosPointingAngle(XiCosPoinAn);
    cutXi->SetMaxPseudorapidity(0.8);
    cutXi->SetMinTPCcluster(-1);
    
    AliRsnCutCascade* cutXibar=new AliRsnCutCascade("cutXibar",kOmegaPlusBar);//kXiPlusBar);
    cutXibar->SetPIDCutV0Proton(XiPIDcut);
    cutXibar->SetPIDCutV0Pion(XiPIDcut);
    cutXibar->SetPIDCutBachelor(XiPIDcut);
    cutXibar->SetESDtrackCuts(esdTrackCuts);
    cutXibar->SetV0MaxDaughtersDCA(V0dDCA);
    cutXibar->SetCascadeMaxDaughtersDCA(XidDCA);
    cutXibar->SetV0MaxDCAVertex(1e5); // not using
    cutXibar->SetV0MinDCAVertex(XiMinDCA);
    cutXibar->SetCascadeMaxDCAVertex(1e5); // not using
    cutXibar->SetCascadeMinDCAVertex(-1e5); // not using
    cutXibar->SetV0LowRadius(0); // not using
    cutXibar->SetV0HighRadius(1e5); // not using
    cutXibar->SetCascadeLowRadius(0); // not using
    cutXibar->SetCascadeHighRadius(1e5); // not using
    cutXibar->SetMassTolerance(Xi_massTol);
    cutXibar->SetMassToleranceVeto(Xi_massTolVeto);//Rejection range for Competing Xi Rejection
    cutXibar->SetSwitch(kFALSE); // not using
    cutXibar->SetV0MinCosPointingAngle(V0CosPoinAn);
    cutXibar->SetCascadeMinCosPointingAngle(XiCosPoinAn);
    cutXibar->SetMaxPseudorapidity(0.8);
    cutXibar->SetMinTPCcluster(-1);
    
    AliRsnCutSet* cutSetXi=new AliRsnCutSet("setXi",AliRsnTarget::kDaughter);
    cutSetXi->AddCut(cutXi);
    cutSetXi->SetCutScheme(cutXi->GetName());
    Int_t icutXi=task->AddTrackCuts(cutSetXi);
    
    AliRsnCutSet* cutSetXibar=new AliRsnCutSet("setXibar",AliRsnTarget::kDaughter);
    cutSetXibar->AddCut(cutXibar);
    cutSetXibar->SetCutScheme(cutXibar->GetName());
    Int_t icutXibar=task->AddTrackCuts(cutSetXibar);
    
    // monitoring
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputCascade.C");
        
        AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());
        AddMonitorOutputV0(isMC, cutSetXi->GetMonitorOutput(), "Lambda", "nokine");
        AddMonitorOutputV0(isMC, cutSetXibar->GetMonitorOutput(), "AntiLambda", "nokine");
        AddMonitorOutputCascade(isMC, cutSetXi->GetMonitorOutput(), "Xi");
        AddMonitorOutputCascade(isMC, cutSetXibar->GetMonitorOutput(), "AntiXi");
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
    /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass,    kFALSE);
    /* mother mass      */ Int_t mmID   = task->CreateValue(AliRsnMiniValue::kInvMassMother,kFALSE);
    /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
    /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
    /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
    /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
    /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);
    /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt,kFALSE);
    /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt,kFALSE);
    
    // -- Create all needed outputs -----------------------------------------------------------------
    // use an array for more compact writing, which are different on mixing and charges
    
    Int_t i,k,xID,cut1,pairID,ipdg;
    TString name,comp;
    Char_t charge1, charge2;
    AliRsnMiniOutput* out;
    
    task->SetMotherAcceptanceCutMinPt(0.15);
    task->SetMotherAcceptanceCutMaxEta(0.8);
    task->KeepMotherInAcceptance(true);
    
    for(i=0;i<2;i++) for(j=0;j<2;j++) for(k=0;k<9;k++){
        if(!i){
            name.Form("Xim");//Xi-
            cut1=icutXi;
            charge1='-';
        }else{
            name.Form("Xip");//anti-Xi+
            cut1=icutXibar;
            charge1='+';
        }
        
        if(!j){
            name.Append("Pip");
            charge2='+';
        }else{
            name.Append("Pim");
            charge2='-';
        }
        
        if(!isMC && k>=3) continue;
        pairID=1;
        xID=imID;
        if(!k){
            comp.Form("PAIR");
            pairID=0;
        }else if(k==1){
            name.Append("Mix");
            comp.Form("MIX");
        }else if(k==2){
            name.Append("Rotated");
            comp.Form("ROTATE1");
            pairID=0;
        }else if(k==3){
            name.Append("_gen");
            comp.Form("MOTHER");
        }else if(k==4){
            name.Append("_rec");
            comp.Form("TRUE");
        }else if(k==5){
            name.Append("_recMM");
            comp.Form("TRUE");
            xID=mmID;
        }else if(k==6){
            name.Append("_res");
            comp.Form("TRUE");
            xID=diffID;
        }else if(k==7){
            name.Append("_genPS");
            comp.Form("MOTHER_IN_ACC");
        }else if(k==8){
            name.Append("_recPS");
            comp.Form("TRUE");
        }
        
        ipdg=(i==0)?3324:-3324;
        
        out=task->CreateOutput(Form("Xipi_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kOmega);//kXi);
        out->SetCutID(0,cut1);
        out->SetCharge(0,charge1);
        
        out->SetDaughter(1,AliRsnDaughter::kPion);
        out->SetCutID(1,iCutPi);
        out->SetCharge(1,charge2);
        
        if(!pairID) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(ipdg);
        out->SetMotherMass(mass);
        
        if(k<=6){
            if(xID==imID || xID==mmID) out->AddAxis(xID,300,1.8,3);//200,1.4,1.8);// axis X: invmass or resolution
            else out->AddAxis(diffID,200,-0.02,0.02);
            out->AddAxis(ptID,200,0.0,20.0);// axis Y: transverse momentum
            out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
        }else{//Phase-space histograms
            out->AddAxis(fdpt,100,0.,10.);
            out->AddAxis(sdpt,100,0.,10.);
            out->AddAxis(ptID,40,0.,20.);
        }
    }
    
    if(isMC){
        AliRsnCutMiniPair* cutEta=new AliRsnCutMiniPair("cutPseudorapidity", AliRsnCutMiniPair::kPseudorapidityRangeMC);
        cutEta->SetRangeD(-0.8,0.8);
        
        AliRsnCutSet* cutsPairDaughter=new AliRsnCutSet("pairCutsDaughter", AliRsnTarget::kMother);
        cutsPairDaughter->AddCut(cutEta);
        cutsPairDaughter->SetCutScheme(cutEta->GetName());
        
        //for efficiency of pi+/-
        out = task->CreateOutput("Xipi_Pip_mother", "HIST", "SINGLE");
        out->SetDaughter(0, AliRsnDaughter::kUnknown);
        out->SetDaughter(1, AliRsnDaughter::kUnknown);
        out->SetMotherPDG(211);
        out->SetMotherMass(0.139);
        out->SetPairCuts(cutsPairDaughter);
        out->AddAxis(ptID,200,0.0,20.0);
        
        out = task->CreateOutput("Xipi_Pim_mother", "HIST", "SINGLE");
        out->SetDaughter(0, AliRsnDaughter::kUnknown);
        out->SetDaughter(1, AliRsnDaughter::kUnknown);
        out->SetMotherPDG(-211);
        out->SetMotherMass(0.139);
        out->SetPairCuts(cutsPairDaughter);
        out->AddAxis(ptID,200,0.0,20.0);
        
        //for efficiency of (anti-)Xi
        out = task->CreateOutput("Xipi_Xim_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kLambda);
        out->SetMotherPDG(3312);
        out->SetMotherMass(1.32171);
        out->SetPairCuts(cutsPairDaughter);
        out->AddAxis(ptID,200,0.0,20.0);
        
        out = task->CreateOutput("Xipi_Xip_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kLambda);
        out->SetMotherPDG(-3312);
        out->SetMotherMass(1.32171);
        out->SetPairCuts(cutsPairDaughter);
        out->AddAxis(ptID,200,0.0,20.0);
    }
    
    return kTRUE;
}


//=============================


Bool_t Config_Xikx(
                   AliRsnMiniAnalysisTask *task,
                   TString     lname,
                   Bool_t      isMC,
                   Int_t       system,
                   Int_t       EventCuts,
                   Int_t       TrackCutsXi,
                   Int_t       TrackCutsK
                   ){
    bool isPP=false;
    if(!system) isPP=true;
    int trigger=EventCuts%10;
    int MultBins=(EventCuts/10)%10;
    if(system==1 || system==2) MultBins=1;
    
    char suffix[1000];
    sprintf(suffix,"_%s",lname.Data());
    Bool_t enableMonitor=kTRUE;
    
    Double_t mass= 0.493677+1672.43;//1.32171;
    
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
    
    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts("qualityTrackcut");
    esdTrackCuts->SetEtaRange(-0.8,0.8);
    esdTrackCuts->SetRequireTPCRefit();
    esdTrackCuts->SetAcceptKinkDaughters(0);
    esdTrackCuts->SetMinNCrossedRowsTPC(70);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetMinDCAToVertexXY(0.06);
    
    // selections for Xi
    Float_t XiPIDcut=5.;
    Float_t V0dDCA=1.6;
    Float_t XidDCA=1.6;
    Float_t XiMinDCA=0.07;
    Float_t Xi_massTol=0.007;
    Float_t Xi_massTolVeto=0.007;
    Float_t V0CosPoinAn=0.97;
    Float_t XiCosPoinAn=0.97;
    
    AliRsnCutCascade* cutXi=new AliRsnCutCascade("cutXi",kOmegaMinus);//kXiMinus);
    cutXi->SetPIDCutV0Proton(XiPIDcut);
    cutXi->SetPIDCutV0Pion(XiPIDcut);
    cutXi->SetPIDCutBachelor(XiPIDcut);
    cutXi->SetESDtrackCuts(esdTrackCuts);
    cutXi->SetV0MaxDaughtersDCA(V0dDCA);
    cutXi->SetCascadeMaxDaughtersDCA(XidDCA);
    cutXi->SetV0MaxDCAVertex(1e5); // not using
    cutXi->SetV0MinDCAVertex(XiMinDCA);
    cutXi->SetCascadeMaxDCAVertex(1e5); // not using
    cutXi->SetCascadeMinDCAVertex(-1e5); // not using
    cutXi->SetV0LowRadius(0); // not using
    cutXi->SetV0HighRadius(1e5); // not using
    cutXi->SetCascadeLowRadius(0); // not using
    cutXi->SetCascadeHighRadius(1e5); // not using
    cutXi->SetMassTolerance(Xi_massTol);
    cutXi->SetMassToleranceVeto(Xi_massTolVeto);//Rejection range for Competing Xi Rejection
    cutXi->SetSwitch(kFALSE); // not using
    cutXi->SetV0MinCosPointingAngle(V0CosPoinAn);
    cutXi->SetCascadeMinCosPointingAngle(XiCosPoinAn);
    cutXi->SetMaxPseudorapidity(0.8);
    cutXi->SetMinTPCcluster(-1);
    
    AliRsnCutCascade* cutXibar=new AliRsnCutCascade("cutXibar",kOmegaPlusBar);//kXiPlusBar);
    cutXibar->SetPIDCutV0Proton(XiPIDcut);
    cutXibar->SetPIDCutV0Pion(XiPIDcut);
    cutXibar->SetPIDCutBachelor(XiPIDcut);
    cutXibar->SetESDtrackCuts(esdTrackCuts);
    cutXibar->SetV0MaxDaughtersDCA(V0dDCA);
    cutXibar->SetCascadeMaxDaughtersDCA(XidDCA);
    cutXibar->SetV0MaxDCAVertex(1e5); // not using
    cutXibar->SetV0MinDCAVertex(XiMinDCA);
    cutXibar->SetCascadeMaxDCAVertex(1e5); // not using
    cutXibar->SetCascadeMinDCAVertex(-1e5); // not using
    cutXibar->SetV0LowRadius(0); // not using
    cutXibar->SetV0HighRadius(1e5); // not using
    cutXibar->SetCascadeLowRadius(0); // not using
    cutXibar->SetCascadeHighRadius(1e5); // not using
    cutXibar->SetMassTolerance(Xi_massTol);
    cutXibar->SetMassToleranceVeto(Xi_massTolVeto);//Rejection range for Competing Xi Rejection
    cutXibar->SetSwitch(kFALSE); // not using
    cutXibar->SetV0MinCosPointingAngle(V0CosPoinAn);
    cutXibar->SetCascadeMinCosPointingAngle(XiCosPoinAn);
    cutXibar->SetMaxPseudorapidity(0.8);
    cutXibar->SetMinTPCcluster(-1);
    
    AliRsnCutSet* cutSetXi=new AliRsnCutSet("setXi",AliRsnTarget::kDaughter);
    cutSetXi->AddCut(cutXi);
    cutSetXi->SetCutScheme(cutXi->GetName());
    Int_t icutXi=task->AddTrackCuts(cutSetXi);
    
    AliRsnCutSet* cutSetXibar=new AliRsnCutSet("setXibar",AliRsnTarget::kDaughter);
    cutSetXibar->AddCut(cutXibar);
    cutSetXibar->SetCutScheme(cutXibar->GetName());
    Int_t icutXibar=task->AddTrackCuts(cutSetXibar);
    
    // monitoring
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputCascade.C");
        
        AddMonitorOutput(isMC,cutSetK->GetMonitorOutput());
        AddMonitorOutputV0(isMC, cutSetXi->GetMonitorOutput(), "Lambda", "nokine");
        AddMonitorOutputV0(isMC, cutSetXibar->GetMonitorOutput(), "AntiLambda", "nokine");
        AddMonitorOutputCascade(isMC, cutSetXi->GetMonitorOutput(), "Omega");//"Xi");
        AddMonitorOutputCascade(isMC, cutSetXibar->GetMonitorOutput(), "AntiOmega");//"AntiXi");
    }
    
    // pair cuts
    AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
    if(system!=1) cutY->SetRangeD(-0.8,0.8);
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
    /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass,    kFALSE);
    /* mother mass      */ Int_t mmID   = task->CreateValue(AliRsnMiniValue::kInvMassMother,kFALSE);
    /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
    /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
    /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
    /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
    /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);
    /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt,kFALSE);
    /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt,kFALSE);
    
    // -- Create all needed outputs -----------------------------------------------------------------
    // use an array for more compact writing, which are different on mixing and charges
    
    Int_t i,k,xID,cut1,pairID,ipdg;
    TString name,comp;
    Char_t charge1, charge2;
    AliRsnMiniOutput* out;
    
    task->SetMotherAcceptanceCutMinPt(0.15);
    task->SetMotherAcceptanceCutMaxEta(0.8);
    task->KeepMotherInAcceptance(true);
    
    for(i=0;i<2;i++) for(j=0;j<2;j++) for(k=0;k<9;k++){
        if(!i){
            name.Form("Xim");//Xi-
            cut1=icutXi;
            charge1='-';
        }else{
            name.Form("Xip");//anti-Xi+
            cut1=icutXibar;
            charge1='+';
        }
        
        if(!j){
            name.Append("Kp");
            charge2='+';
        }else{
            name.Append("Km");
            charge2='-';
        }
        
        if(!isMC && k>=3) continue;
        pairID=1;
        xID=imID;
        if(!k){
            comp.Form("PAIR");
            pairID=0;
        }else if(k==1){
            name.Append("Mix");
            comp.Form("MIX");
        }else if(k==2){
            name.Append("Rotated");
            comp.Form("ROTATE1");
            pairID=0;
        }else if(k==3){
            name.Append("_gen");
            comp.Form("MOTHER");
        }else if(k==4){
            name.Append("_rec");
            comp.Form("TRUE");
        }else if(k==5){
            name.Append("_recMM");
            comp.Form("TRUE");
            xID=mmID;
        }else if(k==6){
            name.Append("_res");
            comp.Form("TRUE");
            xID=diffID;
        }else if(k==7){
            name.Append("_genPS");
            comp.Form("MOTHER_IN_ACC");
        }else if(k==8){
            name.Append("_recPS");
            comp.Form("TRUE");
        }
        
        ipdg=(i==0)?3324:-3324;
        
        out=task->CreateOutput(Form("Xikx_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kOmega);//kXi);
        out->SetCutID(0,cut1);
        out->SetCharge(0,charge1);
        
        out->SetDaughter(1,AliRsnDaughter::kKaon);
        out->SetCutID(1,iCutK);
        out->SetCharge(1,charge2);
        
        if(!pairID) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(ipdg);
        out->SetMotherMass(mass);
        
        if(k<=6){
            if(xID==imID || xID==mmID) out->AddAxis(xID,240,1.8,3.);// axis X: invmass or resolution
            else out->AddAxis(diffID,200,-0.02,0.02);
            out->AddAxis(ptID,200,0.0,20.0);// axis Y: transverse momentum
            out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
        }else{//Phase-space histograms
            out->AddAxis(fdpt,100,0.,10.);
            out->AddAxis(sdpt,100,0.,10.);
            out->AddAxis(ptID,40,0.,20.);
        }
    }
    
    if(isMC){
        AliRsnCutMiniPair* cutEta=new AliRsnCutMiniPair("cutPseudorapidity", AliRsnCutMiniPair::kPseudorapidityRangeMC);
        cutEta->SetRangeD(-0.8,0.8);
        
        AliRsnCutSet* cutsPairDaughter=new AliRsnCutSet("pairCutsDaughter", AliRsnTarget::kMother);
        cutsPairDaughter->AddCut(cutEta);
        cutsPairDaughter->SetCutScheme(cutEta->GetName());
        
        //for efficiency of K+/-
        out = task->CreateOutput("Xikx_Kp_mother", "HIST", "SINGLE");
        out->SetDaughter(0, AliRsnDaughter::kUnknown);
        out->SetDaughter(1, AliRsnDaughter::kUnknown);
        out->SetMotherPDG(321);
        out->SetMotherMass(0.493677);
        out->SetPairCuts(cutsPairDaughter);
        out->AddAxis(ptID,200,0.0,20.0);
        
        out = task->CreateOutput("Xikx_Km_mother", "HIST", "SINGLE");
        out->SetDaughter(0, AliRsnDaughter::kUnknown);
        out->SetDaughter(1, AliRsnDaughter::kUnknown);
        out->SetMotherPDG(-321);
        out->SetMotherMass(0.493677);
        out->SetPairCuts(cutsPairDaughter);
        out->AddAxis(ptID,200,0.0,20.0);
        
        //for efficiency of (anti-)Xi
        out = task->CreateOutput("Xikx_Xim_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kLambda);
        out->SetMotherPDG(3312);
        out->SetMotherMass(1.32171);
        out->SetPairCuts(cutsPairDaughter);
        out->AddAxis(ptID,200,0.0,20.0);
        
        out = task->CreateOutput("Xikx_Xip_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kLambda);
        out->SetMotherPDG(-3312);
        out->SetMotherMass(1.32171);
        out->SetPairCuts(cutsPairDaughter);
        out->AddAxis(ptID,200,0.0,20.0);
    }
    
    return kTRUE;
}


//=============================


Bool_t Config_Xik0(
                   AliRsnMiniAnalysisTask *task,
                   TString     lname,
                   Bool_t      isMC,
                   Int_t       system,
                   Int_t       EventCuts,
                   Int_t       TrackCutsXi,
                   Int_t       TrackCutsK
                   ){
    bool isPP=false;
    if(!system) isPP=true;
    int trigger=EventCuts%10;
    int MultBins=(EventCuts/10)%10;
    if(system==1 || system==2) MultBins=1;
    
    char suffix[1000];
    sprintf(suffix,"_%s",lname.Data());
    Bool_t enableMonitor=kTRUE;
    
    Double_t mass= 0.497611+1.32171;
    
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
    Float_t k0s_radiuslow=0.9; // 0.5;
    Float_t k0s_radiushigh=200.;
    Float_t k0s_massTolSigma=4;
    Int_t   k0s_massTolID=0;
    Float_t k0s_massTol=0.03;
    Float_t k0s_massTolVeto=0.004;
    Bool_t  k0sSwitch=kTRUE;
    Float_t k0sCosPoinAn=0.97;
    
    if(TrackCutsK/100000) k0s_massTolID=1;//use pT-dependent mass tolerance cut
    if((TrackCutsK/10000)%10==1) k0sSwitch=kFALSE;//no competing V0 rejection
    if((TrackCutsK/100)%100) k0s_pLife=(TrackCutsK/100)%100;
    if(TrackCutsK%100) k0s_piPIDCut=((float) (TrackCutsK%100))*0.1;

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
    cutK0s->SetMaxPseudorapidity(0.8);
    cutK0s->SetMinTPCcluster(-1);
    
    AliRsnCutSet* cutSetK0s=new AliRsnCutSet("setK0s",AliRsnTarget::kDaughter);
    cutSetK0s->AddCut(cutK0s);
    cutSetK0s->SetCutScheme(cutK0s->GetName());
    Int_t iCutK0s=task->AddTrackCuts(cutSetK0s);
    
    // selections for Xi
    Float_t XiPPIDcut=3.;
    Float_t XiPiPIDcut=5.;
    Float_t V0dDCA=1.6;
    Float_t XidDCA=1.4; // 1.6;
    Float_t XiMinDCA=0.07;
    Float_t Xi_massTol=0.015;
    Float_t Xi_massTolVeto=0.007;
    Float_t Xi_V0massTol=0.006;
    Float_t V0CosPoinAn=0.97;
    Float_t XiCosPoinAn=0.97;
    Float_t V0lifetime=40.;

    Bool_t CheckOOBP=false;//true;
    if((TrackCutsXi/1000000)==1) CheckOOBP=true;
    if((TrackCutsXi/10000)%100) V0lifetime=(TrackCutsXi/10000)%100;
    if(V0lifetime>98.5) V0lifetime=1e20;//turn off
    if((TrackCutsXi/100)%100) XiPPIDcut=((float) ((TrackCutsXi/100)%100))*0.1;
    if(TrackCutsXi%100) XiPiPIDcut=((float) (TrackCutsXi%100))*0.1;

    AliRsnCutCascade* cutXi=new AliRsnCutCascade("cutXi",kXiMinus);
    cutXi->SetPIDCutV0Proton(XiPPIDcut);
    cutXi->SetPIDCutV0Pion(XiPiPIDcut);
    cutXi->SetPIDCutBachelor(XiPiPIDcut);
    cutXi->SetESDtrackCuts(esdTrackCuts);
    cutXi->SetV0MaxDaughtersDCA(V0dDCA);
    cutXi->SetCascadeMaxDaughtersDCA(XidDCA);
    cutXi->SetV0MaxDCAVertex(1e5); // not using
    cutXi->SetV0MinDCAVertex(XiMinDCA);
    cutXi->SetCascadeMaxDCAVertex(1e5); // not using
    cutXi->SetCascadeMinDCAVertex(-1e5); // not using
    cutXi->SetV0LowRadius(0.9); // 0
    cutXi->SetV0HighRadius(1e5); // not using
    cutXi->SetCascadeLowRadius(0.5); // 0
    cutXi->SetCascadeHighRadius(1e5); // not using
    cutXi->SetV0Life(V0lifetime);
    cutXi->SetMassTolerance(Xi_massTol);
    cutXi->SetMassToleranceVeto(Xi_massTolVeto);//Rejection range for Competing Xi Rejection
    cutXi->SetV0MassTolerance(Xi_V0massTol);
    cutXi->SetSwitch(kFALSE); // not using
    cutXi->SetV0MinCosPointingAngle(V0CosPoinAn);
    cutXi->SetCascadeMinCosPointingAngle(XiCosPoinAn);
    cutXi->SetMaxPseudorapidity(0.8);
    cutXi->SetMinTPCcluster(-1);
    
    AliRsnCutCascade* cutXibar=new AliRsnCutCascade("cutXibar",kXiPlusBar);
    cutXibar->SetPIDCutV0Proton(XiPPIDcut);
    cutXibar->SetPIDCutV0Pion(XiPiPIDcut);
    cutXibar->SetPIDCutBachelor(XiPiPIDcut);
    cutXibar->SetESDtrackCuts(esdTrackCuts);
    cutXibar->SetV0MaxDaughtersDCA(V0dDCA);
    cutXibar->SetCascadeMaxDaughtersDCA(XidDCA);
    cutXibar->SetV0MaxDCAVertex(1e5); // not using
    cutXibar->SetV0MinDCAVertex(XiMinDCA);
    cutXibar->SetCascadeMaxDCAVertex(1e5); // not using
    cutXibar->SetCascadeMinDCAVertex(-1e5); // not using
    cutXibar->SetV0LowRadius(0.9); // 0
    cutXibar->SetV0HighRadius(1e5); // not using
    cutXibar->SetCascadeLowRadius(0.5); // 0
    cutXibar->SetCascadeHighRadius(1e5); // not using
    cutXibar->SetV0Life(V0lifetime);
    cutXibar->SetMassTolerance(Xi_massTol);
    cutXibar->SetMassToleranceVeto(Xi_massTolVeto);//Rejection range for Competing Xi Rejection
    cutXibar->SetV0MassTolerance(Xi_V0massTol);
    cutXibar->SetSwitch(kFALSE); // not using
    cutXibar->SetV0MinCosPointingAngle(V0CosPoinAn);
    cutXibar->SetCascadeMinCosPointingAngle(XiCosPoinAn);
    cutXibar->SetMaxPseudorapidity(0.8);
    cutXibar->SetMinTPCcluster(-1);
    
    AliRsnCutSet* cutSetXi=new AliRsnCutSet("setXi",AliRsnTarget::kDaughter);
    cutSetXi->AddCut(cutXi);
    cutSetXi->SetCutScheme(cutXi->GetName());
    Int_t icutXi=task->AddTrackCuts(cutSetXi);
    
    AliRsnCutSet* cutSetXibar=new AliRsnCutSet("setXibar",AliRsnTarget::kDaughter);
    cutSetXibar->AddCut(cutXibar);
    cutSetXibar->SetCutScheme(cutXibar->GetName());
    Int_t icutXibar=task->AddTrackCuts(cutSetXibar);
    
    // monitoring
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputCascade.C");
        
        AddMonitorOutputV0(isMC, cutSetK0s->GetMonitorOutput(), "K0S");
        AddMonitorOutputV0(isMC, cutSetXi->GetMonitorOutput(), "Lambda", "nokine");
        AddMonitorOutputV0(isMC, cutSetXibar->GetMonitorOutput(), "AntiLambda", "nokine");
        AddMonitorOutputCascade(isMC, cutSetXi->GetMonitorOutput(), "Xi");
        AddMonitorOutputCascade(isMC, cutSetXibar->GetMonitorOutput(), "AntiXi");
    }
    
    // pair cuts
    AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
    if(system!=1) cutY->SetRangeD(-0.8,0.8);
    else cutY->SetRangeD(-0.465,0.035);
    
    AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
    AliRsnCutMiniPair* cutOOBP=new AliRsnCutMiniPair("cutOOBP",AliRsnCutMiniPair::kPassesOOBPileupCut);
    
    AliRsnCutSet* cutsPairSame=new AliRsnCutSet("pairCutsSame",AliRsnTarget::kMother);
    cutsPairSame->AddCut(cutY);
    cutsPairSame->AddCut(cutV0);
    
    AliRsnCutSet* cutsPairMix=new AliRsnCutSet("pairCutsMix", AliRsnTarget::kMother);
    cutsPairMix->AddCut(cutY);
    
    if(CheckOOBP){
        cutsPairSame->AddCut(cutOOBP);
        cutsPairSame->SetCutScheme(TString::Format("%s&(!%s)&%s",cutY->GetName(),cutV0->GetName(),cutOOBP->GetName()).Data());
        cutsPairMix->AddCut(cutOOBP);
        cutsPairMix->SetCutScheme(TString::Format("%s&%s",cutY->GetName(),cutOOBP->GetName()).Data());
    }else{
        cutsPairSame->SetCutScheme(TString::Format("%s&(!%s)",cutY->GetName(),cutV0->GetName()).Data());
        cutsPairMix->SetCutScheme(cutY->GetName());
    }
    
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
    /* mother mass      */ Int_t mmID   = task->CreateValue(AliRsnMiniValue::kInvMassMother,kFALSE);
    /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
    /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
    /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
    /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
    /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);
    /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt,kFALSE);
    /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt,kFALSE);
    /* 1st daughter pt  */ Int_t fdptmc = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt,kTRUE);
    /* 2nd daughter pt  */ Int_t sdptmc = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt,kTRUE);
    
    // -- Create all needed outputs -----------------------------------------------------------------
    // use an array for more compact writing, which are different on mixing and charges
    
    Int_t i,xID,cut1,pairID,ipdg;
    TString name,comp;
    Char_t charge1;
    AliRsnMiniOutput* out;
    
    task->SetMotherAcceptanceCutMinPt(0.15);
    task->SetMotherAcceptanceCutMaxEta(0.8);
    task->KeepMotherInAcceptance(true);
    
    for(i=0;i<2;i++) for(j=0;j<9;j++){
        if(!i){
            name.Form("XimK0");
            cut1=icutXi;
            charge1='-';
            ipdg=3335;
        }else{
            name.Form("XipK0");
            cut1=icutXibar;
            charge1='+';
            ipdg=-3335;
        }
        
        if(!isMC && j>=3) continue;
        pairID=1;
        xID=imID;
        if(!j){
            comp.Form("PAIR");
            pairID=0;
        }else if(j==1){
            name.Append("Mix");
            comp.Form("MIX");
        }else if(j==2){
            name.Append("Rotated");
            comp.Form("ROTATE1");
            pairID=0;
        }else if(j==3){
            name.Append("_gen");
            comp.Form("MOTHER");
        }else if(j==4){
            name.Append("_rec");
            comp.Form("TRUE");
        }else if(j==5){
            name.Append("_recMM");
            comp.Form("TRUE");
            xID=mmID;
        }else if(j==6){
            name.Append("_res");
            comp.Form("TRUE");
            xID=diffID;
        }else if(j==7){
            name.Append("_genPS");
            comp.Form("MOTHER_IN_ACC");
        }else if(j==8){
            name.Append("_recPS");
            comp.Form("TRUE");
        }
        
        out=task->CreateOutput(Form("Xik0_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kXi);
        out->SetCutID(0,cut1);
        out->SetCharge(0,charge1);
        
        out->SetDaughter(1,AliRsnDaughter::kKaon0);
        out->SetCutID(1,iCutK0s);
        out->SetCharge(1,'0');
        
        if(!pairID) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(ipdg);
        out->SetMotherMass(mass);
        
        if(j<=6){
            //if(xID==imID) out->AddAxis(imID,240,1.8,3);// axis X: invmass
            if(xID==imID || xID==mmID) out->AddAxis(xID,400,1.8,2.2);// axis X: invmass
            else out->AddAxis(diffID,200,-0.02,0.02);// axis X: resolution
            out->AddAxis(ptID,200,0.0,20.0);// axis Y: transverse momentum
            out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
        }else if(j==7){//Phase-space histograms
            out->AddAxis(fdptmc,100,0.,10.);
            out->AddAxis(sdptmc,100,0.,10.);
            out->AddAxis(ptID,40,0.,20.);
        }else{
            out->AddAxis(fdpt,100,0.,10.);
            out->AddAxis(sdpt,100,0.,10.);
            out->AddAxis(ptID,40,0.,20.);
        }
    }
    
    if(isMC){
        AliRsnCutMiniPair* cutEta=new AliRsnCutMiniPair("cutPseudorapidity", AliRsnCutMiniPair::kPseudorapidityRangeMC);
        cutEta->SetRangeD(-0.8,0.8);
        
        AliRsnCutSet* cutsPairDaughter=new AliRsnCutSet("pairCutsDaughter", AliRsnTarget::kMother);
        cutsPairDaughter->AddCut(cutEta);
        cutsPairDaughter->SetCutScheme(cutEta->GetName());
        
        //for efficiency of K0S
        out = task->CreateOutput("Xik0_K0S_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kPion);
        out->SetMotherPDG(310);
        out->SetMotherMass(0.497611);
        out->SetPairCuts(cutsPairDaughter);
        out->AddAxis(ptID,200,0.0,20.0);
        
        //for efficiency of (anti-)Xi
        out = task->CreateOutput("Xik0_Xim_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kLambda);
        out->SetMotherPDG(3312);
        out->SetMotherMass(1.32171);
        out->SetPairCuts(cutsPairDaughter);
        out->AddAxis(ptID,200,0.0,20.0);
        
        out = task->CreateOutput("Xik0_Xip_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kLambda);
        out->SetMotherPDG(-3312);
        out->SetMotherMass(1.32171);
        out->SetPairCuts(cutsPairDaughter);
        out->AddAxis(ptID,200,0.0,20.0);
    }
    
    return kTRUE;
}


//=============================


Bool_t Config_XiP(
                  AliRsnMiniAnalysisTask *task,
                  TString     lname,
                  Bool_t      isMC,
                  Int_t       system,
                  Int_t       EventCuts,
                  Int_t       TrackCutsXi,
                  Int_t       TrackCutsP
                  ){
    bool isPP=false;
    if(!system) isPP=true;
    int trigger=EventCuts%10;
    int MultBins=(EventCuts/10)%10;
    if(system==1 || system==2) MultBins=1;
    
    char suffix[1000];
    sprintf(suffix,"_%s",lname.Data());
    Bool_t enableMonitor=kTRUE;
    
    Double_t mass= 0.938272+1.32171;
    
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
    Int_t fineBin=(TrackCutsP%100000)/10000; //option for 2 MeV bins
    
    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts("qualityTrackcut");
    esdTrackCuts->SetEtaRange(-0.8,0.8);
    esdTrackCuts->SetRequireTPCRefit();
    esdTrackCuts->SetAcceptKinkDaughters(0);
    esdTrackCuts->SetMinNCrossedRowsTPC(70);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetMinDCAToVertexXY(0.06);
    
    // selections for Xi
    Float_t XiPIDcut=3.;
    Float_t V0dDCA=1.6;
    Float_t XidDCA=1.6;
    Float_t XiMinDCA=0.07;
    Float_t Xi_massTol=0.007;
    Float_t Xi_massTolVeto=0.007;
    Float_t V0CosPoinAn=0.97;
    Float_t XiCosPoinAn=0.97;
    
    AliRsnCutCascade* cutXi=new AliRsnCutCascade("cutXi",kXiMinus);
    cutXi->SetPIDCutV0Proton(XiPIDcut);
    cutXi->SetPIDCutV0Pion(XiPIDcut);
    cutXi->SetPIDCutBachelor(XiPIDcut);
    cutXi->SetESDtrackCuts(esdTrackCuts);
    cutXi->SetV0MaxDaughtersDCA(V0dDCA);
    cutXi->SetCascadeMaxDaughtersDCA(XidDCA);
    cutXi->SetV0MaxDCAVertex(1e5); // not using
    cutXi->SetV0MinDCAVertex(XiMinDCA);
    cutXi->SetCascadeMaxDCAVertex(1e5); // not using
    cutXi->SetCascadeMinDCAVertex(-1e5); // not using
    cutXi->SetV0LowRadius(0); // not using
    cutXi->SetV0HighRadius(1e5); // not using
    cutXi->SetCascadeLowRadius(0); // not using
    cutXi->SetCascadeHighRadius(1e5); // not using
    cutXi->SetMassTolerance(Xi_massTol);
    cutXi->SetMassToleranceVeto(Xi_massTolVeto);//Rejection range for Competing Xi Rejection
    cutXi->SetSwitch(kFALSE); // not using
    cutXi->SetV0MinCosPointingAngle(V0CosPoinAn);
    cutXi->SetCascadeMinCosPointingAngle(XiCosPoinAn);
    cutXi->SetMaxPseudorapidity(0.8);
    cutXi->SetMinTPCcluster(-1);
    
    AliRsnCutCascade* cutXibar=new AliRsnCutCascade("cutXibar",kXiPlusBar);
    cutXibar->SetPIDCutV0Proton(XiPIDcut);
    cutXibar->SetPIDCutV0Pion(XiPIDcut);
    cutXibar->SetPIDCutBachelor(XiPIDcut);
    cutXibar->SetESDtrackCuts(esdTrackCuts);
    cutXibar->SetV0MaxDaughtersDCA(V0dDCA);
    cutXibar->SetCascadeMaxDaughtersDCA(XidDCA);
    cutXibar->SetV0MaxDCAVertex(1e5); // not using
    cutXibar->SetV0MinDCAVertex(XiMinDCA);
    cutXibar->SetCascadeMaxDCAVertex(1e5); // not using
    cutXibar->SetCascadeMinDCAVertex(-1e5); // not using
    cutXibar->SetV0LowRadius(0); // not using
    cutXibar->SetV0HighRadius(1e5); // not using
    cutXibar->SetCascadeLowRadius(0); // not using
    cutXibar->SetCascadeHighRadius(1e5); // not using
    cutXibar->SetMassTolerance(Xi_massTol);
    cutXibar->SetMassToleranceVeto(Xi_massTolVeto);//Rejection range for Competing Xi Rejection
    cutXibar->SetSwitch(kFALSE); // not using
    cutXibar->SetV0MinCosPointingAngle(V0CosPoinAn);
    cutXibar->SetCascadeMinCosPointingAngle(XiCosPoinAn);
    cutXibar->SetMaxPseudorapidity(0.8);
    cutXibar->SetMinTPCcluster(-1);
    
    AliRsnCutSet* cutSetXi=new AliRsnCutSet("setXi",AliRsnTarget::kDaughter);
    cutSetXi->AddCut(cutXi);
    cutSetXi->SetCutScheme(cutXi->GetName());
    Int_t icutXi=task->AddTrackCuts(cutSetXi);
    
    AliRsnCutSet* cutSetXibar=new AliRsnCutSet("setXibar",AliRsnTarget::kDaughter);
    cutSetXibar->AddCut(cutXibar);
    cutSetXibar->SetCutScheme(cutXibar->GetName());
    Int_t icutXibar=task->AddTrackCuts(cutSetXibar);
    
    // monitoring
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputCascade.C");
        
        AddMonitorOutput(isMC,cutSetP->GetMonitorOutput());
        AddMonitorOutputV0(isMC, cutSetXi->GetMonitorOutput(), "Lambda", "nokine");
        AddMonitorOutputV0(isMC, cutSetXibar->GetMonitorOutput(), "AntiLambda", "nokine");
        AddMonitorOutputCascade(isMC, cutSetXi->GetMonitorOutput(), "Xi");
        AddMonitorOutputCascade(isMC, cutSetXibar->GetMonitorOutput(), "AntiXi");
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
    /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass,    kFALSE);
    /* mother mass      */ Int_t mmID   = task->CreateValue(AliRsnMiniValue::kInvMassMother,kFALSE);
    /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
    /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
    /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
    /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
    /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);
    /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt,kFALSE);
    /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt,kFALSE);
    
    // -- Create all needed outputs -----------------------------------------------------------------
    // use an array for more compact writing, which are different on mixing and charges
    
    Int_t i,k,xID,cut1,pairID,ipdg;
    TString name,comp;
    Char_t charge1, charge2;
    AliRsnMiniOutput* out;
    
    task->SetMotherAcceptanceCutMinPt(0.15);
    task->SetMotherAcceptanceCutMaxEta(0.8);
    task->KeepMotherInAcceptance(true);
    
    for(i=0;i<2;i++) for(j=0;j<2;j++) for(k=0;k<9;k++){
        if(!i){
            name.Form("Xim");//Xi-
            cut1=icutXi;
            charge1='-';
            ipdg=3324;
        }else{
            name.Form("Xip");//anti-Xi+
            cut1=icutXibar;
            charge1='+';
            ipdg=-3324;
        }
        
        if(!j){
            name.Append("Pp");
            charge2='+';
        }else{
            name.Append("Pm");
            charge2='-';
        }
        
        if(!isMC && k>=3) continue;
        pairID=1;
        xID=imID;
        if(!k){
            comp.Form("PAIR");
            pairID=0;
        }else if(k==1){
            name.Append("Mix");
            comp.Form("MIX");
        }else if(k==2){
            name.Append("Rotated");
            comp.Form("ROTATE1");
            pairID=0;
        }else if(k==3){
            name.Append("_gen");
            comp.Form("MOTHER");
        }else if(k==4){
            name.Append("_rec");
            comp.Form("TRUE");
        }else if(k==5){
            name.Append("_recMM");
            comp.Form("TRUE");
            xID=mmID;
        }else if(k==6){
            name.Append("_res");
            comp.Form("TRUE");
            xID=diffID;
        }else if(k==7){
            name.Append("_genPS");
            comp.Form("MOTHER_IN_ACC");
        }else if(k==8){
            name.Append("_recPS");
            comp.Form("TRUE");
        }
        
        out=task->CreateOutput(Form("Xip_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kXi);
        out->SetCutID(0,cut1);
        out->SetCharge(0,charge1);
        
        out->SetDaughter(1,AliRsnDaughter::kProton);
        out->SetCutID(1,iCutP);
        out->SetCharge(1,charge2);
        
        if(!pairID) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(ipdg);
        out->SetMotherMass(mass);
        
        if(k<=6){
            if(xID==imID || xID==mmID) out->AddAxis(xID,200,2.1,3.1);// axis X: invmass or resolution
            else out->AddAxis(diffID,200,-0.02,0.02);
            out->AddAxis(ptID,200,0.0,20.0);// axis Y: transverse momentum
            out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
        }else{//Phase-space histograms
            out->AddAxis(fdpt,100,0.,10.);
            out->AddAxis(sdpt,100,0.,10.);
            out->AddAxis(ptID,40,0.,20.);
        }
    }
    
    if(isMC){
        AliRsnCutMiniPair* cutEta=new AliRsnCutMiniPair("cutPseudorapidity", AliRsnCutMiniPair::kPseudorapidityRangeMC);
        cutEta->SetRangeD(-0.8,0.8);
        
        AliRsnCutSet* cutsPairDaughter=new AliRsnCutSet("pairCutsDaughter", AliRsnTarget::kMother);
        cutsPairDaughter->AddCut(cutEta);
        cutsPairDaughter->SetCutScheme(cutEta->GetName());
        
        //for efficiency of (anti)proton
        out = task->CreateOutput("Xip_Pp_mother", "HIST", "SINGLE");
        out->SetDaughter(0, AliRsnDaughter::kUnknown);
        out->SetDaughter(1, AliRsnDaughter::kUnknown);
        out->SetMotherPDG(2212);
        out->SetMotherMass(0.938272);
        out->SetPairCuts(cutsPairDaughter);
        out->AddAxis(ptID,200,0.0,20.0);
        
        out = task->CreateOutput("Xip_Pm_mother", "HIST", "SINGLE");
        out->SetDaughter(0, AliRsnDaughter::kUnknown);
        out->SetDaughter(1, AliRsnDaughter::kUnknown);
        out->SetMotherPDG(-2212);
        out->SetMotherMass(0.938272);
        out->SetPairCuts(cutsPairDaughter);
        out->AddAxis(ptID,200,0.0,20.0);
        
        //for efficiency of (anti-)Xi
        out = task->CreateOutput("Xip_Xim_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kLambda);
        out->SetMotherPDG(3312);
        out->SetMotherMass(1.32171);
        out->SetPairCuts(cutsPairDaughter);
        out->AddAxis(ptID,200,0.0,20.0);
        
        out = task->CreateOutput("Xip_Xip_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kLambda);
        out->SetMotherPDG(-3312);
        out->SetMotherMass(1.32171);
        out->SetPairCuts(cutsPairDaughter);
        out->AddAxis(ptID,200,0.0,20.0);
    }
    
    return kTRUE;
}


//=============================


Bool_t Config_XiLambda(
                  AliRsnMiniAnalysisTask *task,
                  TString     lname,
                  Bool_t      isMC,
                  Int_t       system,
                  Int_t       EventCuts,
                  Int_t       TrackCutsXi,
                  Int_t       TrackCutsLambda
                  ){
    bool isPP=false;
    if(!system) isPP=true;
    int trigger=EventCuts%10;
    int MultBins=(EventCuts/10)%10;
    if(system==1 || system==2) MultBins=1;
    
    char suffix[1000];
    sprintf(suffix,"_%s",lname.Data());
    Bool_t enableMonitor=kTRUE;
    
    Double_t mass=1.115683+1.32171;
    
    Bool_t CheckOOBP=true;
    if(TrackCutsXi==1) CheckOOBP=false;
    
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
    cutLambda->SetMaxPseudorapidity(0.8);
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
    cutAntiLambda->SetMaxPseudorapidity(0.8);
    cutAntiLambda->SetMinTPCcluster(-1);
    
    AliRsnCutSet* cutSetAntiLambda=new AliRsnCutSet("setAntiLambda",AliRsnTarget::kDaughter);
    cutSetAntiLambda->AddCut(cutAntiLambda);
    cutSetAntiLambda->SetCutScheme(cutAntiLambda->GetName());
    Int_t iCutAntiLambda=task->AddTrackCuts(cutSetAntiLambda);
    
    // selections for Xi
    Float_t XiPIDcut=3.;
    Float_t V0dDCA=1.6;
    Float_t XidDCA=1.6;
    Float_t XiMinDCA=0.07;
    Float_t Xi_massTol=0.007;
    Float_t Xi_massTolVeto=0.007;
    Float_t V0CosPoinAn=0.97;
    Float_t XiCosPoinAn=0.97;
    
    AliRsnCutCascade* cutXi=new AliRsnCutCascade("cutXi",kXiMinus);
    cutXi->SetPIDCutV0Proton(XiPIDcut);
    cutXi->SetPIDCutV0Pion(XiPIDcut);
    cutXi->SetPIDCutBachelor(XiPIDcut);
    cutXi->SetESDtrackCuts(esdTrackCuts);
    cutXi->SetV0MaxDaughtersDCA(V0dDCA);
    cutXi->SetCascadeMaxDaughtersDCA(XidDCA);
    cutXi->SetV0MaxDCAVertex(1e5); // not using
    cutXi->SetV0MinDCAVertex(XiMinDCA);
    cutXi->SetCascadeMaxDCAVertex(1e5); // not using
    cutXi->SetCascadeMinDCAVertex(-1e5); // not using
    cutXi->SetV0LowRadius(0); // not using
    cutXi->SetV0HighRadius(1e5); // not using
    cutXi->SetCascadeLowRadius(0); // not using
    cutXi->SetCascadeHighRadius(1e5); // not using
    cutXi->SetMassTolerance(Xi_massTol);
    cutXi->SetMassToleranceVeto(Xi_massTolVeto);//Rejection range for Competing Xi Rejection
    cutXi->SetSwitch(kFALSE); // not using
    cutXi->SetV0MinCosPointingAngle(V0CosPoinAn);
    cutXi->SetCascadeMinCosPointingAngle(XiCosPoinAn);
    cutXi->SetMaxPseudorapidity(0.8);
    cutXi->SetMinTPCcluster(-1);
    
    AliRsnCutCascade* cutXibar=new AliRsnCutCascade("cutXibar",kXiPlusBar);
    cutXibar->SetPIDCutV0Proton(XiPIDcut);
    cutXibar->SetPIDCutV0Pion(XiPIDcut);
    cutXibar->SetPIDCutBachelor(XiPIDcut);
    cutXibar->SetESDtrackCuts(esdTrackCuts);
    cutXibar->SetV0MaxDaughtersDCA(V0dDCA);
    cutXibar->SetCascadeMaxDaughtersDCA(XidDCA);
    cutXibar->SetV0MaxDCAVertex(1e5); // not using
    cutXibar->SetV0MinDCAVertex(XiMinDCA);
    cutXibar->SetCascadeMaxDCAVertex(1e5); // not using
    cutXibar->SetCascadeMinDCAVertex(-1e5); // not using
    cutXibar->SetV0LowRadius(0); // not using
    cutXibar->SetV0HighRadius(1e5); // not using
    cutXibar->SetCascadeLowRadius(0); // not using
    cutXibar->SetCascadeHighRadius(1e5); // not using
    cutXibar->SetMassTolerance(Xi_massTol);
    cutXibar->SetMassToleranceVeto(Xi_massTolVeto);//Rejection range for Competing Xi Rejection
    cutXibar->SetSwitch(kFALSE); // not using
    cutXibar->SetV0MinCosPointingAngle(V0CosPoinAn);
    cutXibar->SetCascadeMinCosPointingAngle(XiCosPoinAn);
    cutXibar->SetMaxPseudorapidity(0.8);
    cutXibar->SetMinTPCcluster(-1);
    
    AliRsnCutSet* cutSetXi=new AliRsnCutSet("setXi",AliRsnTarget::kDaughter);
    cutSetXi->AddCut(cutXi);
    cutSetXi->SetCutScheme(cutXi->GetName());
    Int_t icutXi=task->AddTrackCuts(cutSetXi);
    
    AliRsnCutSet* cutSetXibar=new AliRsnCutSet("setXibar",AliRsnTarget::kDaughter);
    cutSetXibar->AddCut(cutXibar);
    cutSetXibar->SetCutScheme(cutXibar->GetName());
    Int_t icutXibar=task->AddTrackCuts(cutSetXibar);
    
    // monitoring
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputCascade.C");
        
        AddMonitorOutputV0(isMC,cutSetLambda->GetMonitorOutput(), "lambda");
        AddMonitorOutputV0(isMC,cutSetAntiLambda->GetMonitorOutput(), "antilambda");
        
        AddMonitorOutputV0(isMC, cutSetXi->GetMonitorOutput(), "Lambda", "nokine");
        AddMonitorOutputV0(isMC, cutSetXibar->GetMonitorOutput(), "AntiLambda", "nokine");
        AddMonitorOutputCascade(isMC, cutSetXi->GetMonitorOutput(), "Xi");
        AddMonitorOutputCascade(isMC, cutSetXibar->GetMonitorOutput(), "AntiXi");
    }
    
    // pair cuts
    AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
    if(system!=1) cutY->SetRangeD(-0.5,0.5);
    else cutY->SetRangeD(-0.465,0.035);
    
    AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
    AliRsnCutMiniPair* cutOOBP=new AliRsnCutMiniPair("cutOOBP",AliRsnCutMiniPair::kPassesOOBPileupCut);
    
    AliRsnCutSet* cutsPairSame=new AliRsnCutSet("pairCutsSame",AliRsnTarget::kMother);
    cutsPairSame->AddCut(cutY);
    cutsPairSame->AddCut(cutV0);
    
    AliRsnCutSet* cutsPairMix=new AliRsnCutSet("pairCutsMix", AliRsnTarget::kMother);
    cutsPairMix->AddCut(cutY);
    
    if(CheckOOBP){
        cutsPairSame->AddCut(cutOOBP);
        cutsPairSame->SetCutScheme(TString::Format("%s&(!%s)&%s",cutY->GetName(),cutV0->GetName(),cutOOBP->GetName()).Data());
        cutsPairMix->AddCut(cutOOBP);
        cutsPairMix->SetCutScheme(TString::Format("%s&%s",cutY->GetName(),cutOOBP->GetName()).Data());
    }else{
        cutsPairSame->SetCutScheme(TString::Format("%s&(!%s)",cutY->GetName(),cutV0->GetName()).Data());
        cutsPairMix->SetCutScheme(cutY->GetName());
    }
    
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
    /* mother mass      */ Int_t mmID   = task->CreateValue(AliRsnMiniValue::kInvMassMother,kFALSE);
    /* IM difference    */ Int_t diffID = task->CreateValue(AliRsnMiniValue::kInvMassDiff,kTRUE);
    /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
    /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
    /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
    /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);
    /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt,kFALSE);
    /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt,kFALSE);
    
    // -- Create all needed outputs -----------------------------------------------------------------
    // use an array for more compact writing, which are different on mixing and charges
    
    Int_t i,k,xID,cut1,cut2,pairID,ipdg;
    TString name,comp;
    Char_t charge1;
    AliRsnMiniOutput* out;
    
    task->SetMotherAcceptanceCutMinPt(0.15);
    task->SetMotherAcceptanceCutMaxEta(0.8);
    task->KeepMotherInAcceptance(true);
    
    for(i=0;i<2;i++) for(j=0;j<2;j++) for(k=0;k<9;k++){
        if(!i){
            name.Form("Xim");//Xi-
            cut1=icutXi;
            charge1='-';
        }else{
            name.Form("Xip");//anti-Xi+
            cut1=icutXibar;
            charge1='+';
        }
        
        if(!j){
            name.Append("Lambdap");//Lambda
            cut2=iCutLambda;
        }else{
            name.Append("Lambdaa");//anti-Lambda
            cut2=iCutAntiLambda;
        }
        
        if(!isMC && k>=3) continue;
        pairID=1;
        xID=imID;
        if(!k){
            comp.Form("PAIR");
            pairID=0;
        }else if(k==1){
            name.Append("Mix");
            comp.Form("MIX");
        }else if(k==2){
            name.Append("Rotated");
            comp.Form("ROTATE1");
            pairID=0;
        }else if(k==3){
            name.Append("_gen");
            comp.Form("MOTHER");
        }else if(k==4){
            name.Append("_rec");
            comp.Form("TRUE");
        }else if(k==5){
            name.Append("_recMM");
            comp.Form("TRUE");
            xID=mmID;
        }else if(k==6){
            name.Append("_res");
            comp.Form("TRUE");
            xID=diffID;
        }else if(k==7){
            name.Append("_genPS");
            comp.Form("MOTHER_IN_ACC");
        }else if(k==8){
            name.Append("_recPS");
            comp.Form("TRUE");
        }
        
        ipdg=(i==0)?3324:-3324;
        
        out=task->CreateOutput(Form("XiLambda_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kXi);
        out->SetCutID(0,cut1);
        out->SetCharge(0,charge1);
        
        out->SetDaughter(1,AliRsnDaughter::kLambda);
        out->SetCutID(1,cut2);
        out->SetCharge(1,'0');
        
        if(!pairID) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(ipdg);
        out->SetMotherMass(mass);
        
        if(k<=6){
            if(xID==imID || xID==mmID) out->AddAxis(xID,220,2.4,3.5);// axis X: invmass or resolution
            else out->AddAxis(diffID,200,-0.02,0.02);
            out->AddAxis(ptID,200,0.0,20.0);// axis Y: transverse momentum
            out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
        }else{//Phase-space histograms
            out->AddAxis(fdpt,100,0.,10.);
            out->AddAxis(sdpt,100,0.,10.);
            out->AddAxis(ptID,40,0.,20.);
        }
    }
    
    if(isMC){
        AliRsnCutMiniPair* cutEta=new AliRsnCutMiniPair("cutPseudorapidity", AliRsnCutMiniPair::kPseudorapidityRangeMC);
        cutEta->SetRangeD(-0.8,0.8);
        
        AliRsnCutSet* cutsPairDaughter=new AliRsnCutSet("pairCutsDaughter", AliRsnTarget::kMother);
        cutsPairDaughter->AddCut(cutEta);
        cutsPairDaughter->SetCutScheme(cutEta->GetName());
        
        //for efficiency of (anti-)Lambda
        out = task->CreateOutput("XiLambda_Lambdap_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kProton);
        out->SetMotherPDG(3122);
        out->SetMotherMass(1.115683);
        out->SetPairCuts(cutsPairDaughter);
        out->AddAxis(ptID,200,0.0,20.0);
        
        out = task->CreateOutput("XiLambda_Lambdaa_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kProton);
        out->SetMotherPDG(-3122);
        out->SetMotherMass(1.115683);
        out->SetPairCuts(cutsPairDaughter);
        out->AddAxis(ptID,200,0.0,20.0);
        
        //for efficiency of (anti-)Xi
        out = task->CreateOutput("XiLambda_Xim_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kLambda);
        out->SetMotherPDG(3312);
        out->SetMotherMass(1.32171);
        out->SetPairCuts(cutsPairDaughter);
        out->AddAxis(ptID,200,0.0,20.0);
        
        out = task->CreateOutput("XiLambda_Xip_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kLambda);
        out->SetMotherPDG(-3312);
        out->SetMotherMass(1.32171);
        out->SetPairCuts(cutsPairDaughter);
        out->AddAxis(ptID,200,0.0,20.0);
    }
    
    return kTRUE;
}
