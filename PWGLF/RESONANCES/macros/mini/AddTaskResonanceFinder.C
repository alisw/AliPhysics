/***************************************************************************
 Anders Knospe: anders.knospe@cern.ch
 Macro to configure the resonance package for analyses using
 AliRsnMiniResonanceFinder (where one decay product is itself a resonance).
 
 ****************************************************************************/

#ifdef __CLING__
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C>
#include <PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C>
//#include <PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputCascade.C>

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

Bool_t Config_pikstar0(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
Bool_t Config_kxkstar0(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
Bool_t Config_k0kstar0(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
Bool_t Config_pkstar0(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);

Bool_t Config_pikstarx(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
Bool_t Config_kxkstarx(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
Bool_t Config_k0kstarx(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
Bool_t Config_pkstarx(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);

Bool_t Config_Kstar0Lambda(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
Bool_t Config_KstarxLambda(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);

Bool_t Config_k0kxpi(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);

Bool_t Config_KxLambdastar(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);
Bool_t Config_K0Lambdastar(AliRsnMiniAnalysisTask*,TString,Bool_t,Int_t,Int_t,Int_t,Int_t);

#endif

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
    
    //TH2F* hmc=new TH2F("MultiVsCent","", nmult,multbins, 401,ybins);
    //hmc->GetYaxis()->SetTitle("QUALITY");
    TH2F* hmc=new TH2F("TrackletsVsCent","", nmult,multbins, 401,ybins);
    hmc->GetYaxis()->SetTitle("TRACKLETS");
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
        
    }else if(d1==AliRsnDaughter::kPion && d2==AliRsnDaughter::kKstar0){
        Config_pikstar0(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
    }else if(d2==AliRsnDaughter::kPion && d1==AliRsnDaughter::kKstar0){
        Config_pikstar0(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
        
    }else if(d1==AliRsnDaughter::kKaon && d2==AliRsnDaughter::kKstar0){
        Config_kxkstar0(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
    }else if(d2==AliRsnDaughter::kKaon && d1==AliRsnDaughter::kKstar0){
        Config_kxkstar0(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
        
    }else if(d1==AliRsnDaughter::kKaon0 && d2==AliRsnDaughter::kKstar0){
        Config_k0kstar0(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
    }else if(d2==AliRsnDaughter::kKaon0 && d1==AliRsnDaughter::kKstar0){
        Config_k0kstar0(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
        
    }else if(d1==AliRsnDaughter::kProton && d2==AliRsnDaughter::kKstar0){
        Config_pkstar0(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
    }else if(d2==AliRsnDaughter::kProton && d1==AliRsnDaughter::kKstar0){
        Config_pkstar0(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
        
    }else if(d1==AliRsnDaughter::kPion && d2==AliRsnDaughter::kKstarpm){
        Config_pikstarx(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
    }else if(d2==AliRsnDaughter::kPion && d1==AliRsnDaughter::kKstarpm){
        Config_pikstarx(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
        
    }else if(d1==AliRsnDaughter::kKaon && d2==AliRsnDaughter::kKstarpm){
        Config_kxkstarx(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
    }else if(d2==AliRsnDaughter::kKaon && d1==AliRsnDaughter::kKstarpm){
        Config_kxkstarx(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
        
    }else if(d1==AliRsnDaughter::kKaon0 && d2==AliRsnDaughter::kKstarpm){
        Config_k0kstarx(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
    }else if(d2==AliRsnDaughter::kKaon0 && d1==AliRsnDaughter::kKstarpm){
        Config_k0kstarx(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
        
    }else if(d1==AliRsnDaughter::kProton && d2==AliRsnDaughter::kKstarpm){
        Config_pkstarx(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
    }else if(d2==AliRsnDaughter::kProton && d1==AliRsnDaughter::kKstarpm){
        Config_pkstarx(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
        
    }else if(d1==AliRsnDaughter::kLambda && d2==AliRsnDaughter::kKstar0){
        Config_Kstar0Lambda(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
    }else if(d2==AliRsnDaughter::kLambda && d1==AliRsnDaughter::kKstar0){
        Config_Kstar0Lambda(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
        
    }else if(d1==AliRsnDaughter::kLambda && d2==AliRsnDaughter::kKstarpm){
        Config_KstarxLambda(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
    }else if(d2==AliRsnDaughter::kLambda && d1==AliRsnDaughter::kKstarpm){
        Config_KstarxLambda(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
        
    }else if(d1==AliRsnDaughter::kKaon0 && d2==AliRsnDaughter::kKaon){
        Config_k0kxpi(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
    }else if(d2==AliRsnDaughter::kKaon0 && d1==AliRsnDaughter::kKaon){
        Config_k0kxpi(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
     
     
    }else if(d1==AliRsnDaughter::kLambdastar && d2==AliRsnDaughter::kKaon0){
        Config_K0Lambdastar(task,lname,isMC,system,EventCuts,TrackCuts1,TrackCuts2);
    }else if(d2==AliRsnDaughter::kLambdastar && d1==AliRsnDaughter::kKaon0){
        Config_K0Lambdastar(task,lname,isMC,system,EventCuts,TrackCuts2,TrackCuts1);
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
    Float_t k0s_massTolSigma=4;
    Int_t   k0s_massTolID=0;
    Float_t k0s_massTol=0.03;
    Float_t k0s_massTolVeto=0.004;
    Bool_t  k0sSwitch=kFALSE;
    Float_t k0sCosPoinAn=0.97;
    
    if(TrackCutsK0==1) k0s_massTolID=1;//use pT-dependent mass tolerance cut
    
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
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetK0s->GetMonitorOutput(), "K0S");
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
    AliRsnCutMiniPair* cutYPhi=new AliRsnCutMiniPair("cutPseudorapidityPhi",AliRsnCutMiniPair::kPseudorapidityRange);
    cutYPhi->SetRangeD(-0.8,0.8);
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
    if(isMC) finder[i]->SetPairMode(1);
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
    
    Int_t xID,cut2,pairID,ipdg;
    TString name,comp;
    Char_t charge1;
    AliRsnMiniOutput* out;
    
    task->SetMotherAcceptanceCutMinPt(0.15);
    task->SetMotherAcceptanceCutMaxEta(0.8);
    task->KeepMotherInAcceptance(true);
    
    for(i=0;i<2;i++) for(j=0;j<12;j++){
        if(!i){
            name.Form("Pp");
            charge1='+';
            ipdg=9322132;
        }else{
            name.Form("Pm");
            charge1='-';
            ipdg=-9322132;
        }
        
        xID=imID;
        cut2=0;
        pairID=0;
        if(!j){
            name.Append("Phi");
            comp.Form("PAIR");
        }else if(j==1){
            name.Append("KpKp");
            comp.Form("PAIR");
            cut2=1;
        }else if(j==2){
            name.Append("KmKm");
            comp.Form("PAIR");
            cut2=2;
        }else if(j==3){
            name.Append("SB");
            comp.Form("PAIR");
            cut2=3;
        }else if(j==4){
            name.Append("PhiMix");
            comp.Form("MIX");
            pairID=1;
        }else if(j==5){
            name.Append("PhiRotated");
            comp.Form("ROTATE1");
        }else if(j==6){
            if(!isMC) continue;
            name.Append("Phi_gen");
            comp.Form("MOTHER");
            pairID=1;
        }else if(j==7){
            if(!isMC) continue;
            name.Append("Phi_rec");
            comp.Form("TRUE");
            pairID=1;
        }else if(j==8){
            if(!isMC) continue;
            name.Append("Phi_recMM");
            comp.Form("TRUE");
            xID=mmID;
            pairID=1;
        }else if(j==9){
            if(!isMC) continue;
            name.Append("Phi_res");
            comp.Form("TRUE");
            xID=diffID;
            pairID=1;
        }else if(j==10){
            if(!isMC) continue;
            name.Append("Phi_genPS");
            comp.Form("MOTHER_IN_ACC");
            pairID=1;
        }else if(j==11){
            if(!isMC) continue;
            name.Append("Phi_recPS");
            comp.Form("TRUE");
            pairID=1;
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
        
        if(j<=9){
            if(xID==imID || xID==mmID) out->AddAxis(xID,210,1.95,3.);// axis X: invmass or resolution
            else out->AddAxis(diffID,200,-0.02,0.02);
            out->AddAxis(ptID,200,0.0,20.0);// axis Y: transverse momentum
            out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
        }else{//Phase-space histograms
            out->AddAxis(fdpt,100,0.,10.);
            out->AddAxis(sdpt,100,0.,10.);
            out->AddAxis(ptID,40,0.,20.);
        }
    }
    
    // fill monitoring histogram for the resonance (phi)
    for(i=0;i<7;i++){
        xID=imID;
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
        }else if(i==6){
            if(!isMC) continue;
            name.Form("phimass_res");
            comp.Form("TRUE");
            cut2=0;
            xID=diffID;
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
        
        if(xID==imID) out->AddAxis(imID,70,1.,1.07);
        else out->AddAxis(diffID,200,-0.02,0.02);
        out->AddAxis(ptID,200,0.0,20.0);
        out->AddAxis(centID,nmult,multbins);
    }
    
    if(isMC){
        //for efficiency of (anti)proton
        AliRsnCutMiniPair* cutetaP=new AliRsnCutMiniPair("cutPseudorapidityP", AliRsnCutMiniPair::kPseudorapidityRangeMC);
        cutetaP->SetRangeD(-0.8,0.8);
        
        AliRsnCutSet* cutsPairP=new AliRsnCutSet("pairCutsP", AliRsnTarget::kMother);
        cutsPairP->AddCut(cutetaP);
        cutsPairP->SetCutScheme(cutetaP->GetName());
        
        out = task->CreateOutput("pphi_Pp_mother", "HIST", "SINGLE");
        out->SetDaughter(0, AliRsnDaughter::kUnknown);
        out->SetDaughter(1, AliRsnDaughter::kUnknown);
        out->SetMotherPDG(2212);
        out->SetMotherMass(0.938272);
        out->SetPairCuts(cutsPairP);
        out->AddAxis(ptID,200,0.0,20.0);
        
        out = task->CreateOutput("pPhi_pm_mother", "HIST", "SINGLE");
        out->SetDaughter(0, AliRsnDaughter::kUnknown);
        out->SetDaughter(1, AliRsnDaughter::kUnknown);
        out->SetMotherPDG(-2212);
        out->SetMotherMass(0.938272);
        out->SetPairCuts(cutsPairP);
        out->AddAxis(ptID,200,0.0,20.0);
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
    if(!cutSetK){cerr<<"Error in AddTaskResonanceFinder::Config_phiphi(): missing cutSetK"<<endl; return kFALSE;}
    
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

        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetLambda->GetMonitorOutput(), "lambda");
        AddMonitorOutputV0(isMC,cutSetAntiLambda->GetMonitorOutput(), "antilambda");
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
    Int_t pairRotate=(TrackCutsK/100000)%10;
    
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
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetLambda->GetMonitorOutput(), "lambda");
        AddMonitorOutputV0(isMC,cutSetAntiLambda->GetMonitorOutput(), "antilambda");
    }
    
    // AliRsnMiniResonanceFinder
    
    AliRsnMiniResonanceFinder* finder[8];
    Int_t i,iCutSig[8];
    
    AliRsnCutMiniPair* cutMassS=new AliRsnCutMiniPair("cutMassSigmastar",AliRsnCutMiniPair::kMassRange);
    cutMassS->SetRangeD(1.346,1.427);
    AliRsnCutMiniPair* cutYS=new AliRsnCutMiniPair("cutPseudorapiditySigmastar",AliRsnCutMiniPair::kPseudorapidityRange);
    cutYS->SetRangeD(-0.8,0.8);
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
    if(isMC) finder[i]->SetPairMode(1);
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
    if(isMC) finder[i]->SetPairMode(1);
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
    if(isMC) finder[i]->SetPairMode(1);
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
    if(isMC) finder[i]->SetPairMode(1);
    iCutSig[i]=task->AddResonanceFinder(finder[i]);
    
    // sidebands
    AliRsnCutMiniPair* cutMassSB=new AliRsnCutMiniPair("cutMassSB",AliRsnCutMiniPair::kMassRange);
    cutMassSB->SetRangeD(1.446,1.486);
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
    if(system!=1) cutY->SetRangeD(-0.5,0.5);
    else cutY->SetRangeD(-0.465,0.035);
    
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
    
    Int_t k,xID,cut2,pairID,ipdg=3124;
    AliRsnDaughter::ESpecies d2=AliRsnDaughter::kUnknown;
    TString name,comp;
    Char_t charge1,charge2;
    Double_t rmass;
    AliRsnMiniOutput* out;
    
    task->SetMotherAcceptanceCutMinPt(0.15);
    task->SetMotherAcceptanceCutMaxEta(0.8);
    task->KeepMotherInAcceptance(true);
    
    for(i=0;i<2;i++) for(j=0;j<4;j++) for(k=0;k<10;k++){
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
        
        xID=imID;
        pairID=1;
        if(!k){
            comp.Form("PAIR");
            pairID=0;
        }else if(k==1){
            name.Append("SB");
            comp.Form("PAIR");
            cut2+=4;
            pairID=0;
        }else if(k==2){
            name.Append("Mix");
            comp.Form("MIX");
        }else if(k==3){
            name.Append("ROTATE");
            if(!pairRotate) {comp.Form("ROTATE1");}
            else {comp.Form("ROTATE2");}
            pairID=0;
        }else if(k==4){
            if(!isMC) continue;
            name.Append("_gen");
            comp.Form("MOTHER");
        }else if(k==5){
            if(!isMC) continue;
            name.Append("_rec");
            comp.Form("TRUE");
        }else if(k==6){
            if(!isMC) continue;
            name.Append("_recMM");
            comp.Form("TRUE");
            xID=mmID;
        }else if(k==7){
            if(!isMC) continue;
            name.Append("_res");
            comp.Form("TRUE");
            xID=diffID;
        }else if(k==8){
            if(!isMC) continue;
            name.Append("_genPS");
            comp.Form("MOTHER_IN_ACC");
        }else if(k==9){
            if(!isMC) continue;
            name.Append("_recPS");
            comp.Form("TRUE");
        }
        
        out=task->CreateOutput(Form("kxSigmastar_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kKaon);
        out->SetCutID(0,iCutK);
        out->SetCharge(0,charge1);
        
        out->SetDaughter(1,d2);
        out->SetCutID(1,iCutSig[cut2]);
        out->SetCharge(1,charge2);
        if(k!=1) out->SetUseStoredMass(1);
            
        if(!pairID) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(ipdg);
        out->SetMotherMass(mass);
            
        if(k<=7){
            if(xID==imID || xID==mmID) out->AddAxis(xID,240,1.8,3);// axis X: invmass or resolution
            else out->AddAxis(diffID,200,-0.02,0.02);
            out->AddAxis(ptID,50,0.0,20.0);// axis Y: transverse momentum
            out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
        }else{//Phase-space histograms
            out->AddAxis(fdpt,100,0.,10.);
            out->AddAxis(sdpt,100,0.,10.);
            out->AddAxis(ptID,40,0.,20.);
        }
    }
    
    // fill monitoring histogram for the resonance (Sigma*)
    for(i=0;i<4;i++) for(j=0;j<5;j++){
        if(!i) name.Form("Sigmastarpp");
        else if(i==1) name.Form("Sigmastarma");
        else if(i==2) name.Form("Sigmastarmp");
        else if(i==3) name.Form("Sigmastarpa");
        
        xID=imID;
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
        }else if(j==4){
            if(!isMC) continue;
            name.Append("_recres");
            comp.Form("TRUE");
            xID=diffID;
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
        
        if(xID==imID) out->AddAxis(imID,150,1.3,1.6);
        else out->AddAxis(diffID,200,-0.02,0.02);
        out->AddAxis(ptID,50,0.0,20.0);
        out->AddAxis(centID,nmult,multbins);
    }
    
    if(isMC){
        //for efficiency of K+/-
        AliRsnCutMiniPair* cutetaK=new AliRsnCutMiniPair("cutPseudorapidityK", AliRsnCutMiniPair::kPseudorapidityRangeMC);
        cutetaK->SetRangeD(-0.8,0.8);
        
        AliRsnCutSet* cutsPairK=new AliRsnCutSet("pairCutsK", AliRsnTarget::kMother);
        cutsPairK->AddCut(cutetaK);
        cutsPairK->SetCutScheme(cutetaK->GetName());
        
        out = task->CreateOutput("kxSigmastar_Kp_mother", "HIST", "SINGLE");
        out->SetDaughter(0, AliRsnDaughter::kUnknown);
        out->SetDaughter(1, AliRsnDaughter::kUnknown);
        out->SetMotherPDG(321);
        out->SetMotherMass(0.493677);
        out->SetPairCuts(cutsPairK);
        out->AddAxis(ptID,200,0.0,20.0);
        
        out = task->CreateOutput("kxSigmastar_Km_mother", "HIST", "SINGLE");
        out->SetDaughter(0, AliRsnDaughter::kUnknown);
        out->SetDaughter(1, AliRsnDaughter::kUnknown);
        out->SetMotherPDG(-321);
        out->SetMotherMass(0.493677);
        out->SetPairCuts(cutsPairK);
        out->AddAxis(ptID,200,0.0,20.0);
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
    
    Int_t K0sCuts=TrackCutsK;
    
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
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetK0s->GetMonitorOutput(), "K0S");
        AddMonitorOutputV0(isMC,cutSetLambda->GetMonitorOutput(), "lambda");
        AddMonitorOutputV0(isMC,cutSetAntiLambda->GetMonitorOutput(), "antilambda");
    }
    
    // AliRsnMiniResonanceFinder
    
    AliRsnMiniResonanceFinder* finder[8];
    Int_t i,iCutSig[8];
    
    AliRsnCutMiniPair* cutMassS=new AliRsnCutMiniPair("cutMassSigmastar",AliRsnCutMiniPair::kMassRange);
    cutMassS->SetRangeD(1.346,1.427);
    AliRsnCutMiniPair* cutYS=new AliRsnCutMiniPair("cutPseudorapiditySigmastar",AliRsnCutMiniPair::kPseudorapidityRange);
    cutYS->SetRangeD(-0.8,0.8);
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
    if(isMC) finder[i]->SetPairMode(1);
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
    if(isMC) finder[i]->SetPairMode(1);
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
    if(isMC) finder[i]->SetPairMode(1);
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
    if(isMC) finder[i]->SetPairMode(1);
    iCutSig[i]=task->AddResonanceFinder(finder[i]);
    
    // sidebands
    AliRsnCutMiniPair* cutMassSB=new AliRsnCutMiniPair("cutMassSB",AliRsnCutMiniPair::kMassRange);
    cutMassSB->SetRangeD(1.446,1.486);
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
    if(system!=1) cutY->SetRangeD(-0.5,0.5);
    else cutY->SetRangeD(-0.465,0.035);
    
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
    
    Int_t k,xID,cut2,pairID,ipdg;
    AliRsnDaughter::ESpecies d2=AliRsnDaughter::kUnknown;
    TString name,comp;
    Char_t charge1='0',charge2;
    AliRsnMiniOutput* out;
    
    task->SetMotherAcceptanceCutMinPt(0.15);
    task->SetMotherAcceptanceCutMaxEta(0.8);
    task->KeepMotherInAcceptance(true);
    
    for(j=0;j<4;j++) for(k=0;k<10;k++){
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
        
        xID=imID;
        pairID=1;
        if(!k){
            comp.Form("PAIR");
            pairID=0;
        }else if(k==1){
            name.Append("SB");
            comp.Form("PAIR");
            cut2+=4;
            pairID=0;
        }else if(k==2){
            name.Append("Mix");
            comp.Form("MIX");
        }else if(k==3){
            name.Append("ROTATE");
            comp.Form("ROTATE1");
            pairID=0;
        }else if(k==4){
            if(!isMC) continue;
            name.Append("_gen");
            comp.Form("MOTHER");
        }else if(k==5){
            if(!isMC) continue;
            name.Append("_rec");
            comp.Form("TRUE");
        }else if(k==6){
            if(!isMC) continue;
            name.Append("_recMM");
            comp.Form("TRUE");
            xID=mmID;
        }else if(k==7){
            if(!isMC) continue;
            name.Append("_res");
            comp.Form("TRUE");
            xID=diffID;
        }else if(k==8){
            if(!isMC) continue;
            name.Append("_genPS");
            comp.Form("MOTHER_IN_ACC");
        }else if(k==9){
            if(!isMC) continue;
            name.Append("_recPS");
            comp.Form("TRUE");
        }
        
        out=task->CreateOutput(Form("k0Sigmastar_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kKaon0);
        out->SetCutID(0,iCutK0s);
        out->SetCharge(0,charge1);
        
        out->SetDaughter(1,d2);
        out->SetCutID(1,iCutSig[cut2]);
        out->SetCharge(1,charge2);
        if(k!=1) out->SetUseStoredMass(1);
        
        if(!pairID) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(3124);
        out->SetMotherMass(mass);
        
        if(k<=7){
            if(xID==imID || xID==mmID) out->AddAxis(xID,240,1.8,3);// axis X: invmass or resolution
            else out->AddAxis(diffID,200,-0.02,0.02);
            out->AddAxis(ptID,50,0.0,20.0);// axis Y: transverse momentum
            out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
        }else{//Phase-space histograms
            out->AddAxis(fdpt,100,0.,10.);
            out->AddAxis(sdpt,100,0.,10.);
            out->AddAxis(ptID,40,0.,20.);
        }
    }
    
    // fill monitoring histogram for the resonance (Sigma*)
    for(i=0;i<4;i++) for(j=0;j<5;j++){
        if(!i) name.Form("Sigmastarpp");
        else if(i==1) name.Form("Sigmastarma");
        else if(i==2) name.Form("Sigmastarmp");
        else if(i==3) name.Form("Sigmastarpa");
        
        xID=imID;
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
        }else if(j==4){
            if(!isMC) continue;
            name.Append("_recres");
            comp.Form("TRUE");
            xID=diffID;
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
        
        if(xID==imID) out->AddAxis(imID,150,1.3,1.6);
        else out->AddAxis(diffID,200,-0.02,0.02);
        out->AddAxis(ptID,50,0.0,20.0);
        out->AddAxis(centID,nmult,multbins);
    }
    
    //for efficiency of K0S
    if(isMC){
        AliRsnCutMiniPair* cutEtaV0=new AliRsnCutMiniPair("cutPseudorapidityV0", AliRsnCutMiniPair::kPseudorapidityRange);
        cutEtaV0->SetRangeD(-0.8,0.8);
        
        AliRsnCutSet* cutsPairV0=new AliRsnCutSet("pairCutsV0", AliRsnTarget::kMother);
        cutsPairV0->AddCut(cutEtaV0);
        cutsPairV0->SetCutScheme(cutEtaV0->GetName());
        
        out = task->CreateOutput("k0Sigmastar_K0S_mother", "HIST", "MOTHER");
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
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetLambda->GetMonitorOutput(), "lambda");
        AddMonitorOutputV0(isMC,cutSetAntiLambda->GetMonitorOutput(), "antilambda");
    }
    
    // AliRsnMiniResonanceFinder - K*0
    AliRsnMiniResonanceFinder* Kfinder[4];
    Int_t i,iCutKstar[4];
    
    AliRsnCutMiniPair* cutYRes=new AliRsnCutMiniPair("cutRapidityResonance",AliRsnCutMiniPair::kPseudorapidityRange);
    cutYRes->SetRangeD(-0.8,0.8);
    
    AliRsnCutMiniPair* cutMassKstar=new AliRsnCutMiniPair("cutMassKstar",AliRsnCutMiniPair::kMassRange);
    if(!SidebandKstar) cutMassKstar->SetRangeD(0.841,0.943);
    else if(SidebandKstar==1) cutMassKstar->SetRangeD(0.766,0.816);
    else cutMassKstar->SetRangeD(0.966,1.014);
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
    if(isMC && !SidebandKstar) Kfinder[i]->SetPairMode(1);
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
    if(isMC && !SidebandKstar) Kfinder[i]->SetPairMode(1);
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
    cutMassS->SetRangeD(1.346,1.427);
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
    if(isMC) Sfinder[i]->SetPairMode(1);
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
    if(isMC) Sfinder[i]->SetPairMode(1);
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
    if(isMC) Sfinder[i]->SetPairMode(1);
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
    if(isMC) Sfinder[i]->SetPairMode(1);
    iCutSig[i]=task->AddResonanceFinder(Sfinder[i]);
    
    // sidebands
    AliRsnCutMiniPair* cutMassSigmastarSB=new AliRsnCutMiniPair("cutMassSigmastarSB",AliRsnCutMiniPair::kMassRange);
    cutMassSigmastarSB->SetRangeD(1.446,1.486);
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
    if(system!=1) cutY->SetRangeD(-0.5,0.5);
    else cutY->SetRangeD(-0.465,0.035);
    
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
    
    Int_t k,xID,cut1,cut2,pairID,ipdg;
    AliRsnDaughter::ESpecies d2=AliRsnDaughter::kUnknown;
    TString name,comp;
    Char_t charge2;
    AliRsnMiniOutput* out;
    
    task->SetMotherAcceptanceCutMinPt(0.15);
    task->SetMotherAcceptanceCutMaxEta(0.8);
    task->KeepMotherInAcceptance(true);
    
    for(i=0;i<4;i++) for(j=0;j<4;j++) for(k=0;k<10;k++){
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
        
        xID=imID;
        pairID=1;
        if(!k){
            comp.Form("PAIR");
            pairID=0;
        }else if(k==1){
            name.Append("SB");
            comp.Form("PAIR");
            cut2+=4;
            pairID=0;
        }else if(k==2){
            name.Append("Mix");
            comp.Form("MIX");
        }else if(k==3){
            name.Append("ROTATE");
            comp.Form("ROTATE1");
            pairID=0;
        }else if(k==4){
            if(!isMC) continue;
            name.Append("_gen");
            comp.Form("MOTHER");
        }else if(k==5){
            if(!isMC) continue;
            name.Append("_rec");
            comp.Form("TRUE");
        }else if(k==6){
            if(!isMC) continue;
            name.Append("_recMM");
            comp.Form("TRUE");
            xID=mmID;
        }else if(k==7){
            if(!isMC) continue;
            name.Append("_res");
            comp.Form("TRUE");
            xID=diffID;
        }else if(k==8){
            if(!isMC) continue;
            name.Append("_genPS");
            comp.Form("MOTHER_IN_ACC");
        }else if(k==9){
            if(!isMC) continue;
            name.Append("_recPS");
            comp.Form("TRUE");
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
        
        if(!pairID) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(3124);
        out->SetMotherMass(mass);
        
        if(k<=7){
            if(xID==imID || xID==mmID) out->AddAxis(xID,280,2.1,3.5);// axis X: invmass or resolution
            else out->AddAxis(diffID,200,-0.02,0.02);
            out->AddAxis(ptID,50,0.0,20.0);// axis Y: transverse momentum
            out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
        }else{//Phase-space histograms
            out->AddAxis(fdpt,100,0.,10.);
            out->AddAxis(sdpt,100,0.,10.);
            out->AddAxis(ptID,40,0.,20.);
        }
    }
    
    // fill monitoring histogram for the K*0
    for(i=0;i<4;i++) for(j=0;j<4;j++){
        if(!i) name.Form("Kstar0p");
        else if(i==1) name.Form("Kstar0a");
        else if(i==2) name.Form("KpPip");
        else if(i==3) name.Form("KmPim");
        
        xID=imID;
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
        }else if(j==3){
            if(!isMC) continue;
            name.Append("_recres");
            comp.Form("TRUE");
            xID=diffID;
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
        
        if(xID==imID) out->AddAxis(imID,170,0.75,1.09);
        else out->AddAxis(diffID,200,-0.02,0.02);
        out->AddAxis(ptID,50,0.0,20.0);
        out->AddAxis(centID,nmult,multbins);
    }
    
    // fill monitoring histogram for the Sigma*
    for(i=0;i<4;i++) for(j=0;j<5;j++){
        if(!i) name.Form("Sigmastarpp");
        else if(i==1) name.Form("Sigmastarma");
        else if(i==2) name.Form("Sigmastarmp");
        else if(i==3) name.Form("Sigmastarpa");
        
        xID=imID;
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
        }else if(j==4){
            if(!isMC) continue;
            name.Append("_recres");
            comp.Form("TRUE");
            xID=diffID;
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
        
        if(xID==imID) out->AddAxis(imID,150,1.3,1.6);
        else out->AddAxis(diffID,200,-0.02,0.02);
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
    Int_t K0sCuts=TrackCutsK%10000;
    Int_t pairRotate=(TrackCutsK/100000)%10;
    
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
    cutK0s->SetMaxRapidity(2.);
    cutK0s->SetMinTPCcluster(-1);
    
    AliRsnCutSet* cutSetK0s=new AliRsnCutSet("setK0s",AliRsnTarget::kDaughter);
    cutSetK0s->AddCut(cutK0s);
    cutSetK0s->SetCutScheme(cutK0s->GetName());
    Int_t iCutK0s=task->AddTrackCuts(cutSetK0s);
    
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
    TString pname="k0s";
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        //AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetK0s->GetMonitorOutput(), "K0S");
        AddMonitorOutputV0(isMC,cutSetLambda->GetMonitorOutput(), "lambda");
        AddMonitorOutputV0(isMC,cutSetAntiLambda->GetMonitorOutput(), "antilambda");
    }
    
    // AliRsnMiniResonanceFinder - K*
    AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
    AliRsnCutMiniPair* cutYRes=new AliRsnCutMiniPair("cutPseudorapidityResonance",AliRsnCutMiniPair::kPseudorapidityRange);
    cutYRes->SetRangeD(-0.8,0.8);
    
    AliRsnMiniResonanceFinder* Kfinder[4];
    Int_t i,iCutKstar[4];
    
    AliRsnCutMiniPair* cutMassKstar=new AliRsnCutMiniPair("cutMassKstar",AliRsnCutMiniPair::kMassRange);
    cutMassKstar->SetRangeD(0.841,0.943);
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
    if(isMC) Kfinder[i]->SetPairMode(1);
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
    if(isMC) Kfinder[i]->SetPairMode(1);
    iCutKstar[i]=task->AddResonanceFinder(Kfinder[i]);
    
    AliRsnCutMiniPair* cutMassKstarSB=new AliRsnCutMiniPair("cutMassKstarSB",AliRsnCutMiniPair::kMassRange);
    if(SidebandKstar<=1) cutMassKstarSB->SetRangeD(0.766,0.816);
    else cutMassKstarSB->SetRangeD(0.966,1.014);
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
    cutMassS->SetRangeD(1.346,1.427);
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
    if(isMC) Sfinder[i]->SetPairMode(1);
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
    if(isMC) Sfinder[i]->SetPairMode(1);
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
    if(isMC) Sfinder[i]->SetPairMode(1);
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
    if(isMC) Sfinder[i]->SetPairMode(1);
    iCutSig[i]=task->AddResonanceFinder(Sfinder[i]);
    
    // sidebands
    AliRsnCutMiniPair* cutMassSigmastarSB=new AliRsnCutMiniPair("cutMassSigmastarSB",AliRsnCutMiniPair::kMassRange);
    cutMassSigmastarSB->SetRangeD(1.446,1.486);
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
    if(system!=1) cutY->SetRangeD(-0.5,0.5);
    else cutY->SetRangeD(-0.465,0.035);
    
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
    
    Int_t k,xID,cut1,cut2,pairID,ipdg;
    AliRsnDaughter::ESpecies d2=AliRsnDaughter::kUnknown;
    TString name,comp;
    Char_t charge1,charge2;
    AliRsnMiniOutput* out;
    
    task->SetMotherAcceptanceCutMinPt(0.15);
    task->SetMotherAcceptanceCutMaxEta(0.8);
    task->KeepMotherInAcceptance(true);
    
    for(i=0;i<4;i++) for(j=0;j<4;j++) for(k=0;k<10;k++){
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
        
        xID=imID;
        pairID=1;
        if(!k){
            comp.Form("PAIR");
            pairID=0;
        }else if(k==1){
            name.Append("SB");
            comp.Form("PAIR");
            cut2+=4;
            pairID=0;
        }else if(k==2){
            name.Append("Mix");
            comp.Form("MIX");
        }else if(k==3){
            name.Append("ROTATE");
            if(!pairRotate) {comp.Form("ROTATE1");}
            else {comp.Form("ROTATE2");}
            pairID=0;
        }else if(k==4){
            if(!isMC) continue;
            name.Append("_gen");
            comp.Form("MOTHER");
        }else if(k==5){
            if(!isMC) continue;
            name.Append("_rec");
            comp.Form("TRUE");
        }else if(k==6){
            if(!isMC) continue;
            name.Append("_recMM");
            comp.Form("TRUE");
            xID=mmID;
        }else if(k==7){
            if(!isMC) continue;
            name.Append("_res");
            comp.Form("TRUE");
            xID=diffID;
        }else if(k==8){
            if(!isMC) continue;
            name.Append("_genPS");
            comp.Form("MOTHER_IN_ACC");
        }else if(k==9){
            if(!isMC) continue;
            name.Append("_recPS");
            comp.Form("TRUE");
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
        
        if(!pairID) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(3124);
        out->SetMotherMass(mass);
        
        if(k<=7){
            if(xID==imID || xID==mmID) out->AddAxis(xID,280,2.1,3.5);// axis X: invmass or resolution
            else out->AddAxis(diffID,200,-0.02,0.02);
            out->AddAxis(ptID,50,0.0,20.0);// axis Y: transverse momentum
            out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
        }else{//Phase-space histograms
            out->AddAxis(fdpt,100,0.,10.);
            out->AddAxis(sdpt,100,0.,10.);
            out->AddAxis(ptID,40,0.,20.);
        }
    }
    
    // fill monitoring histogram for the K*
    for(i=0;i<4;i++) for(j=0;j<4;j++){
        if(!i) name.Form("Kstarp");
        else if(i==1) name.Form("Kstarm");
        else if(i==2) name.Form("KstarpSB");
        else if(i==3) name.Form("KstarmSB");
        
        xID=imID;
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
        }else if(j==3){
            if(!isMC) continue;
            name.Append("_recres");
            comp.Form("TRUE");
            xID=diffID;
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
        
        if(xID==imID) out->AddAxis(imID,170,0.75,1.09);
        else out->AddAxis(diffID,200,-0.02,0.02);
        out->AddAxis(ptID,50,0.0,20.0);
        out->AddAxis(centID,nmult,multbins);
    }
    
    // fill monitoring histogram for the Sigma*
    for(i=0;i<4;i++) for(j=0;j<5;j++){
        if(!i) name.Form("Sigmastarpp");
        else if(i==1) name.Form("Sigmastarma");
        else if(i==2) name.Form("Sigmastarmp");
        else if(i==3) name.Form("Sigmastarpa");
        
        xID=imID;
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
        }else if(j==4){
            if(!isMC) continue;
            name.Append("_recres");
            comp.Form("TRUE");
            xID=diffID;
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
        
        if(xID==imID) out->AddAxis(imID,150,1.3,1.6);
        else out->AddAxis(diffID,200,-0.02,0.02);
        out->AddAxis(ptID,50,0.0,20.0);
        out->AddAxis(centID,nmult,multbins);
    }
    
    return kTRUE;
}


//=============================


Bool_t Config_pikstar0(
                       AliRsnMiniAnalysisTask *task,
                       TString     lname,
                       Bool_t      isMC,
                       Int_t       system,
                       Int_t       EventCuts,
                       Int_t       TrackCutsPi,
                       Int_t       TrackCutsKstar0 //kstar0
){
    bool isPP=false;
    if(!system) isPP=true;
    int trigger=EventCuts%10;
    int MultBins=(EventCuts/10)%10;
    if(system==1 || system==2) MultBins=1;
    
    char suffix[1000];
    sprintf(suffix,"_%s",lname.Data());
    Bool_t enableMonitor=kTRUE;
    
    Double_t mass=0.139571+0.895810;//895.81 from PDG 2012 //modified
    
    // set daughter cuts :
    //for all pions
    if(!(TrackCutsPi%10000)) TrackCutsPi+=3020;//default settings
    Float_t nsigmaPiTPC=0.1*(TrackCutsPi%100);
    Float_t nsigmaPiTOF=0.1*((TrackCutsPi/100)%100);
    Int_t CutTypePi=(TrackCutsPi/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only
    
    //for all Kaons
    if(!(TrackCutsKstar0%10000)) TrackCutsKstar0+=3020;//default settings
    Float_t nsigmaKTPC=0.1*(TrackCutsKstar0%100);
    Float_t nsigmaKTOF=0.1*((TrackCutsKstar0/100)%100);
    Int_t CutTypeK=(TrackCutsKstar0/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only
    ////////////////////////////////////
    
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
    if(!cutSetPi){cerr<<"Error in AddTaskResonanceFinder::Config_pikstar0(): missing cutSetPi"<<endl; return kFALSE;}
    
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
    if(!cutSetK){cerr<<"Error in AddTaskResonanceFinder::Config_pikstar0(): missing cutSetK"<<endl; return kFALSE;}
    
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
    
    // AliRsnMiniResonanceFinder//kstar0
    AliRsnCutMiniPair* cutMassKstar0=new AliRsnCutMiniPair("cutMassKstar0",AliRsnCutMiniPair::kMassRange);
    //cutMassKstar0->SetRangeD(0.85,0.95);
    cutMassKstar0->SetRangeD(0.841,0.943);
    AliRsnCutMiniPair* cutYKstar0=new AliRsnCutMiniPair("cutRapidityKstar0",AliRsnCutMiniPair::kRapidityRange);
    cutYKstar0->SetRangeD(-0.6,0.6);
    AliRsnCutSet* cutsKstar0=new AliRsnCutSet("pairCutsKstar0",AliRsnTarget::kMother);
    cutsKstar0->AddCut(cutMassKstar0);
    cutsKstar0->AddCut(cutYKstar0);
    cutsKstar0->SetCutScheme(TString::Format("%s&%s",cutMassKstar0->GetName(),cutYKstar0->GetName()).Data());
    
    //kstar0
    AliRsnMiniResonanceFinder* rsnfinder=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder",task->GetName()));
    rsnfinder->SetCutID(0,iCutK);
    rsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
    rsnfinder->SetCharge(0,'+');
    rsnfinder->SetCutID(1,iCutPi);
    rsnfinder->SetDaughter(1,AliRsnDaughter::kPion);
    rsnfinder->SetCharge(1,'-');
    rsnfinder->SetResonanceMass(0.895810);//kstar0
    rsnfinder->SetResonancePDG(313);
    rsnfinder->SetPairCuts(cutsKstar0);
    Int_t iCutKstar0=task->AddResonanceFinder(rsnfinder);
    
    //anti-kstar0
    AliRsnMiniResonanceFinder* antirsnfinder=new AliRsnMiniResonanceFinder(Form("%s_AntiResonanceFinder",task->GetName()));
    antirsnfinder->SetCutID(0,iCutK);
    antirsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
    antirsnfinder->SetCharge(0,'-');
    antirsnfinder->SetCutID(1,iCutPi);
    antirsnfinder->SetDaughter(1,AliRsnDaughter::kPion);
    antirsnfinder->SetCharge(1,'+');
    antirsnfinder->SetResonanceMass(0.895810);
    antirsnfinder->SetResonancePDG(-313);
    antirsnfinder->SetPairCuts(cutsKstar0);
    Int_t iCutantiKstar0=task->AddResonanceFinder(antirsnfinder);
    
    /////////////////////////////////
    //like sign K+ pi+
    AliRsnMiniResonanceFinder* likepprsnfinder=new AliRsnMiniResonanceFinder(Form("%s_LikeppResonanceFinder",task->GetName()));//likesign k plus pion plus
    likepprsnfinder->SetCutID(0,iCutK);
    likepprsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
    likepprsnfinder->SetCharge(0,'+');
    likepprsnfinder->SetCutID(1,iCutPi);
    likepprsnfinder->SetDaughter(1,AliRsnDaughter::kPion);
    likepprsnfinder->SetCharge(1,'+');
    likepprsnfinder->SetResonanceMass(0.895810);
    likepprsnfinder->SetResonancePDG(313);
    likepprsnfinder->SetPairCuts(cutsKstar0);
    Int_t iCutlikeppKstar0=task->AddResonanceFinder(likepprsnfinder);
    
    //like sign K- pi-
    AliRsnMiniResonanceFinder* likemmrsnfinder=new AliRsnMiniResonanceFinder(Form("%s_LikemmResonanceFinder",task->GetName()));//likesign k minus pion minus
    likemmrsnfinder->SetCutID(0,iCutK);
    likemmrsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
    likemmrsnfinder->SetCharge(0,'-');
    likemmrsnfinder->SetCutID(1,iCutPi);
    likemmrsnfinder->SetDaughter(1,AliRsnDaughter::kPion);
    likemmrsnfinder->SetCharge(1,'-');
    likemmrsnfinder->SetResonanceMass(0.895810);
    likemmrsnfinder->SetResonancePDG(313);
    likemmrsnfinder->SetPairCuts(cutsKstar0);
    Int_t iCutlikemmKstar0=task->AddResonanceFinder(likemmrsnfinder);
    
    ////////////////////////////////
    
    //side bands
    AliRsnCutMiniPair* cutMassSB=new AliRsnCutMiniPair("cutMassSB",AliRsnCutMiniPair::kMassRange);
    //cutMassSB->SetRangeD(0.95,1.05);
    cutMassSB->SetRangeD(0.966,1.014);
    AliRsnCutSet* cutsSB=new AliRsnCutSet("pairCutsSB",AliRsnTarget::kMother);
    cutsSB->AddCut(cutMassSB);
    cutsSB->AddCut(cutYKstar0);
    cutsSB->SetCutScheme(TString::Format("%s&%s",cutMassSB->GetName(),cutYKstar0->GetName()).Data());
    
    AliRsnMiniResonanceFinder* rsnfinderSB=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderSB",task->GetName()));
    rsnfinderSB->SetCutID(0,iCutK);
    rsnfinderSB->SetDaughter(0,AliRsnDaughter::kKaon);
    rsnfinderSB->SetCharge(0,'+');
    rsnfinderSB->SetCutID(1,iCutPi);
    rsnfinderSB->SetDaughter(1,AliRsnDaughter::kPion);
    rsnfinderSB->SetCharge(1,'-');
    rsnfinder->SetResonanceMass(0.895810);
    rsnfinderSB->SetResonancePDG(313);
    rsnfinderSB->SetPairCuts(cutsSB);
    Int_t iCutSB=task->AddResonanceFinder(rsnfinderSB);
    
    AliRsnMiniResonanceFinder* antirsnfinderSB=new AliRsnMiniResonanceFinder(Form("%s_AntiResonanceFinderSB",task->GetName()));
    antirsnfinderSB->SetCutID(0,iCutK);
    antirsnfinderSB->SetDaughter(0,AliRsnDaughter::kKaon);
    antirsnfinderSB->SetCharge(0,'-');
    antirsnfinderSB->SetCutID(1,iCutPi);
    antirsnfinderSB->SetDaughter(1,AliRsnDaughter::kPion);
    antirsnfinderSB->SetCharge(1,'+');
    antirsnfinder->SetResonanceMass(0.895810);
    antirsnfinderSB->SetResonancePDG(-313);
    antirsnfinderSB->SetPairCuts(cutsSB);
    Int_t iCutantiSB=task->AddResonanceFinder(antirsnfinderSB);
    
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
    /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
    /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
    /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
    /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
    /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);
    
    // -- Create all needed outputs -----------------------------------------------------------------
    // use an array for more compact writing, which are different on mixing and charges
    
    Int_t i,xID,cut2,pairID,ipdg;
    TString name,comp;
    Char_t charge1;
    AliRsnMiniOutput* out;
    
    for(i=0;i<20;i++){
        if(!i){
            xID=imID;
            name.Form("PipKstar0");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==1){
            xID=imID;
            name.Form("PimKstar0");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutKstar0;
            pairID=0;
            ipdg=-3124;
        }else if(i==2){
            xID=imID;
            name.Form("PipSB");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutSB;
            pairID=0;
            ipdg=3124;
        }else if(i==3){
            xID=imID;
            name.Form("PimSB");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutSB;
            pairID=0;
            ipdg=-3124;
        }else if(i==4){
            xID=imID;
            name.Form("PipKstar0Mix");
            comp.Form("MIX");
            charge1='+';
            cut2=iCutKstar0;
            pairID=1;
            ipdg=3124;
        }else if(i==5){
            xID=imID;
            name.Form("PimKstar0Mix");
            comp.Form("MIX");
            charge1='-';
            cut2=iCutKstar0;
            pairID=1;
            ipdg=-3124;
        }else if(i==6){
            xID=imID;
            name.Form("PipantiKstar0");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutantiKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==7){
            xID=imID;
            name.Form("PimantiKstar0");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutantiKstar0;
            pairID=0;
            ipdg=-3124;
        }else if(i==8){
            xID=imID;
            name.Form("PipantiSB");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutantiSB;
            pairID=0;
            ipdg=3124;
        }else if(i==9){
            xID=imID;
            name.Form("PimantiSB");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutantiSB;
            pairID=0;
            ipdg=-3124;
        }else if(i==10){
            xID=imID;
            name.Form("PipantiKstar0Mix");
            comp.Form("MIX");
            charge1='+';
            cut2=iCutantiKstar0;
            pairID=1;
            ipdg=3124;
        }else if(i==11){
            xID=imID;
            name.Form("PimantiKstar0Mix");
            comp.Form("MIX");
            charge1='-';
            cut2=iCutantiKstar0;
            pairID=1;
            ipdg=3124;
        }else if(i==12){
            xID=imID;
            name.Form("PiplikeppKstar0");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutlikeppKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==13){
            xID=imID;
            name.Form("PimlikeppKstar0");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutlikeppKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==14){
            xID=imID;
            name.Form("PiplikemmKstar0");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutlikemmKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==15){
            xID=imID;
            name.Form("PimlikemmKstar0");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutlikemmKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==16){
            xID=imID;
            name.Form("PipKstar0Rotated");
            comp.Form("ROTATE1");
            charge1='+';
            cut2=iCutKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==17){
            xID=imID;
            name.Form("PimKstar0Rotated");
            comp.Form("ROTATE1");
            charge1='-';
            cut2=iCutKstar0;
            pairID=0;
            ipdg=-3124;
        }else if(i==18){
            xID=imID;
            name.Form("PipantiKstar0Rotated");
            comp.Form("ROTATE1");
            charge1='+';
            cut2=iCutantiKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==19){
            xID=imID;
            name.Form("PimantiKstar0Rotated");
            comp.Form("ROTATE1");
            charge1='-';
            cut2=iCutantiKstar0;
            pairID=0;
            ipdg=-3124;
        }
        
        out=task->CreateOutput(Form("pikstar0_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kPion);
        out->SetCutID(0,iCutPi);
        out->SetCharge(0,charge1);
        
        out->SetDaughter(1,AliRsnDaughter::kKstar0);
        out->SetCutID(1,cut2);
        out->SetCharge(1,'0');
        if(cut2!=iCutSB && cut2!=iCutantiSB) out->SetUseStoredMass(1);
        
        if(!pairID) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(ipdg);
        out->SetMotherMass(mass);
        if(xID==imID) out->AddAxis(imID,220,0.9,2.0);
        else out->AddAxis(resID,200,-0.02,0.02);
        out->AddAxis(ptID,200,0.0,20.0);// axis Y: transverse momentum
        out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
    }
    
    // fill monitoring histogram for the resonance (K*0)
    for(i=0;i<8;i++){
        if(!i){
            name.Form("kstar0mass");
            comp.Form("PAIR");
            cut2=iCutKstar0;
        }else if(i==1){
            name.Form("SBmass");
            comp.Form("PAIR");
            cut2=iCutSB;
        }else if(i==2){
            if(!isMC) continue;
            name.Form("kstar0mass_gen");
            comp.Form("MOTHER");
            cut2=iCutKstar0;
        }else if(i==3){
            if(!isMC) continue;
            name.Form("kstar0mass_rec");
            comp.Form("TRUE");
            cut2=iCutKstar0;
        }else if(i==4){
            name.Form("antikstar0mass");
            comp.Form("PAIR");
            cut2=iCutantiKstar0;
        }else if(i==5){
            name.Form("antiSBmass");
            comp.Form("PAIR");
            cut2=iCutantiSB;
        }else if(i==6){
            name.Form("likeppkstar0mass");
            comp.Form("PAIR");
            cut2=iCutlikeppKstar0;
        }else if(i==7){
            name.Form("likemmkstar0mass");
            comp.Form("PAIR");
            cut2=iCutlikemmKstar0;
        }
        
        out=task->CreateOutput(Form("pikstar0_%s",name.Data()),"HIST",comp.Data());
        out->SetMotherPDG(313);//fix this
        
        if(cut2==iCutKstar0){
            for(j=0;j<2;j++){
                out->SetDaughter(j,rsnfinder->GetDaughter(j));
                out->SetCutID(j,rsnfinder->GetCutID(j));
                out->SetCharge(j,rsnfinder->GetCharge(j));
            }
            out->SetMotherMass(rsnfinder->GetResonanceMass());
            out->SetPairCuts(cutsKstar0);
        }else if(cut2==iCutSB){
            for(j=0;j<2;j++){
                out->SetDaughter(j,rsnfinderSB->GetDaughter(j));
                out->SetCutID(j,rsnfinderSB->GetCutID(j));
                out->SetCharge(j,rsnfinderSB->GetCharge(j));
            }
            out->SetMotherMass(rsnfinderSB->GetResonanceMass());
            out->SetPairCuts(cutsSB);
        }
        //for anti
        
        else if(cut2==iCutantiKstar0){
            for(j=0;j<2;j++){
                out->SetDaughter(j,antirsnfinder->GetDaughter(j));
                out->SetCutID(j,antirsnfinder->GetCutID(j));
                out->SetCharge(j,antirsnfinder->GetCharge(j));
            }
            out->SetMotherMass(antirsnfinder->GetResonanceMass());
            out->SetPairCuts(cutsKstar0);
        }else if(cut2==iCutantiSB){
            for(j=0;j<2;j++){
                out->SetDaughter(j,antirsnfinderSB->GetDaughter(j));
                out->SetCutID(j,antirsnfinderSB->GetCutID(j));
                out->SetCharge(j,antirsnfinderSB->GetCharge(j));
            }
            out->SetMotherMass(antirsnfinderSB->GetResonanceMass());
            out->SetPairCuts(cutsSB);
        }
        
        else if(cut2==iCutlikeppKstar0){
            for(j=0;j<2;j++){
                out->SetDaughter(j,likepprsnfinder->GetDaughter(j));
                out->SetCutID(j,likepprsnfinder->GetCutID(j));
                out->SetCharge(j,likepprsnfinder->GetCharge(j));
            }
            out->SetMotherMass(likepprsnfinder->GetResonanceMass());
            out->SetPairCuts(cutsKstar0);
        }
        
        else if(cut2==iCutlikemmKstar0){
            for(j=0;j<2;j++){
                out->SetDaughter(j,likemmrsnfinder->GetDaughter(j));
                out->SetCutID(j,likemmrsnfinder->GetCutID(j));
                out->SetCharge(j,likemmrsnfinder->GetCharge(j));
            }
            out->SetMotherMass(likemmrsnfinder->GetResonanceMass());
            out->SetPairCuts(cutsKstar0);
        }
        
        //out->AddAxis(imID,200,0.85,1.05);//including all
        out->AddAxis(imID,170,0.75,1.09);
        out->AddAxis(ptID,200,0.0,20.0);
        out->AddAxis(centID,nmult,multbins);
    }
    
    return kTRUE;
}


//=============================


Bool_t Config_kxkstar0(
                       AliRsnMiniAnalysisTask *task,
                       TString     lname,
                       Bool_t      isMC,
                       Int_t       system,
                       Int_t       EventCuts,
                       Int_t       TrackCutsPi,
                       Int_t       TrackCutsKstar0 //kstar0
){
    bool isPP=false;
    if(!system) isPP=true;
    int trigger=EventCuts%10;
    int MultBins=(EventCuts/10)%10;
    if(system==1 || system==2) MultBins=1;
    
    char suffix[1000];
    sprintf(suffix,"_%s",lname.Data());
    Bool_t enableMonitor=kTRUE;
    
    Double_t mass=0.493677+0.895810;
    
    // set daughter cuts :
    //for all pions
    if(!(TrackCutsPi%10000)) TrackCutsPi+=3020;//default settings
    Float_t nsigmaPiTPC=0.1*(TrackCutsPi%100);
    Float_t nsigmaPiTOF=0.1*((TrackCutsPi/100)%100);
    Int_t CutTypePi=(TrackCutsPi/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only
    
    //for all Kaons
    if(!(TrackCutsKstar0%10000)) TrackCutsKstar0+=3020;//default settings
    Float_t nsigmaKTPC=0.1*(TrackCutsKstar0%100);
    Float_t nsigmaKTOF=0.1*((TrackCutsKstar0/100)%100);
    Int_t CutTypeK=(TrackCutsKstar0/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only
    ////////////////////////////////////
    
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
    if(!cutSetPi){cerr<<"Error in AddTaskResonanceFinder::Config_kxKstar0(): missing cutSetPi"<<endl; return kFALSE;}
    
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
    if(!cutSetK){cerr<<"Error in AddTaskResonanceFinder::Config_kxKstar0(): missing cutSetK"<<endl; return kFALSE;}
    
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
    
    // AliRsnMiniResonanceFinder//kstar0
    AliRsnCutMiniPair* cutMassKstar0=new AliRsnCutMiniPair("cutMassKstar0",AliRsnCutMiniPair::kMassRange);
    //cutMassKstar0->SetRangeD(0.85,0.95);
    cutMassKstar0->SetRangeD(0.841,0.943);
    AliRsnCutMiniPair* cutYKstar0=new AliRsnCutMiniPair("cutRapidityKstar0",AliRsnCutMiniPair::kRapidityRange);
    cutYKstar0->SetRangeD(-0.6,0.6);
    AliRsnCutSet* cutsKstar0=new AliRsnCutSet("pairCutsKstar0",AliRsnTarget::kMother);
    cutsKstar0->AddCut(cutMassKstar0);
    cutsKstar0->AddCut(cutYKstar0);
    cutsKstar0->SetCutScheme(TString::Format("%s&%s",cutMassKstar0->GetName(),cutYKstar0->GetName()).Data());
    
    //kstar0
    AliRsnMiniResonanceFinder* rsnfinder=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder",task->GetName()));
    rsnfinder->SetCutID(0,iCutK);
    rsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
    rsnfinder->SetCharge(0,'+');
    rsnfinder->SetCutID(1,iCutPi);
    rsnfinder->SetDaughter(1,AliRsnDaughter::kPion);
    rsnfinder->SetCharge(1,'-');
    rsnfinder->SetResonanceMass(0.895810);//kstar0
    rsnfinder->SetResonancePDG(313);
    rsnfinder->SetPairCuts(cutsKstar0);
    Int_t iCutKstar0=task->AddResonanceFinder(rsnfinder);
    
    //anti-kstar0
    AliRsnMiniResonanceFinder* antirsnfinder=new AliRsnMiniResonanceFinder(Form("%s_AntiResonanceFinder",task->GetName()));
    antirsnfinder->SetCutID(0,iCutK);
    antirsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
    antirsnfinder->SetCharge(0,'-');
    antirsnfinder->SetCutID(1,iCutPi);
    antirsnfinder->SetDaughter(1,AliRsnDaughter::kPion);
    antirsnfinder->SetCharge(1,'+');
    antirsnfinder->SetResonanceMass(0.895810);
    antirsnfinder->SetResonancePDG(-313);
    antirsnfinder->SetPairCuts(cutsKstar0);
    Int_t iCutantiKstar0=task->AddResonanceFinder(antirsnfinder);
    
    /////////////////////////////////
    //like sign K+ pi+
    AliRsnMiniResonanceFinder* likepprsnfinder=new AliRsnMiniResonanceFinder(Form("%s_LikeppResonanceFinder",task->GetName()));//likesign k plus pion plus
    likepprsnfinder->SetCutID(0,iCutK);
    likepprsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
    likepprsnfinder->SetCharge(0,'+');
    likepprsnfinder->SetCutID(1,iCutPi);
    likepprsnfinder->SetDaughter(1,AliRsnDaughter::kPion);
    likepprsnfinder->SetCharge(1,'+');
    likepprsnfinder->SetResonanceMass(0.895810);
    likepprsnfinder->SetResonancePDG(313);
    likepprsnfinder->SetPairCuts(cutsKstar0);
    Int_t iCutlikeppKstar0=task->AddResonanceFinder(likepprsnfinder);
    
    //like sign K- pi-
    AliRsnMiniResonanceFinder* likemmrsnfinder=new AliRsnMiniResonanceFinder(Form("%s_LikemmResonanceFinder",task->GetName()));//likesign k minus pion minus
    likemmrsnfinder->SetCutID(0,iCutK);
    likemmrsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
    likemmrsnfinder->SetCharge(0,'-');
    likemmrsnfinder->SetCutID(1,iCutPi);
    likemmrsnfinder->SetDaughter(1,AliRsnDaughter::kPion);
    likemmrsnfinder->SetCharge(1,'-');
    likemmrsnfinder->SetResonanceMass(0.895810);
    likemmrsnfinder->SetResonancePDG(313);
    likemmrsnfinder->SetPairCuts(cutsKstar0);
    Int_t iCutlikemmKstar0=task->AddResonanceFinder(likemmrsnfinder);
    
    ////////////////////////////////
    
    //side bands
    AliRsnCutMiniPair* cutMassSB=new AliRsnCutMiniPair("cutMassSB",AliRsnCutMiniPair::kMassRange);
    //cutMassSB->SetRangeD(0.95,1.05);
    cutMassSB->SetRangeD(0.966,1.014);
    AliRsnCutSet* cutsSB=new AliRsnCutSet("pairCutsSB",AliRsnTarget::kMother);
    cutsSB->AddCut(cutMassSB);
    cutsSB->AddCut(cutYKstar0);
    cutsSB->SetCutScheme(TString::Format("%s&%s",cutMassSB->GetName(),cutYKstar0->GetName()).Data());
    
    AliRsnMiniResonanceFinder* rsnfinderSB=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderSB",task->GetName()));
    rsnfinderSB->SetCutID(0,iCutK);
    rsnfinderSB->SetDaughter(0,AliRsnDaughter::kKaon);
    rsnfinderSB->SetCharge(0,'+');
    rsnfinderSB->SetCutID(1,iCutPi);
    rsnfinderSB->SetDaughter(1,AliRsnDaughter::kPion);
    rsnfinderSB->SetCharge(1,'-');
    rsnfinder->SetResonanceMass(0.895810);
    rsnfinderSB->SetResonancePDG(313);
    rsnfinderSB->SetPairCuts(cutsSB);
    Int_t iCutSB=task->AddResonanceFinder(rsnfinderSB);
    
    AliRsnMiniResonanceFinder* antirsnfinderSB=new AliRsnMiniResonanceFinder(Form("%s_AntiResonanceFinderSB",task->GetName()));
    antirsnfinderSB->SetCutID(0,iCutK);
    antirsnfinderSB->SetDaughter(0,AliRsnDaughter::kKaon);
    antirsnfinderSB->SetCharge(0,'-');
    antirsnfinderSB->SetCutID(1,iCutPi);
    antirsnfinderSB->SetDaughter(1,AliRsnDaughter::kPion);
    antirsnfinderSB->SetCharge(1,'+');
    antirsnfinder->SetResonanceMass(0.895810);
    antirsnfinderSB->SetResonancePDG(-313);
    antirsnfinderSB->SetPairCuts(cutsSB);
    Int_t iCutantiSB=task->AddResonanceFinder(antirsnfinderSB);
    
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
    /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
    /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
    /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
    /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
    /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);
    
    // -- Create all needed outputs -----------------------------------------------------------------
    // use an array for more compact writing, which are different on mixing and charges
    
    Int_t i,xID,cut2,pairID,ipdg;
    TString name,comp;
    Char_t charge1;
    AliRsnMiniOutput* out;
    
    for(i=0;i<20;i++){
        if(!i){
            xID=imID;
            name.Form("KpKstar0");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==1){
            xID=imID;
            name.Form("KmKstar0");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutKstar0;
            pairID=0;
            ipdg=-3124;
        }else if(i==2){
            xID=imID;
            name.Form("KpSB");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutSB;
            pairID=0;
            ipdg=3124;
        }else if(i==3){
            xID=imID;
            name.Form("KmSB");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutSB;
            pairID=0;
            ipdg=-3124;
        }else if(i==4){
            xID=imID;
            name.Form("KpKstar0Mix");
            comp.Form("MIX");
            charge1='+';
            cut2=iCutKstar0;
            pairID=1;
            ipdg=3124;
        }else if(i==5){
            xID=imID;
            name.Form("KmKstar0Mix");
            comp.Form("MIX");
            charge1='-';
            cut2=iCutKstar0;
            pairID=1;
            ipdg=-3124;
        }else if(i==6){
            xID=imID;
            name.Form("KpantiKstar0");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutantiKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==7){
            xID=imID;
            name.Form("KmantiKstar0");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutantiKstar0;
            pairID=0;
            ipdg=-3124;
        }else if(i==8){
            xID=imID;
            name.Form("KpantiSB");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutantiSB;
            pairID=0;
            ipdg=3124;
        }else if(i==9){
            xID=imID;
            name.Form("KmantiSB");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutantiSB;
            pairID=0;
            ipdg=-3124;
        }else if(i==10){
            xID=imID;
            name.Form("KpantiKstar0Mix");
            comp.Form("MIX");
            charge1='+';
            cut2=iCutantiKstar0;
            pairID=1;
            ipdg=3124;
        }else if(i==11){
            xID=imID;
            name.Form("KmantiKstar0Mix");
            comp.Form("MIX");
            charge1='-';
            cut2=iCutantiKstar0;
            pairID=1;
            ipdg=3124;
        }else if(i==12){
            xID=imID;
            name.Form("KplikeppKstar0");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutlikeppKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==13){
            xID=imID;
            name.Form("KmlikeppKstar0");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutlikeppKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==14){
            xID=imID;
            name.Form("KplikemmKstar0");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutlikemmKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==15){
            xID=imID;
            name.Form("KmlikemmKstar0");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutlikemmKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==16){
            xID=imID;
            name.Form("KpKstar0Rotated");
            comp.Form("ROTATE1");
            charge1='+';
            cut2=iCutKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==17){
            xID=imID;
            name.Form("KmKstar0Rotated");
            comp.Form("ROTATE1");
            charge1='-';
            cut2=iCutKstar0;
            pairID=0;
            ipdg=-3124;
        }else if(i==18){
            xID=imID;
            name.Form("KpantiKstar0Rotated");
            comp.Form("ROTATE1");
            charge1='+';
            cut2=iCutantiKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==19){
            xID=imID;
            name.Form("KmantiKstar0Rotated");
            comp.Form("ROTATE1");
            charge1='-';
            cut2=iCutantiKstar0;
            pairID=0;
            ipdg=-3124;
        }
        
        //modification in here
        out=task->CreateOutput(Form("kxkstar0_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kKaon);
        out->SetCutID(0,iCutK);
        out->SetCharge(0,charge1);
        
        out->SetDaughter(1,AliRsnDaughter::kKstar0);
        out->SetCutID(1,cut2);
        out->SetCharge(1,'0');
        if(cut2!=iCutSB && cut2!=iCutantiSB) out->SetUseStoredMass(1);
        
        if(!pairID) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(ipdg);
        out->SetMotherMass(mass);
        if(xID==imID) out->AddAxis(imID,250,1.25,2.5);
        else out->AddAxis(resID,200,-0.02,0.02);
        out->AddAxis(ptID,200,0.0,20.0);// axis Y: transverse momentum
        out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
    }
    
    // fill monitoring histogram for the resonance (K*0)
    for(i=0;i<8;i++){
        if(!i){
            name.Form("kstar0mass");
            comp.Form("PAIR");
            cut2=iCutKstar0;
        }else if(i==1){
            name.Form("SBmass");
            comp.Form("PAIR");
            cut2=iCutSB;
        }else if(i==2){
            if(!isMC) continue;
            name.Form("kstar0mass_gen");
            comp.Form("MOTHER");
            cut2=iCutKstar0;
        }else if(i==3){
            if(!isMC) continue;
            name.Form("kstar0mass_rec");
            comp.Form("TRUE");
            cut2=iCutKstar0;
        }else if(i==4){
            name.Form("antikstar0mass");
            comp.Form("PAIR");
            cut2=iCutantiKstar0;
        }else if(i==5){
            name.Form("antiSBmass");
            comp.Form("PAIR");
            cut2=iCutantiSB;
        }else if(i==6){
            name.Form("likeppkstar0mass");
            comp.Form("PAIR");
            cut2=iCutlikeppKstar0;
        }else if(i==7){
            name.Form("likemmkstar0mass");
            comp.Form("PAIR");
            cut2=iCutlikemmKstar0;
        }
        //modified
        out=task->CreateOutput(Form("kxkstar0_%s",name.Data()),"HIST",comp.Data());
        out->SetMotherPDG(313);//
        
        if(cut2==iCutKstar0){
            for(j=0;j<2;j++){
                out->SetDaughter(j,rsnfinder->GetDaughter(j));
                out->SetCutID(j,rsnfinder->GetCutID(j));
                out->SetCharge(j,rsnfinder->GetCharge(j));
            }
            out->SetMotherMass(rsnfinder->GetResonanceMass());
            out->SetPairCuts(cutsKstar0);
        }else if(cut2==iCutSB){
            for(j=0;j<2;j++){
                out->SetDaughter(j,rsnfinderSB->GetDaughter(j));
                out->SetCutID(j,rsnfinderSB->GetCutID(j));
                out->SetCharge(j,rsnfinderSB->GetCharge(j));
            }
            out->SetMotherMass(rsnfinderSB->GetResonanceMass());
            out->SetPairCuts(cutsSB);
        }
        //for anti
        
        else if(cut2==iCutantiKstar0){
            for(j=0;j<2;j++){
                out->SetDaughter(j,antirsnfinder->GetDaughter(j));
                out->SetCutID(j,antirsnfinder->GetCutID(j));
                out->SetCharge(j,antirsnfinder->GetCharge(j));
            }
            out->SetMotherMass(antirsnfinder->GetResonanceMass());
            out->SetPairCuts(cutsKstar0);
        }else if(cut2==iCutantiSB){
            for(j=0;j<2;j++){
                out->SetDaughter(j,antirsnfinderSB->GetDaughter(j));
                out->SetCutID(j,antirsnfinderSB->GetCutID(j));
                out->SetCharge(j,antirsnfinderSB->GetCharge(j));
            }
            out->SetMotherMass(antirsnfinderSB->GetResonanceMass());
            out->SetPairCuts(cutsSB);
        }
        
        else if(cut2==iCutlikeppKstar0){
            for(j=0;j<2;j++){
                out->SetDaughter(j,likepprsnfinder->GetDaughter(j));
                out->SetCutID(j,likepprsnfinder->GetCutID(j));
                out->SetCharge(j,likepprsnfinder->GetCharge(j));
            }
            out->SetMotherMass(likepprsnfinder->GetResonanceMass());
            out->SetPairCuts(cutsKstar0);
        }
        
        else if(cut2==iCutlikemmKstar0){
            for(j=0;j<2;j++){
                out->SetDaughter(j,likemmrsnfinder->GetDaughter(j));
                out->SetCutID(j,likemmrsnfinder->GetCutID(j));
                out->SetCharge(j,likemmrsnfinder->GetCharge(j));
            }
            out->SetMotherMass(likemmrsnfinder->GetResonanceMass());
            out->SetPairCuts(cutsKstar0);
        }
        
        //out->AddAxis(imID,200,0.85,1.05);//including all
        out->AddAxis(imID,170,0.75,1.09);
        out->AddAxis(ptID,200,0.0,20.0);
        out->AddAxis(centID,nmult,multbins);
    }
    
    return kTRUE;
}


//=============================


Bool_t Config_k0kstar0(
                       AliRsnMiniAnalysisTask *task,
                       TString     lname,
                       Bool_t      isMC,
                       Int_t       system,
                       Int_t       EventCuts,
                       Int_t       TrackCutsK0,
                       Int_t       TrackCutsKstar0 //kstar0
){
    bool isPP=false;
    if(!system) isPP=true;
    int trigger=EventCuts%10;
    int MultBins=(EventCuts/10)%10;
    if(system==1 || system==2) MultBins=1;
    
    char suffix[1000];
    sprintf(suffix,"_%s",lname.Data());
    Bool_t enableMonitor=kTRUE;
    
    Double_t mass=0.497611+0.895810;//895.81 from PDG 2012 //modified
    
    // set daughter cuts :
    //for all pions
    Int_t TrackCutsPi=0;
    if(!(TrackCutsPi%10000)) TrackCutsPi+=3020;//default settings
    Float_t nsigmaPiTPC=0.1*(TrackCutsPi%100);
    Float_t nsigmaPiTOF=0.1*((TrackCutsPi/100)%100);
    Int_t CutTypePi=(TrackCutsPi/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only
    
    //for all Kaons
    if(!(TrackCutsKstar0%10000)) TrackCutsKstar0+=3020;//default settings
    Float_t nsigmaKTPC=0.1*(TrackCutsKstar0%100);
    Float_t nsigmaKTOF=0.1*((TrackCutsKstar0/100)%100);
    Int_t CutTypeK=(TrackCutsKstar0/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only
    ////////////////////////////////////
    
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
    if(!cutSetPi){cerr<<"Error in AddTaskResonanceFinder::Config_k0Kstar0(): missing cutSetPi"<<endl; return kFALSE;}
    
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
    Float_t k0s_massTolSigma=4;
    Int_t   k0s_massTolID=0;
    Float_t k0s_massTol=0.03;
    Float_t k0s_massTolVeto=0.004;
    Bool_t  k0sSwitch=kFALSE;
    Float_t k0sCosPoinAn=0.97;
    
    if(TrackCutsK0==1) k0s_massTolID=1;//use pT-dependent mass tolerance cut
    
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
    if(!cutSetK){cerr<<"Error in AddTaskResonanceFinder::Config_pikstar0(): missing cutSetK"<<endl; return kFALSE;}
    
    Int_t iCutQ=task->AddTrackCuts(cutSetQ);
    Int_t iCutPi=task->AddTrackCuts(cutSetPi);
    // Int_t iCutP=task->AddTrackCuts(cutSetP);
    Int_t iCutK=task->AddTrackCuts(cutSetK);
    
    // monitoring
    TString pname="k0";
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());
        // AddMonitorOutput(isMC,cutSetP->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetK->GetMonitorOutput());
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetK0s->GetMonitorOutput(), "K0S");
    }
    
    // AliRsnMiniResonanceFinder//kstar0
    AliRsnCutMiniPair* cutMassKstar0=new AliRsnCutMiniPair("cutMassKstar0",AliRsnCutMiniPair::kMassRange);
    //cutMassKstar0->SetRangeD(0.85,0.95);
    cutMassKstar0->SetRangeD(0.841,0.943);
    AliRsnCutMiniPair* cutYKstar0=new AliRsnCutMiniPair("cutRapidityKstar0",AliRsnCutMiniPair::kRapidityRange);
    cutYKstar0->SetRangeD(-0.6,0.6);
    AliRsnCutSet* cutsKstar0=new AliRsnCutSet("pairCutsKstar0",AliRsnTarget::kMother);
    cutsKstar0->AddCut(cutMassKstar0);
    cutsKstar0->AddCut(cutYKstar0);
    cutsKstar0->SetCutScheme(TString::Format("%s&%s",cutMassKstar0->GetName(),cutYKstar0->GetName()).Data());
    
    //kstar0
    AliRsnMiniResonanceFinder* rsnfinder=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder",task->GetName()));
    rsnfinder->SetCutID(0,iCutK);
    rsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
    rsnfinder->SetCharge(0,'+');
    rsnfinder->SetCutID(1,iCutPi);
    rsnfinder->SetDaughter(1,AliRsnDaughter::kPion);
    rsnfinder->SetCharge(1,'-');
    rsnfinder->SetResonanceMass(0.895810);//kstar0
    rsnfinder->SetResonancePDG(313);
    rsnfinder->SetPairCuts(cutsKstar0);
    Int_t iCutKstar0=task->AddResonanceFinder(rsnfinder);
    
    //anti-kstar0
    AliRsnMiniResonanceFinder* antirsnfinder=new AliRsnMiniResonanceFinder(Form("%s_AntiResonanceFinder",task->GetName()));
    antirsnfinder->SetCutID(0,iCutK);
    antirsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
    antirsnfinder->SetCharge(0,'-');
    antirsnfinder->SetCutID(1,iCutPi);
    antirsnfinder->SetDaughter(1,AliRsnDaughter::kPion);
    antirsnfinder->SetCharge(1,'+');
    antirsnfinder->SetResonanceMass(0.895810);
    antirsnfinder->SetResonancePDG(-313);
    antirsnfinder->SetPairCuts(cutsKstar0);
    Int_t iCutantiKstar0=task->AddResonanceFinder(antirsnfinder);
    
    /////////////////////////////////
    //like sign K+ pi+
    AliRsnMiniResonanceFinder* likepprsnfinder=new AliRsnMiniResonanceFinder(Form("%s_LikeppResonanceFinder",task->GetName()));//likesign k plus pion plus
    likepprsnfinder->SetCutID(0,iCutK);
    likepprsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
    likepprsnfinder->SetCharge(0,'+');
    likepprsnfinder->SetCutID(1,iCutPi);
    likepprsnfinder->SetDaughter(1,AliRsnDaughter::kPion);
    likepprsnfinder->SetCharge(1,'+');
    likepprsnfinder->SetResonanceMass(0.895810);
    likepprsnfinder->SetResonancePDG(313);
    likepprsnfinder->SetPairCuts(cutsKstar0);
    Int_t iCutlikeppKstar0=task->AddResonanceFinder(likepprsnfinder);
    
    //like sign K- pi-
    AliRsnMiniResonanceFinder* likemmrsnfinder=new AliRsnMiniResonanceFinder(Form("%s_LikemmResonanceFinder",task->GetName()));//likesign k minus pion minus
    likemmrsnfinder->SetCutID(0,iCutK);
    likemmrsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
    likemmrsnfinder->SetCharge(0,'-');
    likemmrsnfinder->SetCutID(1,iCutPi);
    likemmrsnfinder->SetDaughter(1,AliRsnDaughter::kPion);
    likemmrsnfinder->SetCharge(1,'-');
    likemmrsnfinder->SetResonanceMass(0.895810);
    likemmrsnfinder->SetResonancePDG(313);
    likemmrsnfinder->SetPairCuts(cutsKstar0);
    Int_t iCutlikemmKstar0=task->AddResonanceFinder(likemmrsnfinder);
    
    ////////////////////////////////
    
    //side bands
    AliRsnCutMiniPair* cutMassSB=new AliRsnCutMiniPair("cutMassSB",AliRsnCutMiniPair::kMassRange);
    //cutMassSB->SetRangeD(0.95,1.05);
    cutMassSB->SetRangeD(0.966,1.014);
    AliRsnCutSet* cutsSB=new AliRsnCutSet("pairCutsSB",AliRsnTarget::kMother);
    cutsSB->AddCut(cutMassSB);
    cutsSB->AddCut(cutYKstar0);
    cutsSB->SetCutScheme(TString::Format("%s&%s",cutMassSB->GetName(),cutYKstar0->GetName()).Data());
    
    AliRsnMiniResonanceFinder* rsnfinderSB=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderSB",task->GetName()));
    rsnfinderSB->SetCutID(0,iCutK);
    rsnfinderSB->SetDaughter(0,AliRsnDaughter::kKaon);
    rsnfinderSB->SetCharge(0,'+');
    rsnfinderSB->SetCutID(1,iCutPi);
    rsnfinderSB->SetDaughter(1,AliRsnDaughter::kPion);
    rsnfinderSB->SetCharge(1,'-');
    rsnfinder->SetResonanceMass(0.895810);
    rsnfinderSB->SetResonancePDG(313);
    rsnfinderSB->SetPairCuts(cutsSB);
    Int_t iCutSB=task->AddResonanceFinder(rsnfinderSB);
    
    AliRsnMiniResonanceFinder* antirsnfinderSB=new AliRsnMiniResonanceFinder(Form("%s_AntiResonanceFinderSB",task->GetName()));
    antirsnfinderSB->SetCutID(0,iCutK);
    antirsnfinderSB->SetDaughter(0,AliRsnDaughter::kKaon);
    antirsnfinderSB->SetCharge(0,'-');
    antirsnfinderSB->SetCutID(1,iCutPi);
    antirsnfinderSB->SetDaughter(1,AliRsnDaughter::kPion);
    antirsnfinderSB->SetCharge(1,'+');
    antirsnfinder->SetResonanceMass(0.895810);
    antirsnfinderSB->SetResonancePDG(-313);
    antirsnfinderSB->SetPairCuts(cutsSB);
    Int_t iCutantiSB=task->AddResonanceFinder(antirsnfinderSB);
    
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
    /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
    /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
    /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
    /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
    /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);
    
    // -- Create all needed outputs -----------------------------------------------------------------
    // use an array for more compact writing, which are different on mixing and charges
    
    Int_t i,xID,cut2,pairID,ipdg;
    TString name,comp;
    Char_t charge1;
    AliRsnMiniOutput* out;
    
    for(i=0;i<10;i++){
        if(!i){
            xID=imID;
            name.Form("K0Kstar0");
            comp.Form("PAIR");
            charge1='0';
            cut2=iCutKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==1){
            xID=imID;
            name.Form("K0SB");
            comp.Form("PAIR");
            charge1='0';
            cut2=iCutSB;
            pairID=0;
            ipdg=3124;
        }else if(i==2){
            xID=imID;
            name.Form("K0Kstar0Mix");
            comp.Form("MIX");
            charge1='0';
            cut2=iCutKstar0;
            pairID=1;
            ipdg=3124;
        }else if(i==3){
            xID=imID;
            name.Form("K0antiKstar0");
            comp.Form("PAIR");
            charge1='0';
            cut2=iCutantiKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==4){
            xID=imID;
            name.Form("K0antiSB");
            comp.Form("PAIR");
            charge1='0';
            cut2=iCutantiSB;
            pairID=0;
            ipdg=3124;
        }else if(i==5){
            xID=imID;
            name.Form("K0antiKstar0Mix");
            comp.Form("MIX");
            charge1='0';
            cut2=iCutantiKstar0;
            pairID=1;
            ipdg=3124;
        }else if(i==6){
            xID=imID;
            name.Form("K0likeppKstar0");
            comp.Form("PAIR");
            charge1='0';
            cut2=iCutlikeppKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==7){
            xID=imID;
            name.Form("K0likemmKstar0");
            comp.Form("PAIR");
            charge1='0';
            cut2=iCutlikemmKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==8){
            xID=imID;
            name.Form("K0Kstar0Rotated");
            comp.Form("ROTATE1");
            charge1='0';
            cut2=iCutKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==9){
            xID=imID;
            name.Form("K0antiKstar0Rotated");
            comp.Form("ROTATE1");
            charge1='0';
            cut2=iCutantiKstar0;
            pairID=0;
            ipdg=3124;
        }
        
        //modification in here
        out=task->CreateOutput(Form("k0kstar0_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kKaon0);
        out->SetCutID(0,iCutK0s);
        out->SetCharge(0,charge1);
        
        out->SetDaughter(1,AliRsnDaughter::kKstar0);
        out->SetCutID(1,cut2);
        out->SetCharge(1,'0');
        if(cut2!=iCutSB && cut2!=iCutantiSB) out->SetUseStoredMass(1);
        
        if(!pairID) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(ipdg);
        out->SetMotherMass(mass);
        if(xID==imID) out->AddAxis(imID,250,1.25,2.5);// if(xID==imID) out->AddAxis(imID,246,1.77,3.0); //    if(xID==imID) out->AddAxis(imID,200,1.4,2.4);
        else out->AddAxis(resID,200,-0.02,0.02);
        out->AddAxis(ptID,200,0.0,20.0);// axis Y: transverse momentum
        out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
    }
    
    // fill monitoring histogram for the resonance (K*0)
    for(i=0;i<8;i++){
        if(!i){
            name.Form("kstar0mass");
            comp.Form("PAIR");
            cut2=iCutKstar0;
        }else if(i==1){
            name.Form("SBmass");
            comp.Form("PAIR");
            cut2=iCutSB;
        }else if(i==2){
            if(!isMC) continue;
            name.Form("kstar0mass_gen");
            comp.Form("MOTHER");
            cut2=iCutKstar0;
        }else if(i==3){
            if(!isMC) continue;
            name.Form("kstar0mass_rec");
            comp.Form("TRUE");
            cut2=iCutKstar0;
        }else if(i==4){
            name.Form("antikstar0mass");
            comp.Form("PAIR");
            cut2=iCutantiKstar0;
        }else if(i==5){
            name.Form("antiSBmass");
            comp.Form("PAIR");
            cut2=iCutantiSB;
        }else if(i==6){
            name.Form("likeppkstar0mass");
            comp.Form("PAIR");
            cut2=iCutlikeppKstar0;
        }else if(i==7){
            name.Form("likemmkstar0mass");
            comp.Form("PAIR");
            cut2=iCutlikemmKstar0;
        }
        //modified
        out=task->CreateOutput(Form("pkstar0_%s",name.Data()),"HIST",comp.Data());
        out->SetMotherPDG(313);//
        
        if(cut2==iCutKstar0){
            for(j=0;j<2;j++){
                out->SetDaughter(j,rsnfinder->GetDaughter(j));
                out->SetCutID(j,rsnfinder->GetCutID(j));
                out->SetCharge(j,rsnfinder->GetCharge(j));
            }
            out->SetMotherMass(rsnfinder->GetResonanceMass());
            out->SetPairCuts(cutsKstar0);
        }else if(cut2==iCutSB){
            for(j=0;j<2;j++){
                out->SetDaughter(j,rsnfinderSB->GetDaughter(j));
                out->SetCutID(j,rsnfinderSB->GetCutID(j));
                out->SetCharge(j,rsnfinderSB->GetCharge(j));
            }
            out->SetMotherMass(rsnfinderSB->GetResonanceMass());
            out->SetPairCuts(cutsSB);
        }
        //for anti
        
        else if(cut2==iCutantiKstar0){
            for(j=0;j<2;j++){
                out->SetDaughter(j,antirsnfinder->GetDaughter(j));
                out->SetCutID(j,antirsnfinder->GetCutID(j));
                out->SetCharge(j,antirsnfinder->GetCharge(j));
            }
            out->SetMotherMass(antirsnfinder->GetResonanceMass());
            out->SetPairCuts(cutsKstar0);
        }else if(cut2==iCutantiSB){
            for(j=0;j<2;j++){
                out->SetDaughter(j,antirsnfinderSB->GetDaughter(j));
                out->SetCutID(j,antirsnfinderSB->GetCutID(j));
                out->SetCharge(j,antirsnfinderSB->GetCharge(j));
            }
            out->SetMotherMass(antirsnfinderSB->GetResonanceMass());
            out->SetPairCuts(cutsSB);
        }
        
        else if(cut2==iCutlikeppKstar0){
            for(j=0;j<2;j++){
                out->SetDaughter(j,likepprsnfinder->GetDaughter(j));
                out->SetCutID(j,likepprsnfinder->GetCutID(j));
                out->SetCharge(j,likepprsnfinder->GetCharge(j));
            }
            out->SetMotherMass(likepprsnfinder->GetResonanceMass());
            out->SetPairCuts(cutsKstar0);
        }
        
        else if(cut2==iCutlikemmKstar0){
            for(j=0;j<2;j++){
                out->SetDaughter(j,likemmrsnfinder->GetDaughter(j));
                out->SetCutID(j,likemmrsnfinder->GetCutID(j));
                out->SetCharge(j,likemmrsnfinder->GetCharge(j));
            }
            out->SetMotherMass(likemmrsnfinder->GetResonanceMass());
            out->SetPairCuts(cutsKstar0);
        }
        
        //out->AddAxis(imID,200,0.85,1.05);//including all
        out->AddAxis(imID,170,0.75,1.09);
        out->AddAxis(ptID,200,0.0,20.0);
        out->AddAxis(centID,nmult,multbins);
    }
    
    return kTRUE;
}


//=============================


Bool_t Config_pkstar0(
                      AliRsnMiniAnalysisTask *task,
                      TString     lname,
                      Bool_t      isMC,
                      Int_t       system,
                      Int_t       EventCuts,
                      Int_t       TrackCutsP,
                      Int_t       TrackCutsKstar0 //kstar0
){
    bool isPP=false;
    if(!system) isPP=true;
    int trigger=EventCuts%10;
    int MultBins=(EventCuts/10)%10;
    if(system==1 || system==2) MultBins=1;
    
    char suffix[1000];
    sprintf(suffix,"_%s",lname.Data());
    Bool_t enableMonitor=kTRUE;
    
    Double_t mass=0.938272+0.895810;//895.81 from PDG 2012 //modified
    
    // set daughter cuts :
    //for all pions
    Int_t TrackCutsPi=0;
    if(!(TrackCutsPi%10000)) TrackCutsPi+=3020;//default settings
    Float_t nsigmaPiTPC=0.1*(TrackCutsPi%100);
    Float_t nsigmaPiTOF=0.1*((TrackCutsPi/100)%100);
    Int_t CutTypePi=(TrackCutsPi/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only
    
    //proton
    if(!(TrackCutsP%10000)) TrackCutsP+=3020;//default settings
    Float_t nsigmaPTPC=0.1*(TrackCutsP%100);
    Float_t nsigmaPTOF=0.1*((TrackCutsP/100)%100);
    Int_t CutTypeP=(TrackCutsP/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only
    
    //for all Kaons
    if(!(TrackCutsKstar0%10000)) TrackCutsKstar0+=3020;//default settings
    Float_t nsigmaKTPC=0.1*(TrackCutsKstar0%100);
    Float_t nsigmaKTOF=0.1*((TrackCutsKstar0/100)%100);
    Int_t CutTypeK=(TrackCutsKstar0/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only
    ////////////////////////////////////
    
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
    if(!cutSetPi){cerr<<"Error in AddTaskResonanceFinder::Config_pKstar0(): missing cutSetPi"<<endl; return kFALSE;}
    
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
    if(!cutSetP){cerr<<"Error in AddTaskResonanceFinder::Config_pKstar0(): missing cutSetP"<<endl; return kFALSE;}
    
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
    if(!cutSetK){cerr<<"Error in AddTaskResonanceFinder::Config_pKstar0(): missing cutSetK"<<endl; return kFALSE;}
    
    Int_t iCutQ=task->AddTrackCuts(cutSetQ);
    Int_t iCutPi=task->AddTrackCuts(cutSetPi);
    Int_t iCutP=task->AddTrackCuts(cutSetP);
    Int_t iCutK=task->AddTrackCuts(cutSetK);
    
    // monitoring
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetP->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetK->GetMonitorOutput());
    }
    
    // AliRsnMiniResonanceFinder//kstar0
    AliRsnCutMiniPair* cutMassKstar0=new AliRsnCutMiniPair("cutMassKstar0",AliRsnCutMiniPair::kMassRange);
    //cutMassKstar0->SetRangeD(0.85,0.95);
    cutMassKstar0->SetRangeD(0.841,0.943);
    AliRsnCutMiniPair* cutYKstar0=new AliRsnCutMiniPair("cutRapidityKstar0",AliRsnCutMiniPair::kRapidityRange);
    cutYKstar0->SetRangeD(-0.6,0.6);
    AliRsnCutSet* cutsKstar0=new AliRsnCutSet("pairCutsKstar0",AliRsnTarget::kMother);
    cutsKstar0->AddCut(cutMassKstar0);
    cutsKstar0->AddCut(cutYKstar0);
    cutsKstar0->SetCutScheme(TString::Format("%s&%s",cutMassKstar0->GetName(),cutYKstar0->GetName()).Data());
    
    //kstar0
    AliRsnMiniResonanceFinder* rsnfinder=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder",task->GetName()));
    rsnfinder->SetCutID(0,iCutK);
    rsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
    rsnfinder->SetCharge(0,'+');
    rsnfinder->SetCutID(1,iCutPi);
    rsnfinder->SetDaughter(1,AliRsnDaughter::kPion);
    rsnfinder->SetCharge(1,'-');
    rsnfinder->SetResonanceMass(0.895810);//kstar0
    rsnfinder->SetResonancePDG(313);
    rsnfinder->SetPairCuts(cutsKstar0);
    Int_t iCutKstar0=task->AddResonanceFinder(rsnfinder);
    
    //anti-kstar0
    AliRsnMiniResonanceFinder* antirsnfinder=new AliRsnMiniResonanceFinder(Form("%s_AntiResonanceFinder",task->GetName()));
    antirsnfinder->SetCutID(0,iCutK);
    antirsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
    antirsnfinder->SetCharge(0,'-');
    antirsnfinder->SetCutID(1,iCutPi);
    antirsnfinder->SetDaughter(1,AliRsnDaughter::kPion);
    antirsnfinder->SetCharge(1,'+');
    antirsnfinder->SetResonanceMass(0.895810);
    antirsnfinder->SetResonancePDG(-313);
    antirsnfinder->SetPairCuts(cutsKstar0);
    Int_t iCutantiKstar0=task->AddResonanceFinder(antirsnfinder);
    
    /////////////////////////////////
    //like sign K+ pi+
    AliRsnMiniResonanceFinder* likepprsnfinder=new AliRsnMiniResonanceFinder(Form("%s_LikeppResonanceFinder",task->GetName()));//likesign k plus pion plus
    likepprsnfinder->SetCutID(0,iCutK);
    likepprsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
    likepprsnfinder->SetCharge(0,'+');
    likepprsnfinder->SetCutID(1,iCutPi);
    likepprsnfinder->SetDaughter(1,AliRsnDaughter::kPion);
    likepprsnfinder->SetCharge(1,'+');
    likepprsnfinder->SetResonanceMass(0.895810);
    likepprsnfinder->SetResonancePDG(313);
    likepprsnfinder->SetPairCuts(cutsKstar0);
    Int_t iCutlikeppKstar0=task->AddResonanceFinder(likepprsnfinder);
    
    //like sign K- pi-
    AliRsnMiniResonanceFinder* likemmrsnfinder=new AliRsnMiniResonanceFinder(Form("%s_LikemmResonanceFinder",task->GetName()));//likesign k minus pion minus
    likemmrsnfinder->SetCutID(0,iCutK);
    likemmrsnfinder->SetDaughter(0,AliRsnDaughter::kKaon);
    likemmrsnfinder->SetCharge(0,'-');
    likemmrsnfinder->SetCutID(1,iCutPi);
    likemmrsnfinder->SetDaughter(1,AliRsnDaughter::kPion);
    likemmrsnfinder->SetCharge(1,'-');
    likemmrsnfinder->SetResonanceMass(0.895810);
    likemmrsnfinder->SetResonancePDG(313);
    likemmrsnfinder->SetPairCuts(cutsKstar0);
    Int_t iCutlikemmKstar0=task->AddResonanceFinder(likemmrsnfinder);
    
    ////////////////////////////////
    
    //side bands
    AliRsnCutMiniPair* cutMassSB=new AliRsnCutMiniPair("cutMassSB",AliRsnCutMiniPair::kMassRange);
    //cutMassSB->SetRangeD(0.95,1.05);
    cutMassSB->SetRangeD(0.966,1.014);
    AliRsnCutSet* cutsSB=new AliRsnCutSet("pairCutsSB",AliRsnTarget::kMother);
    cutsSB->AddCut(cutMassSB);
    cutsSB->AddCut(cutYKstar0);
    cutsSB->SetCutScheme(TString::Format("%s&%s",cutMassSB->GetName(),cutYKstar0->GetName()).Data());
    
    AliRsnMiniResonanceFinder* rsnfinderSB=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderSB",task->GetName()));
    rsnfinderSB->SetCutID(0,iCutK);
    rsnfinderSB->SetDaughter(0,AliRsnDaughter::kKaon);
    rsnfinderSB->SetCharge(0,'+');
    rsnfinderSB->SetCutID(1,iCutPi);
    rsnfinderSB->SetDaughter(1,AliRsnDaughter::kPion);
    rsnfinderSB->SetCharge(1,'-');
    rsnfinder->SetResonanceMass(0.895810);
    rsnfinderSB->SetResonancePDG(313);
    rsnfinderSB->SetPairCuts(cutsSB);
    Int_t iCutSB=task->AddResonanceFinder(rsnfinderSB);
    
    AliRsnMiniResonanceFinder* antirsnfinderSB=new AliRsnMiniResonanceFinder(Form("%s_AntiResonanceFinderSB",task->GetName()));
    antirsnfinderSB->SetCutID(0,iCutK);
    antirsnfinderSB->SetDaughter(0,AliRsnDaughter::kKaon);
    antirsnfinderSB->SetCharge(0,'-');
    antirsnfinderSB->SetCutID(1,iCutPi);
    antirsnfinderSB->SetDaughter(1,AliRsnDaughter::kPion);
    antirsnfinderSB->SetCharge(1,'+');
    antirsnfinder->SetResonanceMass(0.895810);
    antirsnfinderSB->SetResonancePDG(-313);
    antirsnfinderSB->SetPairCuts(cutsSB);
    Int_t iCutantiSB=task->AddResonanceFinder(antirsnfinderSB);
    
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
    /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
    /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
    /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
    /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
    /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);
    
    // -- Create all needed outputs -----------------------------------------------------------------
    // use an array for more compact writing, which are different on mixing and charges
    
    Int_t i,xID,cut2,pairID,ipdg;
    TString name,comp;
    Char_t charge1;
    AliRsnMiniOutput* out;
    
    for(i=0;i<20;i++){
        if(!i){
            xID=imID;
            name.Form("PpKstar0");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==1){
            xID=imID;
            name.Form("PmKstar0");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutKstar0;
            pairID=0;
            ipdg=-3124;
        }else if(i==2){
            xID=imID;
            name.Form("PpSB");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutSB;
            pairID=0;
            ipdg=3124;
        }else if(i==3){
            xID=imID;
            name.Form("PmSB");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutSB;
            pairID=0;
            ipdg=-3124;
        }else if(i==4){
            xID=imID;
            name.Form("PpKstar0Mix");
            comp.Form("MIX");
            charge1='+';
            cut2=iCutKstar0;
            pairID=1;
            ipdg=3124;
        }else if(i==5){
            xID=imID;
            name.Form("PmKstar0Mix");
            comp.Form("MIX");
            charge1='-';
            cut2=iCutKstar0;
            pairID=1;
            ipdg=-3124;
        }else if(i==6){
            xID=imID;
            name.Form("PpantiKstar0");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutantiKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==7){
            xID=imID;
            name.Form("PmantiKstar0");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutantiKstar0;
            pairID=0;
            ipdg=-3124;
        }else if(i==8){
            xID=imID;
            name.Form("PpantiSB");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutantiSB;
            pairID=0;
            ipdg=3124;
        }else if(i==9){
            xID=imID;
            name.Form("PmantiSB");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutantiSB;
            pairID=0;
            ipdg=-3124;
        }else if(i==10){
            xID=imID;
            name.Form("PpantiKstar0Mix");
            comp.Form("MIX");
            charge1='+';
            cut2=iCutantiKstar0;
            pairID=1;
            ipdg=3124;
        }else if(i==11){
            xID=imID;
            name.Form("PmantiKstar0Mix");
            comp.Form("MIX");
            charge1='-';
            cut2=iCutantiKstar0;
            pairID=1;
            ipdg=3124;
        }else if(i==12){
            xID=imID;
            name.Form("PplikeppKstar0");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutlikeppKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==13){
            xID=imID;
            name.Form("PmlikeppKstar0");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutlikeppKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==14){
            xID=imID;
            name.Form("PplikemmKstar0");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutlikemmKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==15){
            xID=imID;
            name.Form("PmlikemmKstar0");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutlikemmKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==16){
            xID=imID;
            name.Form("PpKstar0Rotated");
            comp.Form("ROTATE1");
            charge1='+';
            cut2=iCutKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==17){
            xID=imID;
            name.Form("PmKstar0Rotated");
            comp.Form("ROTATE1");
            charge1='-';
            cut2=iCutKstar0;
            pairID=0;
            ipdg=-3124;
        }else if(i==18){
            xID=imID;
            name.Form("PpantiKstar0Rotated");
            comp.Form("ROTATE1");
            charge1='+';
            cut2=iCutantiKstar0;
            pairID=0;
            ipdg=3124;
        }else if(i==19){
            xID=imID;
            name.Form("PmantiKstar0Rotated");
            comp.Form("ROTATE1");
            charge1='-';
            cut2=iCutantiKstar0;
            pairID=0;
            ipdg=-3124;
        }
        
        //modification in here
        out=task->CreateOutput(Form("pkstar0_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kProton);
        out->SetCutID(0,iCutP);
        out->SetCharge(0,charge1);
        
        out->SetDaughter(1,AliRsnDaughter::kKstar0);
        out->SetCutID(1,cut2);
        out->SetCharge(1,'0');
        if(cut2!=iCutSB && cut2!=iCutantiSB) out->SetUseStoredMass(1);
        
        if(!pairID) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(ipdg);
        out->SetMotherMass(mass);
        if(xID==imID) out->AddAxis(imID,260,1.7,3.0); //    if(xID==imID) out->AddAxis(imID,200,1.4,2.4);
        else out->AddAxis(resID,200,-0.02,0.02);
        out->AddAxis(ptID,200,0.0,20.0);// axis Y: transverse momentum
        out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
    }
    
    // fill monitoring histogram for the resonance (K*0)
    for(i=0;i<8;i++){
        if(!i){
            name.Form("kstar0mass");
            comp.Form("PAIR");
            cut2=iCutKstar0;
        }else if(i==1){
            name.Form("SBmass");
            comp.Form("PAIR");
            cut2=iCutSB;
        }else if(i==2){
            if(!isMC) continue;
            name.Form("kstar0mass_gen");
            comp.Form("MOTHER");
            cut2=iCutKstar0;
        }else if(i==3){
            if(!isMC) continue;
            name.Form("kstar0mass_rec");
            comp.Form("TRUE");
            cut2=iCutKstar0;
        }else if(i==4){
            name.Form("antikstar0mass");
            comp.Form("PAIR");
            cut2=iCutantiKstar0;
        }else if(i==5){
            name.Form("antiSBmass");
            comp.Form("PAIR");
            cut2=iCutantiSB;
        }else if(i==6){
            name.Form("likeppkstar0mass");
            comp.Form("PAIR");
            cut2=iCutlikeppKstar0;
        }else if(i==7){
            name.Form("likemmkstar0mass");
            comp.Form("PAIR");
            cut2=iCutlikemmKstar0;
        }
        //modified
        out=task->CreateOutput(Form("pkstar0_%s",name.Data()),"HIST",comp.Data());
        out->SetMotherPDG(313);//
        
        if(cut2==iCutKstar0){
            for(j=0;j<2;j++){
                out->SetDaughter(j,rsnfinder->GetDaughter(j));
                out->SetCutID(j,rsnfinder->GetCutID(j));
                out->SetCharge(j,rsnfinder->GetCharge(j));
            }
            out->SetMotherMass(rsnfinder->GetResonanceMass());
            out->SetPairCuts(cutsKstar0);
        }else if(cut2==iCutSB){
            for(j=0;j<2;j++){
                out->SetDaughter(j,rsnfinderSB->GetDaughter(j));
                out->SetCutID(j,rsnfinderSB->GetCutID(j));
                out->SetCharge(j,rsnfinderSB->GetCharge(j));
            }
            out->SetMotherMass(rsnfinderSB->GetResonanceMass());
            out->SetPairCuts(cutsSB);
        }
        //for anti
        
        else if(cut2==iCutantiKstar0){
            for(j=0;j<2;j++){
                out->SetDaughter(j,antirsnfinder->GetDaughter(j));
                out->SetCutID(j,antirsnfinder->GetCutID(j));
                out->SetCharge(j,antirsnfinder->GetCharge(j));
            }
            out->SetMotherMass(antirsnfinder->GetResonanceMass());
            out->SetPairCuts(cutsKstar0);
        }else if(cut2==iCutantiSB){
            for(j=0;j<2;j++){
                out->SetDaughter(j,antirsnfinderSB->GetDaughter(j));
                out->SetCutID(j,antirsnfinderSB->GetCutID(j));
                out->SetCharge(j,antirsnfinderSB->GetCharge(j));
            }
            out->SetMotherMass(antirsnfinderSB->GetResonanceMass());
            out->SetPairCuts(cutsSB);
        }
        
        else if(cut2==iCutlikeppKstar0){
            for(j=0;j<2;j++){
                out->SetDaughter(j,likepprsnfinder->GetDaughter(j));
                out->SetCutID(j,likepprsnfinder->GetCutID(j));
                out->SetCharge(j,likepprsnfinder->GetCharge(j));
            }
            out->SetMotherMass(likepprsnfinder->GetResonanceMass());
            out->SetPairCuts(cutsKstar0);
        }
        
        else if(cut2==iCutlikemmKstar0){
            for(j=0;j<2;j++){
                out->SetDaughter(j,likemmrsnfinder->GetDaughter(j));
                out->SetCutID(j,likemmrsnfinder->GetCutID(j));
                out->SetCharge(j,likemmrsnfinder->GetCharge(j));
            }
            out->SetMotherMass(likemmrsnfinder->GetResonanceMass());
            out->SetPairCuts(cutsKstar0);
        }
        
        //out->AddAxis(imID,200,0.85,1.05);//including all
        out->AddAxis(imID,170,0.75,1.09);
        out->AddAxis(ptID,200,0.0,20.0);
        out->AddAxis(centID,nmult,multbins);
    }
    
    return kTRUE;
}


//=============================


Bool_t Config_pikstarx(
                       AliRsnMiniAnalysisTask *task,
                       TString     lname,
                       Bool_t      isMC,
                       Int_t       system,
                       Int_t       EventCuts,
                       Int_t       TrackCutsK,
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
    
    Double_t mass=0.139571+0.89166;
    
    // set cuts for pions
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
    if(!cutSetPi){cerr<<"Error in AddTaskRare_pp13::Config_pikstarx(): missing cutSetPi"<<endl; return kFALSE;}
    
    //Int_t iCutQ=task->AddTrackCuts(cutSetQ);
    Int_t iCutPi=task->AddTrackCuts(cutSetPi);
    
    // set cuts for K0S from K*
    Int_t SidebandKstar=(TrackCutsK/10000)%10;
    Int_t K0sCuts=TrackCutsK%10000;
    
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
    cutK0s->SetMaxRapidity(2.);
    cutK0s->SetMinTPCcluster(-1);
    
    AliRsnCutSet* cutSetK0s=new AliRsnCutSet("setK0s",AliRsnTarget::kDaughter);
    cutSetK0s->AddCut(cutK0s);
    cutSetK0s->SetCutScheme(cutK0s->GetName());
    Int_t iCutK0s=task->AddTrackCuts(cutSetK0s);
    
    // monitoring
    TString pname="k0s";
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        //AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetK0s->GetMonitorOutput(), "K0S");
    }
    // AliRsnMiniResonanceFinder - K*
    AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
    AliRsnCutMiniPair* cutYRes=new AliRsnCutMiniPair("cutRapidityResonance",AliRsnCutMiniPair::kRapidityRange);
    cutYRes->SetRangeD(-0.6,0.6);
    AliRsnMiniResonanceFinder* Kfinder[4];
    Int_t i,iCutKstar[4];
    
    AliRsnCutMiniPair* cutMassKstar=new AliRsnCutMiniPair("cutMassKstar",AliRsnCutMiniPair::kMassRange);
    cutMassKstar->SetRangeD(0.841,0.943);
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
    cutMassKstarSB->SetRangeD(0.966,1.014);
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
    
    // pair cuts
    AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
    if(system!=1) cutY->SetRangeD(-0.5,0.5);
    else cutY->SetRangeD(-0.465,0.035);
    
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
    
    Int_t k,cut1;
    TString name,comp;
    Char_t charge1,charge2;
    AliRsnMiniOutput* out;
    
    for(i=0;i<4;i++) for(j=0;j<2;j++) for(k=0;k<3;k++){
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
            name.Append("Pip");
            charge2='+';
            
        }else if(j==1){
            name.Append("Pim");
            charge2='-';
        }
        
        if(!k){
            comp.Form("PAIR");
        }else if(k==1){
            name.Append("Mix");
            comp.Form("MIX");
        }else if(k==2){
            name.Append("Rotate");
            comp.Form("ROTATE1");
        }
        
        if(i>=2 && k) continue;
        
        out=task->CreateOutput(Form("pikstarx_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kKstarpm);
        out->SetCutID(0,iCutKstar[cut1]);
        out->SetCharge(0,charge1);
        if(i<=1) out->SetUseStoredMass(0);
        
        out->SetDaughter(1,AliRsnDaughter::kPion);
        out->SetCutID(1,iCutPi);
        out->SetCharge(1,charge2);
        
        if(k<=1) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(3124);
        out->SetMotherMass(mass);
        
        out->AddAxis(imID,220,0.9,2.0);// axis X: invmass or resolution
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
        
        out=task->CreateOutput(Form("pikstarx_%s",name.Data()),"HIST",comp.Data());
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
    
    return kTRUE;
}


//=============================


Bool_t Config_kxkstarx(
                       AliRsnMiniAnalysisTask *task,
                       TString     lname,
                       Bool_t      isMC,
                       Int_t       system,
                       Int_t       EventCuts,
                       Int_t       TrackCutsK,
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
    
    Double_t mass=0.493677+0.89166;
    
    // set cuts for pions
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
    if(!cutSetPi){cerr<<"Error in AddTaskRare_pp13::Config_kxkstar0(): missing cutSetPi"<<endl; return kFALSE;}
    
    // set cuts for kaons
    int TrackCutsKx=0;
    if(!(TrackCutsKx%10000)) TrackCutsKx+=3020;//default settings
    Float_t nsigmaKTPC=0.1*(TrackCutsKx%100);
    Float_t nsigmaKTOF=0.1*((TrackCutsKx/100)%100);
    Int_t CutTypeK=(TrackCutsKx/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only
    
    AliRsnCutSetDaughterParticle* cutSetKx=0;
    if(!CutTypeK) cutSetKx=new AliRsnCutSetDaughterParticle(
                                                            Form("cutKx%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaKTPC),
                                                            trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kKaon,nsigmaKTPC,nsigmaKTOF);
    else if(CutTypeK==1) cutSetKx=new AliRsnCutSetDaughterParticle(
                                                                   Form("cutKx%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,nsigmaKTPC),
                                                                   trkQualityCut,AliRsnCutSetDaughterParticle::kFastTPCpidNsigma,AliPID::kKaon,nsigmaKTPC,-1.);
    else if(CutTypeK==2) cutSetKx=new AliRsnCutSetDaughterParticle(
                                                                   Form("cutKx%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,nsigmaKTOF),
                                                                   trkQualityCut,AliRsnCutSetDaughterParticle::kFastTOFpidNsigma,AliPID::kKaon,-1.,nsigmaKTOF);
    if(!cutSetKx){cerr<<"Error in AddTaskResonanceFinder::Config_kxkstar0(): missing cutSetKx"<<endl; return kFALSE;}
    
    //Int_t iCutQ=task->AddTrackCuts(cutSetQ);
    Int_t iCutPi=task->AddTrackCuts(cutSetPi);
    Int_t iCutKx=task->AddTrackCuts(cutSetKx);
    
    // set cuts for K0S from K*
    Int_t SidebandKstar=(TrackCutsK/10000)%10;
    Int_t K0sCuts=TrackCutsK%10000;
    
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
    cutK0s->SetMaxRapidity(2.);
    cutK0s->SetMinTPCcluster(-1);
    
    AliRsnCutSet* cutSetK0s=new AliRsnCutSet("setK0s",AliRsnTarget::kDaughter);
    cutSetK0s->AddCut(cutK0s);
    cutSetK0s->SetCutScheme(cutK0s->GetName());
    Int_t iCutK0s=task->AddTrackCuts(cutSetK0s);
    
    // monitoring
    TString pname="k0s";
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        //AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetKx->GetMonitorOutput());
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetK0s->GetMonitorOutput(), "K0S");
    }
    
    AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
    AliRsnCutMiniPair* cutYRes=new AliRsnCutMiniPair("cutRapidityResonance",AliRsnCutMiniPair::kRapidityRange);
    cutYRes->SetRangeD(-0.6,0.6);
    AliRsnMiniResonanceFinder* Kfinder[4];
    Int_t i,iCutKstar[4];
    
    AliRsnCutMiniPair* cutMassKstar=new AliRsnCutMiniPair("cutMassKstar",AliRsnCutMiniPair::kMassRange);
    cutMassKstar->SetRangeD(0.841,0.943);
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
    cutMassKstarSB->SetRangeD(0.966,1.014);
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
    
    // pair cuts
    AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
    if(system!=1) cutY->SetRangeD(-0.5,0.5);
    else cutY->SetRangeD(-0.465,0.035);
    
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
    
    Int_t k,cut1;
    TString name,comp;
    Char_t charge1,charge2;
    AliRsnMiniOutput* out;
    
    for(i=0;i<4;i++) for(j=0;j<2;j++) for(k=0;k<3;k++){
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
            name.Append("Kp");
            charge2='+';
            
        }else if(j==1){
            name.Append("Km");
            charge2='-';
        }
        
        if(!k){
            comp.Form("PAIR");
        }else if(k==1){
            name.Append("Mix");
            comp.Form("MIX");
        }else if(k==2){
            name.Append("Rotate");
            comp.Form("ROTATE1");
        }
        
        if(i>=2 && k) continue;
        
        out=task->CreateOutput(Form("kxkstarx_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kKstarpm);
        out->SetCutID(0,iCutKstar[cut1]);
        out->SetCharge(0,charge1);
        if(i<=1) out->SetUseStoredMass(0);
        
        out->SetDaughter(1,AliRsnDaughter::kKaon);
        out->SetCutID(1,iCutKx);
        out->SetCharge(1,charge2);
        
        if(k<=1) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(3124);
        out->SetMotherMass(mass);
        
        out->AddAxis(imID,250,1.25,2.5);// axis X: invmass or resolution
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
        
        out=task->CreateOutput(Form("kxkstarx_%s",name.Data()),"HIST",comp.Data());
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
    
    return kTRUE;
}


//=============================


Bool_t Config_k0kstarx(
                       AliRsnMiniAnalysisTask *task,
                       TString     lname,
                       Bool_t      isMC,
                       Int_t       system,
                       Int_t       EventCuts,
                       Int_t       TrackCutsK,
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
    
    Double_t mass=0.497611+0.89166;
    
    // set cuts for pions
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
    if(!cutSetPi){cerr<<"Error in AddTaskRare_pp13::Config_k0kstarx(): missing cutSetPi"<<endl; return kFALSE;}
    
    //Int_t iCutQ=task->AddTrackCuts(cutSetQ);
    Int_t iCutPi=task->AddTrackCuts(cutSetPi);
    
    // set cuts for K0S from K*
    Int_t SidebandKstar=(TrackCutsK/10000)%10;
    Int_t K0sCuts=TrackCutsK%10000;
    
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
    cutK0s->SetMaxRapidity(2.);
    cutK0s->SetMinTPCcluster(-1);
    
    AliRsnCutSet* cutSetK0s=new AliRsnCutSet("setK0s",AliRsnTarget::kDaughter);
    cutSetK0s->AddCut(cutK0s);
    cutSetK0s->SetCutScheme(cutK0s->GetName());
    Int_t iCutK0s=task->AddTrackCuts(cutSetK0s);
    
    // monitoring
    TString pname="k0s";
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        //AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetK0s->GetMonitorOutput(), "K0S");
    }
    
    AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
    AliRsnCutMiniPair* cutYRes=new AliRsnCutMiniPair("cutRapidityResonance",AliRsnCutMiniPair::kRapidityRange);
    cutYRes->SetRangeD(-0.6,0.6);
    AliRsnMiniResonanceFinder* Kfinder[4];
    Int_t i,iCutKstar[4];
    
    AliRsnCutMiniPair* cutMassKstar=new AliRsnCutMiniPair("cutMassKstar",AliRsnCutMiniPair::kMassRange);
    cutMassKstar->SetRangeD(0.841,0.943);
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
    cutMassKstarSB->SetRangeD(0.966,1.014);
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
    
    // pair cuts
    AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
    if(system!=1) cutY->SetRangeD(-0.5,0.5);
    else cutY->SetRangeD(-0.465,0.035);
    
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
    
    Int_t k,cut1;
    TString name,comp;
    Char_t charge1;
    AliRsnMiniOutput* out;
    
    for(i=0;i<4;i++) for(k=0;k<3;k++){
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
        
        name.Append("K0");
        
        if(!k){
            comp.Form("PAIR");
        }else if(k==1){
            name.Append("Mix");
            comp.Form("MIX");
        }else if(k==2){
            name.Append("Rotate");
            comp.Form("ROTATE1");
        }
        
        if(i>=2 && k) continue;
        
        out=task->CreateOutput(Form("k0kstarx_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kKstarpm);
        out->SetCutID(0,iCutKstar[cut1]);
        out->SetCharge(0,charge1);
        if(i<=1) out->SetUseStoredMass(0);
        
        out->SetDaughter(1,AliRsnDaughter::kKaon0);
        out->SetCutID(1,iCutK0s);
        out->SetCharge(1,'0');
        
        if(k<=1) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(3124);
        out->SetMotherMass(mass);
        
        out->AddAxis(imID,250,1.25,2.5);// axis X: invmass or resolution
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
        
        out=task->CreateOutput(Form("k0kstarx_%s",name.Data()),"HIST",comp.Data());
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
    
    return kTRUE;
}


//=============================


Bool_t Config_pkstarx(
                       AliRsnMiniAnalysisTask *task,
                       TString     lname,
                       Bool_t      isMC,
                       Int_t       system,
                       Int_t       EventCuts,
                       Int_t       TrackCutsK,
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
    
    Double_t mass=0.938272+0.89166;
    
    // set cuts for pions
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
    if(!cutSetPi){cerr<<"Error in AddTaskRare_pp13::Config_pkstarx(): missing cutSetPi"<<endl; return kFALSE;}
    
    // set cuts for protons
    int TrackCutsP=0;
    if(!(TrackCutsP%10000)) TrackCutsP+=3020;//default settings
    Float_t nsigmaPTPC=0.1*(TrackCutsP%100);
    Float_t nsigmaPTOF=0.1*((TrackCutsP/100)%100);
    Int_t CutTypeP=(TrackCutsP/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only
    
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
    if(!cutSetP){cerr<<"Error in AddTaskResonanceFinder::Config_pkstarx(): missing cutSetP"<<endl; return kFALSE;}
    
    //Int_t iCutQ=task->AddTrackCuts(cutSetQ);
    Int_t iCutPi=task->AddTrackCuts(cutSetPi);
    Int_t iCutP=task->AddTrackCuts(cutSetP);
    
    // set cuts for K0S from K*
    Int_t SidebandKstar=(TrackCutsK/10000)%10;
    Int_t K0sCuts=TrackCutsK%10000;
    
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
    cutK0s->SetMaxRapidity(2.);
    cutK0s->SetMinTPCcluster(-1);
    
    AliRsnCutSet* cutSetK0s=new AliRsnCutSet("setK0s",AliRsnTarget::kDaughter);
    cutSetK0s->AddCut(cutK0s);
    cutSetK0s->SetCutScheme(cutK0s->GetName());
    Int_t iCutK0s=task->AddTrackCuts(cutSetK0s);
    
    // monitoring
    TString pname="k0s";
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        //AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetP->GetMonitorOutput());
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetK0s->GetMonitorOutput(), "K0S");
    }
    
    AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
    AliRsnCutMiniPair* cutYRes=new AliRsnCutMiniPair("cutRapidityResonance",AliRsnCutMiniPair::kRapidityRange);
    cutYRes->SetRangeD(-0.6,0.6);
    AliRsnMiniResonanceFinder* Kfinder[4];
    Int_t i,iCutKstar[4];
    
    AliRsnCutMiniPair* cutMassKstar=new AliRsnCutMiniPair("cutMassKstar",AliRsnCutMiniPair::kMassRange);
    cutMassKstar->SetRangeD(0.841,0.943);
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
    cutMassKstarSB->SetRangeD(0.966,1.014);
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
    
    // pair cuts
    AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
    if(system!=1) cutY->SetRangeD(-0.5,0.5);
    else cutY->SetRangeD(-0.465,0.035);
    
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
    
    Int_t k,cut1;
    TString name,comp;
    Char_t charge1,charge2;
    AliRsnMiniOutput* out;
    
    for(i=0;i<4;i++) for(j=0;j<2;j++) for(k=0;k<3;k++){
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
            name.Append("Pp");
            charge2='+';
            
        }else if(j==1){
            name.Append("Pm");
            charge2='-';
        }
        
        if(!k){
            comp.Form("PAIR");
        }else if(k==1){
            name.Append("Mix");
            comp.Form("MIX");
        }else if(k==2){
            name.Append("Rotate");
            comp.Form("ROTATE1");
        }
        
        if(i>=2 && k) continue;
        
        out=task->CreateOutput(Form("pkstarx_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kKstarpm);
        out->SetCutID(0,iCutKstar[cut1]);
        out->SetCharge(0,charge1);
        if(i<=1) out->SetUseStoredMass(0);
        
        out->SetDaughter(1,AliRsnDaughter::kProton);
        out->SetCutID(1,iCutP);
        out->SetCharge(1,charge2);
        
        if(k<=1) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(3124);
        out->SetMotherMass(mass);
        
        out->AddAxis(imID,260,1.7,3);// axis X: invmass or resolution
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
        
        out=task->CreateOutput(Form("pkstarx_%s",name.Data()),"HIST",comp.Data());
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
    
    return kTRUE;
}


//=============================


Bool_t Config_Kstar0Lambda(
                           AliRsnMiniAnalysisTask *task,
                           TString     lname,
                           Bool_t      isMC,
                           Int_t       system,
                           Int_t       EventCuts,
                           Int_t       TrackCutsKstar,
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
    
    Double_t mass= 0.896+1.115683;
    
    // set cuts for pions and kaons
    Int_t TrackCutsPi=TrackCutsKstar%1000000;//Changed from TrackCutsS?
    if(!(TrackCutsPi%10000)) TrackCutsPi+=3020;//default settings
    Float_t nsigmaPiTPC=0.1*(TrackCutsPi%100);
    Float_t nsigmaPiTOF=0.1*((TrackCutsPi/100)%100);
    Int_t CutTypePi=(TrackCutsPi/10000)%100;//0=TPC+TOFveto (default), 1=TPC only, 2=TOF only, 3 TPC+TOFcut
    Int_t SidebandKstar=(TrackCutsKstar/1000000)%10;
    
    Int_t TrackCutsKx=TrackCutsK%10000;
    if(!(TrackCutsKx)) TrackCutsKx+=3020;//default settings
    Float_t nsigmaKTPC=0.1*(TrackCutsKx%100);
    Float_t nsigmaKTOF=0.1*((TrackCutsKx/100)%100);
    Int_t pairRotate=((TrackCutsK/10000)%10);
    
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
    
    
    // selections for V0 daughters
    Int_t v0d_xrows=70;
    Float_t v0d_rtpc=0.8;
    Float_t v0d_dcaxy=0.06;
    
    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts("qualityDaughterK0s");  //not sure about this line
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
    TString pname="lambdap";
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        //AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetK->GetMonitorOutput());
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetLambda->GetMonitorOutput(), "lambda");
        AddMonitorOutputV0(isMC,cutSetAntiLambda->GetMonitorOutput(), "antilambda");
    }
    
    // AliRsnMiniResonanceFinder - K*0
    AliRsnMiniResonanceFinder* Kfinder[6];
    Int_t i,iCutKstar[6];
    
    AliRsnCutMiniPair* cutYRes=new AliRsnCutMiniPair("cutPseudorapidityResonance",AliRsnCutMiniPair::kPseudorapidityRange);
    cutYRes->SetRangeD(-0.8,0.8);
    
    AliRsnCutMiniPair* cutMassKstar=new AliRsnCutMiniPair("cutMassKstar",AliRsnCutMiniPair::kMassRange);
    cutMassKstar->SetRangeD(0.841,0.943);
    AliRsnCutSet* cutsKstar=new AliRsnCutSet("pairCutsKstar",AliRsnTarget::kMother);
    cutsKstar->AddCut(cutMassKstar);
    cutsKstar->AddCut(cutYRes);
    cutsKstar->SetCutScheme(TString::Format("%s&%s",cutMassKstar->GetName(),cutYRes->GetName()).Data());
    
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
    if(isMC) Kfinder[i]->SetPairMode(1);
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
    if(isMC) Kfinder[i]->SetPairMode(1);
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
    
    AliRsnCutMiniPair* cutMassKstarSB=new AliRsnCutMiniPair("cutMassKstarSB",AliRsnCutMiniPair::kMassRange);
    if(SidebandKstar<=1) cutMassKstarSB->SetRangeD(0.766,0.816);
    else cutMassKstarSB->SetRangeD(0.966,1.014);
    AliRsnCutSet* cutsKstarSB=new AliRsnCutSet("pairCutsKstarSB",AliRsnTarget::kMother);
    cutsKstarSB->AddCut(cutMassKstarSB);
    cutsKstarSB->AddCut(cutYRes);
    cutsKstarSB->SetCutScheme(TString::Format("%s&%s",cutMassKstarSB->GetName(),cutYRes->GetName()).Data());
    
    i=4; // K*0 Sideband
    Kfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinderSB",task->GetName()));
    Kfinder[i]->SetCutID(0,iCutK);
    Kfinder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
    Kfinder[i]->SetCharge(0,'+');
    Kfinder[i]->SetCutID(1,iCutPi);
    Kfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
    Kfinder[i]->SetCharge(1,'-');
    Kfinder[i]->SetResonanceMass(0.895810);
    Kfinder[i]->SetResonancePDG(313);
    Kfinder[i]->SetPairCuts(cutsKstarSB);
    iCutKstar[i]=task->AddResonanceFinder(Kfinder[i]);
    
    i=5; //Anti-K*0 Sideband
    Kfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_AntiResonanceFinderSB",task->GetName()));
    Kfinder[i]->SetCutID(0,iCutK);
    Kfinder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
    Kfinder[i]->SetCharge(0,'-');
    Kfinder[i]->SetCutID(1,iCutPi);
    Kfinder[i]->SetDaughter(1,AliRsnDaughter::kPion);
    Kfinder[i]->SetCharge(1,'+');
    Kfinder[i]->SetResonanceMass(0.895810);
    Kfinder[i]->SetResonancePDG(-313);
    Kfinder[i]->SetPairCuts(cutsKstarSB);
    iCutKstar[i]=task->AddResonanceFinder(Kfinder[i]);
    
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
    
    Int_t k,xID,cut1,cut2,pairID,ipdg;
    TString name,comp;
    AliRsnMiniOutput* out;
    
    task->SetMotherAcceptanceCutMinPt(0.15);
    task->SetMotherAcceptanceCutMaxEta(0.8);
    task->KeepMotherInAcceptance(true);
    
    for(i=0;i<6;i++) for(j=0;j<2;j++) for(k=0;k<9;k++){
        if(!i){
            name.Form("Kstar0p");
            cut2=0;
        }else if(i==1){
            name.Form("Kstar0a");
            cut2=1;
        }else if(i==2){
            name.Form("KpPip");
            cut2=2;
        }else if(i==3){
            name.Form("KmPim");
            cut2=3;
        }else if(i==4){
            name.Form("Kstar0pSB");
            cut2=4;
        }else if(i==5){
            name.Form("Kstar0aSB");
            cut2=5;
        }
        
        if(!j){
            name.Append("Lambdap");
            cut1=iCutLambda;
            ipdg=3124;
        }else if(j==1){
            name.Append("Lambdaa");
            cut1=iCutAntiLambda;
            ipdg=-3124;
        }
        
        if(!isMC && k>2) continue;
        if(i>1 && k) continue;
        
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
            if(!isMC) continue;
            name.Append("_genPS");
            comp.Form("MOTHER_IN_ACC");
        }else if(k==8){
            if(!isMC) continue;
            name.Append("_recPS");
            comp.Form("TRUE");
        }
        
        out=task->CreateOutput(Form("Kstar0Lambda_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kLambda);
        out->SetCutID(0,cut1);
        out->SetCharge(0,'0');
        
        out->SetDaughter(1,AliRsnDaughter::kKstar0);
        out->SetCutID(1,iCutKstar[cut2]);
        out->SetCharge(1,'0');
        if(cut2!=4 && cut2!=5) out->SetUseStoredMass(1);
        
        if(!pairID) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(ipdg);
        out->SetMotherMass(mass);
        
        if(k<=6){
            if(xID==imID || xID==mmID) out->AddAxis(xID,240,1.8,3.0);// axis X: invmass or resolution
            else out->AddAxis(diffID,200,-0.02,0.02);
            out->AddAxis(ptID,200,0.0,20.0);// axis Y: transverse momentum
            out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
        }else{
            out->AddAxis(fdpt,100,0.,10.);
            out->AddAxis(sdpt,100,0.,10.);
            out->AddAxis(ptID,40,0.,20.);
        }
    }
    
    // fill monitoring histogram for the K*0
    for(i=0;i<6;i++) for(j=0;j<4;j++){
        if(!i) name.Form("Kstar0p");
        else if(i==1) name.Form("Kstar0a");
        else if(i==2) name.Form("KpPip");
        else if(i==3) name.Form("KmPim");
        else if(i==4) name.Form("Kstar0pSB");
        else if(i==5) name.Form("Kstar0aSB");
        
        xID=imID;
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
        }else if(j==3){
            if(!isMC) continue;
            name.Append("_recres");
            comp.Form("TRUE");
            xID=diffID;
        }
        
        if(i>=2 && (j || SidebandKstar)) continue;
        
        out=task->CreateOutput(Form("Kstar0Lambda_%s",name.Data()),"HIST",comp.Data());
        out->SetMotherPDG(Kfinder[i]->GetResonancePDG());
        
        out->SetDaughter(0,Kfinder[i]->GetDaughter(0));
        out->SetCutID(0,Kfinder[i]->GetCutID(0));
        out->SetCharge(0,Kfinder[i]->GetCharge(0));
        
        out->SetDaughter(1,Kfinder[i]->GetDaughter(1));
        out->SetCutID(1,Kfinder[i]->GetCutID(1));
        out->SetCharge(1,Kfinder[i]->GetCharge(1));
        
        out->SetMotherMass(Kfinder[i]->GetResonanceMass());
        if(i<=3) out->SetPairCuts(cutsKstar);
        else out->SetPairCuts(cutsKstarSB);
        
        if(xID==imID) out->AddAxis(imID,170,0.75,1.09);
        else out->AddAxis(diffID,200,-0.02,0.02);
        out->AddAxis(ptID,50,0.0,20.0);
        out->AddAxis(centID,nmult,multbins);
    }
    
    //for efficiency of (anti-)Lambda
    if(isMC){
        AliRsnCutMiniPair* cutYV0=new AliRsnCutMiniPair("cutPseudorapidityV0", AliRsnCutMiniPair::kPseudorapidityRange);
        cutYV0->SetRangeD(-0.8,0.8);
        
        AliRsnCutSet* cutsPairV0=new AliRsnCutSet("pairCutsV0", AliRsnTarget::kMother);
        cutsPairV0->AddCut(cutYV0);
        cutsPairV0->SetCutScheme(cutYV0->GetName());
        
        out = task->CreateOutput("Kstar0Lambda_Lambdap_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kProton);
        out->SetMotherPDG(3122);
        out->SetMotherMass(1.115683);
        out->SetPairCuts(cutsPairV0);
        out->AddAxis(ptID,200,0.0,20.0);
        
        out = task->CreateOutput("Kstar0Lambda_Lambdaa_mother", "HIST", "MOTHER");
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


Bool_t Config_KstarxLambda(
                           AliRsnMiniAnalysisTask *task,
                           TString     lname,
                           Bool_t      isMC,
                           Int_t       system,
                           Int_t       EventCuts,
                           Int_t       TrackCutsLambda,
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
    
    Double_t mass=0.892+1.115683;
    
    // set cuts for pions
    Int_t TrackCutsPi=TrackCutsK%1000000;
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
    if(!cutSetPi){cerr<<"Error in AddTaskRare_pp13::Config_kstarxLambda(): missing cutSetPi"<<endl; return kFALSE;}
    
    //Int_t iCutQ=task->AddTrackCuts(cutSetQ);
    Int_t iCutPi=task->AddTrackCuts(cutSetPi);
    
    // set cuts for K0S from K*
    
    Int_t SidebandKstar=(TrackCutsK/1000000)%10;
    Int_t K0sCuts=(TrackCutsK/10000000)%10;
    Int_t pairRotate=(TrackCutsLambda%10);
    
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
    TString pname="k0s";
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        //AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());
        
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutputV0(isMC,cutSetK0s->GetMonitorOutput(), "K0S");
        AddMonitorOutputV0(isMC,cutSetLambda->GetMonitorOutput(), "lambda");
        AddMonitorOutputV0(isMC,cutSetAntiLambda->GetMonitorOutput(), "antilambda");
    }
    
    // AliRsnMiniResonanceFinder - K*
    AliRsnCutMiniPair* cutV0=new AliRsnCutMiniPair("cutV0",AliRsnCutMiniPair::kContainsV0Daughter);
    AliRsnCutMiniPair* cutYRes=new AliRsnCutMiniPair("cutPseudorapidityResonance",AliRsnCutMiniPair::kPseudorapidityRange);
    cutYRes->SetRangeD(-0.8,0.8);
    
    AliRsnMiniResonanceFinder* Kfinder[4];
    Int_t i,iCutKstar[4];
    
    AliRsnCutMiniPair* cutMassKstar=new AliRsnCutMiniPair("cutMassKstar",AliRsnCutMiniPair::kMassRange);
    cutMassKstar->SetRangeD(0.841,0.943);
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
    if(isMC) Kfinder[i]->SetPairMode(1);
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
    if(isMC) Kfinder[i]->SetPairMode(1);
    iCutKstar[i]=task->AddResonanceFinder(Kfinder[i]);
    
    AliRsnCutMiniPair* cutMassKstarSB=new AliRsnCutMiniPair("cutMassKstarSB",AliRsnCutMiniPair::kMassRange);
    if(SidebandKstar<=1) cutMassKstarSB->SetRangeD(0.766,0.816);
    else cutMassKstarSB->SetRangeD(0.966,1.014);
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
    
    // pair cuts
    AliRsnCutMiniPair* cutY=new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
    if(system!=1) cutY->SetRangeD(-0.5,0.5);
    else cutY->SetRangeD(-0.465,0.035);
    
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
    
    Int_t k,xID,cut1,cut2,pairID,ipdg;
    AliRsnDaughter::ESpecies d2=AliRsnDaughter::kUnknown;
    TString name,comp;
    Char_t charge1;
    AliRsnMiniOutput* out;
    
    task->SetMotherAcceptanceCutMinPt(0.15);
    task->SetMotherAcceptanceCutMaxEta(0.8);
    task->KeepMotherInAcceptance(true);
    
    for(i=0;i<4;i++) for(j=0;j<2;j++) for(k=0;k<9;k++){
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
            name.Append("Lambdap");
            d2=AliRsnDaughter::kLambda;
            cut2=iCutLambda;
        }else if(j==1){
            name.Append("Lambdaa");
            d2=AliRsnDaughter::kLambda;
            cut2=iCutAntiLambda;
        }
        
        xID=imID;
        pairID=1;
        if(!k){
            comp.Form("PAIR");
            pairID=0;
        }else if(k==1){
            name.Append("Mix");
            comp.Form("MIX");
        }else if (k==2){
            name.Append("ROTATE");
            if(!pairRotate) {comp.Form("ROTATE1");}
            else {comp.Form("ROTATE2");}
            pairID=0;
        }else if(k==3){
            if(!isMC) continue;
            name.Append("_gen");
            comp.Form("MOTHER");
        }else if(k==4){
            if(!isMC) continue;
            name.Append("_rec");
            comp.Form("TRUE");
        }else if(k==5){
            if(!isMC) continue;
            name.Append("_recMM");
            comp.Form("TRUE");
            xID=mmID;
        }else if(k==6){
            if(!isMC) continue;
            name.Append("_res");
            comp.Form("TRUE");
            xID=diffID;
        }else if(k==7){
            if(!isMC) continue;
            name.Append("_genPS");
            comp.Form("MOTHER_IN_ACC");
        }else if(k==8){
            if(!isMC) continue;
            name.Append("_recPS");
            comp.Form("TRUE");
        }
        
        if(i>=2 && k) continue;
        
        out=task->CreateOutput(Form("KstarxLambda_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kKstarpm);
        out->SetCutID(0,iCutKstar[cut1]);
        out->SetCharge(0,charge1);
        if(i<=1) out->SetUseStoredMass(0);
        
        out->SetDaughter(1,d2);
        out->SetCutID(1,cut2);
        out->SetCharge(1,'0');
        
        if(k!=1) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(3124);
        out->SetMotherMass(mass);
        
        if(k<=6){
            if(xID==imID || xID==mmID) out->AddAxis(xID,240,1.8,3.0);// axis X: invmass or resolution
            else out->AddAxis(diffID,200,-0.02,0.02);
            out->AddAxis(ptID,50,0.0,20.0);// axis Y: transverse momentum
            out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
        }else{//Phase-space histograms
            out->AddAxis(fdpt,100,0.,10.);
            out->AddAxis(sdpt,100,0.,10.);
            out->AddAxis(ptID,40,0.,20.);
        }
    }
    
    // fill monitoring histogram for the K*
    for(i=0;i<4;i++) for(j=0;j<4;j++){
        if(!i) name.Form("Kstarp");
        else if(i==1) name.Form("Kstarm");
        else if(i==2) name.Form("KstarpSB");
        else if(i==3) name.Form("KstarmSB");

        xID=imID;
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
        }else if(j==3){
            if(!isMC) continue;
            name.Append("_recres");
            comp.Form("TRUE");
            xID=diffID;
        }
        
        if(i>=2 && j) continue;
        
        out=task->CreateOutput(Form("KstarxLambda_%s",name.Data()),"HIST",comp.Data());
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
        
        if(xID==imID) out->AddAxis(imID,170,0.75,1.09);
        else out->AddAxis(diffID,200,-0.02,0.02);
        out->AddAxis(ptID,50,0.0,20.0);
        out->AddAxis(centID,nmult,multbins);
    }
    
    //for efficiency of (anti-)Lambda
    if(isMC){
        AliRsnCutMiniPair* cutYV0=new AliRsnCutMiniPair("cutPseudorapidityV0", AliRsnCutMiniPair::kPseudorapidityRange);
        cutYV0->SetRangeD(-0.8,0.8);
        
        AliRsnCutSet* cutsPairV0=new AliRsnCutSet("pairCutsV0", AliRsnTarget::kMother);
        cutsPairV0->AddCut(cutYV0);
        cutsPairV0->SetCutScheme(cutYV0->GetName());
        
        out = task->CreateOutput("KstarxLambda_Lambdap_mother", "HIST", "MOTHER");
        out->SetDaughter(0, AliRsnDaughter::kPion);
        out->SetDaughter(1, AliRsnDaughter::kProton);
        out->SetMotherPDG(3122);
        out->SetMotherMass(1.115683);
        out->SetPairCuts(cutsPairV0);
        out->AddAxis(ptID,200,0.0,20.0);
        
        out = task->CreateOutput("KstarxLambda_Lambdaa_mother", "HIST", "MOTHER");
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


Bool_t Config_k0kxpi(
                       AliRsnMiniAnalysisTask *task,
                       TString     lname,
                       Bool_t      isMC,
                       Int_t       system,
                       Int_t       EventCuts,
                       Int_t       TrackCutsK0,
                       Int_t       TrackCutsKstar0 //K+/- pi+/- pair
){
    bool isPP=false;
    if(!system) isPP=true;
    int trigger=EventCuts%10;
    int MultBins=(EventCuts/10)%10;
    if(system==1 || system==2) MultBins=1;

    char suffix[1000];
    sprintf(suffix,"_%s",lname.Data());
    Bool_t enableMonitor=kTRUE;

    Double_t mass=1.2819;//f1(1285)

    // set daughter cuts :
    Int_t TrackCutsPi=0;
    if(!(TrackCutsPi%10000)) TrackCutsPi+=3020;//default settings
    Float_t nsigmaPiTPC=0.1*(TrackCutsPi%100);
    Float_t nsigmaPiTOF=0.1*((TrackCutsPi/100)%100);
    Int_t CutTypePi=(TrackCutsPi/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only

    if(!(TrackCutsKstar0%10000)) TrackCutsKstar0+=3020;//default settings
    Float_t nsigmaKTPC=0.1*(TrackCutsKstar0%100);
    Float_t nsigmaKTOF=0.1*((TrackCutsKstar0/100)%100);
    Int_t CutTypeK=(TrackCutsKstar0/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only

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
    if(!cutSetPi){cerr<<"Error in AddTaskResonanceFinder::Config_k0kxpi(): missing cutSetPi"<<endl; return kFALSE;}

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
    Float_t k0s_massTolSigma=4;
    Int_t   k0s_massTolID=0;
    Float_t k0s_massTol=0.03;
    Float_t k0s_massTolVeto=0.004;
    Bool_t  k0sSwitch=kFALSE;
    Float_t k0sCosPoinAn=0.97;

    if(TrackCutsK0==1) k0s_massTolID=1;//use pT-dependent mass tolerance cut

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
    if(!cutSetK){cerr<<"Error in AddTaskResonanceFinder::Config_k0kxpi(): missing cutSetK"<<endl; return kFALSE;}

    Int_t iCutQ=task->AddTrackCuts(cutSetQ);
    Int_t iCutPi=task->AddTrackCuts(cutSetPi);
    Int_t iCutK=task->AddTrackCuts(cutSetK);

    // monitoring
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddMonitorOutputV0.C");
        AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetPi->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetK->GetMonitorOutput());
        AddMonitorOutputV0(isMC,cutSetK0s->GetMonitorOutput(),"K0S");
    }

    // AliRsnMiniResonanceFinder
    AliRsnCutMiniPair* cutMassKxK0=new AliRsnCutMiniPair("cutMassKxK0",AliRsnCutMiniPair::kMassRange);
    cutMassKxK0->SetRangeD(0,1.07);
    AliRsnCutSet* cutsKxK0=new AliRsnCutSet("pairCutsKxK0",AliRsnTarget::kMother);
    cutsKxK0->AddCut(cutMassKxK0);
    cutsKxK0->SetCutScheme(TString::Format("%s",cutMassKxK0->GetName()));

    AliRsnMiniResonanceFinder* pairfinder[2];

    // K+ K0S pairs
    Int_t i=0,iCutPair[2];
    pairfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_KpK0Finder",task->GetName()));
    pairfinder[i]->SetCutID(0,iCutK);
    pairfinder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
    pairfinder[i]->SetCharge(0,'+');
    pairfinder[i]->SetCutID(1,iCutK0s);
    pairfinder[i]->SetDaughter(1,AliRsnDaughter::kKaon0);
    pairfinder[i]->SetCharge(1,'0');
    pairfinder[i]->SetResonanceMass(0.98);
    pairfinder[i]->SetResonancePDG(9000211);
    pairfinder[i]->SetPairCuts(cutsKxK0);
    iCutPair[i]=task->AddResonanceFinder(pairfinder[i]);

    // K- K0S pairs
    i=1;
    pairfinder[i]=new AliRsnMiniResonanceFinder(Form("%s_KmK0Finder",task->GetName()));
    pairfinder[i]->SetCutID(0,iCutK);
    pairfinder[i]->SetDaughter(0,AliRsnDaughter::kKaon);
    pairfinder[i]->SetCharge(0,'-');
    pairfinder[i]->SetCutID(1,iCutK0s);
    pairfinder[i]->SetDaughter(1,AliRsnDaughter::kKaon0);
    pairfinder[i]->SetCharge(1,'0');
    pairfinder[i]->SetResonanceMass(0.98);
    pairfinder[i]->SetResonancePDG(-9000211);
    pairfinder[i]->SetPairCuts(cutsKxK0);
    iCutPair[i]=task->AddResonanceFinder(pairfinder[i]);

    // triplet cuts
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
    /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
    /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt,         kFALSE);
    /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult,       kFALSE);
    /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta,        kFALSE);
    /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY,          kFALSE);

    // -- Create all needed outputs -----------------------------------------------------------------
    // use an array for more compact writing, which are different on mixing and charges

    Int_t xID,cut2,pairID;
    Char_t charge1;
    TString name,comp;
    AliRsnMiniOutput* out;

    for(i=0;i<12;i++){
        if(!i){
            xID=imID;
            name.Form("Pim_KpK0");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutPair[0];
            pairID=0;
        }else if(i==1){
            xID=imID;
            name.Form("Pip_KmK0");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutPair[1];
            pairID=0;
        }else if(i==2){
            xID=imID;
            name.Form("Pip_KpK0");
            comp.Form("PAIR");
            charge1='+';
            cut2=iCutPair[0];
            pairID=0;
        }else if(i==3){
            xID=imID;
            name.Form("Pim_KmK0");
            comp.Form("PAIR");
            charge1='-';
            cut2=iCutPair[1];
            pairID=0;
        }else if(i==4){
            xID=imID;
            name.Form("Pim_KpK0Mix");
            comp.Form("MIX");
            charge1='-';
            cut2=iCutPair[0];
            pairID=1;
        }else if(i==5){
            xID=imID;
            name.Form("Pip_KmK0Mix");
            comp.Form("MIX");
            charge1='+';
            cut2=iCutPair[1];
            pairID=1;
        }else if(i==6){
            xID=imID;
            name.Form("Pip_KpK0Mix");
            comp.Form("MIX");
            charge1='+';
            cut2=iCutPair[0];
            pairID=1;
        }else if(i==7){
            xID=imID;
            name.Form("Pim_KmK0Mix");
            comp.Form("MIX");
            charge1='-';
            cut2=iCutPair[1];
            pairID=1;
        }else if(i==8){
            xID=imID;
            name.Form("Pim_KpK0Rotated");
            comp.Form("ROTATE1");
            charge1='-';
            cut2=iCutPair[0];
            pairID=0;
        }else if(i==9){
            xID=imID;
            name.Form("Pip_KmK0Rotated");
            comp.Form("ROTATE1");
            charge1='+';
            cut2=iCutPair[1];
            pairID=0;
        }else if(i==10){
            xID=imID;
            name.Form("Pip_KpK0Rotated");
            comp.Form("ROTATE1");
            charge1='+';
            cut2=iCutPair[0];
            pairID=0;
        }else if(i==11){
            xID=imID;
            name.Form("Pim_KmK0Rotated");
            comp.Form("ROTATE1");
            charge1='-';
            cut2=iCutPair[1];
            pairID=0;
        }

        out=task->CreateOutput(Form("k0kxpi_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kPion);
        out->SetCutID(0,iCutPi);
        out->SetCharge(0,charge1);

        out->SetDaughter(1,AliRsnDaughter::kKstarpm);//dummy
        out->SetCutID(1,cut2);
        if(cut2==iCutPair[0]) out->SetCharge(1,'+');
        else out->SetCharge(1,'-');
        out->SetUseStoredMass(1);

        if(!pairID) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(20223);
        out->SetMotherMass(mass);
        
        if(xID==imID) out->AddAxis(imID,280,1.1,2.5);
        else out->AddAxis(resID,200,-0.02,0.02);
        out->AddAxis(ptID,200,0.0,20.0);// axis Y: transverse momentum
        out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
    }

    return kTRUE;
}

//=============================


Bool_t Config_KxLambdastar(
                           AliRsnMiniAnalysisTask *task,
                           TString     lname,
                           Bool_t      isMC,
                           Int_t       system,
                           Int_t       EventCuts,
                           Int_t       TrackCutsK,
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
    
    Double_t mass=0.493677+1.5195;
    
    // set cuts for primary bachelor kaon
    if(!(TrackCutsK%10000)) TrackCutsK+=3020;//default settings
    Float_t nsigmaKTPC=0.1*(TrackCutsK%100);
    Float_t nsigmaKTOF=0.1*((TrackCutsK/100)%100);
    Int_t pairRotate=(TrackCutsK/100000)%10;
    
    AliRsnCutTrackQuality* trkQualityCut=new AliRsnCutTrackQuality("myQualityCut");
    trkQualityCut->SetDefaults2011(kTRUE,kTRUE);
    
    AliRsnCutSetDaughterParticle* cutSetQ=new AliRsnCutSetDaughterParticle("cutQ",trkQualityCut,AliRsnCutSetDaughterParticle::kQualityStd2010,AliPID::kPion,-1.);
    AliRsnCutSetDaughterParticle* cutSetK=new AliRsnCutSetDaughterParticle(Form("cutK%i_%2.1fsigma",AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,nsigmaKTPC),
                                                                           trkQualityCut,AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015,AliPID::kKaon,nsigmaKTPC,nsigmaKTOF);
    
    
    
    // set cuts for primary proton from Lambda(1520)
    Int_t TrackCutsP=TrackCutsLambda/10000;
    if(!(TrackCutsP%10000)) TrackCutsP+=3020;//default settings
    Float_t nsigmaPTPC=0.1*(TrackCutsP%100);
    Float_t nsigmaPTOF=0.1*((TrackCutsP/100)%100);
    Int_t CutTypeP=(TrackCutsP/10000)%100;//0=TPC+TOFveto (default), 1=TPC only, 2=TOF only, 3 TPC+TOFcut
    
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
    if(!cutSetP){cerr<<"Error in AddTaskResonanceFinder::Config_KxLambdastar(): missing cutSetP"<<endl; return kFALSE;}
    
    Int_t iCutP=task->AddTrackCuts(cutSetP);
    Int_t iCutQ=task->AddTrackCuts(cutSetQ);
    Int_t iCutK=task->AddTrackCuts(cutSetK);
    
    
    // monitoring
    if(enableMonitor){
        Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
        AddMonitorOutput(isMC,cutSetQ->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetK->GetMonitorOutput());
        AddMonitorOutput(isMC,cutSetP->GetMonitorOutput());
    }
    
    // AliRsnMiniResonanceFinder
    
    AliRsnMiniResonanceFinder* finder[6];
    Int_t i,iCutLambdaStar[6];
    
    AliRsnCutMiniPair* cutMassL=new AliRsnCutMiniPair("cutMassLambdaStar",AliRsnCutMiniPair::kMassRange);
    cutMassL->SetRangeD(1.503,1.536);   //Mass +/- one gamma
    AliRsnCutMiniPair* cutYL=new AliRsnCutMiniPair("cutRapidityLambdaStar",AliRsnCutMiniPair::kRapidityRange);
    cutYL->SetRangeD(-0.6,0.6);
    AliRsnCutSet* cutsL=new AliRsnCutSet("pairCutsLambdaStar",AliRsnTarget::kMother);
    cutsL->AddCut(cutMassL);
    cutsL->AddCut(cutYL);
    cutsL->SetCutScheme(TString::Format("%s&%s&(!%s)",cutMassL->GetName(),cutYL->GetName()).Data());
    
    i=0; // Lambda(1520)
    finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_LambdaStar",task->GetName()));
    finder[i]->SetCutID(0,iCutP);
    finder[i]->SetDaughter(0,AliRsnDaughter::kProton);
    finder[i]->SetCharge(0,'+');
    finder[i]->SetCutID(1,iCutK);
    finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
    finder[i]->SetCharge(1,'-');
    finder[i]->SetResonanceMass(1.5195);
    finder[i]->SetResonancePDG(3124);
    finder[i]->SetPairCuts(cutsL);
    iCutLambdaStar[i]=task->AddResonanceFinder(finder[i]);
    
    i=1; // anti-Lambda1520
    finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_LambdaStar",task->GetName()));
    finder[i]->SetCutID(0,iCutP);
    finder[i]->SetDaughter(0,AliRsnDaughter::kProton);
    finder[i]->SetCharge(0,'+');
    finder[i]->SetCutID(1,iCutK);
    finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
    finder[i]->SetCharge(1,'-');
    finder[i]->SetResonanceMass(1.5195);
    finder[i]->SetResonancePDG(-3124);
    finder[i]->SetPairCuts(cutsL);
    iCutLambdaStar[i]=task->AddResonanceFinder(finder[i]);
    
    i=2; // P+ K+
    finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_PpKp",task->GetName()));
    finder[i]->SetCutID(0,iCutP);
    finder[i]->SetDaughter(0,AliRsnDaughter::kProton);
    finder[i]->SetCharge(0,'+');
    finder[i]->SetCutID(1,iCutK);
    finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
    finder[i]->SetCharge(1,'+');
    finder[i]->SetResonanceMass(1.5195);
    finder[i]->SetResonancePDG(3124);
    finder[i]->SetPairCuts(cutsL);
    iCutLambdaStar[i]=task->AddResonanceFinder(finder[i]);
    
    i=3; //P- K-
    finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_KmPm",task->GetName()));
    finder[i]->SetCutID(0,iCutP);
    finder[i]->SetDaughter(0,AliRsnDaughter::kProton);
    finder[i]->SetCharge(0,'-');
    finder[i]->SetCutID(1,iCutK);
    finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
    finder[i]->SetCharge(1,'-');
    finder[i]->SetResonanceMass(1.5195);
    finder[i]->SetResonancePDG(-3124);
    finder[i]->SetPairCuts(cutsL);
    iCutLambdaStar[i]=task->AddResonanceFinder(finder[i]);
    
    // sidebands
    AliRsnCutMiniPair* cutMassSB=new AliRsnCutMiniPair("cutMassSB",AliRsnCutMiniPair::kMassRange);
    cutMassSB->SetRangeD(1.480,1.496);
    AliRsnCutSet* cutsSB=new AliRsnCutSet("pairCutsSB",AliRsnTarget::kMother);
    cutsSB->AddCut(cutMassSB);
    cutsSB->AddCut(cutYS);
    cutsSB->SetCutScheme(TString::Format("%s&%s&(!%s)",cutMassSB->GetName(),cutYS->GetName()).Data());
    
    i=4; // Lambda(1520) sideband
    finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBLambdaStar",task->GetName()));
    finder[i]->SetCutID(0,iCutP);
    finder[i]->SetDaughter(0,AliRsnDaughter::kProton);
    finder[i]->SetCharge(0,'+');
    finder[i]->SetCutID(1,iCutK);
    finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
    finder[i]->SetCharge(1,'-');
    finder[i]->SetResonanceMass(1.5195);
    finder[i]->SetResonancePDG(3124);
    finder[i]->SetPairCuts(cutsSB);
    iCutLambdaStar[i]=task->AddResonanceFinder(finder[i]);
    
    i=5; // anti-Lambda(1520) sideband
    finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBAntiLambdaStar",task->GetName()));
    finder[i]->SetCutID(0,iCutP);
    finder[i]->SetDaughter(0,AliRsnDaughter::kProton);
    finder[i]->SetCharge(0,'+');
    finder[i]->SetCutID(1,iCutK);
    finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
    finder[i]->SetCharge(1,'-');
    finder[i]->SetResonanceMass(1.5195);
    finder[i]->SetResonancePDG(-3124);
    finder[i]->SetPairCuts(cutsSB);
    iCutLambdaStar[i]=task->AddResonanceFinder(finder[i]);
    
    
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
    
    Int_t k,xID,cut2,pairID,ipdg;
    AliRsnDaughter::ESpecies d2=AliRsnDaughter::kUnknown;
    TString name,comp;
    Char_t charge1,charge2;
    AliRsnMiniOutput* out;
    
    for(i=0;i<2;i++) for(j=0;j<2;j++) for(k=0;k<4;k++){
        if(!i){
            name.Form("Kp");
            charge1='+';
        }else{
            name.Form("Km");
            charge1='-';
        }
        
        if(!j){
            name.Append("Lambdastarp");
            cut2=0;
        }else if(j==1){
            name.Append("Lambdastara");
            cut2=1;
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
        }else if (k==3){
            name.Append("ROTATE");
            if(!pairRotate) {comp.Form("ROTATE1");}
            else {comp.Form("ROTATE2");}
        }
        
        out=task->CreateOutput(Form("kxLambdastar_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kKaon);
        out->SetCutID(0,iCutK);
        out->SetCharge(0,charge1);
        
        out->SetDaughter(1,AliRsnDaughter::kLambdastar);
        out->SetCutID(1,iCutLambdaStar[cut2]);
        out->SetCharge(1,0);
        if(k!=1) out->SetUseStoredMass(1);
        
        if(k!=2) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(3124);
        out->SetMotherMass(mass);
        
        out->AddAxis(imID,240,1.8,3);// axis X: invmass or resolution
        //out->AddAxis(resID,200,-0.02,0.02);
        out->AddAxis(ptID,50,0.0,20.0);// axis Y: transverse momentum
        out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
    }
    
    // fill monitoring histogram for the resonance (Lambda(1520))
    for(i=0;i<2;i++) for(j=0;j<4;j++){
        if(!i) name.Form("Lambdastarp");
        else if(i==1) name.Form("Lambdastarm");
        
        
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
        
        out=task->CreateOutput(Form("kxLambdastar_%s",name.Data()),"HIST",comp.Data());
        out->SetMotherPDG(finder[k]->GetResonancePDG());
        
        out->SetDaughter(0,finder[k]->GetDaughter(0));
        out->SetCutID(0,finder[k]->GetCutID(0));
        out->SetCharge(0,finder[k]->GetCharge(0));
        
        out->SetDaughter(1,finder[k]->GetDaughter(1));
        out->SetCutID(1,finder[k]->GetCutID(1));
        out->SetCharge(1,finder[k]->GetCharge(1));
        
        out->SetMotherMass(finder[k]->GetResonanceMass());
        if(j!=1) out->SetPairCuts(cutsL);
        else out->SetPairCuts(cutsSB);
        
        out->AddAxis(imID,150,1.3,1.6);
        out->AddAxis(ptID,50,0.0,20.0);
        out->AddAxis(centID,nmult,multbins);
    }
    
    return kTRUE;
}



//=============================


Bool_t Config_K0Lambdastar(
                           AliRsnMiniAnalysisTask *task,
                           TString     lname,
                           Bool_t      isMC,
                           Int_t       system,
                           Int_t       EventCuts,
                           Int_t       TrackCutsP,
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
    
    Double_t mass=0.497611+1.5195;
    
    // set cuts for kaon from Lambda(1520) decay
    if(!(TrackCutsK%10000)) TrackCutsK+=3020;//default settings
    Float_t nsigmaKTPC=0.1*(TrackCutsK%100);
    Float_t nsigmaKTOF=0.1*((TrackCutsK/100)%100);
    Int_t CutTypeK=(TrackCutsK/10000)%10;//0=TPC+TOF (default), 1=TPC only, 2=TOF only
    Int_t K0sCuts = (TrackCutsK/1000000)%10;
    
    //Integer for Rotated background, 0=ROTATE1, 1=ROTATE2
    Int_t pairRotate=(TrackCutsK/100000)%10;
    
    // set cuts for proton from Lambda(1520) decay
    if(!(TrackCutsP%10000)) TrackCutsP+=3020;//default settings
    Float_t nsigmaPTPC=0.1*(TrackCutsP%100);
    Float_t nsigmaPTOF=0.1*((TrackCutsP/100)%100);
    Int_t CutTypeP=(TrackCutsP/10000)%100;//0=TPC+TOFveto (default), 1=TPC only, 2=TOF only, 3 TPC+TOFcut
    
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
    if(!cutSetP){cerr<<"Error in AddTaskResonanceFinder::Config_k0Lambdastar(): missing cutSetP"<<endl; return kFALSE;}
    
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
        AddMonitorOutput(isMC,cutSetK->GetMonitorOutput());
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
    
    // AliRsnMiniResonanceFinder
    
    AliRsnMiniResonanceFinder* finder[6];
    Int_t i,iCutLambdaStar[6];
    
    AliRsnCutMiniPair* cutMassL=new AliRsnCutMiniPair("cutMassLambdaStar",AliRsnCutMiniPair::kMassRange);
    cutMassL->SetRangeD(1.503,1.536);   //Mass +/- one gamma
    AliRsnCutMiniPair* cutYL=new AliRsnCutMiniPair("cutRapidityLambdaStar",AliRsnCutMiniPair::kRapidityRange);
    cutYL->SetRangeD(-0.6,0.6);
    AliRsnCutSet* cutsL=new AliRsnCutSet("pairCutsLambdaStar",AliRsnTarget::kMother);
    cutsL->AddCut(cutMassL);
    cutsL->AddCut(cutYL);
    cutsL->SetCutScheme(TString::Format("%s&%s",cutMassL->GetName(),cutYL->GetName()).Data());
    
    i=0; // Lambda(1520)
    finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_LambdaStarp",task->GetName()));
    finder[i]->SetCutID(0,iCutP);
    finder[i]->SetDaughter(0,AliRsnDaughter::kProton);
    finder[i]->SetCharge(0,'+');
    finder[i]->SetCutID(1,iCutK);
    finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
    finder[i]->SetCharge(1,'-');
    finder[i]->SetResonanceMass(1.5195);
    finder[i]->SetResonancePDG(3124);
    finder[i]->SetPairCuts(cutsL);
    iCutLambdaStar[i]=task->AddResonanceFinder(finder[i]);
    
    i=1; // anti-Lambda1520
    finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_LambdaStara",task->GetName()));
    finder[i]->SetCutID(0,iCutP);
    finder[i]->SetDaughter(0,AliRsnDaughter::kProton);
    finder[i]->SetCharge(0,'-');
    finder[i]->SetCutID(1,iCutK);
    finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
    finder[i]->SetCharge(1,'+');
    finder[i]->SetResonanceMass(1.5195);
    finder[i]->SetResonancePDG(-3124);
    finder[i]->SetPairCuts(cutsL);
    iCutLambdaStar[i]=task->AddResonanceFinder(finder[i]);
    
    i=2; // P+ K+
    finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_PpKp",task->GetName()));
    finder[i]->SetCutID(0,iCutP);
    finder[i]->SetDaughter(0,AliRsnDaughter::kProton);
    finder[i]->SetCharge(0,'+');
    finder[i]->SetCutID(1,iCutK);
    finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
    finder[i]->SetCharge(1,'+');
    finder[i]->SetResonanceMass(1.5195);
    finder[i]->SetResonancePDG(3124);
    finder[i]->SetPairCuts(cutsL);
    iCutLambdaStar[i]=task->AddResonanceFinder(finder[i]);
    
    i=3; //P- K-
    finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_KmPm",task->GetName()));
    finder[i]->SetCutID(0,iCutP);
    finder[i]->SetDaughter(0,AliRsnDaughter::kProton);
    finder[i]->SetCharge(0,'-');
    finder[i]->SetCutID(1,iCutK);
    finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
    finder[i]->SetCharge(1,'-');
    finder[i]->SetResonanceMass(1.5195);
    finder[i]->SetResonancePDG(-3124);
    finder[i]->SetPairCuts(cutsL);
    iCutLambdaStar[i]=task->AddResonanceFinder(finder[i]);
    
    // sidebands
    AliRsnCutMiniPair* cutMassSB=new AliRsnCutMiniPair("cutMassSB",AliRsnCutMiniPair::kMassRange);
    cutMassSB->SetRangeD(1.480,1.496);
    AliRsnCutSet* cutsSB=new AliRsnCutSet("pairCutsSB",AliRsnTarget::kMother);
    cutsSB->AddCut(cutMassSB);
    cutsSB->AddCut(cutYL);
    cutsSB->SetCutScheme(TString::Format("%s&%s",cutMassSB->GetName(),cutYL->GetName()).Data());
    
    i=4; // Lambda(1520) sideband
    finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBLambdaStar",task->GetName()));
    finder[i]->SetCutID(0,iCutP);
    finder[i]->SetDaughter(0,AliRsnDaughter::kProton);
    finder[i]->SetCharge(0,'+');
    finder[i]->SetCutID(1,iCutK);
    finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
    finder[i]->SetCharge(1,'-');
    finder[i]->SetResonanceMass(1.5195);
    finder[i]->SetResonancePDG(3124);
    finder[i]->SetPairCuts(cutsSB);
    iCutLambdaStar[i]=task->AddResonanceFinder(finder[i]);
    
    i=5; // anti-Lambda(1520) sideband
    finder[i]=new AliRsnMiniResonanceFinder(Form("%s_ResonanceFinder_SBAntiLambdaStar",task->GetName()));
    finder[i]->SetCutID(0,iCutP);
    finder[i]->SetDaughter(0,AliRsnDaughter::kProton);
    finder[i]->SetCharge(0,'-');
    finder[i]->SetCutID(1,iCutK);
    finder[i]->SetDaughter(1,AliRsnDaughter::kKaon);
    finder[i]->SetCharge(1,'+');
    finder[i]->SetResonanceMass(1.5195);
    finder[i]->SetResonancePDG(-3124);
    finder[i]->SetPairCuts(cutsSB);
    iCutLambdaStar[i]=task->AddResonanceFinder(finder[i]);
    
    
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
    
    Int_t k,xID,cut1,cut2,pairID,ipdg;
    AliRsnDaughter::ESpecies d2=AliRsnDaughter::kUnknown;
    TString name,comp;
    Char_t charge1,charge2;
    AliRsnMiniOutput* out;
    
    for(i=0;i<10;i++){
        if(!i){
            xID=imID;
            name.Form("K0Lambdastarp");
            comp.Form("PAIR");
            cut1=iCutK0s;
            cut2=0;
            pairID=0;
            ipdg=3224;
        }else if(i==1){
            xID=imID;
            name.Form("K0Lambdastara");
            comp.Form("PAIR");
            cut1=iCutK0s;
            cut2=1;
            pairID=0;
            ipdg=-3224;
        }else if(i==2){
            xID=imID;
            name.Form("K0KpPp");
            comp.Form("PAIR");
            cut1=iCutK0s;
            cut2=2;
            pairID=0;
            ipdg=3224;
        }else if(i==3){
            xID=imID;
            name.Form("K0KmPm");
            comp.Form("PAIR");
            cut1=iCutK0s;
            cut2=3;
            pairID=0;
            ipdg=-3224;
        }else if(i==4){
            xID=imID;
            name.Form("K0LambdastarpSB");
            comp.Form("PAIR");
            cut1=iCutK0s;
            cut2=4;
            pairID=0;
            ipdg=3224;
        }else if(i==5){
            xID=imID;
            name.Form("K0LambdastaraSB");
            comp.Form("PAIR");
            cut1=iCutK0s;
            cut2=5;
            pairID=0;
            ipdg=-3224;
        }else if(i==6){
            xID=imID;
            name.Form("K0LambdastarpROTATE");
            if(!pairRotate) {comp.Form("ROTATE1");}
            else {comp.Form("ROTATE2");}
            cut1=iCutK0s;
            cut2=0;
            pairID=0;
            ipdg=3224;
        }else if(i==7){
            xID=imID;
            name.Form("K0LambdastaraROTATE");
            if(!pairRotate) {comp.Form("ROTATE1");}
            else {comp.Form("ROTATE2");}
            cut1=iCutK0s;
            cut2=1;
            pairID=0;
            ipdg=-3224;
        }else if(i==8){
            xID=imID;
            name.Form("K0LambdastarpMIX");
            comp.Form("MIX");
            cut1=iCutK0s;
            cut2=0;
            pairID=1;
            ipdg=3224;
        }else if(i==9){
            xID=imID;
            name.Form("K0LambdastaraMIX");
            comp.Form("MIX");
            cut1=iCutK0s;
            cut2=1;
            pairID=1;
            ipdg=-3224;
        }
        
        
        
        out=task->CreateOutput(Form("k0Lambdastar_%s%s",name.Data(),suffix),"HIST",comp.Data());
        out->SetDaughter(0,AliRsnDaughter::kKaon0);
        out->SetCutID(0,iCutK0s);
        out->SetCharge(0,0);
        
        out->SetDaughter(1,AliRsnDaughter::kLambdastar);
        out->SetCutID(1,iCutLambdaStar[cut2]);
        out->SetCharge(1,0);
        if(cut2!=4 && cut2!=5) out->SetUseStoredMass(1);
        
        if(!pairID) out->SetPairCuts(cutsPairSame);
        else out->SetPairCuts(cutsPairMix);
        out->SetMotherPDG(3224);
        out->SetMotherMass(mass);
        
        out->AddAxis(imID,240,1.8,3);// axis X: invmass or resolution
        //out->AddAxis(resID,200,-0.02,0.02);
        out->AddAxis(ptID,50,0.0,20.0);// axis Y: transverse momentum
        out->AddAxis(centID,nmult,multbins);// axis Z: centrality-multiplicity
    }
    
    // fill monitoring histogram for the resonance (Lambda(1520))
    for(i=0;i<2;i++) for(j=0;j<4;j++){
        if(!i) name.Form("Lambdastarp");
        else if(i==1) name.Form("Lambdastara");
        
        
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
        
        out=task->CreateOutput(Form("kxLambdastar_%s",name.Data()),"HIST",comp.Data());
        out->SetMotherPDG(finder[k]->GetResonancePDG());
        
        out->SetDaughter(0,finder[k]->GetDaughter(0));
        out->SetCutID(0,finder[k]->GetCutID(0));
        out->SetCharge(0,finder[k]->GetCharge(0));
        
        out->SetDaughter(1,finder[k]->GetDaughter(1));
        out->SetCutID(1,finder[k]->GetCutID(1));
        out->SetCharge(1,finder[k]->GetCharge(1));
        
        out->SetMotherMass(finder[k]->GetResonanceMass());
        if(j!=1) out->SetPairCuts(cutsL);
        else out->SetPairCuts(cutsSB);
        
        out->AddAxis(imID,150,1.3,1.6);
        out->AddAxis(ptID,50,0.0,20.0);
        out->AddAxis(centID,nmult,multbins);
    }
    
    return kTRUE;
}
