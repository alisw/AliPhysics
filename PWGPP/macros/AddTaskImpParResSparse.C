#if !(defined(__CINT__) || defined(__MAKECINT__))

#include "AliAnalysisTaskSEImpParResSparse.h"
#include "AliAnalysisManager.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisDataContainer.h"
#include "TChain.h"


#endif

AliAnalysisTaskSEImpParResSparse *AddTaskImpParResSparse(Bool_t readMC=kFALSE,
                                                         Int_t selPdg=-1,
                                                         Bool_t diamond=kTRUE,
                                                         Bool_t skipTrack=kTRUE,
                                                         Int_t minmult=0,
                                                         Int_t maxmult=1000000,
                                                         Int_t checkSDDIsIn=1,
                                                         TString dirName="",
                                                         Bool_t isAOD=kFALSE,
                                                         Int_t tracktype=0,
                                                         Int_t SPDreq=7,
                                                         Bool_t useCutGeoNcrNcl=kFALSE,
                                                         Bool_t fillexpertSparse=kFALSE,
                                                         Bool_t usetriggersel=kFALSE,
                                                         Int_t species=-1,
                                                         Bool_t finebinphi = kFALSE,
                                                         Bool_t checkSPDmod = kFALSE
                                                         ,Bool_t switchImpParrphiSparsePtBchargePhi = kTRUE
                                                         ,Bool_t switchImpParrphiSparsePtzVtxEtaPhi = kTRUE
                                                         ,Bool_t switchImpParrphiSparsePtEtaPhi = kTRUE 
                                                         ,Bool_t switchImpParPullrphiSparsePtEtaPhi = kTRUE
                                                         ,Bool_t switchImpParPullrphiSparsePtBchargePhi = kTRUE
                                                         ,Bool_t switchImpParzSparsePtBchargePhi = kTRUE
                                                         ,Bool_t switchImpParzSparsePtzVtxEtaPhi = kTRUE
                                                         ,Bool_t switchImpParzSparsePtEtaPhi = kTRUE
                                                         ,Bool_t switchImpParPullzSparsePtEtaPhi = kTRUE
                                                         ,Bool_t switchImpParPullzSparsePtBchargePhi = kTRUE
                                                         ,Bool_t useITSrefitTrks_highMev = kFALSE
                                          				 )

{
    
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskImpParRes", "No analysis manager to connect to.");
        return NULL;
    }
    
    
    //
    
    
    // Aanalysis task
    AliAnalysisTaskSEImpParResSparse *d0ResTask = new AliAnalysisTaskSEImpParResSparse("d0ResAnalysis");
    d0ResTask->SetDebugLevel(2);
    d0ResTask->SetReadMC(readMC);
    d0ResTask->SetIsAOD(isAOD);
    d0ResTask->SetSelectedPdg(selPdg);
    d0ResTask->SetUseDiamond(diamond);
    d0ResTask->SetSkipTrack(skipTrack);
    d0ResTask->SetMultiplicityRange(minmult,maxmult);
    //d0ResTask->SetFillSparseForExpert(fillexpertSparse);
    d0ResTask->SetUsePtWeights(0,1.);
    d0ResTask->SetCheckSDDIsIn(checkSDDIsIn);
    d0ResTask->SetParticleSpecies(species);
    if(usetriggersel){
        d0ResTask->SetUseTriggerSelection(usetriggersel);
        d0ResTask->SetTriggerClass("");
        d0ResTask->SetTriggerMask(AliVEvent::kINT7);
        //d0ResTask->SetTriggerMask(0);
    }
    d0ResTask->SetTrackType(tracktype);
    // mfaggin
    d0ResTask->SetUseFinerPhiBins(finebinphi);
    d0ResTask->SetStoreSPDmodulesInfo(checkSPDmod);
    d0ResTask->SetFillSparse_ImpParrphiSparsePtBchargePhi(switchImpParrphiSparsePtBchargePhi);
    d0ResTask->SetFillSparse_ImpParrphiSparsePtzVtxEtaPhi(switchImpParrphiSparsePtzVtxEtaPhi);
    d0ResTask->SetFillSparse_ImpParrphiSparsePtEtaPhi(switchImpParrphiSparsePtEtaPhi);
    d0ResTask->SetFillSparse_ImpParPullrphiSparsePtEtaPhi(switchImpParPullrphiSparsePtEtaPhi);
    d0ResTask->SetFillSparse_ImpParPullrphiSparsePtBchargePhi(switchImpParPullrphiSparsePtBchargePhi);
    d0ResTask->SetFillSparse_ImpParzSparsePtBchargePhi(switchImpParzSparsePtBchargePhi);
    d0ResTask->SetFillSparse_ImpParzSparsePtzVtxEtaPhi(switchImpParzSparsePtzVtxEtaPhi);
    d0ResTask->SetFillSparse_ImpParzSparsePtEtaPhi(switchImpParzSparsePtEtaPhi);
    d0ResTask->SetFillSparse_ImpParPullzSparsePtEtaPhi(switchImpParPullzSparsePtEtaPhi);
    d0ResTask->SetFillSparse_ImpParPullzSparsePtBchargePhi(switchImpParPullzSparsePtBchargePhi);

    d0ResTask->SetFillSparseForExpert(fillexpertSparse);

    // mfaggin
    d0ResTask->SetEstimatorHighMultEv(useITSrefitTrks_highMev);
    
    AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("d0ResAnalysisESDTrackCuts");
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    switch(SPDreq){
        case(1):
            esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOnlyFirst);
            break;
        case(2):
            esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOnlySecond);
            break;
        case(3):
            esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
            break;
        case(4):
            esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
            break;
        case(5):
            esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
            break;
        case(6):
            esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kNone);
            break;
        default:
            esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
            if(checkSDDIsIn==1)esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSDD,AliESDtrackCuts::kBoth);
            esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSSD,AliESDtrackCuts::kBoth);
            break;
    }
    esdTrackCuts->SetMaxChi2PerClusterTPC(4);
    esdTrackCuts->SetMinNClustersTPC(70);
    if(useCutGeoNcrNcl && !isAOD) esdTrackCuts->SetCutGeoNcrNcl(3,130,1.5,0.85,0.7);
    d0ResTask->SetESDtrackCuts(esdTrackCuts);
    
    
    if(useCutGeoNcrNcl && isAOD){
        d0ResTask->SetUseCutGeoNcrNcl(useCutGeoNcrNcl);
        d0ResTask->ConfigureCutGeoNcrNcl(3,130,1.5,0.85,0.7);
    }
    
    mgr->AddTask(d0ResTask);
    
    TString fname=Form("%s:ImpParRes_Performance",mgr->GetCommonFileName());
    if(selPdg>0) {fname+=selPdg;}
    fname += dirName.Data();
    
    TString name=dirName;
    //
    // Create containers for input/output
    AliAnalysisDataContainer *cinputd0distr = mgr->CreateContainer(Form("cinputd0distr%s",name.Data()),TChain::Class(),
                                                                   AliAnalysisManager::kInputContainer);
    
    AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("coutputd0_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             fname.Data());
    
    
    // Attach input
    mgr->ConnectInput(d0ResTask,0,mgr->GetCommonInputContainer());
    // Attack output
    mgr->ConnectOutput(d0ResTask,1,coutput);
    
    return d0ResTask;
}
