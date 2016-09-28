

/* $Id: AddTaskDplusCorrelations.C 58712 2012-09-20 08:38:36Z prino $ */
//AddTask for the Dplus - Hadron (or Kaon/K0) Corelation with same/mixed event
//Jitendra Kumar (Last updated on 31.01.2016)

AliAnalysisTaskSEDplusCorrelations *AddTaskDplusCorrelations(TString suffix="",
                                                             Int_t  fOption = 1,
                                                             Bool_t fMixing = kFALSE,
                                                             Bool_t readMC  = kFALSE,
                                                             Bool_t genMC   = kFALSE,
                                                             Bool_t tracks  = kTRUE,
                                                             TString fileDplusCuts="",
                                                             TString DplusCutsObjName="",
                                                             TString fileTrackCuts="",
                                                             TString TrackCutsObjName="",
                                                             Bool_t isTrackEff = kFALSE,
                                                             TString fileTrackeff="",
                                                             Bool_t isDplusEff = kFALSE,
                                                             TString fileDplusEff="",
                                                             Bool_t PoolbyPool=kFALSE,
                                                             Bool_t useCentrality = kFALSE,
			  				     Int_t AODprot=1)
{
    
    
    const Int_t centralityEstimator = 7; // enum from AliRDHFCuts.h
    Double_t etacorr  =  0.9;
    
    //1. Dplus Cuts either from PPStd OR from cut file
    Bool_t PPstdcuts=kFALSE;
    TFile* fDplusCuts;
    if(fileDplusCuts.EqualTo("")){
        PPstdcuts=kTRUE;
    } else {
        fDplusCuts=TFile::Open(fileDplusCuts.Data());
        if(!fDplusCuts ||(fDplusCuts && !fDplusCuts->IsOpen())){
            AliFatal("Input Dplus cuts file(object) is not found: EXIT");
            return;
        }
    }
    AliRDHFCutsDplustoKpipi* DplusCorrCuts=new AliRDHFCutsDplustoKpipi();
    if(PPstdcuts)  DplusCorrCuts->SetStandardCutsPP2010();
    else DplusCorrCuts = (AliRDHFCutsDplustoKpipi*)fDplusCuts->Get(DplusCutsObjName.Data());
    if(!DplusCorrCuts){
        Printf("Dplus cut objects not found");
        return;
    }
    //DplusCorrCuts->PrintAll();
    
    
    
    //2.Track cuts AliHFAssociated from cut file
    TFile* fAssoTrackCuts;
    if(fileTrackCuts.EqualTo("") ) {
        AliFatal("Input associated track cuts file not loaded");
        return;
    } else {
        fAssoTrackCuts=TFile::Open(fileTrackCuts.Data());
        if(!fAssoTrackCuts ||(fAssoTrackCuts&& !fAssoTrackCuts->IsOpen())){
            AliFatal("Input associated cuts file object are not found");
            return;
        }
    }
    AliHFAssociatedTrackCuts* HFAssoTrackCuts=new AliHFAssociatedTrackCuts();
    HFAssoTrackCuts = (AliHFAssociatedTrackCuts*)fAssoTrackCuts->Get(TrackCutsObjName.Data());
    if(!HFAssoTrackCuts){
        AliFatal("Specific AliHFAssociatedTrackCuts not found");
        return;
    }
    HFAssoTrackCuts->SetTrackCutsNames();
    HFAssoTrackCuts->SetvZeroCutsNames();
    
    
    
    //3. Track efficeincy map
    TFile* fAssoTracksEffMap;
    if(isTrackEff){
        if(fileTrackeff.EqualTo("") ) {
            AliFatal("Input Associated Track Efficeincy Map not loaded ");
            return;
        } else {
            fAssoTracksEffMap=TFile::Open(fileTrackeff.Data());
            if(!fAssoTracksEffMap ||(fAssoTracksEffMap&& !fAssoTracksEffMap->IsOpen())){
                AliFatal("Input Associated Track Efficeincy Map object not found");
                return;
            }
        }
        TCanvas *cTrackEffMap = (TCanvas*)fAssoTracksEffMap->Get("c");
        TH3D *hEffTrackMap = (TH3D*)cTrackEffMap->FindObject("heff_rebin");
        HFAssoTrackCuts->SetEfficiencyWeightMap(hEffTrackMap);
    }
    
    
    
    //4. Dplus efficeincy map
    TFile* fDplusEffMap;
    if(isDplusEff){
        if(fileDplusEff.EqualTo("") ){
            AliFatal("Input Dplus Efficeincy Map not loaded");
            return;
        }else{
            fDplusEffMap=TFile::Open(fileDplusEff.Data());
            if(!fDplusEffMap ||(fDplusEffMap&& !fDplusEffMap->IsOpen())){
                AliFatal("Input Dplus Track Efficeincy Map object not found");
                return;
            }
        }
        TCanvas *cDplusEffMap = (TCanvas*)fDplusEffMap->Get("c1");
        TH2D *hEffDplusMap_Frmc = (TH2D*)cDplusEffMap->FindObject("heff2D_rebin");
        HFAssoTrackCuts->SetTriggerEffWeightMap(hEffDplusMap_Frmc);
    }
    
    //Priting all cuts
    HFAssoTrackCuts->PrintAll();
    
    
    
    
    //6. Correlation Setter
    AliAnalysisTaskSEDplusCorrelations *dpluscorrTask = new AliAnalysisTaskSEDplusCorrelations("DplusCorrelation",DplusCorrCuts,HFAssoTrackCuts);
    dpluscorrTask->SetEventMixing(fMixing);
    dpluscorrTask->SetCorrelator(fOption);
    dpluscorrTask->SetDebugLevel(0);
    dpluscorrTask->SetEtaRagne(etacorr);
    dpluscorrTask->SetCorrFormPart(genMC);
    dpluscorrTask->SetCorrFormTrack(tracks);
    dpluscorrTask->SetDataOrMC(readMC);
    dpluscorrTask->SetUseBit(kTRUE);
    dpluscorrTask->SetTCConfig(kTRUE);
    dpluscorrTask->SetTrackEffActive(isTrackEff);
    dpluscorrTask->SetDplusEffActive(isDplusEff);
    dpluscorrTask->SetSystem(useCentrality); //TRUE means pbpb Or pA
    dpluscorrTask->SetPoolByPoolCorr(PoolbyPool); //TRUE means pbpb Or pA
    if(useCentrality)dpluscorrTask->SetUseCentrality(useCentrality, centralityEstimator);
    dpluscorrTask->SetCheckCutDist(kTRUE);
    dpluscorrTask->SetAODMismatchProtection(AODprot);    
    
    //7. Create container for input/output
    TString finDirname = "";
    if(PPstdcuts)finDirname+="PP_Dplus";
    else if(!PPstdcuts)finDirname+="pPb_Dplus";
    
    if(fOption==1)finDirname+="HadCorr";
    else if(fOption==2)finDirname+="KaonCorr";
    else if(fOption==3)finDirname+="KZeroCorr";
    
    if(fMixing)finDirname+="_ME";
    else if(!fMixing)finDirname+="_SE";
    
    if(isTrackEff)finDirname+="_wTrkEff";
    else if(!isTrackEff)finDirname+="_woTrkEff";
    
    if(isDplusEff)finDirname+="_wDkEff_";
    else if(!isDplusEff)finDirname+="_woDkEff_";
    finDirname += suffix.Data();
    
    TString inname             = "cin_";
    TString outBasicname       = "coutBasicPlots_";
    TString outCorrname        = "coutHistos_";
    TString outDcutsname       = "coutDplusCuts_";
    TString outTrcutsname      = "coutTrackCuts_";
    TString outNormname        = "coutNorm_";
    inname            +=   finDirname.Data();
    outBasicname      +=   finDirname.Data();
    outCorrname       +=   finDirname.Data();
    outDcutsname      +=   finDirname.Data();
    outTrcutsname     +=   finDirname.Data();
    outNormname       +=   finDirname.Data();
    
    
    
    //8. Get the pointer to the existing analysis manager via the static access method.
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        AliFatal("AddTaskDplusCorrelations: No analysis manager to connect to");
        return;
    }
    mgr->AddTask(dpluscorrTask);
    
    //Input and Output Slots:
    AliAnalysisDataContainer *cinputDplusCorrelations = mgr->CreateContainer(inname,TChain::Class(), AliAnalysisManager::kInputContainer);
    TString outputfile = AliAnalysisManager::GetCommonFileName();
    outputfile += Form(":%s",finDirname.Data());
    
    AliAnalysisDataContainer *coutputCorrDplus1 = mgr->CreateContainer(outBasicname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutputCorrDplus2 = mgr->CreateContainer(outCorrname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutputCorrDplus3 = mgr->CreateContainer(outDcutsname,AliRDHFCutsDplustoKpipi::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutputCorrDplus4 = mgr->CreateContainer(outTrcutsname,AliHFAssociatedTrackCuts::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());
    AliAnalysisDataContainer *coutputCorrDplus5 = mgr->CreateContainer(outNormname,AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());
    mgr->ConnectInput(dpluscorrTask,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(dpluscorrTask,1,coutputCorrDplus1);
    mgr->ConnectOutput(dpluscorrTask,2,coutputCorrDplus2);
    mgr->ConnectOutput(dpluscorrTask,3,coutputCorrDplus3);  
    mgr->ConnectOutput(dpluscorrTask,4,coutputCorrDplus4);
    mgr->ConnectOutput(dpluscorrTask,5,coutputCorrDplus5);
    
    
    return dpluscorrTask;
    
}


//EOF
