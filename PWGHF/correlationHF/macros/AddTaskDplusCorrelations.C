

/* $Id: AddTaskDplusCorrelations.C 58712 2012-09-20 08:38:36Z prino $ */
//AddTask for the Dplus - Hadron (or Kaon/K0) Corelation with same/mixed event
//Jitendra Kumar(Last updated on 02.10.2016) //AOD production setting:
AliAnalysisTaskSEDplusCorrelations *AddTaskDplusCorrelations(TString suffix="",
                                                             Bool_t useCutFileSBRanges=kFALSE,
                                                             TString fileSignalSBRange="",
                                                             TString  fSys = "PP",
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
                                                             Int_t centralityEstimator=AliRDHFCuts::kCentZNA,
                                                             Int_t AODproduction=1,
                                                             Bool_t IncDCutQA=kTRUE,
                                                             Bool_t IncDCutQABefore=kFALSE,
                                                             Bool_t IsfillTrees=kTRUE,
                                                             Double_t fractAccME=100.,
                                                             Double_t DPtThrs=3,
                                                             Bool_t LeadPartCorr=kFALSE,
                                                             Double_t massWidth=0.004,
                                                             Bool_t PoolCent = kFALSE,
                                                             Bool_t autosignalSBrange)
{
    
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
    dpluscorrTask->SetBinWidth(massWidth);
    dpluscorrTask->SetDebugLevel(0);
    //dpluscorrTask->SetEtaRagne(etacorr);
    dpluscorrTask->SetCorrFormPart(genMC);
    dpluscorrTask->SetAutoSignalSBRange(autosignalSBrange);
    dpluscorrTask->SetCorrFormTrack(tracks);
    dpluscorrTask->SetLeadPartCorrelation(LeadPartCorr);
    dpluscorrTask->SetDataOrMC(readMC);
    dpluscorrTask->SetUseBit(kTRUE);
    dpluscorrTask->SetTCConfig(kTRUE);
    dpluscorrTask->SetTrackEffActive(isTrackEff);
    dpluscorrTask->SetDplusEffActive(isDplusEff);
    dpluscorrTask->SetSystem(useCentrality); //TRUE means pbpb Or pA
    dpluscorrTask->SetPoolByPoolCorr(PoolbyPool); //TRUE means pbpb Or pA
    if(useCentrality)dpluscorrTask->SetUseCentrality(useCentrality, centralityEstimator, PoolCent);
    dpluscorrTask->SetCheckCutDistandChoice(IncDCutQA, IncDCutQABefore);
    dpluscorrTask->SetAODMismatchProtection(AODproduction);
    dpluscorrTask->SetMinDPt(DPtThrs);
    dpluscorrTask->SetFillTrees(IsfillTrees,fractAccME);
    
    //7. Create container for input/output
    TString finDirname = "";
    
    if(useCutFileSBRanges) { //use SB ranges from cut file
        
        filerange=TFile::Open(fileSignalSBRange.Data());
        if(!filerange && !filerange->IsOpen()){
            AliFatal("Signal and Side Band Ranges are  not found: EXIT");
            return;
        }
        
        TVectorD *LSBLow = (TVectorD*)filerange->Get("vLSBLow");
        TVectorD *LSBUpp = (TVectorD*)filerange->Get("vLSBUpp");
        TVectorD *RSBLow = (TVectorD*)filerange->Get("vRSBLow");
        TVectorD *RSBUpp = (TVectorD*)filerange->Get("vRSBUpp");
        TVectorD *SandBLow = (TVectorD*)filerange->Get("vSandBLow");
        TVectorD *SandBUpp = (TVectorD*)filerange->Get("vSandBUpp");
        
        if(!LSBLow||!LSBUpp||!RSBLow||!RSBUpp ||!SandBLow ||!SandBUpp) {printf("Error! No SB ranges found in the Associated track cut file, but useCutFileSBRanges==kTRUE! Exiting...\n"); return;}
        
        dpluscorrTask->SetLSBLowerUpperLim(LSBLow->GetMatrixArray(), LSBUpp->GetMatrixArray());
        dpluscorrTask->SetSandBLowerUpperLim(SandBLow->GetMatrixArray(), SandBUpp->GetMatrixArray());
        dpluscorrTask->SetRSBLowerUpperLim(RSBLow->GetMatrixArray(), RSBUpp->GetMatrixArray());
        
    }
    if (!useCutFileSBRanges){
        if(fSys=="pp" || fSys=="PP"|| fSys=="p-p"){
            
            //Setting up the mass ranges for LSB,S+B and RSB region (Check values from fits)
            Double_t LSBLowLim[15] = {0.,0.,0.,1.7930,1.7770,1.7930,1.7530,1.7690,1.7690,1.7690,1.7690,1.7690,1.7690,0.,0.};
            Double_t LSBUppLim[15] = {0.,0.,0.,1.8410,1.8330,1.8330,1.8170,1.8250,1.8250,1.8250,1.8250,1.8250,1.8250,0.,0.};
            Double_t RSBLowLim[15] = {0.,0.,0.,1.8970,1.9050,1.9050,1.9210,1.9130,1.9130,1.9130,1.9130,1.9130,1.9130,0.,0.};
            Double_t RSBUppLim[15] = {0.,0.,0.,1.9450,1.9610,1.9450,1.9770,1.9690,1.9690,1.9690,1.9690,1.9690,1.9690,0.,0.};
            Double_t SandBLowLim[15] = {0.,0.,0.,1.8490,1.8410,1.8490,1.8410,1.8410,1.8410,1.8410,1.8410,1.8410,1.8410,0.,0.};
            Double_t SandBUppLim[15] = {0.,0.,0.,1.8890,1.8970,1.8890,1.8970,1.8970,1.8970,1.8970,1.8970,1.8970,1.8970,0.,0.};
            
            
        }else if(fSys=="pPb" || fSys=="p-Pb" || fSys=="PPb"){
            
            //Setting up the mass ranges for LSB,S+B and RSB region (Check values from fits)
            Double_t LSBLowLim[15]   = {0.,0.,0.,1.7930,1.7930,1.7850,1.7850,1.7770,1.7770,1.7690,1.7690,1.7530,1.7610,0.,0.};
            Double_t LSBUppLim[15]   = {0.,0.,0.,1.8410,1.8410,1.8330,1.8330,1.8330,1.8330,1.8250,1.8250,1.8170,1.8250,0.,0.};
            Double_t SandBLowLim[15] = {0.,0.,0.,1.8490,1.8490,1.8490,1.8490,1.8490,1.8410,1.8410,1.8410,1.8410,1.8410,0.,0.};
            Double_t SandBUppLim[15] = {0.,0.,0.,1.8970,1.8970,1.8970,1.8970,1.8970,1.8970,1.8970,1.8970,1.9050,1.8970,0.,0.};
            Double_t RSBLowLim[15]   = {0.,0.,0.,1.9050,1.9050,1.9050,1.9050,1.9130,1.9130,1.9130,1.9130,1.9210,1.9130,0.,0.};
            Double_t RSBUppLim[15]   = {0.,0.,0.,1.9530,1.9530,1.9610,1.9530,1.9610,1.9610,1.9690,1.9770,1.9850,1.9770,0.,0.};
            
            
        }else if(fSys=="PbPb" || fSys=="Pb-Pb" || fSys=="PbPb"){
            
            //Setting up the mass ranges for LSB,S+B and RSB region (Check values from fits)
            Double_t LSBLowLim[15]   = {0.,0.,0.,1.797,1.775,1.791,1.759,1.767,1.767,1.767,1.767,1.767,1.767,1.767,1.767};
            Double_t LSBUppLim[15]   = {0.,0.,0.,1.837,1.829,1.837,1.821,1.821,1.821,1.821,1.821,1.821,1.821,1.821,1.821};
            Double_t SandBLowLim[15] = {0.,0.,0.,1.843,1.843,1.843,1.843,1.843,1.843,1.843,1.843,1.843,1.843,1.843,1.843};
            Double_t SandBUppLim[15] = {0.,0.,0.,1.885,1.885,1.885,1.885,1.885,1.885,1.885,1.885,1.885,1.885,1.885,1.885};
            Double_t RSBLowLim[15]   = {0.,0.,0.,1.901,1.907,1.901,1.917,1.917,1.917,1.917,1.917,1.917,1.917,1.917,1.917};
            Double_t RSBUppLim[15]   = {0.,0.,0.,1.939,1.963,1.947,1.979,1.971,1.971,1.971,1.971,1.971,1.971,1.971,1.971};
            
            
        }else{
            cout << "EROOR in AddTask: Please check system name" << endl;
            return;
        }
        dpluscorrTask->SetLSBLowerUpperLim(LSBLowLim, LSBUppLim);
        dpluscorrTask->SetSandBLowerUpperLim(SandBLowLim, SandBUppLim);
        dpluscorrTask->SetRSBLowerUpperLim(RSBLowLim, RSBUppLim);
    }
    
    if(fSys=="pp" || fSys=="PP"|| fSys=="p-p") finDirname+="PP";
    else if(fSys=="pPb" || fSys=="p-Pb" || fSys=="PPb") finDirname+="pPb";
    else if(fSys=="PbPb" || fSys=="Pb-Pb" || fSys=="PbPb") finDirname+="PbPb";
    
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
    
    if(PPstdcuts)finDirname+="DefCuts";
    else if(!PPstdcuts)finDirname+="FileCuts";
    
    if(PoolbyPool)finDirname+="PlbyPl";
    else if(!PoolbyPool)finDirname+="IntPl";
    
    TString inname             = "cin_";
    TString outBasicname       = "coutBasicPlots_";
    TString outCorrname        = "coutHistos_";
    TString outDcutsname       = "coutDplusCuts_";
    TString outTrcutsname      = "coutTrackCuts_";
    TString outNormname        = "coutNorm_";
    TString Dplustree          = "TreeDplus_";
    TString Trackstree         = "TreeTracks_";
    
    inname            +=   finDirname.Data();
    outBasicname      +=   finDirname.Data();
    outCorrname       +=   finDirname.Data();
    outDcutsname      +=   finDirname.Data();
    outTrcutsname     +=   finDirname.Data();
    outNormname       +=   finDirname.Data();
    Dplustree         +=   finDirname.Data();
    Trackstree        +=   finDirname.Data();
    
    
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
    AliAnalysisDataContainer *coutputCorrDplus6 = mgr->CreateContainer(Dplustree,TTree::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());
    AliAnalysisDataContainer *coutputCorrDplus7 = mgr->CreateContainer(Trackstree,TTree::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());
    mgr->ConnectInput(dpluscorrTask,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(dpluscorrTask,1,coutputCorrDplus1);
    mgr->ConnectOutput(dpluscorrTask,2,coutputCorrDplus2);
    mgr->ConnectOutput(dpluscorrTask,3,coutputCorrDplus3);
    mgr->ConnectOutput(dpluscorrTask,4,coutputCorrDplus4);
    mgr->ConnectOutput(dpluscorrTask,5,coutputCorrDplus5);
    if(IsfillTrees) mgr->ConnectOutput(dpluscorrTask,6,coutputCorrDplus6);
    if(IsfillTrees) mgr->ConnectOutput(dpluscorrTask,7,coutputCorrDplus7);
    
    return dpluscorrTask;
}
//EOF
