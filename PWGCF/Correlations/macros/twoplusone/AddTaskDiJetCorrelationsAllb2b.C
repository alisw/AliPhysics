
// $Id: AddTaskDiJetCorrelationsAllb2b.C

AliAnalysisTaskDiJetCorrelationsAllb2b *AddTaskDiJetCorrelationsAllb2b(TString suffixName="",
                                                                               Bool_t ppOrPbPb = kTRUE,
                                                                               Bool_t SEorME = kTRUE,
                                                                               Bool_t twoplus1 = kTRUE,
                                                                               Bool_t bkgSE = kFALSE,
                                                                               Double_t pTrg1min = 12.0,
                                                                               Double_t pTrg1max = 16.0,
                                                                               Double_t pTrg2min = 5.0,
                                                                               Double_t pTrg2max = 8.0,
                                                                               Bool_t T1T2Equal = kFALSE,
                                                                           
                                         Bool_t fineBinsME = kTRUE,
                                                                               Double_t bit = 272,
                                                                               TString fileTrackeff = "",
                                                                               Bool_t resCut = kFALSE,
                                                                               Bool_t conversionCut = kFALSE,
                                                                               Bool_t TTRcut = kFALSE
                                                                               )
{
    
    Bool_t UseFbits = kTRUE;
    
    
    
    
    
    
    //____________________________________| Correlation class setting..
    AliAnalysisTaskDiJetCorrelationsAllb2b *dijetcorrelations = new AliAnalysisTaskDiJetCorrelationsAllb2b("");
    dijetcorrelations->SelectCollisionCandidates(AliVEvent::kMB);
    dijetcorrelations->SetCorr2plus1or1plus1(twoplus1);
    dijetcorrelations->SetSystem(ppOrPbPb); //PbPb = kTRUE
    dijetcorrelations->SetSEorME(SEorME); //kTRUE for mixed events
    if(SEorME)dijetcorrelations->SetMESettings(1000, 50000, 5); //evt,track,minMixEvents
    //if(SEorME)dijetcorrelations->SetMESettings(2, 10, 2); //evt,track,minMixEvents
    dijetcorrelations->SetTrigger1PTValue(pTrg1min, pTrg1max); //GeV/c
    dijetcorrelations->SetTrigger2PTValue(pTrg2min, pTrg2max); //GeV/c
    dijetcorrelations->SetFilterBit(UseFbits);
    if(UseFbits)dijetcorrelations->SetFilterType(bit);
    if(ppOrPbPb)dijetcorrelations->SetCentralityRange(0., 100); // 0-100%
    dijetcorrelations->SetDataType(kTRUE); //track Data/MC tracks=1 or MC Part=0?
    dijetcorrelations->SetVarCentBin(kTRUE);
    dijetcorrelations->SetVarPtBin(kTRUE);
    dijetcorrelations->SetResonanceCut(resCut);
    dijetcorrelations->SetConversionCut(conversionCut);
    dijetcorrelations->SetTwoTrackEfficiencyCut(TTRcut);
    dijetcorrelations->SetEqualT1T2Demand(T1T2Equal);
    dijetcorrelations->SetFinerBinsME(fineBinsME);
    
   // dijetcorrelations->SetDiJetAlphaAngle(TMath::Pi())/8);
    dijetcorrelations->SetDiJetAlphaAngle(0.392699);
    //if(effLoc!="")dijetcorrelations->SetEffCorrection(GetEfficiencyCorr(effLoc));
    
    dijetcorrelations->SetBkgSE(bkgSE);
    dijetcorrelations->SetBkgSEBothSide(kTRUE);
    
    
    
    
    if(fileTrackeff!="")
        
    {
        
        TFile* fAssoTracksEffMap;
        THnF* hEff;
        fAssoTracksEffMap=TFile::Open(fileTrackeff.Data());
        
        if(!fAssoTracksEffMap ||(fAssoTracksEffMap&& !fAssoTracksEffMap->IsOpen())){
            
            Printf("Input Associated Track Efficeincy Map object not found");
            
            return;
            
        }
        
        THnF *hEff = (THnF*)fAssoTracksEffMap->Get("correction");
        if(!hEff){
            Printf("%s%d Couldn't find correction",(char*)__FILE__,__LINE__);
            return;
        }
        
        
        
        dijetcorrelations->SetEfficiencyWeightMap(hEff);
        
    }
    
    
    
    // Create containers for input/output
    TString finDirname         = "_DiJet";
    TString inname             = "cinputDiJetCorrelations";
    TString outBasicname       = "coutputDiJetBase";
    TString outCorrname        = "coutputDiJetCorrHistos";
    
    finDirname += suffixName.Data();
    inname            +=   finDirname.Data();
    outBasicname      +=   finDirname.Data();
    outCorrname       +=   finDirname.Data();
    
    
    // Get the pointer to the existing analysis manager via the static access method.
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        cout<<"AddTaskDiJetCorrelationsAllb2b", "No analysis manager to connect to."<<endl;
    }
    
    mgr->AddTask(dijetcorrelations);
    
    //Input and Output Slots:
    AliAnalysisDataContainer *cinputDiJetCorrelations = mgr->CreateContainer(inname,TChain::Class(), AliAnalysisManager::kInputContainer);
    //TString outputfile = AliAnalysisManager::GetCommonFileName();
    //outputfile += ":PWGCF_Di_Jet_Corr";
    TString outputfile = AliAnalysisManager::GetCommonFileName();//"resultsDiJetCorrelationsT112to16T25to8Dec8.root";
    
    
    AliAnalysisDataContainer *coutputDiJetCorrelations1 = mgr->CreateContainer(outBasicname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutputDiJetCorrelations2 = mgr->CreateContainer(outCorrname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    
    
    mgr->ConnectInput(dijetcorrelations,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(dijetcorrelations,1,coutputDiJetCorrelations1);
    mgr->ConnectOutput(dijetcorrelations,2,coutputDiJetCorrelations2);
    
    return dijetcorrelations;
    
}
