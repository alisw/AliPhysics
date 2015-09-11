// $Id: AddTaskDiJetCorr1plus1Bkg.C

AliAnalysisTaskDiJetCorr1plus1Bkg *AddTaskDiJetCorr1plus1Bkg(TString suffixName = " ",
                                                             TString fileTrackeff = "",
                                                             Bool_t ppOrPbPb = kTRUE,
                                                             Int_t bit = 272,
                                                             Double_t pTrgmin = 12.0,
                                                             Double_t pTrgmax = 16.0,
                                                             Bool_t resCut = kFALSE,
                                                             Bool_t conversionCut = kFALSE,
                                                             Bool_t TTRcut = kFALSE
                                                             
                                                             )
{
    
    
    AliAnalysisTaskDiJetCorr1plus1Bkg *dijetcorr = new AliAnalysisTaskDiJetCorr1plus1Bkg("DiJetCorr");
    
    dijetcorr->SelectCollisionCandidates(AliVEvent::kCentral|AliVEvent::kSemiCentral|AliVEvent::kMB);
    
    if(fileTrackeff!="")
        
    {
        
        TFile* fAssoTracksEffMap;
        THnF* hEff;
        fAssoTracksEffMap=TFile::Open(fileTrackeff.Data());
        
        if(!fAssoTracksEffMap ||(fAssoTracksEffMap&& !fAssoTracksEffMap->IsOpen())){
            
            AliFatal("Input Associated Track Efficeincy Map object not found");
            
            return;
            
        }
        
        THnF *hEff = (THnF*)fAssoTracksEffMap->Get("correction");
        if(!hEff){
            Printf("%s%d Couldn't find correction",(char*)__FILE__,__LINE__);
            return;
        }
        
        
        // hEff = dynamic_cast<THnF*>tmp1->Clone("hEff");
        
        dijetcorr->SetEfficiencyWeightMap(hEff);
        
    }
    
    
    dijetcorr->SetSystemType(ppOrPbPb);
    dijetcorr->SetMixingTracks(50000);
    dijetcorr->SetFilterBit(bit);
    dijetcorr->SetTriggerpTValue(pTrgmin, pTrgmax); //GeV/c
    dijetcorr->SetResonanceCut(resCut);
    dijetcorr->SetConversionCut(conversionCut);
    dijetcorr->SetTwoTrackEfficiencyCut(TTRcut);
    
    
    
    
    // Create containers for input/output
    TString finDirname         = "_DiJetBkg";
    TString inname             = "cinputDiJetCorrelations";
    TString outCorrname        = "coutput1plus1";
    
    finDirname += suffixName.Data();
    inname            +=   finDirname.Data();
    outCorrname       +=   finDirname.Data();
    
    
    // Get the pointer to the existing analysis manager via the static access method.
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        cout<<"AddTaskDiJetCorr1plus1Bkg", "No analysis manager to connect to."<<endl;
    }
    
    mgr->AddTask(dijetcorr);
    
    //Input and Output Slots:
    AliAnalysisDataContainer *cinputDiJetCorrelations = mgr->CreateContainer(inname,TChain::Class(), AliAnalysisManager::kInputContainer);
    //TString outputfile = AliAnalysisManager::GetCommonFileName();
    //outputfile += ":PWGCF_Di_Jet_Corr";
    TString outputfile = AliAnalysisManager::GetCommonFileName();//"resultsDiJetCorrelationsT112to16T25to8Dec8.root";
    
    
    AliAnalysisDataContainer *coutputDiJetCorr = mgr->CreateContainer(outCorrname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    
    
    mgr->ConnectInput(dijetcorr,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(dijetcorr,1,coutputDiJetCorr);
    
    return dijetcorr;
    
}



