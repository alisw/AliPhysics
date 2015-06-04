
void AddTaskPIDconfig(Int_t CentralityTriggerSelection = AliVEvent::kMB, Int_t centralityMinlimit=0, Int_t centralityMaxlimit=5 ,Double_t FilterBit=1, Bool_t PIDcuts=kFALSE,Bool_t isLHC11h=kFALSE,TString useroutputfile="output"){
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskPID", "No analysis manager to connect to.");
        return 0x0;
    }
    
    // standard with task
    printf("========================================================================================\n");
    printf("PID: Initialising AliAnalysisTaskPIDconfig\n");
    printf("========================================================================================\n");
    
    Double_t centrMin[8] = {0,0,10,20,30,40,60,60};
    Double_t centrMax[8] = {1,2,20,30,40,50,70,80};
    
    const int ncentr = centralityMaxlimit - centralityMinlimit ;
    TString outputfile[ncentr];
    AliAnalysisDataContainer *coutput1[ncentr];
    AliAnalysisTaskPIDconfig *pidTask[ncentr];
    int icentr = 0;
    
    for(int i=0;i<ncentr;i++){
        icentr = i + centralityMinlimit;
        outputfile[i] = useroutputfile;
        outputfile[i].Append(".root");
        pidTask[i] = new AliAnalysisTaskPIDconfig(Form("pidTask_%.f-%.f",centrMin[icentr],centrMax[icentr]));
        pidTask[i]->SelectCollisionCandidates(CentralityTriggerSelection);
        pidTask[i]->SetCutTPCmultiplicityOutliersAOD(kTRUE);
        pidTask[i]->SetData2011(isLHC11h);
        pidTask[i]->SetFilterBit(FilterBit);
        pidTask[i]->SetUseCentrality(kTRUE);
        pidTask[i]->SetCentralityPercentileMin(centrMin[icentr]);
        pidTask[i]->SetCentralityPercentileMax(centrMax[icentr]);
        pidTask[i]->SetCentralityEstimator("V0M");
        pidTask[i]->SetDCAxyCut(10);
        pidTask[i]->SetDCAzCut(10);
        pidTask[i]->SetCuts(PIDcuts);
        pidTask[i]->SetPIDPurityFunctions(0.8);
        
        
        mgr->AddTask(pidTask[i]);
        
        coutput1[i] = mgr->CreateContainer(Form("PID_%.f-%.f",centrMin[icentr],centrMax[icentr]), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile[i]);
        
        //connect containers
        mgr->ConnectInput  (pidTask[i],  0, mgr->GetCommonInputContainer());
        mgr->ConnectOutput (pidTask[i],  1, coutput1[i]);
        
        //return pidTask[icentr];
        
    }
}