
void AddTaskPIDconfig(Int_t CentralityTriggerSelection = AliVEvent::kMB, Double_t centralityMin=0, Double_t centralityMax=5 ,Double_t FilterBit=1, Bool_t PIDcuts=kFALSE,TString useroutputfile="output.root"){

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskPID", "No analysis manager to connect to.");
        return 0x0;
    }
    
    // standard with task
    printf("========================================================================================\n");
    printf("PID: Initialising AliAnalysisTaskPIDconfig\n");
    printf("========================================================================================\n");
    
    Double_t centrMin[9] = {0,5,10,20,30,40,50,60,70};
    Double_t centrMax[9] = {5,10,20,30,40,50,60,70,80};
    Bool_t Pass_Min = kFALSE;
    Bool_t Pass_Max = kFALSE;

    for(int i=0;i<9;i++)
    {
        if(centralityMin == centrMin[i]){
            const int iMin = i;
            Pass_Min = kTRUE;
        }
        if(centralityMax == centrMax[i]){
            const int iMax = i;
            Pass_Max = kTRUE;
        }
    }
    if(!Pass_Min || !Pass_Max){
        ::Error("centrality Min and Max don't match the defined ranges");
        return 0x0;
    }
    
    const int ncentr = iMax - iMin +1;
    TString outputfile[ncentr];
    AliAnalysisDataContainer *coutput1[ncentr];
    AliAnalysisTaskPIDconfig *pidTask[ncentr];
    int icentr = 0;
    
    
    
    for(int i=0;i<iMax-iMin+1;i++){
        icentr = i + iMin;
        outputfile[i] = useroutputfile;
        pidTask[i] = new AliAnalysisTaskPIDconfig(Form("pidTask_%.f-%.f",centrMin[icentr],centrMax[icentr]));
        pidTask[i]->SelectCollisionCandidates(CentralityTriggerSelection);
        pidTask[i]->SetCutTPCmultiplicityOutliersAOD(kTRUE);
        pidTask[i]->SetData2011(kFALSE);
        pidTask[i]->SetFilterBit(FilterBit);
        pidTask[i]->SetUseCentrality(kTRUE);
        pidTask[i]->SetCentralityPercentileMin(centrMin[icentr]);
        pidTask[i]->SetCentralityPercentileMax(centrMax[icentr]);
        pidTask[i]->SetCentralityEstimator("V0M");
        pidTask[i]->SetDCAxyCut(10);
        pidTask[i]->SetDCAzCut(10);
        pidTask[i]->SetCuts(PIDcuts);
        if(PIDcuts){
            TFile *ContoursFile = new TFile(Form("PurityHistContours_%.f-%.f.root",centrMin[icentr],centrMax[icentr]));

            Contourlist = new TDirectory;
            Contourlist=(TDirectory*)ContoursFile->Get("Filterbit1");
            if(!Contourlist){printf("The contour file is empty"); continue;}
    
            pidTask[i]->SetPIDcontoursList(Contourlist);
            
        }
        
        mgr->AddTask(pidTask[i]);
        
        coutput1[i] = mgr->CreateContainer(Form("PID_%.f-%.f",centrMin[icentr],centrMax[icentr]), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile[i]);
  
        //connect containers
        mgr->ConnectInput  (pidTask[i],  0, mgr->GetCommonInputContainer());
        mgr->ConnectOutput (pidTask[i],  1, coutput1[i]);
        
        //return pidTask[icentr];

    }
}