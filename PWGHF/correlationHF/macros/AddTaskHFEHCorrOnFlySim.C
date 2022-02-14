
AliAnalysisTask *AddTaskHFEHCorrOnFlySim(TString suffixName ="")
{
    
    AliAnalysisHFEHCorrOnFlySim* clus = new  AliAnalysisHFEHCorrOnFlySim("");
    clus->SetEtaRange(-20.0, 20.0);
    clus->SetPtRange(0.3, 1000.0);
    clus->SetYRange(-20., 20.);
    clus->SetMultRange(0,5000);
    clus->SetEventProperties(kTRUE);
    clus->SetPartProperties(kTRUE);
    clus->SetHFCorrelations(kTRUE);
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
      Printf("AliAnalysisHFEHCorrOnFlySim, No analysis manager to connect to.");
      return NULL;
    }
   
    if(!mgr->GetMCtruthEventHandler()){
      Printf("AliAnalysisHFEHCorrOnFlySim; This task requires an input MC event handler");
      return NULL;
    }
    
    mgr->AddTask(clus);
    
    // Create containers for input/output
    TString finDirname   = "";
    TString inname       = "cinput";
    TString outBasic     = "BasicPlots";
    TString Specific     = "Specific";
    
    finDirname += suffixName.Data();
    inname           +=   finDirname.Data();
    outBasic         +=   finDirname.Data();
    Specific         +=   finDirname.Data();

    //Input and Output Slots:
    AliAnalysisDataContainer *cinputSim = mgr->CreateContainer(inname,TChain::Class(), AliAnalysisManager::kInputContainer);
    TString outputfile = AliAnalysisManager::GetCommonFileName();
    outputfile += ":KineSimulations";
    //TString outputfile = "AnaKineResults.root";

    AliAnalysisDataContainer *coutputSim1 = mgr->CreateContainer(outBasic,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutputSim2 = mgr->CreateContainer(Specific,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

    mgr->ConnectInput(clus,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(clus,1,coutputSim1);
    mgr->ConnectOutput(clus,2,coutputSim2);

    return clus;
}
