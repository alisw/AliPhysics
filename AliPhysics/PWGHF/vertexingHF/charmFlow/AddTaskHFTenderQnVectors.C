AliAnalysisTaskSEHFTenderQnVectors* AddTaskHFTenderQnVectors(TString taskname = "HFTenderQnVectors",
                                                             TString outputSuffix = "", 
                                                             int harmonic = 2, 
                                                             int normmethod = 1,//AliHFQnVectorHandler::kQoverM,
                                                             int calibType = 0,//AliHFQnVectorHandler::kQnCalib, 
                                                             TString AODBfileName = "",
                                                             TString qnSplineFileName = "",
                                                             bool enableTPCphivscentrhistos=false,
                                                             bool enableTPCQvecvscentrhistos=false)
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AliAnalysisTaskSEHFTenderQnVectors", "No analysis manager found.");
        return NULL;
    }
    if (!mgr->GetInputEventHandler()) {
        ::Error("AliAnalysisTaskSEHFTenderQnVectors", "This task requires an input event handler");
        return NULL;
    }

    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    if(type.Contains("ESD")){
        ::Error("AliAnalysisTaskSEHFTenderQnVectors", "This task requires to run on AOD");
        return NULL;
    }

    //========= Add task for standard analysis to the ANALYSIS manager ====
    AliAnalysisTaskSEHFTenderQnVectors *task = new AliAnalysisTaskSEHFTenderQnVectors(taskname.Data(),harmonic,calibType,AODBfileName);
    task->SetNormalisationMethod(normmethod);
    if(qnSplineFileName!="") task->LoadSplinesForqnPercentile(qnSplineFileName);
    if(enableTPCphivscentrhistos) task->EnableTPCPhiVsCentrDistrHistosVsRun();
    if(enableTPCQvecvscentrhistos) task->EnableQVecTPCVsCentrDistrHistosVsRun();
    mgr->AddTask(task);

    TString outputfile = AliAnalysisManager::GetCommonFileName();
    outputfile += ":PWGHF_D2H_QnVectorTender";

    //define input container
    AliAnalysisDataContainer *cinput = mgr->CreateContainer(Form("cinputQnVectorTender%s",outputSuffix.Data()),TChain::Class(),AliAnalysisManager::kInputContainer);
    //define output containers
    AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("coutputQnVectorTender%s",outputSuffix.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutputphidistr[2], *coutputQvecdistr[3];
    if(calibType==0 && enableTPCphivscentrhistos) {
        coutputphidistr[0] = mgr->CreateContainer("coutputQnVectorTenderPhiDistrPosEta", TH2F::Class(),AliAnalysisManager::kOutputContainer,Form("%sPhiDistr",outputfile.Data()));
        coutputphidistr[1] = mgr->CreateContainer("coutputQnVectorTenderPhiDistrNegEta", TH2F::Class(),AliAnalysisManager::kOutputContainer,Form("%sPhiDistr",outputfile.Data()));
    }
    if(calibType==0 && enableTPCQvecvscentrhistos) {
        coutputQvecdistr[0] = mgr->CreateContainer("coutputQnVectorTenderQvecDistrFullEta", TH2F::Class(),AliAnalysisManager::kOutputContainer,Form("%sPhiDistr",outputfile.Data()));
        coutputQvecdistr[1] = mgr->CreateContainer("coutputQnVectorTenderQvecDistrPosEta", TH2F::Class(),AliAnalysisManager::kOutputContainer,Form("%sPhiDistr",outputfile.Data()));
        coutputQvecdistr[2] = mgr->CreateContainer("coutputQnVectorTenderQvecDistrNegEta", TH2F::Class(),AliAnalysisManager::kOutputContainer,Form("%sPhiDistr",outputfile.Data()));
    }

    //connect containers
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, coutput);
    if(calibType==0 && enableTPCphivscentrhistos) {
        mgr->ConnectOutput(task, 2, coutputphidistr[0]);
        mgr->ConnectOutput(task, 3, coutputphidistr[1]);
    }
    if(calibType==0 && enableTPCQvecvscentrhistos) {
        mgr->ConnectOutput(task, 4, coutputQvecdistr[0]);
        mgr->ConnectOutput(task, 5, coutputQvecdistr[1]);
        mgr->ConnectOutput(task, 6, coutputQvecdistr[2]);
    }
    return task;
}
