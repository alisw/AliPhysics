
AliAnalysisTaskRecEff *AddTaskRecEff(TString  centralityEstimator="V0M",
						Double_t centrMin=0.,
						Double_t centrMax=90.,
						Double_t vertexZ=10.,
						Double_t ptMin = 0.1,
						Double_t ptMax = 10,
						Double_t etaMin = -0.8,
						Double_t etaMax = 0.8,
						Int_t NclusMin = 70,
						Int_t AODfilterBit = 768,
						TString fileNameBase="AnalysisResults",
						UInt_t triggerSel = AliVEvent::kINT7,
						TString extraString = ""
						) {

    // Creates a balance function analysis task and adds it to the analysis manager.
    // Get the pointer to the existing analysis manager via the static access method.
    TString centralityName("");
    centralityName+=Form("%.0f-%.0f_%.0f",centrMin,centrMax,vertexZ);
    TString outputFileName(fileNameBase);
    outputFileName.Append(".root");

    //===========================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskTriggeredBF", "No analysis manager to connect to.");
        return NULL;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //===========================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskTriggeredBF", "This task requires an input event handler");
        return NULL;
    }
    TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    if(dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())) analysisType = "MC";


    // Create the task, add it to manager and configure it.
    //===========================================================================

    AliAnalysisTaskRecEff *taskRecEff = new AliAnalysisTaskRecEff(Form("TaskEffContBF,FB%i_%s",AODfilterBit, centralityEstimator.Data()));

    // centrality
    if(centralityEstimator) {
        taskRecEff->SetCentralityEstimator(centralityEstimator);
        taskRecEff->SetCentralityPercentileMin(centrMin);
        taskRecEff->SetCentralityPercentileMax(centrMax);
    }

    taskRecEff->SelectCollisionCandidates(triggerSel);

    // vertex
	taskRecEff->SetVzRangeMin(-vertexZ);
	taskRecEff->SetVzRangeMax(vertexZ);
	
    //analysis kinematic cuts
    
    //taskRecEff->SetEtaRange(-0.8,0.8,100,0.0,1.6, 64); //acceptance cuts
    //taskRecEff->SetPtRange(0.1, 20.0, 100);  //acceptance cuts //5.0,49
    //taskRecEff->SetNSigmaCutTPC(4);
    taskRecEff->SetTrackCutdEdxMin(10);
    taskRecEff->SetTrackCutChi2Min(0); // max fix to 4.0
    taskRecEff->SetEtaRangeMin(etaMin);
    taskRecEff->SetEtaRangeMax(etaMax);
    taskRecEff->SetPtRangeMin(ptMin);
    taskRecEff->SetPtRangeMax(ptMax); //5.0
    taskRecEff->SetFilterBit(AODfilterBit);
    taskRecEff->SetTrackCutNclusterMin(NclusMin);
    
    TString pidsuffix ="ch";
    
    //AODs
    mgr->AddTask(taskRecEff);


    // Create ONLY the output containers for the data produced by the task.
    // Get and connect other common input/output containers via the manager as below
    //==============================================================================
    //TString outputFileName = AliAnalysisManager::GetCommonFileName();
    outputFileName += ":PWGCFEbyE.outputBalanceFunctionEffContAnalysis";
    AliAnalysisDataContainer *coutQA = mgr->CreateContainer(Form("listQA_FB%i_%s_%s_%s_%s",AODfilterBit,centralityEstimator.Data(),centralityName.Data(), pidsuffix.Data(), extraString.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
    AliAnalysisDataContainer *coutEffContBF = mgr->CreateContainer(Form("listEffContBF_FB%i_%s_%s_%s_%s",AODfilterBit,centralityEstimator.Data(),centralityName.Data(), pidsuffix.Data(), extraString.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());

    mgr->ConnectInput(taskRecEff, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskRecEff, 1, coutQA);
    mgr->ConnectOutput(taskRecEff, 2, coutEffContBF);

    return taskRecEff;
}
