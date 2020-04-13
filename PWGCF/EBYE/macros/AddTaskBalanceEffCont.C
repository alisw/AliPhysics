
AliAnalysisTaskEffContBF *AddTaskBalanceEffCont(TString  centralityEstimator="V0M",
						Double_t centrMin=10.,
						Double_t centrMax=80.,
						Double_t vertexZ=10.,
						Int_t AODfilterBit = 96,
						Bool_t bUseElectronRejection = kFALSE,
						TString fileNameBase="AnalysisResults",
						UInt_t triggerSel = AliVEvent::kINT7,
						Bool_t usePIDfromPDG=kTRUE,
						AliPID::EParticleType particleType = AliPID::kPion,
						Bool_t usePIDstrategy=kFALSE,
						Bool_t usePIDnSigmaComb=kTRUE,
						Double_t BayesThr = 0.8,
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

    AliAnalysisTaskEffContBF *taskEffContBF = new AliAnalysisTaskEffContBF(Form("TaskEffContBF,FB%i_%s",AODfilterBit, centralityEstimator.Data()));

    // centrality
    if(centralityEstimator) {
        taskEffContBF->UseCentrality();
        taskEffContBF->SetCentralityEstimator(centralityEstimator);
        taskEffContBF->SetCentralityPercentileRange(centrMin,centrMax);
    }

    taskEffContBF->SelectCollisionCandidates(triggerSel);

    // vertex
    taskEffContBF->SetVertexDiamond(1.,1.,vertexZ);

    //analysis kinematic cuts
    taskEffContBF->SetMinPt(0.0);
    taskEffContBF->SetMaxPt(20.0); //5.0
    //taskEffContBF->SetEtaRange(-0.8,0.8,100,0.0,1.6, 64); //acceptance cuts
    //taskEffContBF->SetPtRange(0.1, 20.0, 100);  //acceptance cuts //5.0,49
    taskEffContBF->SetEtaRange(-0.8,0.8,100,0.0,1.6, 64); //acceptance cuts
    taskEffContBF->SetPtRange(0.0, 20.0, 100);  //acceptance cuts //5.0,49

    TString pidsuffix ="ch";
    if (usePIDfromPDG) {
        pidsuffix = AliPID::ParticleShortName(particleType);
        taskEffContBF->SetUsePIDfromPDG(usePIDfromPDG, particleType);

    }

    if (usePIDstrategy) {
      pidsuffix = AliPID::ParticleShortName(particleType);
      taskEffContBF->SetUseParticleID(usePIDstrategy, particleType);
      taskEffContBF->SetUsePIDnSigmaComb(usePIDnSigmaComb);
      taskEffContBF->SetBayesPIDThr(BayesThr);
    }
    //AODs
    taskEffContBF->SetAODtrackCutBit(AODfilterBit);
    mgr->AddTask(taskEffContBF);


    // Create ONLY the output containers for the data produced by the task.
    // Get and connect other common input/output containers via the manager as below
    //==============================================================================
    //TString outputFileName = AliAnalysisManager::GetCommonFileName();
    outputFileName += ":PWGCFEbyE.outputBalanceFunctionEffContAnalysis";
    AliAnalysisDataContainer *coutQA = mgr->CreateContainer(Form("listQA_FB%i_%s_%s_%s_%s",AODfilterBit,centralityEstimator.Data(),centralityName.Data(), pidsuffix.Data(), extraString.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
    AliAnalysisDataContainer *coutEffContBF = mgr->CreateContainer(Form("listEffContBF_FB%i_%s_%s_%s_%s",AODfilterBit,centralityEstimator.Data(),centralityName.Data(), pidsuffix.Data(), extraString.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());

    mgr->ConnectInput(taskEffContBF, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskEffContBF, 1, coutQA);
    mgr->ConnectOutput(taskEffContBF, 2, coutEffContBF);

    return taskEffContBF;
}
