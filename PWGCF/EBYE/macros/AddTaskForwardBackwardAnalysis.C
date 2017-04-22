AliForwardBackwardAnalysis *AddTaskForwardBackwardAnalysis(
        TString  centralityEstimator="V0M",
//        Double_t centrMin=0.,
//        Double_t centrMax=80.,
//        Double_t vertexZ=10.,
//        Double_t DCAxy=-1,
//        Double_t DCAz=-1,
//        Double_t maxTPCchi2 = -1,
//        Int_t minNClustersTPC = -1,
//        Double_t ptmin = 0.2,
//        Double_t ptmax = 2.0,
//        Double_t pttpcmax = 0.6,
//        Double_t nsigmapid = 3.0,
//        int DetectorType = 0,
//        Bool_t tofMistMatch=kFALSE,
//        Double_t tofMisValue=.01,
        Int_t AODfilterBit = 128,
//        Bool_t RapidityUse=kFALSE,
//        Bool_t bUseElectronRejection = kFALSE,
//        TString fileNameBase = "AnalysisResults",
//        TString dirNameExtra = ""
        ) {

    // Creates a Forward-Backward correlations analysis task and adds it to the analysis manager.
    // Get the pointer to the existing analysis manager via the static access method.

    //    TString outputFileName(fileNameBase);
//    outputFileName.Append(".root");

    //===========================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskAliForwardBackwardAnalysis", "No analysis manager to connect to.");
        return NULL;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //===========================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskAliForwardBackwardAnalysis", "This task requires an input event handler");
        return NULL;
    }
    TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    if(dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())) analysisType = "MC";


    // Create the task, add it to manager and configure it.
    //===========================================================================
    AliForwardBackwardAnalysis *taskFB = new AliForwardBackwardAnalysis( Form( "AliForwardBackwardAnalysis_bit%d_%s"
                                                                                         , AODfilterBit, centralityEstimator.Data() ) );

    // centrality
    if( centralityEstimator )
        taskFB->SetCentralityEstimator(centralityEstimator);

    taskFB->SelectCollisionCandidates(AliVEvent::kMB);



    // vertex
//    taskFB->SetVertexDiamond(.3,.3,vertexZ);

    //analysis kinematic cuts
//    taskFB->SetMinPt(ptmin);
//    taskFB->SetMaxPt(ptmax); //5.0
//    taskFB->SetTPCPtMax(pttpcmax); // For TPC Max
//    //taskFB->SetEtaRange(-0.8,0.8,100,0.0,1.6, 64); //acceptance cuts
//    //taskFB->SetPtRange(0.1, 20.0, 100);  //acceptance cuts //5.0,49
//    taskFB->SetEtaRange(-0.8,0.8,100,0.0,1.6, 64); //acceptance cuts
//    taskFB->SetPtRange(0.2, 20.0, 100);  //acceptance cuts //5.0,49
//    taskFB->SetNSigmaCut(nsigmapid); // NSigma Cut;
//    taskFB->SetDetectorPID(DetectorType);

//    taskFB->SetMisMatchTOFProb(tofMisValue,tofMistMatch);
//    taskFB->SetRapidityUse(RapidityUse);


//    taskFB->SetExtraDCACutsAOD(DCAxy,DCAz);

    // set extra TPC chi2 / nr of clusters cut
//    taskFB->SetExtraTPCCutsAOD(maxTPCchi2, minNClustersTPC);

    // electron rejection
//    if(bUseElectronRejection){
//        taskFB->SetElectronOnlyRejection(3.); // no other particle in nsigma (this is what we use standard in BF code)
//    }

    //AODs
    taskFB->SetAODtrackCutBit(AODfilterBit);
    mgr->AddTask(taskFB);

    // Create ONLY the output containers for the data produced by the task.
    // Get and connect other common input/output containers via the manager as below
    //==============================================================================
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    outputFileName += ":PWGCFEbyE.outputFBanalysis";

    AliAnalysisDataContainer *coutput = mgr->CreateContainer( Form("list_cEst%s_bit%d",centralityEstimator.Data(),AODfilterBit), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data() );
//    AliAnalysisDataContainer *coutEffContBF = mgr->CreateContainer(Form("listEffContBF_%s_Bit%d_%s%s",centralityName.Data(),AODfilterBit,centralityEstimator.Data(),dirNameExtra.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
//    AliAnalysisDataContainer *coutQATruthReco = mgr->CreateContainer(Form("listQARecoTruth_%s_Bit%d_%s%s",centralityName.Data(),AODfilterBit,centralityEstimator.Data(),dirNameExtra.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());

    mgr->ConnectInput( taskFB, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput( taskFB, 1, coutput);
//    mgr->ConnectOutput(taskFB, 2, coutEffContBF);
//    mgr->ConnectOutput(taskFB, 3, coutQATruthReco);

    return taskFB;
}
