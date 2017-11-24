AliAnalysisTaskPIDMCEffCont *AddTaskPIDMCEffCont( TString  centralityEstimator="V0M",
                                                UInt_t triggerSelectionString = AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral,
                                                Double_t vertexZ=10., Int_t nsigma = 3,
                                                Double_t ptMin=0, Double_t ptMax=10,
                                                Double_t etaMin=-0.8, Double_t etaMax=0.8,
                                                Bool_t PID = kFALSE,
                                                AliAnalysisTaskPIDMCEffCont::kParticleOfInterest particleType = AliAnalysisTaskPIDMCEffCont::kMuon,
                                                Int_t AODfilterBit = 768,
                                                TString fileNameBase="AnalysisResults"
                                                ) {
    
    // Creates a balance function analysis task and adds it to the analysis manager.
    // Get the pointer to the existing analysis manager via the static access method.
    
    TString centralityName("");
    
    static const Int_t nCentralities = 9;
    Double_t gCentrality[nCentralities];
    
    for(Int_t i = 1; i < 10; i++)
        gCentrality[i-1] = 0 + (i-1)*10;

    TString suffixName[nCentralities];
    TString outputFileName(fileNameBase);
    outputFileName.Append(".root");
    AliAnalysisTaskPIDMCEffCont *task[nCentralities];
    AliAnalysisDataContainer *coutQA[nCentralities];
    AliAnalysisDataContainer *coutEff[nCentralities];

  //===========================================================================
    for(Int_t iCentralityBin = 0; iCentralityBin < nCentralities-1 ; iCentralityBin++) {
    
    suffixName[iCentralityBin] = "Centrality_";
    suffixName[iCentralityBin] += gCentrality[iCentralityBin];
    suffixName[iCentralityBin] += "To";
    suffixName[iCentralityBin] += gCentrality[iCentralityBin+1];

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
    //if(dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())) analysisType = "MC";
    
    
    // Create the task, add it to manager and configure it.
    //===========================================================================
    
    task[iCentralityBin] = new AliAnalysisTaskPIDMCEffCont(Form("Task_%s",suffixName[iCentralityBin].Data()));
    
    // centrality
    
    task[iCentralityBin]->UseCentrality();
    task[iCentralityBin]->SetCentralityEstimator(centralityEstimator);

    task[iCentralityBin]->SetCentralityPercentileRange(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1]);
    task[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
    
    // vertex
    task[iCentralityBin]->SetVertexDiamond(3.,3.,vertexZ);
    task[iCentralityBin]->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax);
    task[iCentralityBin]->SetAODtrackCutBit(AODfilterBit);
    
    if(PID){
            task[iCentralityBin]->UsePID();
            task[iCentralityBin]->SetNSigmaPID(nsigma);
            task[iCentralityBin]->setParticleType(particleType);
    }
        
    mgr->AddTask(task[iCentralityBin]);
    
    // Create ONLY the output containers for the data produced by the task.
    // Get and connect other common input/output containers via the manager as below
    //==============================================================================
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    outputFileName += ":BalanceFunctionEffContAnalysis";
    coutQA[iCentralityBin] = mgr->CreateContainer(Form("listQA_%s",suffixName[iCentralityBin].Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());//centralityName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
    coutEff[iCentralityBin] = mgr->CreateContainer(Form("listEff_%s",suffixName[iCentralityBin].Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());//centralityName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
    mgr->ConnectInput(task[iCentralityBin], 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task[iCentralityBin], 1, coutQA[iCentralityBin]);
    mgr->ConnectOutput(task[iCentralityBin], 2, coutEff[iCentralityBin]);

    }
}
