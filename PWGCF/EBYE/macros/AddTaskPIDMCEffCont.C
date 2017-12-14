AliAnalysisTaskPIDMCEffCont *AddTaskPIDMCEffCont( TString  centralityEstimator="V0M",
                                                UInt_t triggerSelectionString = AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral,
						Double_t partData = 1,
					        Double_t vertexZ=10., Int_t nsigma = 3,
                                                Double_t ptMin=0, Double_t ptMax=10,
                                                Double_t etaMin=-0.8, Double_t etaMax=0.8,
                                                Bool_t PIDNsigma = kFALSE,
						Bool_t PIDPDG = kFALSE,
						Bool_t DCAext = kFALSE,
                                                AliAnalysisTaskPIDMCEffCont::kParticleOfInterest particleType = AliAnalysisTaskPIDMCEffCont::kPion,
                                                Int_t AODfilterBit = 768,
                                                TString fileNameBase="AnalysisResults"
                                                ) {
    
    // Creates a balance function analysis task and adds it to the analysis manager.
    // Get the pointer to the existing analysis manager via the static access method.
    
    TString centralityName("");
    
    static const Int_t nCentralitiespart1 = 2;
    static const Int_t nCentralitiespart2 = 5;
    static const Int_t nCentralitiespart3 = 4;

    Double_t gCentrality[nCentralitiespart2];

    if(partData == 1) {

             gCentrality[0]=0;
             gCentrality[1]=10;
        }

    else if(partData == 2) {

            gCentrality[0] = 10;
            gCentrality[1] = 20;
            gCentrality[2] = 30;
            gCentrality[3] = 40;
            gCentrality[4] = 50;
        }

     else if(partData == 3) {

            gCentrality[0] = 50;
            gCentrality[1] = 60;
            gCentrality[2] = 70;
            gCentrality[3] = 80;
        }

    Int_t nCentralities = 0;
    
    if(partData == 1) nCentralities = nCentralitiespart1;
    else if(partData == 2) nCentralities = nCentralitiespart2;
    else if(partData == 3) nCentralities = nCentralitiespart3;    


    TString suffixName[nCentralitiespart2];
    TString outputFileName(fileNameBase);
    outputFileName.Append(".root");
    AliAnalysisTaskPIDMCEffCont *task[nCentralitiespart2];
    AliAnalysisDataContainer *coutQA[nCentralitiespart2];
    AliAnalysisDataContainer *coutEff[nCentralitiespart2];

  //===========================================================================
    for(Int_t iCentralityBin = 0; iCentralityBin < nCentralities-1 ; iCentralityBin++) {
    
    suffixName[iCentralityBin] = "Centrality_";
    suffixName[iCentralityBin] += gCentrality[iCentralityBin];
    suffixName[iCentralityBin] += "To";
    suffixName[iCentralityBin] += gCentrality[iCentralityBin+1];
    suffixName[iCentralityBin] += "_fb";
    suffixName[iCentralityBin] += AODfilterBit;

    if(PIDNsigma)
    suffixName[iCentralityBin] += "_nsigmaPID";

    if(PIDPDG)
    suffixName[iCentralityBin] += "_nsigmaPDG";

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
    task[iCentralityBin]->SetRejectInjectedSignals();
    task[iCentralityBin]->SetCentralityPercentileRange(gCentrality[iCentralityBin],gCentrality[iCentralityBin+1]);
    task[iCentralityBin]->SelectCollisionCandidates(triggerSelectionString);
    
    // vertex
    task[iCentralityBin]->SetVertexDiamond(3.,3.,vertexZ);
    task[iCentralityBin]->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax);
    task[iCentralityBin]->SetAODtrackCutBit(AODfilterBit);
        
    if(DCAext)
	task[iCentralityBin]->USEextendedDCA();	    

    if(PIDNsigma) {
	task[iCentralityBin]->UsePIDNsigma();
        task[iCentralityBin]->SetNSigmaPID(nsigma);
        task[iCentralityBin]->setParticleType(particleType);
	mgr->AddTask(task[iCentralityBin]); 
   }
    
    if(PIDPDG) {  
	task[iCentralityBin]->UsePIDPDG();
    	task[iCentralityBin]->setParticleType(particleType);
        mgr->AddTask(task[iCentralityBin]); 
   }

    if(!PIDNsigma && !PIDPDG){
  	mgr->AddTask(task[iCentralityBin]);
   }	
       
    // Create ONLY the output containers for the data produced by the task.
    // Get and connect other common input/output containers via the manager as below
    //==============================================================================
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    outputFileName += ":BalanceFunctionEffContAnalysis";
    coutQA[iCentralityBin] = mgr->CreateContainer(Form("listQA_%s",suffixName[iCentralityBin].Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
    coutEff[iCentralityBin] = mgr->CreateContainer(Form("listEff_%s",suffixName[iCentralityBin].Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
    mgr->ConnectInput(task[iCentralityBin], 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task[iCentralityBin], 1, coutQA[iCentralityBin]);
    mgr->ConnectOutput(task[iCentralityBin], 2, coutEff[iCentralityBin]);

    }
}
