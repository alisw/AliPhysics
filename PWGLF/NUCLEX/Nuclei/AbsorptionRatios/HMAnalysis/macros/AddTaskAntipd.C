AliAnalysisTaskAntipd* AddTaskAntipd( TString name = "TaskAntiPD",
                                     int filterBit = 256,
                                     float lowp = 0.1,
                                     float eta  = 0.8,
                                     int minITScl = 2,
                                     float maxDCAxy = 1.0,
                                     float maxDCAz = 0.2,
                                     float maxTPCnSig = 3.0,
                                     float maxTOFnSig = 3.0,
                                     float momTOFanaProt = 0.7,
                                     float momTOFanaDeut = 1.4,
                                     float minITSnSigmaDeut = -2.0,
                                     float maxITSnSigmaDeut = 1e30,
                                     const char* suffix =""
                                     )
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }

    // create an instance of the task
	AliAnalysisTaskAntipd *taskantipd = new AliAnalysisTaskAntipd(name.Data());
    if(!taskantipd) return 0x0;

    // event selection
     //taskantipd->SelectCollisionCandidates(triggerMask); //  trigger

    // set AOD filter bit and further track cuts
    taskantipd->SetFilterBit(filterBit);
    taskantipd->SetLowPCut(lowp);
    taskantipd->SetHighPCut(1e30);
    taskantipd->SetEtaCut(eta);
    taskantipd->SetMinNITSCl(minITScl);
    taskantipd->SetMaxDCAxy(maxDCAxy);
    taskantipd->SetMaxDCAz(maxDCAz);

    // set PID cuts
    taskantipd->SetMaxTPCnSigma(maxTPCnSig);
    //taskantipd->SetUseTOFPidCut(UseTOFcut);
    //if (UseTOFcut) taskantipd->SetMaxTOFnSigma(maxTOFnSig);
    // momentum p from which a hit/cut in TOF is required
    taskantipd->SetMomForTOFanaProt(momTOFanaProt);
    taskantipd->SetMomForTOFanaDeut(momTOFanaDeut);
    // ITS n sigma cut up to TOF momentum
    taskantipd->SetITSnSigmaRange(minITSnSigmaDeut,maxITSnSigmaDeut);

    // add task to the manager

    mgr->AddTask(taskantipd);
    //TString file = "AnalysisResults.root";
    //TString OutputFile = Form("Output_%s",suffix);


   //mgr->ConnectInput(taskantipd,0,mgr->GetCommonInputContainer());
   //mgr->ConnectOutput(taskantipd,1,mgr->CreateContainer(OutputFile.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s%s",file.Data(),OutputFile.Data())));

   TString list = Form("Output_%s",suffix);
   TString file = "AnalysisResults.root";
   // connect the manager to the task
   mgr->ConnectInput(taskantipd,0,mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskantipd,1,mgr->CreateContainer(list, TList::Class(), AliAnalysisManager::kOutputContainer, file));

    return taskantipd;
}
