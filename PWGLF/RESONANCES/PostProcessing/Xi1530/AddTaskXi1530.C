AliAnalysisTaskXi1530* AddTaskXi1530(const char *taskname = "Xi1530"
                                     , const char *option = "SYS"
                                     , int nmix=20
                                     , const char* suffix = "MB")
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    TString foption = option;
    AliAnalysisTaskXi1530 *taskXi1530 = new AliAnalysisTaskXi1530(Form("%s%s", taskname,suffix), option);
    taskXi1530->fEventCuts.fCentralityFramework = 1;
    taskXi1530->fEventCuts.SetMaxVertexZposition(10);
    taskXi1530->SetnMix(nmix);

    std::cout << "AliAnaylsisTaskXi1530:: Option: " << option << std::endl;
    if(foption.Contains("MC")){
        taskXi1530->SetIsMC(kTRUE); // default: kFALSE
        taskXi1530->SetMixing(kFALSE);  // default: kTRUE
        std::cout << "AliAnaylsisTaskXi1530:: MC mode " << std::endl;
        if (foption.Contains("Gen")){
            taskXi1530->SetIsPrimaryMC(kFALSE); // default: kTRUE
            std::cout << "<GENERAL PURPOSE MC>" << std::endl;
        }
    }
    if(foption.Contains("AA")){
        taskXi1530->SetIsAA(kTRUE); // default: kFALSE
        std::cout << "AliAnaylsisTaskXi1530:: AA mode " << std::endl;
    }
    if(foption.Contains("pA")){
        taskXi1530->SetIsAA(kTRUE); // default: kFALSE
        //taskXi1530->SetXi1530RapidityCut_high(0.5); // default: 0.5
        taskXi1530->SetXi1530RapidityCut_low(0); // default: -0.5
        std::cout << "AliAnaylsisTaskXi1530:: pA mode " << std::endl;
    }
    if(foption.Contains("Ap")){
        taskXi1530->SetIsAA(kTRUE); // default: kFALSE
        taskXi1530->SetXi1530RapidityCut_high(0); // default: 0.5
        //taskXi1530->SetXi1530RapidityCut_low(-0.5); // default: -0.5
        std::cout << "AliAnaylsisTaskXi1530:: Ap mode " << std::endl;
    }
    if(foption.Contains("NoMix")){
        taskXi1530->SetMixing(kFALSE); // default: kTRUE
        std::cout << "AliAnaylsisTaskXi1530:: Event Mix disabled " << std::endl;
    } 
    if(foption.Contains("HM")){
        taskXi1530->SetHighMult(kTRUE); // default: kFALSE
        taskXi1530->fEventCuts.fTriggerMask =
            AliVEvent::kHighMultV0;  // default: kINT7
        std::cout << "AliAnaylsisTaskXi1530:: HighMultV0 mode " << std::endl;
    }  
    if(foption.Contains("SYS")){
        taskXi1530->SetSystematics(kTRUE); // default: kFALSE
        std::cout << "AliAnaylsisTaskXi1530:: Systematic Study mode " << std::endl;
    }
    if (foption.Contains("NoQA")) {
        taskXi1530->SetQA(kFALSE);  // default: kTRUE
        std::cout << "AliAnaylsisTaskXi1530:: NoQA mode " << std::endl;
    }
    if (foption.Contains("EXO")) {
        taskXi1530->SetExoticFinder(kTRUE);  // default: kFALSE
        std::cout << "AliAnaylsisTaskXi1530:: ExoticFinder mode " << std::endl;
    }
    if (foption.Contains("HF")) {
        taskXi1530->SetExoticFinder2(kTRUE);  // default: kFALSE
        std::cout << "AliAnaylsisTaskXi1530:: ExoticFinder2 mode " << std::endl;
    }
    if (foption.Contains("INEL")) {
        taskXi1530->fEventCuts.fCentralityFramework = 0;
        taskXi1530->fEventCuts.SelectOnlyInelGt0(false);
        taskXi1530->SetINEL(kTRUE);  // default: kFALSE
        std::cout << "AliAnaylsisTaskXi1530:: Inelastic mode " << std::endl;
    }

    if(!taskXi1530) return 0x0;
    mgr->AddTask(taskXi1530);
    
    
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutputXi1530 = mgr->CreateContainer(Form("%s%s", taskname, suffix), TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

    mgr->ConnectInput(taskXi1530, 0, cinput);
    mgr->ConnectOutput(taskXi1530, 1, coutputXi1530);

    return taskXi1530;
}

