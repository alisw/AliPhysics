AliAnalysisTaskSECheckCharmHadronBkg *AddTaskCheckCharmHadronBkg(double ptmin = 0.,
                                                                 double ptmax = 50.,
                                                                 double ptbinwidth = 0.5,
                                                                 double massmin = 1.6,
                                                                 double massmax = 2.1,
                                                                 double massbinwidth = 0.002,
                                                                 TString cutfilename = "",
                                                                 TString suffix = "",
                                                                 bool applyML = false,
                                                                 TString configFileML = "")
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
        ::Error("AddTaskDplusBkg", "No analysis manager to connect to.");

    bool stdcuts = false;
    TFile *filecuts = NULL;
    if (cutfilename.EqualTo(""))
        stdcuts = true;
    else
    {
        filecuts = TFile::Open(cutfilename.Data());
        if (!filecuts || (filecuts && !filecuts->IsOpen()))
            std::cerr << "Input file not found : check your cut object" << std::endl;
    }

    //Analysis Task
    AliRDHFCutsDplustoKpipi *analysiscuts = new AliRDHFCutsDplustoKpipi();
    if (stdcuts)
        analysiscuts->SetStandardCutsPP2010();
    else
        analysiscuts = (AliRDHFCutsDplustoKpipi *)filecuts->Get("AnalysisCuts");

    AliAnalysisTaskSECheckCharmHadronBkg *dplusTaskBkg = new AliAnalysisTaskSECheckCharmHadronBkg("DplusBkgStudy", analysiscuts);
    dplusTaskBkg->SetMassLimits(massmin, massmax);
    dplusTaskBkg->SetMassBinWidth(massbinwidth);
    dplusTaskBkg->SetPtLimits(ptmin, ptmax);
    dplusTaskBkg->SetPtBinWidth(ptbinwidth);
    if (applyML)
    {
        dplusTaskBkg->SetDoMLApplication();
        dplusTaskBkg->SetMLConfigFile(configFileML);
    }
    mgr->AddTask(dplusTaskBkg);

    // Create containers for input/output
    TString inname = "cinputDplus";
    TString outname = "coutputDplus";
    TString cutsname = "coutputDplusCuts";
    inname += suffix.Data();
    outname += suffix.Data();
    cutsname += suffix.Data();

    AliAnalysisDataContainer *cinputDplus = mgr->CreateContainer(inname, TChain::Class(), AliAnalysisManager::kInputContainer);
    TString outputfile = AliAnalysisManager::GetCommonFileName();
    outputfile += ":PWGHF_D2H_BkgDplus";

    AliAnalysisDataContainer *coutputDplusCuts = mgr->CreateContainer(cutsname, TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
    AliAnalysisDataContainer *coutputDplus = mgr->CreateContainer(outname, TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());

    mgr->ConnectInput(dplusTaskBkg, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(dplusTaskBkg, 1, coutputDplus);
    mgr->ConnectOutput(dplusTaskBkg, 2, coutputDplusCuts);

    return dplusTaskBkg;
}
