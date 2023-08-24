#include "AliVEvent.h"

AliAnalysisTaskSECharmHadronMLSelector *AddTaskCharmHadronMLSelector(TString fileName = "alien:///alice/cern.ch/user/f/fgrosa/CutObjects/pp13TeVHM/DplustoKpipiCuts_pp_femto_loose_kHighMultV0.root",
                                                                     TString MLconfigFileName = "alien:///alice/cern.ch/user/f/fgrosa/MLconfigfiles/config_DplusFemto_applyML_pp13TeVHM.yml",
                                                                     int decCh = AliAnalysisTaskSECharmHadronMLSelector::kDplustoKpipi,
                                                                     TString cutObjName = "analysisCut",
                                                                     TString trigClass = "",
                                                                     unsigned long long trigMask = AliVEvent::kINT7,
                                                                     TString suffix = "",
                                                                     int useAODProtection = 0)
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        ::Error("AddTaskCharmHadronMLSelector", "No analysis manager to connect to.");
        return NULL;
    }
    if (!mgr->GetInputEventHandler())
    {
        ::Error("AddTaskCharmHadronMLSelector", "This task requires an input event handler");
        return NULL;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    if (type.Contains("ESD"))
    {
        ::Error("AddTaskCharmHadronMLSelector", "This task requires to run on AOD");
        return NULL;
    }

    TFile *filecuts = TFile::Open(fileName.Data());
    if (!filecuts || (filecuts && !filecuts->IsOpen()))
    {
        Printf("FATAL: Input file not found : check your cut object");
        return NULL;
    }

    //Analysis cuts
    AliRDHFCuts *analysisCut = NULL;
    TString taskName = "MLSelector";
    if (decCh == AliAnalysisTaskSECharmHadronMLSelector::kDplustoKpipi)
    {
        analysisCut = (AliRDHFCutsDplustoKpipi *)filecuts->Get(cutObjName);
        taskName += "Dplus";
        suffix.Prepend("Dplus");
    }
    else if (decCh == AliAnalysisTaskSECharmHadronMLSelector::kDstoKKpi)
    {
        analysisCut = (AliRDHFCutsDstoKKpi *)filecuts->Get(cutObjName);
        taskName += "Ds";
        suffix.Prepend("Ds");
    }
    else if (decCh == AliAnalysisTaskSECharmHadronMLSelector::kD0toKpi)
    {
        analysisCut = (AliRDHFCutsD0toKpi *)filecuts->Get(cutObjName);
        taskName += "D0";
        suffix.Prepend("D0");
    }
    else if (decCh == AliAnalysisTaskSECharmHadronMLSelector::kDstartoD0pi)
    {
        analysisCut = (AliRDHFCutsDStartoKpipi *)filecuts->Get(cutObjName);
        taskName += "DStar";
        suffix.Prepend("DStar");
    }
    if (!analysisCut)
    {
        Printf("FATAL: Specific AliRDHFCuts not found");
        return NULL;
    }


    AliAnalysisTaskSECharmHadronMLSelector *task = new AliAnalysisTaskSECharmHadronMLSelector(taskName, decCh, analysisCut);
    task->SetAODMismatchProtection(useAODProtection);
    task->SetTriggerInfo(trigClass, trigMask);
    task->SetMLConfigFile(MLconfigFileName);
    mgr->AddTask(task);

    TString outFileName = AliAnalysisManager::GetCommonFileName();
    outFileName += ":PWGHF_D2H_CharmMLSelector";

    //define input container
    AliAnalysisDataContainer *cinput = mgr->CreateContainer(Form("cinputCharmMLSelector%s", suffix.Data()), TChain::Class(), AliAnalysisManager::kInputContainer);
    //define output containers
    AliAnalysisDataContainer *coutPut = mgr->CreateContainer(Form("coutputCharmMLSelector%s", suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, outFileName.Data());
    AliAnalysisDataContainer *coutPutCuts = mgr->CreateContainer(Form("coutputCutsCharmMLSelector%s", suffix.Data()), AliRDHFCuts::Class(), AliAnalysisManager::kOutputContainer, outFileName.Data());

    //connect containers
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, coutPut);
    mgr->ConnectOutput(task, 2, coutPutCuts);

    return task;
}
