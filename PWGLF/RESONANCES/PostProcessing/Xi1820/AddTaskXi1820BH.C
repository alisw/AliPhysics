#include "AliAnalysisTaskXi1820BH.h"
AliAnalysisTaskXi1820BH* AddTaskXi1820BH(
    const char* taskname = "Xi1820",
    const char* option = "MB_Mix",
    int nmix = 10,
    const char* suffix = "MB") {
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    TString foption = option;
    Bool_t IsMC = kFALSE;
    if (foption.Contains("MC"))
        IsMC = kTRUE;
    AliAnalysisTaskXi1820BH* taskXi1820BH = new AliAnalysisTaskXi1820BH(Form("%s%s", taskname, suffix), IsMC);

    taskXi1820BH->fEventCuts.fCentralityFramework = 1;
    taskXi1820BH->fEventCuts.SetMaxVertexZposition(10);
    taskXi1820BH->fEventCuts.SelectOnlyInelGt0(false);
    std::cout << "AddtaskXi1820BH:: Option: " << option << std::endl;
    if (foption.Contains("MC")) {
        std::cout << "AliAnalysisTaskXi1820BH:: MC mode " << std::endl;
        if (foption.Contains("Gen")) {
            taskXi1820BH->SetIsPrimaryMC(kFALSE);  // default: kTRUE
            std::cout << "<GENERAL PURPOSE MC>" << std::endl;
        }
    }
    if (foption.Contains("pA")) {
        taskXi1820BH->SetXi1820RapidityCutLow(0);  // default: -0.5
        std::cout << "AliAnalysisTaskXi1820BH:: pA mode " << std::endl;
    }
    if (foption.Contains("Ap")) {
        taskXi1820BH->SetXi1820RapidityCutHigh(0);  // default: 0.5
        std::cout << "AliAnalysisTaskXi1820BH:: Ap mode " << std::endl;
    }
    if (foption.Contains("HM")) {
        taskXi1820BH->SetHighMult(kTRUE); // default: kFALSE
        taskXi1820BH->fEventCuts.fTriggerMask =
            AliVEvent::kHighMultV0;  // default: kINT7
        std::cout << "AliAnalysisTaskXi1820BH:: HighMultV0 mode "
                  << std::endl;
    }
    // Mixing
    if (foption.Contains("Mix")) {
        taskXi1820BH->SetnMix(nmix);
        std::cout << "Event Mix mode: " << nmix << "times" << std::endl;
    }

    if (!taskXi1820BH)
        return 0x0;
    mgr->AddTask(taskXi1820BH);

    AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
    mgr->ConnectInput(taskXi1820BH, 0, cinput);

    AliAnalysisDataContainer* coutputXi1820 = mgr->CreateContainer(
        Form("%s%s", taskname, suffix), TList::Class(),
        AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
    mgr->ConnectOutput(taskXi1820BH, 1, coutputXi1820);
    return taskXi1820BH;
}
