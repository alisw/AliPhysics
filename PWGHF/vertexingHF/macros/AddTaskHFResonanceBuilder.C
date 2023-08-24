#if !defined (__CINT__) || defined (__CLING__)
#include <vector>
#include <TFile.h>
#include <TString.h>
#include <TList.h>

#include "AliRDHFCuts.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliAnalysisTaskSEHFResonanceBuilder.h"
#endif

AliAnalysisTaskSEHFResonanceBuilder *AddTaskHFResonanceBuilder(int decCh = AliAnalysisTaskSEHFResonanceBuilder::kDplustoKpipi,
                                                               bool readMC = false,
                                                               TString fileName = "",
                                                               TString suffix = "",
                                                               TString cutObjName = "AnalysisCuts",
                                                               bool enablePi = true,
                                                               bool enableKa = false,
                                                               bool enablePr = false,
                                                               bool enableDe = false,
                                                               bool enableKz = false,
                                                               bool enableLa = false,
                                                               float nSigmaTPCPi = 3.,
                                                               float nSigmaTPCKa = 0.,
                                                               float nSigmaTPCPr = 0.,
                                                               float nSigmaTPCDe = 0.,
                                                               float nSigmaTOFPi = 3.,
                                                               float nSigmaTOFKa = 0.,
                                                               float nSigmaTOFPr = 0.,
                                                               float nSigmaTOFDe = 0.,
                                                               float ptMinBach = 0.1,
                                                               std::vector<float> massMinPi = {1.9},
                                                               std::vector<float> massMaxPi = {3.2},
                                                               std::vector<float> massMinKa = {-1.},
                                                               std::vector<float> massMaxKa = {-1.},
                                                               std::vector<float> massMinPr = {-1.},
                                                               std::vector<float> massMaxPr = {-1.},
                                                               std::vector<float> massMinDe = {-1.},
                                                               std::vector<float> massMaxDe = {-1.},
                                                               std::vector<float> massMinKz = {-1.},
                                                               std::vector<float> massMaxKz = {-1.},
                                                               std::vector<float> massMinLa = {-1.},
                                                               std::vector<float> massMaxLa = {-1.},
                                                               int AODProtection = 0,
                                                               bool applyMultWeights = false,
                                                               std::string fileNameMultWeights = "",
                                                               std::string histoMultWeights = "")
{
    // \brief: AddTask for AliAnalysisTaskSEHFResonanceBuilder
    // \authors:
    // F. Grosa, fabrizio.grosa@cern.ch
    //==============================================================================

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
        ::Error("AddTaskHFResonanceBuilder", "No analysis manager to connect to.");

    TFile *fileCuts = TFile::Open(fileName.Data());
    if (!fileCuts || (fileCuts && !fileCuts->IsOpen()))
        ::Fatal("AddTaskHFResonanceBuilder", "Cut file not found on Grid: analysis will not start!\n");
    else
        printf("Cut file correctly found\n");

    AliRDHFCuts *analysisCuts = NULL;
    if(decCh == AliAnalysisTaskSEHFResonanceBuilder::kD0toKpi)
    {
        analysisCuts = (AliRDHFCutsD0toKpi*)fileCuts->Get(cutObjName.Data());
        suffix.Prepend("D0toKpi");
    }
    else if(decCh == AliAnalysisTaskSEHFResonanceBuilder::kDplustoKpipi)
    {
        analysisCuts = (AliRDHFCutsDplustoKpipi*)fileCuts->Get(cutObjName.Data());
        suffix.Prepend("DplustoKpipi");
    }
    else if(decCh == AliAnalysisTaskSEHFResonanceBuilder::kDstartoD0pi)
    {
        analysisCuts = (AliRDHFCutsDStartoKpipi*)fileCuts->Get(cutObjName.Data());
        suffix.Prepend("DstartoD0pi");
    }
    else
    {
        ::Fatal("AddTaskHFResonanceBuilder", "Specific AliRDHFCuts not found!\n");
        return NULL;
    }

    TH1F* hMultWeights = NULL;
    if (applyMultWeights) {
        TFile *fileMultWeights = TFile::Open(fileNameMultWeights.data());
        if (!fileMultWeights || (fileMultWeights && !fileMultWeights->IsOpen()))
            ::Fatal("AddTaskHFResonanceBuilder", "Multiplicity weight file not found!\n");
        hMultWeights = (TH1F*)fileMultWeights->Get(histoMultWeights.data());
        if (!hMultWeights)
            ::Fatal("AddTaskHFResonanceBuilder", "Multiplicity weight histo not found!\n");
    }

    // Analysis Task
    AliAnalysisTaskSEHFResonanceBuilder *hfResoTask = new AliAnalysisTaskSEHFResonanceBuilder("HFResonanceBuilderAnalysis", decCh, analysisCuts);
    hfResoTask->SetPtBachelorSelection(ptMinBach);
    hfResoTask->SetNsigmaBachelorSelection(nSigmaTPCPi, nSigmaTPCKa, nSigmaTPCPr, nSigmaTPCDe, nSigmaTOFPi, nSigmaTOFKa, nSigmaTOFPr, nSigmaTOFDe);
    hfResoTask->SetCharmResoMassWindows(massMinPi, massMaxPi, massMinKa, massMaxKa, massMinPr, massMaxPr, massMinDe, massMaxDe, massMinKz, massMaxKz, massMinLa, massMaxLa);
    hfResoTask->SetAODMismatchProtection(AODProtection);
    hfResoTask->SetReadMC(readMC);
    hfResoTask->EnableBachelors(enablePi, enableKa, enablePr, enableDe);
    hfResoTask->EnableV0s(enableKz, enableLa);
    if (applyMultWeights) {
        hfResoTask->SetMultiplicityWeights(hMultWeights);
    }
    mgr->AddTask(hfResoTask);

    // Create containers for input/output
    TString name = Form("cinputTree%s", suffix.Data());
    AliAnalysisDataContainer *cinput = mgr->CreateContainer(name, TChain::Class(), AliAnalysisManager::kInputContainer);
    TString outputfile = AliAnalysisManager::GetCommonFileName();
    outputfile += ":PWGHF_D2H_HFResoBuilder";
    outputfile += suffix.Data();

    name = Form("coutputCuts%s", suffix.Data());
    AliAnalysisDataContainer *coutputCuts = NULL;
    if(decCh == AliAnalysisTaskSEHFResonanceBuilder::kD0toKpi)
    {
        coutputCuts = mgr->CreateContainer(name, AliRDHFCutsD0toKpi::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
    }
    else if(decCh == AliAnalysisTaskSEHFResonanceBuilder::kDplustoKpipi)
    {
        coutputCuts = mgr->CreateContainer(name, AliRDHFCutsDplustoKpipi::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
    }
    else if(decCh == AliAnalysisTaskSEHFResonanceBuilder::kDstartoD0pi)
    {
        coutputCuts = mgr->CreateContainer(name, AliRDHFCutsDStartoKpipi::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
    }
    name = Form("coutput%s", suffix.Data());
    AliAnalysisDataContainer *coutput = mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
    name = Form("coutputNorm%s", suffix.Data());
    AliAnalysisDataContainer *coutputNorm = mgr->CreateContainer(name, AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
    name = Form("coutputTree%s", suffix.Data());
    AliAnalysisDataContainer *coutputTree = mgr->CreateContainer(name, TNtuple::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
    coutputTree->SetSpecialOutput();

    mgr->ConnectInput(hfResoTask, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(hfResoTask, 1, coutput);
    mgr->ConnectOutput(hfResoTask, 2, coutputCuts);
    mgr->ConnectOutput(hfResoTask, 3, coutputNorm);
    mgr->ConnectOutput(hfResoTask, 4, coutputTree);

    return hfResoTask;
}
