AliAnalysisTask *AddTaskRsnVsLeading(TString taskName = "phi", Bool_t isMC = kFALSE, UInt_t triggerMask = AliVEvent::kMB, Double_t nSigmaKaon = 3.0, TString configName = "ConfigPhiLeading.C", TString path = "")
{

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        ::Error("AddTaskRsnVsLeading", "No analysis manager to connect to.");
        return 0;
    }

    AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName, isMC);

    // Trigger setting
    //UInt_t triggerMask = AliVEvent::kINT7; // for pp
   // triggerMask=AliVEvent::kAny; // else
    task->SelectCollisionCandidates(triggerMask); //AOD

    // Mixing setings
    Int_t nmix = 10;
    Float_t maxDiffVzMix = 1.;
    Float_t maxDiffMultMix = 10.;
    task->UseContinuousMix();
    //task->UseBinnedMix();
    task->SetNMix(nmix);
    task->SetMaxDiffVz(maxDiffVzMix);
    task->SetMaxDiffMult(maxDiffMultMix);

    Double_t minYlab = -0.5;
    Double_t maxYlab = 0.5;
    // Pair cuts (common to all resonances)
    AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
    cutY->SetRangeD(minYlab, maxYlab);

    AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
    cutsPair->AddCut(cutY);
    cutsPair->SetCutScheme(cutY->GetName());

    // We will add RSN config file
    TString macroArgs = TString::Format("(AliRsnMiniAnalysisTask *)%p,%d,%f", task, isMC, nSigmaKaon);
    TMacro cfg(gSystem->ExpandPathName(TString::Format("%s%s", path.Data(), configName.Data()).Data()));
    Long_t rc = reinterpret_cast<Long_t>(cfg.Exec(macroArgs.Data()));
    if (!rc)
        return 0;

    mgr->AddTask(task);

    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    AliAnalysisDataContainer *output = mgr->CreateContainer(TString::Format("%sOut", taskName.Data()).Data(),
                                                            TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            outputFileName);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, output);

    return task;
}
