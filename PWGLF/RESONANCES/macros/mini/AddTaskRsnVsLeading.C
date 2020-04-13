AliAnalysisTask *AddTaskRsnVsLeading(TString taskName = "phi", Bool_t isMC = kFALSE, Bool_t isPP = kTRUE, UInt_t triggerMask = AliVEvent::kAny, Double_t nSigmaPart1 = 3.0, Double_t nSigmaPart2 = -1, TString configName = "ConfigPhiLeading.C", TString path = "")
{

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        ::Error("AddTaskRsnVsLeading", "No analysis manager to connect to.");
        return 0;
    }

    AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask(taskName, isMC);

    // Trigger setting
    task->SelectCollisionCandidates(triggerMask); //AOD
    task->UseESDTriggerMask(triggerMask);

    if(isPP){
        task->UseMultiplicity("QUALITY");
    } else {
        task->UseCentrality("V0M");
    }

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
    TString macroArgs = TString::Format("(AliRsnMiniAnalysisTask *)%p,%d,%d,%f,%f", task, isMC, isPP, nSigmaPart1, nSigmaPart2);
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
