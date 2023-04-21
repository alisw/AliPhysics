AliJHFTagging *AddTaskJHFTagging(
    const char *ntracks = "tracks",
    const char *njets = "Jets",
    const char *nrho = "",
    double jetradius = 0.4,
    bool isMC = false,
    const char *type = "TPC",
    const char *taskname = "AliJHFTagging",
    const char *njetsMC = "Jets",
    const char *nrhoMC = "",
    const char *suffix = "")
{

    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        ::Error("AliJHFTagging", "No analysis manager to connect to.");
        return NULL;
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler())
    {
        ::Error("AliJHFTagging", "This task requires an input event handler");
        return NULL;
    }

    TString name(taskname);

    TString combinedName;
    combinedName.Form("%s%s", name.Data(), suffix);

    AliJHFTagging *jetTask = new AliJHFTagging(combinedName);

    AliTrackContainer *trackCont = jetTask->AddTrackContainer(ntracks);

    TString strType(type);

    AliJetContainer *jetCont = jetTask->AddJetContainer(njets, strType, jetradius);

    if (jetCont)
    {
        jetCont->SetRhoName(nrho);
        jetCont->ConnectParticleContainer(trackCont);
        jetCont->SetJetEtaLimits(-0.5, 0.5);
        jetCont->SetJetPtCut(0.0);
        jetCont->SetMaxTrackPt(1000);
        jetCont->SetPercAreaCut(0.6);
    }

    if (isMC)
    {
        AliJetContainer *jetContMC = jetTask->AddJetContainer(njetsMC, strType, jetradius);

        if (jetContMC)
        {

            jetContMC->SetRhoName(nrhoMC);
            jetContMC->SetIsParticleLevel(true);
            jetContMC->SetJetEtaLimits(-0.5, 0.5);
            jetContMC->SetJetPtCut(0.0);
            jetContMC->SetMaxTrackPt(1000);
            jetContMC->SetPercAreaCut(0.6);
        }
    }
    //-------------------------------------------------------
    //  Configure analysis task
    //-------------------------------------------------------
    jetTask->SetIsPythia(isMC);

    //-------------------------------------------------------
    // Final settings, pass to manager and set the containers
    //-------------------------------------------------------
    mgr->AddTask(jetTask);

    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();

    TString contname(combinedName);
    contname += "_histos";
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
                                                              TList::Class(), AliAnalysisManager::kOutputContainer,
                                                              Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectInput(jetTask, 0, cinput1);
    mgr->ConnectOutput(jetTask, 1, coutput1);

    return jetTask;
}
