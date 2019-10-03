AliAnalysisTaskSE* AddTaskJetMultiPreparation(
  const char*    periodstr          = "LHC13b",
  TString        usedTracks         = "PicoTracks",
  const char*    usedMCParticles    = "MCParticlesSelected",
  const char*    usedClusters       = "CaloClusters",
  TString        outClusName        = "CaloClustersCorr",
  const Double_t hadCorr            = 2.0,
  const Double_t Eexcl              = 0.00,
  const Double_t phiMatch           = 0.03,
  const Double_t etaMatch           = 0.015,
  const Double_t minPtEt            = 0.15,
  const UInt_t   pSel               = AliVEvent::kAny,
  const Bool_t   trackclus          = kTRUE,
  const Bool_t   doHistos           = kFALSE,
  const Bool_t   makePicoTracks     = kTRUE,
  const Bool_t   makeTrigger        = kFALSE,
  const Bool_t   isEmcalTrain       = kFALSE,
  const Double_t trackEff           = 1.0,
  const Bool_t   doAODTrackProp     = kFALSE
)
{
    // Add task macros for all jet related helper tasks.
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        Error("AddTaskJetPreparation","No analysis manager found.");
        return 0;
    }

    AliVEventHandler *evhand = mgr->GetInputEventHandler();
    if (!evhand)
    {
        Error("AddTaskJetPreparation", "This task requires an input event handler");
        return NULL;
    }

    // Convert hadCorr & trackEff to TStrings to be used in naming conventions
    TString sTrackEff = "";
    TString sHadCorr = "";
    Int_t tEff = Int_t(trackEff*100);
    Int_t hC = Int_t(hadCorr*100);
    
    if (trackEff<0.1)
    {
        sTrackEff = Form("00%d",tEff);
    }
    else if (trackEff<1.0)
    {
        sTrackEff = Form("0%d",tEff);
    }
    else
    {
        sTrackEff = Form("%d",tEff);
    }
    
    if (hadCorr<0.1)
    {
        sHadCorr = Form("00%d",hC);
    }
    else if (hadCorr<1.0)
    {
        sHadCorr = Form("0%d",hC);
    }
    else
    {
        sHadCorr = Form("%d",hC);
    }
    
    // Establish names for output collections and tasks
    TString sPicoTracks = usedTracks + "_" + sTrackEff;
    TString sEmcalTriggers = "EmcalTriggers_" + sTrackEff;
    TString sEmcalTracks = "EmcalTracks_" + sTrackEff;
    TString sEmcalClusters = "EmcalClusters_" + sTrackEff;
    TString sCaloClustersCorr = outClusName + "_" + sTrackEff + "_" + sHadCorr;
    
    TString sESDTrackFilterName = "AliEmcalEsdTrackFilterTask";
    TString sAODTrackPropagatorName = "AliEmcalTrackPropagatorTaskAOD";
    TString sEmcalPicoTrackMakerName = "EmcalPicoTrackMaker_" + sTrackEff;
    TString sEmcalParticleMakerName = "EmcalParticleMaker_" + sTrackEff;
    TString sEmcalTriggerMakerName = "EmcalTriggerMaker_" + sTrackEff;
    TString sEmcalTriggerMakerSetupOutName = "EmcalTriggerMakerSetupOut_" + sTrackEff;
    TString sMCTrackSelectorName = "AliEmcalMCTrackSelector";
    
    // Set trackcuts according to period. Every period used should be definied here
    TString period(periodstr);
    TString clusterColName(usedClusters);
    TString particleColName(usedMCParticles);
    TString dType("ESD");
    
    if (!evhand->InheritsFrom("AliESDInputHandler"))
    {
        dType = "AOD";
    }

    if ((dType == "AOD") && (clusterColName == "CaloClusters"))
    {
        clusterColName = "caloClusters";
    }
    if ((dType == "ESD") && (clusterColName == "caloClusters"))
    {
        clusterColName = "CaloClusters";
    }

    if (0)
    {
        gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalTrackPropagator.C");
        AliEmcalTrackPropagatorTask *proptask = AddTaskEmcalTrackPropagator();
        proptask->SelectCollisionCandidates(pSel);
    }

    if (makePicoTracks && (dType == "ESD" || dType == "AOD") )
    {
        TString inputTracks = "AODFilterTracks";
        const Double_t edist = 440;
        if (dType == "ESD")
        {
            inputTracks = "ESDFilterTracks";
            TString trackCuts(Form("Hybrid_%s",period.Data()));
            
            // Hybrid tracks maker for ESD
            AliEmcalEsdTrackFilterTask *esdFilter = mgr->GetTask(sESDTrackFilterName.Data());
            if (!esdFilter)
            {
                gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalEsdTrackFilter.C");
                AliEmcalEsdTrackFilterTask *esdFilter = AddTaskEmcalEsdTrackFilter(inputTracks.Data(),trackCuts.Data());
                esdFilter->SetDoPropagation(kTRUE);
                esdFilter->SetDist(edist);
                esdFilter->SelectCollisionCandidates(pSel);
            }
        }
        else if (dType == "AOD")
        {
            // Track propagator to extend track to the EMCal surface
            AliEmcalAodTrackFilterTask *aodFilter = mgr->GetTask(sAODTrackPropagatorName.Data());
            if (!aodFilter)
            {
                gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalAodTrackFilter.C");
                AliEmcalAodTrackFilterTask *aodFilter = AddTaskEmcalAodTrackFilter(inputTracks.Data(),"tracks",period.Data());
                aodFilter->SetDist(edist);
                aodFilter->SelectCollisionCandidates(pSel);
                if (doAODTrackProp == kTRUE)
                {
                    aodFilter->SetDoPropagation(kTRUE);
                }
            }
        }

        // Produce PicoTracks
        AliEmcalPicoTrackMaker *pTrackTask = mgr->GetTask(sEmcalPicoTrackMakerName.Data());
        if (!pTrackTask)
        {
            gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
            pTrackTask = AddTaskEmcalPicoTrackMaker(sPicoTracks.Data(), inputTracks.Data(),0,1000,-10,10,-10,10,trackEff,sEmcalPicoTrackMakerName.Data());
            pTrackTask->SelectCollisionCandidates(pSel);
        }
    }

    // Trigger maker
    if (makeTrigger)
    {
        AliEmcalTriggerMaker *emcalTriggers = mgr->GetTask(sEmcalTriggerMakerName.Data());
        if (!emcalTriggers)
        {
            gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalTriggerMaker.C");
            emcalTriggers = AddTaskEmcalTriggerMaker(sEmcalTriggers.Data(),sEmcalTriggerMakerSetupOutName.Data(),0,0,sEmcalTriggerMakerName.Data(),0,0,0,0,0,0);
            emcalTriggers->SelectCollisionCandidates(pSel);
        }
    }

    // Produce particles used for hadronic correction, then match!
    AliEmcalParticleMaker *emcalParts = mgr->GetTask(sEmcalParticleMakerName.Data());
    if (!emcalParts)
    {
        gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalParticleMaker.C");
        emcalParts = AddTaskEmcalParticleMaker(sPicoTracks.Data(),clusterColName.Data(),sEmcalTracks.Data(),sEmcalClusters.Data(),sEmcalParticleMakerName.Data());
        emcalParts->SelectCollisionCandidates(pSel);

        // Relate tracks and clusters
        gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
        AliEmcalClusTrackMatcherTask *emcalClus = AddTaskEmcalClusTrackMatcher(sEmcalTracks.Data(),sEmcalClusters.Data(),0.1,doHistos);
        emcalClus->SetModifyObjs(kFALSE);
        emcalClus->SelectCollisionCandidates(pSel);
        if (isEmcalTrain)
        {
            RequestMemory(emcalClus,100*1024);
        }
    }

    // Make Corrected CaloClusters
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskHadCorr.C");
    AliHadCorrTask *hCorr = AddTaskHadCorr(sEmcalTracks.Data(),sEmcalClusters.Data(),sCaloClustersCorr.Data(),hadCorr,minPtEt,phiMatch,etaMatch,Eexcl,trackclus,doHistos);
    hCorr->SelectCollisionCandidates(pSel);
    if (isEmcalTrain)
    {
        if (doHistos)
        {
            RequestMemory(hCorr,500*1024);
        }
    }

    // Produce MC particles
    if(particleColName != "")
    {
        AliEmcalMCTrackSelector *mcPartTask = mgr->GetTask(sMCTrackSelectorName.Data());
        if (!mcPartTask)
        {
            gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskMCTrackSelector.C");
            AliEmcalMCTrackSelector *mcPartTask = AddTaskMCTrackSelector(particleColName.Data(), kFALSE, kFALSE);
            mcPartTask->SelectCollisionCandidates(pSel);
        }
    }

    // Return one task that represents the jet preparation on LEGO trains
    return emcalParts;
}
