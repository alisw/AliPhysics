void AddTaskJetPreparationWrapper(const char *configstring){
  AliEMCALConfiguration defaultConfig("defaultConfig");
  defaultConfig.AddParam("periodstr", new AliJSONString("LHC11h"));
  defaultConfig.AddParam("pTracksName", new AliJSONString("PicoTracks"));
  defaultConfig.AddParam("usedMCParticles", new AliJSONString("MCParticlesSelected"));
  defaultConfig.AddParam("usedClusters", new AliJSONString("CaloClusters"));
  defaultConfig.AddParam("outClusName", new AliJSONString("CaloClustersCorr"));
  defaultConfig.AddParam("hadcorr", new AliJSONDouble(2.0));
  defaultConfig.AddParam("Eexcl", new AliJSONDouble(0.00));
  defaultConfig.AddParam("phiMatch", new AliJSONDouble(0.03));
  defaultConfig.AddParam("etaMatch", new AliJSONDouble(0.015));
  defaultConfig.AddParam("minPtEt", new AliJSONDouble(0.15));
  defaultConfig.AddParam("pSel", new AliJSONInt(AliVEvent::kAny));
  defaultConfig.AddParam("trackclus", new AliJSONBool(kTRUE));
  defaultConfig.AddParam("doHistos", new AliJSONBool(kFALSE));
  defaultConfig.AddParam("makePicoTracks", new AliJSONBool(kTRUE));
  defaultConfig.AddParam("makeTrigger", new AliJSONBool(kTRUE));
  defaultConfig.AddParam("isEMCALTrain", new AliJSONBool(kFALSE));
  defaultConfig.AddParam("trackeff", new AliJSONDouble(1.0));
  defaultConfig.AddParam("doAODTrackProp", new AliJSONBool(kTRUE));
  defaultConfig.AddParam("modifyMatchObjects", new AliJSONBool(kTRUE));
  defaultConfig.AddParam("doTriggerQA", new AliJSONBool(kTRUE));
  AliEMCALConfiguration userConfig("userConfig");
  userConfig.Build(configstring);
  AliEMCALConfigMatcher combinedConfiguration(&defaultConfig, &userConfig);

  TString periodstr = (static_cast<AliJSONString *>(combinedConfig.GetValue("periodstr")))->GetValue(),
          pTracksName = (static_cast<AliJSONString *>(combinedConfig.GetValue("pTracksName")))->GetValue(),
          usedMCParticles = (static_cast<AliJSONString *>(combinedConfig.GetValue("usedMCParticles")))->GetValue(),
          usedClusters = (static_cast<AliJSONString *>(combinedConfig.GetValue("usedClusters")))->GetValue(),
          outClusName = (static_cast<AliJSONString *>(combinedConfig.GetValue("outClusName")))->GetValue();
  Double_t hadcorr = (static_cast<AliJSONDouble *>(combinedConfig.GetValue("hadcorr")))->GetValue(),
          Eexcl = (static_cast<AliJSONDouble *>(combinedConfig.GetValue("Eexcl")))->GetValue(),
          phiMatch = (static_cast<AliJSONDouble *>(combinedConfig.GetValue("phiMatch")))->GetValue(),
          etaMatch = (static_cast<AliJSONDouble *>(combinedConfig.GetValue("etaMatch")))->GetValue(),
          minPtEt = (static_cast<AliJSONDouble *>(combinedConfig.GetValue("minPtEt")))->GetValue(),
          trackeff = (static_cast<AliJSONDouble *>(combinedConfig.GetValue("trackeff")))->GetValue();
  UInt_t pSel = static_cast<UInt_t>((static_cast<AliJSONInt *>(combinedConfig.GetValue("pSel")))->GetValue());
  Bool_t trackclus = (static_cast<AliJSONBool *>(combinedConfig.GetValue("trackclus")))->GetValue(),
         doHistos = (static_cast<AliJSONBool *>(combinedConfig.GetValue("doHistos")))->GetValue(),
         makePicoTracks = (static_cast<AliJSONBool *>(combinedConfig.GetValue("makePicoTracks")))->GetValue(),
         makeTrigger = (static_cast<AliJSONBool *>(combinedConfig.GetValue("makeTrigger")))->GetValue(),
         isEMCALTrain = (static_cast<AliJSONBool *>(combinedConfig.GetValue("isEMCALTrain")))->GetValue(),
         doAODTrackProp = (static_cast<AliJSONBool *>(combinedConfig.GetValue("doAODTrackProp")))->GetValue(),
         modifyMatchObjects = (static_cast<AliJSONBool *>(combinedConfig.GetValue("modifyMatchObjects")))->GetValue(),
         doTriggerQA = (static_cast<AliJSONBool *>(combinedConfig.GetValue("doTriggerQA")))->GetValue();

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCLAJetTasks/macros/AddTaskJetPreparation.C");
  AddTaskJetPreparation(periodstr.Data(), pTracksName.Data(), usedMCParticles.Data(), usedClusters.Data(), outClusName.Data(),
      hadcorr, Eexcl, phiMatch, etaMatch, minPtEt, pSel, trackclus, doHistos, makePicoTracks, makeTrigger, isEMCALTrain, trackeff,
      doAODTrackProp, modifyMatchObjects, doTriggerQA);
}
