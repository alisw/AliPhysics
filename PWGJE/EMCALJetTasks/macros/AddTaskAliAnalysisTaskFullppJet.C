AliAnalysisTaskFullppJet *AddTaskAliAnalysisTaskFullppJet(const char *name = "Baseline",
							   const char *period = "lhc11a",
							   const Bool_t isMC = kFALSE,
							   const Bool_t IsPhySelForMC = kFALSE,
							   const Bool_t offlineTrig = kFALSE,
							   const Double_t minTrkPt = 0.15,
							   const Double_t minClsEt = 0.30,
							   const Bool_t hc = kTRUE,
							   const Double_t fraction = 1)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr)
    {
      AliError("No analysise manager is availabe !");
      return NULL;
    }

  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/macros/CreateTrackCutsPWGJE.C");
  AliESDtrackCuts *esdTrackCuts = 0x0;
  AliESDtrackCuts *hybridTrackCuts1 = 0x0;
  AliESDtrackCuts *hybridTrackCuts2 = 0x0;
  printf("\n===== Use Hybrid track cuts =====\n");
  esdTrackCuts = CreateTrackCutsPWGJE(10001006);
  hybridTrackCuts1 = CreateTrackCutsPWGJE(1006);
  hybridTrackCuts2 = CreateTrackCutsPWGJE(10041006);

  AliAnalysisTaskFullppJet *jetTask = new AliAnalysisTaskFullppJet(Form("ppJet_%s_%s",period,name));
  jetTask->SetNonStdBranch(name);
  jetTask->SetRunPeriod(period);
  jetTask->SetAnaType(1);
  jetTask->SetCheckTriggerMask(kFALSE);
  jetTask->SetMCAna(isMC);
  jetTask->SetPhySelForMC(IsPhySelForMC);
  jetTask->SetRejectSPDPileup(kTRUE);
  jetTask->SetZvtx(10);
  jetTask->SetTrackCutsType(AliAnalysisTaskFullppJet::kHybrid);      
  jetTask->SetEsdTrackCuts(esdTrackCuts);
  jetTask->SetHybridTrackCuts1(hybridTrackCuts1);
  jetTask->SetHybridTrackCuts2(hybridTrackCuts2);
  jetTask->SetOfflineTrigger(offlineTrig);
  jetTask->SetEtaMax(1);
  jetTask->SetdEdxRange(75,95);
  jetTask->SetEoverPRange(0.8,1.2);
  jetTask->SetPtRange(minTrkPt,1e4,minTrkPt,1e4);
  jetTask->SetEtRange(minClsEt,1e4,minClsEt,1e4);
  jetTask->SetRejectExoticCluster(kTRUE);
  jetTask->SetRemoveProblematicSM4(kTRUE);  
  jetTask->SetStudySubEInHC(kFALSE);
  jetTask->SetJetNEFCut(0.02,0.98);
  jetTask->SetRejectElectron(hc);
  jetTask->SetCorrectHadron(hc);
  jetTask->SetHCFraction(fraction);
  jetTask->SetRadius("0.4 0.2 0.3");
  jetTask->SetCheckTrkEffCorr(kFALSE);

  mgr->AddTask(jetTask);
  TString outfileName = "ppJetOutput.root";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("JetOutputList_%s",name), TList::Class(), AliAnalysisManager::kOutputContainer,outfileName.Data());
  mgr->ConnectInput(jetTask,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(jetTask, 1, coutput1);

  return jetTask;
}
