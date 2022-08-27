AliAnalysisTaskAODTrackPair *AddTaskAODTrackPair(
    UInt_t offlineTriggerMask = AliVEvent::kINT7,
    string spd_multi_correction_file_path =
        "/Users/syano_mbp2021/analysis/run2/Glueball/SPDMultiCorrection.root",
    string dimuon_ds_file_path =
        "/Users/syano_mbp2021/analysis/run2/Glueball/DownScale_Run2_CTP.root",
    float min_vtxz = -10, float max_vtxz = 10, int min_vtx_cont = 1,
    float min_pair_rap = -0.5, float max_pair_rap = 0.5, float alpha = 0.2,
    float pangle = 0.998, float v0Dca = 0.1, float trackDca = 1.0,
    float min_dlength = 5.0, float max_dlength = 100., string period = "LHC18c",
    string multi_method = "SPDTracklets", bool onPURej = true,
    bool onLBcut = true, bool onMuEtaCut = true, bool onMuThetaAbsCut = true,
    bool onMuMatchAptCut = true, bool onMuMatchLptCut = true,
    bool onMuMatchHptCut = false, bool onMuChi2Cut = true,
    bool onMuPdcaCut = true, bool isMC = false, bool isSelectEvt = true,
    int paircuttype = 1, double min_pairtrackptcut = 0.5,
    bool onMixingAnalysis = false, bool isMidTrackAnalysis = false,
    bool isPrimTrackAnalysis = false, bool isV0TrackAnalysis = true,
    float min_track_pt = 0.0, float max_track_pt = 999.,
    float min_track_eta = -0.8, float max_track_eta = 0.8,
    float min_track_p = 0.3, float max_track_p = 2.5,
    float min_pion_sigma_tpc = -5, float max_pion_sigma_tpc = 5,
    float min_pion_sigma_tof = -5, float max_pion_sigma_tof = 5,
    float min_kaon_sigma_tpc = -2, float max_kaon_sigma_tpc = 2,
    float min_kaon_sigma_tof = -2, float max_kaon_sigma_tof = 2,
    float min_proton_sigma_tpc = -2, float max_proton_sigma_tpc = 2,
    float min_proton_sigma_tof = -2, float max_proton_sigma_tof = 2,
    float findable = 0.8, string dcaxy = "0.0105+0.035/pow(x,1.1)",
    float dcaz = 2.0, float chi2tpc = 4., float chi2its = 36.,
    int nclusttpc = 70, int nclustits = 1, int pid1 = 211, int pid2 = 211) {

  AliMuonTrackCuts *fMuonTrackCuts =
      new AliMuonTrackCuts("StandardMuonTrackCuts", "StandardMuonTrackCuts");
  fMuonTrackCuts->SetIsMC(isMC);
  fMuonTrackCuts->SetAllowDefaultParams(true);

  int selectionMask = 0;
  if (onMuEtaCut) {
    selectionMask |= AliMuonTrackCuts::kMuEta;
  }
  if (onMuThetaAbsCut) {
    selectionMask |= AliMuonTrackCuts::kMuThetaAbs;
  }
  if (onMuMatchAptCut) {
    selectionMask |= AliMuonTrackCuts::kMuMatchApt;
  }
  if (onMuMatchLptCut) {
    selectionMask |= AliMuonTrackCuts::kMuMatchLpt;
  }
  if (onMuMatchHptCut) {
    selectionMask |= AliMuonTrackCuts::kMuMatchHpt;
  }
  if (onMuPdcaCut) {
    selectionMask |= AliMuonTrackCuts::kMuPdca;
  }
  if (onMuChi2Cut) {
    selectionMask |= AliMuonTrackCuts::kMuTrackChiSquare;
  }
  fMuonTrackCuts->SetFilterMask(selectionMask);

  TFile *input = TFile::Open(dimuon_ds_file_path.c_str());
  TFile *input2 = TFile::Open(spd_multi_correction_file_path.c_str());

  AliAnalysisTaskAODTrackPairUtils *utils =
      new AliAnalysisTaskAODTrackPairUtils();
  utils->setMC(isMC);
  utils->setEvtSelection(isSelectEvt);
  utils->setDownScalingHist(input);
  utils->setSPDTrkCorrHist(input2, period);
  utils->setVertexCut(min_vtxz, max_vtxz, min_vtx_cont);
  utils->setPairRapidityCut(min_pair_rap, max_pair_rap);
  utils->setV0SelectCuts(alpha, pangle, v0Dca, trackDca, min_dlength,
                         max_dlength);
  utils->setPileupRejectionCut(onPURej);
  utils->setLocalBoardCut(onLBcut);
  utils->setMultiEstimateMethod(multi_method);
  utils->setMuonTrackCut(fMuonTrackCuts);
  utils->setPairKinematicCut(paircuttype, min_pairtrackptcut);
  utils->setMidTrackAna(isMidTrackAnalysis);
  utils->setPeriod(period);
  utils->setTrackKinematicCut(min_track_pt, max_track_pt, min_track_eta,
                              max_track_eta, min_track_p, max_track_p);
  utils->setPionSelectSigmaTPC(min_pion_sigma_tpc, max_pion_sigma_tpc);
  utils->setKaonSelectSigmaTPC(min_kaon_sigma_tpc, max_kaon_sigma_tpc);
  utils->setProtonSelectSigmaTPC(min_proton_sigma_tpc, max_proton_sigma_tpc);
  utils->setPionSelectSigmaTOF(min_pion_sigma_tof, max_pion_sigma_tof);
  utils->setKaonSelectSigmaTOF(min_kaon_sigma_tof, max_kaon_sigma_tof);
  utils->setProtonSelectSigmaTOF(min_proton_sigma_tof, max_proton_sigma_tof);
  utils->setTrackQualities(findable, dcaxy.c_str(), dcaz, chi2tpc, chi2its,
                           nclusttpc, nclustits);
  utils->setPairTargetPIDs(pid1, pid2);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskAODMuonEventSelection",
            "No analysis manager to connect to");
    return NULL;
  }
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskAODMuonEventSelection",
            "This task requires an input event handler");
    return NULL;
  }

  AliAnalysisTaskAODTrackPair *task = new AliAnalysisTaskAODTrackPair("dimuon");
  if (isSelectEvt) {
    task->SelectCollisionCandidates(offlineTriggerMask);
  }
  task->setUtils(utils);
  task->setEvtMixingTrackDepth(100);
  task->setEvtMixingPoolSize(100);
  task->setEvtMixingReadyFraction(0.1);
  task->setEvtMixingPoolVtxZ(true);
  task->setEvtMixingPoolCent(true);
  task->setEvtMixingPoolPsi(true);
  task->setMixingAnalysis(onMixingAnalysis);
  task->setMC(isMC);
  task->setMidTrackAna(isMidTrackAnalysis);
  task->setPrimTrackAna(isPrimTrackAnalysis);
  task->setV0TrackAna(isV0TrackAnalysis);
  mgr->AddTask(task);

  cout << "min_vtxz=" << min_vtxz << endl;
  cout << "max_vtxz=" << max_vtxz << endl;
  cout << "min_vtx_cont=" << min_vtx_cont << endl;
  cout << "min_pair_rap=" << min_pair_rap << endl;
  cout << "max_pair_rap=" << max_pair_rap << endl;
  cout << "alpha=" << alpha << endl;
  cout << "pangle=" << pangle << endl;
  cout << "v0Dca=" << v0Dca << endl;
  cout << "trackDca=" << trackDca << endl;
  cout << "min_dlength=" << min_dlength << endl;
  cout << "max_dlength=" << max_dlength << endl;
  cout << "period=" << period << endl;
  cout << "multi_method=" << multi_method << endl;
  cout << "onPURej=" << onPURej << endl;
  cout << "onLBcut=" << onLBcut << endl;
  cout << "onMuEtaCut=" << onMuEtaCut << endl;
  cout << "onMuThetaAbsCut=" << onMuThetaAbsCut << endl;
  cout << "onMuMatchAptCut=" << onMuMatchAptCut << endl;
  cout << "onMuMatchLptCut=" << onMuMatchLptCut << endl;
  cout << "onMuMatchHptCut=" << onMuMatchHptCut << endl;
  cout << "onMuChi2Cut=" << onMuChi2Cut << endl;
  cout << "onMuPdcaCut=" << onMuPdcaCut << endl;
  cout << "isMC=" << isMC << endl;
  cout << "isSelectEvt=" << isSelectEvt << endl;
  cout << "paircuttype=" << paircuttype << endl;
  cout << "min_pairtrackptcut=" << min_pairtrackptcut << endl;
  cout << "onMixingAnalysis=" << onMixingAnalysis << endl;
  cout << "isMidTrackAnalysis=" << isMidTrackAnalysis << endl;
  cout << "isPrimTrackAnalysis=" << isPrimTrackAnalysis << endl;
  cout << "isV0TrackAnalysis=" << isV0TrackAnalysis << endl;
  cout << "min_track_pt=" << min_track_pt << endl;
  cout << "max_track_pt=" << max_track_pt << endl;
  cout << "min_track_eta=" << min_track_eta << endl;
  cout << "max_track_eta=" << max_track_eta << endl;
  cout << "min_track_p=" << min_track_p << endl;
  cout << "max_track_p=" << max_track_p << endl;
  cout << "min_pion_sigma_tpc=" << min_pion_sigma_tpc << endl;
  cout << "max_pion_sigma_tpc=" << max_pion_sigma_tpc << endl;
  cout << "min_pion_sigma_tpc=" << min_pion_sigma_tof << endl;
  cout << "max_pion_sigma_tpc=" << max_pion_sigma_tof << endl;
  cout << "min_kaon_sigma_tpc=" << min_kaon_sigma_tpc << endl;
  cout << "max_kaon_sigma_tpc=" << max_kaon_sigma_tpc << endl;
  cout << "min_kaon_sigma_tpc=" << min_kaon_sigma_tof << endl;
  cout << "max_kaon_sigma_tpc=" << max_kaon_sigma_tof << endl;
  cout << "min_proton_sigma_tpc=" << min_proton_sigma_tpc << endl;
  cout << "max_proton_sigma_tpc=" << max_proton_sigma_tpc << endl;
  cout << "min_proton_sigma_tpc=" << min_proton_sigma_tof << endl;
  cout << "max_proton_sigma_tpc=" << max_proton_sigma_tof << endl;
  cout << "findable=" << findable << endl;
  cout << "dcaxy=" << dcaxy.c_str() << endl;
  cout << "dcaz=" << dcaz << endl;
  cout << "chi2tpc=" << chi2tpc << endl;
  cout << "chi2its=" << chi2its << endl;
  cout << "nclusttpc=" << nclusttpc << endl;
  cout << "nclustits=" << nclustits << endl;
  cout << "pid1=" << pid1 << endl;
  cout << "pid2=" << pid2 << endl;

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 =
      mgr->CreateContainer("dimuon", TList::Class(),
                           AliAnalysisManager::kOutputContainer, "Dimuon.root");

  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);

  return task;
}
