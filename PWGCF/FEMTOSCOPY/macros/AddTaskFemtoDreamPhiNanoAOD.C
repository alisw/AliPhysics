AliAnalysisTaskSE *AddTaskFemtoDreamPhiNanoAOD(bool isMC = false,
                                               TString CentEst = "kInt7",
                                               const char *cutVariation = "0") {
  TString suffix = TString::Format("%s", cutVariation);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    printf("No analysis manager to connect to!\n");
    return nullptr;
  }
  if (!mgr->GetInputEventHandler()) {
    printf("This task requires an input event handler!\n");
    return nullptr;
  }

  AliFemtoDreamEventCuts *evtCuts = AliFemtoDreamEventCuts::StandardCutsRun2();
  evtCuts->CleanUpMult(false, false, false, true);
  evtCuts->SetSphericityCuts(0.7, 1);

  //    if (suffix == "1") {
  //        evtCuts->SetSpherocityCuts(0.5,1);
  //    }

  //    if (suffix == "2") {
  //        evtCuts->SetSpherocityCuts(0.6,1);
  //    }

  //    if (suffix == "3") {
  //        evtCuts->SetSpherocityCuts(0.7,1);
  //    }

  //    if (suffix == "4") {
  //        evtCuts->SetSpherocityCuts(0.8,1);
  //    }

  //    if (suffix == "5") {
  //        evtCuts->SetSpherocityCuts(0.9,1);
  //    }

  // Proton cuts
  const float ProtonPtlow = 0.45;
  const float ProtonPtup = 0.55;
  const float ProtonEtaLow = 0.78;
  const float ProtonEtaUp = 0.82;
  const float ProtonNsigmaLow = 2.7;
  const float ProtonNsigmaUp = 3.3;
  const float ProtonNClsLow = 70;
  const float ProtonNClsUp = 90;

  AliFemtoDreamTrackCuts *TrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  TrackCuts->SetCutCharge(1);

  AliFemtoDreamTrackCuts *AntiTrackCuts =
      AliFemtoDreamTrackCuts::PrimProtonCuts(isMC, true, false, false);
  AntiTrackCuts->SetCutCharge(-1);

  if (suffix != "0" && suffix != "999") {
    TrackCuts->SetMinimalBooking(true);
    AntiTrackCuts->SetMinimalBooking(true);
  }
  //  if (suffix == "1") {
  //    TrackCuts->SetPtRange(ProtonPtlow, 4.05);
  //    AntiTrackCuts->SetPtRange(ProtonPtlow, 4.05);
  //    TrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
  //    AntiTrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
  //  } else if (suffix == "2") {
  //    TrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  //    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  //  } else if (suffix == "3") {
  //    TrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
  //    AntiTrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
  //    TrackCuts->SetNClsTPC(ProtonNClsUp);
  //    AntiTrackCuts->SetNClsTPC(ProtonNClsUp);
  //  } else if (suffix == "4") {
  //    TrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaUp);
  //    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaUp);
  //  } else if (suffix == "5") {
  //    TrackCuts->SetEtaRange(-ProtonEtaLow, ProtonEtaLow);
  //    AntiTrackCuts->SetEtaRange(-ProtonEtaLow, ProtonEtaLow);
  //    TrackCuts->SetNClsTPC(ProtonNClsLow);
  //    AntiTrackCuts->SetNClsTPC(ProtonNClsLow);
  //  } else if (suffix == "6") {
  //    TrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaUp);
  //    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaUp);
  //  } else if (suffix == "7") {
  //    TrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
  //    AntiTrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
  //    TrackCuts->SetNClsTPC(ProtonNClsLow);
  //    AntiTrackCuts->SetNClsTPC(ProtonNClsLow);
  //  } else if (suffix == "8") {
  //    TrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  //    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  //  } else if (suffix == "9") {
  //    TrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    TrackCuts->SetNClsTPC(ProtonNClsUp);
  //    AntiTrackCuts->SetNClsTPC(ProtonNClsUp);
  //  } else if (suffix == "10") {
  //    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  //    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  //    TrackCuts->SetNClsTPC(ProtonNClsLow);
  //    AntiTrackCuts->SetNClsTPC(ProtonNClsLow);
  //  } else if (suffix == "11") {
  //    TrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  //    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  //  } else if (suffix == "12") {
  //    TrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    TrackCuts->SetNClsTPC(ProtonNClsUp);
  //    AntiTrackCuts->SetNClsTPC(ProtonNClsUp);
  //  } else if (suffix == "13") {
  //    TrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
  //    AntiTrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
  //    TrackCuts->SetPtRange(ProtonPtlow, 4.05);
  //    AntiTrackCuts->SetPtRange(ProtonPtlow, 4.05);
  //  } else if (suffix == "14") {
  //    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  //    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  //    TrackCuts->SetEtaRange(-ProtonEtaLow, ProtonEtaLow);
  //    AntiTrackCuts->SetEtaRange(-ProtonEtaLow, ProtonEtaLow);
  //  } else if (suffix == "15") {
  //    TrackCuts->SetNClsTPC(ProtonNClsUp);
  //    AntiTrackCuts->SetNClsTPC(ProtonNClsUp);
  //    TrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
  //  } else if (suffix == "16") {
  //    TrackCuts->SetEtaRange(-ProtonEtaLow, ProtonEtaLow);
  //    AntiTrackCuts->SetEtaRange(-ProtonEtaLow, ProtonEtaLow);
  //    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  //    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  //  } else if (suffix == "17") {
  //    TrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    TrackCuts->SetNClsTPC(ProtonNClsLow);
  //    AntiTrackCuts->SetNClsTPC(ProtonNClsLow);
  //  } else if (suffix == "18") {
  //    TrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaUp);
  //    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaUp);
  //  } else if (suffix == "19") {
  //    TrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  //    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  //  } else if (suffix == "20") {
  //    TrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  //    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  //  } else if (suffix == "21") {
  //    TrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    TrackCuts->SetNClsTPC(ProtonNClsLow);
  //    AntiTrackCuts->SetNClsTPC(ProtonNClsLow);
  //  } else if (suffix == "22") {
  //    TrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    AntiTrackCuts->SetPtRange(ProtonPtup, 4.05);
  //    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  //    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaLow);
  //  } else if (suffix == "23") {
  //    TrackCuts->SetPtRange(ProtonPtlow, 4.05);
  //    AntiTrackCuts->SetPtRange(ProtonPtlow, 4.05);
  //    TrackCuts->SetNClsTPC(ProtonNClsUp);
  //    AntiTrackCuts->SetNClsTPC(ProtonNClsUp);
  //  } else if (suffix == "24") {
  //    TrackCuts->SetEtaRange(-ProtonEtaLow, ProtonEtaLow);
  //    AntiTrackCuts->SetEtaRange(-ProtonEtaLow, ProtonEtaLow);
  //    TrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaUp);
  //    AntiTrackCuts->SetPID(AliPID::kProton, 0.75, ProtonNsigmaUp);
  //  } else if (suffix == "25") {
  //    TrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
  //    AntiTrackCuts->SetEtaRange(-ProtonEtaUp, ProtonEtaUp);
  //    TrackCuts->SetNClsTPC(ProtonNClsLow);
  //    AntiTrackCuts->SetNClsTPC(ProtonNClsLow);
  //  }

  // Kaon cuts
  const float KaonPtlow = 0.1;
  const float KaonPtup = 0.2;
  const float KaonEtaLow = 0.78;
  const float KaonEtaUp = 0.82;
  const float KaonNsigmaLow = 4.5;
  const float KaonNsigmaUp = 5.5;
  const float KaonNClsLow = 70;
  const float KaonNClsUp = 90;
  const float SpheriLow = 0.68;
  const float SpheriUp = 0.72;

  AliFemtoDreamTrackCuts *TrackPosKaonCuts =
      AliFemtoDreamTrackCuts::PrimKaonCuts(isMC, true, false, false);
  TrackPosKaonCuts->SetCutCharge(1);
  TrackPosKaonCuts->SetFilterBit(128);

  AliFemtoDreamTrackCuts *TrackNegKaonCuts =
      AliFemtoDreamTrackCuts::PrimKaonCuts(isMC, true, false, false);
  TrackNegKaonCuts->SetCutCharge(-1);
  TrackNegKaonCuts->SetFilterBit(128);

  TrackPosKaonCuts->SetDCAVtxZ(0.4);
  TrackNegKaonCuts->SetDCAVtxZ(0.4);
  TrackPosKaonCuts->SetDCAVtxXY(0.8);
  TrackNegKaonCuts->SetDCAVtxXY(0.8);

  if (suffix != "0" && suffix != "999") {
    TrackPosKaonCuts->SetMinimalBooking(true);
    TrackNegKaonCuts->SetMinimalBooking(true);
  }
  //  if (suffix == "1") {
  //    TrackPosKaonCuts->SetPtRange(KaonPtlow, 999);
  //    TrackNegKaonCuts->SetPtRange(KaonPtlow, 999);
  //    TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
  //    TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
  //  } else if (suffix == "2") {
  //    evtCuts->SetSphericityCuts(SpheriLow,1);
  //    TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
  //    TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);
  //  } else if (suffix == "3") {
  //    evtCuts->SetSphericityCuts(SpheriUp,1);
  //    TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
  //    TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
  //  } else if (suffix == "4") {
  //    evtCuts->SetSphericityCuts(SpheriLow,1);
  //    TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
  //     TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
  //  } else if (suffix == "5") {
  //    evtCuts->SetSphericityCuts(SpheriUp,1);
  //    TrackPosKaonCuts->SetPtRange(KaonPtlow, 999);
  //    TrackNegKaonCuts->SetPtRange(KaonPtlow, 999);
  //  } else if (suffix == "6") {
  //    TrackPosKaonCuts->SetPtRange(KaonPtlow, 999);
  //    TrackNegKaonCuts->SetPtRange(KaonPtlow, 999);
  //    TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
  //     TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
  //  } else if (suffix == "7") {
  //    evtCuts->SetSphericityCuts(SpheriLow,1);
  //    TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
  //    TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
  //  } else if (suffix == "8") {
  //    TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
  //    TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
  //    TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
  //    TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);
  //  } else if (suffix == "9") {
  //    TrackPosKaonCuts->SetPtRange(KaonPtup, 999);
  //    TrackNegKaonCuts->SetPtRange(KaonPtup, 999);
  //    TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
  //    TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
  //  } else if (suffix == "10") {
  //    TrackPosKaonCuts->SetPtRange(KaonPtup, 999);
  //    TrackNegKaonCuts->SetPtRange(KaonPtup, 999);
  //    TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
  //    TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
  //  } else if (suffix == "11") {
  //    TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
  //    TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
  //    TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
  //    TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
  //  } else if (suffix == "12") {
  //    TrackPosKaonCuts->SetPtRange(KaonPtlow, 999);
  //    TrackNegKaonCuts->SetPtRange(KaonPtlow, 999);
  //    TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
  //    TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
  //  } else if (suffix == "13") {
  //    TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
  //    TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
  //    TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
  //    TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);
  //  } else if (suffix == "14") {
  //    TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
  //    TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);
  //    TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
  //    TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaLow);
  //  } else if (suffix == "15") {
  //    TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
  //    TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
  //    TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
  //    TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);
  //  } else if (suffix == "16") {
  //    TrackPosKaonCuts->SetPtRange(KaonPtup, 999);
  //    TrackNegKaonCuts->SetPtRange(KaonPtup, 999);
  //    TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
  //    TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
  //  } else if (suffix == "17") {
  //    TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
  //    TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
  //    TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
  //    TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
  //    TrackPosKaonCuts->SetPtRange(KaonPtlow, 999);
  //    TrackNegKaonCuts->SetPtRange(KaonPtlow, 999);
  //  } else if (suffix == "18") {
  //    TrackPosKaonCuts->SetPtRange(KaonPtup, 999);
  //    TrackNegKaonCuts->SetPtRange(KaonPtup, 999);
  //    TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
  //    TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);
  //  } else if (suffix == "19") {
  //    TrackPosKaonCuts->SetPtRange(KaonPtup, 999);
  //    TrackNegKaonCuts->SetPtRange(KaonPtup, 999);
  //    TrackPosKaonCuts->SetNClsTPC(KaonNClsUp);
  //    TrackNegKaonCuts->SetNClsTPC(KaonNClsUp);
  //  } else if (suffix == "20") {
  //    TrackPosKaonCuts->SetPtRange(KaonPtlow, 999);
  //    TrackNegKaonCuts->SetPtRange(KaonPtlow, 999);
  //    TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
  //    TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
  //    TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
  //    TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
  //  } else if (suffix == "21") {
  //    TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
  //    TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
  //    TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
  //    TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);
  //  } else if (suffix == "22") {
  //    TrackPosKaonCuts->SetPtRange(KaonPtlow, 999);
  //    TrackNegKaonCuts->SetPtRange(KaonPtlow, 999);
  //    TrackPosKaonCuts->SetNClsTPC(KaonNClsLow);
  //    TrackNegKaonCuts->SetNClsTPC(KaonNClsLow);
  //  } else if (suffix == "23") {
  //    evtCuts->SetSphericityCuts(SpheriUp,1);
  //    TrackPosKaonCuts->SetPtRange(KaonPtup, 999);
  //    TrackNegKaonCuts->SetPtRange(KaonPtup, 999);
  //    TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
  //    TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
  //  } else if (suffix == "24") {
  //    TrackPosKaonCuts->SetPtRange(KaonPtup, 999);
  //    TrackNegKaonCuts->SetPtRange(KaonPtup, 999);
  //    TrackPosKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
  //    TrackNegKaonCuts->SetEtaRange(-KaonEtaLow, KaonEtaLow);
  //  } else if (suffix == "25") {
  //    TrackPosKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
  //    TrackNegKaonCuts->SetEtaRange(-KaonEtaUp, KaonEtaUp);
  //    TrackPosKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
  //    TrackNegKaonCuts->SetPID(AliPID::kKaon, 0.4, KaonNsigmaUp);
  //  }

  AliFemtoDreamv0Cuts *TrackCutsPhi = new AliFemtoDreamv0Cuts();
  TrackCutsPhi->SetIsMonteCarlo(isMC);
  TrackCutsPhi->SetAxisInvMassPlots(400, 0.95, 2);
  TrackCutsPhi->SetCutInvMass(0.008);
  AliFemtoDreamTrackCuts *dummyCutsPos = new AliFemtoDreamTrackCuts();
  dummyCutsPos->SetIsMonteCarlo(isMC);
  AliFemtoDreamTrackCuts *dummyCutsNeg = new AliFemtoDreamTrackCuts();
  dummyCutsNeg->SetIsMonteCarlo(isMC);
  TrackCutsPhi->SetPosDaugterTrackCuts(dummyCutsPos);
  TrackCutsPhi->SetNegDaugterTrackCuts(dummyCutsNeg);
  TrackCutsPhi->SetPDGCodePosDaug(321);
  TrackCutsPhi->SetPDGCodeNegDaug(321);
  TrackCutsPhi->SetPDGCodev0(333);
  double Phimass = TDatabasePDG::Instance()->GetParticle(333)->Mass();

  //  if (suffix != "0") {
  //    TrackCutsPhi->SetMinimalBooking(true);
  //  }

    if (suffix == "1") {
      TrackCutsPhi->SetCutWindow(0.987, 1.011);
    }
    if (suffix == "2") {
      TrackCutsPhi->SetCutWindow(1.027, 1.1);
    }
    if (suffix == "3") {
      evtCuts->SetSphericityCuts(0, 1);
    }
    if (suffix == "4") {
      TrackCutsPhi->SetCutWindow(0.987, 1.011);
      evtCuts->SetSphericityCuts(0, 1);
    }
    if (suffix == "5") {
      TrackCutsPhi->SetCutWindow(1.027, 1.1);
      evtCuts->SetSphericityCuts(0, 1);
    }
    if (suffix == "7") {
      TrackCutsPhi->SetCutWindow(0.987, 1.011);
    }
    if (suffix == "8") {
      TrackCutsPhi->SetCutWindow(1.027, 1.1);
    }
    if (suffix == "9") {
      evtCuts->SetSphericityCuts(0, 1);
    }
    if (suffix == "10") {
      TrackCutsPhi->SetCutWindow(0.987, 1.011);
      evtCuts->SetSphericityCuts(0, 1);
    }
    if (suffix == "11") {
      TrackCutsPhi->SetCutWindow(1.027, 1.1);
      evtCuts->SetSphericityCuts(0, 1);
    }

  //  if (suffix == "1") {
  //    TrackCutsPhi->SetCutWindow(0.987, 1.011);
  //  }
  //  if (suffix == "2") {
  //    TrackCutsPhi->SetCutWindow(1.027, 1.1);
  //  }
  //  if (suffix == "3") {
  //    TrackCutsPhi->SetCutWindow(1.1, 1.2);
  //  }
  //  if (suffix == "4") {
  //    TrackCutsPhi->SetCutWindow(1.2, 1.3);
  //  }
  //  if (suffix == "5") {
  //    TrackCutsPhi->SetCutWindow(1.3, 1.4);
  //  }
  //  if (suffix == "6") {
  //    TrackCutsPhi->SetCutWindow(1.4, 1.5);
  //  }
  //  if (suffix == "7") {
  //    TrackCutsPhi->SetCutWindow(1.5, 1.6);
  //  }
  //  if (suffix == "8") {
  //    TrackCutsPhi->SetCutWindow(1.6, 1.7);
  //  }
  //  if (suffix == "9") {
  //    TrackCutsPhi->SetCutWindow(1.7, 1.8);
  //  }
  //  if (suffix == "10") {
  //    TrackCutsPhi->SetCutWindow(1.8, 1.9);
  //  }
  //  if (suffix == "11") {
  //    TrackCutsPhi->SetCutWindow(1.9, 2);
  //  }
  //  if (suffix == "12") {
  //    evtCuts->SetSphericityCuts(0, 1);
  //  }
  //  if (suffix == "13") {
  //    TrackCutsPhi->SetCutWindow(0.987, 1.011);
  //    evtCuts->SetSphericityCuts(0, 1);
  //  }
  //  if (suffix == "14") {
  //    TrackCutsPhi->SetCutWindow(1.027, 1.1);
  //    evtCuts->SetSphericityCuts(0, 1);
  //  }
  //  if (suffix == "15") {
  //    TrackCutsPhi->SetCutWindow(1.1, 1.2);
  //    evtCuts->SetSphericityCuts(0, 1);
  //  }
  //  if (suffix == "16") {
  //    TrackCutsPhi->SetCutWindow(1.2, 1.3);
  //    evtCuts->SetSphericityCuts(0, 1);
  //  }
  //  if (suffix == "17") {
  //    TrackCutsPhi->SetCutWindow(1.3, 1.4);
  //    evtCuts->SetSphericityCuts(0, 1);
  //  }
  //  if (suffix == "18") {
  //    TrackCutsPhi->SetCutWindow(1.4, 1.5);
  //    evtCuts->SetSphericityCuts(0, 1);
  //  }
  //  if (suffix == "19") {
  //    TrackCutsPhi->SetCutWindow(1.5, 1.6);
  //    evtCuts->SetSphericityCuts(0, 1);
  //  }
  //  if (suffix == "20") {
  //    TrackCutsPhi->SetCutWindow(1.6, 1.7);
  //    evtCuts->SetSphericityCuts(0, 1);
  //  }
  //  if (suffix == "21") {
  //    TrackCutsPhi->SetCutWindow(1.7, 1.8);
  //    evtCuts->SetSphericityCuts(0, 1);
  //  }
  //  if (suffix == "22") {
  //    TrackCutsPhi->SetCutWindow(1.8, 1.9);
  //    evtCuts->SetSphericityCuts(0, 1);
  //  }
  //  if (suffix == "23") {
  //    TrackCutsPhi->SetCutWindow(1.9, 2);
  //    evtCuts->SetSphericityCuts(0, 1);
  //  }

  // Now we define stuff we want for our Particle collection
  // Thanks, CINT - will not compile due to an illegal constructor
  // std::vector<int> PDGParticles ={2212,2212,3122,3122,3312,3312};
  // First we need to tell him about the particles we mix, from the
  // PDG code the mass is obtained.
  std::vector<int> PDGParticles;
  PDGParticles.push_back(2212);  // 0 protons
  PDGParticles.push_back(2212);  // 1 antiprot
  PDGParticles.push_back(333);   // 2 v0 particle
  PDGParticles.push_back(2212);  // 3 proton truth
  PDGParticles.push_back(2212);  // 4 antiprot truth
  PDGParticles.push_back(333);   // 5 phi truth
  PDGParticles.push_back(333);   // 6 phiall
//  PDGParticles.push_back(2212);  // 7 proton common
//  PDGParticles.push_back(2212);  // 8 aproton common
//  PDGParticles.push_back(333);   // 9 phi common
//  PDGParticles.push_back(2212);  // 10 proton no prim
//  PDGParticles.push_back(2212);  // 11 antiprot no prim

  // We need to set the ZVtx bins
  std::vector<float> ZVtxBins;
  ZVtxBins.push_back(-10);
  ZVtxBins.push_back(-8);
  ZVtxBins.push_back(-6);
  ZVtxBins.push_back(-4);
  ZVtxBins.push_back(-2);
  ZVtxBins.push_back(0);
  ZVtxBins.push_back(2);
  ZVtxBins.push_back(4);
  ZVtxBins.push_back(6);
  ZVtxBins.push_back(8);
  ZVtxBins.push_back(10);
  // The Multiplicity bins are set here
  std::vector<int> MultBins;
  MultBins.push_back(0);
  MultBins.push_back(4);
  MultBins.push_back(8);
  MultBins.push_back(12);
  MultBins.push_back(16);
  MultBins.push_back(20);
  MultBins.push_back(24);
  MultBins.push_back(28);
  MultBins.push_back(32);
  MultBins.push_back(36);
  MultBins.push_back(40);
  MultBins.push_back(44);
  MultBins.push_back(48);
  MultBins.push_back(52);
  MultBins.push_back(56);
  MultBins.push_back(60);
  MultBins.push_back(64);
  MultBins.push_back(68);
  MultBins.push_back(72);
  MultBins.push_back(76);
  MultBins.push_back(80);

  // Number of bins
  std::vector<int> NBins;
  //  NBins.push_back(750);
  //  NBins.push_back(750);
  //  NBins.push_back(750);
  //  NBins.push_back(750);
  //  NBins.push_back(750);
  //  NBins.push_back(750);
  std::vector<float> kMin;
  // minimum k* value
  //  kMin.push_back(0.);
  //  kMin.push_back(0.);
  //  kMin.push_back(0.);
  //  kMin.push_back(0.);
  //  kMin.push_back(0.);
  //  kMin.push_back(0.);
  // maximum k* value
  std::vector<float> kMax;
  //  kMax.push_back(3.);
  //  kMax.push_back(3.);
  //  kMax.push_back(3.);
  //  kMax.push_back(3.);
  //  kMax.push_back(3.);
  //  kMax.push_back(3.);
  // pair QA extended
  std::vector<int> pairQA;
  //  pairQA.push_back(11); //pp
  //  pairQA.push_back(0); //pap
  //  pairQA.push_back(12); //pphi
  //  pairQA.push_back(11); //apap
  //  pairQA.push_back(12); //apphi
  //  pairQA.push_back(22); //phiphi

 // for (int i = 0; i < (1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 + 10 + 11 + 12); i++) {


  for (int i = 0; i < (1 + 2 + 3 + 4 + 5 + 6 + 7); i++) {
    NBins.push_back(750);
    kMin.push_back(0.);
    kMax.push_back(3.);
    pairQA.push_back(0);
  }

  pairQA[0] = 10;   // pp
  pairQA[1] = 10;    // pap
  pairQA[2] = 10;   // pphi
  pairQA[7] = 10;  // apap
  pairQA[8] = 10;  // apphi
  pairQA[13] = 10;  // phiphi

  if (isMC) {
    pairQA[18] = 10;  // TRUE
    pairQA[19] = 10;
    pairQA[20] = 10;
    pairQA[22] = 10;
    pairQA[23] = 10;
    pairQA[25] = 10;

//    pairQA[63] = 10;  // COMMON
//    pairQA[64] = 10;
//    pairQA[65] = 10;
//    pairQA[68] = 10;
//    pairQA[69] = 10;
//    pairQA[72] = 10;

    pairQA[21] = 10;  // phiALL pTRUE
    pairQA[24] = 10;

//    pairQA[55] = 10;  // phiTURE pNOprim
//    pairQA[56] = 10;
  }

  AliFemtoDreamCollConfig *config =
      new AliFemtoDreamCollConfig("Femto", "Femto");
  config->SetPtQA(true);
  config->SetMassQA(true);
  config->SetExtendedQAPairs(pairQA);
  config->SetZBins(ZVtxBins);
  config->SetMultBins(MultBins);
  config->SetMultBinning(true);
  config->SetPDGCodes(PDGParticles);
  config->SetNBinsHist(NBins);
  config->SetMinKRel(kMin);
  config->SetMaxKRel(kMax);

  if (suffix == "0") {
    config->SetUseEventMixing(true);
    config->SetMixingDepth(10);
  }
  if (suffix == "1") {
    config->SetUseEventMixing(true);
    config->SetMixingDepth(10);
  }
  if (suffix == "2") {
    config->SetUseEventMixing(true);
    config->SetMixingDepth(10);
  }
  if (suffix == "3") {
    config->SetUseEventMixing(true);
    config->SetMixingDepth(10);
  }
  if (suffix == "4") {
    config->SetUseEventMixing(true);
    config->SetMixingDepth(10);
  }
  if (suffix == "5") {
    config->SetUseEventMixing(true);
    config->SetMixingDepth(10);
  }
  if (suffix == "6") {
    config->SetUseEventMixing(false);
    config->SetUsePhiSpinning(true);
    config->SetControlMethod(AliFemtoDreamCollConfig::kPhiSpin);
    config->SetSpinningDepth(1);
  }
  if (suffix == "7") {
    config->SetUseEventMixing(false);
    config->SetUsePhiSpinning(true);
    config->SetControlMethod(AliFemtoDreamCollConfig::kPhiSpin);
    config->SetSpinningDepth(1);
  }
  if (suffix == "8") {
    config->SetUseEventMixing(false);
    config->SetUsePhiSpinning(true);
    config->SetControlMethod(AliFemtoDreamCollConfig::kPhiSpin);
    config->SetSpinningDepth(1);
  }
  if (suffix == "9") {
    config->SetUseEventMixing(false);
    config->SetUsePhiSpinning(true);
    config->SetControlMethod(AliFemtoDreamCollConfig::kPhiSpin);
    config->SetSpinningDepth(1);
  }
  if (suffix == "10") {
    config->SetUseEventMixing(false);
    config->SetUsePhiSpinning(true);
    config->SetControlMethod(AliFemtoDreamCollConfig::kPhiSpin);
    config->SetSpinningDepth(1);
  }
  if (suffix == "11") {
    config->SetUseEventMixing(false);
    config->SetUsePhiSpinning(true);
    config->SetControlMethod(AliFemtoDreamCollConfig::kPhiSpin);
    config->SetSpinningDepth(1);
  }


  //-------MIXED EVENTS---------------------------
  // config->SetUseEventMixing(true);
  // config->SetMixingDepth(10);

  //-------ROTATED SAMPLE-------------------------
  //  config->SetUseEventMixing(false);
  //  config->SetUsePhiSpinning(true);
  //  config->SetControlMethod(AliFemtoDreamCollConfig::kPhiSpin);
  //  config->SetSpinningDepth(1);

  /*
  //This is just to show off what would be possible in case you are interested,
  don't be confused by this at the beginning
  //you can just ignore it!
  if (false) {
    config->SetkTBinning(false);
    config->SetmTBinning(false);
    config->SetkTCentralityBinning(false);
    std::vector<float> centBins;
    centBins.push_back(20);
    centBins.push_back(40);
    centBins.push_back(90);
    config->SetCentBins(centBins);
    config->SetZBins(ZVtxBins);

    if (isMC) {
      config->SetMomentumResolution(true);
    } else {
      std::cout << "You are trying to request the Momentum Resolution without MC
  Info; fix it wont work! \n";
    }
    if (isMC) {
      config->SetPhiEtaBinnign(true);
    } else {
      std::cout << "You are trying to request the Eta Phi Plots without MC Info;
  fix it wont work! \n";
    }
  }
   */
  // now we create the task
  AliAnalysisTaskNanoAODFemtoDreamPhi *task =
      new AliAnalysisTaskNanoAODFemtoDreamPhi(
          "AliAnalysisTaskNanoAODFemtoDreamPhi", isMC);
  // THIS IS VERY IMPORTANT ELSE YOU DONT PROCESS ANY EVENTS
  // kINT7 == Minimum bias
  // kHighMultV0 high multiplicity triggered by the V0 detector
  if (CentEst == "kInt7") {
    task->SetTrigger(AliVEvent::kINT7);
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    std::cout << "Added kINT7 Trigger \n";
  } else if (CentEst == "kHM") {
    task->SetTrigger(AliVEvent::kHighMultV0);
    task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    std::cout << "Added kHighMultV0 Trigger \n";
  } else {
    std::cout << "============================================================="
                 "========"
              << std::endl;
    std::cout << "============================================================="
                 "========"
              << std::endl;
    std::cout << "Centrality Estimator not set, fix it else your Results will "
                 "be empty!"
              << std::endl;
    std::cout << "============================================================="
                 "========"
              << std::endl;
    std::cout << "============================================================="
                 "========"
              << std::endl;
  }

  // Throw all our settings to the task
  task->SetEventCuts(evtCuts);
  task->SetProtonCuts(TrackCuts);
  task->SetAntiProtonCuts(AntiTrackCuts);
  task->SetPosKaonCuts(TrackPosKaonCuts);
  task->SetNegKaonCuts(TrackNegKaonCuts);
  task->SetCollectionConfig(config);
  task->SetPhiCuts(TrackCutsPhi);
  task->SetUseDumpster(true);
  mgr->AddTask(task);

  TString file = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  mgr->ConnectInput(task, 0, cinput);

  AliAnalysisDataContainer *coutputQA;
  TString addon = "";
  if (CentEst == "kInt7") {
    addon += "MB";
  } else if (CentEst == "kHM") {
    addon += "HM";
  }
  TString QAName = Form("%sPhiResults%s", addon.Data(), suffix.Data());
  coutputQA = mgr->CreateContainer(QAName.Data(), TList::Class(),
                                   AliAnalysisManager::kOutputContainer,
                                   Form("%s:%s", file.Data(), QAName.Data()));
  mgr->ConnectOutput(task, 1, coutputQA);

  AliAnalysisDataContainer *coutputDumpsterQA;
  TString DumpsterName = Form("%sDumpster%s", addon.Data(), suffix.Data());
  coutputDumpsterQA = mgr->CreateContainer(
      DumpsterName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s:%s", file.Data(), DumpsterName.Data()));
  mgr->ConnectOutput(task, 2, coutputDumpsterQA);

  return task;
}
