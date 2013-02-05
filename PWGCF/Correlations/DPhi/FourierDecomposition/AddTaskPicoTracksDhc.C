void AddTaskPicoTracksDhc(TString chNOutTracks = "PicoTracks") {
  TString chNIntermTracks = "HybridTracks";
  
  // Get the analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskPicoTracksDhc", "No analysis manager found.");
    return;
  }
  
  // ESD or AOD? Create track cuts with pico track maker
  AliEmcalPicoTrackMaker *pTrackTask = 0x0;
  TString chIsESD("ESD");
  
  if (chIsESD.EqualTo(mgr->GetInputEventHandler()->GetDataType())) {
    Info("AddTaskPicoTracksDhc","adding ESD track selection tasks ...");
    // ESD Track Cuts
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalEsdTpcTrack.C");
    AliEmcalEsdTpcTrackTask *hybTask = AddTaskEmcalEsdTpcTrack(chNIntermTracks.Data(),"Hybrid_LHC11h",kFALSE);
    // Pico Tracks
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
    pTrackTask = AddTaskEmcalPicoTrackMaker(chNOutTracks.Data(), chNIntermTracks.Data(), "LHC11h");
  }
  else {
    Info("AddTaskPicoTracksDhc","AOD analysis, adding picotrackmaker ...");
    pTrackTask = AddTaskEmcalPicoTrackMaker(chNOutTracks.Data(),"tracks","lhc11h",kFALSE);
  }

}