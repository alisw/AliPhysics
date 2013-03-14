// $Id$

void AddTaskPicoTracksDhc(
  TString chNOutTracks   = "PicoTracks",
  TString period         = "LHC11h",
  Bool_t includeNoITS    = kFALSE
) 
{
  // Get the analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskPicoTracksDhc", "No analysis manager found.");
    return;
  }
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
  
  // ESD or AOD? Create track cuts with pico track maker
  AliEmcalPicoTrackMaker *pTrackTask = 0x0;
  TString chIsESD("ESD");
  
  if (chIsESD.EqualTo(mgr->GetInputEventHandler()->GetDataType())) {
    TString cuts("Hybrid_");
    cuts += period;
    Info("AddTaskPicoTracksDhc","adding ESD track selection task ...");
    // ESD Track Cuts
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalEsdTpcTrack.C");
    AliEmcalEsdTpcTrackTask *hybTask = AddTaskEmcalEsdTpcTrack("HybridTracks", cuts.Data(), includeNoITS);
    hybTask->SelectCollisionCandidates(AliVEvent::kAny);
    // Pico Tracks
    pTrackTask = AddTaskEmcalPicoTrackMaker(chNOutTracks.Data(), "HybridTracks", period);
    pTrackTask->SelectCollisionCandidates(AliVEvent::kAny);
  } else {
    Info("AddTaskPicoTracksDhc","AOD analysis, adding PicoTrack maker ...");
    pTrackTask = AddTaskEmcalPicoTrackMaker(chNOutTracks.Data(),"tracks", period, includeNoITS);
    pTrackTask->SelectCollisionCandidates(AliVEvent::kAny);
  }
}
