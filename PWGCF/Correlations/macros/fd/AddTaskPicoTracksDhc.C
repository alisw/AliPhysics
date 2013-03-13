// $Id$

void AddTaskPicoTracksDhc(
  TString chNOutTracks = "PicoTracks"
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
    Info("AddTaskPicoTracksDhc","adding ESD track selection task ...");
    // ESD Track Cuts
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalEsdTpcTrack.C");
    AliEmcalEsdTpcTrackTask *hybTask = AddTaskEmcalEsdTpcTrack(chNIntermTracks.Data(),"Hybrid_LHC11h",kFALSE);
    hybTask->SelectCollisionCandidates(AliVEvent::kAny);
    // Pico Tracks
    pTrackTask = AddTaskEmcalPicoTrackMaker(chNOutTracks.Data(), "Hybrid_LHC11", "LHC11h");
    pTrackTask->SelectCollisionCandidates(AliVEvent::kAny);
  } else {
    Info("AddTaskPicoTracksDhc","AOD analysis, adding PicoTrack maker ...");
    pTrackTask = AddTaskEmcalPicoTrackMaker(chNOutTracks.Data(),"tracks","LHC11h",kFALSE);
    pTrackTask->SelectCollisionCandidates(AliVEvent::kAny);
  }
}
