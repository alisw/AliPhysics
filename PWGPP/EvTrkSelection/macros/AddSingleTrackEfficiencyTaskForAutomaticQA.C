//
// Add task to add the instances of the single track efficiency calculation
//  required for the automatic Montecarlo QA trending. This concerns:
//   charged particles with filter bit 0 or 4,
//   and pion/proton/kaon/electron with filter bit 0
//
// Authors: Jitendra Kumar, Zaida Conesa del Valle
//

class AliAnalysisTask;

AliAnalysisTask *AddSingleTrackEfficiencyTaskForAutomaticQA(const Bool_t readAOD = 0, // Flag to read AOD:1 or ESD:0
                                           ULong64_t triggerMask=AliVEvent::kAnyINT,
                                           Bool_t useCentrality = kFALSE)
{
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EvTrkSelection/macros/AddSingleTrackEfficiencyTask.C");
    gROOT->Macro(TString::Format("$ALICE_PHYSICS/PWGPP/EvTrkSelection/main_AddSingleTrackEfficiencyTaskForAutomaticQA.C(%d, %llu, %d)", readAOD, triggerMask, useCentrality));
    return NULL;
}
