//
// Add task to add the instances of the single track efficiency calculation
//  required for the automatic Montecarlo QA trending. This concerns:
//   charged particles with filter bit 0 or 4,
//   and pion/proton/kaon/electron with filter bit 0
//
// Authors: Jitendra Kumar, Zaida Conesa del Valle
//

#if (!defined(__CLING__) && !defined(__CINT__)) || defined(__ROOTCLING__) || defined(__ROOTCINT__)
#include "AliCFSingleTrackEfficiencyTask.h"
#endif

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
#include "PWGPP/EvTrkSelection/macros/AddSingleTrackEfficiencyTask.C"
#endif

AliAnalysisTask *AddSingleTrackEfficiencyTaskForAutomaticQA(const Bool_t readAOD = 0, // Flag to read AOD:1 or ESD:0
                                           ULong64_t triggerMask=AliVEvent::kAnyINT,
                                           Bool_t useCentrality = kFALSE)
{

    Info("AliCFSingleTrackEfficiencyTaskForAutomaticQA","Setting up instances");
    
    // The AddTask
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
#else
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EvTrkSelection/macros/AddSingleTrackEfficiencyTask.C");
#endif
    
    // Charged particles, filter bit 0
    AliCFSingleTrackEfficiencyTask * nchFB0 = AddSingleTrackEfficiencyTask(readAOD,"NchFbit0",AliPID::kPion,0,triggerMask,useCentrality);
    nchFB0->SetFilterType(0);
    
    // Charged particles, filter bit 4
    AliCFSingleTrackEfficiencyTask * nchFB4 = AddSingleTrackEfficiencyTask(readAOD,"NchFbit4",AliPID::kPion,0,triggerMask,useCentrality);
    nchFB4->SetFilterType(4);
    
    // Pions, filter bit 0
    AliCFSingleTrackEfficiencyTask * pions = AddSingleTrackEfficiencyTask(readAOD,"PionFbit0",AliPID::kPion,211,triggerMask,useCentrality);
    
    // Kaons, filter bit 0
    AliCFSingleTrackEfficiencyTask * kaons = AddSingleTrackEfficiencyTask(readAOD,"KaonFbit0",AliPID::kKaon,321,triggerMask,useCentrality);
    
    // Protons, filter bit 0
    AliCFSingleTrackEfficiencyTask * protons = AddSingleTrackEfficiencyTask(readAOD,"ProtonFbit0",AliPID::kProton,2212,triggerMask,useCentrality);
    
    // Electrons, filter bit 0
    AliCFSingleTrackEfficiencyTask * electrons = AddSingleTrackEfficiencyTask(readAOD,"ElectronFbit0",AliPID::kElectron,11,triggerMask,useCentrality);

    
    return NULL;
}
