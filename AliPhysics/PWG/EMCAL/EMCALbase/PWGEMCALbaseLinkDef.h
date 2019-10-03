#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AliEMCALRecoUtils+;

#pragma link C++ class AliAnalysisTaskEmcalLight+;
#pragma link C++ class AliAnalysisTaskEmcal+;
#pragma link C++ class AliAnalysisTaskEmcalEmbeddingHelper+;
#pragma link C++ class AliAnalysisTaskEmcalEmbeddingHelperData+;
#pragma link C++ class AliEmcalEmbeddingQA+;
#pragma link C++ class AliClusterContainer+;
#pragma link C++ class AliEmcalContainer+;
#pragma link C++ class AliEmcalContainerUtils+;
#pragma link C++ class AliEmcalParticle+;
#pragma link C++ class AliEmcalPhysicsSelection+;
#pragma link C++ class AliEmcalPythiaInfo+;
#pragma link C++ class AliEmcalManagedObject+;
#pragma link C++ class AliEmcalTrackSelection+;
#pragma link C++ class AliEmcalTrackSelectionESD+;
#pragma link C++ class AliEmcalTrackSelectionAOD+;
#pragma link C++ class AliParticleContainer+;
#pragma link C++ class AliPicoTrack+;
#pragma link C++ class AliMCParticleContainer+;
#pragma link C++ class AliTrackContainer+;
#pragma link C++ class AliTrackContainer::TrackOwnerHandler+;
#pragma link C++ class AliEmcalList+;
#pragma link C++ class std::map<std::string, AliParticleContainer*>+;
#pragma link C++ class std::pair<std::string, AliParticleContainer*>+;
#pragma link C++ class std::map<std::string, AliClusterContainer*>+;
#pragma link C++ class std::pair<std::string, AliClusterContainer*>+;

#pragma link C++ namespace PWG;
#pragma link C++ namespace PWG::EMCAL;
#pragma link C++ class PWG::EMCAL::AliEmcalDownscaleFactorsOCDB+;
#pragma link C++ class PWG::EMCAL::AliEmcalTrackSelResultPtr+;
#pragma link C++ class PWG::EMCAL::AliEmcalTrackSelResultUserPtr+;
#pragma link C++ class PWG::EMCAL::AliEmcalTrackSelResultUserStorage+;
#pragma link C++ class PWG::EMCAL::AliEmcalTrackSelResultCombined+;
#pragma link C++ class PWG::EMCAL::AliEmcalTrackSelResultHybrid+;
#pragma link C++ class PWG::EMCAL::AliEmcalAODFilterBitCuts+;
#pragma link C++ class PWG::EMCAL::AliEmcalCutBase+;
#pragma link C++ class PWG::EMCAL::AliEmcalVCutsWrapper+;
#pragma link C++ class PWG::EMCAL::AliEmcalAODHybridTrackCuts+;
#pragma link C++ class PWG::EMCAL::AliEmcalAODTPCOnlyTrackCuts+;
#pragma link C++ class PWG::EMCAL::AliEmcalESDHybridTrackCuts+;
#pragma link C++ class PWG::EMCAL::AliEmcalESDTrackCutsGenerator+;
#pragma link C++ class PWG::EMCAL::AliEmcalESDtrackCutsWrapper+;
#pragma link C++ class PWG::EMCAL::TestAliEmcalTrackSelResultPtr+;
#pragma link C++ class PWG::EMCAL::TestAliEmcalAODHybridTrackCuts+;
#pragma link C++ class PWG::EMCAL::TestAliEmcalTrackSelectionAOD+;
#pragma link C++ class std::vector<PWG::EMCAL::AliEmcalTrackSelResultPtr>+;
#endif
