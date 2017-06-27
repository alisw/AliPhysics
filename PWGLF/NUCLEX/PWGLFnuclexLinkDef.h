#ifdef __CINT__

#pragma link off all glols;
#pragma link off all classes;
#pragma link off all functions;

/// Exotica
/// * Hdibaryon
#pragma link C++ class AliAnalysisTaskHdibaryonLPpi+;
/// * LambdaN
#pragma link C++ class AliAnalysisTaskLambdaNAOD+;
/// * NOmega
#pragma link C++ class AliAnalysisTaskNOmegaLPK+;
#pragma link C++ class AliAnalysisTaskOmegaOmegaOX+;
#pragma link C++ class AliAnalysisTaskNOmegaLX+;
/// * Ps
#pragma link C++ class AliAnalysisTaskPsEfficiency+;

/// Nuclei
/// * DeltaMasses
#pragma link C++ class AliAnalysisNucleiMass+;
#pragma link C++ class AliAnalysisNuclMult+;
/// * DeuteronpA
#pragma link C++ class AliAnalysisDeuteronpA+;
#pragma link C++ class AliAnalysisDeuteronTree+;
/// * He4
#pragma link C++ class AliAnalysisTaskAntiHe4+;
/// * NucleiPbPb
#pragma link C++ class AliAnalysisTaskNucleiYield+;
#pragma link C++ class AliAnalysisTaskNucleiYieldESD+;
#pragma link C++ class AliAnalysisTaskNucleiPIDqa+;
/// * EventCuts
#pragma link C++ class AliNuclexEventCuts+;
#pragma link C++ class AliNuclexEventCutsContainer+;

/// * Nucleipp
#pragma link C++ class AliLnID+;
#pragma link C++ class AliLnHistoMap+;
#pragma link C++ class AliLnAODtrackCuts+;
#pragma link C++ class AliAnalysisTaskB2+;
#pragma link C++ class AliAnalysisTaskB2AOD+;
/// * v2
#pragma link C++ class AliAnalysisTaskNucleiv2+;
#pragma link C++ class AliAnalysisTaskNucleiv2SP+;
#pragma link C++ class AliAnalysisTaskNucleiv2pPb+;
#pragma link C++ class AliAnalysisTaskAllPtcv2+;
/// * NucleiKine
#pragma link C++ class AliAnalysisTaskNucleiKine+;

/// Hypernuclei
/// * Hyp2body
#pragma link C++ class AliAnalysisTaskHelium3Pi+;
#pragma link C++ class AliAnalysisTaskHelium3PiAOD+;
#pragma link C++ class AliAnalysisTaskHelium3PiMC+;
#pragma link C++ class AliAnalysisTaskHypTritEventTree+;
#pragma link C++ class AliReducedHypTritV0+;
#pragma link C++ class AliReducedHypTritTrack+;
#pragma link C++ class AliReducedHypTritEvent+;
#pragma link C++ class AliAnalysisTaskHypCrossCheck+;

/// * Hyp3body
#pragma link C++ class AliAnalysisTaskHypertriton3+;
#pragma link C++ class AliAnalysisTaskHypertriton3Dev+;
#pragma link C++ class AliAnalysisTaskHypertriton3AOD+;

/// Utils
/// * RecoDecay
#pragma link C++ class AliAODRecoDecayLF+;
#pragma link C++ class AliAODRecoDecayLF2Prong+;
/// * NuclexFilter
#pragma link C++ class AliAODNuclExReplicator+;
#pragma link C++ class AliAnalysisTaskESDNuclExFilter+;
#pragma link C++ class AliAODMCNuclExReplicator+;
#pragma link C++ class AliAnalysisTaskESDNuclExFilterMC+;
#pragma link C++ class AliAnalysisTaskReadNuclexAOD+;
/// * CODEX
#pragma link C++ class AliAnalysisCODEX::Header+;
#pragma link C++ class AliAnalysisCODEX::Track+;
#pragma link C++ class std::vector<AliAnalysisCODEX::Track>+;
#pragma link C++ class AliAnalysisCODEXtask+;
#endif

