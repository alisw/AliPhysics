#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

/// Exotica
/// * Hdibaryon
#pragma link C++ class AliAnalysisTaskHdibaryonLPpi+;
/// * LambdaN
#pragma link C++ class AliAnalysisTaskLambdaNAOD+;
#pragma link C++ class AliAnalysisTaskLambdaNRun2+;
#pragma link C++ class AliAnalysisTaskLambdaNRun2::AnalysisV0+;
#pragma link C++ class AliAnalysisTaskLambdaNRun2::AnalysisEvent+;
/// * LambdaNN
#pragma link C++ class AliAnalysisTaskLNNntuple+;
#pragma link C++ class AliAnalysisTaskLNNv0Bkg+;
/// * Ps
#pragma link C++ class AliAnalysisTaskPsEfficiency+;
/// * dStar
#pragma link C++ class AliAnalysisTaskdStar+;
#if defined __CINT__ && !defined __CLING__
#pragma link C++ class ROOT::Math::PtEtaPhiM4D<float>+;
#pragma link C++ class ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>+;
#endif
#pragma link C++ class daughter_struct+;
#pragma link C++ class std::vector<daughter_struct>+;
#pragma link C++ typedef FourVector_t;

/// Nuclei
/// * Absorption
#pragma link C++ class AliAnalysisTaskDeuteronAbsorption+;
/// * AbsorptionRatios
#pragma link C++ class AliAnalysisTaskLightN+;
#pragma link C++ class AliLightNEventHist+;
#pragma link C++ class AliLightNAnalysis+;
#pragma link C++ class AliLightNTrack+;
#pragma link C++ class AliLightNBasePart+;
#pragma link C++ class AliLightNTrackCuts+;
#pragma link C++ class AliLightNEvent+;
#pragma link C++ class AliLightNTrackHist+;
#pragma link C++ class AliLightNEventCuts+;
#pragma link C++ class AliLightNTrackMCHist+;
#pragma link C++ class AliAnalysisTaskAntipd+;
#pragma link C++ class AliAnalysisTaskHe3+;
#pragma link C++ class AliAnalysisTaskHe3_ESD+;
/// * DeltaMasses
#pragma link C++ class AliAnalysisNucleiMass+;
#pragma link C++ class AliAnalysisNuclMult+;
/// * He4
#pragma link C++ class AliAnalysisTaskAntiHe4+;
/// * He4pp
#pragma link C++ class AliAnalysisHe4+;
/// * NucleiPbPb
#pragma link C++ class AliAnalysisTaskNucleiYield+;
#pragma link C++ class AliAnalysisTaskNucleiPIDqa+;
#pragma link C++ class AliAnalysisTaskSignalLoss+;
#pragma link C++ class RLightNucleus+;
#pragma link C++ class SLightNucleus+;
/// * Triton
#pragma link C++ class AliAnalysisTaskTritonVsMultiplicity_PbPb+;
/// * ReducedTreeNuclei
#pragma link C++ class AliAnalysisTaskReducedTreeNuclei+;
#pragma link C++ class AliAnalysisTaskReducedTreeHypertriton+;
/// * v2
#pragma link C++ class AliAnalysisTaskNucleiv2+;
#pragma link C++ class AliAnalysisTaskNucleiv2SP+;
#pragma link C++ class AliAnalysisTaskNucleiv2pPb+;
#pragma link C++ class AliAnalysisTaskAllPtcv2+;
#pragma link C++ class AliAnalysishDEventCollection+;  
#pragma link C++ class AliReconstructed2pcFirst+;
#pragma link C++ class AliReconstructed2pcSecond+;
#pragma link C++ class AliAnalysishDEvent+;
#pragma link C++ class AliAnalysishDEventCollection+;
#pragma link C++ class AliAnalysisTaskDeuFlow2PC+;
#pragma link C++ class AliAnalysisTaskHypv2PbPb18+;
#pragma link C++ class AliAnalysisTaskNucleiv2PbPb18+;

/// * NucleiKine
#pragma link C++ class AliAnalysisTaskNucleiKine+;
#pragma link C++ class AliAnalysisTaskNucleiKineCor+;
/// Hypernuclei
/// * Hyp2body
#pragma link C++ class AliAnalysisTaskHelium3Pi+;
#pragma link C++ class AliAnalysisTaskHelium3PiAOD+;
#pragma link C++ class AliAnalysisTaskHelium3PiMC+;
#pragma link C++ class AliAnalysisTaskHypTritEventTree+;
#pragma link C++ class AliReducedHypTritV0+;
#pragma link C++ class AliReducedHypTritTrack+;
#pragma link C++ class AliReducedHypTritEvent+;
#pragma link C++ class AliAnalysisTaskS3ParticleYields+;
#pragma link C++ class AliAnalysisTaskHe3EffTree+;
#pragma link C++ class AliAnalysisTaskHypCrossCheck+;
#pragma link C++ class AliAnalysisTaskHyperTriton2He3piML+;
#pragma link C++ class AliAnalysisTaskReducedTreeHypertritonBindingEnergy+;
#pragma link C++ class RHyperTritonHe3pi+;
#pragma link C++ class SHyperTritonHe3pi+;
#pragma link C++ class std::vector<RHyperTritonHe3pi>+;
#pragma link C++ class std::vector<SHyperTritonHe3pi>+;
#pragma link C++ class RCollision+;
#pragma link C++ class RTracklet+;
#pragma link C++ class std::vector<RTracklet>+;
#pragma link C++ class SGenericV0+;
#pragma link C++ class SGenericTracklet+;
#pragma link C++ class std::vector<SGenericV0>+;
#pragma link C++ class std::vector<SGenericTracklet>+;

/// * Hyp3body
#pragma link C++ class AliAnalysisTaskFindableHypertriton3+;
#pragma link C++ class AliAnalysisTaskHypertriton3+;

/// * KF2Body
#pragma link C++ class AliAnalysisTaskHypertritonKFTree+;

/// Hyp3body Selector
#ifdef __CLING__
#pragma link C++ class AliSelectorFindableHyperTriton3Body+;
#pragma link C++ class AliAnalysisTaskHypertriton3ML+;
#pragma link C++ class AliAnalysisTaskHyperTriton3KF+;
#pragma link C++ class REvent3KF+;
#pragma link C++ class RHyperTriton3KF+;
#pragma link C++ class SHyperTriton3KF+;
#pragma link C++ class std::vector<RHyperTriton3KF>+;
#pragma link C++ class std::vector<SHyperTriton3KF>+;
#pragma link C++ class AliAnalysisTaskHypertritonO2+;
#pragma link C++ class RHyperTritonO2+;
#pragma link C++ class SHyperTritonO2+;
#pragma link C++ class std::vector<RHyperTritonO2>+;
#pragma link C++ class std::vector<SHyperTritonO2>+;
#pragma link C++ class RHypertriton3+;
#pragma link C++ class SHypertriton3+;
#pragma link C++ class MLSelected+;
#pragma link C++ class std::vector<RHypertriton3>+;
#pragma link C++ class std::vector<SHypertriton3>+;
#pragma link C++ class std::vector<MLSelected>+;
#pragma link C++ class REvent+;
#pragma link C++ class o2::vertexing::DCAFitter2+;
#pragma link C++ class o2::vertexing::DCAFitter3+;
#endif

/// * VertexerHyp3Body
#pragma link C++ class AliVertexerHyperTriton2Body+;
#pragma link C++ class AliVertexerHyperTriton3Body+;

// * DoubleHypNuc
#pragma link C++ class AliAnalysisTaskDoubleHypNucTree+;

/// Utils
/// * CODEX
#pragma link C++ class AliAnalysisCODEX::Header+;
#pragma link C++ class AliAnalysisCODEX::Track+;
#pragma link C++ class std::vector<AliAnalysisCODEX::Track>+;
#pragma link C++ class AliAnalysisCODEXtask+;
/// * NanoAOD
#pragma link C++ class AliNanoFilterPID+;
#pragma link C++ class AliNanoSkimmingPID+;
/// * ChunkFilter
#pragma link C++ class AliAnalysisTaskFilterHe3+;
#endif
