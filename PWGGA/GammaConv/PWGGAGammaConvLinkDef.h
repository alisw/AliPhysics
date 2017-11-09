#ifdef __CINT__
 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// Base classes
#pragma link C++ class AliConversionPhotonBase+;
#pragma link C++ class AliAODConversionParticle+;
#pragma link C++ class AliAODConversionMother+;
#pragma link C++ class AliAODConversionPhoton+;
#pragma link C++ class AliKFConversionPhoton+;
#pragma link C++ class AliKFConversionMother+;
#pragma link C++ class AliCaloPhotonCuts+;
#pragma link C++ class AliConvEventCuts+;
#pragma link C++ class AliConversionPhotonCuts+;
#pragma link C++ class AliConversionCuts+;
#pragma link C++ class AliConversionSelection+;
#pragma link C++ class AliV0ReaderV1+;
#pragma link C++ class AliConversionAODBGHandlerRP+;
#pragma link C++ class AliConversionTrackCuts+;
#pragma link C++ class AliConversionMesonCuts+;
#pragma link C++ class AliDalitzElectronCuts+;
#pragma link C++ class AliDalitzElectronSelector+;
#pragma link C++ class AliCaloTrackMatcher+;

// User tasks
#pragma link C++ class AliAnalysisTaskPi0v2+;
#pragma link C++ class AliGammaConversionAODBGHandler+;
#pragma link C++ class AliAnalysisTaskGammaConvV1+;
#pragma link C++ class AliAnalysisTaskGammaConvDalitzV1+;
#pragma link C++ class AliAnalysisTaskConversionQA+;
#pragma link C++ class AliAnalysisTaskMaterial+;
#pragma link C++ class AliAnalysisTaskMaterialHistos+;
#pragma link C++ class AliAnalysisTaskResolution+;
#pragma link C++ class AliAnalysisTaskGammaCaloDalitzV1+;
#pragma link C++ class AliAnalysisTaskGammaPureMC+;
#pragma link C++ class AliAnalysisTaskGammaCocktailMC+;
#pragma link C++ class AliAnalysisTaskHadronicCocktailMC+;

#pragma link C++ class AliAnaConvIsolation+;
#pragma link C++ class AliAnaConvCorrBase+;
#pragma link C++ class AliAnaConvCorrPion+;
#pragma link C++ class AliAnaConvCorrPhoton+;
#pragma link C++ class AliAnalysisTaskdPhi+;
#pragma link C++ class AliAnalysisTaskOmegaToPiZeroGamma+;
#pragma link C++ class AliAnalysisTaskK0toPi0Pi0+;

// Old tasks
#pragma link C++ class AliAnalysisTaskGCPartToPWG4Part+;

#pragma link C++ class AliPrimaryPionSelector+;
#pragma link C++ class AliPrimaryPionCuts+;
#pragma link C++ class AliAnalysisTaskEtaToPiPlPiMiGamma+;
#pragma link C++ class AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero+;
#pragma link C++ class AliAnalysisTaskGammaConvCalo+;
#pragma link C++ class AliAnalysisTaskGammaCalo+;
#pragma link C++ class AliAnalysisTaskGammaCaloMerged+;
#pragma link C++ class AliAnalysisTaskGammaConvFlow+;

#pragma link C++ class AliV0ReaderStrange+;
#pragma link C++ class AliV0CutsStrange+;
#pragma link C++ class AliV0ParticleStrange+;

#endif
