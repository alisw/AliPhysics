#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// helper classes
#pragma link C++ class AliCutHandlerPCM+;
#pragma link C++ class AliIsoInfoHelper+;
#pragma link C++ class AliExtraClusterInfoHelper+;
#pragma link C++ class AliIdentifiedPrimarySelector+;
#pragma link C++ class AliIdentifiedPrimaryCuts+;
#pragma link C++ class AliGammaConversionAODBGHandler+;
#pragma link C++ class AliAnalysisTaskJetOutlierRemoval+;

// User tasks
#pragma link C++ class AliAnalysisTaskPi0v2+;
#pragma link C++ class AliAnalysisTaskConvJet+;
#pragma link C++ class AliConversionAodSkimTask+;
#pragma link C++ class AliAnalysisTaskGammaIsoTree+;
#pragma link C++ class AliAnalysisTaskGammaCaloMix+;
#pragma link C++ class AliAnalysisTaskGammaConvFlow+;

// QA tasks
#pragma link C++ class AliAnalysisTaskConversionQA+;
#pragma link C++ class AliAnalysisTaskConversionTree+;
#pragma link C++ class AliAnalysisTaskClusterQA+;
#pragma link C++ class AliAnalysisTaskMaterial+;
#pragma link C++ class AliAnalysisTaskMaterialHistos+;
#pragma link C++ class AliAnalysisTaskResolution+;
#pragma link C++ class AliAnalysisTaskConvCaloTree+;
#pragma link C++ class AliAnalysisTaskQA+;
#pragma link C++ class AliAnalysisTaskGammaTriggerQA+;
#pragma link C++ class AliAnalysisTaskConvCaloCalibration+;
#pragma link C++ class AliAnalysisTaskTrackQA+;
#pragma link C++ class AliAnalysisTRDEfficiency+;
#pragma link C++ class AliAnalysisTaskElectronStudies+;

// MC and cocktail tasks
#pragma link C++ class AliAnalysisTaskGammaPureMC+;
#pragma link C++ class AliAnalysisTaskGammaCocktailMC+;
#pragma link C++ class AliAnalysisTaskHadronicCocktailMC+;
#pragma link C++ class AliAnalysisTaskOmegaMCStudies+;

// main pion, eta and photon tasks
#pragma link C++ class AliAnalysisTaskGammaConvV1+;
#pragma link C++ class AliAnalysisTaskGammaConvDalitzV1+;
#pragma link C++ class AliAnalysisTaskGammaCaloMerged+;
#pragma link C++ class AliAnalysisTaskGammaCaloMergedML+;
#pragma link C++ class AliAnalysisTaskGammaConvCalo+;
#pragma link C++ class AliAnalysisTaskGammaCalo+;
#pragma link C++ class AliAnalysisTaskGammaCaloDalitzV1+;
#pragma link C++ class AliAnalysisTaskGammaConvCaloIso+;
#pragma link C++ class AliAnalysisTaskGammaCaloIso+;

// heavier resonances
#pragma link C++ class AliAnalysisTaskOmegaToPiZeroGamma+;
#pragma link C++ class AliAnalysisTaskHeavyNeutralMesonToGG+;
#pragma link C++ class AliAnalysisTaskK0toPi0Pi0+;
#pragma link C++ class AliAnalysisTaskGammaConvDtrue+;
#pragma link C++ class AliAnalysisTaskSigmaPlToProtonPiZero+;
#pragma link C++ class AliAnalysisTaskSigmaPlToProtonPiZeroAOD+;
#pragma link C++ class AliAnalysisTaskEtaToPiPlPiMiGamma+;
#pragma link C++ class AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson+;

// Old tasks
// #pragma link C++ class AliAnaConvIsolation+;
// #pragma link C++ class AliAnaConvCorrBase+;
// #pragma link C++ class AliAnaConvCorrPion+;
// #pragma link C++ class AliAnaConvCorrPhoton+;
// #pragma link C++ class AliAnalysisTaskdPhi+;
// #pragma link C++ class AliAnalysisTaskGCPartToPWG4Part+;
#pragma link C++ class AliPrimaryPionSelector+;
#pragma link C++ class AliPrimaryPionCuts+;
// #pragma link C++ class AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero+;

#endif
