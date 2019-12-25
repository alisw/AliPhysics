#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// ClusterSelection
#pragma link C++ class AliPHOSTriggerUtils+;

// CorrectionFW
#pragma link C++ class AliPHOSCorrectionFW+;

// PHOS_Run2
#pragma link C++ class AliAnalysisTaskPHOSObjectCreator+;
#pragma link C++ class AliPHOSEventCuts+;
#pragma link C++ class AliPHOSClusterCuts+;
#pragma link C++ class AliPHOSTriggerHelper+;
#pragma link C++ class AliPHOSJetJetMC+;
#pragma link C++ class AliAnalysisTaskPHOSPi0EtaToGammaGamma+;
#pragma link C++ class AliAnalysisTaskPHOSEmbeddedDiffObjectCreator+;
#pragma link C++ class AliAnalysisTaskPHOSEmbedding+;
#pragma link C++ class AliAnalysisTaskPHOSEmbeddingEfficiency+;
#pragma link C++ class AliAnalysisTaskPHOSSingleSim+;


// PHOS_pp_8TeV_2012 
#pragma link C++ class AliCaloPhoton+;
#pragma link C++ class AliAnalysisTaskPHOSTrigPi0+;
#pragma link C++ class AliCaloTriggerSimulator+;

// PHOS_pp_pi0
#pragma link C++ class AliCaloPhoton+;
#pragma link C++ class AliAnalysisTaskPi0+;
#pragma link C++ class AliAnalysisTaskPi0Conversion+;
#pragma link C++ class AliAnalysisTaskPi0PP+;

//PHOS_PbPb
#pragma link C++ class AliAnalysisTaskPi0Flow+;
#pragma link C++ class AliAnalysisTaskGammaFlow+;
#pragma link C++ class AliAnalysisTaskgg+;
#pragma link C++ class AliAnalysisTaskggMC+;
#pragma link C++ class AliAnalysisTaskPHOSPCMgg+;
#pragma link C++ class AliAnalysisTaskPi0FlowMC+;
#pragma link C++ class AliAnalysisTaskEtaPhigg+;
#pragma link C++ class AliAnalysisTaskEtaPhiMultgg+;
#pragma link C++ class AliAnalysisTaskPi0FlowMCAOD+;
#pragma link C++ class AliAnalysisTaskPi0FlowMCHijing+;
#pragma link C++ class AliAnalysisTaskPi0FlowMCParamWeights+;
#pragma link C++ class AliPHOSTenderTask+;

//PHOS_EpRatio
#pragma link C++ class AliAnalysisTaskEpRatio+;

//PHOS_PbPb_MC
#pragma link C++ class AliPHOSHijingEfficiency+;

// PHOS_embedding
#pragma link C++ class AliPHOSEmbedding+;
#pragma link C++ class AliAnalysisTaskPi0Efficiency+;
#pragma link C++ class AliAnalysisTaskPi0DiffEfficiency+;

// PHOS_Run2embedding
#pragma link C++ class AliPHOSEmbeddingRun2+;
#pragma link C++ class AliPHOSEventCounter+;

//PHOS_GAFlow
#pragma link C++ class AliAnalysisTaskThermalGAFlow+;
#pragma link C++ class AliAnalysisTaskThermalGAFlowMC+;

//PHOS_NeutralMeson
#pragma link C++ class AliAnalysisTaskPHOSNeutralMeson+;

// PHOS_PbPbQA
#pragma link C++ class AliAnalysisTaskPHOSPbPbQA+;

// PHOS_Tagging
#pragma link C++ class AliAnalysisTaskTaggedPhotons+;

// PHOS_TriggerQA
#pragma link C++ class AliAnalysisTaskPHOSTriggerQA+;

// CaloCellQA
#pragma link C++ class AliCaloCellsQA+;
#pragma link C++ class AliAnalysisTaskCaloCellsQA+;

// CaloCellPhysQA
#pragma link C++ class AliCaloCellsPhysQA+;
#pragma link C++ class AliAnalysisTaskCaloCellsPhysQA+;

// Omega3pi
#pragma link C++ class AliAnalysisTaskOmegaPi0PiPi+;

// UserTasks
#pragma link C++ class AliCaloClusterInfo+;
#pragma link C++ class AliPHOSpPbPi0Header+;
#pragma link C++ class AliAnalysisTaskSEPHOSpPbPi0+;
#pragma link C++ class AliAnalysisTaskPHOSCluster+;
#pragma link C++ class AliCaloClusterContent+;

//PHOS_Correlations
#pragma link C++ class AliPHOSCorrelations+;

//PHOSCalib
#pragma link C++ class AliAnalysisTaskEmeanCalib+;
#pragma link C++ class AliAnalysisTaskPHOSTimeCalib+;

//CPV performance
#pragma link C++ class AliAnalysisTaskCPV+;

// LHC16_pp
#pragma link C++ class AliPP13AnalysisCluster+;
#pragma link C++ class AliPP13ClusterCuts+;
#pragma link C++ class AliPP13SelectionWeights+;
#pragma link C++ class AliPP13SelectionWeightsTOF+;
#pragma link C++ class AliPP13SelectionWeightsMC+;
#pragma link C++ class AliPP13SelectionWeightsFeeddown+;
#pragma link C++ class AliPP13SelectionWeightsSPMC+;
#pragma link C++ class AliPP13SelectionWeightsScan+;
#pragma link C++ class AliPP13DetectorHistogram+;
#pragma link C++ class SelectionLimits+;
#pragma link C++ class AliPP13PhysicsSelection+;
#pragma link C++ class AliPP13PhotonSpectrumSelection+;
#pragma link C++ class AliPP13QualityPhotonSelection+;
#pragma link C++ class AliPP13ParticlesHistogram+;
#pragma link C++ class AliPP13PhotonTimecutStudySelection+;
#pragma link C++ class AliPP13SpectrumSelection+;
#pragma link C++ class AliPP13SpectrumSelectionSimple+;
#pragma link C++ class AliPP13TagAndProbeSelection+;
#pragma link C++ class AliPP13MesonSelectionMC+;
#pragma link C++ class AliPP13EfficiencySelectionMC+;
#pragma link C++ class AliPP13EfficiencySelectionSPMC+;
#pragma link C++ class AliPP13PythiaInfoSelection+;
#pragma link C++ class AliPP13SpectrumSelectionMC+;
#pragma link C++ class AliPP13NonlinearityScanSelection+;
#pragma link C++ class AliPP13NonlinearitySelection+;
#pragma link C++ class AliPP13KaonToPionRatioMC+;
#pragma link C++ class AliPP13EpRatioSelection+;
#pragma link C++ class AliPP13FeeddownSelection+;
#pragma link C++ class AliPP13TriggerProperties+;
#pragma link C++ class AliPP13TriggerEfficiency+;
#pragma link C++ class AliPP13MixingSample+;
#pragma link C++ class AliAnalysisTaskPP13+;

//Resonances
#pragma link C++ class AliAnalysisPHOSResonances+ ;
#endif
