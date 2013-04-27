#ifdef __CINT__
 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// PHOS_pp_pi0
#pragma link C++ class AliCaloPhoton+;
#pragma link C++ class AliAnalysisTaskPi0+;

//PHOS_PbPb
#pragma link C++ class AliAnalysisTaskPi0Flow+;
#pragma link C++ class AliAnalysisTaskPi0FlowMC+;
#pragma link C++ class AliPHOSTenderTask+;

//PHOS_PbPb_MC
#pragma link C++ class AliPHOSHijingEfficiency+;

// PHOS_embedding
#pragma link C++ class AliPHOSEmbedding+;
#pragma link C++ class AliAnalysisTaskPi0Efficiency+;
#pragma link C++ class AliAnalysisTaskPi0DiffEfficiency+;

// PHOS_PbPbQA
#pragma link C++ class AliAnalysisTaskPHOSPbPbQA+;

// PHOS_TriggerQA
#pragma link C++ class AliAnalysisTaskPHOSTriggerQA+;

// CaloCellQA
#pragma link C++ class AliCaloCellsQA+;
#pragma link C++ class AliAnalysisTaskCaloCellsQA+;

// Omega3pi
#pragma link C++ class AliAnalysisTaskOmegaPi0PiPi+;

// UserTasks
#pragma link C++ class AliCaloClusterInfo+;
#pragma link C++ class AliPHOSpPbPi0Header+;
#pragma link C++ class AliAnalysisTaskSEPHOSpPbPi0+;

#endif
