#ifdef __CINT__
 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// PHOS_pp_pi0
#pragma link C++ class AliCaloPhoton+;
#pragma link C++ class AliAnalysisTaskPi0+;

// PHOS_embedding
#pragma link C++ class AliPHOSEmbedding+;
#pragma link C++ class AliAnalysisTaskPi0Efficiency+;
#pragma link C++ class AliAnalysisTaskPi0DiffEfficiency+;
#pragma link C++ class AliPHOSDigitDecalibrate+;

// CaloCellQA
#pragma link C++ class AliCaloCellsQA+;
#pragma link C++ class AliAnalysisTaskCaloCellsQA+;

// EmcalTasks
#pragma link C++ class AliEmcalPhysicsSelection+;
#pragma link C++ class AliEmcalPhysicsSelectionTask+;
#pragma link C++ class AliEmcalEsdTpcTrackTask+;
#pragma link C++ class AliAnalysisTaskEMCALPi0PbPb+;
#pragma link C++ class AliAnalysisTaskEMCALClusterizeFast+;
#pragma link C++ class AliAnalysisTaskEMCALPi0PbPb+;
#pragma link C++ class AliStaHeader+;
#pragma link C++ class AliStaCluster+;
#pragma link C++ class AliStaVertex+;
#pragma link C++ class AliStaTrigger+;
#pragma link C++ class AliStaPart+;
#pragma link C++ class AliAnalysisTaskEMCALTriggerQA+;

#endif
