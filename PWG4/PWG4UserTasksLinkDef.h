#ifdef __CINT__
 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// PHOS_pp_pi0
#pragma link C++ class AliCaloPhoton+;
#pragma link C++ class AliAnalysisTaskPi0+;

// CaloCellQA
#pragma link C++ class AliCaloCellsQA+;
#pragma link C++ class AliAnalysisTaskCaloCellsQA+;

// EmcalTasks
#pragma link C++ class AliAnalysisTaskEMCALClusterizeFast+;
#pragma link C++ class AliAnalysisTaskEMCALPi0PbPb+;
#pragma link C++ class AliStaHeader+;
#pragma link C++ class AliStaCluster+;
#pragma link C++ class AliStaVertex+;
#pragma link C++ class AliStaTrigger+;
#pragma link C++ class AliStaPart+;

#endif
