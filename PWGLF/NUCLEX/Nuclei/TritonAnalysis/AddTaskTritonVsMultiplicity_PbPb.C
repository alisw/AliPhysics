#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Rtypes.h>
#include <TString.h>
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliPID.h"
#endif

//____________________________________________________________________________________
AliAnalysisTask *AddTaskTritonVsMultiplicity_PbPb(unsigned long fTriggerMask, const char *filename, Double_t CentralityMin, Double_t CentralityMax, Double_t VertexZmin, Double_t VertexZmax, Int_t NumberOfVertexContributorsMin, const char *CentralityEstimator, Double_t PtMin, Double_t PtMax, Double_t EtaMax, Double_t YMax, Int_t NumberClustersITSMin, Int_t NumberClustersTPCMin, Int_t NumberCrossedRowsTPCMin, Double_t CrossedRowsFindableClsMin, Int_t NumberClustersTPCdEdxMin, Double_t ChiSquarePerNDFMax, const char *ITSrequirement, Double_t DCAzMax, Double_t DCAxyMax, Double_t nSigmaTOFmax, Double_t nSigmaTPCmax, Int_t TRDntracklets, TF1 *MeanTPC, TF1 *WidthTPC, TF1 *MeanTOF, TF1 *WidthTOF )  {

  //Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskTritonVsMultiplicity_PbPb","No analysis manager found.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskTritonVsMultiplicity_PbPb", "This task requires an input event handler");
    return NULL;
  }   
  
  //Input container
    AliAnalysisDataContainer *input = mgr -> GetCommonInputContainer();
  
  //Analysis Task
    AliAnalysisTaskTritonVsMultiplicity_PbPb *task = new AliAnalysisTaskTritonVsMultiplicity_PbPb ("task");
    task -> SelectCollisionCandidates (fTriggerMask);
    task -> SetCentralityRange (CentralityMin,CentralityMax);
    task -> SetEventCuts       (VertexZmin,VertexZmax,NumberOfVertexContributorsMin,CentralityEstimator);
    task -> SetTrackCuts (PtMin, PtMax, EtaMax, YMax, NumberClustersITSMin,NumberClustersTPCMin, NumberCrossedRowsTPCMin,CrossedRowsFindableClsMin, NumberClustersTPCdEdxMin, ChiSquarePerNDFMax, ITSrequirement,DCAzMax, DCAxyMax, nSigmaTOFmax, nSigmaTPCmax, TRDntracklets);
    task -> RecalibrationTPCandTOF (fMeanTPC_central, fWidthTPC_central, fMeanTOF_central, fWidthTOF_central );
    //task -> RecalibrationTPCandTOF (fMeanTPC_semicentral, fWidthTPC_semicentral, fMeanTOF_semicentral, fWidthTOF_semicentral );
    
       
    mgr -> AddTask(task);
    
    AliAnalysisDataContainer *outputData = mgr -> CreateContainer("Input",TList::Class(),AliAnalysisManager::kOutputContainer,filename);
    AliAnalysisDataContainer *outputQA   = mgr -> CreateContainer("InputQA",TList::Class(),AliAnalysisManager::kOutputContainer,"QA.root");
    mgr -> ConnectInput (task,0,input);
    mgr -> ConnectOutput(task,1,outputData);
    mgr -> ConnectOutput(task,2,outputQA);
    
    

    return task;
  
}
//____________________________________________________________________________________
