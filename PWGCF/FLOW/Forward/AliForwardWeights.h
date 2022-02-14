#ifndef AliForwardWeights_cxx
#define AliForwardWeights_cxx
/**
 * @file AliForwardWeights.h
 * @author Freja Thoresen <freja.thoresen@cern.ch>
 *
 * @brief
 *
 * @ingroup pwgcf_forward_flow
 */
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTaskSE.h"
#include "TFile.h"
#include "TList.h"
#include "AliAnalysisManager.h"
#include "AliForwardSettings.h"
class AliForwardWeights : public TObject {

 public:
  AliForwardWeights();
  AliForwardSettings fSettings;

  void connectNUA();
  void connectNUE();
  void connectSec();
  void connectSecCent();
  AliAnalysisDataContainer* makeWeightContainer(TString nua_file, TString containerName);
  AliAnalysisDataContainer* makeWeightContainerNUE(TString nue_file, TString containerName);
  AliAnalysisDataContainer* makeWeightContainerSec(TString sec_file, TString containerName);
  AliAnalysisDataContainer* makeWeightContainerSecCent(TString sec_file, TString containerName);
  void connectContainer(AliAnalysisDataContainer* container);
  void connectSecContainer(AliAnalysisDataContainer* container);
  void connectNUEContainer(AliAnalysisDataContainer* container);
  void connectSecCentContainer(AliAnalysisDataContainer* container);

  private:
    ClassDef(AliForwardWeights, 1);  
 };
#endif
