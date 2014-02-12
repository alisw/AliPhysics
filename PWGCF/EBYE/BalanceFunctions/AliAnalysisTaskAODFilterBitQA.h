#ifndef ALIANALYSISTASKAODFILTERBITQA_H
#define ALIANALYSISTASKAODFILTERBITQA_H

// Analysis task for the QA of AOD track filter bits
// Authors: m.weber@cern.ch

#define gBitMax 16// number of maximum filter bits

class TList;
class TH2D;
class TH3D;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskAODFilterBitQA : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskAODFilterBitQA(const char *name = "AliAnalysisTaskAODFilterBitQA");
  virtual ~AliAnalysisTaskAODFilterBitQA(); 
   
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   FinishTaskOutput();
  virtual void   Terminate(Option_t *);

 private:
  Double_t IsEventAccepted(AliVEvent *event);
  void GetAcceptedTracks(AliVEvent *event, Double_t gCentrality);
  
  AliAnalysisTaskAODFilterBitQA(const AliAnalysisTaskAODFilterBitQA&); // not implemented
  AliAnalysisTaskAODFilterBitQA& operator=(const AliAnalysisTaskAODFilterBitQA&); // not implemented


  TList *fListQA;//output list for QA histograms

  TH2D* fHistTrackStats;//QA histogram for track filter bit statistics vs. centrality
  TH3D* fHistKinematics[gBitMax];//QA histograms for kinematics (eta, phi, pT) for different filter bits
  TH2D* fHistDCA[gBitMax];//QA histograms for DCA (xy,z) for different filter bits
  TH2D* fHistDCAprop[gBitMax];//QA histograms for DCA (xy,z) for different filter bits after PropagateToDCA
  TH2D* fHistChiClus[gBitMax];//QA histograms for Chi2 and number of TPC clusters for different filter bits


  
  ClassDef(AliAnalysisTaskAODFilterBitQA, 0); //
};



#endif
