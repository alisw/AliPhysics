#ifndef AliAnalysisTaskPHOSTriggerQA_cxx
#define AliAnalysisTaskPHOSTriggerQA_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

// QA of PHOS Trigger data.
// Author: Boris Polishchuk

#include "AliAnalysisTaskSE.h"

class AliPHOSGeometry;
class AliESDCaloCells;
class AliESDCaloCluster;

class AliAnalysisTaskPHOSTriggerQA : public AliAnalysisTaskSE {

public:

  AliAnalysisTaskPHOSTriggerQA();
  AliAnalysisTaskPHOSTriggerQA(const char *name);
  virtual ~AliAnalysisTaskPHOSTriggerQA() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  
private:

  AliAnalysisTaskPHOSTriggerQA(const AliAnalysisTaskPHOSTriggerQA&); // not implemented
  AliAnalysisTaskPHOSTriggerQA& operator=(const AliAnalysisTaskPHOSTriggerQA&); // not implemented

  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key

  void   MaxEnergyCellPos(AliESDCaloCells *cells, AliESDCaloCluster* clu, Int_t& maxId);
  Bool_t Matched(Int_t *trig_relid,Int_t *cluster_relid); //is cluster position coincides with 4x4 position?

private:

  TList * fOutputContainer;   //final histogram container
  AliPHOSGeometry  *fPHOSGeo; //! PHOS geometry
  Int_t fEventCounter;        // number of analyzed events

  ClassDef(AliAnalysisTaskPHOSTriggerQA, 1); // PHOS analysis task
};

#endif
