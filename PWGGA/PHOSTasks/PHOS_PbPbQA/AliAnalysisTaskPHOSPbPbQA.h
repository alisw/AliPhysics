#ifndef AliAnalysisTaskPHOSPbPbQA_cxx
#define AliAnalysisTaskPHOSPbPbQA_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

// QA of PbPb data.
// Author: Boris Polishchuk

#include "AliAnalysisTaskSE.h"

class AliPHOSGeometry;

class AliAnalysisTaskPHOSPbPbQA : public AliAnalysisTaskSE {

public:

  AliAnalysisTaskPHOSPbPbQA();
  AliAnalysisTaskPHOSPbPbQA(const char *name);
  virtual ~AliAnalysisTaskPHOSPbPbQA() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  
private:

  AliAnalysisTaskPHOSPbPbQA(const AliAnalysisTaskPHOSPbPbQA&); // not implemented
  AliAnalysisTaskPHOSPbPbQA& operator=(const AliAnalysisTaskPHOSPbPbQA&); // not implemented

  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key

private:

  TList * fOutputContainer;   //final histogram container
  TList * fPHOSEvents[1][2];  //Containers for events with PHOS photons
  TClonesArray * fPHOSEvent ; //PHOS photons in current event
  Float_t fCentrality ;       //!Centrality of the currecnt event
  Int_t fCenBin ;             //! Current centrality bin
  AliPHOSGeometry  *fPHOSGeo; //! PHOS geometry
  Int_t fEventCounter;        // number of analyzed events

  ClassDef(AliAnalysisTaskPHOSPbPbQA, 1); // PHOS analysis task
};

#endif
