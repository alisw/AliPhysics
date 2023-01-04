#ifndef AliAnalysisTaskPHOSPbPbQARun2_cxx
#define AliAnalysisTaskPHOSPbPbQARun2_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

// QA of PbPb data.
// Author: Boris Polishchuk

#include "AliAnalysisTaskSE.h"

class AliPHOSGeometry;

class AliAnalysisTaskPHOSPbPbQARun2 : public AliAnalysisTaskSE {

public:

  AliAnalysisTaskPHOSPbPbQARun2();
  AliAnalysisTaskPHOSPbPbQARun2(const char *name);
  virtual ~AliAnalysisTaskPHOSPbPbQARun2() {}

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);

private:

  AliAnalysisTaskPHOSPbPbQARun2(const AliAnalysisTaskPHOSPbPbQARun2&); // not implemented
  AliAnalysisTaskPHOSPbPbQARun2& operator=(const AliAnalysisTaskPHOSPbPbQARun2&); // not implemented

  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key

private:

  TList * fOutputContainer;   //final histogram container
  TList * fPHOSEvents[1][5];  //Containers for events with PHOS photons
  TClonesArray * fPHOSEvent ; //PHOS photons in current event
  Float_t fCentrality ;       //!Centrality of the currecnt event
  Int_t fCenBin ;             //! Current centrality bin
  AliPHOSGeometry  *fPHOSGeo; //! PHOS geometry
  Int_t fEventCounter;        // number of analyzed events
  TClonesArray *       fMCArray;  // MC array

  ClassDef(AliAnalysisTaskPHOSPbPbQARun2, 1); // PHOS analysis task
};

#endif
