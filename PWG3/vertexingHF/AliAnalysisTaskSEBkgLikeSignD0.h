#ifndef AliAnalysisTaskSEBkgLikeSignD0_H
#define AliAnalysisTaskSEBkgLikeSignD0_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
// Class AliAnalysisTaskSEBkgLikeSignD0
// AliAnalysisTaskSE for reading both reconstructed D0 -> Kpi candidates
// and like sign pairs and for drawing corresponding distributions
// Author: C.Di Giglio, carmelo.digiglio@ba.infn.it
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"

class AliAnalysisTaskSEBkgLikeSignD0 : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSEBkgLikeSignD0();
  AliAnalysisTaskSEBkgLikeSignD0(const char *name);
  virtual ~AliAnalysisTaskSEBkgLikeSignD0();


  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

 private:

  AliAnalysisTaskSEBkgLikeSignD0(const AliAnalysisTaskSEBkgLikeSignD0 &source);
  AliAnalysisTaskSEBkgLikeSignD0& operator=(const AliAnalysisTaskSEBkgLikeSignD0& source); 

  TList   *fOutput;                //! list send on output slot 0
  TH1F    *fHistMassD0;            //! output histograms
  TH1F    *fHistMassLS;            //!
  TH1F    *fHistCtsD0;             //! Cosine of decay angle
  TH1F    *fHistCtsLS;             //!
  TH1F    *fHistCtsLSpos;          //!
  TH1F    *fHistCtsLSneg;          //!
  TH1F    *fHistCPtaD0;            //! Cosine of pointing angle
  TH1F    *fHistCPtaLS;            //!
  TH1F    *fHistd0d0D0;            //! Product of impact parameters
  TH1F    *fHistd0d0LS;            //!
  TH1F    *fHistDCAD0;             //! Distance of closest approach
  TH1F    *fHistDCALS;             //! like-sign
  AliAnalysisVertexingHF *fVHF;    // Vertexer heavy flavour (used to pass the cuts)
  TH1F *fNentries;                 //! histogram with number of events
  Int_t fTotPosPairs;              //
  Int_t fTotNegPairs;              // normalization
  Double_t fLsNormalization;       //
 
  ClassDef(AliAnalysisTaskSEBkgLikeSignD0,1); // comparison of unlike-sign and like-sign background for D0->Kpi
};

#endif
