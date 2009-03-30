#ifndef AliAnalysisTaskSEBkgLikeSignJPSI_H
#define AliAnalysisTaskSEBkgLikeSignJPSI_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskSEBkgLikeSignJPSI
// AliAnalysisTaskSE for reading both reconstructed JPSI -> ee candidates
// and like sign pairs and for drawing corresponding distributions
// Author: C.Di Giglio, carmelo.digiglio@ba.infn.it
//*************************************************************************

#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"

class AliAnalysisTaskSEBkgLikeSignJPSI : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSEBkgLikeSignJPSI();
  AliAnalysisTaskSEBkgLikeSignJPSI(const char *name);
  virtual ~AliAnalysisTaskSEBkgLikeSignJPSI();


  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

 private:

  AliAnalysisTaskSEBkgLikeSignJPSI(const AliAnalysisTaskSEBkgLikeSignJPSI &source);
  AliAnalysisTaskSEBkgLikeSignJPSI& operator=(const AliAnalysisTaskSEBkgLikeSignJPSI& source); 

  TList   *fOutput;                //! list send on output slot 0
  TH1F    *fHistMassJPSI;          // output histograms
  TH1F    *fHistMassLS;            //
  TH1F    *fHistCtsJPSI;           // Cosine of decay angle
  TH1F    *fHistCtsLS;             //
  TH1F    *fHistCtsLSpos;          //
  TH1F    *fHistCtsLSneg;          //
  TH1F    *fHistCPtaJPSI;          // Cosine of pointing angle
  TH1F    *fHistCPtaLS;            //
  TH1F    *fHistd0d0JPSI;          // Product of impact parameters
  TH1F    *fHistd0d0LS;            //
  TH1F    *fHistDCAJPSI;           // Distance of closest approach
  TH1F    *fHistDCALS;             //
  AliAnalysisVertexingHF *fVHF;    // Vertexer heavy flavour (used to pass the cuts)

  Int_t fTotPosPairs;              //
  Int_t fTotNegPairs;              // normalization
  Double_t fLsNormalization;       //
 
  ClassDef(AliAnalysisTaskSEBkgLikeSignJPSI,0); // comparison of unlike-sign and like-sign background for J/psi->ee
};

#endif

