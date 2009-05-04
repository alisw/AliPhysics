#ifndef ALIANALYSISTASKSEBTOJPSITOELE_H
#define ALIANALYSISTASKSEBTOJPSITOELE_H
/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
//               Class AliAnalysisTaskSEBtoJPSItoEle
//                AliAnalysisTaskSE for the reconstruction 
//                   of heavy-flavour decay candidates
//            Author: C.Di Giglio, carmelo.digiglio@ba.infn.it
//*************************************************************************

class TH1F;
class TList;
class AliAnalysisBtoJPSItoEle;
#include <TClonesArray.h>
#include <TNtuple.h>
#include <TH2F.h> 
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSEBtoJPSItoEle : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSEBtoJPSItoEle();
  AliAnalysisTaskSEBtoJPSItoEle(const char *name);
  virtual ~AliAnalysisTaskSEBtoJPSItoEle();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void SetCutsD0(const Double_t cutsD0[9]);
  void SetPtCuts(const Double_t ptcuts[2]);
  void SetAODMCInfo(Bool_t OkMCInfo) { fOkAODMC = OkMCInfo;}
  void SetMinimize(Bool_t OkMinimize) { fOkMinimize = OkMinimize;}
  void ReadAODMCInfo(const AliAODEvent* aodEv, const TClonesArray* inArray);
  void SetPathMCfile (TString mcfilename="CsiMCfromKine_10PtBins.root") {fNameMCfile = mcfilename;}
  TH1F* OpenLocalMCFile(TString dataDir, Int_t nbin);  

 private:

  AliAnalysisTaskSEBtoJPSItoEle(const AliAnalysisTaskSEBtoJPSItoEle &source);
  AliAnalysisTaskSEBtoJPSItoEle& operator=(const AliAnalysisTaskSEBtoJPSItoEle& source); 
  //
  TList *fOutput;                            //! list send on output slot 0
  AliAnalysisBtoJPSItoEle *fCdfFit;          //! Unbinned log-likelihood minimizer 
  TNtuple *fNtupleJPSI;                      //! Ntuple of pseudo-proper decay time & invariant mass values
  TH1F *fhDecayTimeMCjpsifromB;              //! Pseudo-proper decay time distribution used as template for JPSIs from B
  TH1F *fhDecayTime;                         //! Pseudo-proper decay time distribution
  TH1F *fhInvMass;                           //! Invariant mass distribution
  TH1F *fhD0;                                //! Impact parameter distribution
  TH1F *fhD0D0;                              //! Product of impact parameters distributions
  TH1F *fhCosThetaStar;                      //! Cosine of decay angle distribution
  TH1F *fhCosThetaPointing;                  //! Cosine of pointing angle distribution
  TH2F *fhCtsVsD0D0;                         //! Cos theta star Vs. D0D0 distribution
  Bool_t fOkAODMC;                           // Flag to read AOD monte carlo information
  TString fNameMCfile;                       // Name of the MC file with X template
  Bool_t fOkMinimize;                        // Flag to enable unbinned log-likelihood minimization
   
  Double_t fCuts[9];                         // cuts for N-tuple values selection
  Double_t fPtCuts[2];                       // Max and min pt of the candidates

  ClassDef(AliAnalysisTaskSEBtoJPSItoEle,0); // AliAnalysisTaskSE for the reconstruction of heavy-flavour decay candidates
};
#endif
