#ifndef ALIANALYSISTASKSEDSTAR_H
#define ALIANALYSISTASKSEDSTAR_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskSEDStar
// AliAnalysisTaskSE for D* candidates invariant mass histogram
// and comparison to MC truth (kinematics stored in the AOD)
// Authors: Y.wang, yifei@physi.uni-heidelberg.de
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TH1F.h>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"

class AliAnalysisTaskSEDStar : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSEDStar();
  AliAnalysisTaskSEDStar(const char *name);
  virtual ~AliAnalysisTaskSEDStar();


  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void	       SetReadMC(Bool_t readMC=kFALSE){fReadMC=readMC;}
  void         PIDon() {fPID=kTRUE;}
  void         PIDoff() {fPID=kFALSE;}
  const Bool_t       GetPIDStatus() {return fPID;}
  void         SetNSigmaTPC(Double_t nsigma) {fNSigma=nsigma;}
  const Double_t     GetNSigmaTPC() {return fNSigma;}
 private:

  AliAnalysisTaskSEDStar(const AliAnalysisTaskSEDStar &source);
  AliAnalysisTaskSEDStar& operator=(const AliAnalysisTaskSEDStar& source); 
  TList   *fOutput; //! list send on output slot 0

  TH1F    *fhistMass;  //! output total invariant mass histogram - no MC truth
  TH1F    *fhistSgn;  //! output signal invariant mass histogram - use cuts
  TH1F    *fhistBkg;  //! output background invariant mass histogram - use cuts
  Bool_t   fReadMC;  // flag for MC array: kTRUE = read it, kFALSE = do not read it
  Bool_t   fPID;  //flag of TPCpid selection
  Double_t fNSigma;   //nsigma selection for TPCpid

  AliAnalysisVertexingHF *fVHF;  // Vertexer heavy flavour (used to pass the cuts)
  Bool_t SelectTPCPID(AliAODTrack *trk, Int_t pid, Double_t nsig);

  ClassDef(AliAnalysisTaskSEDStar,1); // AliAnalysisTaskSE for the MC association of heavy-flavour decay candidates
};

#endif

