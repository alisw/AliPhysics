#ifndef ALIANALYSISTASKSED0MASS_H
#define ALIANALYSISTASKSED0MASS_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskSED0Mass
// AliAnalysisTaskSE for D0 candidates invariant mass histogram
// and comparison to MC truth (kinematics stored in the AOD) and cut variables
// distributions
// Authors: A.Dainese, andrea.dainese@ln.infn.it
// and C.Bianchin, chiara.bianchin@pd.infn.it
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1F.h>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"

class AliAnalysisTaskSED0Mass : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSED0Mass();
  AliAnalysisTaskSED0Mass(const char *name);
  virtual ~AliAnalysisTaskSED0Mass();


  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  
 private:

  AliAnalysisTaskSED0Mass(const AliAnalysisTaskSED0Mass &source);
  AliAnalysisTaskSED0Mass& operator=(const AliAnalysisTaskSED0Mass& source); 
  void    FillHists(Int_t ptbin, AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliAnalysisVertexingHF *vhf, TList *listout);
  TList   *fOutputtight; //! list send on output slot 1
  TList   *fOutputloose; //! list send on output slot 2
  TList   *fDistr;       //! list send on output slot 3
  AliAnalysisVertexingHF *fVHFtight;  // Vertexer heavy flavour (used to pass the cuts)
  AliAnalysisVertexingHF *fVHFloose;  // Vertexer heavy flavour (used to pass the cuts)
  TH1F *fNentries; //! histogram with number of events
  ClassDef(AliAnalysisTaskSED0Mass,2); // AliAnalysisTaskSE for the MC association of heavy-flavour decay candidates
};

#endif
