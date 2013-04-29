#ifndef ALIANALYSISTASKSEDMESONSFILTERCJ_H
#define ALIANALYSISTASKSEDMESONSFILTERCJ_H

// $Id$

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskSEDmesonsCJ
// AliAnalysisTaskSE for Dmesons - jet correlations analysis
// Author: Xiaoming Zhang, xmzhang@lbl.gov
//*************************************************************************

#include "AliAnalysisTaskSE.h"

class TList;
class AliRDHFCutsD0toKpi;
class AliRDHFCutsDStartoKpipi;

class AliAnalysisTaskSEDmesonsFilterCJ : public AliAnalysisTaskSE {

 public :

  AliAnalysisTaskSEDmesonsFilterCJ();
  AliAnalysisTaskSEDmesonsFilterCJ(const char *name);
  virtual ~AliAnalysisTaskSEDmesonsFilterCJ();

  virtual void Init();
  virtual void LocalInit() { Init(); }
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *opt);
  virtual void Terminate(Option_t *opt);
  virtual void NotifyRun();

  Bool_t SetCutDzero(AliRDHFCutsD0toKpi      *cut);
  Bool_t SetCutDstar(AliRDHFCutsDStartoKpipi *cut);

 private :

  AliAnalysisTaskSEDmesonsFilterCJ(const AliAnalysisTaskSEDmesonsFilterCJ &);
  AliAnalysisTaskSEDmesonsFilterCJ& operator=(const AliAnalysisTaskSEDmesonsFilterCJ &);

  Bool_t ExecOnce();
  Bool_t ExecEach();

  void ExecAnas();
  void ExecAnasDzero();
  void ExecAnasDstar();

  void CreateDzeroHistograms();
  void CreateDstarHistograms();
//=============================================================================

  AliAODEvent *fEventAOD;  //! in put AOD event

  AliRDHFCutsD0toKpi      *fCutDzero;  //! Dzero selection cut
  AliRDHFCutsDStartoKpipi *fCutDstar;  //! Dstar selection cut

  TClonesArray *fDzeroClArr;  //! intput Dzero candidates array
  TClonesArray *fDstarClArr;  //! intput Dstar candidates array

  TClonesArray *fUsedDzeros;  //! selected Dzero candidates array
  TClonesArray *fUsedD0bars;  //! selected D0bar candidates array
  TClonesArray *fUsedDstars;  //! selected Dstar candidates array

  Bool_t fIsNotExecOnce;  // flag for ExecOnce

  TList *fListCutsDmesons;  //! list of Dmeson selection cuts
  TList *fListDzeroHistos;  //! list of Dzero output control histrograms
  TList *fListDstarHistos;  //! list of Dstar output control histrograms

  ClassDef(AliAnalysisTaskSEDmesonsFilterCJ, 1);
};

#endif
