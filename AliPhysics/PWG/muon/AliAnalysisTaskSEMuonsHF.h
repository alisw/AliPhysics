#ifndef ALIANALYSISTASKSEMUONSHF_H
#define ALIANALYSISTASKSEMUONSHF_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
// Class AliAnalysisTaskSEMuonsHF
// AliAnalysisTaskSE for the single muon and dimuon from HF analysis
// Author: X-M. Zhang, zhang@clermont.in2p3.fr
//                     zhangxm@iopp.ccnu.edu.cn
//*************************************************************************

#include "AliAnalysisTaskSE.h"

class TString;
class TList;
class TClonesArray;
class AliMuonsHFHeader;
class AliMuonTrackCuts;
class AliMuonPairCuts;

class AliAnalysisTaskSEMuonsHF : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskSEMuonsHF();
  AliAnalysisTaskSEMuonsHF(const char *name, const AliMuonTrackCuts& cutsMuon, const AliMuonPairCuts& cutsDimu);
  virtual ~AliAnalysisTaskSEMuonsHF();

  virtual void Init();
  virtual void LocalInit() { Init(); }
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *opt);
  virtual void Terminate(Option_t *opt);
  virtual void NotifyRun();

  void SetAnaMode(Int_t mode)      { fAnaMode      = ((mode>=0 && mode<=2) ? mode : 0); }
  void SetIsOutputTree(Bool_t ist) { fIsOutputTree = ist;                               }
  void SetUseMC(Bool_t isMC)       { fIsMC         = isMC;                              }
  void SetIsFull(Bool_t isFull)    { fIsFull       = isFull;                            }

  void SetEvsHCuts(Double_t cuts[5]) const;

 private:

  AliAnalysisTaskSEMuonsHF(const AliAnalysisTaskSEMuonsHF&);
  AliAnalysisTaskSEMuonsHF& operator=(const AliAnalysisTaskSEMuonsHF&);

  Int_t fAnaMode;       // = 0, ana both single muon and dimuon
                        // = 1, ana single muon
                        // = 2, ana dimuon
  Bool_t fIsOutputTree; // flag used to switch on/off tree output
  Bool_t fIsMC;         // flag to use MC
  Bool_t fIsFull;       // flag to save the parton info in MC (PYTHIA)

  AliMuonTrackCuts *fCutsMuon; // single muon selection cuts
  AliMuonPairCuts  *fCutsDimu; // dimuon selection cuts

  AliMuonsHFHeader *fHeader;     // output for info at ev level
  TClonesArray     *fMuonClArr;  // output clones array for single mu
  TClonesArray     *fDimuClArr;  // output clones array for dimu
  TList            *fListOutput; // output list of histos

  ClassDef(AliAnalysisTaskSEMuonsHF, 8);
};

#endif
