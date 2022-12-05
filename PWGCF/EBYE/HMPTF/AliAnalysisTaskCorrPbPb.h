/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#ifndef AliAnalysisTaskCorrPbPb_H
#define AliAnalysisTaskCorrPbPb_H

#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "TTree.h"
#include "AliEventCuts.h"

class TList;
class TTree;

class AliESDEvent;
//class AliAODEvent;
class AliESDtrackCuts;
class AliMCEvent;
class AliStack;
class AliVTrack;
class AliAODTrack;

class TH1D;
class TH2D;
class TH3D;
class TProfile;
class AliPIDResponse;
class AliMultSelection;


class AliAnalysisTaskCorrPbPb : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskCorrPbPb();
  AliAnalysisTaskCorrPbPb(const char *name);
  virtual ~AliAnalysisTaskCorrPbPb();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  Bool_t GetEvent();
  Bool_t PassedTrackQualityCuts (AliAODTrack *track);
  Bool_t KaonSelector (AliVTrack *track); 
  Bool_t ProtonSelector (AliVTrack *track); 
  Bool_t PionSelector (AliVTrack *track); 
  Bool_t PassedPIDSelection (AliAODTrack *track, AliPID::EParticleType type);
  Bool_t PassedSingleParticlePileUpCuts(AliAODTrack *track);

  /*
  void SetVzRangeMax(Double_t VzMax)
  {
    this->fVertexZMax = VzMax;
  }
  void SetDCAXYRangeMax(Double_t dcaxy)          
  {
    fDCAxyMax = dcaxy;
  }
  void SetDCAZRangeMax(Double_t dcaz)
  {
    fDCAzMax = dcaz;
  }
  void SetMaxChi2PerTPCClusterRange(Double_t chi2tpc)
  {
    fChi2TPC = chi2tpc;
  }
  void SetMaxChi2PerITSClusterRange(Double_t chi2its)
  {
    fChi2ITS = chi2its;
  }
  void SetMinNCrossedRowsTPCRange(Double_t nCrossedRow)
  {
    fNCrossedRowsTPC = nCrossedRow;
  }
  void SetTreeName(TString TreeName)
  {
    fTreeName = TreeName;
  }
  void SetEtaCut(Double_t EtaMax)
  {
    fEtaMax=EtaMax;
  }
  */

  
 private:
  
  
  AliESDEvent *fESDevent;
  AliAODEvent *fAODevent;
  AliVEvent *fInputEvent;
  AliPIDResponse   *fPIDResponse;
  AliESDtrackCuts  *fESDtrackCuts;
  AliESDtrackCuts  *fESDtrackCuts_primary;
  AliAnalysisUtils *fUtils;
  AliEventCuts fAODeventCuts;
  TList *fOutputList;
  TList *fQAList;
  Double_t fMultLow;
  Double_t fMultHigh;
  TH1D *hNumberOfEvents;
  TH1D *hNumberOfKaonPlus;
  TH1D *hNumberOfKaonMinus;
  TH1D *hNumberOfPionPlus;
  TH1D *hNumberOfPionMinus;
  TH1D *hNumberOfProtonPlus;
  TH1D *hNumberOfProtonMinus;

  //ftreeEvent object variables
  TTree *fTreeEvent;
  Float_t fTreeVariableCentrality;
  //reconstruced
  Float_t fNoKaonPlus_ptmax2;
  Float_t fNoKaonMinus_ptmax2;
  Float_t fNoKaonPlus_ptmax3;
  Float_t fNoKaonMinus_ptmax3;
  Float_t fNoPionPlus_ptmax2;
  Float_t fNoPionMinus_ptmax2;
  Float_t fNoPionPlus_ptmax3;
  Float_t fNoPionMinus_ptmax3;
  Float_t fNoProtonPlus_ptmax2;
  Float_t fNoProtonMinus_ptmax2;
  Float_t fNoProtonPlus_ptmax3;
  Float_t fNoProtonMinus_ptmax3;


  /*
  //Custom Functions:
  TF1 *fCenCutLowPU;
  TF1 *fCenCutHighPU;
  TF1 *fSPDCutPU;
  TF1 *fV0CutPU;
  TF1 *fMultCutPU;

  //Argument variables
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  TString fTreeName;
  Double_t fEtaMax;
  Double_t fVertexZMax;
  Double_t fDCAxyMax;
  Double_t fDCAzMax;
  Double_t fChi2TPC;
  Double_t fChi2ITS;
  Double_t fNCrossedRowsTPC;
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
 
  
  AliAnalysisTaskCorrPbPb(const AliAnalysisTaskCorrPbPb&);
  AliAnalysisTaskCorrPbPb& operator=(const AliAnalysisTaskCorrPbPb&);  
  ClassDef(AliAnalysisTaskCorrPbPb, 1);
};

#endif
