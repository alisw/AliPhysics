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
#ifndef AliAnalysisTaskResonanceDP_H
#define AliAnalysisTaskResonanceDP_H

#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "TTree.h"
#include "AliEventCuts.h"

class TList;
class TTree;

class AliESDEvent;
class AliESDtrackCuts;
class AliMCEvent;

class TH1F;
class TH2F;
class TH3F;
class TProfile;
class AliPIDResponse;
class AliMultSelection;


class AliAnalysisTaskResonanceDP : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskResonanceDP();
  AliAnalysisTaskResonanceDP(const char *name);
  virtual ~AliAnalysisTaskResonanceDP();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  void           SetPersonalESDtrackCuts(AliESDtrackCuts* trackcuts);
  Double_t       MyRapidity(Double_t rE, Double_t rPz) const;
  //---------------------------------------------------------------------------------------
  Double_t      GetTOFBeta(AliVTrack *esdtrack);
  Bool_t        MatchTOF(AliVTrack *vtrack);
  void setTriggerType(UInt_t type) {fTriggerMask=type;}
  
 private:
  enum
  {
    kMaxTrack=500
  };
  
  TTree  *fTreeEvent;//!
  AliPIDResponse   *fPIDResponse;//!
  AliESDtrackCuts  *fESDtrackCuts;//!
  AliEventCuts fEventCuts; //!
  UInt_t fTriggerMask;
  Float_t fTreeTrackVariableCentrality;
  Float_t fTreeTrackVariableVtxz;
  Int_t fTreeTrackVariableNTrack;//!
  Float_t fTreeTrackVariabledeuteronnsigmaTPC[kMaxTrack];//!
  Float_t fTreeTrackVariableprotonnsigmaTPC[kMaxTrack];//!
  Float_t fTreeTrackVariablemasssquare[kMaxTrack];//!
  Int_t fTreeTrackVariableCharge[kMaxTrack];//!
  Float_t fTreeTrackVariableDCAXY[kMaxTrack];//!
  Float_t fTreeTrackVariableDCAZ[kMaxTrack];//!
  Float_t fTreeTrackVariableITSchi2[kMaxTrack];//!
  Float_t fTreeTrackVariableTPCchi2[kMaxTrack];//!
  Float_t fTreeTrackVariableNCR[kMaxTrack];//!
  Float_t fTreeTrackVariableMomentumPx[kMaxTrack];//!
  Float_t fTreeTrackVariableMomentumPy[kMaxTrack];//!
  Float_t fTreeTrackVariableMomentumPz[kMaxTrack];//!
  AliAnalysisTaskResonanceDP(const AliAnalysisTaskResonanceDP&);
  AliAnalysisTaskResonanceDP& operator=(const AliAnalysisTaskResonanceDP&);  
  ClassDef(AliAnalysisTaskResonanceDP, 1);
};
#endif
