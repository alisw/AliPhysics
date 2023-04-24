/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskFlatenicityLambdaK0sMC_H
#define AliAnalysisTaskFlatenicityLambdaK0sMC_H
#include "AliEventCuts.h"
#include "AliAnalysisTaskSE.h"
class AliESDtrackCuts;
class AliMCEventHandler;
class AliMCEvent;
class AliStack;
class AliESDAD;
class AliESDEvent;
class AliAnalysisUtils;
class AliPIDResponse;
class AliAnalysisTaskFlatenicityLambdaK0sMC : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskFlatenicityLambdaK0sMC();
  AliAnalysisTaskFlatenicityLambdaK0sMC(const char *name);
  virtual ~AliAnalysisTaskFlatenicityLambdaK0sMC();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  Double_t MyRapidity(Double_t rE, Double_t rPz) const;
  Float_t GetDCAz(AliESDtrack *lTrack);
  Float_t GetFlatenicityV0();
  void CheckChargeV0(AliESDv0 *v0);
  AliEventCuts fEventCuts;

private:
  AliESDEvent *fESD;  //! input event
  TList *fOutputList; //! output list
  TList *fOutputListMC; //! output list

  AliPIDResponse *fPIDResponse;
  AliESDpid *fESDpid;
  TH1F *hinvmassK0s;        //! dummy histogram
  TH1F *hinvmassLambda;     //! dummy histogram
  TH1F *hinvmassAntiLambda; //! dummy histogram
  TH1F *hflat;              //! dummy histogram
//  AliAnalysisFilter *fTrackFilter;
  Float_t invmK0s;
  Float_t invpK0s;
  Float_t invptK0s;
  Float_t invyK0s;
  Float_t flatenicityK0s;
  AliESDtrackCuts *fESDtrackCuts; // ESD track cuts used for primary track definition

  TTree *treeK0s;

  Float_t invmLambda;
  Float_t invpLambda;
  Float_t invptLambda;
  Float_t invyLambda;
  Float_t flatenicityLambda;

  TTree *treeLambda;

  Float_t invmAntiLambda;
  Float_t invpAntiLambda;
  Float_t invptAntiLambda;
  Float_t invyAntiLambda;
  Float_t flatenicityAntiLambda;

  TTree *treeAntiLambda;

  //MC True
    TH1F *hinvmassK0sMC;        //! dummy histogram
  TH1F *hinvmassLambdaMC;     //! dummy histogram
  TH1F *hinvmassAntiLambdaMC; //! dummy histogram
  TH1F *hNeventMC;              //! dummy histogram
//  AliAnalysisFilter *fTrackFilter;
  Float_t invmK0sMC;
  Float_t invpzK0sMC;
  Float_t invptK0sMC;
  Float_t invyK0sMC;
  Float_t flatenicityK0sMC;
  AliESDtrackCuts *fESDtrackCutsMC; // ESD track cuts used for primary track definition

  TTree *treeK0sMC;

  Float_t invmLambdaMC;
  Float_t invpzLambdaMC;
  Float_t invptLambdaMC;
  Float_t invyLambdaMC;
  Float_t flatenicityLambdaMC;

  TTree *treeLambdaMC;

  Float_t invmAntiLambdaMC;
  Float_t invpzAntiLambdaMC;
  Float_t invptAntiLambdaMC;
  Float_t invyAntiLambdaMC;
  Float_t flatenicityAntiLambdaMC;

  TTree *treeAntiLambdaMC;

  AliAnalysisTaskFlatenicityLambdaK0sMC(
      const AliAnalysisTaskFlatenicityLambdaK0sMC &); // not implemented
  AliAnalysisTaskFlatenicityLambdaK0sMC &
  operator=(const AliAnalysisTaskFlatenicityLambdaK0sMC &); // not implemented

  ClassDef(AliAnalysisTaskFlatenicityLambdaK0sMC, 1);
};

#endif
