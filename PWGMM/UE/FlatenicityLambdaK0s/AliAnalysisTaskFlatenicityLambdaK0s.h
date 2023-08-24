/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskFlatenicityLambdaK0s_H
#define AliAnalysisTaskFlatenicityLambdaK0s_H
#include "AliEventCuts.h"
#include "AliAnalysisTaskSE.h"
class AliESDtrackCuts;
class AliESDAD;
class AliESDEvent;
class AliAnalysisUtils;
class AliPIDResponse;
class AliAnalysisTaskFlatenicityLambdaK0s : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskFlatenicityLambdaK0s();
  AliAnalysisTaskFlatenicityLambdaK0s(const char *name);
  virtual ~AliAnalysisTaskFlatenicityLambdaK0s();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  Float_t GetDCAz(AliESDtrack *lTrack);
  Float_t GetFlatenicityV0();
  void CheckChargeV0(AliESDv0 *v0);
  AliEventCuts fEventCuts;

private:
  AliESDEvent *fESD;  //! input event
  TList *fOutputList; //! output list
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

  AliAnalysisTaskFlatenicityLambdaK0s(
      const AliAnalysisTaskFlatenicityLambdaK0s &); // not implemented
  AliAnalysisTaskFlatenicityLambdaK0s &
  operator=(const AliAnalysisTaskFlatenicityLambdaK0s &); // not implemented

  ClassDef(AliAnalysisTaskFlatenicityLambdaK0s, 1);
};

#endif
