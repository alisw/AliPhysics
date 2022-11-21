/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskFlatenicityLambdaK0s_H
#define AliAnalysisTaskFlatenicityLambdaK0s_H

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskFlatenicityLambdaK0s : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskFlatenicityLambdaK0s();
  AliAnalysisTaskFlatenicityLambdaK0s(const char *name);
  virtual ~AliAnalysisTaskFlatenicityLambdaK0s();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  void SetTrackCuts(AliESDtrackCuts *myTracksCuts) { fTracksCuts = myTracksCuts; }

private:
  AliAODEvent *fAOD;  //! input event
  AliESDEvent *fESD;  //! input event
  TList *fOutputList; //! output list
  AliESDpid *fESDpid;
  TH1F *hinvmassK0s;        //! dummy histogram
  TH1F *hinvmassLambda;     //! dummy histogram
  TH1F *hinvmassAntiLambda; //! dummy histogram
  TH1F *hflat;              //! dummy histogram
  AliPIDResponse *fPIDResponse;
  AliESDtrackCuts *fTracksCuts;
  AliAnalysisTaskFlatenicityLambdaK0s(const AliAnalysisTaskFlatenicityLambdaK0s &);            // not implemented
  AliAnalysisTaskFlatenicityLambdaK0s &operator=(const AliAnalysisTaskFlatenicityLambdaK0s &); // not implemented
  Double_t GetFlatenicityV0();
  ClassDef(AliAnalysisTaskFlatenicityLambdaK0s, 1);
};

#endif
