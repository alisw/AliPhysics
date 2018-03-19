/*
 * AliAnalysisTaskMCFemtoDream.h
 *
 *  Created on: Mar 15, 2018
 *      Author: hohlweger
 */

#ifndef ALIANALYSISTASKMCFEMTODREAM_H_
#define ALIANALYSISTASKMCFEMTODREAM_H_
#include "Rtypes.h"
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamBasePart.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamCollConfig.h"
#include "TList.h"
#include "TH1F.h"
class AliAnalysisTaskMCFemtoDream : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskMCFemtoDream();
  AliAnalysisTaskMCFemtoDream(const char *name);
  virtual ~AliAnalysisTaskMCFemtoDream();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *){};
 private:
  TList *foutput;
  AliFemtoDreamBasePart *fPart;
  AliFemtoDreamCollConfig *fcollConfig;
  AliFemtoDreamPartCollection *fPartCollAll;   //!
  AliFemtoDreamPartCollection *fPartCollAccept;   //!
  TH1F *fCounter;
  TH1F *fProtonEta;
  TH1F *fProtonPt;
  TH1F *fProtonResonances;
  TH1F *fAntiProtonEta;
  TH1F *fAntiProtonPt;
  TH1F *fAntiProtonResonances;
  TH1F *fLambdaEta;
  TH1F *fLambdaPt;
  TH1F *fLambdaResonances;
  TH1F *fAntiLambdaEta;
  TH1F *fAntiLambdaPt;
  TH1F *fAntiLambdaResonances;
  TH1F *fXiEta;
  TH1F *fXiPt;
  TH1F *fXiResonances;
  TH1F *fAntiXiEta;
  TH1F *fAntiXiPt;
  TH1F *fAntiXiResonances;
  ClassDef(AliAnalysisTaskMCFemtoDream,1);
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKMCFEMTODREAM_H_ */
