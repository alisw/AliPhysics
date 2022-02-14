/*
 * AliAnalysisTaskEtaPhiEfficiency.h
 *
 *  Created on: Feb 4, 2016
 *      Author: markus
 */

#ifndef ALIANALYSISTASKETAPHIEFFICIENCY_H
#define ALIANALYSISTASKETAPHIEFFICIENCY_H

#include "AliAnalysisTaskSE.h"

class THistManager;
class AliAnalysisUtils;
class AliESDtrackCuts;

namespace PWGJE {
  
namespace EMCALJetTasks {

class AliAnalysisTaskEtaPhiEfficiency : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEtaPhiEfficiency();
  AliAnalysisTaskEtaPhiEfficiency(const char *name);
  virtual ~AliAnalysisTaskEtaPhiEfficiency();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual Bool_t UserNotify();

private:
  Bool_t PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials, Int_t &pthard) const;
  Bool_t FindPtBin(Double_t ptin, Int_t &ptmin, Int_t &ptmax) const;

  AliAnalysisUtils        *fAnalysisUtil;
  THistManager            *fHistos;
  AliESDtrackCuts         *fTrackCuts;

  AliAnalysisTaskEtaPhiEfficiency(const AliAnalysisTaskEtaPhiEfficiency &);
  AliAnalysisTaskEtaPhiEfficiency &operator=(const AliAnalysisTaskEtaPhiEfficiency &);

  ClassDef(AliAnalysisTaskEtaPhiEfficiency, 1);
};

}

}

#endif /* ALIANALYSISTASKETAPHIEFFICIENCY_H */
