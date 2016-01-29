/*
 * AliAnalysisTaskEventFilter.h
 *
 *  Created on: Jan 29, 2016
 *      Author: markus
 */

#ifndef ALIANALYSISTASKEVENTFILTER_H
#define ALIANALYSISTASKEVENTFILTER_H

#include <AliAnalysisTaskSE.h>

class THistManager;

class AliAnalysisUtils;
class AliVEvent;

namespace EMCalTriggerPtAnalysis {

class AliAnalysisTaskEventFilter: public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEventFilter();
  AliAnalysisTaskEventFilter(const char *name);
  virtual ~AliAnalysisTaskEventFilter();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);

protected:
  Bool_t FakeVertexSelection2013pA(const AliVEvent * const inputevent) const;

private:
  AliAnalysisUtils          *fAnalysisUtils;
  THistManager              *fHistos;

  ClassDef(AliAnalysisTaskEventFilter, 1);
};

}

#endif /* ALIANALYSISTASKEVENTFILTER_H */
