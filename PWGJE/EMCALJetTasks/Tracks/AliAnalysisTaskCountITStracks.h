#ifndef ALIANALYSISTASKCOUNTITSTRACKS_H
#define ALIANALYSISTASKCOUNTITSTRACKS_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class THistManager;
class AliAnalysisUtils;

#include "AliAnalysisTaskSE.h"

namespace PWGJE {
  
namespace EMCALJetTasks {

class AliAnalysisTaskCountITStracks : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskCountITStracks();
  AliAnalysisTaskCountITStracks(const char *name);
  virtual ~AliAnalysisTaskCountITStracks();

  static AliAnalysisTaskCountITStracks *AddTaskCountITStracks(const char *name);

protected:
  void UserCreateOutputObjects();
  void UserExec(Option_t *);

private:
  AliAnalysisUtils          *fAnalysisUtils;          //!<!
  THistManager              *fHistos;                 //!<!

  ClassDef(AliAnalysisTaskCountITStracks, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIANALYSISTASKCOUNTITSTRACKS_H */
