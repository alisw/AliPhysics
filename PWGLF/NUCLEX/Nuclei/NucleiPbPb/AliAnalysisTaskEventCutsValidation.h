#ifndef __AliAnalysisTaskEventCutsValidation__
#define __AliAnalysisTaskEventCutsValidation__

#include <TF1.h>

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliESDtrackCuts.h"

class TList;
class TH2F;

class AliAnalysisTaskEventCutsValidation : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEventCutsValidation(bool storeCuts = false, TString taskname = "EventCutsValidation");
  virtual ~AliAnalysisTaskEventCutsValidation();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);

  AliEventCuts  fEventCut;    ///<  Event cuts

private:
  AliAnalysisTaskEventCutsValidation (const AliAnalysisTaskEventCutsValidation &source);
  AliAnalysisTaskEventCutsValidation &operator=(const AliAnalysisTaskEventCutsValidation &source);

  TList *fList;            //!<! Output list

  bool                fStoreCuts; ///<  If true the EventCuts object is put in the output
  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskEventCutsValidation, 1);
  /// \endcond
};


#endif /* defined(__AliAnalysisTaskEventCutsValidation__) */
