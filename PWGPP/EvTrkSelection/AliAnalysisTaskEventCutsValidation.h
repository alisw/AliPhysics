#ifndef __AliAnalysisTaskEventCutsValidation__
#define __AliAnalysisTaskEventCutsValidation__

#include <TF1.h>

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliESDtrackCuts.h"

class TList;
class TH2F;
class TTree;

class AliAnalysisTaskEventCutsValidation : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEventCutsValidation(bool storeCuts = false, TString taskname = "EventCutsValidation");
  virtual ~AliAnalysisTaskEventCutsValidation();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);

  bool          fFillTree;
  AliEventCuts  fEventCut;    ///<  Event cuts

  struct EventSummary {
    ULong64_t trigger;
    float x,y,z,v0m,cl0,fb32,fb32acc,fb32tof,tpc,tpcOut,esd,multvzero;
  };

private:
  AliAnalysisTaskEventCutsValidation (const AliAnalysisTaskEventCutsValidation &source);
  AliAnalysisTaskEventCutsValidation &operator=(const AliAnalysisTaskEventCutsValidation &source);

  TList *fList;            //!<! Output list
  TTree *fTree;            //!<! Output list

  EventSummary        fCurrentEvent;
  bool                fStoreCuts; ///<  If true the EventCuts object is put in the output
  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskEventCutsValidation, 1);
  /// \endcond
};


#endif /* defined(__AliAnalysisTaskEventCutsValidation__) */
