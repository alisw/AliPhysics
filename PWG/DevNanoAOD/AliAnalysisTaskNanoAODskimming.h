#ifndef ALIANALYSISTASKNanoAODFILTERSELECTION_H
#define ALIANALYSISTASKNanoAODFILTERSELECTION_H

#include "AliAnalysisCuts.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisNanoAODCuts.h"

#include <TList.h>

class AliNanoFilterNormalisation;

/// The NanoAOD selection perform an offline triggering using a combination of event selection
/// and PID triggering on the tracks. If the trigger conditions are not satisfied, then the
/// following tasks are not executed
class AliAnalysisTaskNanoAODskimming : public AliAnalysisTaskSE {
  public:
  AliAnalysisTaskNanoAODskimming(std::string taskName = "NanoAODSkimming");
  virtual ~AliAnalysisTaskNanoAODskimming();

  virtual void UserExec(Option_t* /*option*/);
  virtual void UserCreateOutputObjects();
  virtual void Terminate(Option_t *);
  virtual void FinishTaskOutput();

  void AddEventCut(AliAnalysisCuts* cut) { fOtherEventCuts.Add(cut); }
  static AliAnalysisTaskNanoAODskimming* AddTask(std::string name = "NanoAODSkimming");

  AliAnalysisNanoAODEventCuts fEventCuts;

  private:
  AliNanoFilterNormalisation* fNormalisation;          //!<! Normalisation object
  TList fOtherEventCuts;                              /// List of cuts to be applied

  ClassDef(AliAnalysisTaskNanoAODskimming,1);
};

#endif
