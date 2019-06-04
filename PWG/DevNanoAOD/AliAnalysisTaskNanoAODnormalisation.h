#ifndef AliAnalysisTaskNanoAODnormalisation_H
#define AliAnalysisTaskNanoAODnormalisation_H

#include "AliAnalysisCuts.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisNanoAODCuts.h"

#include <TList.h>

class AliAnalysisTaskNanoAODnormalisation : public AliAnalysisTaskSE {
  public:
  AliAnalysisTaskNanoAODnormalisation(std::string taskName = "AliAnalysisTaskNanoAODnormalisation");
  virtual ~AliAnalysisTaskNanoAODnormalisation();

  virtual void UserCreateOutputObjects();
  virtual Bool_t UserNotify();

  static AliAnalysisTaskNanoAODnormalisation* AddTask(std::string name = "AliAnalysisTaskNanoAODnormalisation");

  private:
  bool fFirst;                   //!<! Current filke name
  TList*  fOutputList;                        //!<! Output list                    

  TH2D* fCandidateEvents[2]; //!<!
  TH2D* fSelectedEvents[2];  //!<!

  ClassDef(AliAnalysisTaskNanoAODnormalisation,1);
};

#endif
