#ifndef AliAnalysisTaskNanoAODnormalisation_H
#define AliAnalysisTaskNanoAODnormalisation_H

#include "AliAnalysisCuts.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisNanoAODCuts.h"

#include <TList.h>
#include <TString.h>

class AliAnalysisTaskNanoAODnormalisation : public AliAnalysisTaskSE {
  public:
  AliAnalysisTaskNanoAODnormalisation(std::string taskName = "AliAnalysisTaskNanoAODnormalisation");
  virtual ~AliAnalysisTaskNanoAODnormalisation();

  virtual void UserCreateOutputObjects();
  virtual Bool_t UserNotify();
  virtual void   UserExec(Option_t *) {}

  void SetMultBinning(int nbins, float min, float max) {
    fNmultBins = nbins;
    fMinMult = min;
    fMaxMult = max;
  }
  static AliAnalysisTaskNanoAODnormalisation* AddTask(std::string name = "AliAnalysisTaskNanoAODnormalisation");

  private:
  void FillHistograms(TH2D* candidate[2], TH2D* selected[2]);

  TList*  fOutputList;                        //!<! Output list                    

  int   fNmultBins; ///
  float fMinMult;   ///
  float fMaxMult;   ///

  TH2D* fCandidateEvents[2]; //!<!
  TH2D* fSelectedEvents[2];  //!<!
  TH2D* fCandidateEventsUE[2]; //!<!
  TH2D* fSelectedEventsUE[2];  //!<!

  ClassDef(AliAnalysisTaskNanoAODnormalisation,1);
};

#endif
