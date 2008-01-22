#ifndef AliCutTask_cxx
#define AliCutTask_cxx

// simple task that runs the esd track cuts to evaluate the basic plots created during the cuts

class TH1F;
class AliESDtrackCuts;
class AliESDEvent;
class TList;

#include "AliAnalysisTask.h"
#include "AliPWG0Helper.h"

class AliCutTask : public AliAnalysisTask {
 public:
  AliCutTask(const char *name = "AliCutTask");
  virtual ~AliCutTask() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetTrackCuts(AliESDtrackCuts* cuts) { fTrackCuts = cuts; }
  void SetAnalysisMode(AliPWG0Helper::AnalysisMode mode) { fAnalysisMode = mode; }
  void EnableSecondaryStudy();
  
 private:
  AliESDEvent *fESD;           //! ESD object
  AliESDtrackCuts* fTrackCuts; // track cuts
  AliPWG0Helper::AnalysisMode fAnalysisMode; // detector that is used for analysis
  
  AliESDtrackCuts* fTrackCutsPrimaries; // cuts for tracks from primary particles
  AliESDtrackCuts* fTrackCutsSecondaries; // cuts for tracks from secondary particles

  TH1F* fVertex;   //! event z vertex distribution
  TH1F* fTriggerStats;  //! triggers

  TList* fOutput;                  //! list send on output slot 0

  TH1F* fPrimStats; //! statistics about primaries, see bin names in CreateOutputData

  AliCutTask(const AliCutTask&); // not implemented
  AliCutTask& operator=(const AliCutTask&); // not implemented
  
  ClassDef(AliCutTask, 1); // example of analysis
};

#endif
