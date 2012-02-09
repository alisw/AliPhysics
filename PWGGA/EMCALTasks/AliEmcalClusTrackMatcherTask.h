#ifndef ALIEMCALCLUSTRACKMATCHERTASK_H
#define ALIEMCALCLUSTRACKMATCHERTASK_H

// $Id$

class TClonesArray;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliEmcalClusTrackMatcherTask : public AliAnalysisTaskSE {
 public:
  AliEmcalClusTrackMatcherTask(const char *name=0);
  virtual ~AliEmcalClusTrackMatcherTask();

  void         UserCreateOutputObjects();
  void         UserExec(Option_t *option);
  void         Terminate(Option_t *option);

  void         SetCaloName(const char *n)        { fCaloName   = n; }
  void         SetTracksName(const char *n)      { fTracksName = n; }

 protected:
  TString      fTracksName;         // name of track collection (if "" use branch)
  TString      fCaloName;           // name of calo collection

 private:
  AliEmcalClusTrackMatcherTask(const AliEmcalClusTrackMatcherTask&);            // not implemented
  AliEmcalClusTrackMatcherTask &operator=(const AliEmcalClusTrackMatcherTask&); // not implemented

  ClassDef(AliEmcalClusTrackMatcherTask, 1) // Cluster-Track matching task
};
#endif
