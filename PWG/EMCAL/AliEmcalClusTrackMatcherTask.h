#ifndef ALIEMCALCLUSTRACKMATCHERTASK_H
#define ALIEMCALCLUSTRACKMATCHERTASK_H

// $Id$

#include "AliAnalysisTaskEmcal.h"

class AliEmcalClusTrackMatcherTask : public AliAnalysisTaskEmcal {
 public:
  AliEmcalClusTrackMatcherTask();
  AliEmcalClusTrackMatcherTask(const char *name);
  virtual ~AliEmcalClusTrackMatcherTask();

  Bool_t       Run();
  void         SetMaxDistance(Double_t d)       { fMaxDistance = d; }

 protected:
  Double_t     fMaxDistance;        // maximum distance to match clusters and tracks

 private:
  AliEmcalClusTrackMatcherTask(const AliEmcalClusTrackMatcherTask&);            // not implemented
  AliEmcalClusTrackMatcherTask &operator=(const AliEmcalClusTrackMatcherTask&); // not implemented

  ClassDef(AliEmcalClusTrackMatcherTask, 3) // Cluster-Track matching task
};
#endif
