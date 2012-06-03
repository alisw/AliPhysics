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

  void         SetDoClusTrackMatching(Bool_t b) { fDoClusTrack = b; }
  void         SetDoTrackClusMatching(Bool_t b) { fDoTrackClus = b; }
  void         SetMaxDistance(Double_t d)       { fMaxDistance = d; }

 protected:
  void         DoMatching(TClonesArray *array1, TClonesArray *array2);

  Bool_t       fDoClusTrack;        // match clusters to tracks (one -> many)
  Bool_t       fDoTrackClus;        // match tracks to clusters (one -> many) 
  Double_t     fMaxDistance;        // maximum distance to match clusters and tracks

 private:
  AliEmcalClusTrackMatcherTask(const AliEmcalClusTrackMatcherTask&);            // not implemented
  AliEmcalClusTrackMatcherTask &operator=(const AliEmcalClusTrackMatcherTask&); // not implemented

  ClassDef(AliEmcalClusTrackMatcherTask, 2) // Cluster-Track matching task
};
#endif
