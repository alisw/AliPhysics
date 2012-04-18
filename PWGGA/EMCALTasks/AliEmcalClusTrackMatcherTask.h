#ifndef ALIEMCALCLUSTRACKMATCHERTASK_H
#define ALIEMCALCLUSTRACKMATCHERTASK_H

// $Id$

class TClonesArray;

#include "AliAnalysisTaskSE.h"

class AliEmcalClusTrackMatcherTask : public AliAnalysisTaskSE {
 public:
  AliEmcalClusTrackMatcherTask(const char *name=0);
  virtual ~AliEmcalClusTrackMatcherTask();

  void         UserCreateOutputObjects();
  void         UserExec(Option_t *option);
  void         Terminate(Option_t *option);

  void         SetDoClusTrackMatching(Bool_t b) { fDoClusTrack = b; }
  void         SetDoTrackClusMatching(Bool_t b) { fDoTrackClus = b; }
  void         SetClusName(const char *n)       { fCaloName    = n; }
  void         SetTracksName(const char *n)     { fTracksName  = n; }

 protected:
  TString      fTracksName;         // name of track collection
  TString      fCaloName;           // name of calo cluster collection
  Bool_t       fDoClusTrack;        // match clusters to tracks (one -> many)
  Bool_t       fDoTrackClus;        // match tracks to clusters (one -> many) 

 private:
  AliEmcalClusTrackMatcherTask(const AliEmcalClusTrackMatcherTask&);            // not implemented
  AliEmcalClusTrackMatcherTask &operator=(const AliEmcalClusTrackMatcherTask&); // not implemented

  ClassDef(AliEmcalClusTrackMatcherTask, 1) // Cluster-Track matching task
};
#endif
