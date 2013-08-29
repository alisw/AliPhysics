#ifndef ALIEMCALCLUSTRACKMATCHERTASK_H
#define ALIEMCALCLUSTRACKMATCHERTASK_H

// $Id$

#include "AliAnalysisTaskEmcalDev.h"

class AliEmcalClusTrackMatcherTask : public AliAnalysisTaskEmcalDev {
 public:
  AliEmcalClusTrackMatcherTask();
  AliEmcalClusTrackMatcherTask(const char *name, Bool_t histo=kFALSE);
  virtual ~AliEmcalClusTrackMatcherTask();

  void         UserCreateOutputObjects();

  void         SetMaxDistance(Double_t d)       { fMaxDistance = d; }

 protected:
  void         ExecOnce()                 ;
  Bool_t       Run();
  Int_t        GetMomBin(Double_t p) const;

  Double_t     fMaxDistance            ;// maximum distance to match clusters and tracks

  TH1         *fHistMatchEta[8][9][2]  ;//!deta distribution
  TH1         *fHistMatchPhi[8][9][2]  ;//!dphi distribution

 private:
  AliEmcalClusTrackMatcherTask(const AliEmcalClusTrackMatcherTask&);            // not implemented
  AliEmcalClusTrackMatcherTask &operator=(const AliEmcalClusTrackMatcherTask&); // not implemented

  ClassDef(AliEmcalClusTrackMatcherTask, 5) // Cluster-Track matching task
};
#endif
