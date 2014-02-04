#ifndef ALIEMCALCLUSTRACKMATCHERTASK_H
#define ALIEMCALCLUSTRACKMATCHERTASK_H

// $Id$

#include "AliAnalysisTaskEmcal.h"

class AliEmcalClusTrackMatcherTask : public AliAnalysisTaskEmcal {
 public:
  AliEmcalClusTrackMatcherTask();
  AliEmcalClusTrackMatcherTask(const char *name, Bool_t histo=kFALSE);
  virtual ~AliEmcalClusTrackMatcherTask();

  void         SetMaxDistance(Double_t d)       { fMaxDistance = d; }
  void         UserCreateOutputObjects();

 protected:
  void         ExecOnce()                 ;
  Int_t        GetMomBin(Double_t p) const;
  Bool_t       Run();

  Double_t     fMaxDistance            ;// maximum distance to match clusters and tracks
  TH1         *fHistMatchEtaAll        ;//!deta distribution
  TH1         *fHistMatchPhiAll        ;//!dphi distribution
  TH1         *fHistMatchEta[8][9][2]  ;//!deta distribution
  TH1         *fHistMatchPhi[8][9][2]  ;//!dphi distribution

 private:
  AliEmcalClusTrackMatcherTask(const AliEmcalClusTrackMatcherTask&);            // not implemented
  AliEmcalClusTrackMatcherTask &operator=(const AliEmcalClusTrackMatcherTask&); // not implemented

  ClassDef(AliEmcalClusTrackMatcherTask, 5) // Cluster-Track matching task
};
#endif
