#ifndef ALIEMCALCLUSTRACKMATCHERTASK_H
#define ALIEMCALCLUSTRACKMATCHERTASK_H

// $Id$

#include "AliAnalysisTaskEmcal.h"

class AliEmcalClusTrackMatcherTask : public AliAnalysisTaskEmcal {
 public:
  AliEmcalClusTrackMatcherTask();
  AliEmcalClusTrackMatcherTask(const char *name, Bool_t histo=kFALSE);
  virtual ~AliEmcalClusTrackMatcherTask();

  void          SetMaxDistance(Double_t d)       { fMaxDistance = d; }
  void          SetModifyObjs(Bool_t b)          { fModifyObjs  = b; }

 protected:
  void          ExecOnce();
  Int_t         GetMomBin(Double_t p) const;
  Bool_t        Run();
  void          UserCreateOutputObjects();

  Double_t      fMaxDistance;           // maximum distance to match clusters and tracks
  Bool_t        fModifyObjs;            // if true then modify original tracks/clusters
  TClonesArray *fOrigTracks;            //!ptr to original tracks (used if fModifyObjs true)
  TClonesArray *fOrigClus;              //!ptr to original clusters (used if fModifyObjs true)
  TH1          *fHistMatchEtaAll;       //!deta distribution
  TH1          *fHistMatchPhiAll;       //!dphi distribution
  TH1          *fHistMatchEta[8][9][2]; //!deta distribution
  TH1          *fHistMatchPhi[8][9][2]; //!dphi distribution
  
 private:
  AliEmcalClusTrackMatcherTask(const AliEmcalClusTrackMatcherTask&);            // not implemented
  AliEmcalClusTrackMatcherTask &operator=(const AliEmcalClusTrackMatcherTask&); // not implemented

  ClassDef(AliEmcalClusTrackMatcherTask, 6) // Cluster-Track matching task
};
#endif
