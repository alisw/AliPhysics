#ifndef ALIEMCALCLUSTRACKMATCHERTASK_H
#define ALIEMCALCLUSTRACKMATCHERTASK_H

#include "AliAnalysisTaskEmcal.h"

class AliEmcalClusTrackMatcherTask : public AliAnalysisTaskEmcal {
 public:
  AliEmcalClusTrackMatcherTask();
  AliEmcalClusTrackMatcherTask(const char *name, Bool_t histo=kFALSE);
  virtual ~AliEmcalClusTrackMatcherTask();

  void          SetPropDist(Double_t d)           { fPropDist              = d; }
  void          SetDoPropagation(Bool_t b)        { fDoPropagation         = b; }
  void          SetAttemptProp(Bool_t b)          { fAttemptProp           = b; }
  void          SetAttemptPropMatch(Bool_t b)     { fAttemptPropMatch      = b; }
  void          SetMaxDistance(Double_t d)        { fMaxDistance           = d; }
  void          SetAttachEmcalParticles(Bool_t b) { fAttachEmcalParticles  = b; }
  void          SetUpdateTracks(Bool_t b)         { fUpdateTracks          = b; }
  void          SetUpdateClusters(Bool_t b)       { fUpdateClusters        = b; }

 protected:
  void          ExecOnce();
  Int_t         GetMomBin(Double_t p) const;
  Bool_t        Run();
  void          UserCreateOutputObjects();

  void          GenerateEmcalParticles();
  void          DoMatching();
  void          UpdateTracks();
  void          UpdateClusters();
  
  Double_t      fPropDist;              // distance to surface (440cm default)
  Bool_t        fDoPropagation;         // if true then propagate all hybrid tracks to EMCal surface
  Bool_t        fAttemptProp;           // if true then attempt to propagate if not done yet
  Bool_t        fAttemptPropMatch;      // if true then attempt to propagate if not done yet but IsEMCAL is true
  Double_t      fMaxDistance;           // maximum distance to match clusters and tracks
  Bool_t        fAttachEmcalParticles;  // attach emcal particles to the event, so that other tasks can use them
  Bool_t        fUpdateTracks;          // update tracks with matching info
  Bool_t        fUpdateClusters;        // update clusters with matching info

  TClonesArray *fEmcalTracks;           //!emcal tracks
  TClonesArray *fEmcalClusters;         //!emcal clusters
  Int_t         fNEmcalTracks;          //!number of emcal tracks
  Int_t         fNEmcalClusters;        //!number of emcal clusters
  TH1          *fHistMatchEtaAll;       //!deta distribution
  TH1          *fHistMatchPhiAll;       //!dphi distribution
  TH1          *fHistMatchEta[10][9][2]; //!deta distribution
  TH1          *fHistMatchPhi[10][9][2]; //!dphi distribution
  
 private:
  AliEmcalClusTrackMatcherTask(const AliEmcalClusTrackMatcherTask&);            // not implemented
  AliEmcalClusTrackMatcherTask &operator=(const AliEmcalClusTrackMatcherTask&); // not implemented

  ClassDef(AliEmcalClusTrackMatcherTask, 8) // Cluster-Track matching task
};
#endif
