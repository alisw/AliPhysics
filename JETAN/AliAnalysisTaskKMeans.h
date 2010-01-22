#ifndef AliAnalysisTaskKMeans_cxx
#define AliAnalysisTaskKMeans_cxx
 /* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     Analysis Task that uses the Soft K-Means Algorithm to find clusters in
//     the eta-phi space of Minimum Bias. No pt information is used for the clustering.
//     
//
//     Author: Andreas Morsch (CERN)
//     andreas.morsch@cern.ch
//-------------------------------------------------------------------------

class TH1F;
class TH2F;
class TList;
class TProfile;

class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;


#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskKMeans : public AliAnalysisTaskSE {
 public:
    AliAnalysisTaskKMeans();
  AliAnalysisTaskKMeans(const char *name);
  virtual ~AliAnalysisTaskKMeans() {}
  virtual void     UserCreateOutputObjects();
  virtual void     UserExec(Option_t *option);
  virtual void     Terminate(Option_t *);
  virtual void     SetCuts(AliESDtrackCuts* cuts) {fCuts = cuts;}
  virtual Double_t DeltaPhi(Double_t phi1, Double_t phi2);
  virtual Double_t DeltaR(Double_t phi1, Double_t eta1, Double_t phi2, Double_t eta2);
  
 private:
  // Others
  TList*           fHists;         // Histograms
  TH1F*            fH1CEta;        // Eta distribution of clusters
  TH1F*            fH1CPhi;        // Phi distribution of clusters  
  TH1F*            fH1CEtaR;       // Eta distribution of clusters for rndm evnt
  TH2F*            fH2N1N2;        // Cluster sizes 
  TH1F*            fH1Pt;          // pt outside clusters
  TH1F*            fH1PtC;         // pt outside clusters
  TH1F*            fH1PtC1;        // pt dr > 0.4
  TH1F*            fH1PtC2;        // pt dr > 0.2 
  TH1F*            fH1SPt;         // sum pt
  TH1F*            fH1SPtC;        // sum pt
  TH1F*            fH1DPhi;        // Dphi wr to cluster
  TH1F*            fH1DR;          // DR   wr to cluster
  TH1F*            fH1DRR;         // DR   wr to cluster from rndm events   
  TH2F*            fH2DPhiEta;     // eta-phi wr to cluster
  TH2F*            fH2DPhiEtaR;    // eta-phi wr to cluster for rndm events 
  TH2F*            fH2DPhiEtaL;    // eta-phi of leading particle
  TH2F*            fH2DPhiEtaLR;   // eta-phi of leading particle
  TH2F*            fH2DPhiEtaC;    // eta-phi of Clusters
  TH2F*            fH2DPhiEtaCR;   // eta-phi of Clusters
  
  AliESDtrackCuts* fCuts;             // List of cuts
  ClassDef(AliAnalysisTaskKMeans, 1); // A k-means clustering analysis
};

#endif
