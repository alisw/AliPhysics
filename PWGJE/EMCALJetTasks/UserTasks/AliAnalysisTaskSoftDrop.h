#ifndef ALIANALYSISTASKSOFTDROP_H
#define ALIANALYSISTASKSOFTDROP_H

// $Id$

class TH1;
class TH2;
class TH3;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskSoftDrop : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskSoftDrop();
  AliAnalysisTaskSoftDrop(const char *name);
  virtual ~AliAnalysisTaskSoftDrop();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

 protected:
  void                        ExecOnce();
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;
  void                        CheckClusTrackMatching();

  // General histograms
  TH1                       **fHistTracksPt;            //!Track pt spectrum
  TH1                       **fHistClustersPt;          //!Cluster pt spectrum
  TH1                       **fHistLeadingJetPt;        //!Leading jet pt spectrum
  TH2                       **fHistJetsPhiEta;          //!Phi-Eta distribution of jets
  TH2                       **fHistJetsPtArea;          //!Jet pt vs. area
  TH2                       **fHistJetsPtLeadHad;       //!Jet pt vs. leading hadron
  TH2                       **fHistJetsCorrPtArea;      //!Jet pt - bkg vs. area
  TH3                        *fHistPtDEtaDPhiTrackClus; //!track pt, delta eta, delta phi to matched cluster
  TH3                        *fHistPtDEtaDPhiClusTrack; //!cluster pt, delta eta, delta phi to matched track
  TH1                       **fHistNTracks;             //! number of tracks per event

  TH1                        *fHistClustDx;             //!
  TH1                        *fHistClustDz;             //!
  TH1                        *fNAccJets;                //! number of jets per event
  TH1                        *fhZg;                     //!<! distribution of zg
  TH2                        *fhCorrPtZg;               //!<! distribution of zg, jet pt-diff

  AliJetContainer            *fJetsCont;                   //!Jets
  AliParticleContainer       *fTracksCont;                 //!Tracks
  AliClusterContainer        *fCaloClustersCont;           //!Clusters  

 private:
  AliAnalysisTaskSoftDrop(const AliAnalysisTaskSoftDrop&);            // not implemented
  AliAnalysisTaskSoftDrop &operator=(const AliAnalysisTaskSoftDrop&); // not implemented

  ClassDef(AliAnalysisTaskSoftDrop, 1) // jet sample analysis task
};
#endif
