#ifndef ALIANALYSISTASKSAQA_H
#define ALIANALYSISTASKSAQA_H

// $Id$

class TClonesArray;
class TString;
class TH1F;
class TH2F;
class TH3F;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskSAQA : public AliAnalysisTaskEmcalJet {
 public:
  AliAnalysisTaskSAQA();
  AliAnalysisTaskSAQA(const char *name);
  virtual ~AliAnalysisTaskSAQA();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  void                        SetCellEnergyCut(Float_t cut)                        { fCellEnergyCut      = cut        ; }
  void                        SetParticleLevel(Bool_t s)                           { fParticleLevel      = s          ; }

 protected:

  Bool_t                      FillHistograms()                                              ;
  Bool_t                      RetrieveEventObjects()                                        ;
  Int_t                       DoCellLoop(Float_t &sum, Float_t &sum_cut)                    ;
  Float_t                     DoTrackLoop()                                                 ;
  Float_t                     DoClusterLoop()                                               ;
  void                        DoJetLoop()                                                   ;
  Double_t                    GetFcross(AliVCluster *cluster, AliVCaloCells *cells)         ;

  Float_t                     fCellEnergyCut;            // Energy cell cut
  Bool_t                      fParticleLevel;            // Set particle level analysis
  Int_t                       fNclusters;                //!Number of accepted clusters in the event
  Int_t                       fNtracks;                  //!Number of accepted tracks in the event
  Int_t                       fNjets;                    //!Number of accepted jets in the event
 
  // General histograms
  TH2F                       *fHistTracksCent;           //!Number of tracks vs. centrality
  TH2F                       *fHistClusCent;             //!Number of clusters vs. centrality
  TH2F                       *fHistJetsCent;             //!Number of jets vs. centrality
  TH2F                       *fHistClusTracks;           //!Number of clusters vs. number of tracks
  TH2F                       *fHistJetsParts;            //!Number of jets vs. number of particles (tracks+clusters)
  TH2F                       *fHistCellsCent;            //!Number of cells vs. centrality
  TH2F                       *fHistCellsTracks;          //!Number of cells vs. number of tracks

  // Tracks
  TH3F                       *fHistTrPhiEtaPt[4][4];     //!Phi-Eta-Pt distribution of tracks
  TH3F                       *fHistTrPhiEtaPtNegLab[4];  //!Phi-Eta-Pt distribution of tracks with negative labels
  TH2F                       *fHistTrEmcPhiEta;          //!Phi-Eta emcal propagated distribution of tracks
  TH2F                       *fHistTrPhiEtaNonProp;      //!Phi-Eta distribution of non emcal propagated tracks
  TH2F                       *fHistDeltaEtaPt;           //!Eta-EtaProp vs. Pt
  TH2F                       *fHistDeltaPhiPt;           //!Phi-PhiProp vs. Pt

  // Clusters
  TH3F                       *fHistClusPhiEtaEnergy[4];  //!Phi-Eta-Energy distribution of clusters
  TH2F                       *fHistNCellsEnergy;         //!Number of cells vs. energy of cluster
  TH2F                       *fHistFcrossEnergy;         //!Fcross vs. energy of cluster
  TH2F                       *fHistClusTimeEnergy;       //!Time vs. energy of cluster

  //Jets
  TH3F                       *fHistJetsPhiEtaPt[4];      //!Phi-Eta distribution of jets
  TH2F                       *fHistJetsPtArea[4];        //!Pt vs. area of jets

  // EMCAL Cells
  TH2F                       *fHistCellsAbsIdEnergy;    //!Energy spectrum of cells

  // Had corr QA
  TH2F                       *fHistChVSneCells;          //!Charged vs. neutral (cells) energy
  TH2F                       *fHistChVSneClus;           //!Charged vs. neutral (clusters) energy
  TH2F                       *fHistChVSneCorrCells;      //!Charged vs. neutral (corrected cells) energy

 private:
  AliAnalysisTaskSAQA(const AliAnalysisTaskSAQA&);            // not implemented
  AliAnalysisTaskSAQA &operator=(const AliAnalysisTaskSAQA&); // not implemented

  ClassDef(AliAnalysisTaskSAQA, 17) // Quality task for Emcal analysis
};
#endif
