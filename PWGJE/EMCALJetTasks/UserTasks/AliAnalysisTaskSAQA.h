#ifndef ALIANALYSISTASKSAQA_H
#define ALIANALYSISTASKSAQA_H

// $Id$

class TH1;
class TH2;
class TH3;

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
  void                        SetMC(Bool_t m)                                      { fIsMC               = m          ; }

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
  Bool_t                      fIsMC;                     // Trigger, MC analysis
  Int_t                       fNclusters;                //!Number of accepted clusters in the event
  Int_t                       fNtracks;                  //!Number of accepted tracks in the event
  Int_t                       fNjets;                    //!Number of accepted jets in the event
 
  // General histograms
  TH2                        *fHistTracksCent;           //!Number of tracks vs. centrality
  TH2                        *fHistClusCent;             //!Number of clusters vs. centrality
  TH2                        *fHistJetsCent;             //!Number of jets vs. centrality
  TH2                        *fHistClusTracks;           //!Number of clusters vs. number of tracks
  TH2                        *fHistJetsParts;            //!Number of jets vs. number of particles (tracks+clusters)
  TH2                        *fHistCellsCent;            //!Number of cells vs. centrality
  TH2                        *fHistCellsTracks;          //!Number of cells vs. number of tracks

  // Tracks
  TH1                        *fHistTrNegativeLabels;     //!Percentage of negative label tracks
  TH1                        *fHistTrZeroLabels;         //!Percentage of tracks with label=0
  TH3                        *fHistTrPhiEtaPt[4][4];     //!Phi-Eta-Pt distribution of tracks
  TH3                        *fHistTrPhiEtaPtZeroLab[4]; //!Phi-Eta-Pt distribution of tracks with label=0
  TH2                        *fHistTrEmcPhiEta[4];       //!Phi-Eta emcal propagated distribution of tracks
  TH1                        *fHistTrEmcPt[4];           //!Phi-Eta emcal propagated distribution of tracks
  TH2                        *fHistTrPhiEtaNonProp[4];   //!Phi-Eta distribution of non emcal propagated tracks
  TH2                        *fHistDeltaEtaPt[4];        //!Eta-EtaProp vs. Pt
  TH2                        *fHistDeltaPhiPt[4];        //!Phi-PhiProp vs. Pt
  TH2                        *fHistDeltaPtvsPt[4];       //!Pt-PtProp vs. Pt

  // Clusters
  TH3                        *fHistClusPhiEtaEnergy[4];  //!Phi-Eta-Energy distribution of clusters
  TH2                        *fHistNCellsEnergy;         //!Number of cells vs. energy of cluster
  TH2                        *fHistFcrossEnergy;         //!Fcross vs. energy of cluster
  TH2                        *fHistClusTimeEnergy;       //!Time vs. energy of cluster

  // Jets
  TH3                        *fHistJetsPhiEtaPt[4];      //!Phi-Eta distribution of jets
  TH2                        *fHistJetsPtArea[4];        //!Pt vs. area of jets

  // EMCAL Cells
  TH2                        *fHistCellsAbsIdEnergy;    //!Energy spectrum of cells

  // Had corr QA
  TH2                        *fHistChVSneCells;          //!Charged vs. neutral (cells) energy
  TH2                        *fHistChVSneClus;           //!Charged vs. neutral (clusters) energy
  TH2                        *fHistChVSneCorrCells;      //!Charged vs. neutral (corrected cells) energy

 private:
  AliAnalysisTaskSAQA(const AliAnalysisTaskSAQA&);            // not implemented
  AliAnalysisTaskSAQA &operator=(const AliAnalysisTaskSAQA&); // not implemented

  ClassDef(AliAnalysisTaskSAQA, 18) // Quality task for Emcal analysis
};
#endif
