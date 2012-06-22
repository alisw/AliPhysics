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

  void                        SetTrgClusName(const char *n)                        { fTrgClusName   = n          ; }
  void                        SetCellEnergyCut(Float_t cut)                        { fCellEnergyCut = cut        ; }
  void                        SetDoTrigger(Bool_t trg = kTRUE)                     { fDoTrigger     = trg        ; }

 protected:

  Bool_t                      FillHistograms()                                          ;
  Bool_t                      RetrieveEventObjects()                                    ;
  void                        DoCellLoop(Float_t &sum, Float_t &sum_cut)                ;
  void                        DoTriggerPrimitives(Int_t &maxL1amp, Int_t &maxL1thr)     ;
  Float_t                     DoTriggerClusLoop()                                       ;
  Float_t                     DoTrackLoop()                                             ;
  Float_t                     DoClusterLoop()                                           ;
  void                        DoJetLoop()                                               ;


  Float_t                     fCellEnergyCut;          // Energy cell cut
  Bool_t                      fDoTrigger;              // Make trigger qa plots
  TString                     fTrgClusName;            // Name of trg clus name
  TClonesArray               *fTrgClusters;            //!Trg Clusters

  TH1F                       *fHistCentrality;         //!Event centrality distribution
  TH2F                       *fHistTracksCent;         //!Number of tracks vs. centrality
  TH2F                       *fHistClusCent;           //!Number of clusters vs. centrality

  TH2F                       *fHistMaxL1FastORCent;    //!Maximum L1 trigger FastOR amplitude vs. centrality
  TH2F                       *fHistMaxL1ClusCent;      //!Maximum L1 trigger cluster amplitude vs. centrality
  TH2F                       *fHistMaxL1ThrCent;       //!Maximum L1 trigger threshold vs. centrality
 
  TH1F                       *fHistTracksPt;             //!Pt spectrum of tracks
  TH2F                       *fHistTrPhiEta;             //!Phi-Eta distribution of tracks
  TH2F                       *fHistTrEmcPhiEta;          //!Phi-Eta emcal distribution of tracks

  TH3F                       *fHistClusPhiEtaEnergy;     //!Phi-Eta-Energy distribution of clusters
  TH2F                       *fHistNCellsEnergy;         //!Number of cells vs. energy of cluster

  TH3F                       *fHistJetsPhiEtaPt[4];      //!Phi-Eta-Pt distribution of jets
  TH1F                       *fHistJetsPtNonBias[4];     //!Non biased inclusive jet pt spectrum
  TH1F                       *fHistJetsPtClus[4];        //!Inclusive jet pt spectrum cluster biased
  TH1F                       *fHistJetsPtTrack[4];       //!Inclusive jet pt spectrum track biased
  TH1F                       *fHistJetsPt[4];            //!Biased inclusive jet pt spectrum
  TH2F                       *fHistJetsPtAreaNonBias[4]; //!Non biased pt vs. area of jets
  TH2F                       *fHistJetsPtArea[4];        //!Biased pt vs. area of jets

  TH1F                       *fHistCellsEnergy;        //!Energy spectrum of cells

  TH2F                       *fHistChVSneCells;        //!Charged vs. neutral (cells) energy
  TH2F                       *fHistChVSneClus;         //!Charged vs. neutral (clusters) energy
  TH2F                       *fHistChVSneCorrCells;    //!Charged vs. neutral (corrected cells) energy

  TH1F                       *fHistTrackPhi[5];        //!Phi distribution of hybrid tracks
  TH1F                       *fHistTrackEta[5];        //!Eta distribution of hybrid tracks

 private:
  AliAnalysisTaskSAQA(const AliAnalysisTaskSAQA&);            // not implemented
  AliAnalysisTaskSAQA &operator=(const AliAnalysisTaskSAQA&); // not implemented

  ClassDef(AliAnalysisTaskSAQA, 8) // Quality task for Emcal analysis
};
#endif
