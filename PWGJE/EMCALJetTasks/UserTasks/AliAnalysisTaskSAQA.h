#ifndef ALIANALYSISTASKSAQA_H
#define ALIANALYSISTASKSAQA_H

// $Id$

class TH1;
class TH2;
class TH3;
class THnSparse;
class AliVVZERO;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskSAQA : public AliAnalysisTaskEmcalJet {
 public:
  AliAnalysisTaskSAQA();
  AliAnalysisTaskSAQA(const char *name);
  virtual ~AliAnalysisTaskSAQA();

  void                        UserCreateOutputObjects();

  void                        SetCellEnergyCut(Float_t cut)                        { fCellEnergyCut      = cut        ; }
  void                        SetParticleLevel(Bool_t s)                           { fParticleLevel      = s          ; }
  void                        SetMC(Bool_t m)                                      { fIsMC               = m          ; }
  void                        SetAdditionalCentEst(const char* meth2, const char* meth3="") { fCentMethod2 = meth2; fCentMethod3 = meth3; }
  void                        SetDoV0QA(Int_t b)                                   { fDoV0QA             = b          ; }
  void                        SetDoEPQA(Int_t b)                                   { fDoEPQA             = b          ; }
  void                        SetMaxCellsInCluster(Int_t b)                        { fMaxCellsInCluster  = b          ; }
  void                        SetDoLeadingObjectPosition(Int_t b)                  { fDoLeadingObjectPosition = b     ; }

 protected:

  void                        ExecOnce()                                                    ;
  Bool_t                      FillHistograms()                                              ;
  void                        FillEventQAHisto(Float_t cent, Float_t cent2, Float_t cent3, Float_t v0a, Float_t v0c, Float_t ep, Float_t rho, 
					       Int_t ntracks, Int_t nclusters, Int_t ncells, Int_t njets, 
					       Float_t maxTrackPt, Float_t maxTrackEta, Float_t maxTrackPhi,
					       Float_t maxClusterE, Float_t maxClusterEta, Float_t maxClusterPhi,
					       Float_t maxJetPt, Float_t maxJetEta, Float_t maxJetPhi);
  Bool_t                      RetrieveEventObjects()                                        ;
  Int_t                       DoCellLoop(Float_t &sum)                    ;
  Int_t                       DoTrackLoop(Float_t &sum, AliVParticle* &leading)             ;
  Int_t                       DoClusterLoop(Float_t &sum, AliVCluster* &leading)            ;
  Int_t                       DoJetLoop(AliEmcalJet* &leading)                              ;
  Double_t                    GetFcross(AliVCluster *cluster, AliVCaloCells *cells)         ;

  Float_t                     fCellEnergyCut;            // Energy cell cut
  Bool_t                      fParticleLevel;            // Set particle level analysis
  Bool_t                      fIsMC;                     // Trigger, MC analysis
  TString                     fCentMethod2;              // Centrality method 2
  TString                     fCentMethod3;              // Centrality method 3
  Int_t                       fDoV0QA;                   // Add V0 QA histograms
  Int_t                       fDoEPQA;                   // Add event plane QA histograms
  Int_t                       fDoLeadingObjectPosition;  // Add axis for leading object position (eta-phi)
  Int_t                       fMaxCellsInCluster;        // Maximum number (approx) of cells in a cluster
  Double_t                    fCent2;                    //!Event centrality with method 2
  Double_t                    fCent3;                    //!Event centrality with method 3
  AliVVZERO                  *fVZERO;                    //!Event V0 object
  Double_t                    fV0ATotMult;               //!Event V0A total multiplicity
  Double_t                    fV0CTotMult;               //!Event V0C total multiplicity
 
  // General histograms
  THnSparse                  *fHistEventQA;              //!Event-wise QA observables

  // Tracks
  TH1                        *fHistTrNegativeLabels[4];  //!Percentage of negative label tracks
  TH1                        *fHistTrZeroLabels[4];      //!Percentage of tracks with label=0
  TH3                        *fHistTrPhiEtaPt[4][4];     //!Phi-Eta-Pt distribution of tracks
  TH2                        *fHistTrPhiEtaZeroLab[4];   //!Phi-Eta distribution of tracks with label=0
  TH1                        *fHistTrPtZeroLab[4];       //!Pt distribution of tracks with label=0
  TH2                        *fHistTrEmcPhiEta[4];       //!Phi-Eta emcal propagated distribution of tracks
  TH1                        *fHistTrEmcPt[4];           //!Pt emcal propagated distribution of tracks
  TH2                        *fHistTrPhiEtaNonProp[4];   //!Phi-Eta distribution of non emcal propagated tracks
  TH1                        *fHistTrPtNonProp[4];       //!Pt distribution of non emcal propagated tracks
  TH2                        *fHistDeltaEtaPt[4];        //!Eta-EtaProp vs. Pt
  TH2                        *fHistDeltaPhiPt[4];        //!Phi-PhiProp vs. Pt
  TH2                        *fHistDeltaPtvsPt[4];       //!Pt-PtProp vs. Pt

  // Clusters
  TH3                        *fHistClusPhiEtaEnergy[4];       //!Phi-Eta-Energy distribution of clusters
  TH2                        *fHistClusDeltaPhiEPEnergy[4];   //!DeltaPhi EP vs Energy of clusters
  TH2                        *fHistNCellsEnergy[4];           //!Number of cells vs. energy of cluster
  TH2                        *fHistFcrossEnergy[4];           //!Fcross vs. energy of cluster
  TH2                        *fHistClusTimeEnergy[4];         //!Time vs. energy of cluster
  TH1                        *fHistClusMCEnergyFraction[4];   //!MC energy fraction (embedding)

  // EMCAL Cells
  TH2                        *fHistCellsAbsIdEnergy[4];  //!Energy spectrum of cells

  // Jets
  TH2                        *fHistJetsPhiEta[4];        //!Phi-Eta distribution of jets
  TH2                        *fHistJetsPtArea[4];        //!Pt vs. area of jets

 private:
  AliAnalysisTaskSAQA(const AliAnalysisTaskSAQA&);            // not implemented
  AliAnalysisTaskSAQA &operator=(const AliAnalysisTaskSAQA&); // not implemented

  ClassDef(AliAnalysisTaskSAQA, 23) // Quality task for Emcal analysis
};
#endif
