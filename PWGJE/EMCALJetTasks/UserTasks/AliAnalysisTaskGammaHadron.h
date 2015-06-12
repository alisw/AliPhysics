#ifndef ALIANALYSISTASKGAMMAHADRON_H
#define ALIANALYSISTASKGAMMAHADRON_H

// $Id$

class TH1;
class TH2;
class TH3;
class THnSparse;
class AliVVZERO;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskGammaHadron : public AliAnalysisTaskEmcalJet {
 public:
	AliAnalysisTaskGammaHadron();
	AliAnalysisTaskGammaHadron(const char *name);
  virtual ~AliAnalysisTaskGammaHadron();

  void                        UserCreateOutputObjects();

  void                        SetCellEnergyCut(Float_t cut)                        { fCellEnergyCut      = cut        ; }
  void                        SetParticleLevel(Bool_t s)                           { fParticleLevel      = s          ; }
  void                        SetMC(Bool_t m)                                      { fIsMC               = m          ; }
  void                        SetAdditionalCentEst(const char* meth2, const char* meth3="") { fCentMethod2 = meth2; fCentMethod3 = meth3; }
  void                        SetDoV0QA(Int_t b)                                   { fDoV0QA             = b          ; }
  void                        SetDoEPQA(Int_t b)                                   { fDoEPQA             = b          ; }
  void                        SetMaxCellsInCluster(Int_t b)                        { fMaxCellsInCluster  = b          ; }
  void                        SetDoLeadingObjectPosition(Int_t b)                  { fDoLeadingObjectPosition = b     ; }
  void                        SetAODfilterBits(Int_t b0 = 0, Int_t b1 = 0)         { fAODfilterBits[0]   = b0  ; fAODfilterBits[1] = b1  ; }

 protected:

  void                        AllocateHistogramArrays()                                     ;
  void                        ExecOnce()        											  ;
  Bool_t                      RetrieveEventObjects()                                        ;
  Bool_t                      FillHistograms()                                              ;
  void                        FillEventQAHisto(Float_t cent, Float_t cent2, Float_t cent3, Float_t v0a, Float_t v0c, Float_t ep, Float_t rho,
					       Int_t ntracks, Int_t nclusters, Int_t ncells,
					       Float_t maxTrackPt, Float_t maxTrackEta, Float_t maxTrackPhi,
					       Float_t maxClusterE, Float_t maxClusterEta, Float_t maxClusterPhi);

  Int_t                       DoCellLoop(Float_t &sum)                                      ;
  Double_t                    GetFcross(AliVCluster *cluster, AliVCaloCells *cells)         ;
  Int_t                       DoClusterLoop(Float_t &sum, AliVCluster* &leading)            ;
  Int_t                       DoTrackLoop(Float_t &sum, AliVParticle* &leading)             ;

  Float_t                     fCellEnergyCut;            // Energy cell cut
  Bool_t                      fParticleLevel;            // Set particle level analysis
  Bool_t                      fIsMC;                     // Trigger, MC analysis
  TString                     fCentMethod2;              // Centrality method 2
  TString                     fCentMethod3;              // Centrality method 3
  Int_t                       fDoV0QA;                   // Add V0 QA histograms
  Int_t                       fDoEPQA;                   // Add event plane QA histograms
  Int_t                       fDoLeadingObjectPosition;  // Add axis for leading object position (eta-phi)
  Int_t                       fMaxCellsInCluster;        // Maximum number (approx) of cells in a cluster
  UInt_t                      fAODfilterBits[2];         // AOD track filter bit map
  Double_t                    fCent2;                    //!Event centrality with method 2
  Double_t                    fCent3;                    //!Event centrality with method 3
  AliVVZERO                  *fVZERO;                    //!Event V0 object
  Double_t                    fV0ATotMult;               //!Event V0A total multiplicity
  Double_t                    fV0CTotMult;               //!Event V0C total multiplicity
 
  // General histograms
  THnSparse                  *fHistEventQA;              //!Event-wise QA observables

  // Tracks
  TH1                       **fHistTrNegativeLabels;  //!Percentage of negative label tracks
  TH1                       **fHistTrZeroLabels;      //!Percentage of tracks with label=0
  TH3                      ***fHistTrPhiEtaPt;        //!Phi-Eta-Pt distribution of tracks
  TH2                       **fHistTrPhiEtaZeroLab;   //!Phi-Eta distribution of tracks with label=0
  TH1                       **fHistTrPtZeroLab;       //!Pt distribution of tracks with label=0
  TH2                       **fHistTrEmcPhiEta;       //!Phi-Eta emcal propagated distribution of tracks
  TH1                       **fHistTrEmcPt;           //!Pt emcal propagated distribution of tracks
  TH2                       **fHistTrPhiEtaNonProp;   //!Phi-Eta distribution of non emcal propagated tracks
  TH1                       **fHistTrPtNonProp;       //!Pt distribution of non emcal propagated tracks
  TH2                       **fHistDeltaEtaPt;        //!Eta-EtaProp vs. Pt
  TH2                       **fHistDeltaPhiPt;        //!Phi-PhiProp vs. Pt
  TH2                       **fHistDeltaPtvsPt;       //!Pt-PtProp vs. Pt

  // Clusters
  TH3                       **fHistClusPhiEtaEnergy;       //!Phi-Eta-Energy distribution of clusters
  TH2                       **fHistClusDeltaPhiEPEnergy;   //!DeltaPhi EP vs Energy of clusters
  TH2                       **fHistNCellsEnergy;           //!Number of cells vs. energy of cluster
  TH2                       **fHistFcrossEnergy;           //!Fcross vs. energy of cluster
  TH2                       **fHistClusTimeEnergy;         //!Time vs. energy of cluster
  TH1                       **fHistClusMCEnergyFraction;   //!MC energy fraction (embedding)

  // EMCAL Cells
  TH2                       **fHistCellsAbsIdEnergy;  //!Energy spectrum of cells

  // ELIANE
  TH1  						*fHistNoClus_pt;           //! No of calorimeter Clusters as a function of p_T
  TH1					    *fHistNoClus_ptH;          //! No of calorimeter Clusters as a function of p_T with a hadron in the second hemisphere
  TH1					    *fHistNoClus_ptLeadH;      //! No of calorimeter Clusters as a function of p_T with a leading hadron in the second hemisphere
  TH1					    *fHistNoClus_Leadpt;        //! No of leading calorimeter Clusters as a function of p_T
  TH1					    *fHistNoClus_LeadptH;       //! No of leading calorimeter Clusters as a function of p_T with a hadron in the second hemisphere
  TH1					    *fHistNoClus_LeadptLeadH;   //! No of leading calorimeter Clusters as a function of p_T with a leading hadron in the second hemisphere

  TH1					    *fHistNoClus_xEH;
  TH1					    *fHistNoClus_LeadxEH;
  TH1					    *fHistNoClus_xELeadH;
  TH1					    *fHistNoClus_LeadxELeadH;
  //
  //

 private:
  AliAnalysisTaskGammaHadron(const AliAnalysisTaskGammaHadron&);            // not implemented
  AliAnalysisTaskGammaHadron &operator=(const AliAnalysisTaskGammaHadron&); // not implemented

  ClassDef(AliAnalysisTaskGammaHadron, 1) // Quality task for Emcal analysis
};
#endif
