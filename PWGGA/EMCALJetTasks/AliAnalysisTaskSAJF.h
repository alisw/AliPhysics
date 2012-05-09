#ifndef ALIANALYSISTASKSAJF_H
#define ALIANALYSISTASKSAJF_H

// $Id$

class TClonesArray;
class TString;
class AliVTrack;
class AliVCluster;
class AliEmcalJet;
class TList;
class TH1F;
class TH2F;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSAJF : public AliAnalysisTaskSE {
 public:
  
  enum SAJFAnaType {
    kFullAcceptance  = 0,     // Full acceptance
    kEMCAL           = 1,     // EMCal acceptance only
    kEMCALFiducial   = 2      // EMCal fiduacial region only
  };

  AliAnalysisTaskSAJF();
  AliAnalysisTaskSAJF(const char *name);
  virtual ~AliAnalysisTaskSAJF();

  void                        UserCreateOutputObjects();
  void                        UserExec(Option_t *option);
  void                        Terminate(Option_t *option);
  void                        Init();

  void                        SetTracksName(const char *n)                         { fTracksName    = n          ; }
  void                        SetClusName(const char *n)                           { fCaloName      = n          ; }
  void                        SetJetsName(const char *n)                           { fJetsName      = n          ; }
  void                        SetKtJetsName(const char *n)                         { fKtJetsName    = n          ; }
  void                        SetEmbJetsName(const char *n)                        { fEmbJetsName   = n          ; }
  void                        SetTrgClusName(const char *n)                        { fTrgClusName   = n          ; }
  void                        SetAnaType(SAJFAnaType type)                         { fAnaType       = type       ; }
  void                        SetJetRadius(Float_t r)                              { fJetRadius     = r          ; } 
  void                        SetHistoBins(Int_t nbins, Float_t min, Float_t max)  { fNbins = nbins; fMinPt = min; fMaxPt = max; }

 protected:

  AliVTrack                  *GetTrack(const Int_t i)                                              const;
  Int_t                       GetNumberOfTracks()                                                  const;
  AliVCluster                *GetCaloCluster(const Int_t i)                                        const;
  Int_t                       GetNumberOfCaloClusters()                                            const;
  AliEmcalJet                *GetJet(const Int_t i)                                                const;
  Int_t                       GetNumberOfJets()                                                    const;
  AliEmcalJet                *GetKtJet(const Int_t i)                                              const;
  Int_t                       GetNumberOfKtJets()                                                  const;
  AliEmcalJet                *GetEmbJet(const Int_t i)                                             const;
  Int_t                       GetNumberOfEmbJets()                                                 const;
  AliVCluster                *GetTrgCluster(const Int_t i)                                         const;
  Int_t                       GetNumberOfTrgClusters()                                             const;
  Bool_t                      AcceptTrack(AliVTrack* track, Bool_t acceptMC = kFALSE)              const;
  Bool_t                      AcceptCluster(AliVCluster* clus, Bool_t acceptMC = kFALSE)           const;
  Bool_t                      AcceptJet(AliEmcalJet* jet)                                          const;
  Bool_t                      AcceptJet(Float_t eta, Float_t phi)                                  const;
  Bool_t                      IsJetTrack(AliEmcalJet* jet, Int_t itrack, Bool_t sorted = kTRUE)    const;
  Bool_t                      IsJetCluster(AliEmcalJet* jet, Int_t iclus, Bool_t sorted = kTRUE)   const;
  Float_t                     GetArea()                                                            const;
  void                        RetrieveEventObjects()                                                    ;
  void                        FillHistograms()                                                          ;
  void                        DoJetLoop(Int_t &maxJetIndex, Int_t &max2JetIndex)                        ;
  Float_t                     DoKtJetLoop(Int_t nLJs = 2)                                                             ;
  Bool_t                      DoEmbJetLoop(Float_t &maxJetPt, Float_t &maxPartPt)                       ;
  Float_t                     DoTrackLoop(Int_t maxJetIndex, Int_t max2JetIndex)                        ;
  Float_t                     DoClusterLoop(Int_t maxJetIndex, Int_t max2JetIndex)                      ;
  Float_t                     GetRigidConePt(AliEmcalJet *jet = 0,  Float_t minD = 0.8)                 ;


  SAJFAnaType                 fAnaType;                    // Analysis type
  Float_t                     fMinEta;                     // Minimum eta accepatance
  Float_t                     fMaxEta;                     // Maximum eta accepatance
  Float_t                     fMinPhi;                     // Minimum phi accepatance
  Float_t                     fMaxPhi;                     // Maximum phi accepatance  
  Float_t                     fJetRadius;                  // jet radius
  TList                      *fOutput;                     // Output list
  TString                     fTracksName;                 // Name of track collection
  TString                     fCaloName;                   // Name of calo cluster collection
  TString                     fJetsName;                   // Name of jet collection
  TString                     fKtJetsName;                 // Name of kt jet collection
  TString                     fEmbJetsName;                // Name of embedded jets collection
  TString                     fTrgClusName;                // Name of trg clus name
  TClonesArray               *fTracks;                     //!Tracks
  TClonesArray               *fCaloClusters;               //!Clusters
  TClonesArray               *fJets;                       //!Jets
  TClonesArray               *fKtJets;                     //!Kt Jets
  TClonesArray               *fEmbJets;                    //!Embedded Jets
  TClonesArray               *fTrgClusters;                //!Trg Clusters
  AliCentrality              *fCent;                       // Event centrality
  Int_t                       fCentBin;                    // Event centrality bin
  TH1F                       *fHistCentrality;             // Event centrality distribution
  TH2F                       *fHistJetPhiEta;              // Phi-Eta distribution of jets
  TH2F                       *fHistRhoPartVSleadJetPt;     // Background et density of particles (clusters+tracks) vs. leading jet pt
  TH2F                       *fHistMedKtVSRhoPart;         // Median of the pt density of kt jets vs. background pt density of particles
  TH2F                       *fHistRCPtVSRhoPart;          // Random cone Pt vs. background pt density of particles
  TH1F                       *fHistJetsPt[4];              // Jet pt spectrum
  TH1F                       *fHistJetsNEF[4];             // Jet neutral energy fraction
  TH1F                       *fHistJetsZ[4];               // Constituent Pt over Jet Pt ratio
  TH1F                       *fHistLeadingJetPt[4];        // Leading jet pt spectrum
  TH1F                       *fHist2LeadingJetPt[4];       // Second leading jet pt spectrum
  TH1F                       *fHistTracksPtLJ[4];          // Pt spectrum of tracks
  TH1F                       *fHistClusEtLJ[4];            // Et spectrum of clusters
  TH1F                       *fHistTracksPtBkg[4];         // Pt spectrum of tracks
  TH1F                       *fHistClusEtBkg[4];           // Et spectrum of clusters
  TH1F                       *fHistBkgClusPhiCorr[4];      // Background clusters phi correlations
  TH1F                       *fHistBkgTracksPhiCorr[4];    // Background tracks phi correlations
  TH1F                       *fHistBkgClusMeanRho[4];      // Background clusters mean et density
  TH1F                       *fHistBkgTracksMeanRho[4];    // Background tracks mean pt density
  TH1F                       *fHistBkgLJetPhiCorr[4];      // Background particles phi correlation with leading jet
  TH1F                       *fHistMedianPtKtJet[4];       // Median of the pt density of kt jets, excluded the 2 most energetic
  TH1F                       *fHistDeltaPtRC[4];           // deltaPt = Pt(RC) - A * rhoKt
  TH1F                       *fHistDeltaPtRCExLJ[4];       // deltaPt = Pt(RC) - A * rhoKt, imposing min distance from leading jet
  TH1F                       *fHistRCPt[4];                // Random cone pt
  TH1F                       *fHistRCPtExLJ[4];            // Random cone pt, imposing min distance from leading jet
  TH1F                       *fHistEmbJets[4];             // Pt distribution of embedded jets
  TH1F                       *fHistEmbPart[4];             // Pt distribution of embedded particle
  TH1F                       *fHistDeltaPtMedKtEmb[4];     // deltaPt = Pt(embjet) - A * rhoKt - Pt(embtrack)
  Int_t                       fNbins;                      // No. of pt bins
  Float_t                     fMinPt;                      // Min pt
  Float_t                     fMaxPt;                      // Max pt

 private:
  AliAnalysisTaskSAJF(const AliAnalysisTaskSAJF&);            // not implemented
  AliAnalysisTaskSAJF &operator=(const AliAnalysisTaskSAJF&); // not implemented

  ClassDef(AliAnalysisTaskSAJF, 2) // jet finder task
};
#endif
