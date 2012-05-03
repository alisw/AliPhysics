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

  void                        SetClusName(const char *n)                           { fCaloName      = n          ; }
  void                        SetJetsName(const char *n)                           { fJetsName      = n          ; }
  void                        SetKtJetsName(const char *n)                         { fKtJetsName    = n          ; }
  void                        SetTracksName(const char *n)                         { fTracksName    = n          ; }
  void                        SetTrgClusName(const char *n)                        { fTrgClusName   = n          ; }
  void                        SetAnaType(SAJFAnaType type)                         { fAnaType       = type       ; }
  void                        SetHistoBins(Int_t nbins, Float_t min, Float_t max)  { fNbins = nbins; fMinPt = min; fMaxPt = max; }

 protected:

  AliVTrack                  *GetTrack(const Int_t i)          const;
  Int_t                       GetNumberOfTracks()              const;
  AliVCluster                *GetCaloCluster(const Int_t i)    const;
  Int_t                       GetNumberOfCaloClusters()        const;
  AliEmcalJet                *GetJet(const Int_t i)            const;
  Int_t                       GetNumberOfJets()                const;
  AliEmcalJet                *GetKtJet(const Int_t i)          const;
  Int_t                       GetNumberOfKtJets()              const;
  AliVCluster                *GetTrgCluster(const Int_t i)     const;
  Int_t                       GetNumberOfTrgClusters()         const;
  void                        FillHistograms()                      ;
  void                        RetrieveEventObjects()                ;
  Bool_t                      AcceptTrack(AliVTrack* track)    const;
  Bool_t                      AcceptJet(AliEmcalJet* jet)      const;
  Bool_t                      IsJetTrack(AliEmcalJet* jet, Int_t itrack, Bool_t sorted = kTRUE)    const;
  Bool_t                      IsJetCluster(AliEmcalJet* jet, Int_t iclus, Bool_t sorted = kTRUE)   const;
  Float_t                     GetArea()                        const;
  void                        DoJetLoop(Int_t &maxJetIndex, Int_t &max2JetIndex);
  Float_t                     DoKtJetLoop();
  Float_t                     DoTrackLoop(Int_t maxJetIndex, Int_t max2JetIndex);
  Float_t                     DoClusterLoop(Int_t maxJetIndex, Int_t max2JetIndex);


  SAJFAnaType                 fAnaType;                    // analysis type
  TList                      *fOutput;                     // Output list
  TString                     fTracksName;                 // name of track collection
  TString                     fCaloName;                   // name of calo cluster collection
  TString                     fJetsName;                   // name of jet collection
  TString                     fKtJetsName;                 // name of kt jet collection
  TString                     fTrgClusName;                // name of trg clus name
  TClonesArray               *fTracks;                     //!Tracks
  TClonesArray               *fCaloClusters;               //!Clusters
  TClonesArray               *fJets;                       //!Jets
  TClonesArray               *fKtJets;                     //!Kt Jets
  TClonesArray               *fTrgClusters;                //!Trg Clusters
  AliCentrality              *fCent;                       // Event centrality
  Int_t                       fCentBin;                    // Event centrality bin
  TH1F                       *fHistCentrality;             // Event centrality distribution
  TH2F                       *fHistJetPhiEta;              // Phi-Eta distribution of jets
  TH2F                       *fHistRhoPartVSleadJetPt;     // Background et density of particles (clusters+tracks) vs. leading jet pt
  TH2F                       *fHistMedKtVSRhoPart;         // Median of the pt density of kt jets vs. background pt density of particles
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
  Int_t                       fNbins;                      // No. of pt bins
  Float_t                     fMinPt;                      // Min pt
  Float_t                     fMaxPt;                      // Max pt

 private:
  AliAnalysisTaskSAJF(const AliAnalysisTaskSAJF&);            // not implemented
  AliAnalysisTaskSAJF &operator=(const AliAnalysisTaskSAJF&); // not implemented

  ClassDef(AliAnalysisTaskSAJF, 1) // Isolated photons task
};
#endif
