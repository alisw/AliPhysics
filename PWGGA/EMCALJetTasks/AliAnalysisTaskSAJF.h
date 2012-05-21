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
    kTPC       = 0,     // TPC only analysis
    kEMCAL     = 1,     // EMCal + TPC analysis
  };

  AliAnalysisTaskSAJF();
  AliAnalysisTaskSAJF(const char *name);
  virtual ~AliAnalysisTaskSAJF();

  void                        UserCreateOutputObjects();
  void                        UserExec(Option_t *option);
  void                        Terminate(Option_t *option);
  void                        Init();

  void                        SetTracksName(const char *n)                         { fTracksName     = n          ; }
  void                        SetClusName(const char *n)                           { fCaloName       = n          ; }
  void                        SetJetsName(const char *n)                           { fJetsName       = n          ; }
  void                        SetEmbJetsName(const char *n)                        { fEmbJetsName    = n          ; }
  void                        SetRhoName(const char *n)                            { fRhoName        = n          ; }
  void                        SetAnaType(SAJFAnaType type)                         { fAnaType        = type       ; }
  void                        SetJetRadius(Float_t r)                              { fJetRadius      = r          ; } 
  void                        SetPtCut(Float_t cut)                                { fPtCut          = cut        ; }
  void                        SetPtCutJetPart(Float_t cut)                         { fPtCutJetPart   = cut        ; }
  void                        SetHistoBins(Int_t nbins, Float_t min, Float_t max)  { fNbins = nbins; fMinPt = min; fMaxPt = max; }

 protected:

  AliVTrack                  *GetTrack(const Int_t i)                                              const;
  Int_t                       GetNumberOfTracks()                                                  const;
  AliVCluster                *GetCaloCluster(const Int_t i)                                        const;
  Int_t                       GetNumberOfCaloClusters()                                            const;
  AliEmcalJet                *GetJet(const Int_t i)                                                const;
  Int_t                       GetNumberOfJets()                                                    const;
  AliEmcalJet                *GetEmbJet(const Int_t i)                                             const;
  Int_t                       GetNumberOfEmbJets()                                                 const;
  Bool_t                      AcceptTrack(AliVTrack* track, Bool_t acceptMC = kFALSE)              const;
  Bool_t                      AcceptCluster(AliVCluster* clus, Bool_t acceptMC = kFALSE)                ;
  Bool_t                      AcceptJet(AliEmcalJet* jet, Bool_t checkPt = kFALSE)                 const;
  Bool_t                      AcceptJet(Float_t eta, Float_t phi)                                  const;
  Bool_t                      IsJetTrack(AliEmcalJet* jet, Int_t itrack, Bool_t sorted = kTRUE)    const;
  Bool_t                      IsJetCluster(AliEmcalJet* jet, Int_t iclus, Bool_t sorted = kTRUE)   const;
  Float_t                     GetArea()                                                            const;

  void                        RetrieveEventObjects()                                                                             ;
  void                        FillHistograms()                                                                                   ;
  void                        DoJetLoop(Int_t &maxJetIndex, Int_t &max2JetIndex)                                                 ;
  Bool_t                      DoEmbJetLoop(AliEmcalJet* &maxJet, TObject* &maxPart)                                              ;
  void                        DoTrackLoop(Float_t &rhoTracksLJex, Float_t &rhoTracks, Int_t maxJetIndex, Int_t max2JetIndex)     ;
  void                        DoClusterLoop(Float_t &rhoClusLJex, Float_t &rhoClus,Int_t maxJetIndex, Int_t max2JetIndex)        ;
  Bool_t                      GetRigidConePt(Float_t &pt, Float_t &eta, Float_t &phi, AliEmcalJet *jet = 0,  Float_t minD = 1.0) ;

  SAJFAnaType                 fAnaType;                    // Analysis type
  Float_t                     fMinEta;                     // Minimum eta accepatance
  Float_t                     fMaxEta;                     // Maximum eta accepatance
  Float_t                     fMinPhi;                     // Minimum phi accepatance
  Float_t                     fMaxPhi;                     // Maximum phi accepatance  
  Float_t                     fJetRadius;                  // Jet radius
  TString                     fTracksName;                 // Name of track collection
  TString                     fCaloName;                   // Name of calo cluster collection
  TString                     fJetsName;                   // Name of jet collection
  TString                     fEmbJetsName;                // Name of embedded jets collection
  TString                     fRhoName;                    // Name of rho object
  Int_t                       fNbins;                      // No. of pt bins
  Float_t                     fMinPt;                      // Min pt in histograms
  Float_t                     fMaxPt;                      // Max pt in histograms
  Float_t                     fPtCut;                      // Cut on particle pt
  Float_t                     fPtCutJetPart;               // Select jets with a particle with minimum pt
  TClonesArray               *fTracks;                     //!Tracks
  TClonesArray               *fCaloClusters;               //!Clusters
  TClonesArray               *fJets;                       //!Jets
  TClonesArray               *fEmbJets;                    //!Embedded Jets
  Float_t                     fCent;                       //!Event centrality
  Int_t                       fCentBin;                    //!Event centrality bin
  Float_t                     fRho;                        //!Event rho
  TList                      *fOutput;                     //!Output list
  TH1F                       *fHistCentrality;             //!Event centrality distribution
  TH2F                       *fHistJetPhiEta;              //!Phi-Eta distribution of jets
  TH2F                       *fHistEmbJetPhiEta;           //!Phi-Eta distribution of embedded jets
  TH2F                       *fHistEmbPartPhiEta;          //!Phi-Eta distribution of embedded particles
  TH2F                       *fHistRhoPartVSleadJetPt;     //!Background et density of particles (clusters+tracks) vs. leading jet pt
  TH2F                       *fHistMedKtVSRhoPart;         //!Median of the pt density of kt jets vs. bkg pt density of particles, excluding the 2 most energetic jets
  TH2F                       *fHistRCPtVSRhoPart;          //!Random cone Pt vs. background pt density of particles
  TH2F                       *fHistMedKtVSEmbBkg;          //!Area(embjet) * rhoKt vs. Pt(embjet) - Pt(embtrack)
  TH2F                       *fHistMedKtVSRCPt;            //!Area(RC) * rhoKt vs. Pt(RC)
  TH2F                       *fHistRCPtExLJVSDPhiLJ;       //!Random cone pt, imposing min distance from leading jet, vs deltaPhi leading jet
  TH2F                       *fHistRCPhiEta;               //!Phi-Eta distribution of embedded jets
  TH1F                       *fHistJetsPt[4];              //!Inclusive jet pt spectrum
  TH1F                       *fHistCorrJetsPt[4];          //!Corrected inclusive jet pt spectrum
  TH1F                       *fHistUnfoldedJetsPt[4];      //!Unfolded inclusive jet pt spectrum
  TH1F                       *fHistJetsPtNonBias[4];       //!Inclusive jet pt spectrum selecting jet with a particle minimum fPtCutJetPart
  TH1F                       *fHistCorrJetsPtNonBias[4];   //!Corrected inclusive jet pt spectrum selecting jet with a particle minimum fPtCutJetPart
  TH2F                       *fHistJetsNEFvsPt[4];         //!Jet neutral energy fraction
  TH2F                       *fHistJetsZvsPt[4];           //!Constituent Pt over Jet Pt ratio
  TH1F                       *fHistLeadingJetPt[4];        //!Leading jet pt spectrum
  TH1F                       *fHistCorrLeadingJetPt[4];    //!Leading jet pt spectrum
  TH1F                       *fHist2LeadingJetPt[4];       //!Second leading jet pt spectrum
  TH1F                       *fHistTracksPtLJ[4];          //!Pt spectrum of tracks
  TH1F                       *fHistClusEtLJ[4];            //!Et spectrum of clusters
  TH1F                       *fHistTracksPtBkg[4];         //!Pt spectrum of tracks
  TH1F                       *fHistClusEtBkg[4];           //!Et spectrum of clusters
  TH1F                       *fHistBkgClusMeanRho[4];      //!Background clusters mean et density
  TH1F                       *fHistBkgTracksMeanRho[4];    //!Background tracks mean pt density
  TH1F                       *fHistMedianPtKtJet[4];       //!Median of the pt density of kt jets, excluded the 2 most energetic
  TH1F                       *fHistDeltaPtRC[4];           //!deltaPt = Pt(RC) - A * rhoKt
  TH1F                       *fHistDeltaPtRCExLJ[4];       //!deltaPt = Pt(RC) - A * rhoKt, imposing min distance from leading jet
  TH1F                       *fHistRCPt[4];                //!Random cone pt
  TH1F                       *fHistRCPtExLJ[4];            //!Random cone pt, imposing min distance from leading jet
  TH1F                       *fHistEmbJets[4];             //!Pt distribution of embedded jets
  TH1F                       *fHistEmbPart[4];             //!Pt distribution of embedded particle
  TH1F                       *fHistDeltaPtMedKtEmb[4];     //!deltaPt = Pt(embjet) - Area(embjet) * rhoKt - Pt(embtrack)

 private:
  AliAnalysisTaskSAJF(const AliAnalysisTaskSAJF&);            // not implemented
  AliAnalysisTaskSAJF &operator=(const AliAnalysisTaskSAJF&); // not implemented

  ClassDef(AliAnalysisTaskSAJF, 3) // jet finder task
};
#endif
