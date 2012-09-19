#ifndef ALIANALYSISTASKDELTAPT_H
#define ALIANALYSISTASKDELTAPT_H

// $Id$

class TClonesArray;
class TString;
class TH1F;
class TH2F;
class AliRhoParameter;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskDeltaPt : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskDeltaPt();
  AliAnalysisTaskDeltaPt(const char *name);
  virtual ~AliAnalysisTaskDeltaPt();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  void                        SetJetMinRC2LJ(Float_t d)                            { fMinRC2LJ                = d          ; }
  void                        SetEmbJetsName(const char *n)                        { fEmbJetsName             = n          ; } 
  void                        SetEmbTracksName(const char *n)                      { fEmbTracksName           = n          ; }
  void                        SetEmbClusName(const char *n)                        { fEmbCaloName             = n          ; }
  void                        SetRandTracksName(const char *n)                     { fRandTracksName          = n          ; }
  void                        SetRandClusName(const char *n)                       { fRandCaloName            = n          ; }
  void                        SetMC(Bool_t m)                                      { fMCAna                   = m          ; }
  void                        SetRCperEvent(Int_t n)                               { fRCperEvent              = n          ; }

 protected:
  void                        ExecOnce()                                                                                    ;
  Bool_t                      FillHistograms()                                                                              ;
  void                        GetLeadingJets(Int_t &maxJetIndex, Int_t &max2JetIndex)                                       ;
  void                        DoEmbJetLoop(AliEmcalJet* &embJet, TObject* &embPart)                                         ;
  void                        DoEmbTrackLoop()                                                                              ;
  void                        DoEmbClusterLoop()                                                                            ;
  void                        GetRandomCone(Float_t &pt, Float_t &eta, Float_t &phi, 
					    AliEmcalJet *jet = 0, TClonesArray* tracks = 0, TClonesArray* clusters = 0) const;

  Bool_t                      fMCAna;                      // =true MC analysis (toy model)
  Float_t                     fMinRC2LJ;                   // Minimum distance random cone to leading jet
  TString                     fEmbJetsName;                // Name of embedded jet collection
  TString                     fEmbTracksName;              // Name of embedded track collection
  TString                     fEmbCaloName;                // Name of embedded calo cluster collection
  TString                     fRandTracksName;             // Name of randomized track collection
  TString                     fRandCaloName;               // Name of randomized calo cluster collection
  Int_t                       fRCperEvent;                 // No. of random cones per event

  TClonesArray               *fEmbJets;                    //!Embedded jets
  TClonesArray               *fEmbTracks;                  //!Embedded tracks
  TClonesArray               *fEmbCaloClusters;            //!Embedded clusters  
  TClonesArray               *fRandTracks;                 //!Randomized tracks
  TClonesArray               *fRandCaloClusters;           //!Randomized clusters
  Int_t                       fEmbeddedClusterId;          //!Embedded cluster id
  Int_t                       fEmbeddedTrackId;            //!Embedded track id

  // Random cones
  TH2F                       *fHistRCPhiEta;               //!Phi-Eta distribution of random cones
  TH1F                       *fHistRCPt[4];                //!Random cone pt
  TH1F                       *fHistRCPtExLJ[4];            //!Random cone pt, imposing min distance from leading jet
  TH1F                       *fHistRCPtRand[4];            //!Random cone pt, randomized particles
  TH2F                       *fHistRCPtExLJVSDPhiLJ;       //!Random cone pt, imposing min distance from leading jet, vs. deltaPhi leading jet
  TH2F                       *fHistRhoVSRCPt[4];           //!Area(RC) * rho vs. Pt(RC)
  TH1F                       *fHistDeltaPtRC[4];           //!deltaPt = Pt(RC) - A * rho
  TH1F                       *fHistDeltaPtRCExLJ[4];       //!deltaPt = Pt(RC) - A * rho, imposing min distance from leading jet
  TH1F                       *fHistDeltaPtRCRand[4];       //!deltaPt = Pt(RC) - A * rho, randomzied particles

  // Jet embedding
  TH2F                       *fHistEmbNotFoundPhiEta[4];   //!Phi-Eta of "not found" embedded particles
  TH1F                       *fHistEmbJetsPt[4];           //!Pt distribution of embedded jets
  TH1F                       *fHistEmbJetsCorrPt[4];       //!Pt-rho*A distribution of embedded jets
  TH1F                       *fHistEmbJetsArea[4];         //!Area distribution of embedded jets
  TH1F                       *fHistEmbPartPt[4];           //!Pt distribution of embedded particle
  TH2F                       *fHistEmbJetPhiEta;           //!Phi-Eta distribution of embedded jets
  TH2F                       *fHistEmbPartPhiEta;          //!Phi-Eta distribution of embedded particles
  TH1F                       *fHistDistEmbPartJetAxis[4];  //!Distance between embedded particle and jet axis
  TH2F                       *fHistRhoVSEmbBkg[4];         //!Area(embjet) * rho vs. Pt(embjet) - Pt(embtrack)
  TH1F                       *fHistDeltaPtEmb[4];          //!deltaPt = Pt(embjet) - Area(embjet) * rho - Pt(embtrack)

 private:
  AliAnalysisTaskDeltaPt(const AliAnalysisTaskDeltaPt&);            // not implemented
  AliAnalysisTaskDeltaPt &operator=(const AliAnalysisTaskDeltaPt&); // not implemented

  ClassDef(AliAnalysisTaskDeltaPt, 1) // deltaPt analysis task
};
#endif
