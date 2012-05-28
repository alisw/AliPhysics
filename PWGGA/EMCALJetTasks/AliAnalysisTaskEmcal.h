#ifndef ALIANALYSISTASKEMCAL_H
#define ALIANALYSISTASKEMCAL_H

// $Id: AliAnalysisTaskEmcal.h 56670 2012-05-24 13:24:04Z loizides $

class TClonesArray;
class TString;
class AliVParticle;
class AliVCluster;
class AliEmcalJet;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEmcal : public AliAnalysisTaskSE {
 public:
  
  enum EmcalAnaType {
    kTPC       = 0,     // TPC only analysis
    kEMCAL     = 1,     // EMCal + TPC analysis
  };

  AliAnalysisTaskEmcal();
  AliAnalysisTaskEmcal(const char *name);
  virtual ~AliAnalysisTaskEmcal();

  virtual void                UserCreateOutputObjects();
  virtual void                UserExec(Option_t *option);
  virtual void                Terminate(Option_t *option);
  virtual void                Init();

  void                        SetTracksName(const char *n)                         { fTracksName     = n          ; }
  void                        SetClusName(const char *n)                           { fCaloName       = n          ; }
  void                        SetJetsName(const char *n)                           { fJetsName       = n          ; }
  void                        SetAnaType(EmcalAnaType type)                        { fAnaType        = type       ; }
  void                        SetJetRadius(Float_t r)                              { fJetRadius      = r          ; } 
  void                        SetPtCut(Float_t cut)                                { fPtCut          = cut        ; }
  void                        SetPtBiasJetTrack(Float_t b)                         { fPtBiasJetTrack = b          ; }
  void                        SetPtBiasJetClus(Float_t b)                          { fPtBiasJetClus  = b          ; }
  void                        SetHistoBins(Int_t nbins, Float_t min, Float_t max)  { fNbins = nbins; fMinPt = min; fMaxPt = max; }
  void                        SetEtaLimits(Float_t min, Float_t max)               { fMinEta = min, fMaxEta = max ; }
  void                        SetPhiLimits(Float_t min, Float_t max)               { fMinPhi = min, fMaxPhi = max ; }
  void                        SetInitialized(Bool_t ini = kTRUE)                   { fInitialized    = ini        ; }
  void                        SetJetPtCut(Float_t cut)                             { fJetPtCut       = cut        ; }
  void                        SetJetAreaCut(Float_t cut)                           { fJetAreaCut     = cut        ; }

 protected:

  Bool_t                      AcceptTrack(AliVParticle* track, Bool_t acceptMC = kFALSE)           const;
  Bool_t                      AcceptCluster(AliVCluster* clus, Bool_t acceptMC = kFALSE)           const;
  Bool_t                      AcceptJet(AliEmcalJet* jet)                                          const;
  Bool_t                      IsJetTrack(AliEmcalJet* jet, Int_t itrack, Bool_t sorted = kTRUE)    const;
  Bool_t                      IsJetCluster(AliEmcalJet* jet, Int_t iclus, Bool_t sorted = kTRUE)   const;

  virtual void                RetrieveEventObjects()        ;
  virtual void                FillHistograms()           = 0; 

  EmcalAnaType                fAnaType;                    // Analysis type
  Bool_t                      fInitialized;                // Whether or not the task has been already initialized
  Float_t                     fMinEta;                     // Minimum eta accepatance
  Float_t                     fMaxEta;                     // Maximum eta accepatance
  Float_t                     fMinPhi;                     // Minimum phi accepatance
  Float_t                     fMaxPhi;                     // Maximum phi accepatance  
  Float_t                     fJetRadius;                  // Jet radius
  TString                     fTracksName;                 // Name of track collection
  TString                     fCaloName;                   // Name of calo cluster collection
  TString                     fJetsName;                   // Name of jet collection
  Int_t                       fNbins;                      // No. of pt bins
  Float_t                     fMinPt;                      // Min pt in histograms
  Float_t                     fMaxPt;                      // Max pt in histograms
  Float_t                     fPtCut;                      // Cut on particle pt
  Float_t                     fPtBiasJetTrack;             // Select jets with a minimum pt track
  Float_t                     fPtBiasJetClus;              // Select jets with a minimum pt cluster
  Float_t                     fJetPtCut;                   // Cut on jet pt
  Float_t                     fJetAreaCut;                 // Cut on jet area

  TClonesArray               *fTracks;                     //!Tracks
  TClonesArray               *fCaloClusters;               //!Clusters
  TClonesArray               *fJets;                       //!Jets
  Float_t                     fCent;                       //!Event centrality
  Int_t                       fCentBin;                    //!Event centrality bin
  Double_t                    fVertex[3];                  //!Event vertex

  TList                      *fOutput;                     //!Output list

 private:
  AliAnalysisTaskEmcal(const AliAnalysisTaskEmcal&);            // not implemented
  AliAnalysisTaskEmcal &operator=(const AliAnalysisTaskEmcal&); // not implemented

  ClassDef(AliAnalysisTaskEmcal, 0) // emcal base analysis task
};
#endif
