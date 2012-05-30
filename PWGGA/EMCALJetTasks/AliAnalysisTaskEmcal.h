#ifndef ALIANALYSISTASKEMCAL_H
#define ALIANALYSISTASKEMCAL_H

// $Id$

class TClonesArray;
class TString;
class TList;
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

  enum BeamType {
    kNA       = -1,
    kpp       = 0,
    kAA       = 1,
    kpA       = 2
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
  Int_t                       GetBeamType()                                                             ;

  virtual void                RetrieveEventObjects()        ;
  virtual void                FillHistograms()             {;} 

  EmcalAnaType                fAnaType;                    // analysis type
  Bool_t                      fInitialized;                // whether or not the task has been already initialized
  Float_t                     fMinEta;                     // minimum eta accepatance
  Float_t                     fMaxEta;                     // maximum eta accepatance
  Float_t                     fMinPhi;                     // minimum phi accepatance
  Float_t                     fMaxPhi;                     // maximum phi accepatance  
  Float_t                     fJetRadius;                  // jet radius
  TString                     fTracksName;                 // name of track collection
  TString                     fCaloName;                   // name of calo cluster collection
  TString                     fJetsName;                   // name of jet collection
  Int_t                       fNbins;                      // no. of pt bins
  Float_t                     fMinPt;                      // min pt in histograms
  Float_t                     fMaxPt;                      // max pt in histograms
  Float_t                     fPtCut;                      // cut on particle pt
  Float_t                     fPtBiasJetTrack;             // select jets with a minimum pt track
  Float_t                     fPtBiasJetClus;              // select jets with a minimum pt cluster
  Float_t                     fJetPtCut;                   // cut on jet pt
  Float_t                     fJetAreaCut;                 // cut on jet area

  TClonesArray               *fTracks;                     //!tracks
  TClonesArray               *fCaloClusters;               //!clusters
  TClonesArray               *fJets;                       //!jets
  Float_t                     fCent;                       //!event centrality
  Int_t                       fCentBin;                    //!event centrality bin
  Double_t                    fVertex[3];                  //!event vertex

  TList                      *fOutput;                     //!output list

 private:
  AliAnalysisTaskEmcal(const AliAnalysisTaskEmcal&);            // not implemented
  AliAnalysisTaskEmcal &operator=(const AliAnalysisTaskEmcal&); // not implemented

  ClassDef(AliAnalysisTaskEmcal, 1) // EMCAL base analysis task
};
#endif
