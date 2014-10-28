#ifndef ALIANALYSISTASKDELTAPT_H
#define ALIANALYSISTASKDELTAPT_H

// $Id$

class TClonesArray;
class TString;
class TH1;
class TH2;
class TH3;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskDeltaPt : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskDeltaPt();
  AliAnalysisTaskDeltaPt(const char *name);
  virtual ~AliAnalysisTaskDeltaPt() {;}

  void                        UserCreateOutputObjects();

  void                        SetJetMinRC2LJ(Float_t d)                            { fMinRC2LJ                = d          ; }
  void                        SetRCperEvent(Int_t n)                               { fRCperEvent              = n          ; }
  void                        SetMCJetPtThreshold(Double_t t)                      { fMCJetPtThreshold        = t          ; }
  void                        SetConeRadius(Double_t r)                            { fConeRadius              = r          ; }
  void                        SetConeEtaPhiEMCAL() ;
  void                        SetConeEtaPhiTPC()   ;
  void                        SetConeEtaLimits(Float_t min, Float_t max)           { fConeMinEta = min, fConeMaxEta = max  ; }
  void                        SetConePhiLimits(Float_t min, Float_t max)           { fConeMinPhi = min, fConeMaxPhi = max  ; }

 protected:
  void                        AllocateHistogramArrays()                                                                     ;
  void                        ExecOnce()                                                                                    ;
  Bool_t                      FillHistograms()                                                                              ;
  void                        GetLeadingJets(Int_t &maxJetIndex, Int_t &max2JetIndex)                                       ;
  AliEmcalJet*                NextEmbeddedJet(Bool_t reset=kFALSE)                                                          ;
  void                        DoEmbTrackLoop()                                                                              ;
  void                        DoEmbClusterLoop()                                                                            ;
  void                        GetRandomCone(Float_t &pt, Float_t &eta, Float_t &phi, AliParticleContainer* tracks, AliClusterContainer* clusters,
					    AliEmcalJet *jet = 0, Bool_t bPartialExclusion = 0) const;
  Double_t                    GetNColl() const;


  Double_t                    fMCJetPtThreshold;           // threshold for MC jets
  Float_t                     fMinRC2LJ;                   // Minimum distance random cone to leading jet
  Int_t                       fRCperEvent;                 // No. of random cones per event
  Double_t                    fConeRadius;                 // Radius of the random cones
  Float_t                     fConeMinEta;                 // Minimum eta of the random cones
  Float_t                     fConeMaxEta;                 // Maximum eta of the random cones
  Float_t                     fConeMinPhi;                 // Minimum phi of the random cones
  Float_t                     fConeMaxPhi;                 // Maximum phi of the random cones

  AliJetContainer            *fJetsCont;                   //!Jets
  AliParticleContainer       *fTracksCont;                 //!Tracks
  AliClusterContainer        *fCaloClustersCont;           //!Clusters  
  AliJetContainer            *fEmbJetsCont;                //!Embedded jets
  AliParticleContainer       *fEmbTracksCont;              //!Embedded tracks
  AliClusterContainer        *fEmbCaloClustersCont;        //!Embedded clusters  
  AliParticleContainer       *fRandTracksCont;             //!Randomized tracks
  AliClusterContainer        *fRandCaloClustersCont;       //!Randomized clusters

  // General
  TH2                        *fHistRhovsCent;              //!Rho vs. centrality

  // Random cones
  TH2                        *fHistRCPhiEta;               //!Phi-Eta distribution of random cones
  TH1                       **fHistRCPt;                   //!Random cone pt
  TH1                       **fHistRCPtExLJ;               //!Random cone pt, imposing min distance from leading jet
  TH1                       **fHistRCPtExPartialLJ;        //!Random cone pt, imposing min distance from leading jet with 1/ncoll probability
  TH1                       **fHistRCPtRand;               //!Random cone pt, randomized particles
  TH2                       **fHistRhoVSRCPt;              //!Area(RC) * rho vs. Pt(RC)
  TH2                       **fHistDeltaPtRCvsEP;          //!deltaPt = Pt(RC) - A * rho vs. event plane
  TH1                       **fHistDeltaPtRCExLJ;          //!deltaPt = Pt(RC) - A * rho, imposing min distance from leading jet
  TH1                       **fHistDeltaPtRCExPartialLJ;   //!deltaPt = Pt(RC) - A * rho, imposing min distance from leading jet with 1/ncoll probability
  TH1                       **fHistDeltaPtRCRand;          //!deltaPt = Pt(RC) - A * rho, randomzied particles

  // Jet embedding
  TH3                       **fHistEmbJetsPtArea;          //!Pt vs. area of embedded jets
  TH3                       **fHistEmbJetsCorrPtArea;      //!Pt-rho*A vs. area of embedded jets
  TH2                       **fHistEmbPartPtvsJetPt;       //!MC jet pt total jet pt
  TH2                       **fHistEmbPartPtvsJetCorrPt;   //!MC jet pt total jet pt - rho*A
  TH2                       **fHistJetPtvsJetCorrPt;       //!Pt vs jet pt - rho*A
  TH1                       **fHistDistLeadPart2JetAxis;   //!Distance between leading particle and jet axis
  TH2                       **fHistEmbBkgArea;             //!Pt(embjet) - Pt(embtrack) vs. area of embedded jets
  TH2                       **fHistRhoVSEmbBkg;            //!Area(embjet) * rho vs. Pt(embjet) - Pt(embtrack)
  TH2                       **fHistDeltaPtEmbArea;         //!deltaPt = Pt(embjet) - Area(embjet) * rho - Pt(embtrack) vs. Area(embjet)
  TH2                       **fHistDeltaPtEmbvsEP;         //!deltaPt = Pt(embjet) - Area(embjet) * rho - Pt(embtrack) vs. event plane
  TH2                        *fHistRCPtExLJVSDPhiLJ;       //!Random cone pt, imposing min distance from leading jet, vs. deltaPhi leading jet
  TH2                        *fHistRCPtExPartialLJVSDPhiLJ;//!Random cone pt, imposing min distance from leading jet, vs. deltaPhi leading jet with 1/ncoll probability
  TH2                        *fHistEmbJetsPhiEta;          //!Phi-Eta distribution of embedded jets
  TH2                        *fHistLeadPartPhiEta;         //!Phi-Eta distribution of the leading particle of embedded jets

 private:
  AliAnalysisTaskDeltaPt(const AliAnalysisTaskDeltaPt&);            // not implemented
  AliAnalysisTaskDeltaPt &operator=(const AliAnalysisTaskDeltaPt&); // not implemented

  ClassDef(AliAnalysisTaskDeltaPt, 5) // deltaPt analysis task
};
#endif
