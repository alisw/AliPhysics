#ifndef ALIANALYSISTASKDELTAPT_H
#define ALIANALYSISTASKDELTAPT_H

// $Id$

class TClonesArray;
class TString;
class TH1;
class TH2;
class TH3;

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
  void                        SetRCperEvent(Int_t n)                               { fRCperEvent              = n          ; }
  void                        SetMCJetPtThreshold(Double_t t)                      { fMCJetPtThreshold        = t          ; }

 protected:
  void                        ExecOnce()                                                                                    ;
  Bool_t                      FillHistograms()                                                                              ;
  void                        GetLeadingJets(Int_t &maxJetIndex, Int_t &max2JetIndex)                                       ;
  AliEmcalJet*                NextEmbeddedJet(Int_t i=-1)                                                                   ;
  void                        DoEmbTrackLoop()                                                                              ;
  void                        DoEmbClusterLoop()                                                                            ;
  void                        GetRandomCone(Float_t &pt, Float_t &eta, Float_t &phi, 
					    AliEmcalJet *jet = 0, TClonesArray* tracks = 0, TClonesArray* clusters = 0, Bool_t bPartialExclusion = 0) const;
  Double_t                    GetNColl() const;


  Double_t                    fMCJetPtThreshold;           // threshold for MC jets
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
  Int_t                       fEmbeddedClusterNIds;        //!Embedded cluster id count
  Int_t                       fEmbeddedClusterIds[999];    //!Embedded cluster ids
  Int_t                       fEmbeddedTrackNIds;          //!Embedded track id count
  Int_t                       fEmbeddedTrackIds[999];      //!Embedded track ids

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
  TH1                       **fHistEmbNotFoundPt;          //!Pt of "not found" embedded particles
  TH2                       **fHistEmbNotFoundPhiEta;      //!Phi-Eta of "not found" embedded particles
  TH2                       **fHistEmbRejectedJetsPhiEta;  //!Phi-Eta of rejected embedded jets
  TH1                       **fHistEmbRejectedJetsPtArea;  //!Pt-area of rejected embedded jets
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
  TH2                        *fHistEmbJetsPhiEta;          //!Phi-Eta distribution of embedded jets<
  TH2                        *fHistLeadPartPhiEta;         //!Phi-Eta distribution of the leading particle of embedded jets

 private:
  AliAnalysisTaskDeltaPt(const AliAnalysisTaskDeltaPt&);            // not implemented
  AliAnalysisTaskDeltaPt &operator=(const AliAnalysisTaskDeltaPt&); // not implemented

  ClassDef(AliAnalysisTaskDeltaPt, 4) // deltaPt analysis task
};
#endif
