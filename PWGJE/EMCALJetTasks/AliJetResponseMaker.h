#ifndef ALIJETRESPONSEMAKER_H
#define ALIJETRESPONSEMAKER_H

// $Id$

class AliGenPythiaEventHeader;
class TClonesArray;
class TH1;
class TH2;
class TH3;

#include "AliAnalysisTaskEmcalJet.h"

class AliJetResponseMaker : public AliAnalysisTaskEmcalJet {
 public:
  AliJetResponseMaker();
  AliJetResponseMaker(const char *name);
  virtual ~AliJetResponseMaker();

  void                        UserCreateOutputObjects();

  void                        SetMCJetsName(const char *n)       { fMCJetsName      = n; }
  void                        SetMCTracksName(const char *n)     { fMCTracksName    = n; }
  void                        SetMaxDistance(Double_t d)         { fMaxDistance     = d; }
  void                        SetDoWeighting(Bool_t d = kTRUE)   { fDoWeighting     = d; }
  void                        SetMCFiducial(Bool_t f = kTRUE)    { fMCFiducial      = f; }
  void                        SetEventWeightHist(Bool_t h)       { fEventWeightHist = h; }
  void                        SetPtHardBin(Int_t b)              { fSelectPtHardBin = b; }
  void                        SetDoMatching(Bool_t b)            { fDoMatching      = b; }

 protected:
  Bool_t                      IsEventSelected();
  Bool_t                      AcceptJet(AliEmcalJet* jet)   const;
  void                        ExecOnce();
  void                        DoJetLoop(TClonesArray *jets1, TClonesArray *jets2, Bool_t mc);
  Bool_t                      FillHistograms();
  Bool_t                      RetrieveEventObjects();
  Bool_t                      Run();

  TString                     fMCTracksName;                  // name of MC particle collection
  TString                     fMCJetsName;                    // name of MC jet collection
  Double_t                    fMaxDistance;                   // maximum distance between matched particle and detector level jets
  Bool_t                      fDoWeighting;                   // = true, weight using trials and given x section
  Bool_t                      fEventWeightHist;               // = true create event weight histogram
  Bool_t                      fMCFiducial;                    // = true MC jets in fiducial acceptance
  Double_t                    fMCminEta;                      // MC jets minimum eta
  Double_t                    fMCmaxEta;                      // MC jets maximum eta
  Double_t                    fMCminPhi;                      // MC jets minimum phi
  Double_t                    fMCmaxPhi;                      // MC jets maximum phi
  Int_t                       fSelectPtHardBin;               // select one pt hard bin for analysis
  Bool_t                      fDoMatching;                    // whether or not it should run the matching between MC and rec jets

  AliGenPythiaEventHeader    *fPythiaHeader;                  //!event Pythia header
  Double_t                    fEventWeight;                   //!event weight
  Int_t                       fPtHardBin;                     //!event pt hard bin
  Int_t                       fNTrials;                       //!event trials
  TClonesArray               *fMCTracks;                      //!MC particles
  TClonesArray               *fMCJets;                        //!MC jets
  // General histograms
  TH1                        *fHistNTrials;                   //!total number of trials per pt hard bin
  TH1                        *fHistEvents;                    //!total number of events per pt hard bin
  TH1                        *fHistEventWeight[11];           //!event weight
  // Particle level jets
  TH2                        *fHistMCJetsPhiEta;              //!phi-eta distribution of jets
  TH2                        *fHistMCJetsPtArea;              //!inclusive jet pt vs area histogram
  TH2                        *fHistMCJetsPhiEtaFiducial;      //!phi-eta distribution of jets in fiducial acceptance (plus lead hadron bias)
  TH2                        *fHistMCJetsPtAreaFiducial;      //!inclusive jet pt spectrum in fiducial acceptance (plus lead hadron bias)
  TH2                        *fHistMCJetsNEFvsPt;             //!jet neutral energy fraction vs. jet pt
  TH2                        *fHistMCJetsZvsPt;               //!constituent Pt over Jet Pt ratio vs. jet pt
  // Detector level jets
  TH2                        *fHistJetsPhiEta;                //!phi-eta distribution of jets
  TH2                        *fHistJetsPtArea;                //!inclusive jet pt vs. area histogram
  TH2                        *fHistJetsCorrPtArea;            //!inclusive jet pt vs. area histogram
  TH2                        *fHistJetsNEFvsPt;               //!jet neutral energy fraction vs. jet pt
  TH2                        *fHistJetsZvsPt;                 //!constituent Pt over Jet Pt ratio vs. jet pt
  // Detector-particle level matching
  TH2                        *fHistMatchingLevelMCPt;         //!matching level vs. particle level pt
  TH3                        *fHistClosestDeltaEtaPhiMCPt;    //!delta eta-phi between closest particle level to detector level jet vs. particle level pt
  TH2                        *fHistClosestDeltaPtMCPt;        //!delta pt between closest particle level to detector level jet vs. particle level pt
  TH2                        *fHistClosestDeltaCorrPtMCPt;    //!delta pt between closest particle level to detector level jet vs. particle level pt
  TH2                        *fHistNonMatchedMCJetsPtArea;    //!non-matched mc jet pt distribution
  TH2                        *fHistNonMatchedJetsPtArea;      //!non-matched jet pt distribution
  TH2                        *fHistNonMatchedJetsCorrPtArea;  //!non-matched jet pt distribution
  TH2                        *fHistPartvsDetecPt;             //!particle vs detector level jet pt
  TH2                        *fHistPartvsDetecCorrPt;         //!particle vs detector level jet pt
  TH2                        *fHistMissedMCJetsPtArea;         //!mc jets not measured

 private:
  AliJetResponseMaker(const AliJetResponseMaker&);            // not implemented
  AliJetResponseMaker &operator=(const AliJetResponseMaker&); // not implemented

  ClassDef(AliJetResponseMaker, 9) // Jet response matrix producing task
};
#endif
