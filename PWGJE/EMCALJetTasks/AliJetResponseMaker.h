#ifndef ALIJETRESPONSEMAKER_H
#define ALIJETRESPONSEMAKER_H

// $Id$

class AliGenPythiaEventHeader;
class TClonesArray;
class TH1F;
class TH2F;

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
  Bool_t                      AcceptJet(AliEmcalJet* jet, Bool_t /*bias*/ = kFALSE, Bool_t /*upCut*/ = kFALSE)   const;
  void                        ExecOnce();
  void                        DoJetLoop(TClonesArray *jets1, TClonesArray *jets2, Bool_t mc);
  Bool_t                      FillHistograms();
  Bool_t                      RetrieveEventObjects();
  Bool_t                      Run();

  TString                     fMCTracksName;              // name of MC particle collection
  TString                     fMCJetsName;                // name of MC jet collection
  Double_t                    fMaxDistance;               // maximum distance between matched particle and detector level jets
  Bool_t                      fDoWeighting;               // = true, weight using trials and given x section
  Bool_t                      fEventWeightHist;           // = true create event weight histogram
  Bool_t                      fMCFiducial;                // = true MC jets in fiducial acceptance
  Double_t                    fMCminEta;                  // MC jets minimum eta
  Double_t                    fMCmaxEta;                  // MC jets maximum eta
  Double_t                    fMCminPhi;                  // MC jets minimum phi
  Double_t                    fMCmaxPhi;                  // MC jets maximum phi
  Int_t                       fSelectPtHardBin;           // select one pt hard bin for analysis
  Bool_t                      fDoMatching;                // whether or not it should run the matching between MC and rec jets

  AliGenPythiaEventHeader    *fPythiaHeader;              //!event Pythia header
  Double_t                    fEventWeight;               //!event weight
  Int_t                       fPtHardBin;                 //!event pt hard bin
  Int_t                       fNTrials;                   //!event trials
  TClonesArray               *fMCTracks;                  //!MC particles
  TClonesArray               *fMCJets;                    //!MC jets
  // General histograms
  TH1F                       *fHistZVertex;               //!Z vertex position
  TH1F                       *fHistNTrials;               //!total number of trials per pt hard bin
  TH1F                       *fHistEvents;                //!total number of events per pt hard bin
  TH1F                       *fHistEventWeight[11];       //!event weight
  // Particle level jets
  TH2F                       *fHistMCJetPhiEta;           //!phi-eta distribution of jets
  TH1F                       *fHistMCJetsPt;              //!inclusive jet pt spectrum
  TH2F                       *fHistMCJetPhiEtaFiducial;   //!phi-eta distribution of jets in fiducial acceptance (plus lead hadron bias)
  TH1F                       *fHistMCJetsPtFiducial;      //!inclusive jet pt spectrum in fiducial acceptance (plus lead hadron bias)
  TH2F                       *fHistMCJetsNEFvsPt;         //!jet neutral energy fraction vs. jet pt
  TH2F                       *fHistMCJetsZvsPt;           //!constituent Pt over Jet Pt ratio vs. jet pt
  // Detector level jets
  TH2F                       *fHistJetPhiEta;             //!phi-eta distribution of jets
  TH1F                       *fHistJetsPt;                //!inclusive jet pt spectrum
  TH2F                       *fHistJetsNEFvsPt;           //!jet neutral energy fraction vs. jet pt
  TH2F                       *fHistJetsZvsPt;             //!constituent Pt over Jet Pt ratio vs. jet pt
  // Detector-particle level matching
  TH1F                       *fHistClosestDistance;       //!distance between closest particle to detector level jet
  TH1F                       *fHistClosestDeltaPhi;       //!delta phi between closest particle to detector level jet
  TH1F                       *fHistClosestDeltaEta;       //!delta eta between closest particle to detector level jet
  TH1F                       *fHistClosestDeltaPt;        //!delta pt between closest particle to detector level jet
  TH1F                       *fHistNonMatchedMCJetPt;     //!non-matched mc jet pt distribution
  TH1F                       *fHistNonMatchedJetPt;       //!non-matched jet pt distribution
  TH2F                       *fHistPartvsDetecPt;         //!particle vs detector level jet pt
  TH1F                       *fHistMissedMCJets;          //!mc jets not measured

 private:
  AliJetResponseMaker(const AliJetResponseMaker&);            // not implemented
  AliJetResponseMaker &operator=(const AliJetResponseMaker&); // not implemented

  ClassDef(AliJetResponseMaker, 7) // Jet response matrix producing task
};
#endif
