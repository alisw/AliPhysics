#ifndef ALIJETRESPONSEMAKER_H
#define ALIJETRESPONSEMAKER_H

// $Id$

class TClonesArray;
class TH1F;
class TH2F;

#include "AliAnalysisTaskEmcal.h"

class AliJetResponseMaker : public AliAnalysisTaskEmcal {
 public:
  AliJetResponseMaker();
  AliJetResponseMaker(const char *name);
  virtual ~AliJetResponseMaker();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  void                        SetMCJetsName(const char *n)       { fMCJetsName    = n; }
  void                        SetMCTracksName(const char *n)     { fMCTracksName  = n; }
  void                        SetMaxDistance(Double_t d)         { fMaxDistance   = d; }

 protected:
  void                        DoJetLoop(TClonesArray *jets1, TClonesArray *jets2);
  void                        FillHistograms();
  void                        RetrieveEventObjects();

  TString                     fMCTracksName;              // name of MC particle collection
  TString                     fMCJetsName;                // name of MC jet collection
  Double_t                    fMaxDistance;               // maximum distance between matched particle and detector level jets
  TClonesArray               *fMCTracks;                  //!MC particles
  TClonesArray               *fMCJets;                    //!MC jets
  // Particle level jets
  TH2F                       *fHistMCJetPhiEta;           //!phi-eta distribution of jets
  TH1F                       *fHistMCJetsPt;              //!inclusive jet pt spectrum
  TH1F                       *fHistMCJetsPtNonBias;       //!non biased inclusive jet pt spectrum
  TH2F                       *fHistMCJetsNEFvsPt;         //!jet neutral energy fraction vs. jet pt
  TH2F                       *fHistMCJetsZvsPt;           //!constituent Pt over Jet Pt ratio vs. jet pt
  // Detector level jets
  TH2F                       *fHistJetPhiEta;             //!phi-eta distribution of jets
  TH1F                       *fHistJetsPt;                //!inclusive jet pt spectrum
  TH1F                       *fHistJetsPtNonBias;         //!non biased inclusive jet pt spectrum
  TH2F                       *fHistJetsNEFvsPt;           //!jet neutral energy fraction vs. jet pt
  TH2F                       *fHistJetsZvsPt;             //!constituent Pt over Jet Pt ratio vs. jet pt
  // Detector-particle level matching
  TH1F                       *fHistClosestDistance;       //!distance between closest particle to detector level jet
  TH1F                       *fHistClosestDeltaPhi;       //!delta phi between closest particle to detector level jet
  TH1F                       *fHistClosestDeltaEta;       //!delta eta between closest particle to detector level jet
  TH1F                       *fHistClosestDeltaPt;        //!delta pt between closest particle to detector level jet
  TH2F                       *fHistPartvsDetecPt;         //!particle vs detector level jet pt

 private:
  AliJetResponseMaker(const AliJetResponseMaker&);            // not implemented
  AliJetResponseMaker &operator=(const AliJetResponseMaker&); // not implemented

  ClassDef(AliJetResponseMaker, 1) // Jet response matrix producing task
};
#endif
