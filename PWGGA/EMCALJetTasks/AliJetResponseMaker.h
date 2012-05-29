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

  void                        SetMCJetsName(const char *n)       { fMCJetsName    = n;  }
  void                        SetMCTracksName(const char *n)     { fMCTracksName  = n;  }

 protected:
  void                        DoJetLoop(Int_t &maxJetIndex, TClonesArray *jets, TClonesArray *tracks, TClonesArray *clusters = 0);
  void                        FillHistograms();
  void                        RetrieveEventObjects();

  TString                     fMCTracksName;              // name of MC particle collection
  TString                     fMCJetsName;                // name of MC jet collection
  TClonesArray               *fMCTracks;                  //!MC particles
  TClonesArray               *fMCJets;                    //!MC jets
  // Particle level jets
  TH2F                       *fHistMCJetPhiEta;           //!phi-eta distribution of jets
  TH1F                       *fHistMCJetsPt;              //!inclusive jet pt spectrum
  TH1F                       *fHistMCJetsPtTrack;         //!inclusive jet pt spectrum track biased
  TH1F                       *fHistMCJetsPtClus;          //!inclusive jet pt spectrum cluster biased
  TH1F                       *fHistMCJetsPtNonBias;       //!non biased inclusive jet pt spectrum
  TH1F                       *fHistMCLeadingJetPt;        //!leading jet pt spectrum
  TH2F                       *fHistMCJetsNEFvsPt;         //!jet neutral energy fraction vs. jet pt
  TH2F                       *fHistMCJetsZvsPt;           //!constituent Pt over Jet Pt ratio vs. jet pt
  // Detector level jets
  TH2F                       *fHistJetPhiEta;             //!phi-eta distribution of jets
  TH1F                       *fHistJetsPt;                //!inclusive jet pt spectrum
  TH1F                       *fHistJetsPtTrack;           //!inclusive jet pt spectrum track biased
  TH1F                       *fHistJetsPtClus;            //!inclusive jet pt spectrum cluster biased
  TH1F                       *fHistJetsPtNonBias;         //!non biased inclusive jet pt spectrum
  TH1F                       *fHistLeadingJetPt;          //!leading jet pt spectrum
  TH2F                       *fHistJetsNEFvsPt;           //!jet neutral energy fraction vs. jet pt
  TH2F                       *fHistJetsZvsPt;             //!constituent Pt over Jet Pt ratio vs. jet pt

 private:
  AliJetResponseMaker(const AliJetResponseMaker&);            // not implemented
  AliJetResponseMaker &operator=(const AliJetResponseMaker&); // not implemented

  ClassDef(AliJetResponseMaker, 1) // Jet response matrix producing task
};
#endif
