#ifndef ALIANALYSISTASKEMCALJETSAMPLE_H
#define ALIANALYSISTASKEMCALJETSAMPLE_H

// $Id$

class TH1;
class TH2;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalJetSample : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskEmcalJetSample();
  AliAnalysisTaskEmcalJetSample(const char *name);
  virtual ~AliAnalysisTaskEmcalJetSample();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

 protected:
  void                        ExecOnce();
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;

  // General histograms
  TH1                        *fHistTracksPt[4];            //!Track pt spectrum
  TH1                        *fHistClustersPt[4];          //!Cluster pt spectrum
  TH1                        *fHistLeadingJetPt[4];        //!Leading jet pt spectrum
  TH2                        *fHistJetsPhiEta[4];          //!Phi-Eta distribution of jets
  TH2                        *fHistJetsPtArea[4];          //!Jet pt vs. area
  TH2                        *fHistJetsPtLeadHad[4];       //!Jet pt vs. leading hadron
  TH2                        *fHistJetsCorrPtArea[4];      //!Jet pt - bkg vs. area

  AliJetContainer            *fJetsCont;                   //!Jets
  AliParticleContainer       *fTracksCont;                 //!Tracks
  AliClusterContainer        *fCaloClustersCont;           //!Clusters  

 private:
  AliAnalysisTaskEmcalJetSample(const AliAnalysisTaskEmcalJetSample&);            // not implemented
  AliAnalysisTaskEmcalJetSample &operator=(const AliAnalysisTaskEmcalJetSample&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetSample, 2) // jet sample analysis task
};
#endif
