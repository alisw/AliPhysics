#ifndef ALIANALYSISTASKEMCALJETSAMPLEDEV_H
#define ALIANALYSISTASKEMCALJETSAMPLEDEV_H

// $Id$

class TH1;
class TH2;

#include "AliAnalysisTaskEmcalJetDev.h"

class AliAnalysisTaskEmcalJetSampleDev : public AliAnalysisTaskEmcalJetDev {
 public:

  AliAnalysisTaskEmcalJetSampleDev();
  AliAnalysisTaskEmcalJetSampleDev(const char *name);
  virtual ~AliAnalysisTaskEmcalJetSampleDev();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

 protected:
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

 private:
  AliAnalysisTaskEmcalJetSampleDev(const AliAnalysisTaskEmcalJetSampleDev&);            // not implemented
  AliAnalysisTaskEmcalJetSampleDev &operator=(const AliAnalysisTaskEmcalJetSampleDev&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetSampleDev, 1) // jet sample analysis task
};
#endif
