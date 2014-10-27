#ifndef ALIANALYSISTASKEMCALHIGHMULTTRIGGER_H
#define ALIANALYSISTASKEMCALHIGHMULTTRIGGER_H

// $Id$

class TH1;
class TH2;
class TH3;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalHighMultTrigger : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskEmcalHighMultTrigger();
  AliAnalysisTaskEmcalHighMultTrigger(const char *name);
  virtual ~AliAnalysisTaskEmcalHighMultTrigger();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  //Setters
  void                        SetNExcludeLeadingPatches(Int_t n)  {fNExLP = n; }

 protected:
  void                        ExecOnce();
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;
  
 private:
  Int_t                       fNExLP;                  //nr of leading patched to exclude from estimate
  Int_t                       fNAccPatches;            //nr of accepted patches
  Double_t                    fMedianEnergy;           //median event energy
  Double_t                    fMedianEnergyExLP;       //median event energy
  Double_t                    fSumEnergy;              //summed energy
  Double_t                    fSumEnergyExLP;          //summed energy


  //Histograms
  TH2F                       *fHistPatchEtaPhi;         //!
  TH1F                       *fHistEnergyMedian;        //!
  TH1F                       *fHistEnergyMedianExLP;    //!
  TH1F                       *fHistEnergySum;           //!
  TH1F                       *fHistEnergySumExLP;       //!

  TH1F                       *fHistTracks;              //!
  TH1F                       *fHistTracklets;           //!
  TH1F                       *fHistV0MultSum;           //!

  TH2F                       *fHistEnergyMedianEst[3];        //!
  TH2F                       *fHistEnergyMedianExLPEst[3];    //!
  TH2F                       *fHistEnergySumEst[3];           //!
  TH2F                       *fHistEnergySumExLPEst[3];       //!
  TH2F                       *fHistEnergySumAvgEst[3];        //!
  TH2F                       *fHistEnergySumAvgExLPEst[3];    //!

  TH2F                       *fHistTracksTracklets;           //!
  TH2F                       *fHistTracksV0MultSum;           //!

  AliAnalysisTaskEmcalHighMultTrigger(const AliAnalysisTaskEmcalHighMultTrigger&);            // not implemented
  AliAnalysisTaskEmcalHighMultTrigger &operator=(const AliAnalysisTaskEmcalHighMultTrigger&); // not implemented

  ClassDef(AliAnalysisTaskEmcalHighMultTrigger, 1) // high multiplicity pp trigger analysis task
};
#endif
