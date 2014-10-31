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
  void                        SetNExcludeLeadingPatches(Int_t n)  { fNExLP             = n; }
  void                        SetTruncateThreshold(Double_t t)    { fTruncateThreshold = t; }

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
  Double_t                    fTruncatedMean;          //truncated mean
  Double_t                    fTruncateThreshold;      //threshold used for truncating

  //Histograms
  TH3F                       *fHistPatchEtaPhiE;              //! patch eta vs phi (center of patch) vs energy
  TH1F                       *fHistEnergyMedian;              //! median energy in EMCal
  TH1F                       *fHistEnergyMedianExLP;          //! median energy in EMCal exclucing N leading patches
  TH1F                       *fHistEnergySum;                 //! total energy in EMCal
  TH1F                       *fHistEnergySumExLP;             //! total energy in EMCal exclucing N leading patches
  TH1F                       *fHistTruncatedMean;             //! truncated mean in EMCal

  TH1F                       *fHistTracks;                    //! N hybrid tracks
  TH1F                       *fHistTracklets;                 //! Ntracklets
  TH1F                       *fHistV0MultSum;                 //! V0A+V0C multiplicity

  TH2F                       *fHistEnergyMedianEst[3];        //! median energy in EMCal vs mult estimator
  TH2F                       *fHistEnergyMedianExLPEst[3];    //! median energy in EMCal excluding N leading patches vs mult estimator
  TH2F                       *fHistEnergySumEst[3];           //! total energy in EMCal vs mult estimator
  TH2F                       *fHistEnergySumExLPEst[3];       //! total energy in EMCal excluding N leading patches vs mult estimator
  TH2F                       *fHistEnergySumAvgEst[3];        //! avg energy in EMCal vs mult estimator
  TH2F                       *fHistEnergySumAvgExLPEst[3];    //! avg energy in EMCal excluding N leading patches vs mult estimator
  TH2F                       *fHistTruncatedMeanEst[3];       //! truncated mean in EMCal vs mult estimator

  TH2F                       *fHistTracksTracklets;           //! Ntracks vs Ntracklets
  TH2F                       *fHistTracksV0MultSum;           //! Ntracks vs V0A+V0C
  TH3F                       *fHistSPDTrkClsSum;              //! correlation between SPD clusters and tracklets and total energy in EMCal
  TH3F                       *fHistSPDTrkClsSumExLP;          //! correlation between SPD clusters and tracklets and total energy in EMCal
  TH3F                       *fHistSPDTrkClsMedian;           //! correlation between SPD clusters and tracklets and median energy in EMCal
  TH3F                       *fHistSPDTrkClsMedianExLP;       //! correlation between SPD clusters and tracklets and median energy in EMCal
  TH3F                       *fHistSPDTrkClsTruncMean;        //! correlation between SPD clusters and tracklets and truncated mean in EMCal

  AliAnalysisTaskEmcalHighMultTrigger(const AliAnalysisTaskEmcalHighMultTrigger&);            // not implemented
  AliAnalysisTaskEmcalHighMultTrigger &operator=(const AliAnalysisTaskEmcalHighMultTrigger&); // not implemented

  ClassDef(AliAnalysisTaskEmcalHighMultTrigger, 4) // high multiplicity pp trigger analysis task
};
#endif
