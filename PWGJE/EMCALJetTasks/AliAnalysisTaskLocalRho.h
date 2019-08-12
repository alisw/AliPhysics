/**************************************************************************************
 * Copyright (C) 2013, Copyright Holders of the ALICE Collaboration                   *
 * All rights reserved.                                                               *
 *                                                                                    *
 * Redistribution and use in source and binary forms, with or without                 *
 * modification, are permitted provided that the following conditions are met:        *
 *     * Redistributions of source code must retain the above copyright               *
 *       notice, this list of conditions and the following disclaimer.                *
 *     * Redistributions in binary form must reproduce the above copyright            *
 *       notice, this list of conditions and the following disclaimer in the          *
 *       documentation and/or other materials provided with the distribution.         *
 *     * Neither the name of the <organization> nor the                               *
 *       names of its contributors may be used to endorse or promote products         *
 *       derived from this software without specific prior written permission.        *
 *                                                                                    *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE             *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY                *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES         *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;       *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND        *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT         *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 **************************************************************************************/
#ifndef ALIANALYSISTASKLOCALRHO_H
#define ALIANALYSISTASKLOCALRHO_H

// uncomment of define externally to enable debug information
//#define ALIANALYSISTASKLOCALRHO_DEBUG_FLAG_0

#include <AliAnalysisTaskEmcalJet.h>
#include <AliEmcalJet.h>
#include <AliVEvent.h>
#include <AliVTrack.h>
#include <AliVCluster.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TRandom3.h>
#include <AliLog.h>
#include <AliJetContainer.h>

class TF1;
class THF1;
class THF2;
class TProfile;
class AliLocalRhoParameter;

/**
 * @class AliAnalysisTaskLocalRho
 * @brief Analysis task to estimate an event's local energy density
 * @ingroup PWGJEBASE
 * @author Redmer Alexander Bertens <rbertens@cern.ch, rbertens@nikhef.nl, r.a.bertens@uu.nl>, Utrecht Univeristy, Utrecht, Netherlands
 * @since July 26, 2013
 * 
 * This task is part of the emcal jet framework and should be run in the emcaljet train
 * The following extensions to an accepted AliVEvent are expected:
 *      - (anti-kt) jets -> necessary if one wants to exclude leading jet contribution to the event plane
 *      - background estimate of rho -> this task estimates modulation, not rho itself
 *      - pico tracks -> a uniform track selection is necessary to estimate the contribution of v_n harmonics
 *        aod's and esd's are handled transparently
 * The task will estimates a phi-dependent background density rho which is added to the event as a 
 * AliLocalRhoParamter object.
 */
class AliAnalysisTaskLocalRho : public AliAnalysisTaskEmcalJet {
 public:
  // enumerators
  enum fitModulationType  { kNoFit, kV2, kV3, kCombined, kFourierSeries, kIntegratedFlow, kQC2, kQC4 }; // fit type
  enum detectorType       { kTPC, kVZEROA, kVZEROC, kVZEROComb};  // detector that was used
  enum qcRecovery         { kFixedRho, kNegativeVn, kTryFit };    // how to deal with negative cn value for qcn value
  enum runModeType        { kLocal, kGrid };                      // run mode type
  // constructors, destructor
  AliAnalysisTaskLocalRho();
  AliAnalysisTaskLocalRho(const char *name, runModeType type);
  virtual                 ~AliAnalysisTaskLocalRho();
  // setting up the task and technical aspects
  void                    ExecOnce();
  Bool_t                  InitializeAnalysis();
  virtual void            UserCreateOutputObjects();
  TH1F*                   BookTH1F(const char* name, const char* x, Int_t bins, Double_t min, Double_t max, Int_t c = -1, Bool_t append = kTRUE);
  TH2F*                   BookTH2F(const char* name, const char* x, const char* y, Int_t binsx, Double_t minx, Double_t maxx, 
				   Int_t binsy, Double_t miny, Double_t maxy, Int_t c = -1, Bool_t append = kTRUE);
  virtual Bool_t          Run();
  /* inline */   Double_t PhaseShift(Double_t x) const {  
    while (x>=TMath::TwoPi())x-=TMath::TwoPi();
    while (x<0.)x+=TMath::TwoPi();
    return x; }
  /* inline */   Double_t PhaseShift(Double_t x, Double_t n) const {
    x = PhaseShift(x);
    if(TMath::Nint(n)==2) while (x>TMath::Pi()) x-=TMath::Pi();
    if(TMath::Nint(n)==3) {
      if(x>2.*TMath::TwoPi()/n) x = TMath::TwoPi() - x;
      if(x>TMath::TwoPi()/n) x = TMath::TwoPi()-(x+TMath::TwoPi()/n);
    }
    return x; }
  /* inline */   Double_t ChiSquarePDF(Int_t ndf, Double_t x) const {
    Double_t n(ndf/2.), denom(TMath::Power(2, n)*TMath::Gamma(n));
    if (denom!=0)  return ((1./denom)*TMath::Power(x, n-1)*TMath::Exp(-x/2.)); 
    return -999; }
  // the cdf of the chisquare distribution is the normalized lower incomplete gamma function
  /* inline */    Double_t ChiSquareCDF(Int_t ndf, Double_t x) const { return TMath::Gamma(ndf/2., x/2.); }
  // setters - setup how to run
  void                    SetCentralityClasses(TArrayI* c)                {fCentralityClasses = c;}
  void                    SetAttachToEvent(Bool_t a)                      {fAttachToEvent = a;}
  void                    SetUseScaledRho(Bool_t s)                       {fUseScaledRho = s;}
  void                    SetFillHistograms(Bool_t b)                     {fFillHistograms = b;}
  // setters - analysis details
  void                    SetNoEventWeightsForQC(Bool_t e)                {fNoEventWeightsForQC = e;}
  void                    SetIntegratedFlow(TH1F* i, TH1F* j)             {fUserSuppliedV2 = i; fUserSuppliedV3 = j; }
  void                    SetOnTheFlyResCorrection(TH1F* r2, TH1F* r3)    {fUserSuppliedR2 = r2; fUserSuppliedR3 = r3; }
  void                    SetModulationFit(TF1* fit);
  void                    SetModulationFitMinMaxP(Float_t m, Float_t n)   {fMinPvalue = m; fMaxPvalue = n; }
  void                    SetModulationFitType(fitModulationType type)    {fFitModulationType = type; }
  void                    SetQCnRecoveryType(qcRecovery type)             {fQCRecovery = type; }
  void                    SetModulationFitOptions(TString opt)            {fFitModulationOptions = opt; }
  void                    SetReferenceDetector(detectorType type)         {fDetectorType = type; }
  void                    SetUsePtWeight(Bool_t w)                        {fUsePtWeight = w; }
  void                    SetUsePtWeightErrorPropagation(Bool_t w)        {fUsePtWeightErrorPropagation = w;}
  void                    SetRunModeType(runModeType type)                {fRunModeType = type; }
  void                    SetForceAbsVnHarmonics(Bool_t f)                {fAbsVnHarmonics = f; }
  void                    SetExcludeLeadingJetsFromFit(Float_t n)         {fExcludeLeadingJetsFromFit = n; }
  void                    SetRebinSwapHistoOnTheFly(Bool_t r)             {fRebinSwapHistoOnTheFly = r; }
  void                    SetSaveThisPercentageOfFits(Float_t p)          {fPercentageOfFits = p; }
  void                    SetUseV0EventPlaneFromHeader(Bool_t h)          {fUseV0EventPlaneFromHeader = h;}
  void                    SetSoftTrackMinMaxPt(Float_t min, Float_t max)  {fSoftTrackMinPt = min; fSoftTrackMaxPt = max;}
  // getters
  TString                 GetLocalRhoName() const                         {return fLocalRhoName; }
  // numerical evaluations
  void                    CalculateEventPlaneVZERO(Double_t vzero[2][2]) const;
  void                    CalculateEventPlaneTPC(Double_t* tpc);
  void                    CalculateEventPlaneCombinedVZERO(Double_t* comb) const;
  Double_t                CalculateQC2(Int_t harm);
  Double_t                CalculateQC4(Int_t harm);
  // helper calculations for the q-cumulant analysis, also used by AliAnalyisTaskJetFlow
  void                    QCnQnk(Int_t n, Int_t k, Double_t &reQ, Double_t &imQ);
  Double_t                QCnS(Int_t i, Int_t j);
  Double_t                QCnM();
  Double_t                QCnM11();
  Double_t                QCnM1111();
  Bool_t                  QCnRecovery(Double_t psi2, Double_t psi3);
  // analysis details
  Bool_t                  CorrectRho(Double_t psi2, Double_t psi3);
  void                    FillEventPlaneHistograms(Double_t psi2, Double_t psi3) const;
  void                    FillAnalysisSummaryHistogram() const;
  // track selection
  /* inline */     Bool_t PassesCuts(AliVTrack* track) const { UInt_t rejectionReason = 0; return GetParticleContainer(0)->AcceptParticle(track, rejectionReason);} // due to framework updates
  /* inline */     Bool_t PassesCuts(AliEmcalJet* jet) { return AcceptJet(jet, 0);}
  virtual void            Terminate(Option_t* option);

 private: 
  Bool_t                  fInitialized;           //!<! is the analysis initialized?
  Bool_t                  fAttachToEvent;         ///< attach local rho to the event
  Bool_t                  fFillHistograms;        ///< fill qa histograms
  Bool_t                  fNoEventWeightsForQC;   ///< don't store event weights for qc analysis
  Bool_t                  fUseScaledRho;          ///< use scaled rho
  TArrayI*                fCentralityClasses;     ///< centrality classes (maximum 10) used for QA
  TH1F*                   fUserSuppliedV2;        ///< histo with integrated v2
  TH1F*                   fUserSuppliedV3;        ///< histo with integrated v3
  TH1F*                   fUserSuppliedR2;        /// correct the extracted v2 with this r
  TH1F*                   fUserSuppliedR3;        ///< correct the extracted v3 with this r
  Int_t                   fNAcceptedTracks;       //!<! number of accepted tracks
  Int_t                   fNAcceptedTracksQCn;    //!<! accepted tracks for QCn
  Int_t                   fInCentralitySelection; //!<! centrality bin, only for QA plots
  fitModulationType       fFitModulationType;     ///< fit modulation type
  qcRecovery              fQCRecovery;            ///< recovery type for e-by-e qc method
  Bool_t                  fUsePtWeight;           ///< use dptdphi instead of dndphi
  Bool_t                  fUsePtWeightErrorPropagation; ///< recalculate the bin error on the dpt dphi histogram
  detectorType            fDetectorType;          ///< type of detector used for modulation fit
  TString                 fFitModulationOptions;  ///< fit options for modulation fit
  runModeType             fRunModeType;           ///< run mode type 
  TF1*                    fFitModulation;         ///< modulation fit for rho
  Float_t                 fMinPvalue;             ///< minimum value of p
  Float_t                 fMaxPvalue;             ///< maximum value of p
  // additional jet cuts (most are inherited)
  Float_t                 fLocalJetMinEta;        ///< local eta cut for jets
  Float_t                 fLocalJetMaxEta;        ///< local eta cut for jets
  Float_t                 fLocalJetMinPhi;        ///< local phi cut for jets
  Float_t                 fLocalJetMaxPhi;        ///< local phi cut for jets
  Float_t                 fSoftTrackMinPt;        ///< min pt for soft tracks
  Float_t                 fSoftTrackMaxPt;        ///< max pt for soft tracks
  // general qa histograms
  TH1F*                   fHistPvalueCDF;         //!<! cdf value of chisquare p
  TH2F*                   fHistRhoStatusCent;     //!<! status of rho vs centrality
  // general settings
  Bool_t                  fAbsVnHarmonics;                ///< force postive local rho
  Float_t                 fExcludeLeadingJetsFromFit;     ///< exclude n leading jets from fit
  Bool_t                  fRebinSwapHistoOnTheFly;        ///< rebin swap histo on the fly
  Float_t                 fPercentageOfFits;              ///< save this percentage of fits
  Bool_t                  fUseV0EventPlaneFromHeader;     ///< use the vzero event plane from the header
  // transient object pointers
  TList*                  fOutputList;            //!<! output list
  TList*                  fOutputListGood;        //!<! output list for local analysis
  TList*                  fOutputListBad;         //!<! output list for local analysis
  TH1F*                   fHistSwap;              //!<! swap histogram
  TH1F*                   fHistAnalysisSummary;   //!<! flags
  TProfile*               fProfV2;                //!<! extracted v2
  TProfile*               fProfV2Cumulant;        //!<! v2 cumulant
  TProfile*               fProfV3;                //!<! extracted v3
  TProfile*               fProfV3Cumulant;        //!<! v3 cumulant
  TH1F*                   fHistPsi2[10];          //!<! psi 2
  TH1F*                   fHistPsi3[10];          //!<! psi 3

  AliAnalysisTaskLocalRho(const AliAnalysisTaskLocalRho&);                  // not implemented
  AliAnalysisTaskLocalRho& operator=(const AliAnalysisTaskLocalRho&);       // not implemented

  ClassDef(AliAnalysisTaskLocalRho, 6);
};
#endif
