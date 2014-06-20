#ifndef AliAnalysisTaskLocalRhoDev_H
#define AliAnalysisTaskLocalRhoDev_H

// $Id$

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
#include <AliParticleContainer.h>
#include <TF1.h>
#include <TH1.h>

class THF1;
class THF2;
class TProfile;
class AliLocalRhoParameter;
class TArrayI;

class AliAnalysisTaskLocalRhoDev : public AliAnalysisTaskEmcalJet {
 public:
  // enumerators
  enum fitModulationType  { kNoFit, kV2, kV3, kCombined, kFourierSeries, kIntegratedFlow, kQC2, kQC4 }; // fit type
  enum fitGoodnessTest    { kChi2ROOT, kChi2Poisson, kLinearFit };
  enum detectorType       { kTPC, kVZEROA, kVZEROC, kVZEROComb};  // detector that was used
  enum qcRecovery         { kFixedRho, kNegativeVn, kTryFit };    // how to deal with negative cn value for qcn value
  enum runModeType        { kLocal, kGrid };                      // run mode type
  // constructors, destructor
  AliAnalysisTaskLocalRhoDev();
  AliAnalysisTaskLocalRhoDev(const char *name, runModeType type);
  virtual                 ~AliAnalysisTaskLocalRhoDev();
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

  /* inline */    Double_t ChiSquarePDF(Int_t ndf, Double_t x) const {
    Double_t n(ndf/2.), denom(TMath::Power(2, n)*TMath::Gamma(n));
    if (denom!=0)  return ((1./denom)*TMath::Power(x, n-1)*TMath::Exp(-x/2.)); 
    return -999; }
  //note that the cdf of the chisquare distribution is the normalized lower incomplete gamma function
  /*inline */    Double_t ChiSquareCDF(Int_t ndf, Double_t x) const { return TMath::Gamma(ndf/2., x/2.); }
  /*inline */    Double_t ChiSquare(TH1& histo, TF1* func) const {
    // evaluate the chi2 using a poissonian error estimate on bins
    Double_t chi2(0.);
    for(Int_t i(0); i < histo.GetXaxis()->GetNbins(); i++) {
      if(histo.GetBinContent(i+1) <= 0.) continue;
      chi2 += TMath::Power((histo.GetBinContent(i+1)-func->Eval(histo.GetXaxis()->GetBinCenter(1+i))), 2)/histo.GetBinContent(i+1);
    }
    return chi2;
  }
  /*inline*/ Double_t KolmogorovTest(TH1F& histo, TF1* func, Bool_t toy = kTRUE) const {
    // return the probability from a Kolmogorov test
    TH1F test(histo);       // stack copy of test statistic
    for(Int_t i(0); i < test.GetXaxis()->GetNbins(); i++) test.SetBinContent(i+1, func->Eval(test.GetXaxis()->GetBinCenter(1+i)));
    if (toy) return histo.TH1::KolmogorovTest((&test), "X");
    else return histo.TH1::KolmogorovTest((&test));
  }
 
  // setters - setup how to run
  void                    SetDebugMode(Int_t d)                           {fDebug = d;}
  void                    SetCentralityClasses(TArrayD* c)                {fCentralityClasses = c;}
  void                    SetAttachToEvent(Bool_t a)                      {fAttachToEvent = a;}
  void                    SetUseScaledRho(Bool_t s)                       {fUseScaledRho = s;}
  void                    SetFillHistograms(Bool_t b)                     {fFillHistograms = b;}
  // setters - analysis details
  void                    SetNoEventWeightsForQC(Bool_t e)                {fNoEventWeightsForQC = e;}
  void                    SetIntegratedFlow(TH1F* i, TH1F* j)             {fUserSuppliedV2 = i; fUserSuppliedV3 = j; }
  void                    SetOnTheFlyResCorrection(TH1F* r2, TH1F* r3)    {fUserSuppliedR2 = r2; fUserSuppliedR3 = r3; }
  void                    SetModulationFit(TF1* fit);
  void                    SetModulationFitMinMaxP(Float_t m, Float_t n)   {fMinPvalue = m; fMaxPvalue = n; }
  void                    SetExpectedRuns(TArrayI* r)                     {fExpectedRuns = r; }
  void                    SetExpectedSemiGoodRuns(TArrayI* r)             {fExpectedSemiGoodRuns = r;}
  void                    SetNameRhoSmall(TString s)                      {fNameSmallRho = s;}
  void                    SetModulationFitType(fitModulationType type)    {fFitModulationType = type; }
  void                    SetControlFunction(TF1* func)                   {fFitControl = func; }
  void                    SetFitGoodnessType(fitGoodnessTest test)        {fFitGoodnessTest = test; }
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
  /* inline */     Bool_t PassesCuts(AliVTrack* track) const { return AcceptTrack(track, 0);}
  /* inline */     Bool_t PassesCuts(AliEmcalJet* jet) { return AcceptJet(jet, 0);}
  /* inline*/      Bool_t PassesSimpleCuts(AliEmcalJet* jet) {
  Float_t minPhi(GetJetContainer()->GetJetPhiMin()), maxPhi(GetJetContainer()->GetJetPhiMax());
  return (jet && jet->Pt() > 1 && jet->Eta() < .9-GetJetRadius() && jet->Eta() > -.9+GetJetRadius() && jet->Phi() > minPhi && jet->Phi() < maxPhi && jet->Area() > .557*GetJetRadius()*GetJetRadius()*TMath::Pi());
  }

  /* inline */     AliEmcalJet* GetLeadingJet() {
      Int_t iJets(fJets->GetEntriesFast());
      Double_t pt(0);
      AliEmcalJet* leadingJet(0x0);
      for(Int_t i(0); i < iJets; i++) {
          AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJets->At(i));
          if(!PassesSimpleCuts(jet)) continue;
          if(jet->Pt() > pt) {
              leadingJet = jet;
              pt = leadingJet->Pt();
          }
      }
      return leadingJet;
  }
  Bool_t                  PassesCuts(AliVEvent* event);
  virtual void            Terminate(Option_t* option);

 private: 
  Int_t                   fDebug;                 // debug level (0 none, 1 fcn calls, 2 verbose)
  Bool_t                  fInitialized;           //! is the analysis initialized?
  Bool_t                  fAttachToEvent;         // attach local rho to the event
  Bool_t                  fFillHistograms;        // fill qa histograms
  Bool_t                  fNoEventWeightsForQC;   // don't store event weights for qc analysis
  Bool_t                  fUseScaledRho;          // use scaled rho
  TArrayD*                fCentralityClasses;     // centrality classes (maximum 10) used for QA
  TH1F*                   fUserSuppliedV2;        // histo with integrated v2
  TH1F*                   fUserSuppliedV3;        // histo with integrated v3
  TH1F*                   fUserSuppliedR2;        // correct the extracted v2 with this r
  TH1F*                   fUserSuppliedR3;        // correct the extracted v3 with this r
  Int_t                   fNAcceptedTracks;       //! number of accepted tracks
  Int_t                   fNAcceptedTracksQCn;    //! accepted tracks for QCn
  Int_t                   fInCentralitySelection; //! centrality bin, only for QA plots
  fitModulationType       fFitModulationType;     // fit modulation type
  fitGoodnessTest         fFitGoodnessTest;       // goodness of fit test
  qcRecovery              fQCRecovery;            // recovery type for e-by-e qc method
  Bool_t                  fUsePtWeight;           // use dptdphi instead of dndphi
  Bool_t                  fUsePtWeightErrorPropagation; // recalculate the bin error on the dpt dphi histogram
  detectorType            fDetectorType;          // type of detector used for modulation fit
  TString                 fFitModulationOptions;  // fit options for modulation fit
  runModeType             fRunModeType;           // run mode type 
  TF1*                    fFitModulation;         // modulation fit for rho
  TF1*                    fFitControl;            // control function
  Float_t                 fMinPvalue;             // minimum value of p
  Float_t                 fMaxPvalue;             // maximum value of p
  TArrayI*                fExpectedRuns;          // list of known run numbers with default cuts
  TArrayI*                fExpectedSemiGoodRuns;  // list of runs that are marked as semi-good in the rct
  Int_t                   fRunNumber;             //! current run number
  AliRhoParameter*        fCachedRho;             //! cached rho object
  TString                 fNameSmallRho;          // name of rho object for semi-good tpc runs
  // additional jet cuts (most are inherited)
  Float_t                 fLocalJetMinEta;        // local eta cut for jets
  Float_t                 fLocalJetMaxEta;        // local eta cut for jets
  Float_t                 fLocalJetMinPhi;        // local phi cut for jets
  Float_t                 fLocalJetMaxPhi;        // local phi cut for jets
  Float_t                 fSoftTrackMinPt;        // min pt for soft tracks
  Float_t                 fSoftTrackMaxPt;        // max pt for soft tracks
  Double_t                fSemiGoodJetMinPhi;     // min phi for semi good tpc runs
  Double_t                fSemiGoodJetMaxPhi;     // max phi for semi good tpc runs
  Double_t                fSemiGoodTrackMinPhi;   // min phi for semi good tpc runs
  Double_t                fSemiGoodTrackMaxPhi;   // max phi for semi good tpc runs

  // general qa histograms
  TH1F*                   fHistPvalueCDFROOT;     //! pdf value of chisquare p
  TH2F*                   fHistPvalueCDFROOTCent; //! p value versus centrlaity from root
  TH2F*                   fHistChi2ROOTCent;      //! reduced chi2 from ROOT, centrality correlation
  TH2F*                   fHistPChi2Root;         //! correlation p value and reduced chi2
  TH1F*                   fHistPvalueCDF;         //! cdf value of chisquare p
  TH2F*                   fHistPvalueCDFCent;     //! p value vs centrality
  TH2F*                   fHistChi2Cent;          //! reduced chi2, centrlaity correlation
  TH2F*                   fHistPChi2;             //! correlation p value and reduced chi2
  TH2F*                   fHistRhoStatusCent;     //! status of rho as function of centrality

  // general settings
  Bool_t                  fAbsVnHarmonics;        // force postive local rho
  Float_t                 fExcludeLeadingJetsFromFit;    // exclude n leading jets from fit
  Bool_t                  fRebinSwapHistoOnTheFly;       // rebin swap histo on the fly
  Float_t                 fPercentageOfFits;      // save this percentage of fits
  Bool_t                  fUseV0EventPlaneFromHeader;    // use the vzero event plane from the header
  // transient object pointers
  TList*                  fOutputList;            //! output list
  TList*                  fOutputListGood;        //! output list for local analysis
  TList*                  fOutputListBad;         //! output list for local analysis
  TH1F*                   fHistSwap;              //! swap histogram
  TH1F*                   fHistAnalysisSummary;   //! flags
  TProfile*               fProfV2;                //! extracted v2
  TProfile*               fProfV2Cumulant;        //! v2 cumulant
  TProfile*               fProfV3;                //! extracted v3
  TProfile*               fProfV3Cumulant;        //! v3 cumulant
  TH1F*                   fHistPsi2[10];          //! psi 2
  TH1F*                   fHistPsi3[10];          //! psi 3

  AliAnalysisTaskLocalRhoDev(const AliAnalysisTaskLocalRhoDev&);                  // not implemented
  AliAnalysisTaskLocalRhoDev& operator=(const AliAnalysisTaskLocalRhoDev&);       // not implemented

  ClassDef(AliAnalysisTaskLocalRhoDev, 6);
};
#endif
