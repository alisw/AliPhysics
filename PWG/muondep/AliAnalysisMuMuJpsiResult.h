#ifndef ALIANALYSISMUMUJPSIRESULT_H
#define ALIANALYSISMUMUJPSIRESULT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/**

@ingroup pwg_muondep_mumu

@class AliAnalysisMuMuJpsiResult

@brief Class to hold results about J/psi

The class is used to hold results like number of of J/psi (before and after Acc x Eff correction),
Acc x Eff correction, Yield, RAB, etc...

A note on "naming conventions" :

-  FitFunctionXY : denotes a function with a prototype double func(double*,double*)
 which can be used in a fit. X = Background, Signal or Total for Background+Signal
Y is the functio name

The FitTypeKey() is the one decoded to select/configure the fit.
This string should be written as : FitType : <key>=<value>:<key>=<value>:<key>=<value>:<key>=<value>...
FitType: func=PSIPSIPRIMECB2VWG:rebin=1:histoType=minv:alJPsi=0.974323:nlJPsi=7.36371:auJPsi=1.84248:nuJPsi=16.0656:range=1.7;4.8

@author Laurent Aphecetche (Subatech)
@author  Benjamin Audurier (Subatech)
*/

#include "TNamed.h"
#include <TString.h>
#include "AliAnalysisMuMuResult.h"
#include "AliAnalysisMuMuBinning.h"

class TH1;
class THashList;
class TF1;
class TMap;
class TFitResultPtr;

class AliAnalysisMuMuJpsiResult : public AliAnalysisMuMuResult
{

public:

  AliAnalysisMuMuJpsiResult(TRootIOCtor* io);

//  AliAnalysisMuMuJpsiResult(const TH1& hminv);

  AliAnalysisMuMuJpsiResult(const char* particle,
                            const TH1& hminv,
                            const char* fitType);

  AliAnalysisMuMuJpsiResult(const char* particle,
                            const TH1& hminv,
                            const char* triggerClass,
                            const char* eventSelection,
                            const char* pairSelection,
                            const char* centSelection,
                            const AliAnalysisMuMuBinning::Range& bin);

  AliAnalysisMuMuJpsiResult(const AliAnalysisMuMuJpsiResult& rhs);
  AliAnalysisMuMuJpsiResult& operator=(const AliAnalysisMuMuJpsiResult& rhs);

  virtual ~AliAnalysisMuMuJpsiResult();

  virtual TObject* Clone(const char* newname = "") const;

  Bool_t Correct(const AliAnalysisMuMuJpsiResult& other, const char* particle, const char* subResultName="");

  TH1* Histo() const { return fHisto; }

  Int_t NofTriggers() const;

  void SetNofTriggers(Int_t n);

  void Print(Option_t* opt="") const;

  Bool_t AddFit(const char* fitType);

  /** All the fit functions should have a prototype starting like :

   AliAnalysisMuMuJpsiResult* FitXXX();

   If extra parameters to the specific FitXXX function are needed, they should be given
   using the SetValue(...) method, and retrieved from the FitXXX method using
   the GetValue(...) method.
   */


  void FitPSICOUNT();
  void FitPSICB2();
  void FitPSINA60NEW();

  // void FitPSICB2VWG();
  void FitPSIPSIPRIMECB2VWG();
  void FitPSIPSIPRIMECB2VWG2();
  void FitPSIPSIPRIMECB2POL1POL2();
  void FitPSIPSIPRIMECB2POL2POL3();
  void FitPSIPSIPRIMECB2POL2POL3V2();
  void FitPSIPSIPRIMECB2POL2EXP();
  void FitPSIPSIPRIMENA60NEWVWG();
  void FitPSIPSIPRIMENA60NEWVWG2();
  void FitPSIPSIPRIMENA60NEWPOL1POL2();
  void FitPSIPSIPRIMENA60NEWPOL2POL3();
  void FitPSIPSIPRIMENA60NEWPOL2EXP();

  void FitPSIPSIPRIMECB2POL4EXP();
  void FitPSIPSIPRIMENA60NEWPOL4EXP();
  void FitPSIPSIPRIMECB2VWGINDEPTAILS();

  //** All the mean pt fit methods MUST contain the corresponding name of the inv mass spectra method
  void FitMPTPSIPSIPRIMECB2VWG_BKGMPTPOL2();
  void FitMPTPSIPSIPRIMECB2POL1POL2_BKGMPTPOL2();
  void FitMPTPSIPSIPRIMECB2VWG_BKGMPTPOL2EXP();
  void FitMPTPSIPSIPRIMECB2POL1POL2_BKGMPTPOL2EXP();
  void FitMPTPSIPSIPRIMECB2POL2EXP_BKGMPTPOL2();
  void FitMPTPSIPSIPRIMECB2POL2EXP_BKGMPTPOL2EXP();

  void FitMPTPSIPSIPRIMECB2VWG_BKGMPTLIN();
  void FitMPTPSIPSIPRIMECB2VWG_BKGMPTPOL3();
  void FitMPTPSIPSIPRIMECB2VWG_BKGMPTPOL4();
  void FitMPTPSIPSIPRIMECB2VWGINDEPTAILS_BKGMPTPOL2();

  void FitMPTPSIPSIPRIMENA60NEWVWG_BKGMPTPOL2();
  void FitMPTPSIPSIPRIMENA60NEWVWG_BKGMPTPOL2EXP();
  void FitMPTPSIPSIPRIMENA60NEWPOL1POL2_BKGMPTPOL2();
  void FitMPTPSIPSIPRIMENA60NEWPOL1POL2_BKGMPTPOL2EXP();
  void FitMPTPSIPSIPRIMENA60NEWPOL2EXP_BKGMPTPOL2();
  void FitMPTPSIPSIPRIMENA60NEWPOL2EXP_BKGMPTPOL2EXP();

  void FitMPTPSI_HFUNCTION();

  Int_t NofRuns() const;

  void SetNofRuns(int n);

  const AliAnalysisMuMuBinning::Range& Bin() const;

  void SetBin(const AliAnalysisMuMuBinning::Range& bin);

  void SetNofInputParticles(const char* particle, int n);

  void SetNofInputParticles(const TH1& hminv);

  void SetHisto(const TH1& h);

  Long64_t Merge(TCollection* list);

  static Double_t CountParticle(const TH1& hminv, const char* particle, Double_t sigma=-1);

  virtual AliAnalysisMuMuJpsiResult* Mother() const { return static_cast<AliAnalysisMuMuJpsiResult*>(AliAnalysisMuMuResult::Mother()); }

  void PrintValue(const char* key, const char* opt, Double_t value, Double_t errorStat, Double_t rms=0.0) const;

  void ProcessMinvFit(TFitResultPtr& fitResult, TF1* fitTotal, TF1* bckInit, const char* fitOption, Int_t iParKPsip, Int_t iLastParBkg);

  void ProcessBkgFit(TFitResultPtr& fitResultInit, TF1* bckInit, const char* bkgFuncName, const char* fitOption);

  TString FitFunctionName() const { return fFitFunction; }

  TString GetFitFunctionMethodName() const;

  void Draw(Option_t* opt="");

  const char* GetParticle() { return fParticle; }

private:

  enum EIndex
  {
    kValue=0,
    kErrorStat=1
  };

  void DecodeFitType(const char* fitType);

  void PrintParticle(const char* particle, const char* opt) const;

  Double_t FitFunctionBackgroundLin(Double_t *x, Double_t *par);

  Double_t FitFunctionBackgroundPol1Pol2(Double_t *x, Double_t *par);

  Double_t FitFunctionBackgroundPol2Pol3(Double_t *x, Double_t *par);

  Double_t FitFunctionBackgroundPol2Pol3V2(Double_t *x, Double_t *par);

  Double_t FitFunctionBackgroundPol2Exp(Double_t* x, Double_t* par);

  Double_t FitFunctionBackgroundPol4Exp(Double_t *x, Double_t *par);

  Double_t FitFunctionBackgroundPol2(Double_t *x, Double_t *par);

  Double_t FitFunctionBackgroundPol3(Double_t *x, Double_t *par);

  Double_t FitFunctionBackgroundPol4(Double_t *x, Double_t *par);

  Double_t FitFunctionBackgroundVWG(Double_t* x, Double_t* par);

  Double_t FitFunctionBackgroundVWG2(Double_t* x, Double_t* par);

  Double_t FitFunctionSignalCrystalBallExtended(Double_t *x,Double_t *par);

  Double_t FitFunctionNA60New(Double_t *x,Double_t *par);

  // Double_t FitFunctionTotalOneCB2VWG(Double_t *x,Double_t *par); // Correction here

  Double_t FitFunctionTotalTwoNA60NewVWG(Double_t *x, Double_t *par);

  Double_t FitFunctionTotalTwoNA60NewVWG2(Double_t *x, Double_t *par);

  Double_t FitFunctionTotalTwoNA60NewPol1Pol2(Double_t *x, Double_t *par);

  Double_t FitFunctionTotalTwoNA60NewPol2Pol3(Double_t *x, Double_t *par);

  Double_t FitFunctionTotalTwoNA60NewPol2Exp(Double_t *x, Double_t *par);

  Double_t FitFunctionTotalTwoNA60NewPol4Exp(Double_t *x, Double_t *par);

  Double_t FitFunctionTotalTwoCB2Pol2Exp(Double_t *x, Double_t *par);

  Double_t FitFunctionTotalTwoCB2Pol1Pol2(Double_t *x, Double_t *par);

  Double_t FitFunctionTotalTwoCB2Pol2Pol3(Double_t *x, Double_t *par);

  Double_t FitFunctionTotalTwoCB2Pol2Pol3V2(Double_t *x, Double_t *par);

  Double_t FitFunctionTotalTwoCB2Pol4Exp(Double_t *x, Double_t *par);

  Double_t FitFunctionTotalTwoCB2VWG(Double_t *x, Double_t *par);

  Double_t FitFunctionTotalTwoCB2VWG2(Double_t *x, Double_t *par);

  Double_t FitFunctionTotalTwoCB2Lin(Double_t *x, Double_t *par);

  Double_t FitFunctionTotalTwoCB2VWGINDEPTAILS(Double_t *x, Double_t *par);

  Double_t hFunction(Double_t*x, Double_t* par);

  Double_t alphaCB2VWG(Double_t*x, Double_t* par);

  Double_t alphaCB2POL1POL2(Double_t*x, Double_t* par);

  Double_t alphaCB2POL2EXP(Double_t*x, Double_t* par);

  Double_t alphaNA60NEWVWG(Double_t*x, Double_t* par);

  Double_t alphaNA60NEWPOL1POL2(Double_t*x, Double_t* par);

  Double_t alphaNA60NEWPOL2EXP(Double_t*x, Double_t* par);

  Double_t FitFunctionMeanPtSCB2Lin(Double_t* x, Double_t* par);

  Double_t FitFunctionMeanPtS2CB2Lin(Double_t *x, Double_t *par);

  Double_t FitFunctionMeanPtSCB2VWGPOL2(Double_t* x, Double_t* par);



  Double_t FitFunctionMeanPtS2CB2VWGPOL2(Double_t *x, Double_t *par);

  Double_t FitFunctionMeanPtS2CB2POL1POL2POL2(Double_t *x, Double_t *par);

  Double_t FitFunctionMeanPtS2CB2VWGPOL2EXP(Double_t *x, Double_t *par);

  Double_t FitFunctionMeanPtS2CB2POL1POL2POL2EXP(Double_t *x, Double_t *par);

  Double_t FitFunctionMeanPtS2CB2POL2EXPPOL2(Double_t *x, Double_t *par);

  Double_t FitFunctionMeanPtS2CB2POL2EXPPOL2EXP(Double_t *x, Double_t *par);

  Double_t FitFunctionMeanPtS2NA60NEWVWGPOL2(Double_t *x, Double_t *par);

  Double_t FitFunctionMeanPtS2NA60NEWPOL1POL2POL2(Double_t *x, Double_t *par);

  Double_t FitFunctionMeanPtS2NA60NEWVWGPOL2EXP(Double_t *x, Double_t *par);

  Double_t FitFunctionMeanPtS2NA60NEWPOL1POL2POL2EXP(Double_t *x, Double_t *par);

  Double_t FitFunctionMeanPtS2NA60NEWPOL2EXPPOL2(Double_t *x, Double_t *par);

  Double_t FitFunctionMeanPtS2NA60NEWPOL2EXPPOL2EXP(Double_t *x, Double_t *par);

  Double_t FitFunctionMeanPtHFunction(Double_t *x, Double_t *par);


  Double_t FitFunctionMeanPtS2CB2VWGPOL3(Double_t *x, Double_t *par);

  Double_t FitFunctionMeanPtS2CB2VWGPOL2INDEPTAILS(Double_t *x, Double_t *par);

  Double_t FitFunctionMeanPtS2CB2VWGPOL4(Double_t *x, Double_t *par);

  void SetFitRejectRange(Double_t a=TMath::Limits<Double_t>::Max(),
                         Double_t b=TMath::Limits<Double_t>::Max());

  void AttachFunctionsToHisto(TF1* signal, TF1* bck, TF1* total, Double_t xmin, Double_t xmax); //Might be unnecesary since we have the one with 2 signals

  void AttachFunctionsToHisto(TF1* signal1, TF1* signal2, TF1* bck, TF1* total,Double_t xmin, Double_t xmax);

  void SetParameter(TF1* func, Int_t npar, Double_t fix, Double_t initialValue,
                    Double_t min, Double_t max) const;

  Bool_t CheckRoots(TFitResultPtr &fitResult, TF1* fitFunction, Int_t deg, Double_t a, Double_t b, Double_t c, Double_t d, const char* fitOption );

  Bool_t WrongParameter(TF1* fitFunction, Int_t npar, Double_t fixValueIfWrong);

  Bool_t StrongCorrelation(TFitResultPtr& fitResult, TF1* fitFunction, Int_t npar1, Int_t npar2, Double_t fixValueIfWrong);


private:
  Int_t fNofRuns; // number of runs used to get this result
  Int_t fNofTriggers; // number of trigger analyzed
  TH1* fHisto; // invariant mass spectrum
  AliAnalysisMuMuBinning::Range fBin; // bin range

  TString fTriggerClass; // trigger class for this result
  TString fEventSelection; // event selection for this result
  TString fPairSelection; // pair selection for this result
  TString fCentralitySelection; // centrality selection for this result

  TString fFitFunction; // fit function used
  Double_t fFitRejectRangeLow; // fit range to reject
  Double_t fFitRejectRangeHigh; // fit range to reject
  Bool_t fRejectFitPoints; // whether or not some fit points should be rejected

  TString fParticle;
  TString fMinvRS; // minv spectra range and sigmaPsiP factor for the mpt fits

  ClassDef(AliAnalysisMuMuJpsiResult,8) // a class to hold invariant mass analysis results (counts, yields, AccxEff, R_AB, etc...)
};

#endif
