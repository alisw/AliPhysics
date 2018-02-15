#ifndef ALIDIELECTRONSIGNALBASE_H
#define ALIDIELECTRONSIGNALBASE_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           # 
//#         Class AliDielectronSignalBase                       #
//#         Manage Cuts on the legs of the pair               #
//#                                                           #
//#  Authors:                                                 #
//#   Anton     Andronic, GSI / A.Andronic@gsi.de             #
//#   Ionut C.  Arsene,   GSI / I.C.Arsene@gsi.de             #
//#   Julian    Book,     Uni Ffm / Julian.Book@cern.ch       #
//#   Frederick Kramer,   Uni Ffm, / Frederick.Kramer@cern.ch #
//#   Magnus    Mager,    CERN / Magnus.Mager@cern.ch         #
//#   WooJin J. Park,     GSI / W.J.Park@gsi.de               #
//#   Jens      Wiechula, Uni HD / Jens.Wiechula@cern.ch      #
//#                                                           #
//#############################################################


#include <TNamed.h>
#include <TVectorT.h>
#include <TMath.h>
#include <TH1F.h>
#include <TF1.h>
#include <TExMap.h>


class TObjArray;
class TPaveText;

class AliDielectronSignalBase : public TNamed {
public:
  enum EBackgroundMethod {
    kFittedMC = 0,
    kFitted,
    kLikeSign,
    kLikeSignArithm,
    kLikeSignRcorr,
    kLikeSignArithmRcorr,
    kLikeSignFit,
    kEventMixing,
    kEventMixingFit,
    kRotation,
    kCombinatorialPlusFit
  };

  enum ESignalExtractionMethod {
    kBinCounting = 0,
    kMCScaledMax,
    kMCScaledInt,
    kMCFitted,
    kCrystalBall,
    kGaus
  };

  AliDielectronSignalBase();
  AliDielectronSignalBase(const char*name, const char* title);
  AliDielectronSignalBase(const char*name, const char* title, bool enummaps);
  
  virtual ~AliDielectronSignalBase();
  
  TExMap MapBackgroundMethod;
  TExMap MapSignalExtractionMethod;

  void SetMCSignalShape(TH1F* hist) { fHistSimPM=hist; }
  void SetParticleOfInterest(Int_t pdgcode) { fPOIpdg=pdgcode; }
  void SetIntegralRange(Double_t min, Double_t max) {fIntMin=min;fIntMax=max;}
  void SetFitRange(Double_t min, Double_t max) {fFitMin=min; fFitMax=max;}
  void SetRebin(Int_t factor) {fRebin = factor;}
  void SetMethod(EBackgroundMethod method) {fMethod = method;}
  Int_t GetMethod() const { return (Int_t)fMethod; }
  void SetExtractionMethod(ESignalExtractionMethod method) {fPeakMethod = method;}
  Int_t GetExtractionMethod() const { return (Int_t)fPeakMethod; }

  const TVectorD& GetValues() const {return fValues;}
  const TVectorD& GetErrors() const {return fErrors;}

  Double_t GetIntegralMin()          const { return fIntMin; }
  Double_t GetIntegralMax()          const { return fIntMax; }
  Double_t GetSignal()               const { return fValues(0);}
  Double_t GetSignalError()          const { return fErrors(0);}
  Double_t GetBackground()           const { return fValues(1);}
  Double_t GetBackgroundError()      const { return fErrors(1);}
  Double_t GetSignificance()         const { return fValues(2);}
  Double_t GetSignificanceError()    const { return fErrors(2);}
  Double_t GetSB()                   const { return fValues(3);}
  Double_t GetSBError()              const { return fErrors(3);}
  Double_t GetMass()                 const { return fValues(4);}
  Double_t GetMassError()            const { return fErrors(4);}
  Double_t GetMassWidth()            const { return fValues(5);}
  Double_t GetMassWidthError()       const { return fErrors(5);}
  static const char* GetValueName(Int_t i) { return (i>=0&&i<6)?fgkValueNames[i]:""; }

  TH1* GetSignalHistogram()      const {return fHistSignal;}
  TH1* GetBackgroundHistogram()  const {return fHistBackground;}
  TH1* GetUnlikeSignHistogram()  const {return fHistDataPM;}
  TH1* GetRfactorHistogram()     const {return fHistRfactor;}
  TObject* GetPeakShape()        const {return fPeakShapeObj;}
  
  void SetScaleRawToBackground(Double_t intMin, Double_t intMax) { fScaleMin=intMin; fScaleMax=intMax; }
  void SetScaleRawToBackground(Double_t intMin, Double_t intMax, Double_t intMin2, Double_t intMax2) { fScaleMin=intMin; fScaleMax=intMax; fScaleMin2=intMin2; fScaleMax2=intMax2; }
  Double_t GetScaleMin()        const { return fScaleMin;    }
  Double_t GetScaleMax()        const { return fScaleMax;    }
  Double_t GetScaleMin2()       const { return fScaleMin2;   }
  Double_t GetScaleMax2()       const { return fScaleMax2;   }
  Double_t GetScaleFactor()     const { return fScaleFactor; }
  Int_t GetParticleOfInterest() const { return fPOIpdg;      }

  void SetMixingCorrection(Bool_t mixcorr=kTRUE) { fMixingCorr=mixcorr; }
  
  static Double_t ScaleHistograms(TH1* histRaw, TH1* histBackground, Double_t intMin, Double_t intMax);
  static Double_t ScaleHistograms(TH1* histRaw, TH1* histBackground, Double_t intMin, Double_t intMax, Double_t intMin2, Double_t intMax2);
  
  virtual void Print(Option_t *option="") const;

  /**
  This function needs to be implemented by the signal extraction classes.
  Here all the work should be done.
  
  The signal extraction is done on the mass spectra.
  The TObjArray should contain the Inv. Mass spectra of the 10 possible combinations
     for single and mixed events defined in AliDielectron.cxx
  */
  virtual void Process(TObjArray * const /*arrhist*/) = 0;
  TObject* DescribePeakShape(ESignalExtractionMethod method, Bool_t replaceValErr=kFALSE,  TH1F *mcShape=0x0);

  TObject* fPeakShapeObj;       // histogram or function used to describe the extracted signal
  TH1F* fHistSimPM;          // simulated peak shape

protected: 

  TH1 *fHistSignal;                  // histogram of pure signal
  TH1 *fHistBackground;              // histogram of background (fitted=0, like-sign=1, event mixing=2)
  TH1 *fHistDataPM;                  // histogram of selected +- pair candidates
  TH1 *fHistDataPP;                  // histogram of selected ++ pair candidates
  TH1 *fHistDataMM;                  // histogram of selected -- pair candidates
  TH1 *fHistDataME;                  // histogram of selected +- pair candidates from mixed event
  TH1 *fHistRfactor;                 // histogram of R factors

  TVectorD fValues;                  // values
  TVectorD fErrors;                  // value errors

  Double_t fIntMin;                  // signal extraction range min
  Double_t fIntMax;                  // signal extraction range max
  Double_t fFitMin;                  // fit range lowest inv. mass
  Double_t fFitMax;                  // fit range highest inv. mass

  Int_t fRebin;                      // histogram rebin factor
  EBackgroundMethod fMethod;         // method for background substraction
  Double_t fScaleMin;                // min for scaling of raw and background histogram
  Double_t fScaleMax;                // max for scaling of raw and background histogram
  Double_t fScaleMin2;                // min for scaling of raw and background histogram
  Double_t fScaleMax2;                // max for scaling of raw and background histogram
  Double_t fScaleFactor;             // scale factor of raw to background histogram scaling
  Bool_t fMixingCorr;                // switch for bin by bin correction with R factor

  ESignalExtractionMethod fPeakMethod;  // method for peak description and signal extraction
  Bool_t fProcessed;                 // flag
  Int_t  fPOIpdg;                    // pdg code particle of interest

  void SetSignificanceAndSOB();      // calculate the significance and S/B
  void SetFWHM();                    // calculate the fwhm
  void SetBackgroundEnumMap();       // set enum string relations
  void SetSignalExtractionEnumMap(); // set enum string relations
  static const char* fgkValueNames[6];  //value names
  static const Double_t fgkErrorZero;  //statistical error if zero entries

  TPaveText* DrawStats(Double_t x1=0., Double_t y1=0., Double_t x2=0., Double_t y2=0.);

  AliDielectronSignalBase(const AliDielectronSignalBase &c);
  AliDielectronSignalBase &operator=(const AliDielectronSignalBase &c);

  ClassDef(AliDielectronSignalBase,6)         // base and abstract class for signal extraction
};

inline void AliDielectronSignalBase::SetSignificanceAndSOB()
{
  //
  // Calculate S/B and significance
  //
  // Signal/Background
  fValues(3) = (fValues(1)>0 ? fValues(0)/fValues(1) : 0);
  Float_t epsSig = (fValues(0)>0 ? fErrors(0)/fValues(0) : 0);
  Float_t epsBknd = (fValues(1)>0 ? fErrors(1)/fValues(1) : 0);
  fErrors(3) = fValues(3)*TMath::Sqrt(epsSig*epsSig + epsBknd*epsBknd);
  // Significance
  fValues(2) = ((fValues(0)+fValues(1))>0 ? fValues(0)/TMath::Sqrt(fValues(0)+fValues(1)) : 0);
  Float_t s = (fValues(0)>0?fValues(0):0); Float_t b = fValues(1);
  Float_t se = fErrors(0); Float_t be = fErrors(1);
  // fErrors(2) = ((s+b)>0 ? TMath::Sqrt((s*(s+2*b)*(s+2*b)+b*s*s)/(4*TMath::Power(s+b,3))) : 0); // old implementation
  fErrors(2) = ((s+b)>0 && s>0? fValues(2)*TMath::Sqrt(be*be + TMath::Power(se*(s+2*b)/s, 2)) / 2 / (s+b) : 0);
}

inline void AliDielectronSignalBase::SetFWHM()
{
  // calculate the fwhm
  if(!fPeakShapeObj) return;

  
  // case for TF1
  if(fPeakShapeObj->IsA() == TF1::Class()) {
    TF1* fit  = (TF1*) fPeakShapeObj->Clone("fit");
    TF1* pfit = (TF1*) fit->Clone("pfit");
    TF1* mfit = (TF1*) fit->Clone("mfit");
    for (Int_t i=0; i<fit->GetNpar(); i++) {
      pfit->SetParameter(i,fit->GetParameter(i) + fit->GetParError(i));
      mfit->SetParameter(i,fit->GetParameter(i) - fit->GetParError(i));
    }
    Double_t maxX   = fit->GetMaximumX();
    Double_t maxY   = fit->GetHistogram()->GetMaximum();
    Double_t xAxMin = fit->GetXmin();
    Double_t xAxMax = fit->GetXmax();
    // fwhms 
    Double_t fwhmMin  = fit->GetX(.5*maxY, xAxMin, maxX);
    Double_t fwhmMax  = fit->GetX(.5*maxY, maxX, xAxMax);
    Double_t pfwhmMin = pfit->GetX(.5*maxY, xAxMin, maxX);
    Double_t pfwhmMax = pfit->GetX(.5*maxY, maxX, xAxMax);
    Double_t mfwhmMin = mfit->GetX(.5*maxY, xAxMin, maxX);
    Double_t mfwhmMax = mfit->GetX(.5*maxY, maxX, xAxMax);
    Double_t pError = TMath::Abs( (fwhmMax-fwhmMin) - (pfwhmMax-pfwhmMin) );
    Double_t mError = TMath::Abs( (fwhmMax-fwhmMin) - (mfwhmMax-mfwhmMin) );
    fValues(5) = (fwhmMax-fwhmMin);
    fErrors(5) = (pError>=mError ? pError : mError);
    delete fit;
    delete pfit;
    delete mfit;
  }
  else if(fPeakShapeObj->IsA() == TH1F::Class()) {
    // th1 calculation
    TH1F *hist = (TH1F*) fPeakShapeObj->Clone("hist");
    Int_t bin1 = hist->FindFirstBinAbove(hist->GetMaximum()/2);
    Int_t bin2 = hist->FindLastBinAbove(hist->GetMaximum()/2);
    fValues(5) = hist->GetBinCenter(bin2) - hist->GetBinCenter(bin1);
    fErrors(5) = 0.0; // not defined
    delete hist;
  }
}

inline void AliDielectronSignalBase::SetBackgroundEnumMap(){

  TString method;
  
  method = "FittedMC";
  MapBackgroundMethod.Add((Long64_t) method.Hash(), (Long64_t) AliDielectronSignalBase::kFittedMC);
  method = "Fitted";
  MapBackgroundMethod.Add((Long64_t) method.Hash(), (Long64_t) AliDielectronSignalBase::kFitted);
  method = "LikeSign";
  MapBackgroundMethod.Add((Long64_t) method.Hash(), (Long64_t) AliDielectronSignalBase::kLikeSign);
  method = "LikeSignArithm";
  MapBackgroundMethod.Add((Long64_t) method.Hash(), (Long64_t) AliDielectronSignalBase::kLikeSignArithm);
  method = "LikeSignRcorr";
  MapBackgroundMethod.Add((Long64_t) method.Hash(), (Long64_t) AliDielectronSignalBase::kLikeSignRcorr);
  method = "LikeSignArithmRcorr";
  MapBackgroundMethod.Add((Long64_t) method.Hash(),(Long64_t) AliDielectronSignalBase::kLikeSignArithmRcorr);
  method = "LikeSignFit";
  MapBackgroundMethod.Add((Long64_t) method.Hash(),(Long64_t) AliDielectronSignalBase::kLikeSignFit);
  method = "EventMixing";
  MapBackgroundMethod.Add((Long64_t) method.Hash(),(Long64_t) AliDielectronSignalBase::kEventMixing);
  method = "EventMixingFit";
  MapBackgroundMethod.Add((Long64_t) method.Hash(),(Long64_t) AliDielectronSignalBase::kEventMixingFit);
  method = "Rotation";
  MapBackgroundMethod.Add((Long64_t) method.Hash(),(Long64_t) AliDielectronSignalBase::kRotation);

}

inline void AliDielectronSignalBase::SetSignalExtractionEnumMap(){

  TString method;

  method = "BinCounting";
  MapSignalExtractionMethod.Add((Long64_t) method.Hash(), (Long64_t) AliDielectronSignalBase::kBinCounting);
  method = "MCScaledMax";
  MapSignalExtractionMethod.Add((Long64_t) method.Hash(), (Long64_t) AliDielectronSignalBase::kMCScaledMax);
  method = "MCScaledInt";
  MapSignalExtractionMethod.Add((Long64_t) method.Hash(), (Long64_t) AliDielectronSignalBase::kMCScaledInt);
  method = "MCFitted";
  MapSignalExtractionMethod.Add((Long64_t) method.Hash(), (Long64_t) AliDielectronSignalBase::kMCFitted);
  method = "CrystalBall";
  MapSignalExtractionMethod.Add((Long64_t) method.Hash(), (Long64_t) AliDielectronSignalBase::kCrystalBall);
  method = "Gaus";
  MapSignalExtractionMethod.Add((Long64_t) method.Hash(), (Long64_t) AliDielectronSignalBase::kGaus);
  
}

#endif
