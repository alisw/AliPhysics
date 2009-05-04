#ifndef ALIPHOSFASTALTROFIT_H
#define ALIPHOSFASTALTROFIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */

/* History of svn commits
 * $Log$
 */


//_________________________________________________________________________
//  Procedure of fast altro fitting     
//                  
//*-- Author:  Aleksei Pavlinov; IHEP, Protvino, Russia

#include <TNamed.h>
#include <TArrayD.h>
// --- ROOT system ---
class TCanvas;
class TVirtualPad;
class TF1;
class TH1F;
//class TLegend;
//class TPaveStats;

class AliPHOSFastAltroFit : public TNamed {

public:

  AliPHOSFastAltroFit();
  AliPHOSFastAltroFit(const char* name, const char* title, const Double_t tau);
  virtual ~AliPHOSFastAltroFit(); 

  void SetTau(const Double_t tau) {fTau = tau;}

  void FastFit(TH1F* h, Double_t sig, Double_t tau, Double_t ped);
  void FastFit(Int_t* t, Int_t* y, Int_t n, Double_t sig, Double_t tau, Double_t ped);
  void Reset();

  void GetFitResult(Double_t &amp, Double_t &eamp, Double_t &t0, Double_t &et0, 
                    Double_t &chi2,Int_t &ndf);
  Double_t GetSig()    const {return fSig;}
  Double_t GetTau()    const {return fTau;}
  Double_t GetN()      const {return fN;}
  Double_t GetPed()    const {return fPed;}

  Double_t GetEnergy() const {return fAmp;}
  Double_t GetAmp()    const {return GetEnergy();}
  Double_t GetAmpErr() const {return fAmpErr;}
  Double_t GetTime()   const {return fT0;}
  Double_t GetT0()     const {return GetTime();}
  Double_t GetT0Err()  const {return fT0Err;}
  Double_t GetChi2()   const {return fChi2;}
  Int_t    GetNDF()    const {return fNDF;}
  Int_t    GetNfit()   const {return fNfit;}
  void     GetFittedPoints(Int_t &nfit, Double_t* ar[2]);

  static void DeductPedestal(Int_t* t, Int_t* y, Int_t n, Double_t ped, Double_t tau, 
  Double_t* tn, Double_t* yn, Int_t &nn);
  static void FastFit(Double_t* t, Double_t* y, Int_t n, Double_t sig, Double_t tau,
  Double_t &amp, Double_t &eamp, Double_t &t0, Double_t &et0, Double_t &chi2);
  static Bool_t QuadraticRoots(Double_t a, Double_t b, Double_t c, Double_t &x1, Double_t &x2);
  static void Amplitude(Double_t* t, Double_t* y, Int_t n, Double_t sig, Double_t tau, 
  Double_t t0, Double_t &amp, Double_t &chi2);
  static void CalculateParsErrors(Double_t* t, Double_t* y, Int_t n, Double_t sig, Double_t tau,
				 Double_t &amp, Double_t &t0, Double_t &eamp, Double_t &et0);

  // Drawing for QA
  TCanvas* DrawFastFunction(); // *MENU*
  static Double_t StdResponseFunction(Double_t *x, Double_t *par); 

private:
  AliPHOSFastAltroFit(const AliPHOSFastAltroFit &obj);
  AliPHOSFastAltroFit& operator= (const AliPHOSFastAltroFit &obj);

protected:
  Double_t fSig;
  Double_t fTau; // filter time response
  Double_t fN;   // order of function (equal 2)
  Double_t fPed; // pedestal

  Double_t fAmp;
  Double_t fAmpErr; 
  Double_t fT0;
  Double_t fT0Err;
  Double_t fChi2;
  Int_t    fNDF;

  // Working variable
  Int_t     fNfit;   //! number points for fit
  Double_t* fTfit;   //! points for fit after selection
  Double_t* fAmpfit; //!
  // 
  TF1*      fStdFun; //! function for drawing

  ClassDef(AliPHOSFastAltroFit,1) // Class for fast altro fitting
};

#endif // ALIPHOSFASTALTROFIT_H
