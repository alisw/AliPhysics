//_________________________________________________________________________
//  Procedure of fast altro fitting     
//                  
//*-- Author:  Aleksei Pavlinov; IHEP, Protvino, Russia & WSU, Detroit, USA

#ifndef ALICALOFASTALTROFITV0_H
#define ALICALOFASTALTROFITV0_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of svn commits
 * $Log$
 */

#include <TNamed.h>
// --- ROOT system ---
class TCanvas;
class TVirtualPad;
class TF1;
class TH1F;

class AliCaloFastAltroFitv0 : public TNamed {

public:

  AliCaloFastAltroFitv0();
  AliCaloFastAltroFitv0(const char* name, const char* title,
  const Double_t sig=1.3, const Double_t tau=2.35, const Double_t n=2.);
  virtual ~AliCaloFastAltroFitv0(); 

  virtual void FastFit(Int_t* t, Int_t* y, Int_t nPoints, Double_t sig, Double_t tau, 
  Double_t n, Double_t ped, Double_t tMax);

  void FastFit(TH1F* h, Double_t sig, Double_t tau, 
  Double_t n, Double_t ped, Double_t tMax);

  void Reset();
  void SetSig(const Double_t sig) {fSig = sig;}
  void SetTau(const Double_t tau) {fTau = tau;}
  void SetN(const Double_t n)     {fN = n;}
  void SetParameters(const Double_t sig, const Double_t tau, const Double_t n) 
  {fSig = sig; fTau = tau; fN = n;}

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
  Int_t    GetNoFit()  const {return fNoFit;}

  void GetFitResult(Double_t &amp, Double_t &eamp, Double_t &t0, Double_t &et0, 
                    Double_t &chi2,Int_t &ndf) const;
  void     GetFittedPoints(Int_t &nfit, Double_t* ar[2]) const;

  // Drawing for QA
  TCanvas* DrawFastFunction(); // *MENU*
  static Double_t StdResponseFunction(const Double_t *x, const  Double_t *par); 

  static void CutRightPart(Int_t *t,Int_t *y,Int_t nPoints, Double_t tMax, Int_t &ii); 
  static void DeductPedestal(Int_t* t, Int_t* y, Int_t nPointsIn, Double_t ped, Double_t tau, 
  Double_t* tn, Double_t* yn, Int_t &nPointsOut);

  static void FastFit(const Double_t* t, const Double_t* y, const Int_t nPoints, 
                      const Double_t sig, const Double_t tau,
                      Double_t &amp, Double_t &eamp, Double_t &t0, Double_t &et0, Double_t &chi2);
  static Bool_t QuadraticRoots(const Double_t a, const Double_t b, const Double_t c, 
                               Double_t &x1, Double_t &x2);
  static void Amplitude(const Double_t* t, const Double_t* y, const Int_t nPoints, 
                        const Double_t sig, const Double_t tau, 
                        Double_t t0, Double_t &amp, Double_t &chi2);
  static void CalculateParsErrors(const Double_t* t, const Double_t* y, const Int_t nPoints, 
                                 const Double_t sig, const Double_t tau,
				 Double_t &amp, Double_t &t0, Double_t &eamp, Double_t &et0);
protected:
  Double_t fSig; // error in amplitude - used in chi2 calculation
  Double_t fTau; // first  fixed parameter od fitting function (should be - filter time response
  Double_t fN;   // second fixed parameter od fitting function (should be positive)
  Double_t fPed; // pedestal

  Double_t fAmp;    // amplitude
  Double_t fAmpErr; // amplitude error
  Double_t fT0;     // time
  Double_t fT0Err;  // time error
  Double_t fChi2;   // chi square
  Int_t    fNDF;    // number degree of freedom
  Int_t    fNoFit;  // no solution for square equation

  // Working variable
  Int_t     fNfit;   //! number points for fit
  Double_t* fTfit;   //! points for fit after selection - time bins
  Double_t* fAmpfit; //!                                - amplitudes
  // 
  TF1*      fStdFun; //! function for drawing

private:
  AliCaloFastAltroFitv0(const AliCaloFastAltroFitv0 &obj);
  AliCaloFastAltroFitv0& operator= (const AliCaloFastAltroFitv0 &obj);

  ClassDef(AliCaloFastAltroFitv0,1) // Class for fast altro fitting
};

#endif // ALICALOFASTALTROFITV0_H
