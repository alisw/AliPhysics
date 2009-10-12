#ifndef ALIBTOJPSITOELECDFFITHANDLER_H
#define ALIBTOJPSITOELECDFFITHANDLER_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//___________________________________________________________________________
//                  Class AliBtoJPSItoEleCDFfitHandler
//            Class to perform unbinned log-likelihood fit
//      
//                         Origin: C. Di Giglio
//     Contact: Carmelo.Digiglio@ba.infn.it; Giuseppe.Bruno@ba.infn.it
//____________________________________________________________________________

#include <TNamed.h>
#include <TBits.h>
#include "TMath.h"
#include "AliBtoJPSItoEleCDFfitFCN.h"
#include "AliLog.h"

class TBits; 

class AliBtoJPSItoEleCDFfitHandler : public TNamed {
 public:
  //
  AliBtoJPSItoEleCDFfitHandler();
  AliBtoJPSItoEleCDFfitHandler& operator= (const AliBtoJPSItoEleCDFfitHandler& c);
  AliBtoJPSItoEleCDFfitHandler(const AliBtoJPSItoEleCDFfitHandler& c);
  AliBtoJPSItoEleCDFfitHandler(Double_t* decaytime, Double_t* invariantmass, Int_t ncand);
  ~AliBtoJPSItoEleCDFfitHandler(); 
  //Double_t Up() const { return fUp*TMath::Sqrt(fLikely->GetIntegral()); }
  Double_t Up() const { return fUp; }
  void SetErrorDef(Double_t up) {fUp = up;}
  void SetPrintStatus(Bool_t okPrintStatus) { fPrintStatus = okPrintStatus; } 
  Bool_t GetPrintStatus() { return fPrintStatus ; }
  void SetParamStartValues (Double_t*);
  Double_t* GetStartParamValues() { return fParamStartValues; }
  TBits GetFixedParamList() { return fIsParamFixed; }
  void FixParam(UInt_t param, Bool_t value) { fIsParamFixed.SetBitNumber(param,value); }
  void FixAllParam(Bool_t value) { for(UInt_t par=0;par<16;par++) fIsParamFixed.SetBitNumber(par,value); }
  Bool_t IsParamFixed(UInt_t param) { return fIsParamFixed.TestBitNumber(param); }
  void SetResolutionConstants();
  void SetCrystalBallFunction(Bool_t okCB);
  void SetMassWndHigh(Double_t limit);
  void SetMassWndLow(Double_t limit);

  Double_t operator()(const Double_t* par) const ;
  void CdfFCN(Int_t & /* npar */, Double_t * /* gin */, Double_t &f, Double_t *par, Int_t /* iflag */);

  Double_t* Decaytime() const         { return fX; }
  Double_t* InvariantMass() const     { return fM; }
  AliBtoJPSItoEleCDFfitFCN* LikelihoodPointer() const { return fLikely; }
  Int_t DoMinimization();

 private:
  //
  TBits fIsParamFixed;                               //array of bits: 0 = param free; 1 = param fixed;
  Bool_t fPrintStatus;                               //flag to enable the prit out of the algorithm at each step
  Double_t fParamStartValues[16];                    //array of parameters input value
  Double_t fUp;                                      //error definition 
  Double_t* fX; 	                     	     //pseudo-proper decay time X
  Double_t* fM;                                      //invariant mass M
  AliBtoJPSItoEleCDFfitFCN* fLikely;                 //Log likelihood function
  Int_t fNcand;                                      //number of candidates
  //
  ClassDef(AliBtoJPSItoEleCDFfitHandler,1);

}; 
#endif
