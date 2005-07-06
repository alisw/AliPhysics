#ifndef ALIGENMUONLIB_H
#define ALIGENMUONLIB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Library class for particle pt and y distributions used for 
// muon spectrometer simulations.
// To be used with AliGenParam.
//
// andreas.morsch@cern.ch

#include "AliGenLib.h"

class AliGenMUONlib :
  public AliGenLib
{
 public:
  enum constants{kPhi, kOmega, kEta, kJpsi, kJpsiFamily, kPsiP, kJpsiFromB, kUpsilon, kUpsilonFamily,
		   kUpsilonP, kUpsilonPP, kCharm, kBeauty, kPion, kKaon};
    
    GenFunc   GetPt(Int_t param, const char* tname=0) const;
    GenFunc   GetY (Int_t param, const char* tname=0) const;
    GenFuncIp GetIp(Int_t param, const char* tname=0) const;
 private:
    
// pions
    static Double_t PtPion(Double_t *px, Double_t *dummy);
    static Double_t PtScal(Double_t pt, Int_t np);
    static Double_t YPion( Double_t *py, Double_t *dummy);
    static Int_t    IpPion(TRandom *ran);
// kaons
    static Double_t PtKaon(Double_t *px, Double_t *dummy);
    static Double_t YKaon( Double_t *py, Double_t *dummy);
    static Int_t    IpKaon(TRandom *ran);
// Phi
    static Double_t PtPhi( Double_t *px, Double_t *dummy);
    static Double_t YPhi( Double_t *px, Double_t *dummy);
    static Int_t    IpPhi(TRandom *ran);
// Omega
    static Double_t PtOmega( Double_t *px, Double_t *dummy);
    static Double_t YOmega( Double_t *px, Double_t *dummy);
    static Int_t    IpOmega(TRandom *ran);
// Eta
    static Double_t PtEta( Double_t *px, Double_t *dummy);
    static Double_t YEta( Double_t *px, Double_t *dummy);
    static Int_t    IpEta(TRandom *ran);
// J/Psi     
    static Double_t PtJpsi( Double_t *px, Double_t *dummy);
    static Double_t PtJpsiCDFscaled( Double_t *px, Double_t *dummy);
    static Double_t YJpsi(Double_t *py, Double_t *dummy);
    static Double_t PtJpsiPbPb( Double_t *px, Double_t *dummy);
    static Double_t PtJpsiBPbPb( Double_t *px, Double_t *dummy);
    static Double_t YJpsiPbPb(Double_t *py, Double_t *dummy);
    static Double_t YJpsiCDFscaled(Double_t *py, Double_t *dummy);
    static Double_t PtJpsiPP( Double_t *px, Double_t *dummy);
    static Double_t YJpsiPP(Double_t *py, Double_t *dummy);
    static Double_t YJpsiBPbPb(Double_t *py, Double_t *dummy);
    static Int_t    IpJpsi(TRandom *ran);
    static Int_t    IpJpsiFamily(TRandom *ran);
    static Int_t    IpPsiP(TRandom *ran);
    static Double_t PtJpsiFlat( Double_t *px, Double_t *dummy );
    static Double_t YJpsiFlat(Double_t *py, Double_t *dummy);

// Upsilon    
    static Double_t PtUpsilon( Double_t *px, Double_t *dummy );
    static Double_t PtUpsilonCDFscaled( Double_t *px, Double_t *dummy );
    static Double_t YUpsilon(Double_t *py, Double_t *dummy);
    static Double_t YUpsilonCDFscaled(Double_t *py, Double_t *dummy);
    static Double_t PtUpsilonPbPb( Double_t *px, Double_t *dummy );
    static Double_t YUpsilonPbPb(Double_t *py, Double_t *dummy);
    static Double_t PtUpsilonPP( Double_t *px, Double_t *dummy );
    static Double_t YUpsilonPP(Double_t *py, Double_t *dummy);
    static Int_t    IpUpsilon(TRandom *ran);
    static Int_t    IpUpsilonFamily(TRandom *ran);
    static Int_t    IpUpsilonP(TRandom *ran);
    static Int_t    IpUpsilonPP(TRandom *ran);
    static Double_t PtUpsilonFlat( Double_t *px, Double_t *dummy );
    static Double_t YUpsilonFlat(Double_t *py, Double_t *dummy);
//
// Charm    
    static Double_t PtCharm( Double_t *px, Double_t *dummy );
    static Double_t PtCharmCentral( Double_t *px, Double_t *dummy );
    static Double_t YCharm(Double_t *py, Double_t *dummy);
    static Int_t    IpCharm(TRandom *ran);
//
// Beauty
    static Double_t PtBeauty( Double_t *px, Double_t *dummy );
    static Double_t PtBeautyCentral( Double_t *px, Double_t *dummy );
    static Double_t YBeauty(Double_t *py, Double_t *dummy);
    static Int_t    IpBeauty(TRandom *ran);
//
    
    static Float_t Interpolate(Float_t x, Float_t* y, Float_t x0, 
			Float_t dx,
			Int_t n, Int_t no);
    
    ClassDef(AliGenMUONlib,0) // Library providing y and pT parameterisations
};
#endif







