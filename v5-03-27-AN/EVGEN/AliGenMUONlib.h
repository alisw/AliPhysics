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
		   kUpsilonP, kUpsilonPP, kCharm, kBeauty, kPion, kKaon, kChic, kChic0, kChic1, kChic2 }; 
    
    GenFunc   GetPt(Int_t param, const char* tname=0) const;
    GenFunc   GetY (Int_t param, const char* tname=0) const;
    GenFuncIp GetIp(Int_t param, const char* tname=0) const;
 private:
    
// pions
    static Double_t PtPion(const Double_t *px, const Double_t *dummy);
    static Double_t PtScal(Double_t pt, Int_t np);
    static Double_t YPion( const Double_t *py, const Double_t *dummy);
    static Int_t    IpPion(TRandom *ran);
// kaons
    static Double_t PtKaon(const Double_t *px, const Double_t *dummy);
    static Double_t YKaon( const Double_t *py, const Double_t *dummy);
    static Int_t    IpKaon(TRandom *ran);
//  XZhang 20110621
    static Double_t PtPionPos2010PP(const Double_t *px, const Double_t *dummy);
    static Double_t PtPionNeg2010PP(const Double_t *px, const Double_t *dummy);
    static Double_t PtKaonPos2010PP(const Double_t *px, const Double_t *dummy);
    static Double_t PtKaonNeg2010PP(const Double_t *px, const Double_t *dummy);
    static Double_t YKaonPion2010PP(const Double_t *px, const Double_t *dummy);
    static Int_t    IpPionPos(TRandom *ran);
    static Int_t    IpPionNeg(TRandom *ran);
    static Int_t    IpKaonPos(TRandom *ran);
    static Int_t    IpKaonNeg(TRandom *ran);
// Phi
    static Double_t PtPhi( const Double_t *px, const Double_t *dummy);
    static Double_t YPhi( const  Double_t *px, const Double_t *dummy);
    static Int_t    IpPhi(TRandom *ran);
// Omega
    static Double_t PtOmega( const Double_t *px, const Double_t *dummy);
    static Double_t YOmega( const Double_t *px, const Double_t *dummy);
    static Int_t    IpOmega(TRandom *ran);
// Eta
    static Double_t PtEta( const Double_t *px, const Double_t *dummy);
    static Double_t YEta( const Double_t *px, const Double_t *dummy);
    static Int_t    IpEta(TRandom *ran);
// J/Psi     
    static Double_t PtJpsiPPdummy(Double_t px, Double_t en);
    static Double_t PtJpsiPP7000(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPP2760(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPP8800(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPbPb2760ShFdummy(Double_t px, Int_t n);
    static Double_t PtJpsiPbPb2760(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPbPb2760c1(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPbPb2760c2(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPbPb2760c3(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPbPb2760c4(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPbPb2760c5(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPbPb2760c6(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPbPb2760c7(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPbPb2760c8(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPbPb2760c9(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPbPb2760c10(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPbPb2760c11(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPPb8800ShFdummy(Double_t px, Int_t n);
    static Double_t PtJpsiPPb8800(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPPb8800c1(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPPb8800c2(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPPb8800c3(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPPb8800c4(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPbP8800ShFdummy(Double_t px, Int_t n);
    static Double_t PtJpsiPbP8800(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPbP8800c1(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPbP8800c2(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPbP8800c3(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPbP8800c4(const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsi( const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiCDFscaled( const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiCDFscaledPP( const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiCDFscaledPP10( const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiCDFscaledPP9( const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiCDFscaledPP7( const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiCDFscaledPP4( const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiCDFscaledPP3( const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiCDFscaledPP2( const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiCDFscaledPPb9( const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiCDFscaledPbP9( const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiCDFscaledPbPb4( const Double_t *px, const Double_t *dummy);
    static Double_t YJpsi(const Double_t *py, const Double_t *dummy);
    static Double_t PtJpsiPbPb( const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiBPbPb( const Double_t *px, const Double_t *dummy);

    static Double_t YJpsiPPdummy(Double_t px, Double_t en);
    static Double_t YJpsiPPpoly(Double_t px, Double_t en);
    static Double_t YJpsiPP7000(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPP2760(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPPpoly7000(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPPpoly2760(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPP8800(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPbPb2760ShFdummy(Double_t px, Int_t n);
    static Double_t YJpsiPbPb2760(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPbPb2760c1(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPbPb2760c2(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPbPb2760c3(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPbPb2760c4(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPbPb2760c5(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPbPb2760c6(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPbPb2760c7(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPbPb2760c8(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPbPb2760c9(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPbPb2760c10(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPbPb2760c11(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPP8800dummy(Double_t px);
    static Double_t YJpsiPPb8800ShFdummy(Double_t px, Int_t n);
    static Double_t YJpsiPPb8800(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPPb8800c1(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPPb8800c2(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPPb8800c3(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPPb8800c4(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPbP8800(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPbP8800c1(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPbP8800c2(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPbP8800c3(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPbP8800c4(const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPbPb(const Double_t *py, const Double_t *dummy);
    static Double_t YJpsiCDFscaled(const Double_t *py, const Double_t *dummy);
    static Double_t YJpsiCDFscaledPP( const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiCDFscaledPP10( const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiCDFscaledPP9( const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiCDFscaledPP9dummy(Double_t px);
    static Double_t YJpsiCDFscaledPP7( const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiCDFscaledPP4( const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiCDFscaledPP3( const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiCDFscaledPP2( const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiCDFscaledPPb9( const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiCDFscaledPbP9( const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiCDFscaledPbPb4( const Double_t *px, const Double_t *dummy);
    static Double_t PtJpsiPP( const Double_t *px, const Double_t *dummy);
    static Double_t YJpsiPP(const Double_t *py, const Double_t *dummy);
    static Double_t YJpsiBPbPb(const Double_t *py, const Double_t *dummy);
    static Int_t    IpJpsi(TRandom *ran);
    static Int_t    IpJpsiFamily(TRandom *ran);
    static Int_t    IpPsiP(TRandom *ran);
    static Double_t PtJpsiFlat( const Double_t *px, const Double_t *dummy );
    static Double_t YJpsiFlat(const Double_t *py, const Double_t *dummy);

// Upsilon    
    static Double_t PtUpsilonPPdummy(Double_t px, Double_t en);
    static Double_t PtUpsilonPP7000(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPP2760(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPP8800(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPbPb2760ShFdummy(Double_t px, Int_t n);
    static Double_t PtUpsilonPbPb2760(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPbPb2760c1(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPbPb2760c2(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPbPb2760c3(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPbPb2760c4(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPbPb2760c5(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPbPb2760c6(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPbPb2760c7(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPbPb2760c8(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPbPb2760c9(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPbPb2760c10(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPbPb2760c11(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPPb8800ShFdummy(Double_t px, Int_t n);
    static Double_t PtUpsilonPPb8800(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPPb8800c1(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPPb8800c2(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPPb8800c3(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPPb8800c4(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPbP8800ShFdummy(Double_t px, Int_t n);
    static Double_t PtUpsilonPbP8800(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPbP8800c1(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPbP8800c2(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPbP8800c3(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPbP8800c4(const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilon( const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonCDFscaled( const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonCDFscaledPP( const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonCDFscaledPP10( const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonCDFscaledPP9( const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonCDFscaledPP7( const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonCDFscaledPP4( const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonCDFscaledPPb9( const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonCDFscaledPbP9( const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonCDFscaledPbPb4( const Double_t *px, const Double_t *dummy );

    static Double_t YUpsilonPPdummy(Double_t px, Double_t en);
    static Double_t YUpsilonPPpoly(Double_t px, Double_t en);
    static Double_t YUpsilonPP7000(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPP2760(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPPpoly7000(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPPpoly2760(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPP8800(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPbPb2760ShFdummy(Double_t px, Int_t n);
    static Double_t YUpsilonPbPb2760(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPbPb2760c1(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPbPb2760c2(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPbPb2760c3(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPbPb2760c4(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPbPb2760c5(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPbPb2760c6(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPbPb2760c7(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPbPb2760c8(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPbPb2760c9(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPbPb2760c10(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPbPb2760c11(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPP8800dummy(Double_t px);
    static Double_t YUpsilonPPb8800ShFdummy(Double_t px, Int_t n);
    static Double_t YUpsilonPPb8800(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPPb8800c1(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPPb8800c2(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPPb8800c3(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPPb8800c4(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPbP8800(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPbP8800c1(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPbP8800c2(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPbP8800c3(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPbP8800c4(const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilon(const Double_t *py, const Double_t *dummy);
    static Double_t YUpsilonCDFscaled(const Double_t *py, const Double_t *dummy);
    static Double_t YUpsilonCDFscaledPP( const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonCDFscaledPP10( const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonCDFscaledPP9( const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonCDFscaledPP9dummy(Double_t px);
    static Double_t YUpsilonCDFscaledPP7( const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonCDFscaledPP4( const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonCDFscaledPPb9( const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonCDFscaledPbP9( const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonCDFscaledPbPb4( const Double_t *px, const Double_t *dummy );
    static Double_t PtUpsilonPbPb( const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPbPb(const Double_t *py, const Double_t *dummy);
    static Double_t PtUpsilonPP( const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonPP(const Double_t *py, const Double_t *dummy);
    static Int_t    IpUpsilon(TRandom *ran);
    static Int_t    IpUpsilonFamily(TRandom *ran);
    static Int_t    IpUpsilonP(TRandom *ran);
    static Int_t    IpUpsilonPP(TRandom *ran);
    static Double_t PtUpsilonFlat( const Double_t *px, const Double_t *dummy );
    static Double_t YUpsilonFlat(const Double_t *py, const Double_t *dummy);
//
// Charm    
    static Double_t PtCharm( const Double_t *px, const Double_t *dummy );
    static Double_t PtCharmCentral( const Double_t *px, const Double_t *dummy );
    static Double_t PtCharmF0M0S0PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtCharmF1M0S0PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtCharmF2M0S0PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtCharmF0M1S0PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtCharmF0M2S0PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtCharmF0M0S1PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtCharmF0M0S2PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtCharmF0M0S3PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtCharmF0M0S4PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtCharmF0M0S5PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtCharmF0M0S6PP( const Double_t *px, const Double_t *dummy );
    static Double_t YCharm(const Double_t *py, const Double_t *dummy);
    static Double_t YCharmF0M0S0PP(const Double_t *py, const Double_t *dummy);
    static Double_t YCharmF1M0S0PP(const Double_t *py, const Double_t *dummy);
    static Double_t YCharmF2M0S0PP(const Double_t *py, const Double_t *dummy);
    static Double_t YCharmF0M1S0PP(const Double_t *py, const Double_t *dummy);
    static Double_t YCharmF0M2S0PP(const Double_t *py, const Double_t *dummy);
    static Double_t YCharmF0M0S1PP(const Double_t *py, const Double_t *dummy);
    static Double_t YCharmF0M0S2PP(const Double_t *py, const Double_t *dummy);
    static Double_t YCharmF0M0S3PP(const Double_t *py, const Double_t *dummy);
    static Double_t YCharmF0M0S4PP(const Double_t *py, const Double_t *dummy);
    static Double_t YCharmF0M0S5PP(const Double_t *py, const Double_t *dummy);
    static Double_t YCharmF0M0S6PP(const Double_t *py, const Double_t *dummy);
    static Int_t    IpCharm(TRandom *ran);
//
// Beauty
    static Double_t PtBeauty( const Double_t *px, const Double_t *dummy );
    static Double_t PtBeautyF0M0S0PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtBeautyF1M0S0PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtBeautyF2M0S0PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtBeautyF0M1S0PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtBeautyF0M2S0PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtBeautyF0M0S1PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtBeautyF0M0S2PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtBeautyF0M0S3PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtBeautyF0M0S4PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtBeautyF0M0S5PP( const Double_t *px, const Double_t *dummy );
    static Double_t PtBeautyF0M0S6PP( const Double_t *px, const Double_t *dummy );
    static Double_t YBeauty(const Double_t *py, const Double_t *dummy);
    static Double_t YBeautyF0M0S0PP(const Double_t *py, const Double_t *dummy);
    static Double_t YBeautyF1M0S0PP(const Double_t *py, const Double_t *dummy);
    static Double_t YBeautyF2M0S0PP(const Double_t *py, const Double_t *dummy);
    static Double_t YBeautyF0M1S0PP(const Double_t *py, const Double_t *dummy);
    static Double_t YBeautyF0M2S0PP(const Double_t *py, const Double_t *dummy);
    static Double_t YBeautyF0M0S1PP(const Double_t *py, const Double_t *dummy);
    static Double_t YBeautyF0M0S2PP(const Double_t *py, const Double_t *dummy);
    static Double_t YBeautyF0M0S3PP(const Double_t *py, const Double_t *dummy);
    static Double_t YBeautyF0M0S4PP(const Double_t *py, const Double_t *dummy);
    static Double_t YBeautyF0M0S5PP(const Double_t *py, const Double_t *dummy);
    static Double_t YBeautyF0M0S6PP(const Double_t *py, const Double_t *dummy);
    static Double_t PtBeautyCentral( const Double_t *px, const Double_t *dummy );
    static Int_t    IpBeauty(TRandom *ran);
//

   // Chi 1c 2c
   static Double_t PtChic0( const Double_t *px, const Double_t *dummy);
   static Double_t YChic0(const Double_t *py, const Double_t *dummy);
   static Int_t    IpChic0(TRandom *ran);

   static Double_t PtChic1( const Double_t *px, const Double_t *dummy);
   static Double_t YChic1(const Double_t *py, const Double_t *dummy);
   static Int_t    IpChic1(TRandom *ran);

   static Double_t PtChic2( const Double_t *px, const Double_t *dummy);
   static Double_t YChic2(const Double_t *py, const Double_t *dummy);
   static Int_t    IpChic2(TRandom *ran);

   static Double_t PtChic( const Double_t *px, const Double_t *dummy);
   static Double_t YChic(const Double_t *py, const Double_t *dummy);
   static Int_t    IpChic(TRandom *ran);

//

    
    static Float_t Interpolate(Float_t x, Float_t* y, Float_t x0, 
			Float_t dx,
			Int_t n, Int_t no);
    
    ClassDef(AliGenMUONlib,0) // Library providing y and pT parameterisations
};
#endif







