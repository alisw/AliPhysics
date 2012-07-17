#ifndef ALIGENEMLIB_H
#define ALIGENEMLIB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliGenEMlib.h 30052 2008-11-25 14:54:18Z morsch $ */

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Implementation of AliGenEMlib for electron, di-electron, and photon     //
// cocktail calculations.                                                  //
// It is based on AliGenGSIlib.                                            //
//                                                                         //
// Responsible: R.Averbeck@gsi.de                                          //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "AliGenLib.h"
class TRandom;

class AliGenEMlib :public AliGenLib {
public:
    GenFunc   GetPt(Int_t param, const char * tname=0) const;
    GenFunc   GetY(Int_t param, const char * tname=0) const;
    GenFuncIp GetIp(Int_t param, const char * tname=0) const;    
  GenFunc   GetV2(Int_t param, const char * tname=0) const;

  typedef enum {k0005=0, k0510, k1020, k2030, k3040, k4050, k5060, kpp} Centbins_t;
  enum particles{kPizero, kEta, kRho, kOmega, kEtaprime, kPhi};

private:

  // Pizero
    static Int_t    IpPizero(TRandom *ran);
    static Double_t PtPizero( const Double_t *px, const Double_t *dummy );
    static Double_t YPizero(const Double_t *py, const Double_t *dummy);
  static Double_t V2Pizero(const Double_t *px, const Double_t *dummy);

  // Eta
    static Int_t    IpEta(TRandom *ran);
    static Double_t PtEta( const Double_t *px, const Double_t *dummy );
    static Double_t YEta(const Double_t *py, const Double_t *dummy);
  static Double_t V2Eta(const Double_t *px, const Double_t *dummy);

  // Rho
    static Int_t    IpRho(TRandom *ran);
    static Double_t PtRho( const Double_t *px, const Double_t *dummy );
    static Double_t YRho(const Double_t *py, const Double_t *dummy);

  // Omega
    static Int_t    IpOmega(TRandom *ran);
    static Double_t PtOmega( const Double_t *px, const Double_t *dummy );
    static Double_t YOmega(const Double_t *py, const Double_t *dummy);

  // Etaprime
    static Int_t    IpEtaprime(TRandom *ran);
    static Double_t PtEtaprime( const Double_t *px, const Double_t *dummy );
    static Double_t YEtaprime(const Double_t *py, const Double_t *dummy);

  // Phi
    static Int_t    IpPhi(TRandom *ran);
    static Double_t PtPhi( const Double_t *px, const Double_t *dummy );
    static Double_t YPhi(const Double_t *py, const Double_t *dummy);

  // General
    static Double_t PtFlat(const Double_t *px, const Double_t *dummy);
    static Double_t YFlat(Double_t y);
  static Double_t MtScal(Double_t pt, Int_t np);
  static Double_t V2Param(const Double_t *px, const Double_t *param);
  static Double_t V2Flat(const Double_t *px, const Double_t *param);

  static Double_t CrossOverLc(const double a, const double b, const double x);
  static Double_t CrossOverRc(const double a, const double b, const double x);

  static const Double_t fv2param[2][8][9];
  static const Double_t fpTparam[8][14];
  static Centbins_t fCentbin;

  ClassDef(AliGenEMlib,0)
};

#endif
