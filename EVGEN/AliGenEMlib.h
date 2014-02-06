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
    
  enum particles{kPizero, kEta, kRho, kOmega, kEtaprime, kPhi, kJpsi};
  enum centrality{kpp=0x0, k0005=0x1, k0510=0x2, k1020=0x3, k2030=0x4, k3040=0x5, k4050=0x6, k5060=0x7, k0010=0x8, k2040=0x9, k4060=0xA, k6080=0xB, k0020=0xC, k0040=0xD, k2080=0xE, k4080=0xF, kCentralities=0x10};
  enum PtParamSet{kPizero7TeVpp=0x10, kPizeroEta7TeVpp=0x20, kPizero7TeVpplow=0x30, kPizeroEta7TeVpplow=0x40, kPizero7TeVpphigh=0x50, kPizeroEta7TeVpphigh=0x60, kPizero2760GeVpp=0x70, kPizeroEta2760GeVpp=0x80, kPizero2760GeVpplow=0x90, kPizeroEta2760GeVpplow=0xA0, kPizero2760GeVpphigh=0xB0, kPizeroEta2760GeVpphigh=0xC0, kPichargedPbPb=0xD0, kPizeroPbPb=0xE0, kPichargedPPb=0xF0 };
  enum v2Sys{kLoV2Sys=-1, kNoV2Sys=0, kUpV2Sys=+1};
  
  AliGenEMlib() { } ;

  static void SelectParams(Int_t ptSelect, Int_t centSelect=kpp, Int_t v2sys=kNoV2Sys) 
  { fgSelectedPtParam=ptSelect; fgSelectedCentrality=centSelect; fgSelectedV2Systematic=v2sys; }

    GenFunc   GetPt(Int_t param, const char * tname=0) const;
    GenFunc   GetY(Int_t param, const char * tname=0) const;
    GenFuncIp GetIp(Int_t param, const char * tname=0) const;    
  GenFunc   GetV2(Int_t param, const char * tname=0) const;

  //private:

  // General functions

  static Int_t fgSelectedPtParam; // selected pT parameter
  static Int_t fgSelectedCentrality; // selected Centrality
  static Int_t fgSelectedV2Systematic; // selected v2 systematics, usefully values: -1,0,1


  static Double_t PtModifiedHagedornThermal(const Double_t pt, 
					    const Double_t c, 
					    const Double_t p0, 
					    const Double_t p1, 
					    const Double_t n,
					    const Double_t cT,
					    const Double_t T);


 
  static Double_t PtModifiedHagedornExp(const Double_t pt,
					const Double_t c,
					const Double_t p0,
					const Double_t p1,
					const Double_t p2,
					const Double_t n); 


  static Double_t PtModifiedHagedornExp2(const Double_t pt,
                                           const Double_t c,
                                           const Double_t a,
                                           const Double_t b,
                                           const Double_t p0,
                                           const Double_t p1,
                                           const Double_t d,
                                           const Double_t n);


  static Double_t PtTsallis(const Double_t pt,
			    const Double_t m,
			    const Double_t c,
			    const Double_t T,
			    const Double_t n);




  // Pizero
    static Int_t    IpPizero(TRandom *ran);
  static Double_t PtPizero(const Double_t *px, const Double_t *dummy);
    static Double_t YPizero(const Double_t *py, const Double_t *dummy);
  static Double_t V2Pizero(const Double_t *px, const Double_t *dummy);

  // Eta
    static Int_t    IpEta(TRandom *ran);
  static Double_t PtEta(const Double_t *px, const Double_t *dummy);
    static Double_t YEta(const Double_t *py, const Double_t *dummy);
  static Double_t V2Eta(const Double_t *px, const Double_t *dummy);

  // Rho
    static Int_t    IpRho(TRandom *ran);
  static Double_t PtRho(const Double_t *px, const Double_t *dummy);
    static Double_t YRho(const Double_t *py, const Double_t *dummy);
  static Double_t V2Rho(const Double_t *py, const Double_t *dummy);

  // Omega
    static Int_t    IpOmega(TRandom *ran);
  static Double_t PtOmega(const Double_t *px, const Double_t *dummy);
    static Double_t YOmega(const Double_t *py, const Double_t *dummy);
  static Double_t V2Omega(const Double_t *py, const Double_t *dummy);

  // Etaprime
    static Int_t    IpEtaprime(TRandom *ran);
  static Double_t PtEtaprime(const Double_t *px, const Double_t *dummy);
    static Double_t YEtaprime(const Double_t *py, const Double_t *dummy);
  static Double_t V2Etaprime(const Double_t *py, const Double_t *dummy);

  // Phi
    static Int_t    IpPhi(TRandom *ran);
  static Double_t PtPhi(const Double_t *px, const Double_t *dummy);
    static Double_t YPhi(const Double_t *py, const Double_t *dummy);
  static Double_t V2Phi(const Double_t *py, const Double_t *dummy);

  // Jpsi
  static Int_t    IpJpsi(TRandom *ran);
  static Double_t PtJpsi(const Double_t *px, const Double_t *dummy);
  static Double_t YJpsi(const Double_t *py, const Double_t *dummy);
  static Double_t V2Jpsi(const Double_t *py, const Double_t *dummy);

  // General
  //static Double_t PtFlat(const Double_t *px, const Double_t *dummy);
    static Double_t YFlat(Double_t y);
  static Double_t MtScal(Double_t pt, Int_t np);
  static Double_t V2Param(const Double_t *px, const Double_t *param);
  static Double_t V2Flat(const Double_t *px, const Double_t *param);
  static Double_t KEtScal(Double_t pt, Int_t np);

  static Double_t CrossOverLc(const double a, const double b, const double x);
  static Double_t CrossOverRc(const double a, const double b, const double x);

  static const Double_t fgkV2param[16][15]; // array of v2 parameter
  static const Double_t fgkHM[7];
  static const Double_t fgkMtFactor[2][7];

  ClassDef(AliGenEMlib,0)
};


#endif
