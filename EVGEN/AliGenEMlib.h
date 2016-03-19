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

  enum Particle_t{ kPizero=0, kEta=1, kRho0=2, kOmega=3, kEtaprime=4, kPhi=5, kJpsi=6,
       kSigma0=7, kK0s=8, kDeltaPlPl=9, kDeltaPl=10, kDeltaMi=11, kDeltaZero=12,
       kRhoPl=13, kRhoMi=14, kK0star=15, kDirectRealGamma=16, kDirectVirtGamma=17 };
  
  enum CollisionSystem_t {kpp900GeV=0x000, kpp2760GeV=0x64, kpp7TeV=0xC8, kpPb=0x12C, kPbPb=0x190};


  
  enum Centrality_t{ kpp = 0x0, k0005=0x1, k0510=0x2, k1020=0x3, k2030=0x4, k3040=0x5, k4050=0x6, k5060=0x7,
         k0010=0x8, k2040=0x9, k4060=0xA, k6080=0xB, k0020=0xC, k0040=0xD, k2080=0xE, k4080=0xF, kCentralities=0x10};
  
  enum PtParamSetPi0_t{ kPizeroParam, kPizeroParamlow, kPizeroParamhigh, kPichargedParam,
      kPichargedParamlow, kPichargedParamhigh, kPizeroParamAlter,
      kPizeroParamAlterlow, kPizeroParamAlterhigh, kNPi0Param,
      kPichargedParamNew, kPichargedParamOld };
  
  enum PtParamSetEta_t{ kEtaMtScal=0, kEtaParampp, kEtaParampplow, kEtaParampphigh,
            kEtaParamRatiopp, kEtaParamRatiopplow, kEtaParamRatiopphigh,
            kEtaParamPbPb, kEtaParamPPb };
  
  enum PtParamSetOmega_t{ kOmegaMtScal=0, kOmegaParampp, kOmegaParampplow, kOmegaParampphigh,
        kOmegaParamRatiopp, kOmegaParamRatiopplow, kOmegaParamRatiopphigh,
        kOmegaParamPbPb, kOmegaParamPPb };
  
  enum PtParamSetPhi_t{ kPhiMtScal=0, kPhiParampp, kPhiParampplow, kPhiParampphigh,
            kPhiParamPbPb, kPhiParamPPb, kPhiParamPPblow, kPhiParamPPbhigh };
  
  enum v2Sys_t{kLoV2Sys=-1, kNoV2Sys=0, kUpV2Sys=+1};
 
  AliGenEMlib() { } ;

  static void SelectParams( Int_t collisionSystem,
                            Int_t ptSelectPi0, 
                            Int_t ptSelectEta     = kEtaMtScal, 
                            Int_t ptSelectOmega   = kOmegaMtScal,
                            Int_t ptSelectPhi     = kPhiMtScal, 
                            Int_t centSelect      = kpp, 
                            Int_t v2sys           = kNoV2Sys) {  
    fgSelectedCollisionsSystem  = collisionSystem;
    fgSelectedPtParamPi0        = ptSelectPi0; 
    fgSelectedPtParamEta        = ptSelectEta; 
    fgSelectedPtParamOmega      = ptSelectOmega; 
    fgSelectedPtParamPhi        = ptSelectPhi; 
    fgSelectedCentrality        = centSelect; 
    fgSelectedV2Systematic      = v2sys; 
  }

  GenFunc   GetPt(Int_t param, const char * tname=0) const;
  GenFunc   GetY(Int_t param, const char * tname=0) const;
  GenFuncIp GetIp(Int_t param, const char * tname=0) const;    
  GenFunc   GetV2(Int_t param, const char * tname=0) const;

  //private:

  // General functions

  // General functions
  static Int_t fgSelectedCollisionsSystem; // selected pT parameter
  static Int_t fgSelectedPtParamPi0; // selected pT parameter
  static Int_t fgSelectedPtParamEta; // selected pT parameter
  static Int_t fgSelectedPtParamOmega; // selected pT parameter
  static Int_t fgSelectedPtParamPhi; // selected pT parameter
  static Int_t fgSelectedCentrality; // selected Centrality
  static Int_t fgSelectedV2Systematic; // selected v2 systematics, usefully values: -1,0,1


  static Double_t PtModifiedHagedornThermal(Double_t pt, 
                                            Double_t c, 
                                            Double_t p0, 
                                            Double_t p1, 
                                            Double_t n,
                                            Double_t cT,
                                            Double_t T);


  
  static Double_t PtModifiedHagedornExp(Double_t pt,
                                        Double_t c,
                                        Double_t p0,
                                        Double_t p1,
                                        Double_t p2,
                                        Double_t n); 


  static Double_t PtModifiedHagedornExp2( Double_t pt,
                                          Double_t c,
                                          Double_t a,
                                          Double_t b,
                                          Double_t p0,
                                          Double_t p1,
                                          Double_t d,
                                          Double_t n);


  static Double_t PtTsallis(Double_t pt,
                            Double_t m,
                            Double_t c,
                            Double_t T,
                            Double_t n);

  static Double_t PtParticleRatiopp(Double_t pt,
                                    Double_t m1,
                                    Double_t m2,
                                    Double_t c1,
                                    Double_t c2,
                                    Double_t T1,
                                    Double_t T2,
                                    Double_t n);
  
  static Double_t PtXQCD( Double_t pt,
                          Double_t a,
                          Double_t b,
                          Double_t c,
                          Double_t d,
                          Double_t e,
                          Double_t f);
  
  static Double_t PtModTsallis( Double_t pt,
                                Double_t a,
                                Double_t b,
                                Double_t c,
                                Double_t d,
                                Double_t e,
                                Double_t f,
                                Double_t g,
                                Double_t mass);
  
  static Double_t PtQCD(  Double_t pt,
                          Double_t a,
                          Double_t b,
                          Double_t c,
                          Double_t d,
                          Double_t e);

  static Double_t PtExponential(const Double_t *pt, const Double_t *param);
  static Double_t PtModifiedHagedornPowerlaw(const Double_t *pt, const Double_t *param);
  static Double_t PtDoublePowerlaw(const Double_t *pt, const Double_t *param);
  static Double_t IntegratedKrollWada(const Double_t *mh, const Double_t *);

  // direct gamma
  static Double_t PtPromptRealGamma(const Double_t *px, const Double_t *dummy);
  static Double_t PtPromptVirtGamma(const Double_t *px, const Double_t *dummy);
  static Double_t PtThermalRealGamma(const Double_t *px, const Double_t *dummy);
  static Double_t PtThermalVirtGamma(const Double_t *px, const Double_t *dummy);

  static Int_t    IpDirectRealGamma(TRandom *ran);
  static Double_t PtDirectRealGamma(const Double_t *px, const Double_t *dummy);
  static Double_t YDirectRealGamma(const Double_t *py, const Double_t *dummy);
  static Double_t V2DirectRealGamma(const Double_t *px, const Double_t *dummy);

  static Int_t    IpDirectVirtGamma(TRandom *ran);
  static Double_t PtDirectVirtGamma(const Double_t *px, const Double_t *dummy);
  static Double_t YDirectVirtGamma(const Double_t *py, const Double_t *dummy);
  static Double_t V2DirectVirtGamma(const Double_t *px, const Double_t *dummy);

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
  static Int_t    IpRho0(TRandom *ran);
  static Double_t PtRho0(const Double_t *px, const Double_t *dummy);
  static Double_t YRho0(const Double_t *py, const Double_t *dummy);
  static Double_t V2Rho0(const Double_t *py, const Double_t *dummy);


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

  // Sigma
  static Int_t    IpSigma(TRandom *ran);
  static Double_t PtSigma( const Double_t *px, const Double_t *dummy );
  static Double_t YSigma(const Double_t *py, const Double_t *dummy);
  static Double_t V2Sigma0( const Double_t *px, const Double_t *dummy );
  
  // K0short
  static Int_t    IpK0short(TRandom *ran);
  static Double_t PtK0short( const Double_t *px, const Double_t *dummy );
  static Double_t YK0short(const Double_t *py, const Double_t *dummy);
  static Double_t V2K0sshort( const Double_t *px, const Double_t *dummy );

  // Delta++
  static Int_t    IpDeltaPlPl(TRandom *ran);
  static Double_t PtDeltaPlPl( const Double_t *px, const Double_t *dummy );
  static Double_t YDeltaPlPl(const Double_t *py, const Double_t *dummy);
  static Double_t V2DeltaPlPl( const Double_t *px, const Double_t *dummy );
  
  // Delta+
  static Int_t    IpDeltaPl(TRandom *ran);
  static Double_t PtDeltaPl( const Double_t *px, const Double_t *dummy );
  static Double_t YDeltaPl(const Double_t *py, const Double_t *dummy);
  static Double_t V2DeltaPl( const Double_t *px, const Double_t *dummy );
  
  // Delta-
  static Int_t    IpDeltaMi(TRandom *ran);
  static Double_t PtDeltaMi( const Double_t *px, const Double_t *dummy );
  static Double_t YDeltaMi(const Double_t *py, const Double_t *dummy);
  static Double_t V2DeltaMi( const Double_t *px, const Double_t *dummy );
  
  // Delta0
  static Int_t    IpDeltaZero(TRandom *ran);
  static Double_t PtDeltaZero( const Double_t *px, const Double_t *dummy );
  static Double_t YDeltaZero(const Double_t *py, const Double_t *dummy);
  static Double_t V2DeltaZero( const Double_t *px, const Double_t *dummy );

  // Rho+
  static Int_t    IpRhoPl(TRandom *ran);
  static Double_t PtRhoPl( const Double_t *px, const Double_t *dummy );
  static Double_t YRhoPl(const Double_t *py, const Double_t *dummy);
  static Double_t V2RhoPl( const Double_t *px, const Double_t *dummy );

  // Rho-
  static Int_t    IpRhoMi(TRandom *ran);
  static Double_t PtRhoMi( const Double_t *px, const Double_t *dummy );
  static Double_t YRhoMi(const Double_t *py, const Double_t *dummy);
  static Double_t V2RhoMi( const Double_t *px, const Double_t *dummy );

  // K0*
  static Int_t    IpK0star(TRandom *ran);
  static Double_t PtK0star( const Double_t *px, const Double_t *dummy );
  static Double_t YK0star(const Double_t *py, const Double_t *dummy);
  static Double_t V2K0star( const Double_t *px, const Double_t *dummy );
  

  // General
  //static Double_t PtFlat(const Double_t *px, const Double_t *dummy);
  static Double_t YFlat(Double_t y);
  static Double_t MtScal(Double_t pt, Int_t np);
  static Double_t V2Param(const Double_t *px, const Double_t *param);
  static Double_t V2Flat(const Double_t *px, const Double_t *param);
  static Double_t KEtScal(Double_t pt, Int_t np, Int_t nq=2);
  static Double_t GetTAA(Int_t cent);

  static Double_t CrossOverLc(double a, double b, double x);
  static Double_t CrossOverRc(double a, double b, double x);

  static const Double_t fgkPtParam[kCentralities][10];                     // parameters of pi pt spectrum
  static const Double_t fgkModTsallisParamPi0PbPb[kCentralities][7];       // parameters for ModTsallis function for pi0 in PbPb 
  static const Double_t fgkModTsallisParamPiChargedPbPb[kCentralities][7]; // parameters for ModTsallis function for pi+- in PbPb 
  static const Double_t fgkV2param[kCentralities][16];                     // parameters of pi v2
  static const Double_t fgkRawPtOfV2Param[kCentralities][10];              // parameters of the raw pt spectrum of v2 analysys
  static const Double_t fgkThermPtParam[kCentralities][2];                 // parameters of thermal gamma pt
  static const Double_t fgkHM[16];                                         // particle masses
  static const Double_t fgkMtFactor[3][16];                                // mt scaling factor
  static const Double_t fgkParamSetPi07TeV[kNPi0Param][7];                 // parameters for pi0 in 7 TeV
  static const Double_t fgkParamSetPi02760GeV[kNPi0Param][7];              // parameters for pi0 in 2.76 TeV
  static const Double_t fgkParamSetPi0900GeV[kNPi0Param][7];               // parameters for pi0 in 0.9 TeV

  ClassDef(AliGenEMlib,2)

};


#endif
