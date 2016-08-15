#ifndef ALIGENEMCOCKTAIL_H
#define ALIGENEMCOCKTAIL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
// Class to create the cocktail for physics with electrons, di-electrons,
// and photons from the decay of the following sources:
// pizero, eta, rho, omega, etaprime, phi
//

#include "AliGenCocktail.h"
#include "AliGenEMlib.h"
#include "AliDecayer.h"
#include "AliGenParam.h"


class AliGenCocktailEntry;

class AliGenEMCocktail : public AliGenCocktail
{
public:

  AliGenEMCocktail();
  enum GeneratorIndex_t { kPizero=0, kEta, kRho0, kOmega, kEtaprime, kPhi, kJpsi, kSigma0, kK0s, kDeltaPlPl, kDeltaPl, kDeltaMi, kDeltaZero, kRhoPl, kRhoMi, kK0star, kDirectRealGamma, kDirectVirtGamma, kGENs };
  enum ParticleGenerator_t { kGenPizero = 0x00001, kGenEta = 0x00002, kGenRho0 = 0x00004, kGenOmega = 0x00008, 
			     kGenEtaprime = 0x00010, kGenPhi = 0x00020, kGenJpsi = 0x00040, 
			     kGenDirectRealGamma = 0x00100, kGenDirectVirtGamma = 0x00200, kGenSigma0 = 0x00400, kGenK0s = 0x00800,
			     kGenDeltaPlPl = 0x01000, kGenDeltaPl = 0x02000, kGenDeltaMi = 0x04000, kGenDeltaZero = 0x08000, 
			     kGenRhoPl = 0x10000, kGenRhoMi = 0x20000, kGenK0star = 0x40000,  
			     kGenHadrons = 0x0007f, kGenGammas = 0x00300 };

  virtual ~AliGenEMCocktail();    
  virtual void Init();
  virtual void CreateCocktail();
  virtual void Generate();    
  Float_t GetDecayMode()                 const  { return fDecayMode   ;}
  Float_t GetWeightingMode()             const  { return fWeightingMode  ;}
  void    SetDecayer(AliDecayer* const decayer)      { fDecayer = decayer  ;}
  void    SetDecayMode(Decay_t decay)         { fDecayMode = decay  ;}
  void    SetWeightingMode(Weighting_t weight)      { fWeightingMode = weight ;}
  void    SetNPart(Int_t npart){ fNPart = npart; }
  void    SetPtParamPi0(AliGenEMlib::PtParamSetPi0_t PtSelect)  { fPtSelectPi0 = PtSelect ;}
  void    SetPtParamEta(AliGenEMlib::PtParamSetEta_t PtSelect)  { fPtSelectEta = PtSelect ;}
  void    SetPtParamOmega(AliGenEMlib::PtParamSetOmega_t PtSelect) { fPtSelectOmega = PtSelect ;}
  void    SetPtParamPhi(AliGenEMlib::PtParamSetPhi_t PtSelect)  { fPtSelectPhi = PtSelect ;}
  void    SetCollisionSystem(AliGenEMlib::CollisionSystem_t col)  { fCollisionSystem = col ;}
  void    SetCentrality(AliGenEMlib::Centrality_t cent)    { fCentrality = cent  ;}
  void    SetV2Systematic(AliGenEMlib::v2Sys_t v2sys)     { fV2Systematic = v2sys  ;}
  void    SetForceGammaConversion(Bool_t force=kTRUE)     { fForceConv=force   ;}
  void    SetHeaviestHadron(ParticleGenerator_t part);
    
  //***********************************************************************************************
  // This function allows to select the particle which should be procude based on 1 Integer value
  // this integer value is then bitwise compare to the values in SelectParticle
  // Examples: 
  // a) you would like to switch on: pi0, eta, rho0, omega, eta', phi, jpsi, sigma0 and all deltas
  //    implies you want the binary number: 00 1111 0100 0111 1111 = 
  //   which translates 62591_10 (in decimal) and F47F_16 (in hexadecimal)
  // b) you would like to switch on: pi0, eta, rho0, omega, eta', phi, sigma0 and 
  //    implies you want the binary number: 00 0000 0100 0011 1111 = 
  //   which translates 1087_10 (in decimal) and 43F_16 (in hexadecimal)
  // c) you would like to switch on: pi0, eta, rho0, omega, eta', phi, sigma0 and all deltas
  //    implies you want the binary number: 00 1111 0100 0011 1111 = 
  //   which translates 62527_10 (in decimal) and F43F_16 (in hexadecimal)
  // d) you would like to switch on: pi0, eta, rho0, omega, eta', phi
  //    implies you want the binary number: 00 0000 0000 0011 1111 = 
  //   which translates 63_10 (in decimal) and 3F_16 (in hexadecimal)
  //***********************************************************************************************
  void    SelectMotherParticles(UInt_t part)       { fSelectedParticles=part ;}

private:
  AliGenEMCocktail(const AliGenEMCocktail &cocktail); 
  AliGenEMCocktail & operator=(const AliGenEMCocktail &cocktail); 

  void AddSource2Generator(Char_t *nameReso, AliGenParam* const genReso);
  
  AliDecayer*      fDecayer;        // External decayer
  Decay_t       fDecayMode;    // decay mode in which resonances are forced to decay, default: kAll
  Weighting_t      fWeightingMode;  // weighting mode: kAnalog or kNonAnalog
  
  Int_t          fNPart;             // multiplicity of each source per event
  Double_t       fYieldArray[kGENs]; // array of dN/dy for each source

  AliGenEMlib::CollisionSystem_t  fCollisionSystem;   // selected collision system
  AliGenEMlib::PtParamSetPi0_t  fPtSelectPi0;   // selected pT parameterization for pi0
  AliGenEMlib::PtParamSetEta_t  fPtSelectEta;   // selected pT parameterization for eta
  AliGenEMlib::PtParamSetOmega_t  fPtSelectOmega;  // selected pT parameterization for omega
  AliGenEMlib::PtParamSetPhi_t  fPtSelectPhi;   // selected pT parameterization for phi
  AliGenEMlib::Centrality_t   fCentrality;   // selected centrality
  AliGenEMlib::v2Sys_t    fV2Systematic;   // selected systematic error for v2 parameters

  Bool_t        fForceConv;   // select whether you want to force all gammas to convert imidediately
  UInt_t        fSelectedParticles; // which particles to simulate, allows to switch on and off 32 different particles



  ClassDef(AliGenEMCocktail,2)       // cocktail for EM physics
};

#endif



