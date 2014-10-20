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
#include "AliPythia.h"
#include "AliDecayer.h"
#include "AliGenParam.h"


class AliGenCocktailEntry;

class AliGenEMCocktail : public AliGenCocktail
{
public:

    AliGenEMCocktail();
  enum GeneratorIndex_t { kPizero=0, kEta, kRho, kOmega, kEtaprime, kPhi, kJpsi, kDirectRealGamma, kDirectVirtGamma, kGENs };
  enum ParticeGenerator_t { kGenPizero=0x001, kGenEta=0x002, kGenRho=0x004, kGenOmega=0x008, kGenEtaprime=0x010, kGenPhi=0x020, kGenJpsi=0x040, kGenDirectRealGamma=0x100, kGenDirectVirtGamma=0x200, kGenHadrons=0x7f, kGenGammas=0x300 };

    virtual ~AliGenEMCocktail();    
    virtual void Init();
    virtual void CreateCocktail();
    virtual void Generate();    
    Float_t GetDecayMode()         const {return fDecayMode;}
    Float_t GetWeightingMode()     const {return fWeightingMode;}
    void    SetDecayer(AliDecayer* const decayer){fDecayer = decayer;}
    void    SetDecayMode(Decay_t decay){ fDecayMode = decay;}
    void    SetWeightingMode(Weighting_t weight){ fWeightingMode = weight;}
  void    SetNPart(Int_t npart){ fNPart = npart; }
  void    SetPtParam(AliGenEMlib::PtParamSet_t PtSelect){ fPtSelect = PtSelect; }
  void    SetCentrality(AliGenEMlib::Centrality_t cent){ fCentrality = cent; }
  void    SetV2Systematic(Int_t v2sys){ fV2Systematic = v2sys; }
  void    SetForceGammaConversion(Bool_t force=kTRUE){ fForceConv=force; }
  void    SetHeaviestHadron(ParticeGenerator_t part);
  void    SelectMotherParticles(ParticeGenerator_t part){ fSelectedParticles=part; }
    
private:
    AliGenEMCocktail(const AliGenEMCocktail &cocktail); 
    AliGenEMCocktail & operator=(const AliGenEMCocktail &cocktail); 

    void AddSource2Generator(Char_t *nameReso, AliGenParam* const genReso);
    
    AliDecayer* fDecayer;        // External decayer
  Decay_t fDecayMode;   //decay mode in which resonances are forced to decay, default: kAll
  Weighting_t fWeightingMode; //weighting mode: kAnalog or kNonAnalog
    
    Int_t    fNPart;             // multiplicity of each source per event
    Double_t fYieldArray[kGENs]; // array of dN/dy for each source

  AliGenEMlib::PtParamSet_t fPtSelect; // selected pT parameterization
  AliGenEMlib::Centrality_t fCentrality; // selected centrality
  Int_t fV2Systematic; //selected systematic error for v2 parameters

  Bool_t fForceConv; //select whether you want to force all gammas to convert imidediately
  Int_t fSelectedParticles; //which particles to simulate

    ClassDef(AliGenEMCocktail,1)  //  cocktail for EM physics
};

#endif



