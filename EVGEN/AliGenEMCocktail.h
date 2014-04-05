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
  enum GeneratorCode { kGenPizero=0, kGenEta, kGenRho, kGenOmega, kGenEtaprime, kGenPhi, kGenJpsi, kGenPromptRealGamma, kGenThermRealGamma, kGenPromptVirtGamma, kGenThermVirtGamma, kGENs };

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
  void    SetPtParam(Int_t PtSelect){ fPtSelect = PtSelect; }
  void    SetCentrality(Int_t cent){ fCentrality = cent; }
  void    SetV2Systematic(Int_t v2sys){ fV2Systematic = v2sys; }
  void    SetForceGammaConversion(Bool_t force=kTRUE){ fForceConv=force; }
  void    SetHeaviestParticle(Int_t part){ fHeaviestParticle=part; }
    
private:
    AliGenEMCocktail(const AliGenEMCocktail &cocktail); 
    AliGenEMCocktail & operator=(const AliGenEMCocktail &cocktail); 

    void AddSource2Generator(Char_t *nameReso, AliGenParam* const genReso);
    
    AliDecayer* fDecayer;        // External decayer
  Decay_t fDecayMode;   //decay mode in which resonances are forced to decay, default: kAll
  Weighting_t fWeightingMode; //weighting mode: kAnalog or kNonAnalog
    
    Int_t    fNPart;             // multiplicity of each source per event
    Double_t fYieldArray[kGENs]; // array of dN/dy for each source

  Int_t fPtSelect; // selected pT parameterization
  Int_t fCentrality; // selected centrality
  Int_t fV2Systematic; //selected systematic error for v2 parameters

  Bool_t fForceConv; //select whether you want to force all gammas to convert imidediately
  Int_t fHeaviestParticle; //select up to which particle to simulate

    ClassDef(AliGenEMCocktail,1)  //  cocktail for EM physics
};

#endif



