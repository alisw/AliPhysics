#ifndef TFLUKACERENKOV
#define TFLUKACERENKOV

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
// Class that gives access to optical properties for Cerenkov photon         //
// production and transport                                                  //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TObject.h>
#include <TMath.h>

const Double_t khc = 2. * TMath::Pi() * 0.1973269602e-13; // GeV cm 


class TFlukaCerenkov : public TObject
{

public:
   // constructors
    TFlukaCerenkov();
    TFlukaCerenkov(Int_t npckov, Float_t *ppckov, Float_t *absco, Float_t *effic, Float_t *rindex);
    TFlukaCerenkov(Int_t npckov, Float_t *ppckov, Float_t *absco, Float_t *effic, Float_t *rindex, Float_t* rfl);
    virtual Float_t   GetAbsorptionCoefficient(Float_t energy);
    virtual Float_t   GetQuantumEfficiency(Float_t energy);
    virtual Float_t   GetRefractionIndex(Float_t energy);
    virtual Float_t   GetReflectivity(Float_t energy);    
    virtual Float_t   GetAbsorptionCoefficientByWaveLength(Float_t energy);
    virtual Float_t   GetQuantumEfficiencyByWaveLength(Float_t energy);
    virtual Float_t   GetRefractionIndexByWaveLength(Float_t energy);
    virtual Float_t   GetReflectivityByWaveLength(Float_t energy);
    virtual Float_t   GetMinimumEnergy()     {return fEnergy[0];}
    virtual Float_t   GetMaximumEnergy()     {return fEnergy[fSamples-1];}
    virtual Float_t   GetMinimumWavelength() {return khc / fEnergy[fSamples-1];}
    virtual Float_t   GetMaximumWavelength() {return khc / fEnergy[0];}
    virtual Int_t     GetNSamples()          {return fSamples;}
    virtual Bool_t    IsMetal()              {return fIsMetal;}
    virtual Bool_t    IsSensitive()          {return fIsSensitive;}
    virtual Double_t  GetMaximumEfficiency() const             {return fMaximumEfficiency;}
    static  Double_t  GetGlobalMaximumEfficiency()             {return fgGlobalMaximumEfficiency;}
    static  void      SetGlobalMaximumEfficiency(Double_t eff) {fgGlobalMaximumEfficiency = eff;}

 protected:

    virtual Float_t  Interpolate(Float_t energy, Float_t* array1, Float_t* array2);

 protected:
    Int_t        fSamples;                  // Number of sampling points
    Bool_t       fIsMetal;                  // Flag for metals
    Bool_t       fIsSensitive;              // Flag for metals  
    Float_t*     fEnergy;                   // [fSamples] Energy                 (GeV) 
    Float_t*     fWaveLength;               // [fSamples] Wafelength             (cm)
    Float_t*     fAbsorptionCoefficient;    // [fSamples] Absorption Coefficient (1/cm)
    Float_t*     fQuantumEfficiency;        // [fSamples] Quantum efficiency
    Float_t*     fRefractionIndex;          // [fSamples] Refraction Index
    Float_t*     fReflectivity;             // [fSamples] Reflectivity
    Double_t     fMaximumEfficiency;        // Local maximum quantum efficiency
    // static 
    static Double_t fgGlobalMaximumEfficiency; // Global maximum quantum efficiency

 private:
    // Copy constructor and operator= declared but not implemented (-Weff++ flag)
    TFlukaCerenkov(const TFlukaCerenkov&);
    TFlukaCerenkov& operator=(const TFlukaCerenkov&);

    ClassDef(TFlukaCerenkov, 1)          // CerenkovProperties
};
	
#endif
	
