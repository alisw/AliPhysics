/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$*/

//
// Stores user defined cerenkov properties of media like 
// absorption coefficient, refraction index and quantum efficiency. 
// The properties are stored in arrays. The array index corresponds to discrete 
// optical photon energies defined in fEnergy.
// Access to the properties is provided by interpolation.
// 
// Author:
// A. Morsch 
// andreas.morsch@cern.ch
//

#include "TFlukaCerenkov.h"

Double_t TFlukaCerenkov::fgGlobalMaximumEfficiency = 0.;
   
ClassImp(TFlukaCerenkov)


TFlukaCerenkov::TFlukaCerenkov()
    : fSamples(0), 
      fIsMetal(kFALSE),
      fIsSensitive(kFALSE),
      fEnergy(0),
      fWaveLength(0),
      fAbsorptionCoefficient(0),
      fQuantumEfficiency(0),
      fRefractionIndex(0),
      fReflectivity(0),
      fMaximumEfficiency(0.)
{
// Default constructor
}


TFlukaCerenkov::TFlukaCerenkov(Int_t npckov, Float_t *ppckov, Float_t *absco, Float_t *effic, Float_t *rindex)
    : fSamples(npckov),
      fIsMetal(kFALSE),
      fIsSensitive(kFALSE),
      fEnergy(new Float_t[fSamples]),
      fWaveLength(new Float_t[fSamples]),
      fAbsorptionCoefficient(new Float_t[fSamples]),
      fQuantumEfficiency(new Float_t[fSamples]),
      fRefractionIndex(new Float_t[fSamples]),
      fReflectivity(new Float_t[fSamples]),
      fMaximumEfficiency(0.)
{
// Constructor

//    fSamples = npckov;
//    fEnergy                = new Float_t[fSamples];
//    fWaveLength            = new Float_t[fSamples];
//    fAbsorptionCoefficient = new Float_t[fSamples];
//    fRefractionIndex       = new Float_t[fSamples];
//    fQuantumEfficiency     = new Float_t[fSamples];
//    fReflectivity          = new Float_t[fSamples];    

    for (Int_t i = 0; i < fSamples; i++) {
        fEnergy[i]             = ppckov[i];
        fWaveLength[i]         = khc / ppckov[i];
        if (absco[i] > 0.) {
            fAbsorptionCoefficient[i]   = 1./ absco[i];
        } else {
            fAbsorptionCoefficient[i]   = 1.e15;
        }
        fRefractionIndex[i]    = rindex[i];
        fQuantumEfficiency[i]  = effic[i];
        //
        // Find local maximum quantum efficiency
        if (effic[i] > fMaximumEfficiency) fMaximumEfficiency = effic[i];
        //
        // Flag is sensitive if quantum efficiency 0 < eff < 1 for at least one value.
        if (effic[i] < 1. && effic[i] > 0.) fIsSensitive = 1;
        // G3 way to define metal
        if (rindex[0] == 0.) {
            fIsMetal = kTRUE;
            fReflectivity[i] = absco[i];
        }
    }
    // Find global  maximum quantum efficiency
    if (fMaximumEfficiency > GetGlobalMaximumEfficiency()) {
        SetGlobalMaximumEfficiency(fMaximumEfficiency);
    }
}

TFlukaCerenkov::TFlukaCerenkov(Int_t npckov, Float_t *ppckov, Float_t *absco, Float_t *effic, Float_t *rindex, Float_t *refl)
    : fSamples(npckov),
      fIsMetal(kFALSE),
      fIsSensitive(kFALSE),
      fEnergy(new Float_t[fSamples]),
      fWaveLength(new Float_t[fSamples]),
      fAbsorptionCoefficient(new Float_t[fSamples]),
      fQuantumEfficiency(new Float_t[fSamples]),
      fRefractionIndex(new Float_t[fSamples]),
      fReflectivity(new Float_t[fSamples]),
      fMaximumEfficiency(0.)
{
    // Constructor including reflectivity
//    fSamples = npckov;
//    fEnergy                = new Float_t[fSamples];
//    fWaveLength            = new Float_t[fSamples];
//    fAbsorptionCoefficient = new Float_t[fSamples];
//    fRefractionIndex       = new Float_t[fSamples];
//    fQuantumEfficiency     = new Float_t[fSamples];
    
    
    for (Int_t i = 0; i < fSamples; i++) {
        fEnergy[i]             = ppckov[i];
        fWaveLength[i]         = khc / ppckov[i];
        if (absco[i] > 0.) {
            fAbsorptionCoefficient[i]   = 1./ absco[i];
        } else {
            fAbsorptionCoefficient[i]   = 1.e15;
        }
        fRefractionIndex[i]    = rindex[i];
        fQuantumEfficiency[i]  = effic[i];
        //
        // Find local maximum quantum efficiency
        if (effic[i] > fMaximumEfficiency) fMaximumEfficiency = effic[i];
        //
        // Flag is sensitive if quantum efficiency 0 < eff < 1 for at least one value.
        if (effic[i] < 1. && effic[i] > 0.) fIsSensitive = 1;
        //

    }
    // Find global  maximum quantum efficiency
    if (fMaximumEfficiency > GetGlobalMaximumEfficiency()) {
        SetGlobalMaximumEfficiency(fMaximumEfficiency);
    }
//    fReflectivity     = new Float_t[fSamples];
    for (Int_t i = 0; i < fSamples; i++) fReflectivity[i] = refl[i];
    fIsMetal = kTRUE;
}

Float_t TFlukaCerenkov::GetAbsorptionCoefficient(Float_t energy)
{
//
// Get AbsorptionCoefficient for given energy 
//
    return Interpolate(energy, fEnergy, fAbsorptionCoefficient);
    
}

Float_t TFlukaCerenkov::GetRefractionIndex(Float_t energy)
{
//
// Get RefractionIndex for given energy 
//
    return Interpolate(energy, fEnergy, fRefractionIndex);
    
}

Float_t TFlukaCerenkov::GetReflectivity(Float_t energy)
{
//
// Get RefractionIndex for given energy 
//
    return Interpolate(energy, fEnergy, fReflectivity);
    
}

Float_t TFlukaCerenkov::GetQuantumEfficiency(Float_t energy)
{
//
// Get QuantumEfficiency for given energy 
//
    return Interpolate(energy, fEnergy, fQuantumEfficiency);
    
}


Float_t TFlukaCerenkov::GetAbsorptionCoefficientByWaveLength(Float_t wavelength)
{
//
// Get AbsorptionCoefficient for given wavelength 
//
    Float_t energy = khc / wavelength;    
    return Interpolate(energy, fEnergy, fAbsorptionCoefficient);
    
}

Float_t TFlukaCerenkov::GetRefractionIndexByWaveLength(Float_t wavelength)
{
//
// Get RefractionIndex for given wavelenth 
//
    Float_t energy = khc / wavelength;    
    return Interpolate(energy, fEnergy, fRefractionIndex);
}

Float_t TFlukaCerenkov::GetReflectivityByWaveLength(Float_t wavelength)
{
//
// Get RefractionIndex for given wavelenth 
//
    Float_t energy = khc / wavelength;    
    return Interpolate(energy, fEnergy, fReflectivity);
}

Float_t TFlukaCerenkov::GetQuantumEfficiencyByWaveLength(Float_t wavelength)
{
//
// Get QuantumEfficiency for given wavelength 
//
    Float_t energy = khc / wavelength;    
    return Interpolate(energy, fEnergy, fQuantumEfficiency);
}

Float_t TFlukaCerenkov::Interpolate(Float_t value, Float_t* array1, Float_t* array2)
{
//
// Interpolate array values 
//
    if (value < array1[0] && value >= array1[fSamples - 1]) {
        Warning("Interpolate()", "Photon energy out of range. Returning 0.");
        return (0.);
    }
    
    Int_t i = TMath::BinarySearch(fSamples, array1, value);
    if (i == (fSamples-1)) {
        return (array2[i]);
    } else {
        return (array2[i] + (array2[i+1] - array2[i]) / (array1[i+1] - array1[i]) * (value - array1[i]));
    }
}

