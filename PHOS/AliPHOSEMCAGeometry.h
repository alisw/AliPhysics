#ifndef ALIPHOSEMCAGEOMETRY_H
#define ALIPHOSEMCAGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Geometry class  for PHOS : EMCA (Electromagnetic Calirometer)
//                  
//*-- Author:   Yves Schutz (Subatech)
//    Modified: Yuri Kharlov (IHEP, Protvino)
//    15 September 2000

// --- ROOT system ---

#include "TObjArray.h"

// --- AliRoot header files ---

class AliPHOSEMCAGeometry : public TObject {

public: 

           AliPHOSEMCAGeometry();
  virtual ~AliPHOSEMCAGeometry(void) {}

  Float_t    GetAirFilledBoxSize(Int_t index)     const { 
    return fAirFilledBoxSize[index] ;}
  Float_t    GetCrystalHolderThickness(void)      const { 
    return fCrystalHolderThickness ; } 
  Float_t    GetCrystalSize(Int_t index)          const { 
    return fXtlSize[index] ; }
  Float_t    GetCrystalSupportHeight(void)        const { 
    return fCrystalSupportHeight ; } 
  Float_t    GetCrystalWrapThickness(void)        const { 
    return fCrystalWrapThickness;}
  Float_t    GetGapBetweenCrystals(void)          const { 
    return fGapBetweenCrystals ; }
  Float_t    GetIPtoCrystalSurface(void)          const { 
    return fIPtoCrystalSurface ; }
  Float_t    GetIPtoOuterCoverDistance(void)      const { 
    return fIPtoOuterCoverDistance ; }
  Float_t    GetLowerThermoPlateThickness(void)   const { 
    return fLowerThermoPlateThickness ; }
  Float_t    GetLowerTextolitPlateThickness(void) const { 
    return fLowerTextolitPlateThickness ; }
  Float_t    GetModuleBoxThickness(void)          const { 
    return fModuleBoxThickness ; }
  Int_t      GetNPhi(void)                        const { 
    return fNPhi ; }
  Int_t      GetNZ(void)                          const { 
    return fNZ ; }
  Float_t    GetOuterBoxSize(Int_t index)         const { 
    return fOuterBoxSize[index] ;    }
  Float_t    GetOuterBoxThickness(Int_t index)    const { 
    return fOuterBoxThickness[index] ; } 
  Float_t    GetPinDiodeSize(Int_t index)         const { 
    return fPinDiodeSize[index] ; }
  Float_t    GetSecondUpperPlateThickness(void)   const { 
    return fSecondUpperPlateThickness ; }
  Float_t    GetSupportPlateThickness(void)       const { 
    return fSupportPlateThickness ; }    
  Float_t    GetTextolitBoxSize(Int_t index)      const { 
    return fTextolitBoxSize[index] ; }
  Float_t    GetTextolitBoxThickness(Int_t index) const { 
    return fTextolitBoxThickness[index]; } 
  Float_t    GetUpperPlateThickness(void)         const { 
    return fUpperPlateThickness ; }
  Float_t    GetUpperCoolingPlateThickness(void)  const { 
    return fUpperCoolingPlateThickness ; }
 
private:

  Float_t fAirFilledBoxSize[3] ;          // Air filled box containing one module
  Float_t fAirThickness[3] ;              // Space filled with air between the module box and the Textolit box
  Float_t fCrystalSupportHeight ;         // Height of the support of the crystal    
  Float_t fCrystalWrapThickness ;         // Thickness of Tyvek wrapping the crystal
  Float_t fCrystalHolderThickness ;       // Titanium holder of the crystal
  Float_t fGapBetweenCrystals ;           // Total Gap between two adjacent crystals 
  Float_t fIPtoOuterCoverDistance ;       // Distances from interaction point to outer cover 
  Float_t fIPtoCrystalSurface ;           // Distances from interaction point to Xtal surface
  Float_t fModuleBoxThickness ;           // Thickness of the thermo insulating box containing one crystals module 
  Float_t fLowerTextolitPlateThickness ;  // Thickness of lower textolit plate
  Float_t fLowerThermoPlateThickness ;    // Thickness of lower thermo insulating plate
  Int_t   fNPhi ;                         // Number of crystal units in X (phi) direction
  Int_t   fNZ ;                           // Number of crystal units in Z direction
  Float_t fOuterBoxSize[3] ;              // Size of the outer  thermo insulating foam box
  Float_t fOuterBoxThickness[3] ;         // Thickness of the outer thermo insulating foam box
  Float_t fPinDiodeSize[3] ;              // Size of the PIN Diode 
  TObjArray *  fRotMatrixArray ;          // Liste of rotation matrices (one per phos module)
  Float_t fSecondUpperPlateThickness ;    // Thickness of  upper polystyrene foam plate
  Float_t fSupportPlateThickness ;        // Thickness of the Aluminium support plate  
  Float_t fUpperCoolingPlateThickness ;   // Thickness of the upper cooling plate 
  Float_t fUpperPlateThickness ;          // Thickness of the uper thermo insulating foam plate 
  Float_t fTextolitBoxSize[3] ;           // Size of the Textolit box inside the insulating foam box
  Float_t fTextolitBoxThickness[3] ;      // Thicknesses of th Textolit box
  Float_t fXtlSize[3] ;                   // PWO4 crystal dimensions

  ClassDef(AliPHOSEMCAGeometry,1)         // EMCA geometry class 

} ;

#endif // AliPHOSEMCAGEOMETRY_H
