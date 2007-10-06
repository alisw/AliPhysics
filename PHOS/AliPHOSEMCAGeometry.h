#ifndef ALIPHOSEMCAGEOMETRY_H
#define ALIPHOSEMCAGEOMETRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Geometry class  for PHOS : EMCA (Electromagnetic Calirometer)
// Its data members provide geometry parametrization of EMCA
// which can be changed in the constructor only.
// Author:   Yves Schutz (Subatech)
// Modified: Yuri Kharlov (IHEP, Protvino)
// 15 September 2000

// --- ROOT system ---
#include "TObject.h"
class TObjArray ;

// --- AliRoot header files ---

class AliPHOSEMCAGeometry : public TObject {

public: 

  AliPHOSEMCAGeometry();
  //Compiler-generated copy ctor and copy-assignment operator _ARE_ OK.

  virtual ~AliPHOSEMCAGeometry(void) {}

  Float_t * GetStripHalfSize() {return fStripHalfSize ;}
  Float_t   GetStripWallWidthOut() const {return fStripWallWidthOut ;}
  Float_t * GetAirCellHalfSize() {return fAirCellHalfSize ;}
  Float_t * GetWrappedHalfSize() {return fWrappedHalfSize ;}
  Float_t   GetAirGapLed() const {return fAirGapLed ;}
  Float_t * GetCrystalHalfSize() {return fCrystalHalfSize ;}
  Float_t * GetSupportPlateHalfSize() { return fSupportPlateHalfSize ;}
  Float_t * GetSupportPlateInHalfSize()  {return fSupportPlateInHalfSize ;}
  Float_t   GetSupportPlateThickness(void)   const { return fSupportPlateThickness ; }    

  Float_t * GetPreampHalfSize() {return fPreampHalfSize ;}
  Float_t * GetAPDHalfSize(void) {return fPinDiodeHalfSize ; }
  Float_t * GetOuterThermoParams(void) {return  fOuterThermoParams ; }
  Float_t * GetCoolerHalfSize(void)  {return fCoolerHalfSize ;}
  Float_t * GetAirGapHalfSize(void)  {return fAirGapHalfSize; }
  Float_t * GetInnerThermoHalfSize(void) {return  fInnerThermoHalfSize ; }
  Float_t * GetAlCoverParams() {return fAlCoverParams ; }
  Float_t * GetFiberGlassHalfSize() {return fFiberGlassHalfSize ; }
  Float_t * GetWarmAlCoverHalfSize() {return fWarmAlCoverHalfSize ;}
  Float_t * GetWarmThermoHalfSize() {return fWarmThermoHalfSize ;}
  Float_t * GetTSupport1HalfSize() {return fTSupport1HalfSize ;}
  Float_t * GetTSupport2HalfSize() {return fTSupport2HalfSize ;}
  Float_t * GetTCables1HalfSize() {return fTCables1HalfSize ; }
  Float_t * GetTCables2HalfSize() {return fTCables2HalfSize ; }
  Float_t   GetTSupportDist() const {return fTSupportDist ; }
  Float_t * GetFrameXHalfSize() {return fFrameXHalfSize ;}
  Float_t * GetFrameZHalfSize() {return fFrameZHalfSize ;}
  Float_t * GetFrameXPosition() {return fFrameXPosition ;}
  Float_t * GetFrameZPosition() {return fFrameZPosition ;}
  Float_t * GetFGupXHalfSize()  {return fFGupXHalfSize ; }
  Float_t * GetFGupXPosition()  {return fFGupXPosition ; }
  Float_t * GetFGupZHalfSize()  {return fFGupZHalfSize ; }
  Float_t * GetFGupZPosition()  {return fFGupZPosition ; }
  Float_t * GetFGlowXHalfSize()  {return fFGlowXHalfSize ; }
  Float_t * GetFGlowXPosition()  {return fFGlowXPosition ; }
  Float_t * GetFGlowZHalfSize()  {return fFGlowZHalfSize ; }
  Float_t * GetFGlowZPosition()  {return fFGlowZPosition ; }
  Float_t * GetFEEAirHalfSize()  {return fFEEAirHalfSize ; }
  Float_t * GetFEEAirPosition()  {return fFEEAirPosition ; }
  Float_t * GetEMCParams() {return fEMCParams ;}

  Float_t    GetIPtoCrystalSurface(void)          const { 
    return fIPtoCrystalSurface ; }
  Float_t    GetIPtoOuterCoverDistance(void)      const { 
    return fIPtoOuterCoverDistance ; }
  Float_t    GetCrystalSize(Int_t index)  const {return 2.*fCrystalHalfSize[index] ; }

  Int_t     GetNCellsXInStrip() const { return fNCellsXInStrip;}
  Int_t     GetNCellsZInStrip() const { return fNCellsZInStrip;}
  Int_t     GetNStripX()        const { return fNStripX ; }
  Int_t     GetNStripZ()        const { return fNStripZ ; }
  Int_t     GetNTSuppots()      const { return fNTSupports; }
  Int_t     GetNPhi()           const { return fNPhi ; }
  Int_t     GetNZ()             const { return fNZ ; }
 
private:

  Float_t fStripHalfSize[3]   ;        // Strip unit size/2
  Float_t fAirCellHalfSize[3] ;        // geometry parameter
  Float_t fWrappedHalfSize[3] ;        // geometry parameter
  Float_t fSupportPlateHalfSize[3] ;   // geometry parameter
  Float_t fSupportPlateInHalfSize[3] ; // geometry parameter
  Float_t fCrystalHalfSize[3] ;        // crystal size/2
  Float_t fAirGapLed ;                 // geometry parameter
  Float_t fStripWallWidthOut ;         // Side to another strip  
  Float_t fStripWallWidthIn ;          // geometry parameter
  Float_t fTyvecThickness ;            // geometry parameter
  Float_t fTSupport1HalfSize[3] ;      // geometry parameter
  Float_t fTSupport2HalfSize[3] ;      // geometry parameter
  Float_t fPreampHalfSize[3] ;         // geometry parameter
  Float_t fPinDiodeHalfSize[3] ;       // Size of the PIN Diode 

  Float_t fOuterThermoParams[4] ;      // geometry parameter
  Float_t fCoolerHalfSize[3] ;         // geometry parameter
  Float_t fAirGapHalfSize[3] ;         // geometry parameter
  Float_t fInnerThermoHalfSize[3] ;    // geometry parameter
  Float_t fAlCoverParams[4] ;          // geometry parameter
  Float_t fFiberGlassHalfSize[3] ;     // geometry parameter


  Float_t fInnerThermoWidthX ;         // geometry parameter
  Float_t fInnerThermoWidthY ;         // geometry parameter
  Float_t fInnerThermoWidthZ ;         // geometry parameter
  Float_t fAirGapWidthX ;              // geometry parameter
  Float_t fAirGapWidthY ;              // geometry parameter
  Float_t fAirGapWidthZ ;              // geometry parameter
  Float_t fCoolerWidthX ;              // geometry parameter
  Float_t fCoolerWidthY ;              // geometry parameter
  Float_t fCoolerWidthZ ;              // geometry parameter
  Float_t fAlCoverThickness ;          // geometry parameter
  Float_t fOuterThermoWidthXUp ;       // geometry parameter
  Float_t fOuterThermoWidthXLow;       // geometry parameter
  Float_t fOuterThermoWidthY ;         // geometry parameter
  Float_t fOuterThermoWidthZ ;         // geometry parameter
  Float_t fAlFrontCoverX  ;            // geometry parameter
  Float_t fAlFrontCoverZ  ;            // geometry parameter
  Float_t fFiberGlassSup2X ;           // geometry parameter
  Float_t fFiberGlassSup1X ;           // geometry parameter
  Float_t fFrameHeight ;               // geometry parameter
  Float_t fFrameThickness ;            // geometry parameter
  Float_t fAirSpaceFeeX ;              // geometry parameter
  Float_t fAirSpaceFeeZ ;              // geometry parameter
  Float_t fAirSpaceFeeY ;              // geometry parameter
  Float_t fTCables2HalfSize[3] ;       // geometry parameter
  Float_t fTCables1HalfSize[3] ;       // geometry parameter
  Float_t fWarmUpperThickness ;        // geometry parameter
  Float_t fWarmBottomThickness ;       // geometry parameter
  Float_t fWarmAlCoverWidthX ;         // geometry parameter
  Float_t fWarmAlCoverWidthY ;         // geometry parameter
  Float_t fWarmAlCoverWidthZ ;         // geometry parameter
  Float_t fWarmAlCoverHalfSize[3] ;    // geometry parameter
  Float_t fWarmThermoHalfSize[3] ;     // geometry parameter
  Float_t fFiberGlassSup1Y ;           // geometry parameter
  Float_t fFiberGlassSup2Y ;           // geometry parameter
  Float_t fTSupportDist ;              // geometry parameter
  Float_t fTSupport1Thickness ;        // geometry parameter
  Float_t fTSupport2Thickness ;        // geometry parameter
  Float_t fTSupport1Width ;            // geometry parameter
  Float_t fTSupport2Width ;            // geometry parameter
  Float_t fFrameXHalfSize[3] ;         // geometry parameter
  Float_t fFrameZHalfSize[3] ;         // geometry parameter
  Float_t fFrameXPosition[3] ;         // geometry parameter
  Float_t fFrameZPosition[3] ;         // geometry parameter
  Float_t fFGupXHalfSize[3] ;          // geometry parameter
  Float_t fFGupXPosition[3] ;          // geometry parameter
  Float_t fFGupZHalfSize[3] ;          // geometry parameter
  Float_t fFGupZPosition[3] ;          // geometry parameter
  Float_t fFGlowXHalfSize[3] ;         // geometry parameter
  Float_t fFGlowXPosition[3] ;         // geometry parameter
  Float_t fFGlowZHalfSize[3] ;         // geometry parameter
  Float_t fFGlowZPosition[3] ;         // geometry parameter
  Float_t fFEEAirHalfSize[3] ;         // geometry parameter
  Float_t fFEEAirPosition[3] ;         // geometry parameter
  Float_t fEMCParams[4] ;              // geometry parameter
  Float_t fIPtoOuterCoverDistance ;       // Distances from interaction point to outer cover 
  Float_t fIPtoCrystalSurface ;           // Distances from interaction point to Xtal surface

  Float_t fSupportPlateThickness ;        // Thickness of the Aluminium support plate for Strip   

  Int_t  fNCellsXInStrip ;              // Number of cells in a strip unit in X
  Int_t  fNCellsZInStrip ;              // Number of cells in a strip unit in Z
  Int_t  fNStripX ;                    // Number of strip units in X
  Int_t  fNStripZ ;                    // Number of strip units in Z
  Int_t  fNTSupports ;                 // geometry parameter
  Int_t  fNPhi ;                       // Number of crystal units in X (phi) direction
  Int_t  fNZ ;                         // Number of crystal units in Z direction
  ClassDef(AliPHOSEMCAGeometry,1)         // EMCA geometry class 

} ;

#endif // AliPHOSEMCAGEOMETRY_H
