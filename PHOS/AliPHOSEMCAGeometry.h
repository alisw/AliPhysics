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

#include <assert.h> 

// --- ROOT system ---
#include "TObject.h"
class TObjArray ;

// --- AliRoot header files ---

class AliPHOSEMCAGeometry : public TObject {

public: 

  AliPHOSEMCAGeometry();
  AliPHOSEMCAGeometry(const AliPHOSEMCAGeometry & cpv) : TObject(cpv) {
    // cpy ctor requested by Coding Convention but not yet needed
    assert(0==1) ;
  } 
  virtual ~AliPHOSEMCAGeometry(void) {}

  AliPHOSEMCAGeometry & operator = (const AliPHOSEMCAGeometry  & /*rvalue*/) {
    // assignement operator requested by coding convention but not needed
    assert(0==1) ;
    return *this ; 
  }

  Float_t * GetStripHalfSize()   {return fStripHalfSize ;}
  Float_t   GetStripWallWidthOut() {return fStripWallWidthOut ;}
  Float_t * GetAirCellHalfSize() {return fAirCellHalfSize ;}
  Float_t * GetWrappedHalfSize() {return fWrappedHalfSize ;}
  Float_t   GetAirGapLed() const {return fAirGapLed ;}
  Float_t * GetCrystalHalfSize() {return fCrystalHalfSize ;}
  Float_t * GetSupportPlateHalfSize() { return fSupportPlateHalfSize ;}
  Float_t * GetSupportPlateInHalfSize() {return fSupportPlateInHalfSize ;}
  Float_t   GetSupportPlateThickness(void)   const { return fSupportPlateThickness ; }    

  Float_t * GetPreampHalfSize() {return fPreampHalfSize ;}
  Float_t * GetAPDHalfSize(void) {return fPinDiodeHalfSize ; }
  Float_t * GetOuterThermoParams(void) {return  fOuterThermoParams ; }
  Float_t * GetCoolerHalfSize(void) {return fCoolerHalfSize ;}
  Float_t * GetAirGapHalfSize(void) {return fAirGapHalfSize; }
  Float_t * GetInnerThermoHalfSize(void) {return  fInnerThermoHalfSize ; }
  Float_t * GetAlCoverParams() {return fAlCoverParams ; }
  Float_t * GetFiberGlassHalfSize() {return fFiberGlassHalfSize ; }
  Float_t * GetWarmAlCoverHalfSize() {return fWarmAlCoverHalfSize ;}
  Float_t * GetWarmThermoHalfSize() {return fWarmThermoHalfSize ;}
  Float_t * GetTSupport1HalfSize() {return fTSupport1HalfSize ;}
  Float_t * GetTSupport2HalfSize() {return fTSupport2HalfSize ;}
  Float_t * GetTCables1HalfSize() {return fTCables1HalfSize ; }
  Float_t * GetTCables2HalfSize() {return fTCables2HalfSize ; }
  Float_t   GetTSupportDist() {return fTSupportDist ; }
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

  
  Int_t     GetNCellsInStrip() const { return fNCellsInStrip;}
  Int_t     GetNStripX()       const { return fNStripX ; }
  Int_t     GetNStripZ()       const { return fNStripZ ; }
  Int_t     GetNTSuppots()     const { return fNTSupports; }
  Int_t     GetNPhi(void)      const { return fNPhi ; }
  Int_t     GetNZ(void)        const { return fNZ ; }

 
private:

  Float_t fStripHalfSize[3]   ;
  Float_t fAirCellHalfSize[3] ;
  Float_t fWrappedHalfSize[3] ;
  Float_t fSupportPlateHalfSize[3] ;
  Float_t fSupportPlateInHalfSize[3] ;
  Float_t fCrystalHalfSize[3] ;
  Float_t fAirGapLed ;
  Float_t fStripWallWidthOut ; // Side to another strip  
  Float_t fStripWallWidthIn ; 
  Float_t fTyvecThickness ; 
  Float_t fTSupport1HalfSize[3] ;
  Float_t fTSupport2HalfSize[3] ;
  Float_t fPreampHalfSize[3] ;
  Float_t fPinDiodeHalfSize[3] ;              // Size of the PIN Diode 

  Float_t fOuterThermoParams[4] ; 
  Float_t fCoolerHalfSize[3] ; 
  Float_t fAirGapHalfSize[3] ;
  Float_t fInnerThermoHalfSize[3] ;
  Float_t fAlCoverParams[4] ;
  Float_t fFiberGlassHalfSize[3] ;


  Float_t fInnerThermoWidthX ;
  Float_t fInnerThermoWidthY ;
  Float_t fInnerThermoWidthZ ;
  Float_t fAirGapWidthX ;
  Float_t fAirGapWidthY ;
  Float_t fAirGapWidthZ ;
  Float_t fCoolerWidthX ;
  Float_t fCoolerWidthY ;
  Float_t fCoolerWidthZ ;
  Float_t fAlCoverThickness ;
  Float_t fOuterThermoWidthXUp ;
  Float_t fOuterThermoWidthXLow;
  Float_t fOuterThermoWidthY ;
  Float_t fOuterThermoWidthZ ;
  Float_t fAlFrontCoverX  ;
  Float_t fAlFrontCoverZ  ; 
  Float_t fFiberGlassSup2X ;
  Float_t fFiberGlassSup1X ; 
  Float_t fFrameHeight ; 
  Float_t fFrameThickness ; 
  Float_t fAirSpaceFeeX ;
  Float_t fAirSpaceFeeZ ;
  Float_t fAirSpaceFeeY ;
  Float_t fTCables2HalfSize[3] ; 
  Float_t fTCables1HalfSize[3] ; 
  Float_t fWarmUpperThickness ;
  Float_t fWarmBottomThickness ;
  Float_t fWarmAlCoverWidthX ;
  Float_t fWarmAlCoverWidthY ;
  Float_t fWarmAlCoverWidthZ ;
  Float_t fWarmAlCoverHalfSize[3] ;
  Float_t fWarmThermoHalfSize[3] ;
  Float_t fFiberGlassSup1Y ;
  Float_t fFiberGlassSup2Y ;
  Float_t fTSupportDist ;
  Float_t fTSupport1Thickness ;
  Float_t fTSupport2Thickness ;
  Float_t fTSupport1Width ;
  Float_t fTSupport2Width ;
  Float_t fFrameXHalfSize[3] ;
  Float_t fFrameZHalfSize[3] ;
  Float_t fFrameXPosition[3] ;
  Float_t fFrameZPosition[3] ;
  Float_t fFGupXHalfSize[3] ;
  Float_t fFGupXPosition[3] ;
  Float_t fFGupZHalfSize[3] ;
  Float_t fFGupZPosition[3] ;
  Float_t fFGlowXHalfSize[3] ;
  Float_t fFGlowXPosition[3] ;
  Float_t fFGlowZHalfSize[3] ;
  Float_t fFGlowZPosition[3] ;
  Float_t fFEEAirHalfSize[3] ;
  Float_t fFEEAirPosition[3] ;
  Float_t fEMCParams[4] ;
  Float_t fIPtoOuterCoverDistance ;       // Distances from interaction point to outer cover 
  Float_t fIPtoCrystalSurface ;           // Distances from interaction point to Xtal surface

  Float_t fSupportPlateThickness ;        // Thickness of the Aluminium support plate for Strip   

  Int_t  fNCellsInStrip ;
  Int_t  fNStripX ;
  Int_t  fNStripZ ;
  Int_t  fNTSupports ;
  Int_t  fNPhi ;                         // Number of crystal units in X (phi) direction
  Int_t  fNZ ;                           // Number of crystal units in Z direction

  ClassDef(AliPHOSEMCAGeometry,1)         // EMCA geometry class 

} ;

#endif // AliPHOSEMCAGEOMETRY_H
