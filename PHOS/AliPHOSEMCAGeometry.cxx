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

/* $Id$ */

//_________________________________________________________________________
// Geometry class  for PHOS : EMCA (Electromagnetic Calorimeter)  
// Its data members provide geometry parametrization of EMCA
// which can be changed in the constructor only.
// Author   : Yves Schutz (SUBATECH)
// Modified : Yuri Kharlov (IHEP, Protvino)
// 13 September 2000
// Modified : Dmitri Peressounko (RRC "Kurchatov Institute")
// 6 August 2001

// --- AliRoot header files ---

#include "AliPHOSEMCAGeometry.h"

ClassImp(AliPHOSEMCAGeometry)

//____________________________________________________________________________
AliPHOSEMCAGeometry::AliPHOSEMCAGeometry():
                     fAirGapLed(0.f),
                     fStripWallWidthOut(0.f),
                     fStripWallWidthIn(0.f),
                     fTyvecThickness(0.f),
                     fInnerThermoWidthX(0.f),
                     fInnerThermoWidthY(0.f),
                     fInnerThermoWidthZ(0.f),
                     fAirGapWidthX(0.f),
                     fAirGapWidthY(0.f),
                     fAirGapWidthZ(0.f),
                     fCoolerWidthX(0.f),
                     fCoolerWidthY(0.f),
                     fCoolerWidthZ(0.f),
                     fAlCoverThickness(0.f),
                     fOuterThermoWidthXUp(0.f),
                     fOuterThermoWidthXLow(0.f),
                     fOuterThermoWidthY(0.f),
                     fOuterThermoWidthZ(0.f),
                     fAlFrontCoverX(0.f),
                     fAlFrontCoverZ(0.f),
                     fFiberGlassSup2X(0.f),
                     fFiberGlassSup1X(0.f),
                     fFrameHeight(0.f),
                     fFrameThickness(0.f),
                     fAirSpaceFeeX(0.f),
                     fAirSpaceFeeZ(0.f),
                     fAirSpaceFeeY(0.f),
                     fWarmUpperThickness(0.f),
                     fWarmBottomThickness(0.f),
                     fWarmAlCoverWidthX(0.f),
                     fWarmAlCoverWidthY(0.f),
                     fWarmAlCoverWidthZ(0.f),
                     fFiberGlassSup1Y(0.f),
                     fFiberGlassSup2Y(0.f),
                     fTSupportDist(0.f),
                     fTSupport1Thickness(0.f),
                     fTSupport2Thickness(0.f),
                     fTSupport1Width(0.f),
                     fTSupport2Width(0.f),
                     fIPtoOuterCoverDistance(0.f),
                     fIPtoCrystalSurface(0.f),
                     fSupportPlateThickness(0.f),
                     fNCellsXInStrip(0),
                     fNCellsZInStrip(0),
                     fNStripX(0),
                     fNStripZ(0),
                     fNTSupports(0),
                     fNPhi(0),
                     fNZ(0)
{


  // Initializes the EMC parameters
  // Coordinate system chosen: x across beam, z along beam, y out of beam.
  // Reference point for all volumes incide module is 
  // center of module in x,z on the upper surface of support beam

  //Distance from IP to surface of the crystals
  fIPtoCrystalSurface     = 460.0 ;    


  //CRYSTAL

  fCrystalHalfSize[0] =  2.2 /2 ;  //Half-Sizes of crystall
  fCrystalHalfSize[1] = 18.0 /2 ;
  fCrystalHalfSize[2] =  2.2 /2 ;

  //APD + preamplifier

  //fPinDiodeSize[0] = 1.71 ;   //Values of ame PIN diode  
  //fPinDiodeSize[1] = 0.0280 ; // OHO 0.0280 is the depth of active layer
  //fPinDiodeSize[2] = 1.61 ;    
 
  fPinDiodeHalfSize[0] = 0.5000 /2 ;    // APD 5 mm side
  fPinDiodeHalfSize[1] = 0.0100 /2 ;    // APD bulk thickness
  fPinDiodeHalfSize[2] = 0.5000 /2 ;    // APD 5 mm side 

  fPreampHalfSize[0] = 1.5 / 2 ;       // Preamplifier
  fPreampHalfSize[1] = 0.5 / 2 ;
  fPreampHalfSize[2] = 1.5 / 2 ;

  //Strip unit (8x2 crystals)

  fNCellsXInStrip =  8 ;       //Number of crystals in strip unit along x-axis
  fNCellsZInStrip =  2 ;       //Number of crystals in strip unit along z-axis
  fNStripX        =  8 ;       //Number of strip units across along x-axis
  fNStripZ        = 28 ;       //Number of strips along z-axis

  fStripWallWidthOut = 0.01 ;  // Side to another strip  
  fStripWallWidthIn  = 0.02 ;  // Side betveen crystals in one strip

  fTyvecThickness = 0.0175 ;     //Thickness of the tyvec

  fAirGapLed = 1.5 - 2 * fPreampHalfSize[1] - 2 * fPinDiodeHalfSize[1] ; // Air gap before crystalls for LED system
                                           // Note, that Cell in Strip 1.5 longer then crystall

  //---Now calculate thechnical sizes for GEANT implementation

  fWrappedHalfSize[0] = (2*fTyvecThickness + 2*fCrystalHalfSize[0])/2 ;   //This will be size of crystall
  fWrappedHalfSize[1] = fCrystalHalfSize[1] ;                             //wrapped into tyvec
  fWrappedHalfSize[2] = (2*fTyvecThickness + 2*fCrystalHalfSize[2])/2 ;   //

  fAirCellHalfSize[0] = fWrappedHalfSize[0] + 0.01;
  fAirCellHalfSize[1] = (fAirGapLed + 2*fPreampHalfSize[1] + 
                       2*fPinDiodeHalfSize[1] + 2*fWrappedHalfSize[1])/2 ;  //in strip
  fAirCellHalfSize[2] = fWrappedHalfSize[2] + 0.01;

  //  fSupportPlateHalfSize[0] = ( (fNCellsXInStrip-1)*fStripWallWidthIn + 2*fStripWallWidthOut + 
  //			       fNCellsXInStrip * (2*fTyvecThickness + 2*fCrystalHalfSize[0]) )/2 ;
  fSupportPlateHalfSize[0] = 18.04  /2 ;
  fSupportPlateHalfSize[1] =   6.0  /2 ;
//  fSupportPlateHalfSize[2] = ( (fNCellsZInStrip-1)*fStripWallWidthIn + 2*fStripWallWidthOut +
//			       fNCellsZInStrip * (2*fTyvecThickness + 2*fCrystalHalfSize[2]) )/2;
  fSupportPlateHalfSize[2] =  4.51  /2 ;
  fSupportPlateThickness = 0.3 ;  
  fSupportPlateInHalfSize[0] = fSupportPlateHalfSize[0] ;                         //Half-sizes of the air
  fSupportPlateInHalfSize[1] = fSupportPlateHalfSize[1]-fSupportPlateThickness ;  //box in the support plate
  fSupportPlateInHalfSize[2] = fSupportPlateHalfSize[2]-fSupportPlateThickness/2 ;

  fStripHalfSize[0]= fSupportPlateHalfSize[0] ;  
  fStripHalfSize[1]= ( 2*fSupportPlateHalfSize[1] + 2*fAirCellHalfSize[1] )/2;      
  fStripHalfSize[2]= fSupportPlateHalfSize[2] ;

  // ------- Inner hermoinsulation ---------------
  fInnerThermoWidthX = 2.0 ;         // Width of the innerthermoinsulation across the beam
  fInnerThermoWidthY = 2.0 ;         // Width of the upper cover of innerthermoinsulation 
  fInnerThermoWidthZ = 2.0 ;         // Width of the innerthermoinsulation along the beam

  fInnerThermoHalfSize[0] = (2 * fStripHalfSize[0] * fNStripX + 2 * fInnerThermoWidthX ) /2 ;
  fInnerThermoHalfSize[1] = (2 * fStripHalfSize[1] + fInnerThermoWidthY ) /2 ; 
  fInnerThermoHalfSize[2] = (2 * fStripHalfSize[2] * fNStripZ + 2 * fInnerThermoWidthZ ) /2 ;

  // ------- Air gap between inner thermoinsulation and passive coller ---------

  fAirGapWidthX = 0.2 ;         // Width of the air gap across the beam
  fAirGapWidthY = 0.2 ;         // Width of the upper air gap
  fAirGapWidthZ = 0.2 ;         // Width of the air gap along the beam

  fAirGapHalfSize[0] = (2 * fInnerThermoHalfSize[0] + 2 * fAirGapWidthX ) /2 ;
  fAirGapHalfSize[1] = (2 * fInnerThermoHalfSize[1] +     fAirGapWidthY ) /2 ; 
  fAirGapHalfSize[2] = (2 * fInnerThermoHalfSize[2] + 2 * fAirGapWidthZ ) /2 ;
  
  // ------- Passive Cooler ------------------------

  fCoolerWidthX = 2.0 ;         // Width of the passive coller across the beam 
  fCoolerWidthY = 0.3 ;         // Width of the upper cover of cooler
  fCoolerWidthZ = 2.0 ;         // Width of the passive cooler along the beam

  fCoolerHalfSize[0] = (2 * fAirGapHalfSize[0] + 2 * fCoolerWidthX ) /2 ; 
  fCoolerHalfSize[1] = (2 * fAirGapHalfSize[1] +     fCoolerWidthY ) /2 ; 
  fCoolerHalfSize[2] = (2 * fAirGapHalfSize[2] + 2 * fCoolerWidthZ ) /2 ; 

  // ------- Outer thermoinsulation and Al cover -------------------------------
  
  fAlCoverThickness = 0.1 ;   //Thickness of the Al cover of the module  

  fOuterThermoWidthXUp  = 156.0 - fAlCoverThickness ; 
                                  //width of the upper surface of the PHOS module accross the beam 
  fOuterThermoWidthY    = 6.0 ;   // with of the upper cover of outer thermoinsulation
  fOuterThermoWidthZ    = 6.0 ;   //width of the thermoinsulation along the beam 

  fAlFrontCoverX = 6.0 ;   //Width of Al strip around fiberglass window: across
  fAlFrontCoverZ = 6.0 ;   //and along the beam


  // Calculate distance from IP to upper cover
  fIPtoOuterCoverDistance = fIPtoCrystalSurface - fAirGapLed - fInnerThermoWidthY - fAirGapWidthY - 
                            fCoolerWidthY - fOuterThermoWidthY - fAlCoverThickness ; 

  Float_t tanA = fOuterThermoWidthXUp / (2.*fIPtoOuterCoverDistance) ; 
                  // tan(a) where A = angle between IP to center and IP to side across beam

  fOuterThermoWidthXLow = fOuterThermoWidthXUp + 
                          2 * (2* fCoolerHalfSize[1] + fOuterThermoWidthY) * tanA  
                          - fAlCoverThickness ; 
                          //width of the lower surface of the COOL section accross the beam 


  fOuterThermoParams[0] = fOuterThermoWidthXUp / 2 ;   // half-length along x at the z surface positioned at -DZ; 
  fOuterThermoParams[1] = fOuterThermoWidthXLow/ 2 ;   // half-length along x at the z surface positioned at +DZ; 
  fOuterThermoParams[2] = ( 2 * fCoolerHalfSize[2] + 2 * fOuterThermoWidthZ ) / 2 ;
                                                       // `half-length along the y-axis' in out case this is z axis 
  fOuterThermoParams[3] = ( 2* fCoolerHalfSize[1] + fOuterThermoWidthY) /2 ;    
                                                       // `half-length along the z-axis' in our case this is y axis 

  fAlCoverParams[0] = fOuterThermoParams[0] + fAlCoverThickness ;
  fAlCoverParams[1] = fOuterThermoParams[1] + fAlCoverThickness ;
  fAlCoverParams[2] = fOuterThermoParams[2] + fAlCoverThickness ;
  fAlCoverParams[3] = fOuterThermoParams[3] + fAlCoverThickness /2 ;


  fFiberGlassHalfSize[0] = fAlCoverParams[0] - fAlFrontCoverX  ;
  fFiberGlassHalfSize[1] = fAlCoverParams[2] - fAlFrontCoverZ  ; //Note, here other ref. system
  fFiberGlassHalfSize[2] = fAlCoverThickness / 2 ;


  //============Now warm section======================
  //Al Cover 
  fWarmAlCoverWidthX = 2 * fAlCoverParams[1] ;  //Across beam
  fWarmAlCoverWidthY = 159.0 ;                  //along beam

  //T-support
  fTSupport1Thickness = 3.5 ; 
  fTSupport2Thickness = 5.0 ; 
  fTSupport1Width =  10.6 ;
  fTSupport2Width = 3.1 ;
  fNTSupports = fNStripX + 1 ;
  fTSupportDist = 7.48 ;

  //Air space for FEE
  fAirSpaceFeeX = 148.6 ;   //Across beam
  fAirSpaceFeeY = 135.0 ;   //along beam
  fAirSpaceFeeZ =  19.0 ;   //out of beam

  //thermoinsulation
  fWarmBottomThickness = 4.0 ;
  fWarmUpperThickness  = 4.0 ;

  //Frame 
  fFrameThickness = 5.0 ;
  fFrameHeight    = 15.0 ;

  //Fiberglass support
  fFiberGlassSup1X = 6.0 ;
  fFiberGlassSup1Y = 4.0 + fWarmUpperThickness ; 

  fFiberGlassSup2X = 3.0 ;
  fFiberGlassSup2Y = fFrameHeight ;

  //Now calculate Half-sizes

  fWarmAlCoverWidthZ = fAirSpaceFeeZ + fWarmBottomThickness + fWarmUpperThickness +  
                       fTSupport1Thickness + fTSupport2Thickness ;


  fWarmAlCoverHalfSize[0] = fWarmAlCoverWidthX / 2 ;
  fWarmAlCoverHalfSize[1] = fWarmAlCoverWidthY / 2 ;
  fWarmAlCoverHalfSize[2] = fWarmAlCoverWidthZ / 2 ;

  
  fWarmThermoHalfSize[0] = fWarmAlCoverHalfSize[0] - fAlCoverThickness ; 
  fWarmThermoHalfSize[1] = fWarmAlCoverHalfSize[1] - fAlCoverThickness ; 
  fWarmThermoHalfSize[2] = fWarmAlCoverHalfSize[2] - fAlCoverThickness /2 ; 


  //T-support
  fTSupport1HalfSize[0] =  fTSupport1Width /2 ;   //Across beam
  fTSupport1HalfSize[1] =  (fAirSpaceFeeY + 2*fFiberGlassSup1X) /2 ;    //along beam
  fTSupport1HalfSize[2] =  fTSupport1Thickness  /2;    //out of beam

  fTSupport2HalfSize[0] = fTSupport2Width /2;               //Across beam  
  fTSupport2HalfSize[1] = fTSupport1HalfSize[1] ; //along beam
  fTSupport2HalfSize[2] = fTSupport2Thickness /2; //out of beam

  //cables
  fTCables1HalfSize[0] = (2*fTSupport1HalfSize[0]*fNTSupports + (fNTSupports-1)* fTSupportDist) / 2 ; //Across beam
  fTCables1HalfSize[1] = fTSupport1HalfSize[1] ;    //along beam
  fTCables1HalfSize[2] = fTSupport1HalfSize[2] ;    //out of beam

  fTCables2HalfSize[0] = fTCables1HalfSize[0]  ; //Across beam
  fTCables2HalfSize[1] = fTSupport2HalfSize[1] ; //along beam
  fTCables2HalfSize[2] = fTSupport2HalfSize[2] ; //out of beam

  //frame: we define two frames along beam ...Z and across beam ...X
  fFrameXHalfSize[0] = (fAirSpaceFeeX + 2 * fFiberGlassSup2X + 2* fFrameThickness) /2 ;
  fFrameXHalfSize[1] = fFrameThickness /2 ;
  fFrameXHalfSize[2] = fFrameHeight    /2 ;

  fFrameXPosition[0] = 0 ;
  fFrameXPosition[1] = fAirSpaceFeeY /2 + fFiberGlassSup2X + fFrameXHalfSize[1] ;
  fFrameXPosition[2] = fWarmThermoHalfSize[2] - fFrameHeight/ 2 - fWarmBottomThickness ;
    
  fFrameZHalfSize[0] = fFrameThickness /2 ;
  fFrameZHalfSize[1] = (fAirSpaceFeeY + 2 * fFiberGlassSup2X) /2 ;
  fFrameZHalfSize[2] = fFrameHeight    /2 ;

  fFrameZPosition[0] = fAirSpaceFeeX /2 + fFiberGlassSup2X + fFrameZHalfSize[0] ;
  fFrameZPosition[1] = 0 ;
  fFrameZPosition[2] = fWarmThermoHalfSize[2] - fFrameHeight/ 2 - fWarmBottomThickness ;

  //Fiberglass support define 4 fiber glass supports 2 along Z  and 2 along X
  
  fFGupXHalfSize[0] = fFrameXHalfSize[0] ;
  fFGupXHalfSize[1] = fFiberGlassSup1X /2 ;
  fFGupXHalfSize[2] = fFiberGlassSup1Y /2;

  fFGupXPosition[0] = 0 ;
  fFGupXPosition[1] = fAirSpaceFeeY /2 + fFGupXHalfSize[1] ;
  fFGupXPosition[2] = fWarmThermoHalfSize[2] - fFrameHeight - fWarmBottomThickness - fFGupXHalfSize[2] ; 

  fFGupZHalfSize[0] = fFiberGlassSup1X /2 ;
  fFGupZHalfSize[1] = fAirSpaceFeeY /2 ;
  fFGupZHalfSize[2] = fFiberGlassSup1Y /2;

  fFGupZPosition[0] = fAirSpaceFeeX /2 + fFGupZHalfSize[0] ;
  fFGupZPosition[1] = 0 ;
  fFGupZPosition[2] = fWarmThermoHalfSize[2] - fFrameHeight - fWarmBottomThickness - fFGupXHalfSize[2] ; 

  fFGlowXHalfSize[0] = fFrameXHalfSize[0] - 2*fFrameZHalfSize[0] ;
  fFGlowXHalfSize[1] = fFiberGlassSup2X /2 ;
  fFGlowXHalfSize[2] = fFrameXHalfSize[2] ;

  fFGlowXPosition[0] = 0 ;
  fFGlowXPosition[1] = fAirSpaceFeeY /2 + fFGlowXHalfSize[1] ;
  fFGlowXPosition[2] = fWarmThermoHalfSize[2] - fWarmBottomThickness - fFGlowXHalfSize[2] ; 

  fFGlowZHalfSize[0] = fFiberGlassSup2X /2 ;
  fFGlowZHalfSize[1] = fAirSpaceFeeY /2 ;
  fFGlowZHalfSize[2] = fFrameZHalfSize[2] ;

  fFGlowZPosition[0] = fAirSpaceFeeX /2 + fFGlowZHalfSize[0] ;
  fFGlowZPosition[1] = 0 ;
  fFGlowZPosition[2] = fWarmThermoHalfSize[2] - fWarmBottomThickness - fFGlowXHalfSize[2] ; 


  // --- Air Gap for FEE ----

  fFEEAirHalfSize[0] =  fAirSpaceFeeX /2 ;
  fFEEAirHalfSize[1] =  fAirSpaceFeeY /2;
  fFEEAirHalfSize[2] =  fAirSpaceFeeZ /2;

  fFEEAirPosition[0] = 0 ;
  fFEEAirPosition[1] = 0 ;
  fFEEAirPosition[2] = fWarmThermoHalfSize[2] - fWarmBottomThickness - fFEEAirHalfSize[2] ;

  // --- Calculate the oveol dimentions of the EMC module
  
  fEMCParams[3] = fAlCoverParams[3] + fWarmAlCoverHalfSize[2] ; //Size out of beam
  fEMCParams[0] = fAlCoverParams[0] ; //Upper size across the beam
  fEMCParams[1] = (fAlCoverParams[1] - fAlCoverParams[0])*fEMCParams[3]/fAlCoverParams[3] 
                 + fAlCoverParams[0]  ; //Lower size across the beam
  fEMCParams[2] = fWarmAlCoverHalfSize[1] ;                     // Size along the beam

  fNPhi = fNStripX * fNCellsXInStrip ;    //number of crystals across the beam
  fNZ   = fNStripZ * fNCellsZInStrip ;    //number of crystals along the beam
}

