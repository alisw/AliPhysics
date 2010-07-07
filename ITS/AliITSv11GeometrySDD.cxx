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


//*************************************************************************
//
// SDD geometry, based on ROOT geometrical modeler
//
//
// This geometry has no dependence with aliroot, you can run it with root
// only, provided that the AliITSv11GeomCable classes are also compiled
//
// Ludovic Gaudichet                                   gaudichet@to.infn.it
//*************************************************************************


// $Id$


// General Root includes
#include <TMath.h>

// Root Geometry includes
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoCone.h>
#include <TGeoTube.h>
#include <TGeoTrd1.h>
#include <TGeoArb8.h>
#include <TGeoCompositeShape.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>
#include <TGeoPcon.h>

#include "AliITSgeom.h"
#include "AliITSgeomSDD.h"
#include "AliITSv11GeometrySDD.h"
#include "AliITSv11GeomCableFlat.h"
#include "AliITSv11GeomCableRound.h"

const char*    AliITSv11GeometrySDD::fgSDDsensitiveVolName3 = "ITSsddSensitivL3";
const char*    AliITSv11GeometrySDD::fgSDDsensitiveVolName4 = "ITSsddSensitivL4";
const Double_t AliITSv11GeometrySDD::fgkSegmentLength     = 37.2*2*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadderWidth       = 50.0*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadderHeight      = 30.0*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadderSegBoxDW    =  3.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadderSegBoxDH    =  3.*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkLadderBeamRadius  =  0.6*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadderLa          =  3.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadderHa          =  0.721979*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadderLb          =  3.7*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadderHb          =  0.890428*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadderl           =  0.25*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkBottomBeamAngle   = 56.5;
const Double_t AliITSv11GeometrySDD::fgkBeamSidePhi       = 65;

const Double_t AliITSv11GeometrySDD::fgkLadWaferSep       = 2*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkPinSuppWidth      = 2.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkPinSuppHeight     = 2.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkPinSuppRmax       = 2.5/2.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkPinR              = 1.5/2.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkPinSuppLength     = 5.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkPinSuppThickness  = 0.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkPinSuppConeAngle  = 4;
const Double_t AliITSv11GeometrySDD::fgkPinDXminOnSensor  = (39./2.)*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkPinPinDDXOnSensor = 3*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkPinDYOnSensor     = (52.5/2.)*fgkmm;

// parameters from ALR-0752/3
const Double_t AliITSv11GeometrySDD::fgkCoolPipeSuppHeight    =  3.2*fgkmm;  
const Double_t AliITSv11GeometrySDD::fgkCoolPipeSuppMaxLength = 14*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCoolPipeSuppWidthExt  =  0.4*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCoolPipeSuppWidthIn   =  0.65*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCoolPipeSuppHoleDiam  =  2*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCoolPipeSuppFulWidth  =  5.15*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCoolPipeSuppTongW     =  0.8*fgkmm; 
const Double_t AliITSv11GeometrySDD::fgkCoolPipeSuppAngle     = 22.5;
const Double_t AliITSv11GeometrySDD::fgkCoolPipeSuppSlitL     =  4.9*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCoolPipeSuppAxeDist   =  3.05*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCoolPipeInnerDiam     =  1.84*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCoolPipeOuterDiam     =  2.*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkBTBthick           =  0.25 *fgkmm;
const Double_t AliITSv11GeometrySDD::fgkBTBlength          = 55. *fgkmm;
const Double_t AliITSv11GeometrySDD::fgkBTBwidth           = 18*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkBTBaxisAtoBottom   =  4*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkBTBaxisAtoBase     =  1.2*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkRadiusAminBTB      =  1. *fgkmm;
const Double_t AliITSv11GeometrySDD::fgkRadiusBminBTB      =  0.53 *fgkmm;
const Double_t AliITSv11GeometrySDD::fgkBTBHoleLength      = 15 *fgkmm;
const Double_t AliITSv11GeometrySDD::fgkBTBHolewidth       =  6 *fgkmm;
const Double_t AliITSv11GeometrySDD::fgkBTBHoleRefX        = 10 *fgkmm;
const Double_t AliITSv11GeometrySDD::fgkBTBHoleRefY        =  6.5 *fgkmm;

const Double_t AliITSv11GeometrySDD::fgkLay3Rmin           = 145.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLay3Rmax           = 200.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLay3Length         = (524.+0.)*fgkmm; // ladder+supporting rings (length of the virtual tube)
const Double_t AliITSv11GeometrySDD::fgkLay3LadderLength   = 524.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLay3DetShortRadius = 146.0*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLay3DetLongRadius  = 152.0*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLay3LaddTopCornerEnd = 15.6*fgkmm;
const Int_t    AliITSv11GeometrySDD::fgkLay3Ndet           =  6;
const Int_t    AliITSv11GeometrySDD::fgkLay3Nladd          = 14;
const Double_t AliITSv11GeometrySDD::fgkLay3CoolPipeSuppH  =  7.5*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkLay4Rmin           = 234.8*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLay4Rmax           = 286.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLay4Length         = (671.+0.)*fgkmm;    // ladder+supporting rings (length of the virtual tube)
const Double_t AliITSv11GeometrySDD::fgkLay4LadderLength   = 671.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLay4DetShortRadius = 235.0*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLay4DetLongRadius  = 240.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLay4LaddTopCornerEnd = 15.6*fgkmm;
const Int_t    AliITSv11GeometrySDD::fgkLay4Ndet           = 8;
const Int_t    AliITSv11GeometrySDD::fgkLay4Nladd          = 22;
const Double_t AliITSv11GeometrySDD::fgkLay4CoolPipeSuppH  = 7.5*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkEndLaddCardsShortRadiusLay3 = fgkLay3DetShortRadius;
const Double_t AliITSv11GeometrySDD::fgkEndLaddCardsShortRadiusLay4 = fgkLay4DetShortRadius;
const Double_t AliITSv11GeometrySDD::fgkDistEndLaddCardsLadd = 0.*fgkmm;

//hybrid 
const Double_t AliITSv11GeometrySDD::fgkHybridAngle       = 46;           // approx !!!
// Origine taken at the hybrid corner :
const Double_t AliITSv11GeometrySDD::fgkHybridLength      = 65*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHybridWidth       = 41*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHybRndHoleRad     = 1.05*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHybRndHoleZ       = 2.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHybRndHoleX       = fgkHybridWidth-23.599*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkHybFLlowHoleDZ    =   9.698*fgkmm; 
const Double_t AliITSv11GeometrySDD::fgkHybFLlowHolePasDX =  10.754*fgkmm; 
const Double_t AliITSv11GeometrySDD::fgkHybFLlowHoleAmbDX =   9.122*fgkmm;
  // center of ships to the border
const Double_t AliITSv11GeometrySDD::fgkHybFLlowChipZ4    = fgkHybridLength-(4.654      )*fgkmm-fgkHybFLlowHoleDZ/2;
const Double_t AliITSv11GeometrySDD::fgkHybFLlowChipZ3    = fgkHybridLength-(4.654+15.  )*fgkmm-fgkHybFLlowHoleDZ/2;
const Double_t AliITSv11GeometrySDD::fgkHybFLlowChipZ2    = fgkHybridLength-(4.654+15.*2)*fgkmm-fgkHybFLlowHoleDZ/2;
const Double_t AliITSv11GeometrySDD::fgkHybFLlowChipZ1    = fgkHybridLength-(4.654+15.*3)*fgkmm-fgkHybFLlowHoleDZ/2;
const Double_t AliITSv11GeometrySDD::fgkHybFLlowPasX      = fgkHybridWidth-32.775*fgkmm;       
const Double_t AliITSv11GeometrySDD::fgkHybFLlowAmbX      = fgkHybridWidth-20.791*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHybChipsDZ        =  9.221*fgkmm; 
const Double_t AliITSv11GeometrySDD::fgkHybPascalDX       = 10.245*fgkmm; 
const Double_t AliITSv11GeometrySDD::fgkHybAmbraDX        =  8.51*fgkmm; 
const Double_t AliITSv11GeometrySDD::fgkHybFLUpperWidth   = 15.012*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHybFLUpperLength  = 59.878*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHybFLUpperAlDZ    = 11.183*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHybFLUpperAldx    =  2.307*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkmu = 1*fgkmicron; // 1*fgkmicron; // can be increase for checking thin objects
const Double_t AliITSv11GeometrySDD::fgkHybridThBridgeThick =  0.25*fgkmm;               // ???
const Double_t AliITSv11GeometrySDD::fgkHybAlThick         =  30*fgkmu;
const Double_t AliITSv11GeometrySDD::fgkHybUpThick         =  20*fgkmu;
const Double_t AliITSv11GeometrySDD::fgkHybGlueScrnThick   =  50*fgkmu;           // ??? ?????
const Double_t AliITSv11GeometrySDD::fgkHybGlueLowThick    =  90*fgkmu;
const Double_t AliITSv11GeometrySDD::fgkHybGlueUpThick     =  90*fgkmu;           // sur ?????
const Double_t AliITSv11GeometrySDD::fgkHybAlCCThick       =  12*fgkmu;
const Double_t AliITSv11GeometrySDD::fgkHybUpCCThick       =  12*fgkmu;
const Double_t AliITSv11GeometrySDD::fgkHybChipThick       = 150*fgkmu;
const Double_t AliITSv11GeometrySDD::fgkHybGlueAgThick     =  50*fgkmu;          // ??? ????
const Double_t AliITSv11GeometrySDD::fgkHybUnderNiThick    =  20*fgkmu;          // ??? ????
const Int_t    AliITSv11GeometrySDD::fgkNHybSMD = 25;
const Double_t AliITSv11GeometrySDD::fgkHybSMDposX[fgkNHybSMD]     = 
          {2.92*fgkmm,6.5*fgkmm,10.13*fgkmm,13.59*fgkmm,21.40*fgkmm,
	   2.92*fgkmm,6.5*fgkmm,10.13*fgkmm,13.59*fgkmm,19.91*fgkmm,
	   2.92*fgkmm,6.5*fgkmm,10.13*fgkmm,13.59*fgkmm,17.09*fgkmm,21.40*fgkmm,
	   2.92*fgkmm,6.5*fgkmm,10.13*fgkmm,13.59*fgkmm,19.91*fgkmm,
	   1.63*fgkmm,5.22*fgkmm,13.59*fgkmm,21.40*fgkmm};
const Double_t AliITSv11GeometrySDD::fgkHybSMDposZ[fgkNHybSMD]     = 
         { 2.3 *fgkmm, 2.3 *fgkmm, 2.3 *fgkmm, 2.3 *fgkmm, 2.3 *fgkmm,
	   17.315*fgkmm,17.315*fgkmm,17.315*fgkmm,17.315*fgkmm,17.315*fgkmm,
	   32.31*fgkmm,32.31*fgkmm,32.31*fgkmm,32.31*fgkmm,32.31*fgkmm,32.31*fgkmm,
	   47.38*fgkmm,47.38*fgkmm,47.38*fgkmm,47.38*fgkmm,47.38*fgkmm,
	   62.68*fgkmm,62.06*fgkmm,62.06*fgkmm,62.06*fgkmm};
const Double_t AliITSv11GeometrySDD::fgkHybSMDmiddleW      =   0.954*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHybSMDmiddleL      =   0.47 *fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHybSMDendW         =   1.132*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHybSMDendL         =   0.925*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHybSMDheight       = 400.*fgkmu;          // ??? ????!!!!!!!

const Double_t AliITSv11GeometrySDD::fgkWaferThickness     = 300.*fgkmu;
const Double_t AliITSv11GeometrySDD::fgkWaferWidth         =  72.5 *fgkmm;
const Double_t AliITSv11GeometrySDD::fgkWaferLength        =  87.6 *fgkmm;
const Double_t AliITSv11GeometrySDD::fgkWaferThickSens     = 299.8*fgkmu;
const Double_t AliITSv11GeometrySDD::fgkWaferWidthSens     =  70.17*fgkmm;
// 256 anodes times 294 microns of pitch
const Double_t AliITSv11GeometrySDD::fgkWaferLengthSens    =  256*294*fgkmicron;

const Double_t AliITSv11GeometrySDD::fgkDigitCablWidth     = 18.4*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkDigitCablAlThick   = (30+30*8./10.)*fgkmicron; // will probably change
const Double_t AliITSv11GeometrySDD::fgkDigitCablPolyThick = (20+12)*fgkmicron;        // will probably change

const Double_t AliITSv11GeometrySDD::fgkWaHVcableAlThick   = 30*2./10.*fgkmu;  // will probably change // Al ratio is random !!!
const Double_t AliITSv11GeometrySDD::fgkWaHVcablePolyThick = 175*fgkmu;        // will probably change
const Double_t AliITSv11GeometrySDD::fgkWaHVcableLength    = 67.08*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkWaHVcableWitdh     = 17.4*fgkmm;              //  check !!!
const Double_t AliITSv11GeometrySDD::fgkWaHVcableDW        = 5.24*fgkmm; //5.24*fgkmm;//  check !!!

const Double_t AliITSv11GeometrySDD::fgkSensorGlassLX      =   5.  *fgkmm; 
const Double_t AliITSv11GeometrySDD::fgkSensorGlassLZ      =   5.  *fgkmm; 
const Double_t AliITSv11GeometrySDD::fgkSensorGlassLY      = 150.  *fgkmu;
const Double_t AliITSv11GeometrySDD::fgkGlassDXOnSensor    =  26.28*fgkmm;             //  check !!!
const Double_t AliITSv11GeometrySDD::fgkGlassDZOnSensor    =  22.50*fgkmm;             //  check !!!

const Double_t AliITSv11GeometrySDD::fgkTransitHVAlThick    = 30*2./10.*fgkmu; //  check // will probably change //Al ratio is random
const Double_t AliITSv11GeometrySDD::fgkTransitHVPolyThick  = 100*fgkmu;       //  check  // will probably change
const Double_t AliITSv11GeometrySDD::fgkTransitHVHeadLX     =  71.46*fgkmm;           //  check !!!
const Double_t AliITSv11GeometrySDD::fgkTransitHVHeadLZ     =  21.3*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkTransitHVBondingLZ  =   3.6*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkTransitHVtailLength =  27*fgkmm;              // ???, not yet fixed ...
const Double_t AliITSv11GeometrySDD::fgkTransitHVtailWidth  =  26*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkTransitHVtailXpos   =   8*fgkmm;    //8*fgkmm          // ???, a mesurer !!!
const Double_t AliITSv11GeometrySDD::fgkTransitHVsideLZ     =  10.34*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkTransitHVsideLeftZ  =   4.11*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkTransitHVsideRightZ =   3.5*fgkmm;           // ???, a mesurer !!!

const Double_t AliITSv11GeometrySDD::fgkLongHVcablePolyThick= (20+30+125+30+20+30+125+30+20)*fgkmu; //  check  // will probably change
const Double_t AliITSv11GeometrySDD::fgkLongHVcableAlThick  = (30+30*2/10+30)*fgkmu;                //  check  // will probably change
const Double_t AliITSv11GeometrySDD::fgkLongHVcableSeparation = 600*fgkmicron;

const Double_t AliITSv11GeometrySDD::fgkRubyDX              =  14.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkRubyZladd3          = 250*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkRubyZladd4          = 325*fgkmm;

// the stesalite ladder foot at its end
const Double_t AliITSv11GeometrySDD::fgkLadFootX         = 60.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadFootZ         = 20.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadFootY         =  8.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadFootMiddleY    =  4.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadBox1X         = 23.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadFingerPrintX  =  6.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadFingerPrintY  =  1.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadFingerPrintBorder = 4.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkRubyCageHoleZ     =  8.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkRubyCageHoleX     =  9.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkRubyCageHoleY     =  6.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkRubyCageAxisShift =  0.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkScrewM4diam       =  4.*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkRubyScrewShiftToCenterY = 0.1;
const Double_t AliITSv11GeometrySDD::fgkRubyHoleDiam            = 0.5;

// the U cooling pipe and its heat exchanger in end-ladder cards system
const Double_t AliITSv11GeometrySDD::fgkEndLadPipeUlengthLay3 = 138*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkEndLadPipeUlengthLay4 = 150*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkEndLadPipeUwidth      =  59*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkEndLadPipeRadius      =   5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkEndLadPipeInnerDiam   =   2.8*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkEndLadPipeOuterDiam   =   3.*fgkmm;
//--- The al body of the cooling syst.of the heat exchanger :
const Double_t AliITSv11GeometrySDD::fgkEndLadPipeArmZLay3    = 112.*fgkmm;   //
const Double_t AliITSv11GeometrySDD::fgkEndLadPipeArmZLay4    = 125.*fgkmm;   //
const Double_t AliITSv11GeometrySDD::fgkEndLadPipeArmX        =   4.75*fgkmm; // the arms of the U cooling tube
const Double_t AliITSv11GeometrySDD::fgkEndLadPipeArmY        =   6.8*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkEndLadPipeArmBoxDY    =   1.03*fgkmm; // shift in Y of the arms from the axis
const Double_t AliITSv11GeometrySDD::fgkEndLadPipeArmBoxDX    =   0.125*fgkmm;// shift in X of the arms from the axis
const Double_t AliITSv11GeometrySDD::fgkEndLadPipeArmZpos     =   8.9*fgkmm;  // 

// LV card :
const Double_t AliITSv11GeometrySDD::fgkLVcardX     = 26.525*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLVcardY     = 44.95*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLVcardZ     = 1.*fgkmm; // all except Cu layer   //???
const Double_t AliITSv11GeometrySDD::fgkLVcardCuZ   = 0.1*fgkmm;   //???

const Double_t AliITSv11GeometrySDD::fgkLVChip0X    = 16.525*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLVChip0Y    = 10.8*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLVChip0Z    =  3.5*fgkmm; // all except si layer   //???
const Double_t AliITSv11GeometrySDD::fgkLVChip0SiZ  =  0.2*fgkmm; //???????????????????????????????????????????????????
const Double_t AliITSv11GeometrySDD::fgkLVChip0PosX = 13.*fgkmm; //19.95*fgkmm;  ???
const Double_t AliITSv11GeometrySDD::fgkLVChip0PosY = 10.3*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkLVChip1X    = 6.00*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLVChip1Y    = 6.00*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLVChip1Z    = 1*fgkmm;  // ???
const Double_t AliITSv11GeometrySDD::fgkLVChip1SiZ  = 0.2*fgkmm;  // ???
const Double_t AliITSv11GeometrySDD::fgkLVChip1PosX = 18.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLVChip1PosY = 27.6*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkLVChip2X    = 6.00*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLVChip2Y    = 6.00*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLVChip2Z    = 1*fgkmm;    // ???
const Double_t AliITSv11GeometrySDD::fgkLVChip2SiZ  = 0.2*fgkmm;  //???
const Double_t AliITSv11GeometrySDD::fgkLVChip2PosX = 18.0*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLVChip2PosY = 39.0*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkLVChip3X    = 4.01*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLVChip3Y    = 4.01*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLVChip3Z    = 1*fgkmm; // ???
const Double_t AliITSv11GeometrySDD::fgkLVChip3SiZ  = 0.2*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLVChip3PosX = 20.7*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLVChip3PosY = 21.4*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkLVcoolX1    = 17.25*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLVcoolY1    =  8.7*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLVcoolZ1    =  1.*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkLVcoolX2    =  3.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLVcoolY2    =  8.7*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLVcoolZ2    =  2.3*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkLVcoolX3    =  4.75*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLVcoolY3    =  3.1*fgkmm; //+0.1=glue
const Double_t AliITSv11GeometrySDD::fgkLVcoolPosY  = 6.5*fgkmm;

// HV card :
const Double_t AliITSv11GeometrySDD::fgkHVCardCeramX    = 54.01*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCeramY    = 40.89*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCeramZ    =  0.7*fgkmm;  // ???

const Double_t AliITSv11GeometrySDD::fgkHVCardCapa1X    =   6.8*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa1Z    =   1.*fgkmm;  // ???
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa1Ymid =   4.1*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa1Yend =   0.95*fgkmm; // doesn't take into account the soldering
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa1PosX =  13.1*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa1PosY =  14.5*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkHVCardCapa2X    =   6.8*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa2Z    =   1.*fgkmm;  // ???
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa2Ymid =   2.9*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa2Yend =   0.95*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa2PosX = -12.6*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa2PosY =  16.54*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkHVCardCapa3Xmid =   3.0*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa3Xend =   0.91*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa3Z    =   2.*fgkmm;      // ???
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa3Y    =   3.43*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkHVCardCapa3PosX1 =  14.6*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa3PosX2 =   7.2*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa3PosX3 =   2.52*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa3PosX4 =  -4.96*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa3PosX5 = -13.82*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa3PosY1 =   6.27*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa3PosY2 =   0.7*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCapa3PosY3 =   9.1*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkHVCardCool1X     =  14.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCool1Y     =   9.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCool1Z     =   2.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCool2X     =  14.25*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCool2Y     =   3.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCool2Z     =   4.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCool3X     =   4.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCool3Y     =   3.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCool3Z     =   7.2*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkHVCardCoolDY     =   6.*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkCarlosSuppX1     = 19.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosSuppY1     =  2*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosSuppX2     = 35.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosSuppY2     =  3.9*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosSuppZ      = 17.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosSuppAngle  = 45;
const Double_t AliITSv11GeometrySDD::fgkCarlosSuppX3     =  4.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosSuppY3     =  3.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosSuppZ3     = 12.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosSuppTopLen = 8.65*fgkmm;

// screws fixing boards to the end-ladder on the U tube
const Double_t AliITSv11GeometrySDD::fgkLittleScrewHeadR   = 1.85*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLittleScrewHeadH   = 1.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLittleScrewR       = 0.7*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkShiftLittleScrewLV = 3*fgkmm;     // ???
const Double_t AliITSv11GeometrySDD::fgkLittleLVScrewHeadR = 1.2*fgkmm;   // ???

// CARLOS board
const Double_t AliITSv11GeometrySDD::fgkCarlosCardX1          = (25.50+28.50)*fgkmm; // length (first part of Carlos card)
const Double_t AliITSv11GeometrySDD::fgkCarlosCardY1          =    1.6*fgkmm;        // thickness
const Double_t AliITSv11GeometrySDD::fgkCarlosCardZ1          =   40.8*fgkmm;        // width 
const Double_t AliITSv11GeometrySDD::fgkCarlosCardCuY         =    0.1*fgkmm;   // thickness of Cu layer (strips)
const Double_t AliITSv11GeometrySDD::fgkCarlosCardX2          =   25.50*fgkmm;  // length (2nd part of Carlos card)
const Double_t AliITSv11GeometrySDD::fgkCarlosCardZ2          =    8.20*fgkmm;  // width 

const Double_t AliITSv11GeometrySDD::fgkCarlosCardChipSiThick =   0.1*fgkmm;  // ??? idem for all chips ???
const Double_t AliITSv11GeometrySDD::fgkCarlosCardShift       =   9*fgkmm;   // ??? shift in z w.r.t. heat bridge 

// size and position of various chips on carlos end-ladder board
const Double_t AliITSv11GeometrySDD::fgkCarlosU1X             =  13*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosU1Y             =   1.68*fgkmm; 
const Double_t AliITSv11GeometrySDD::fgkCarlosU1Z             =  13*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosU1posX          =  18.4*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosU1posZ          =  -7.2*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkCarlosU2X             =  13.75*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosU2Y             =   1.60*fgkmm; 
const Double_t AliITSv11GeometrySDD::fgkCarlosU2Z             =  13.85*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosU2posX          =  -0.375*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosU2posZ          =  -9.725*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkCarlosU3X             =   5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosU3Y             =   1.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosU3Z             =   5*fgkmm; 
const Double_t AliITSv11GeometrySDD::fgkCarlosU3posX          =   6.4*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosU3posZ          =   9.9*fgkmm;

// U4 like U3  
const Double_t AliITSv11GeometrySDD::fgkCarlosU4posX          = -12*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosU4posZ          =   3.6*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkCarlosU17X            =  16*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosU17Y            =   3.5*fgkmm; 
const Double_t AliITSv11GeometrySDD::fgkCarlosU17Z            =  10.9*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosU17posX         = -17.84*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosU17posZ         = -10.95*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkCarlosU35X            =   4*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosU35Y            =   1.*fgkmm; 
const Double_t AliITSv11GeometrySDD::fgkCarlosU35Z            =   4*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosU35posX         = -21.6*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosU35posZ         =   2.3*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkCarlosU36X            =   6*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosU36Y            =   1.*fgkmm; 
const Double_t AliITSv11GeometrySDD::fgkCarlosU36Z            =   6*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosU36posX         = -21.6*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosU36posZ         =   9.6*fgkmm;
  
const Double_t AliITSv11GeometrySDD::fgkCarlosQZ1X            =   8*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosQZ1Y            =   1.7*fgkmm; // look thicker than design number (0.7) ! ??? 
const Double_t AliITSv11GeometrySDD::fgkCarlosQZ1Z            =   3.7*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosQZ1posX         = -12*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkCarlosQZ1posZ         =  10.6*fgkmm;

// distance from the heat bridge center to the card center :
const Double_t AliITSv11GeometrySDD::fgkCarlosCard2HeatBridge = fgkCarlosSuppY2/2+fgkCarlosCardY1/2+fgkCarlosU1Y+0.1*fgkmm;

// some pieces at the end of the carbon fiber ladder
 const Double_t AliITSv11GeometrySDD::fgkCoolPipeLay3Len   = 467.*fgkmm;  // ???
 const Double_t AliITSv11GeometrySDD::fgkCoolPipeLay4Len   = 616.*fgkmm;  // ???
 const Double_t AliITSv11GeometrySDD::fgkHVguideX1         =  42.5*fgkmm;
 const Double_t AliITSv11GeometrySDD::fgkHVguideY1         =   7.*fgkmm;
 const Double_t AliITSv11GeometrySDD::fgkHVguideZ1         =  10.*fgkmm;
 const Double_t AliITSv11GeometrySDD::fgkHVguideZ2         =   6.*fgkmm;
 const Double_t AliITSv11GeometrySDD::fgkHVguideDX         =  -8.5*fgkmm;
 const Double_t AliITSv11GeometrySDD::fgkHVguideSuppFullZ  = 37.5*fgkmm;

// Cooling connector between phynox and plastic cooling water tubes
const Double_t AliITSv11GeometrySDD::fgkConnectorCoolTubeRmin = 1 *fgkmm;
const Double_t AliITSv11GeometrySDD::fgkConnectorCoolTubeR1   = 2.5*fgkmm; // ???
const Double_t AliITSv11GeometrySDD::fgkConnectorCoolTubeL1   = 3.*fgkmm;  // ???
const Double_t AliITSv11GeometrySDD::fgkConnectorCoolTubeR2   = 3.5*fgkmm;  // ???
const Double_t AliITSv11GeometrySDD::fgkConnectorCoolTubeL2   = 2.*fgkmm;  // ???
const Double_t AliITSv11GeometrySDD::fgkConnectorCoolTubeR3   = 3.*fgkmm;  // ???
const Double_t AliITSv11GeometrySDD::fgkConnectorCoolTubeL3   = 5 *fgkmm;  // ???


// parameters for coding SDD cables on SDD and SSD cones
const Double_t AliITSv11GeometrySDD::fgkSectionCuPerMod = 3*2*0.006 + 3*2*0.0005 + 2*0.002;
//                                             copper :     LV     +  signal  +  HV(HV ???)
const Double_t AliITSv11GeometrySDD::fgkSectionPlastPerMod = (TMath::Pi()*(3*0.36*0.36/4 + 3*0.21*0.21/4 + 2*0.115*0.115/4) 
							      - fgkSectionCuPerMod);

const Double_t AliITSv11GeometrySDD::fgkSectionGlassPerMod = 3*0.006; // ???
const Double_t AliITSv11GeometrySDD::fgkSectionCoolPolyuEL = 0.4672;
const Double_t AliITSv11GeometrySDD::fgkSectionCoolWaterEL = 0.3496;
// (sections are given in cm square)
const Double_t AliITSv11GeometrySDD::fgkCableBendRatio = 1.3; // ??? this factor account for the bending of cables

const Double_t AliITSv11GeometrySDD::fgkConeSDDr1 = 11.87574*fgkcm;
const Double_t AliITSv11GeometrySDD::fgkConeSDDr2 = 26.07574*fgkcm;
const Double_t AliITSv11GeometrySDD::fgkConeSDDz1 =  3.36066*fgkcm + 186.0*fgkmm + 0.5*790.0*fgkmm - 19.18934*fgkcm - 1.6;
const Double_t AliITSv11GeometrySDD::fgkConeSDDz2 = 17.56066*fgkcm + 186.0*fgkmm + 0.5*790.0*fgkmm - 19.18934*fgkcm - 1.6;
// These last parameters come from cone's code and define the slope
// and position of the SDD cone end.  For some unknown reason, this doesn't
// allow to stick on the SDD cone. This has to be checked when a correct
// version of the cone is available ... For now 'm applying some approximative
// corrections

const Double_t AliITSv11GeometrySDD::fgkSDDCableR1    = 16*fgkcm; // ??? // part 1 of "cable cone"
const Double_t AliITSv11GeometrySDD::fgkSDDCableR2    = 23*fgkcm; // ??? // part 1/2 of "cable cone"
const Double_t AliITSv11GeometrySDD::fgkSDDCableR3    = 26*fgkcm; // ??? // part 2 of "cable cone"

const Double_t AliITSv11GeometrySDD::fgkSDDCableDZint =  3.5*fgkcm;
const Double_t AliITSv11GeometrySDD::fgkSDDCableR5    =  37*fgkcm; // third part of "cable cone"
const Double_t AliITSv11GeometrySDD::fgkSDDCableZ5    =  65*fgkcm; // third part of "cable cone"







ClassImp(AliITSv11GeometrySDD)

//________________________________________________________________________
AliITSv11GeometrySDD::AliITSv11GeometrySDD(): 
  AliITSv11Geometry(),
  fPinSupport(0),
  fCoolPipeSupportL(0),
  fCoolPipeSupportR(0),
  fSDDsensor3(0),
  fSDDsensor4(0),
  fBaseThermalBridge(0),
  fHybrid(0),
  fLadderFoot(0),
  fCardLVR(0),
  fCardLVL(0),
  fCardHV(0),
  fCardCarlos(0),
  fRaccordoL(0),
  fDigitCableLay3A(0),
  fDigitCableLay3B(0),
  fDigitCableLay4A(0),
  fDigitCableLay4B(0),
  fMotherVol(0),
  fAddHybrids(kTRUE), 
  fAddSensors(kTRUE),
  fAddHVcables(kTRUE),
  fAddCables(kTRUE), 
  fAddCoolingSyst(kTRUE),
  fCoolingOn(kTRUE),
  fAddOnlyLadder3min(-1),
  fAddOnlyLadder3max(-1),
  fAddOnlyLadder4min(-1),
  fAddOnlyLadder4max(-1),
  fColorCarbonFiber(4),
  fColorRyton(5),
  fColorPhynox(7),
  fColorSilicon(3),
  fColorAl(7),
  fColorPolyhamide(5),
  fColorGlass(2),
  fColorSMD(12),
  fColorSMDweld(17),
  fColorStesalite(20),
  fLay3LadderUnderSegDH(0),
  fLay4LadderUnderSegDH(0),
  fLay3LaddShortRadius(0),
  fLay3LaddLongRadius(0),
  fLay4LaddShortRadius(0),
  fLay4LaddLongRadius(0)
{
  //
  // Standard constructor
  //
  SetParameters();
}


//________________________________________________________________________
AliITSv11GeometrySDD::AliITSv11GeometrySDD(Int_t debug) :
  AliITSv11Geometry(debug),
  fPinSupport(0),
  fCoolPipeSupportL(0),
  fCoolPipeSupportR(0),
  fSDDsensor3(0),
  fSDDsensor4(0),
  fBaseThermalBridge(0),
  fHybrid(0),
  fLadderFoot(0),
  fCardLVR(0),
  fCardLVL(0),
  fCardHV(0),
  fCardCarlos(0),
  fRaccordoL(0),
  fDigitCableLay3A(0),
  fDigitCableLay3B(0),
  fDigitCableLay4A(0),
  fDigitCableLay4B(0),
  fMotherVol(0),
  fAddHybrids(kTRUE), 
  fAddSensors(kTRUE),
  fAddHVcables(kTRUE),
  fAddCables(kTRUE), 
  fAddCoolingSyst(kTRUE),
  fCoolingOn(kTRUE),
  fAddOnlyLadder3min(-1),
  fAddOnlyLadder3max(-1),
  fAddOnlyLadder4min(-1),
  fAddOnlyLadder4max(-1),
  fColorCarbonFiber(4),
  fColorRyton(5),
  fColorPhynox(7),
  fColorSilicon(3),
  fColorAl(7),
  fColorPolyhamide(5),
  fColorGlass(2),
  fColorSMD(12),
  fColorSMDweld(17),
  fColorStesalite(20),
  fLay3LadderUnderSegDH(0),
  fLay4LadderUnderSegDH(0),
  fLay3LaddShortRadius(0),
  fLay3LaddLongRadius(0),
  fLay4LaddShortRadius(0),
  fLay4LaddLongRadius(0)
{
  //
  // Constructor setting debugging level
  //
  SetParameters();
}

//________________________________________________________________________
AliITSv11GeometrySDD::AliITSv11GeometrySDD(const AliITSv11GeometrySDD &s) :
  AliITSv11Geometry(s.GetDebug()),
  fPinSupport(s.fPinSupport),
  fCoolPipeSupportL(s.fCoolPipeSupportL),
  fCoolPipeSupportR(s.fCoolPipeSupportR),
  fSDDsensor3(s.fSDDsensor3),
  fSDDsensor4(s.fSDDsensor4),
  fBaseThermalBridge(s.fBaseThermalBridge),
  fHybrid(s.fHybrid),
  fLadderFoot(s.fLadderFoot),
  fCardLVR(s.fCardLVR),
  fCardLVL(s.fCardLVL),
  fCardHV(s.fCardHV),
  fCardCarlos(s.fCardCarlos),
  fRaccordoL(s.fRaccordoL),
  fDigitCableLay3A(s.fDigitCableLay3A),
  fDigitCableLay3B(s.fDigitCableLay3B),
  fDigitCableLay4A(s.fDigitCableLay4A),
  fDigitCableLay4B(s.fDigitCableLay4B),
  fMotherVol(s.fMotherVol),
  fAddHybrids(s.fAddHybrids), 
  fAddSensors(s.fAddSensors),
  fAddHVcables(s.fAddHVcables),
  fAddCables(s.fAddCables), 
  fAddCoolingSyst(s.fAddCoolingSyst),
  fCoolingOn(s.fCoolingOn),
  fAddOnlyLadder3min(s.fAddOnlyLadder3min),
  fAddOnlyLadder3max(s.fAddOnlyLadder3max),
  fAddOnlyLadder4min(s.fAddOnlyLadder4min),
  fAddOnlyLadder4max(s.fAddOnlyLadder4max),
  fColorCarbonFiber(s.fColorCarbonFiber),
  fColorRyton(s.fColorRyton),
  fColorPhynox(s.fColorPhynox),
  fColorSilicon(s.fColorSilicon),
  fColorAl(s.fColorAl),
  fColorPolyhamide(s.fColorPolyhamide),
  fColorGlass(s.fColorGlass),
  fColorSMD(s.fColorSMD),
  fColorSMDweld(s.fColorSMDweld),
  fColorStesalite(s.fColorStesalite),
  fLay3LadderUnderSegDH(s.fLay3LadderUnderSegDH),
  fLay4LadderUnderSegDH(s.fLay4LadderUnderSegDH),
  fLay3LaddShortRadius(s.fLay3LaddShortRadius),
  fLay3LaddLongRadius(s.fLay3LaddLongRadius),
  fLay4LaddShortRadius(s.fLay4LaddShortRadius),
  fLay4LaddLongRadius(s.fLay4LaddLongRadius)
{
  //     Copy Constructor
  // do only a "shallow copy" ...
  SetParameters();
}

//________________________________________________________________________
AliITSv11GeometrySDD& AliITSv11GeometrySDD::
operator=(const AliITSv11GeometrySDD &s) {
  //     Assignment operator
  if(&s == this) return *this;
  fMotherVol = s.fMotherVol;
  fAddHybrids = s.fAddHybrids;
  fAddSensors = s.fAddSensors;
  fAddHVcables = s.fAddHVcables;
  fAddCables = s.fAddCables;
  fAddCoolingSyst = s.fAddCoolingSyst;
  fCoolingOn = s.fCoolingOn;
  fAddOnlyLadder3min = s.fAddOnlyLadder3min;
  fAddOnlyLadder3max = s.fAddOnlyLadder3max;
  fAddOnlyLadder4min = s.fAddOnlyLadder4min;
  fAddOnlyLadder4max = s.fAddOnlyLadder4max;
  return *this;
}

//________________________________________________________________________
AliITSv11GeometrySDD::~AliITSv11GeometrySDD() {
  // Look like a destructor
  // Smell like a destructor
  // And actually is the destructor
  if (fDigitCableLay3A) delete [] fDigitCableLay3A;
  if (fDigitCableLay3B) delete [] fDigitCableLay3B;
  if (fDigitCableLay4A) delete [] fDigitCableLay4A;
  if (fDigitCableLay4B) delete [] fDigitCableLay4B;
}

//________________________________________________________________________
void AliITSv11GeometrySDD::SetParameters() {
  //
  // Define display colors and the non constant geometry parameters
  //

  Double_t detLadderDist = 8*fgkmm; 

  fLay3LadderUnderSegDH = detLadderDist - (fgkWaHVcableAlThick+fgkWaHVcablePolyThick);
  fLay4LadderUnderSegDH = detLadderDist - (fgkWaHVcableAlThick+fgkWaHVcablePolyThick);

  // radius from the center to the CF ladder :
  fLay3LaddShortRadius = (fgkLay3DetShortRadius
			  + fgkLadWaferSep+2*fgkWaferThickness
			  + detLadderDist); 
  fLay3LaddLongRadius  = (fgkLay3DetLongRadius
			  + fgkLadWaferSep+2*fgkWaferThickness
			  + detLadderDist); 
  fLay4LaddShortRadius = (fgkLay4DetShortRadius
			  + fgkLadWaferSep+2*fgkWaferThickness
			  + detLadderDist); 
  fLay4LaddLongRadius  = (fgkLay4DetLongRadius
			  + fgkLadWaferSep+2*fgkWaferThickness
			  + detLadderDist); 

  fLay3sensorZPos[0]= (  35.8+72.4+75.8 )*fgkmm;
  fLay3sensorZPos[1]= (  35.8+72.4      )*fgkmm;
  fLay3sensorZPos[2]= (  35.8           )*fgkmm;
  fLay3sensorZPos[3]= ( -37.9           )*fgkmm;
  fLay3sensorZPos[4]= ( -37.9-74.9      )*fgkmm;
  fLay3sensorZPos[5]= ( -37.9-74.9-71.1 )*fgkmm;

  fLay4sensorZPos[0] = (  38.5+73.2+75.4+71.6 )*fgkmm;
  fLay4sensorZPos[1] = (  38.5+73.2+75.4      )*fgkmm;
  fLay4sensorZPos[2] = (  38.5+73.2           )*fgkmm;
  fLay4sensorZPos[3] = (  38.5                )*fgkmm;
  fLay4sensorZPos[4] = ( -35.6                )*fgkmm;
  fLay4sensorZPos[5] = ( -35.6-74.8           )*fgkmm;
  fLay4sensorZPos[6] = ( -35.6-74.8-72.4      )*fgkmm;
  fLay4sensorZPos[7] = ( -35.6-74.8-72.4-76.  )*fgkmm;
}


//________________________________________________________________________
TGeoMedium* AliITSv11GeometrySDD::GetMedium(const char* mediumName) {
  //
  // Called to get a medium, checks that it exists.
  // If not, prints an error and returns 0
  //

  char ch[30];
  sprintf(ch, "ITS_%s",mediumName);
  TGeoMedium* medium =  gGeoManager->GetMedium(ch);
  if (! medium)
    printf("Error(AliITSv11GeometrySDD)::medium %s not found !\n", mediumName);
  return medium;
}


//________________________________________________________________________
Int_t AliITSv11GeometrySDD::GetLay3NLadders() const{
  // Get the actual number of ladder in layer 3
  if ( (fAddOnlyLadder3min<0) ||
       (fAddOnlyLadder3min >= fgkLay3Nladd) ||
       (fAddOnlyLadder3max<0) ||
       (fAddOnlyLadder3max >= fgkLay3Nladd) )
    return fgkLay3Nladd;
  else return (fAddOnlyLadder3max-fAddOnlyLadder3min+1);
}


//________________________________________________________________________
Int_t AliITSv11GeometrySDD::GetLay4NLadders() const{
  // Get the actual number of ladder in layer 4
  if ( (fAddOnlyLadder4min<0) ||
       (fAddOnlyLadder4min >= fgkLay4Nladd) ||
       (fAddOnlyLadder4max<0) ||
       (fAddOnlyLadder4max >= fgkLay4Nladd) )
    return fgkLay4Nladd;
  else return (fAddOnlyLadder4max-fAddOnlyLadder4min+1);
}


//________________________________________________________________________
void AliITSv11GeometrySDD::CreateBasicObjects() {
  //
  // Create basics objets which will be assembled together
  // in Layer3 and Layer4 functions
  //


  fDigitCableLay3A = new AliITSv11GeomCableFlat[fgkLay3Ndet];
  fDigitCableLay3B = new AliITSv11GeomCableFlat[fgkLay3Ndet];
  fDigitCableLay4A = new AliITSv11GeomCableFlat[fgkLay4Ndet];
  fDigitCableLay4B = new AliITSv11GeomCableFlat[fgkLay4Ndet];

  fPinSupport = CreatePinSupport();
  fCoolPipeSupportL = CreateCoolPipeSupportL();
  fCoolPipeSupportR = CreateCoolPipeSupportR();
  CreateSDDsensor();
  fBaseThermalBridge = CreateBaseThermalBridge();
  fHybrid = CreateHybrid(0);

  TGeoMedium *carbonFiberLadderStruct = GetMedium("SDD C AL (M55J)$"); //ITSsddCarbonM55J
  TGeoMedium *polyhamideSDD   = GetMedium("SDDKAPTON (POLYCH2)$");//ITSsddKAPTON_POLYCH2
  TGeoMedium *alSDD           = GetMedium("AL$"); //ITSal
  TGeoMedium *stainless       = GetMedium("INOX$"); // for screws, what is the material ???????????
  TGeoMedium *coolerMediumSDD = GetMedium("WATER$");
  TGeoMedium *raccordMedium   = GetMedium("INOX$");  // ??? material of raccordo ???

  //********************************************************************
  // pieces of the carbon fiber structure
  //********************************************************************
  Double_t dy             = fgkLadderSegBoxDH/2;
  Double_t triangleHeight = fgkLadderHeight - fgkLadderBeamRadius;
  Double_t halfTheta      = TMath::ATan( 0.5*fgkLadderWidth/triangleHeight );
  Double_t alpha          = TMath::Pi()*3./4. - halfTheta/2.;
  Double_t beta           = (TMath::Pi() - 2.*halfTheta)/4.;
  Double_t dYTranslation  = (fgkLadderHeight/2.
			     -0.5*fgkLadderWidth*TMath::Tan(beta)
			     -fgkLadderBeamRadius);
  Double_t distCenterSideDown =  0.5*fgkLadderWidth/TMath::Cos(beta);

  //--- the top V of the Carbon Fiber Ladder (segment)
  TGeoArb8 *cfLaddTop1 = CreateLadderSide( "CFladdTopCornerVol1shape",
					   fgkSegmentLength/2., halfTheta, 
			  -1, fgkLadderLa, fgkLadderHa, fgkLadderl);
  TGeoVolume *cfLaddTopVol1 = new TGeoVolume("ITSsddCFladdTopCornerVol1",
				  cfLaddTop1,carbonFiberLadderStruct);
  TGeoArb8 *cfLaddTop2 = CreateLadderSide( "CFladdTopCornerVol2shape",
					   fgkSegmentLength/2., halfTheta,
			   1, fgkLadderLa, fgkLadderHa, fgkLadderl);
  TGeoVolume *cfLaddTopVol2 = new TGeoVolume("ITSsddCFladdTopCornerVol2",
				  cfLaddTop2, carbonFiberLadderStruct);
  cfLaddTopVol1->SetLineColor(fColorCarbonFiber);
  cfLaddTopVol2->SetLineColor(fColorCarbonFiber);
  TGeoTranslation *trTop1 = new TGeoTranslation(0, fgkLadderHeight/2-dy, 0);

  //--- the 2 side V
  TGeoArb8 *cfLaddSide1 = CreateLadderSide( "CFladdSideCornerVol1shape",
					    fgkSegmentLength/2., beta, -1,
					    fgkLadderLb, fgkLadderHb, fgkLadderl);
  TGeoVolume *cfLaddSideVol1 = new TGeoVolume( "ITSsddCFladdSideCornerVol1",
				   cfLaddSide1,carbonFiberLadderStruct);
  TGeoArb8 *cfLaddSide2 = CreateLadderSide( "CFladdSideCornerVol2shape",
					    fgkSegmentLength/2., beta, 1,
					    fgkLadderLb, fgkLadderHb, fgkLadderl);
  TGeoVolume *cfLaddSideVol2 = new TGeoVolume( "ITSsddCFladdSideCornerVol2",
				   cfLaddSide2,carbonFiberLadderStruct);
  cfLaddSideVol1->SetLineColor(fColorCarbonFiber);
  cfLaddSideVol2->SetLineColor(fColorCarbonFiber);
  TGeoCombiTrans *ctSideR = CreateCombiTrans("", distCenterSideDown, 0,
					     alpha*TMath::RadToDeg());
  AddTranslationToCombiTrans(ctSideR, 0, -dYTranslation-dy, 0);
  TGeoCombiTrans *ctSideL = CreateCombiTrans("", distCenterSideDown,0,
					     -alpha*TMath::RadToDeg());
  AddTranslationToCombiTrans(ctSideL, 0, -dYTranslation-dy, 0);

  //--- The beams
  // Beams on the sides
  Double_t beamPhiPrime = TMath::ASin(1./TMath::Sqrt( (1+TMath::Sin(2*beta)*
	   TMath::Sin(2*beta)/(TanD(fgkBeamSidePhi)*TanD(fgkBeamSidePhi))) ));
  //cout<<"Phi prime = "<<beamPhiPrime*TMath::RadToDeg()<<endl;
  Double_t beamLength = TMath::Sqrt( fgkLadderHeight*fgkLadderHeight/
			( TMath::Sin(beamPhiPrime)*TMath::Sin(beamPhiPrime))
			+ fgkLadderWidth*fgkLadderWidth/4.)-fgkLadderLa/2-fgkLadderLb/2;
  TGeoTubeSeg *sideBeamS = new TGeoTubeSeg(0, fgkLadderBeamRadius,beamLength/2.,
					   0, 180);
  TGeoVolume *sideBeam = new TGeoVolume("ITSsddCFSideBeamVol", sideBeamS,
			     carbonFiberLadderStruct);
  sideBeam->SetLineColor(fColorCarbonFiber);

  //Euler rotation : about Z, then new X, then new Z
  TGeoRotation *beamRot1 = new TGeoRotation("", 90-2.*beta*TMath::RadToDeg(),
					    -beamPhiPrime*TMath::RadToDeg(),-90);
  TGeoRotation *beamRot2 = new TGeoRotation("", 90-2.*beta*TMath::RadToDeg(),
					    beamPhiPrime*TMath::RadToDeg(), -90);
  TGeoRotation *beamRot3 = new TGeoRotation("", 90+2.*beta*TMath::RadToDeg(),
					    beamPhiPrime*TMath::RadToDeg(), -90);
  TGeoRotation *beamRot4 = new TGeoRotation("", 90+2.*beta*TMath::RadToDeg(),
					    -beamPhiPrime*TMath::RadToDeg(),-90);

  TGeoCombiTrans *beamTransf[8];
  beamTransf[0] = new TGeoCombiTrans( 0.5*triangleHeight*
				      TMath::Tan(halfTheta),
				      fgkLadderBeamRadius/2. - dy,
				      -3*fgkSegmentLength/8, beamRot1);

  beamTransf[1] = new TGeoCombiTrans( 0.5*triangleHeight*
				      TMath::Tan(halfTheta),
				      fgkLadderBeamRadius/2. - dy,
				      -3*fgkSegmentLength/8, beamRot1);
  AddTranslationToCombiTrans(beamTransf[1], 0, 0, fgkSegmentLength/2);

  beamTransf[2] = new TGeoCombiTrans(0.5*triangleHeight*
				     TMath::Tan(halfTheta),
				     fgkLadderBeamRadius/2. - dy,
				     -fgkSegmentLength/8, beamRot2);

  beamTransf[3] = new TGeoCombiTrans(0.5*triangleHeight*
				     TMath::Tan(halfTheta),
				     fgkLadderBeamRadius/2. - dy,
				     -fgkSegmentLength/8, beamRot2);
  AddTranslationToCombiTrans(beamTransf[3], 0, 0, fgkSegmentLength/2);

  beamTransf[4] = new TGeoCombiTrans(-0.5*triangleHeight*
				     TMath::Tan(halfTheta),
				     fgkLadderBeamRadius/2. - dy,
				     -3*fgkSegmentLength/8, beamRot3);

  beamTransf[5] = new TGeoCombiTrans(-0.5*triangleHeight*
				     TMath::Tan(halfTheta),
				     fgkLadderBeamRadius/2. - dy,
				     -3*fgkSegmentLength/8, beamRot3);
  AddTranslationToCombiTrans(beamTransf[5], 0, 0, fgkSegmentLength/2);

  beamTransf[6] = new TGeoCombiTrans(-0.5*triangleHeight*
  TMath::Tan(halfTheta),fgkLadderBeamRadius/2.-dy, -fgkSegmentLength/8,beamRot4);
  beamTransf[7] = new TGeoCombiTrans(-0.5*triangleHeight*
  TMath::Tan(halfTheta),fgkLadderBeamRadius/2.-dy,3*fgkSegmentLength/8,beamRot4);

  //--- Beams of the bottom
  TGeoTubeSeg *bottomBeam1 = new TGeoTubeSeg(0, fgkLadderBeamRadius,
			         fgkLadderWidth/2.-fgkLadderLb/3, 0, 180);
  TGeoVolume *bottomBeam1Vol = new TGeoVolume("ITSsddBottomBeam1Vol",
				   bottomBeam1, carbonFiberLadderStruct);
  bottomBeam1Vol->SetLineColor(fColorCarbonFiber);
  TGeoTubeSeg *bottomBeam2 = new TGeoTubeSeg(0, fgkLadderBeamRadius,
			         fgkLadderWidth/2.-fgkLadderLb/3, 0, 90);
  TGeoVolume *bottomBeam2Vol = new TGeoVolume("ITSsddBottomBeam2Vol",
                               bottomBeam2, carbonFiberLadderStruct);
  bottomBeam2Vol->SetLineColor(fColorCarbonFiber);
  TGeoTubeSeg *bottomBeam3 = new TGeoTubeSeg(0, fgkLadderBeamRadius,
		             0.5*fgkLadderWidth/SinD(fgkBottomBeamAngle)
		             - fgkLadderLb/3, 0, 180);
  TGeoVolume *bottomBeam3Vol = new TGeoVolume("ITSsddBottomBeam3Vol",
				   bottomBeam3, carbonFiberLadderStruct);
  bottomBeam3Vol->SetLineColor(fColorCarbonFiber);
  //bottomBeam3Vol->SetLineColor(2);
  TGeoRotation *bottomBeamRot1 = new TGeoRotation("", 90, 90,  90);
  TGeoRotation *bottomBeamRot2 = new TGeoRotation("",-90, 90, -90);

  TGeoCombiTrans *bottomBeamTransf1 = new TGeoCombiTrans
    (0,-(fgkLadderHeight/2-fgkLadderBeamRadius)-dy,0, bottomBeamRot1);
  TGeoCombiTrans *bottomBeamTransf2 = new TGeoCombiTrans(0,
				      -(fgkLadderHeight/2-fgkLadderBeamRadius)-dy,
				      -fgkSegmentLength/2, bottomBeamRot1);
  TGeoCombiTrans *bottomBeamTransf3 = new TGeoCombiTrans(0,
                                      -(fgkLadderHeight/2 - fgkLadderBeamRadius)
				      - dy, fgkSegmentLength/2, bottomBeamRot2);
  // be careful for beams #3: when "reading" from -z to +z and 
  // from the bottom of the ladder, it should draw a Lambda, and not a V
  TGeoRotation *bottomBeamRot4 = new TGeoRotation("", -90, fgkBottomBeamAngle, -90);
  TGeoRotation *bottomBeamRot5 = new TGeoRotation("" ,-90,-fgkBottomBeamAngle, -90);
  TGeoCombiTrans *bottomBeamTransf4 = new TGeoCombiTrans
    (0,-(fgkLadderHeight/2-fgkLadderBeamRadius)-dy,-fgkSegmentLength/4,bottomBeamRot4);
  TGeoCombiTrans *bottomBeamTransf5 = new TGeoCombiTrans
    (0,-(fgkLadderHeight/2-fgkLadderBeamRadius)-dy,fgkSegmentLength/4, bottomBeamRot5);

  fLaddSegCommonVol[0] = cfLaddTopVol1;  fLaddSegCommonTr[0] = trTop1;
  fLaddSegCommonVol[1] = cfLaddTopVol2;  fLaddSegCommonTr[1] = trTop1;
  fLaddSegCommonVol[2] = cfLaddSideVol1; fLaddSegCommonTr[2] = ctSideR;
  fLaddSegCommonVol[3] = cfLaddSideVol1; fLaddSegCommonTr[3] = ctSideL;
  fLaddSegCommonVol[4] = cfLaddSideVol2; fLaddSegCommonTr[4] = ctSideR;
  fLaddSegCommonVol[5] = cfLaddSideVol2; fLaddSegCommonTr[5] = ctSideL;
  fLaddSegCommonVol[6] = sideBeam;       fLaddSegCommonTr[6] = beamTransf[0];
  fLaddSegCommonVol[7] = sideBeam;       fLaddSegCommonTr[7] = beamTransf[1];
  fLaddSegCommonVol[8] = sideBeam;       fLaddSegCommonTr[8] = beamTransf[2];
  fLaddSegCommonVol[9] = sideBeam;       fLaddSegCommonTr[9] = beamTransf[3];
  fLaddSegCommonVol[10]= sideBeam;       fLaddSegCommonTr[10]= beamTransf[4];
  fLaddSegCommonVol[11]= sideBeam;       fLaddSegCommonTr[11]= beamTransf[5];
  fLaddSegCommonVol[12]= sideBeam;       fLaddSegCommonTr[12]= beamTransf[6];
  fLaddSegCommonVol[13]= sideBeam;       fLaddSegCommonTr[13]= beamTransf[7];
  fLaddSegCommonVol[14]= bottomBeam1Vol; fLaddSegCommonTr[14]= bottomBeamTransf1;
  fLaddSegCommonVol[15]= bottomBeam2Vol; fLaddSegCommonTr[15]= bottomBeamTransf2;
  fLaddSegCommonVol[16]= bottomBeam2Vol; fLaddSegCommonTr[16]= bottomBeamTransf3;
  fLaddSegCommonVol[17]= bottomBeam3Vol; fLaddSegCommonTr[17]= bottomBeamTransf4;
  fLaddSegCommonVol[18]= bottomBeam3Vol; fLaddSegCommonTr[18]= bottomBeamTransf5;

 
  //********************************************************************
  // cables
  //********************************************************************
  char cableName[30];
  for (Int_t i=0; i<fgkLay3Ndet; i++) {
    sprintf(cableName, "digitCableLay3A_%i",i);
    fDigitCableLay3A[i].SetName(cableName);
    fDigitCableLay3A[i].SetWidth(fgkDigitCablWidth);
    fDigitCableLay3A[i].SetThickness(fgkDigitCablPolyThick+fgkDigitCablAlThick);
    fDigitCableLay3A[i].SetNLayers(2);
    fDigitCableLay3A[i].SetLayer( 0, fgkDigitCablPolyThick, polyhamideSDD,
				  fColorPolyhamide);
    fDigitCableLay3A[i].SetLayer(1, fgkDigitCablAlThick, alSDD, fColorAl);
    sprintf(cableName, "digitCableLay3B_%i",i);
    fDigitCableLay3B[i].SetName(cableName);
    fDigitCableLay3B[i].SetWidth(fgkDigitCablWidth);
    fDigitCableLay3B[i].SetThickness(fgkDigitCablPolyThick+fgkDigitCablAlThick);
    fDigitCableLay3B[i].SetNLayers(2);
    fDigitCableLay3B[i].SetLayer( 0, fgkDigitCablPolyThick, polyhamideSDD,
				  fColorPolyhamide);
    fDigitCableLay3B[i].SetLayer(1, fgkDigitCablAlThick, alSDD, fColorAl);
  };
  for (Int_t i=0; i<fgkLay4Ndet; i++) {
    sprintf(cableName, "digitCableLay4A_%i",i);
    fDigitCableLay4A[i].SetName(cableName);
    fDigitCableLay4A[i].SetWidth(fgkDigitCablWidth);
    fDigitCableLay4A[i].SetThickness(fgkDigitCablPolyThick+fgkDigitCablAlThick);
    fDigitCableLay4A[i].SetNLayers(2);
    fDigitCableLay4A[i].SetLayer( 0, fgkDigitCablPolyThick, polyhamideSDD,
				  fColorPolyhamide);
    fDigitCableLay4A[i].SetLayer(1, fgkDigitCablAlThick, alSDD, fColorAl);
    sprintf(cableName, "digitCableLay4B_%i",i);
    fDigitCableLay4B[i].SetName(cableName);
    fDigitCableLay4B[i].SetWidth(fgkDigitCablWidth);
    fDigitCableLay4B[i].SetThickness(fgkDigitCablPolyThick+fgkDigitCablAlThick);
    fDigitCableLay4B[i].SetNLayers(2);
    fDigitCableLay4B[i].SetLayer( 0, fgkDigitCablPolyThick, polyhamideSDD,
				  fColorPolyhamide);
    fDigitCableLay4B[i].SetLayer(1, fgkDigitCablAlThick, alSDD, fColorAl);
  };
                                        // Well, those digit cables could also include the analog cables
                                        // which have the same width and the same path, at least in the ladder.
                                        // It will gain some computing ressources (less volumes) and some
                                        // coding efforts ... !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                        // The only thing to do is to change the names and put the correct total
                                        // thicknesses

  // some transformations and volumes used in several places
  fCommonTr[0] = new TGeoRotation("CarlosSuppRotN",
				  0, -fgkCarlosSuppAngle, 0);

  TGeoTube *littleScrewHead = new TGeoTube("littleScrewHead", 0, fgkLittleScrewHeadR,
					   fgkLittleScrewHeadH/2);
  fCommonVol[0] = new TGeoVolume("vLittleScrewHead",
				 littleScrewHead, stainless);
  fCommonVol[0]->SetLineColor(kGray);

  fLadderFoot = CreateLadderFoot();
  CreateLVCard();
  fCardHV     = CreateHVCard(0);
  fCardCarlos = CreateCarlosCard(0);

  //==================
  // link beteen phynox and plastic cooling tubes
  //==================

  fRaccordoL = new TGeoVolumeAssembly("RaccordoL");
  Double_t fullRacLen = fgkConnectorCoolTubeL1+fgkConnectorCoolTubeL2+fgkConnectorCoolTubeL3;
  TGeoTube *waterRac = new TGeoTube("waterRac", 0, fgkConnectorCoolTubeRmin, fullRacLen/2);
  TGeoVolume * vwaterRac = new TGeoVolume("vwaterRac", waterRac, coolerMediumSDD);
  vwaterRac->SetLineColor(kBlue);

  TGeoTube *tube1Rac = new TGeoTube("tube1Rac", fgkConnectorCoolTubeRmin,
				    fgkConnectorCoolTubeR1, fgkConnectorCoolTubeL1/2);
  TGeoTube *tube2Rac = new TGeoTube("tube2Rac", fgkConnectorCoolTubeRmin,
				    fgkConnectorCoolTubeR2, fgkConnectorCoolTubeL2/2);
  TGeoTube *tube3Rac = new TGeoTube("tube3Rac", fgkConnectorCoolTubeRmin,
				    fgkConnectorCoolTubeR3, fgkConnectorCoolTubeL3/2);
  TGeoVolume *  vtube1Rac = new TGeoVolume("vtube1Rac", tube1Rac, raccordMedium);
  TGeoVolume *  vtube2Rac = new TGeoVolume("vtube2Rac", tube2Rac, raccordMedium);
  TGeoVolume *  vtube3Rac = new TGeoVolume("vtube3Rac", tube3Rac, raccordMedium);
  vtube1Rac->SetLineColor(kGray);
  vtube2Rac->SetLineColor(kGray);
  vtube3Rac->SetLineColor(kGray);

  TGeoTranslation *trTube1Rac = new TGeoTranslation("trTube1Rac",0,0,
						    -fullRacLen/2+fgkConnectorCoolTubeL1/2);
  TGeoTranslation *trTube2Rac = new TGeoTranslation("trTube2Rac",0,0,
			       (-fullRacLen/2+fgkConnectorCoolTubeL1+fgkConnectorCoolTubeL2/2));
  TGeoTranslation *trTube3Rac = new TGeoTranslation("trTube3Rac",0,0,
				    (-fullRacLen/2+fgkConnectorCoolTubeL1+
				     fgkConnectorCoolTubeL2+fgkConnectorCoolTubeL3/2));
  fRaccordoL->AddNode(vwaterRac, 1,0);
  fRaccordoL->AddNode(vtube1Rac, 1,trTube1Rac);
  fRaccordoL->AddNode(vtube2Rac, 1,trTube2Rac);
  fRaccordoL->AddNode(vtube3Rac, 1,trTube3Rac);
}


//________________________________________________________________________
void AliITSv11GeometrySDD::CheckOverlaps(Double_t precision){
  //
  // a debugging function for checking some possible overlaps
  //
  if (fSDDsensor3)        fSDDsensor3->CheckOverlaps(precision);
  if (fSDDsensor4)        fSDDsensor4->CheckOverlaps(precision);
  if (fHybrid)            fHybrid->CheckOverlaps(precision);
}


//________________________________________________________________________
TGeoCombiTrans *AliITSv11GeometrySDD::
CreateCombiTrans(const char *name, Double_t dy, Double_t dz, Double_t dphi,
		 Bool_t planeSym) {
    //
    // return the TGeoCombiTrans which make a translation in y and z
    // and a rotation in phi in the global coord system
    // If planeSym = true, the rotation places the object symetrically
    // (with respect to the transverse plane) to its position in the
    // case planeSym = false
    //

    TGeoTranslation t1(dy*CosD(90.+dphi),dy*SinD(90.+dphi), dz);
    TGeoRotation r1("",0.,0.,dphi);
    TGeoRotation r2("",90, 180, -90-dphi);

    TGeoCombiTrans *combiTrans1 = new TGeoCombiTrans(name);
    combiTrans1->SetTranslation(t1);
    if (planeSym) combiTrans1->SetRotation(r1);
    else  combiTrans1->SetRotation(r2);
    return combiTrans1;
}


//________________________________________________________________________
void AliITSv11GeometrySDD::AddTranslationToCombiTrans(TGeoCombiTrans* ct,
                                                       Double_t dx,
                                                       Double_t dy,
                                                       Double_t dz) const{
  // Add a dx,dy,dz translation to the initial TGeoCombiTrans
  const Double_t *vect = ct->GetTranslation();
  Double_t newVect[3] = {vect[0]+dx, vect[1]+dy, vect[2]+dz};
  ct->SetTranslation(newVect);
}


//________________________________________________________________________
void AliITSv11GeometrySDD::ShowOnePiece(TGeoVolume *moth) {
// for code developpment and debugging purposes

  if (! fSDDsensor3) CreateBasicObjects();

  //  moth->AddNode(fPinSupport, 1, 0);
  //  moth->AddNode(fCoolPipeSupportL, 1, 0);
  //  moth->AddNode(fSDDsensor3, 1, 0);
  //  moth->AddNode(fSDDsensor4, 1, 0);
  //  moth->AddNode(fBaseThermalBridge, 1, 0);
  //  moth->AddNode(fHybrid,100,0);
  //  moth->AddNode(fLadderFoot,1,0);
  //moth->AddNode(fCardLVL,1,0);
  //moth->AddNode(fCardLVR,1,0);

   TGeoVolume* seg = CreateLadderSegment( 3, 0);
   moth->AddNode(seg, 1, 0);

//   TGeoVolumeAssembly *lay3Ladder = CreateLadder(3);
//   moth->AddNode(lay3Ladder, 1, 0);

//   TGeoVolumeAssembly *lay3Detectors = CreateDetectorsAssembly(3);
//   moth->AddNode(lay3Detectors, 1, 0);

//   TGeoVolumeAssembly *lay3Detectors = CreateDetectorsAssembly(3);
//   moth->AddNode(lay3Detectors, 1, 0);


//   TGeoVolumeAssembly *endLadder = CreateEndLadder( 4 );
//   moth->AddNode(endLadder, 1, 0);

//   TGeoVolumeAssembly *highVCard = CreateHVCard( 4 );
//   moth->AddNode(highVCard, 1, 0);

//   TGeoVolumeAssembly *supportRing = CreateSupportRing( 4 );
//   moth->AddNode(supportRing, 1, 0);

//   TGeoVolume *endLadderCards = CreateEndLadderCardsV( 4 );
//   moth->AddNode(endLadderCards, 1, 0);

//   TGeoVolumeAssembly *carlosCard = CreateCarlosCard( 4 );
//   moth->AddNode(carlosCard, 1, 0);



  /*
  //==================================
  //--- test of flat cable curvature
  //==================================

  double angle = 90;
  AliITSv11GeomCableFlat cable("test", 3, 0.3);
  cable.SetNLayers(1);
  cable.SetNLayers(2);
  cable.SetLayer(0, 0.2, coolerMediumSDD, 2);
  cable.SetLayer(1, 0.1, coolerMediumSDD, 3);
  cable.SetInitialNode(endLadderCards);

  Double_t p1[3], p2[3], vX[3] = {1,0,0},vY[3] = {0,5,0};

  p1[0] = -3;
  p1[1] = 1;
  p1[2] = 10;

  p2[0] = 0;
  p2[1] = 1;
  p2[2] = 10;
  cable.AddCheckPoint(endLadderCards, 0, p1, vX);
  cable.AddCheckPoint(endLadderCards, 1, p2, vX);
  cable.CreateAndInsertBoxCableSegment(1,angle);

  Double_t p3[3], p4[3];

  p3[0] = 2;
  p3[1] = 3;
  p3[2] = 10;
  cable.AddCheckPoint(endLadderCards, 2, p3, vY);
  cable.CreateAndInsertCableCylSegment(2,angle);

  p4[0] = 2;
  p4[1] = 6;
  p4[2] = 10;
  cable.AddCheckPoint(endLadderCards, 3, p4, vY);
  cable.CreateAndInsertCableSegment(3,angle);
  */
}


//________________________________________________________________________
void AliITSv11GeometrySDD::Layer3(TGeoVolume *moth) {
  //
  // Insert the layer 3 in the mother volume. This is a virtual volume
  // containing ladders of layer 3 and the supporting rings
  //

  if (! moth) {
    printf("Error::AliITSv11GeometrySDD: Can't insert layer3, mother is null!\n");
    return;
  };

  TGeoMedium *airSDD = GetMedium("SDD AIR$");

  fMotherVol = moth;
  if (! fSDDsensor3) CreateBasicObjects();
  
  //====================================
  // First we create the central barrel
  //====================================

  TGeoVolumeAssembly *lay3Ladder    = CreateLadder(3);
  TGeoVolumeAssembly *lay3Detectors = CreateDetectorsAssembly(3);
  TGeoVolumeAssembly *lay3Ladd2Det  = CreateDetectorsAssemblyLadd2();
  //TGeoVolume *lay3Detectors = CreateDetectors(3);
  TGeoTube *virtualLayer3Shape = new TGeoTube("ITSsddLayer3Shape",
				     fgkLay3Rmin,fgkLay3Rmax,fgkLay3Length*0.5);
  TGeoVolume *virtualLayer3 = new TGeoVolume("ITSsddLayer3",
					     virtualLayer3Shape, airSDD);

  Double_t dPhi = 360./fgkLay3Nladd;
  Double_t detectorsThick = fgkLadWaferSep + 2*fgkWaferThickness;
  // Placing virtual ladder and detectors volumes following
  // ladder ordering convention
  char rotName[30];
  Int_t iLaddMin = 0;
  Int_t iLaddMax = fgkLay3Nladd;
  if ((fAddOnlyLadder3min>=0)&&(fAddOnlyLadder3max<fgkLay3Nladd)) {
    iLaddMin = fAddOnlyLadder3min;
    iLaddMax = fAddOnlyLadder3max+1;
  };

  for (Int_t iLadd = iLaddMin; iLadd < iLaddMax; iLadd++) {

    Double_t ladderPhi = -3*dPhi+iLadd*dPhi;
    sprintf(rotName, "ITSsddLay3Ladd%i",iLadd);
    Double_t minRadiusLadBox = fLay3LaddShortRadius-fLay3LadderUnderSegDH;
    if (iLadd%2 != 0)
      minRadiusLadBox = fLay3LaddLongRadius-fLay3LadderUnderSegDH;
    minRadiusLadBox += ((TGeoBBox*)lay3Ladder->GetShape())->GetDY();
    TGeoCombiTrans *ctLadd;
    //=============================================================
    //
    //   Special modification  for ladder 2 of layer 3:
    //   It has been inverted (pi rotation around y axis)
    //
    //=============================================================
    if (iLadd != 2)
      ctLadd = CreateCombiTrans(rotName,minRadiusLadBox,
				0, ladderPhi, kTRUE);
    else
      ctLadd = CreateCombiTrans(rotName,minRadiusLadBox,
				0, ladderPhi, kFALSE);
    virtualLayer3->AddNode(lay3Ladder, iLadd, ctLadd);
    ///////////////////////////////////////////////////
    sprintf(rotName, "ITSsddLay3DetBox%i",iLadd);
    Double_t minRadiusDetBox = fgkLay3DetShortRadius;
    if (iLadd%2 != 0) minRadiusDetBox = fgkLay3DetLongRadius;
    minRadiusDetBox += detectorsThick/2;
    TGeoCombiTrans *ctDet;
    ctDet = CreateCombiTrans(rotName, minRadiusDetBox,
			     0, ladderPhi, kTRUE);

    if (iLadd != 2)
      virtualLayer3->AddNode(lay3Detectors, iLadd, ctDet);
    else
      virtualLayer3->AddNode(lay3Ladd2Det , iLadd, ctDet);

    ///////////////////////////////////////////////////
  }

  /*
  //====================================
  // Then the forward rapidity pieces
  // (cooling, Carlos, LV, HV ...)
  //====================================

  Double_t fgkForwardLay3Length = fgkEndLadPipeUlengthLay3+10*fgkmm;  // this has to be tune
  Double_t fgkForwardLay3Rmin = fgkLay3Rmin-7*fgkmm;
  Double_t fgkForwardLay3Rmax = fgkLay3Rmax-5*fgkmm;

  TGeoVolumeAssembly* lay3EndLadder = CreateEndLadderCards(3);
  TGeoTube *virtualForward3Shape = new TGeoTube("virtualForward3Shape",
						fgkForwardLay3Rmin, fgkForwardLay3Rmax,
						fgkForwardLay3Length/2.);

//   TGeoPcon *virtualForward3Shape = new TGeoPcon("virtualForward3Shape",0,360,2);
// // virtualForward3Shape->DefineSection(Int_t snum, Double_t z, Double_t rmin, Double_t rmax);
//   virtualForward3Shape->DefineSection(0, Double_t z, Double_t rmin, Double_t rmax);


  TGeoVolume *virtualForward3Pos = new TGeoVolume("ITSsddForward3Pos",
						  virtualForward3Shape, airSDD);
  TGeoVolume *virtualForward3Neg = new TGeoVolume("ITSsddForward3Neg",
						  virtualForward3Shape, airSDD);
//   TGeoVolume *virtualForward3Neg = new TGeoVolumeAssembly("ITSsddForward3Neg");
//   TGeoVolume *virtualForward3Pos = new TGeoVolumeAssembly("ITSsddForward3Pos");

  TGeoTranslation *virtualForward3TrPos = new TGeoTranslation("virtualForward3TrPos",0,0,
					  fgkLay3Length/2+fgkForwardLay3Length/2);
  TGeoTranslation *virtualForward3TrNeg = new TGeoTranslation("virtualForward3TrNeg",0,0,
					  -fgkLay3Length/2-fgkForwardLay3Length/2);

  for (Int_t iLadd = iLaddMin; iLadd < iLaddMax; iLadd++) {

    Double_t ladderPhi = -3*dPhi+iLadd*dPhi;
    Double_t minRadiusDetBox = fgkLay3DetShortRadius;
    if (iLadd%2 != 0) minRadiusDetBox = fgkLay3DetLongRadius;
    minRadiusDetBox += detectorsThick/2;

    sprintf(rotName, "ITSsddLay3EndLadd%i",iLadd);

    TGeoCombiTrans *ctEndLaddPos = CreateCombiTrans(rotName, minRadiusDetBox-0.1,
				   -fgkForwardLay3Length/2, ladderPhi, kTRUE);
    TGeoCombiTrans *ctEndLaddNeg = CreateCombiTrans(rotName, minRadiusDetBox-0.1,
				   fgkForwardLay3Length/2, ladderPhi, kFALSE);

    virtualForward3Pos->AddNode(lay3EndLadder, iLadd*2, ctEndLaddPos);
    virtualForward3Neg->AddNode(lay3EndLadder, iLadd*2, ctEndLaddNeg);
  }

  */


  if(GetDebug(1)) {
    virtualLayer3->CheckOverlaps(0.01);
    //virtualForward3Pos->CheckOverlaps(0.01);
    //virtualForward3Neg->CheckOverlaps(0.01);
  }

  virtualLayer3->SetVisibility(kFALSE);
  //virtualForward3Pos->SetVisibility(kFALSE);
  //virtualForward3Neg->SetVisibility(kFALSE);


  moth->AddNode(virtualLayer3, 1, 0);
  //moth->AddNode(virtualForward3Pos, 1, virtualForward3TrPos);
  //moth->AddNode(virtualForward3Neg, 1, virtualForward3TrNeg);
}


// //________________________________________________________________________
// void AliITSv11GeometrySDD::ForwardLayer3(TGeoVolume *moth) {
//   //
//   // Insert the forward pieces of layer 3 in the mother volume. 
//   // (cooling, Carlos, LV, HV ...)
//   //

//   if (! moth) {
//     printf("Error::AliITSv11GeometrySDD: Can't insert layer3, mother is null!\n");
//     return;
//   };

//   TGeoMedium *airSDD = GetMedium("SDD AIR$");

//   if (! fSDDsensor3) CreateBasicObjects();

//   Double_t dPhi = 360./fgkLay3Nladd;
//   Double_t detectorsThick = fgkLadWaferSep + 2*fgkWaferThickness;
//   Int_t iLaddMin = 0;
//   Int_t iLaddMax = fgkLay3Nladd;
//   if ((fAddOnlyLadder3min>=0)&&(fAddOnlyLadder3max<fgkLay3Nladd)) {
//     iLaddMin = fAddOnlyLadder3min;
//     iLaddMax = fAddOnlyLadder3max+1;
//   };
//   char rotName[30];


//   //=================

//   Double_t fgkForwardLay3Length = fgkEndLadPipeUlengthLay3+10*fgkmm;  // this has to be tune
//   Double_t fgkForwardLay3Rmin = fgkLay3Rmin-7*fgkmm;
//   Double_t fgkForwardLay3Rmax = fgkLay3Rmax-5*fgkmm;

//   TGeoVolumeAssembly* lay3EndLadder = CreateEndLadderCards(3);
//   TGeoTube *virtualForward3Shape = new TGeoTube("virtualForward3Shape",
// 						fgkForwardLay3Rmin, fgkForwardLay3Rmax,
// 						fgkForwardLay3Length/2.);

// //   TGeoPcon *virtualForward3Shape = new TGeoPcon("virtualForward3Shape",0,360,2);
// // // virtualForward3Shape->DefineSection(Int_t snum, Double_t z, Double_t rmin, Double_t rmax);
// //   virtualForward3Shape->DefineSection(0, Double_t z, Double_t rmin, Double_t rmax);


//   TGeoVolume *virtualForward3Pos = new TGeoVolume("ITSsddForward3Pos",
// 						  virtualForward3Shape, airSDD);
//   TGeoVolume *virtualForward3Neg = new TGeoVolume("ITSsddForward3Neg",
// 						  virtualForward3Shape, airSDD);
// //   TGeoVolume *virtualForward3Neg = new TGeoVolumeAssembly("ITSsddForward3Neg");
// //   TGeoVolume *virtualForward3Pos = new TGeoVolumeAssembly("ITSsddForward3Pos");

//   TGeoTranslation *virtualForward3TrPos = new TGeoTranslation("virtualForward3TrPos",0,0,
// 					  fgkLay3Length/2+fgkForwardLay3Length/2);
//   TGeoTranslation *virtualForward3TrNeg = new TGeoTranslation("virtualForward3TrNeg",0,0,
// 					  -fgkLay3Length/2-fgkForwardLay3Length/2);

//   for (Int_t iLadd = iLaddMin; iLadd < iLaddMax; iLadd++) {

//     Double_t ladderPhi = -3*dPhi+iLadd*dPhi;
//     Double_t minRadiusDetBox = fgkLay3DetShortRadius;
//     if (iLadd%2 != 0) minRadiusDetBox = fgkLay3DetLongRadius;
//     minRadiusDetBox += detectorsThick/2;

//     sprintf(rotName, "ITSsddLay3EndLadd%i",iLadd);

//     TGeoCombiTrans *ctEndLaddPos = CreateCombiTrans(rotName, minRadiusDetBox-0.1,
// 				   -fgkForwardLay3Length/2, ladderPhi, kTRUE);
//     TGeoCombiTrans *ctEndLaddNeg = CreateCombiTrans(rotName, minRadiusDetBox-0.1,
// 				   fgkForwardLay3Length/2, ladderPhi, kFALSE);

//     virtualForward3Pos->AddNode(lay3EndLadder, iLadd*2, ctEndLaddPos);
//     virtualForward3Neg->AddNode(lay3EndLadder, iLadd*2, ctEndLaddNeg);
//   }

//   if(GetDebug(1)) {
//     virtualForward3Pos->CheckOverlaps(0.01);
//     virtualForward3Neg->CheckOverlaps(0.01);
//   }

//   virtualForward3Pos->SetVisibility(kFALSE);
//   virtualForward3Neg->SetVisibility(kFALSE);

//   moth->AddNode(virtualForward3Pos, 1, virtualForward3TrPos);
//   moth->AddNode(virtualForward3Neg, 1, virtualForward3TrNeg);
// }



//________________________________________________________________________
void AliITSv11GeometrySDD::ForwardLayer3(TGeoVolume *moth) {
  //
  // Insert the end-ladder of layer 3 in the mother volume. 
  // (cooling, Carlos, LV, HV ...)
  //

  if (! moth) {
    printf("Error::AliITSv11GeometrySDD: Can't insert layer3, mother is null!\n");
    return;
  };

  if (! fSDDsensor3) CreateBasicObjects();

  Int_t iLaddMin = 0;
  Int_t iLaddMax = fgkLay3Nladd;
  if ((fAddOnlyLadder3min>=0)&&(fAddOnlyLadder3max<fgkLay3Nladd)) {
    iLaddMin = fAddOnlyLadder3min;
    iLaddMax = fAddOnlyLadder3max+1;
  };

  TGeoVolume *virtualForward3Neg = new TGeoVolumeAssembly("ITSsddForward3Neg");
  TGeoVolume *virtualForward3Pos = new TGeoVolumeAssembly("ITSsddForward3Pos");

  char rotName[30];
  Double_t dPhi = 360./fgkLay3Nladd;
  TGeoVolume* lay3EndLadder = CreateEndLadderCardsV(3);

  for (Int_t iLadd = iLaddMin; iLadd < iLaddMax; iLadd++) {

    Double_t ladderPhi = -3*dPhi+iLadd*dPhi;
    Double_t dR = 0;
    if (iLadd%2 != 0) dR = fgkLay3DetLongRadius-fgkLay3DetShortRadius;

    sprintf(rotName, "ITSsddLay3EndLadd%i",iLadd);

    TGeoCombiTrans *ctEndLaddPos = CreateCombiTrans(rotName, dR,
				    fgkLay3Length/2, ladderPhi, kTRUE);
    TGeoCombiTrans *ctEndLaddNeg = CreateCombiTrans(rotName, dR,
				   -fgkLay3Length/2, ladderPhi, kFALSE);

    virtualForward3Pos->AddNode(lay3EndLadder, iLadd*2, ctEndLaddPos);
    virtualForward3Neg->AddNode(lay3EndLadder, iLadd*2, ctEndLaddNeg);
  }

  if(GetDebug(1)) {
    virtualForward3Pos->CheckOverlaps(0.01);
    virtualForward3Neg->CheckOverlaps(0.01);
  }

  // 180deg Y rotation to compensate the cancellation of ITSD volume
  // (idortm[199] in AliITSv11Hybrid : z--->  -z;   x ---> -x;   y ---> y)
  TGeoRotation *y180 = new TGeoRotation();
  y180->SetAngles( 90.,180., 90., 90.,180.,  0.);
  moth->AddNode(virtualForward3Pos, 1, y180);
  moth->AddNode(virtualForward3Neg, 1, y180);
}

//________________________________________________________________________
void AliITSv11GeometrySDD::Layer4(TGeoVolume *moth) {
  //
  // Insert the layer 4 in the mother volume. This is a virtual volume
  // containing ladders of layer 4 and the supporting rings
  //

  if (! moth) {
    printf("Error::AliITSv11GeometrySDD: Can't insert layer4, mother is null!\n");
    return;
  };

  fMotherVol = moth;

  if (! fSDDsensor3) CreateBasicObjects();

  TGeoTube *virtualLayer4Shape =new TGeoTube("ITSsddLayer4Shape",
				    fgkLay4Rmin,fgkLay4Rmax,fgkLay4Length*0.5);
  TGeoMedium *airSDD = GetMedium("SDD AIR$");
  TGeoVolume *virtualLayer4 = new TGeoVolume("ITSsddLayer4",
					     virtualLayer4Shape, airSDD);

  //====================================
  // First we create the central barrel
  //====================================

   TGeoVolumeAssembly *lay4Ladder    = CreateLadder(4);
  //TGeoVolume *lay4Detectors = CreateDetectors(4);
  TGeoVolumeAssembly *lay4Detectors = CreateDetectorsAssembly(4);

  Double_t dPhi = 360./fgkLay4Nladd;
  Double_t detBoxThickness = fgkLadWaferSep + 2*fgkWaferThickness;

  // placing virtual ladder and detectors volumes following ladder 
  // ordering convention
  char rotName[20];
  Int_t iLaddMin = 0;
  Int_t iLaddMax = fgkLay4Nladd;
  if ((fAddOnlyLadder4min >= 0)&&(fAddOnlyLadder4max < fgkLay4Nladd)) {
    iLaddMin = fAddOnlyLadder4min;
    iLaddMax = fAddOnlyLadder4max+1;
  }
  for (Int_t iLadd = iLaddMin; iLadd < iLaddMax; iLadd++) {

    Double_t ladderPhi = -5*dPhi + iLadd*dPhi;
    sprintf(rotName, "ITSsddLay4Ladd%i",iLadd);
    Double_t minRadiusLadBox = fLay4LaddShortRadius-fLay4LadderUnderSegDH;
    if (iLadd%2 != 0)
      minRadiusLadBox = fLay4LaddLongRadius-fLay4LadderUnderSegDH;
    minRadiusLadBox += ((TGeoBBox*)lay4Ladder->GetShape())->GetDY();
    TGeoCombiTrans *ctLadd = CreateCombiTrans(rotName, minRadiusLadBox,
					      0, ladderPhi, kTRUE);
    virtualLayer4->AddNode(lay4Ladder, iLadd, ctLadd);
    ///////////////////////////////////////////////////
    sprintf(rotName, "ITSsddLay4DetBox%i",iLadd);
    Double_t minRadiusDetBox = fgkLay4DetShortRadius;
    if (iLadd%2 != 0)
      minRadiusDetBox = fgkLay4DetLongRadius;
    minRadiusDetBox += detBoxThickness/2;
    TGeoCombiTrans *ctDet = CreateCombiTrans(rotName, minRadiusDetBox,
					     0, ladderPhi, kTRUE);
    virtualLayer4->AddNode(lay4Detectors, iLadd, ctDet);
    ///////////////////////////////////////////////////
  }

  /*
  //====================================
  // Then the pieces at forward rapidity
  // (cooling, Carlos, LV, HV ...)
  //====================================

  Double_t fgkForwardLay4Length = fgkEndLadPipeUlengthLay4+10*fgkmm;  // this has to be tuned
  Double_t fgkForwardLay4Rmin = fgkLay4Rmin-9*fgkmm;
  Double_t fgkForwardLay4Rmax = fgkLay4Rmax-5*fgkmm;

  TGeoVolumeAssembly* lay4EndLadder = CreateEndLadderCards(4);
  TGeoTube *virtualForward4Shape = new TGeoTube("virtualForward3Shape",
						fgkForwardLay4Rmin, fgkForwardLay4Rmax,
						fgkForwardLay4Length/2.);
  TGeoVolume *virtualForward4Pos = new TGeoVolume("ITSsddForward4Pos",
						  virtualForward4Shape, airSDD);
  TGeoVolume *virtualForward4Neg = new TGeoVolume("ITSsddForward4Neg",
						  virtualForward4Shape, airSDD);
//   TGeoVolume *virtualForward4Pos = new TGeoVolumeAssembly("ITSsddForward4Pos");
//   TGeoVolume *virtualForward4Neg = new TGeoVolumeAssembly("ITSsddForward4Neg");

  TGeoTranslation *virtualForward4TrPos = new TGeoTranslation("virtualForward4TrPos",0,0,
					  fgkLay4Length/2+fgkForwardLay4Length/2);
  TGeoTranslation *virtualForward4TrNeg = new TGeoTranslation("virtualForward4TrNeg",0,0,
					  -fgkLay4Length/2-fgkForwardLay4Length/2);

  for (Int_t iLadd = iLaddMin; iLadd < iLaddMax; iLadd++) {

    Double_t ladderPhi = -5*dPhi + iLadd*dPhi;
    Double_t minRadiusDetBox = fgkLay4DetShortRadius;
    if (iLadd%2 != 0)
      minRadiusDetBox = fgkLay4DetLongRadius;
    minRadiusDetBox += detBoxThickness/2;

    sprintf(rotName, "ITSsddLay4EndLadd%i",iLadd);

    TGeoCombiTrans *ctEndLaddPos = CreateCombiTrans(rotName, minRadiusDetBox-0.1,
				   -fgkForwardLay4Length/2, ladderPhi, kTRUE);
    TGeoCombiTrans *ctEndLaddNeg = CreateCombiTrans(rotName, minRadiusDetBox-0.1,
				    fgkForwardLay4Length/2, ladderPhi, kFALSE);
    virtualForward4Pos->AddNode(lay4EndLadder, iLadd*2, ctEndLaddPos);
    virtualForward4Neg->AddNode(lay4EndLadder, iLadd*2+1, ctEndLaddNeg);
  }
  */

  if(GetDebug(1)) virtualLayer4->CheckOverlaps(0.01);

  virtualLayer4->SetVisibility(kFALSE);
  //virtualForward4Pos->SetVisibility(kFALSE);
  //virtualForward4Neg->SetVisibility(kFALSE);

  moth->AddNode(virtualLayer4,1,0);
  //moth->AddNode(virtualForward4Pos, 1, virtualForward4TrPos);
  //moth->AddNode(virtualForward4Neg, 1, virtualForward4TrNeg);
}


// //________________________________________________________________________
// void AliITSv11GeometrySDD::ForwardLayer4(TGeoVolume *moth) {
//   //
//   // Insert the layer 4 in the mother volume. This is a virtual volume
//   // containing ladders of layer 4 and the supporting rings
//   // (cooling, Carlos, LV, HV ...)
//   //

//   if (! moth) {
//     printf("Error::AliITSv11GeometrySDD: Can't insert layer4, mother is null!\n");
//     return;
//   };

//   TGeoMedium *airSDD = GetMedium("SDD AIR$");

//   if (! fSDDsensor3) CreateBasicObjects();

//   Double_t dPhi = 360./fgkLay4Nladd;
//   Double_t detBoxThickness = fgkLadWaferSep + 2*fgkWaferThickness;

//   // placing virtual ladder and detectors volumes following ladder 
//   // ordering convention
//   char rotName[20];
//   Int_t iLaddMin = 0;
//   Int_t iLaddMax = fgkLay4Nladd;
//   if ((fAddOnlyLadder4min >= 0)&&(fAddOnlyLadder4max < fgkLay4Nladd)) {
//     iLaddMin = fAddOnlyLadder4min;
//     iLaddMax = fAddOnlyLadder4max+1;
//   }

//   //=================
//   Double_t fgkForwardLay4Length = fgkEndLadPipeUlengthLay4+10*fgkmm;  // this has to be tuned
//   Double_t fgkForwardLay4Rmin = fgkLay4Rmin-9*fgkmm;
//   Double_t fgkForwardLay4Rmax = fgkLay4Rmax-5*fgkmm;

//   TGeoVolumeAssembly* lay4EndLadder = CreateEndLadderCards(4);
//   TGeoTube *virtualForward4Shape = new TGeoTube("virtualForward3Shape",
// 						fgkForwardLay4Rmin, fgkForwardLay4Rmax,
// 						fgkForwardLay4Length/2.);
//   TGeoVolume *virtualForward4Pos = new TGeoVolume("ITSsddForward4Pos",
// 						  virtualForward4Shape, airSDD);
//   TGeoVolume *virtualForward4Neg = new TGeoVolume("ITSsddForward4Neg",
// 						  virtualForward4Shape, airSDD);
// //   TGeoVolume *virtualForward4Pos = new TGeoVolumeAssembly("ITSsddForward4Pos");
// //   TGeoVolume *virtualForward4Neg = new TGeoVolumeAssembly("ITSsddForward4Neg");

//   TGeoTranslation *virtualForward4TrPos = new TGeoTranslation("virtualForward4TrPos",0,0,
// 					  fgkLay4Length/2+fgkForwardLay4Length/2);
//   TGeoTranslation *virtualForward4TrNeg = new TGeoTranslation("virtualForward4TrNeg",0,0,
// 					  -fgkLay4Length/2-fgkForwardLay4Length/2);

//   for (Int_t iLadd = iLaddMin; iLadd < iLaddMax; iLadd++) {

//     Double_t ladderPhi = -5*dPhi + iLadd*dPhi;
//     Double_t minRadiusDetBox = fgkLay4DetShortRadius;
//     if (iLadd%2 != 0)
//       minRadiusDetBox = fgkLay4DetLongRadius;
//     minRadiusDetBox += detBoxThickness/2;

//     sprintf(rotName, "ITSsddLay4EndLadd%i",iLadd);

//     TGeoCombiTrans *ctEndLaddPos = CreateCombiTrans(rotName, minRadiusDetBox-0.1,
// 				   -fgkForwardLay4Length/2, ladderPhi, kTRUE);
//     TGeoCombiTrans *ctEndLaddNeg = CreateCombiTrans(rotName, minRadiusDetBox-0.1,
// 				    fgkForwardLay4Length/2, ladderPhi, kFALSE);
//     virtualForward4Pos->AddNode(lay4EndLadder, iLadd*2, ctEndLaddPos);
//     virtualForward4Neg->AddNode(lay4EndLadder, iLadd*2+1, ctEndLaddNeg);
//   }

//   virtualForward4Pos->SetVisibility(kFALSE);
//   virtualForward4Neg->SetVisibility(kFALSE);

//   moth->AddNode(virtualForward4Pos, 1, virtualForward4TrPos);
//   moth->AddNode(virtualForward4Neg, 1, virtualForward4TrNeg);
// }


//________________________________________________________________________
void AliITSv11GeometrySDD::ForwardLayer4(TGeoVolume *moth) {
  //
  // Insert the end-ladder of layer 4 in the mother volume.
  // (cooling, Carlos, LV, HV ...)
  //

  if (! moth) {
    printf("Error::AliITSv11GeometrySDD: Can't insert layer4, mother is null!\n");
    return;
  };

  if (! fSDDsensor3) CreateBasicObjects();

  // placing virtual ladder and detectors volumes following ladder 
  // ordering convention
  Int_t iLaddMin = 0;
  Int_t iLaddMax = fgkLay4Nladd;
  if ((fAddOnlyLadder4min >= 0)&&(fAddOnlyLadder4max < fgkLay4Nladd)) {
    iLaddMin = fAddOnlyLadder4min;
    iLaddMax = fAddOnlyLadder4max+1;
  }

  TGeoVolume *virtualForward4Pos = new TGeoVolumeAssembly("ITSsddForward4Pos");
  TGeoVolume *virtualForward4Neg = new TGeoVolumeAssembly("ITSsddForward4Neg");

  char rotName[30];
  Double_t dPhi = 360./fgkLay4Nladd;
  TGeoVolume* lay4EndLadder = CreateEndLadderCardsV(4);

  for (Int_t iLadd = iLaddMin; iLadd < iLaddMax; iLadd++) {

    Double_t ladderPhi = -5*dPhi + iLadd*dPhi;
    Double_t dR = 0;
    if (iLadd%2 != 0)
      dR = fgkLay4DetLongRadius-fgkLay4DetShortRadius;

    sprintf(rotName, "ITSsddLay4EndLadd%i",iLadd);

    TGeoCombiTrans *ctEndLaddPos = CreateCombiTrans(rotName, dR,
				   fgkLay4Length/2, ladderPhi, kTRUE);
    TGeoCombiTrans *ctEndLaddNeg = CreateCombiTrans(rotName, dR,
				   -fgkLay4Length/2, ladderPhi, kFALSE);
    virtualForward4Pos->AddNode(lay4EndLadder, iLadd*2, ctEndLaddPos);
    virtualForward4Neg->AddNode(lay4EndLadder, iLadd*2, ctEndLaddNeg);
  }

  // 180deg Y rotation to compensate the cancellation of ITSD volume
  // (idortm[199] in AliITSv11Hybrid : z--->  -z;   x ---> -x;   y ---> y)
  TGeoRotation *y180 = new TGeoRotation();
  y180->SetAngles( 90.,180., 90., 90.,180.,  0.);
  moth->AddNode(virtualForward4Pos, 1, y180);
  moth->AddNode(virtualForward4Neg, 1, y180);
}


//________________________________________________________________________
TGeoVolumeAssembly *AliITSv11GeometrySDD::CreateLadder(Int_t iLay) {
  //
  // return an assembly volume containing the CF ladder
  //

  Int_t    nDetectors   = fgkLay3Ndet;
  Double_t ladderLength = fgkLay3LadderLength;
  Double_t underSegDH   = fLay3LadderUnderSegDH;
  Double_t *sensorZPos  = fLay3sensorZPos;
  AliITSv11GeomCableFlat *digitCableA = fDigitCableLay3A;
  AliITSv11GeomCableFlat *digitCableB = fDigitCableLay3B;

  if (iLay==3) {}
  else if (iLay==4) {
    nDetectors   = fgkLay4Ndet;
    ladderLength = fgkLay4LadderLength;
    digitCableA  = fDigitCableLay4A;
    digitCableB  = fDigitCableLay4B;
    underSegDH   = fLay4LadderUnderSegDH;
    sensorZPos   = fLay4sensorZPos;
  }
  else {
    printf("AliITSv11GeometrySDD::CreateLadder : error=wrong layer\n");
  };
  Double_t ladderBoxDH = fgkLadderHeight+fgkLadderSegBoxDH+underSegDH;
  TGeoVolumeAssembly *virtualLadder = new TGeoVolumeAssembly("ITSsddLadder");
  
  // placing virtual ladder segment following detector ordering convention
  //=======================================================================
  char transName[30];

  // adding segment this way to create cable points in the correct order ...
  for (Int_t iSegment = nDetectors/2-1; iSegment >= 0; iSegment-- ) {

    //TGeoVolumeAssembly *laddSegment = CreateLadderSegment(iLay, iSegment);
    TGeoVolume *laddSegment = CreateLadderSegment(iLay, iSegment);
    sprintf(transName, "ITSsddLay%iLaddSeg%i", iLay, iSegment);
    Double_t segmentPos = fgkSegmentLength*(nDetectors/2-1-iSegment) 
                          + fgkSegmentLength/2;
    TGeoTranslation *segTr = new TGeoTranslation(transName, 0,
						 underSegDH/2,segmentPos);
    ////
    virtualLadder->AddNode(laddSegment, iSegment, segTr);
  };
  for (Int_t iSegment = nDetectors/2; iSegment < nDetectors; iSegment++ ) {

    TGeoVolume *laddSegment = CreateLadderSegment(iLay, iSegment);
    //TGeoVolumeAssembly *laddSegment = CreateLadderSegment(iLay, iSegment);
    sprintf(transName, "ITSsddLay%iLaddSeg%i", iLay, iSegment);
    Double_t segmentPos = fgkSegmentLength*(nDetectors/2-1-iSegment) 
                          + fgkSegmentLength/2;
    TGeoTranslation *segTr = new TGeoTranslation(transName, 0,
						 underSegDH/2,segmentPos);
    ////
    virtualLadder->AddNode(laddSegment, iSegment, segTr);
  };

  // putting virtual volume corresponding to the end of ladder
  //=======================================================================
  TGeoVolumeAssembly *endLadder = CreateEndLadder( iLay );
  Double_t endLength = (ladderLength - nDetectors*fgkSegmentLength)/2.;
  TGeoTranslation *endTrZPos = new TGeoTranslation("ITSsddEndTrZPos",0,0,
				   fgkSegmentLength*(nDetectors/2)+endLength/2.);
  // Euler rotation : about Z, then new X, then new Z
  TGeoRotation *endZNegRot = new TGeoRotation("",90, 180, -90);
  TGeoCombiTrans *endTrZNeg = new TGeoCombiTrans(0,0,
		  -fgkSegmentLength*(nDetectors/2)-endLength/2.,endZNegRot);
  virtualLadder->AddNode(endLadder, 1, endTrZPos);
  virtualLadder->AddNode(endLadder, 2, endTrZNeg);

  // creating and inserting cable segments
  //   (check points are placed while creating segments)
  //=======================================================================
  if (fAddCables)
  for (Int_t iSegment = 0; iSegment < nDetectors; iSegment++ ) {
    
    digitCableA[iSegment].SetInitialNode(virtualLadder);
    digitCableB[iSegment].SetInitialNode(virtualLadder);
    
    for (Int_t iPt=1; iPt<digitCableA[iSegment].GetNCheckPoints(); iPt++ ) {
      Double_t rotation = 0;
      if (iPt>1) {
	rotation = 90-fgkHybridAngle; 
	digitCableA[iSegment].CreateAndInsertCableSegment(iPt, rotation);
      } else
	digitCableA[iSegment].CreateAndInsertCableSegment(iPt);

    };
    
    for (Int_t iPt=1; iPt<digitCableB[iSegment].GetNCheckPoints(); iPt++ ) {
      Double_t rotation = 0;
      if (iPt>1) {
	rotation = fgkHybridAngle-90; 
	digitCableB[iSegment].CreateAndInsertCableSegment(iPt, rotation);
      } else
	digitCableB[iSegment].CreateAndInsertCableSegment(iPt);
    };
  };
  
  // HV cable
  //=======================================================================
  if (fAddHVcables) {
    TGeoMedium *polyhamideSDD = GetMedium("SDDKAPTON (POLYCH2)$");  //ITSsddKAPTON_POLYCH2
    TGeoMedium *alSDD         = GetMedium("AL$");   //ITSal

  AliITSv11GeomCableFlat cableHV[fgkLay4Ndet];  // temp !!!
  char cableHVname[30];
  for (Int_t iSegment = 0; iSegment<nDetectors; iSegment++) {
    sprintf(cableHVname,"ITSsddHVcable%i", iSegment);
    cableHV[iSegment].SetName(cableHVname);
    cableHV[iSegment].SetThickness(fgkLongHVcablePolyThick+fgkLongHVcableAlThick);
    cableHV[iSegment].SetWidth(fgkTransitHVtailWidth);
    cableHV[iSegment].SetNLayers(2);
    cableHV[iSegment].SetLayer(0, fgkLongHVcablePolyThick, polyhamideSDD,
			       fColorPolyhamide);
    cableHV[iSegment].SetLayer(1, fgkLongHVcableAlThick, alSDD, fColorAl);
    cableHV[iSegment].SetInitialNode(virtualLadder);
  };
  Double_t x1[3], x2[3], x3[3],
	   vY[3] = {0,1,0}, vZ[3] = {0,0,1}, vYZ[3]={0,1,1};

  x1[0] = -fgkTransitHVtailXpos;
  x2[0] = -fgkTransitHVtailXpos;
  x3[0] = -fgkTransitHVtailXpos;
  for (Int_t iSegment  = nDetectors/2-1; iSegment >= 0; iSegment-- ) {
    Double_t cableSeparation = TMath::Abs(iSegment - (nDetectors/2-1))
                               *fgkLongHVcableSeparation;
    // adjust where HV long cable starts in Y
    // useful if you want to let some space for alignment
    x1[1] = - ladderBoxDH/2 + 2*fgkmm;
    x2[1] = - ladderBoxDH/2 + underSegDH - cableSeparation
            - (fgkLongHVcablePolyThick+fgkLongHVcableAlThick)/2;
    x3[1] = x2[1];
    x1[2] = sensorZPos[iSegment]+fgkTransitHVtailLength-5*fgkmm;
    x2[2] =  x1[2]+5*fgkmm;
    x3[2] = ladderLength/2-endLength;
    cableHV[iSegment].AddCheckPoint( virtualLadder, 0, x1, vY );
    cableHV[iSegment].AddCheckPoint( virtualLadder, 1, x2, vZ ); // vYZ
    cableHV[iSegment].AddCheckPoint( virtualLadder, 2, x3, vZ );

    //cableHV[iSegment].CreateAndInsertCableSegment(1,0);
    cableHV[iSegment].CreateAndInsertCableCylSegment(1, -45+180);
    //cableHV[iSegment].CreateAndInsertCableSegment(2,0);
    cableHV[iSegment].CreateAndInsertBoxCableSegment(2,0);
  };

  vYZ[2] = -1;
  x1[0] = fgkTransitHVtailXpos;
  x2[0] = fgkTransitHVtailXpos;
  x3[0] = fgkTransitHVtailXpos;

  for (Int_t iSegment = nDetectors/2; iSegment < nDetectors; iSegment++ ) { 
    Double_t cableSeparation = TMath::Abs(iSegment - (nDetectors/2-1))
                               *fgkLongHVcableSeparation;
    x1[1] = - ladderBoxDH/2 + 2*fgkmm; // adjust where HV long cable starts in Y
    x2[1] = - ladderBoxDH/2 + underSegDH - cableSeparation
            - (fgkLongHVcablePolyThick+fgkLongHVcableAlThick)/2;
    x3[1] = x2[1];
    x1[2] = sensorZPos[iSegment]-fgkTransitHVtailLength+5*fgkmm;
    x2[2] =  x1[2]-5*fgkmm;
    x3[2] = -ladderLength/2+endLength;
    cableHV[iSegment].AddCheckPoint( virtualLadder, 0, x1, vY );
    cableHV[iSegment].AddCheckPoint( virtualLadder, 1, x2, vZ ); // vYZ
    cableHV[iSegment].AddCheckPoint( virtualLadder, 2, x3, vZ );

    cableHV[iSegment].CreateAndInsertCableCylSegment(1, -45);
    cableHV[iSegment].CreateAndInsertBoxCableSegment(2,0);
  };
  };

  //**********************************
  if(GetDebug(1)) virtualLadder->CheckOverlaps(0.01);
  return virtualLadder;
}


//________________________________________________________________________
TGeoArb8 *AliITSv11GeometrySDD::CreateLadderSide(const char *name,
                         Double_t dz, Double_t angle, Double_t xSign,
                         Double_t L, Double_t H, Double_t l) {
    // Create one half of the V shape corner of CF ladder
  
    TGeoArb8 *cfLaddSide = new TGeoArb8(dz);
    cfLaddSide->SetName(name);

    // Points must be in clockwise order
    cfLaddSide->SetVertex(0, 0,  0);
    cfLaddSide->SetVertex(2, xSign*(L*TMath::Sin(angle)-l*TMath::Cos(angle)),
			  -L*TMath::Cos(angle)-l*TMath::Sin(angle));
    cfLaddSide->SetVertex(4, 0,  0);
    cfLaddSide->SetVertex(6, xSign*(L*TMath::Sin(angle)-l*TMath::Cos(angle)),
			  -L*TMath::Cos(angle)-l*TMath::Sin(angle));
    if (xSign < 0) {
     cfLaddSide->SetVertex(1, 0, -H);
     cfLaddSide->SetVertex(3, xSign*L*TMath::Sin(angle), -L*TMath::Cos(angle));
     cfLaddSide->SetVertex(5, 0, -H);
     cfLaddSide->SetVertex(7, xSign*L*TMath::Sin(angle), -L*TMath::Cos(angle));
    } else {
     cfLaddSide->SetVertex(1, xSign*L*TMath::Sin(angle), -L*TMath::Cos(angle));
     cfLaddSide->SetVertex(3, 0, -H);
     cfLaddSide->SetVertex(5, xSign*L*TMath::Sin(angle), -L*TMath::Cos(angle));
     cfLaddSide->SetVertex(7, 0, -H);
    }
    return cfLaddSide;
}


//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateHybrid(Int_t iLRSide) {
  //
  // return a box containing the front-end hybrid
  //

  Double_t roundHoleX       = -fgkHybridWidth/2+fgkHybRndHoleX;

  Double_t screenTotalThick = fgkHybGlueScrnThick+fgkHybUpThick+fgkHybAlThick;
  Double_t lowFLTotalThick  = fgkHybGlueLowThick+fgkHybUpThick+fgkHybAlThick;
//   Double_t upFLTotalThick   = fgkHybGlueUpThick +fgkHybUpThick+fgkHybAlThick;
  Double_t chipsCCTotThick  = fgkHybUnderNiThick+fgkHybGlueAgThick
                              +fgkHybChipThick+2*(fgkHybUpCCThick+fgkHybAlCCThick);
  Double_t ccUpLayerTotThick = fgkHybUpCCThick+fgkHybAlCCThick+fgkHybUpCCThick;
//   Double_t volumeThick = (fgkHybridThBridgeThick+screenTotalThick+lowFLTotalThick
// 			  + upFLTotalThick + ccUpLayerTotThick);
  Double_t volumeThick = (fgkHybridThBridgeThick+screenTotalThick+lowFLTotalThick
			  +fgkHybSMDheight);
  Double_t lowLayerYmin     = -volumeThick/2+fgkHybridThBridgeThick
                              +screenTotalThick;
  Double_t flUpThick        = fgkHybGlueUpThick+fgkHybUpThick;

  //**************************************************** media :
  TGeoMedium *airSDD                  = GetMedium("SDD AIR$");
  TGeoMedium *carbonFiberLadderStruct = GetMedium("SDDKAPTON (POLYCH2)$"); //ITSsddCarbonM55J
  TGeoMedium *alSDD                   = GetMedium("AL$"); //ITSal
  TGeoMedium *alSDD80p100             = GetMedium("AL$");                 // to code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoMedium *alSDD50p100             = GetMedium("AL$");                 // to code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoMedium *polyhamideSDD           = GetMedium("SDDKAPTON (POLYCH2)$"); //ITSsddKAPTON_POLYCH2
  TGeoMedium *niSDD                   = GetMedium("COPPER$");                // to code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoMedium *glueAG                  = GetMedium("SDDKAPTON (POLYCH2)$");  // to code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoMedium *siliconSDD              = GetMedium("SDD SI CHIP$"); //ITSsddSiChip
  TGeoMedium *medSMD                  = GetMedium("SDD X7R capacitors$");      //  SDDX7Rcapacitors   TO CHECK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoMedium *medSMDweld              = GetMedium("SDD X7R capacitors$");      //  SDDX7Rcapacitors  TO CHECK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  //**************************************************** main volume :
  TGeoBBox *hybridBox = new TGeoBBox("",fgkHybridWidth/2, volumeThick/2,
				     (fgkHybridLength)/2);
  TGeoVolume *hybrid = new TGeoVolume("ITSsddHybridVol", hybridBox,
				      airSDD);
 
  TGeoBBox *sThermalBridge = new TGeoBBox( "", fgkHybridWidth/2,
					   fgkHybridThBridgeThick/2,
					   fgkHybridLength/2);

  //**************************************************** Thermal bridge :
  TGeoVolume *vThermalBridge = new TGeoVolume("ITSsddHybridThBridge",
					      sThermalBridge,
					      carbonFiberLadderStruct);
  vThermalBridge->SetLineColor(fColorCarbonFiber);
  TGeoTranslation *thBridgeTr = new TGeoTranslation(0, -volumeThick/2
				    +fgkHybridThBridgeThick/2, 0);
  hybrid->AddNode(vThermalBridge, 1, thBridgeTr);

  //**************************************************** Screen layer :
  TGeoBBox *sAlScreenLayer = new TGeoBBox("sAlScreenLayer", fgkHybridWidth/2,
					  fgkHybAlThick/2, fgkHybridLength/2);
  //here the upedex and glue layers are both assumed to be polyimide
  TGeoBBox *sUpGlueScreenLayer = new TGeoBBox("sUpGlueScreenLayer",
					      fgkHybridWidth/2,
				     (fgkHybUpThick+fgkHybGlueScrnThick)/2,
					      fgkHybridLength/2);
  TGeoTube *sRoundHole = new TGeoTube("sRoundHole", 0, fgkHybRndHoleRad,
				      (screenTotalThick+lowFLTotalThick)/2);

  TGeoTranslation *upGlueScreenTr = new TGeoTranslation("upGlueScreenTr",0,
   -volumeThick/2+fgkHybridThBridgeThick+(fgkHybUpThick+fgkHybGlueScrnThick)/2,0);

  TGeoTranslation *alScreenTr = new TGeoTranslation("AlScreenTr", 0,
   -volumeThick/2+fgkHybridThBridgeThick+fgkHybUpThick+fgkHybGlueScrnThick
   +fgkHybAlThick/2, 0);

  TGeoTranslation hybHolePos1Tr(roundHoleX,
   -volumeThick/2+fgkHybridThBridgeThick+(screenTotalThick+lowFLTotalThick)/2,
				-fgkHybridLength/2+fgkHybRndHoleZ);
  TGeoTranslation hybHolePos2Tr(roundHoleX, 
   -volumeThick/2+fgkHybridThBridgeThick+(screenTotalThick+lowFLTotalThick)/2,
				fgkHybridLength/2-fgkHybRndHoleZ);

  TGeoRotation *rotHole = new TGeoRotation("", 0, 90, 0);
  TGeoCombiTrans *hybHolePos1 = new TGeoCombiTrans(hybHolePos1Tr, *rotHole);
  hybHolePos1->SetName("hybHolePos1");
  TGeoCombiTrans *hybHolePos2 = new TGeoCombiTrans(hybHolePos2Tr, *rotHole);
  hybHolePos2->SetName("hybHolePos2");

  upGlueScreenTr->RegisterYourself();
  alScreenTr->RegisterYourself();
  hybHolePos1->RegisterYourself();
  hybHolePos2->RegisterYourself();
  delete rotHole;

  TGeoCompositeShape *sScreenAl = new TGeoCompositeShape(
		      "sAlScreenLayer:AlScreenTr-(sRoundHole:hybHolePos1"
		      "+sRoundHole:hybHolePos2)");
  TGeoVolume *vScreenAl = new TGeoVolume("vScreenAl",sScreenAl, alSDD);
  vScreenAl->SetLineColor(fColorAl);
  TGeoCompositeShape *sScreenUpGlue = new TGeoCompositeShape(
	"sUpGlueScreenLayer:upGlueScreenTr-(sRoundHole:hybHolePos1"
        "+sRoundHole:hybHolePos2)");
  TGeoVolume *vScreenUpGlue = new TGeoVolume("vScreenUpGlue",
					     sScreenUpGlue,polyhamideSDD);
  vScreenUpGlue->SetLineColor(fColorPolyhamide);

  hybrid->AddNode(vScreenUpGlue, 1, 0);
  hybrid->AddNode(vScreenAl, 1, 0);

  //****************************************************  FL low layer :
  Double_t sideWidth1 = fgkHybFLlowChipZ1 - fgkHybFLlowHoleDZ/2;
  Double_t sideWidth2 = fgkHybridLength - fgkHybFLlowChipZ4 - fgkHybFLlowHoleDZ/2;

  //here the upedex and glue layers are both assumed to be polyimide
  TGeoBBox *sUpGlueBar1 = new TGeoBBox("sUpGlueBar1", fgkHybridWidth/2,
				      (fgkHybGlueLowThick+fgkHybUpThick)/2,
				       sideWidth1/2);
  TGeoBBox *sAlBar1 = new TGeoBBox("sAlBar1", fgkHybridWidth/2,
				  fgkHybAlThick/2, sideWidth1/2);
 
  TGeoTranslation *upGlueBarTr1 = new TGeoTranslation("upGlueBarTr1", 0,
			      lowLayerYmin+(fgkHybGlueLowThick+fgkHybUpThick)/2,
					      -(fgkHybridLength-sideWidth1)/2);
  TGeoTranslation *alBarTr1 = new TGeoTranslation("alBarTr1", 0,
		    lowLayerYmin+fgkHybGlueLowThick+fgkHybUpThick+fgkHybAlThick/2,
					      -(fgkHybridLength-sideWidth1)/2);
  upGlueBarTr1->RegisterYourself();
  alBarTr1->RegisterYourself();

  TGeoCompositeShape *sLowUpGlueBar1 = new TGeoCompositeShape(
		     "sUpGlueBar1:upGlueBarTr1-sRoundHole:hybHolePos1");
  TGeoCompositeShape *sLowAlBar1 = new TGeoCompositeShape(
		      "sAlBar1:alBarTr1-sRoundHole:hybHolePos1");
  TGeoVolume *vLowUpGlueBar1 = new TGeoVolume("vLowUpGlueBar1",
					      sLowUpGlueBar1, polyhamideSDD);
  TGeoVolume *vLowAlBar1 = new TGeoVolume("vLowAlBar1",
					  sLowAlBar1, alSDD);
  vLowUpGlueBar1->SetLineColor(fColorPolyhamide);
  vLowAlBar1->SetLineColor(fColorAl);
  hybrid->AddNode(vLowUpGlueBar1,1,0);
  hybrid->AddNode(vLowAlBar1,1,0);

  //---
  //here the upedex and glue layers are both assumed to be polyimide
  TGeoBBox *sUpGlueBar2 = new TGeoBBox("sUpGlueBar2", fgkHybridWidth/2,
				       (fgkHybGlueLowThick+fgkHybUpThick)/2,
				       sideWidth2/2);
  TGeoBBox *sAlBar2 = new TGeoBBox("sAlBar2", fgkHybridWidth/2,
				   fgkHybAlThick/2, sideWidth2/2);

 TGeoTranslation *upGlueBarTr2 = new TGeoTranslation("upGlueBarTr2", 0,
			      lowLayerYmin+(fgkHybGlueLowThick+fgkHybUpThick)/2,
					       (fgkHybridLength-sideWidth2)/2);
  TGeoTranslation *alBarTr2 = new TGeoTranslation("alBarTr2", 0,
		    lowLayerYmin+fgkHybGlueLowThick+fgkHybUpThick+fgkHybAlThick/2,
					       (fgkHybridLength-sideWidth2)/2);
  upGlueBarTr2->RegisterYourself();
  alBarTr2->RegisterYourself();

  TGeoCompositeShape *sLowUpGlueBar2 = new TGeoCompositeShape(
		     "sUpGlueBar2:upGlueBarTr2-sRoundHole:hybHolePos2");
  TGeoCompositeShape *sLowAlBar2 = new TGeoCompositeShape(
		      "sAlBar2:alBarTr2-sRoundHole:hybHolePos2");
  TGeoVolume *vLowUpGlueBar2 = new TGeoVolume("vLowUpGlueBar2",sLowUpGlueBar2,
					      polyhamideSDD);
  TGeoVolume *vLowAlBar2 = new TGeoVolume("vLowAlBar2",sLowAlBar2,
					  alSDD);
  vLowUpGlueBar2->SetLineColor(fColorPolyhamide);
  vLowAlBar2->SetLineColor(fColorAl);
  hybrid->AddNode(vLowUpGlueBar2, 1, 0);
  hybrid->AddNode(vLowAlBar2, 1, 0);

  if(GetDebug(3)) { // Remove compiler warning.
    sAlScreenLayer->InspectShape();
    sUpGlueScreenLayer->InspectShape();
    sRoundHole->InspectShape();
    sUpGlueBar1->InspectShape();
    sUpGlueBar2->InspectShape();
    sAlBar1->InspectShape();
    sAlBar2->InspectShape();
  };
  //---
  //using class AliITSv11GeomCableFlat to add 2-layer segments ... 
  Double_t piece1width = fgkHybFLlowPasX-fgkHybFLlowHolePasDX/2;
  AliITSv11GeomCableFlat lowFLpiece("lowFLpiece1",piece1width,
				       lowFLTotalThick);
  lowFLpiece.SetNLayers(2);
  lowFLpiece.SetLayer(0, fgkHybGlueLowThick+fgkHybUpThick, polyhamideSDD,
		       fColorPolyhamide);
  lowFLpiece.SetLayer(1, fgkHybAlThick, alSDD80p100, fColorAl);
  // alSDD at 80% : mostly to take into account strips of piece 3

  Double_t x1[3] = { -fgkHybridWidth/2 + piece1width/2,
		     lowLayerYmin + lowFLTotalThick/2,
		     -fgkHybridLength/2 + sideWidth1 };
  Double_t x2[3] ={ x1[0], x1[1], fgkHybridLength/2 - sideWidth2 };
  Double_t vZ[3] = {0,0,1};
  lowFLpiece.AddCheckPoint( hybrid, 0, x2, vZ );
  lowFLpiece.AddCheckPoint( hybrid, 1, x1, vZ );
  lowFLpiece.SetInitialNode(hybrid);
  lowFLpiece.CreateAndInsertBoxCableSegment(1);
  lowFLpiece.ResetPoints();

  Double_t piece2width = fgkHybFLlowAmbX-fgkHybFLlowPasX
                         -fgkHybFLlowHolePasDX/2-fgkHybFLlowHoleAmbDX/2;

  lowFLpiece.SetWidth(piece2width);
  lowFLpiece.SetName("lowFLpiece2");
  x1[0] = piece2width/2+fgkHybFLlowPasX+fgkHybFLlowHolePasDX/2-fgkHybridWidth/2;
  x2[0] = x1[0];
  lowFLpiece.AddCheckPoint( hybrid, 0, x2, vZ );
  lowFLpiece.AddCheckPoint( hybrid, 1, x1, vZ );
  lowFLpiece.CreateAndInsertBoxCableSegment(1);
  lowFLpiece.ResetPoints();

  Double_t piece3width = fgkHybridWidth - fgkHybFLlowAmbX 
                         - fgkHybFLlowHoleAmbDX/2;

  lowFLpiece.SetWidth(piece3width);
  lowFLpiece.SetName("lowFLpiece3");
  x1[0] = fgkHybridWidth/2-piece3width/2;
  x2[0] = x1[0];
  lowFLpiece.AddCheckPoint( hybrid, 0, x2, vZ );
  lowFLpiece.AddCheckPoint( hybrid, 1, x1, vZ );
  lowFLpiece.CreateAndInsertBoxCableSegment(1);

  Double_t zChips[4] = {fgkHybFLlowChipZ1,fgkHybFLlowChipZ2,
			fgkHybFLlowChipZ3,fgkHybFLlowChipZ4};
  Double_t vX[3] = {1,0,0};
  for (Int_t i=0; i<3; i++) {
    char ch[20];
    sprintf(ch, "lowFLpieceA%i", i+4);
    lowFLpiece.SetName(ch);
    lowFLpiece.SetWidth(zChips[i+1]-zChips[i]-fgkHybFLlowHoleDZ);

    lowFLpiece.SetLayer(1, fgkHybAlThick, alSDD, fColorAl);
    x1[0] = -fgkHybridWidth/2 + piece1width;
    x2[0] = x1[0] + fgkHybFLlowHolePasDX;
    Double_t zPiece = (zChips[i+1]+zChips[i])/2 - fgkHybridLength/2;
    x1[2] = zPiece; x2[2] = zPiece; 
    lowFLpiece.AddCheckPoint( hybrid, 0, x2, vX );
    lowFLpiece.AddCheckPoint( hybrid, 1, x1, vX );
    lowFLpiece.CreateAndInsertBoxCableSegment(1,90);
    lowFLpiece.ResetPoints();

    sprintf(ch, "lowFLpieceB%i", i+4);
    lowFLpiece.SetName(ch);
    x1[0] = fgkHybridWidth/2 - piece3width;
    x2[0] = x1[0] - fgkHybFLlowHoleAmbDX;
    lowFLpiece.AddCheckPoint( hybrid, 0, x1, vX );
    lowFLpiece.AddCheckPoint( hybrid, 1, x2, vX );
    lowFLpiece.CreateAndInsertBoxCableSegment(1,90);
  };
  
  //**************************************************** chips+CC:
  AliITSv11GeomCableFlat chip("", fgkHybChipsDZ, chipsCCTotThick);
  chip.SetInitialNode(hybrid);
  chip.SetNLayers(5);
  chip.SetLayer(0, fgkHybUnderNiThick, niSDD, 2);
  chip.SetLayer(1, fgkHybGlueAgThick, glueAG, 4);
  chip.SetLayer(2, fgkHybChipThick, siliconSDD, fColorSilicon);
  chip.SetLayer(3, fgkHybUpCCThick+fgkHybUpCCThick, polyhamideSDD,
		fColorPolyhamide);
  chip.SetLayer(4, fgkHybAlCCThick+fgkHybAlCCThick, alSDD80p100, fColorAl);
  // Here the tho CC (low+up) are merged
  // In fact, the last layer has a smaller surface of Al -> I put 80%
  
  x1[1] = lowLayerYmin + chipsCCTotThick/2;
  x2[1] = x1[1];
  char ch[20];

  for (Int_t i=0; i<4; i++) {
    sprintf(ch, "pascalCC%i", i);
    chip.SetName(ch);
    x1[0] = fgkHybFLlowPasX - fgkHybridWidth/2 - fgkHybPascalDX/2;
    x2[0] = x1[0] + fgkHybPascalDX;
    x1[2] =  zChips[i] - fgkHybridLength/2;
    x2[2] = x1[2];
    chip.AddCheckPoint( hybrid, 0, x1, vX );
    chip.AddCheckPoint( hybrid, 1, x2, vX );
    chip.CreateAndInsertBoxCableSegment(1,-90);
    chip.ResetPoints();

    sprintf(ch, "ambraCC%i", i);
    chip.SetName(ch);
    x1[0] = fgkHybFLlowAmbX - fgkHybridWidth/2 - fgkHybAmbraDX/2;
    x2[0] = x1[0] + fgkHybAmbraDX;
    chip.AddCheckPoint( hybrid, 0, x1, vX );
    chip.AddCheckPoint( hybrid, 1, x2, vX );
    chip.CreateAndInsertBoxCableSegment(1,-90);
    chip.ResetPoints();
  };

  //**************************************************** CC outside chips:
  // I don't think there is a second aluminium layer here ...
  for (Int_t i = 0; i<4; i++) {
    sprintf(ch, "ccLayerA%i", i);

    AliITSv11GeomCableFlat ccLayer1(ch, 6.6*fgkmm, ccUpLayerTotThick);
    ccLayer1.SetInitialNode(hybrid);
    ccLayer1.SetNLayers(2);
    ccLayer1.SetLayer(0, 2*fgkHybUpCCThick, polyhamideSDD, fColorPolyhamide);
    ccLayer1.SetLayer(1, fgkHybAlCCThick, alSDD50p100, fColorAl);
    // Al at ~50%

    x1[0] = -fgkHybridWidth/2;
    x2[0] = fgkHybFLlowPasX - fgkHybridWidth/2 - fgkHybPascalDX/2;
    x1[1] = lowLayerYmin + fgkHybUnderNiThick + fgkHybGlueAgThick
            + fgkHybChipThick + ccUpLayerTotThick/2;
    x2[1] = x1[1];
    x1[2] = zChips[i] - fgkHybridLength/2;
    x2[2] = x1[2];
    ccLayer1.AddCheckPoint( hybrid, 0, x1, vX );
    ccLayer1.AddCheckPoint( hybrid, 1, x2, vX );
    ccLayer1.CreateAndInsertBoxCableSegment(1,-90);

    sprintf(ch, "ccLayerB%i", i);
    AliITSv11GeomCableFlat ccLayer2(ch, fgkHybChipsDZ, ccUpLayerTotThick);
    ccLayer2.SetInitialNode(hybrid);
    ccLayer2.SetNLayers(2);
    ccLayer2.SetLayer(0, 2*fgkHybUpCCThick, polyhamideSDD, fColorPolyhamide);
    ccLayer2.SetLayer(1, fgkHybAlCCThick, alSDD50p100, fColorAl);
     // Al at ~50%

    x1[0] = -fgkHybridWidth/2 + fgkHybFLlowPasX + fgkHybPascalDX/2;
    x2[0] = -fgkHybridWidth/2 + fgkHybFLlowAmbX - fgkHybAmbraDX/2;
    ccLayer2.AddCheckPoint( hybrid, 0, x1, vX );
    ccLayer2.AddCheckPoint( hybrid, 1, x2, vX );
    ccLayer2.CreateAndInsertBoxCableSegment(1,-90);
    ccLayer2.ResetPoints();
    sprintf(ch, "ccLayerC%i", i);
    ccLayer2.SetName(ch);
    x1[0] =  -fgkHybridWidth/2 + fgkHybFLlowAmbX + fgkHybAmbraDX/2;
    x2[0] = fgkHybridWidth/2 - fgkHybFLUpperWidth + 3*fgkmm;
    x1[1] = lowLayerYmin + lowFLTotalThick + flUpThick + fgkHybAlThick
             + ccUpLayerTotThick/2;
    x2[1] = x1[1];

    ccLayer2.AddCheckPoint( hybrid, 0, x1, vX );
    ccLayer2.AddCheckPoint( hybrid, 1, x2, vX );
    ccLayer2.CreateAndInsertBoxCableSegment(1,-90);
  };

  //**************************************************** FL UP:
  // (last Al layer will be a special triangular shape)
  TGeoBBox *sFLupPolyhamide = new TGeoBBox("sFLupPolyhamide",
			      fgkHybFLUpperWidth/2, flUpThick/2,
			      fgkHybFLUpperLength/2);
  TGeoVolume *vFLupPolyhamide = new TGeoVolume("vFLupPolyhamide",
		    		    sFLupPolyhamide, polyhamideSDD);
  vFLupPolyhamide->SetLineColor(fColorPolyhamide);
  TGeoTranslation *trFLupPolyhamide = 
    new TGeoTranslation(fgkHybridWidth/2-fgkHybFLUpperWidth/2,
			lowLayerYmin+lowFLTotalThick+flUpThick/2,0);

  hybrid->AddNode(vFLupPolyhamide, 1, trFLupPolyhamide);

  TGeoArb8 *aluStrip = new TGeoArb8(fgkHybAlThick/2);
  aluStrip->SetVertex( 0,-fgkHybFLUpperAlDZ/2, fgkHybFLUpperWidth);
  aluStrip->SetVertex( 1, fgkHybFLUpperAlDZ/2, fgkHybFLUpperWidth);
  aluStrip->SetVertex( 2, fgkHybFLUpperAlDZ/2, fgkHybFLUpperWidth-fgkHybFLUpperAldx);
  aluStrip->SetVertex( 3,-fgkHybFLUpperAlDZ/2, 0);
  aluStrip->SetVertex( 4,-fgkHybFLUpperAlDZ/2, fgkHybFLUpperWidth);
  aluStrip->SetVertex( 5, fgkHybFLUpperAlDZ/2, fgkHybFLUpperWidth);
  aluStrip->SetVertex( 6, fgkHybFLUpperAlDZ/2, fgkHybFLUpperWidth-fgkHybFLUpperAldx);
  aluStrip->SetVertex( 7,-fgkHybFLUpperAlDZ/2, 0);
  TGeoVolume *vAluStrip = new TGeoVolume("vAluStrip",aluStrip, alSDD50p100);
  // Al at ~50%

  vAluStrip->SetLineColor(fColorAl);
  //TGeoRotation rotAluStrip("rotAluStrip",0, -90, 90);
  TGeoRotation *rotAluStrip = new TGeoRotation("rotAluStrip",0, -90, 90);

  Double_t yRotAluStrip = lowLayerYmin+lowFLTotalThick
                          +flUpThick+fgkHybAlThick/2;
  TGeoCombiTrans *aluStripTr1 = new TGeoCombiTrans(
				    fgkHybridWidth/2,yRotAluStrip,
		    fgkHybridLength/2-fgkHybFLlowChipZ1+1*fgkmm, rotAluStrip);
  TGeoCombiTrans *aluStripTr2 = new TGeoCombiTrans(*aluStripTr1);
  AddTranslationToCombiTrans(aluStripTr2,0,0,
			     fgkHybFLlowChipZ1-fgkHybFLlowChipZ2);
  TGeoCombiTrans *aluStripTr3 = new TGeoCombiTrans(*aluStripTr2);
  AddTranslationToCombiTrans(aluStripTr3,0,0,
			     fgkHybFLlowChipZ2-fgkHybFLlowChipZ3);
  TGeoCombiTrans *aluStripTr4 = new TGeoCombiTrans(*aluStripTr3);
  AddTranslationToCombiTrans(aluStripTr4,0,0,
			     fgkHybFLlowChipZ3-fgkHybFLlowChipZ4);

  hybrid->AddNode(vAluStrip, 1, aluStripTr1); 
  hybrid->AddNode(vAluStrip, 2, aluStripTr2); 
  hybrid->AddNode(vAluStrip, 3, aluStripTr3); 
  hybrid->AddNode(vAluStrip, 4, aluStripTr4); 
  //**************************************************** SMD:
  TGeoBBox *hybSMD = new TGeoBBox("ITSsddSMDshape",
				  fgkHybSMDmiddleL/2+fgkHybSMDendL,
				  fgkHybSMDheight/2,fgkHybSMDendW/2);
  TGeoVolume *vHybSMD = new TGeoVolume("ITSsddSMD",hybSMD,airSDD);

  TGeoBBox *hybSMDmiddle = new TGeoBBox("ITSsddSMDmiddleShape",
			       fgkHybSMDmiddleL/2,fgkHybSMDheight/2,
					fgkHybSMDmiddleW/2);
  TGeoVolume *vHybSMDmiddle = new TGeoVolume("ITSsddSMDmiddle",
					    hybSMDmiddle,medSMD);
  vHybSMDmiddle->SetLineColor(fColorSMD);
  TGeoBBox *hybSMDend = new TGeoBBox("ITSsddSMDendShape",
			    fgkHybSMDendL/2,fgkHybSMDheight/2,fgkHybSMDendW/2);
  TGeoVolume *vHybSMDend = new TGeoVolume("ITSsddSMDend",
					  hybSMDend,medSMDweld);
  vHybSMDend->SetLineColor(fColorSMDweld);
  TGeoTranslation *vHybSMDendTr1 = new TGeoTranslation("",
				       (fgkHybSMDmiddleL+fgkHybSMDendL)/2,0,0);
  TGeoTranslation *vHybSMDendTr2 = new TGeoTranslation("",
				       -(fgkHybSMDmiddleL+fgkHybSMDendL)/2,0,0);
  vHybSMD->AddNode(vHybSMDmiddle,1,0);
  vHybSMD->AddNode(vHybSMDend,1,vHybSMDendTr1);
  vHybSMD->AddNode(vHybSMDend,2,vHybSMDendTr2);
  for (Int_t i=0; i<fgkNHybSMD; i++) {
    TGeoTranslation *vHybSMDtr = new TGeoTranslation("",
				 -fgkHybridWidth/2+fgkHybSMDposX[i],
			         lowLayerYmin+lowFLTotalThick+fgkHybSMDheight/2,
			         -fgkHybridLength/2+fgkHybSMDposZ[i]);
    hybrid->AddNode(vHybSMD, i+1, vHybSMDtr);
  };

  if (iLRSide == 0) {
  };

  if(GetDebug(1)) hybrid->CheckOverlaps(0.01);
  hybrid->SetVisibility(kFALSE);
  return hybrid;
}

//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateLadderSegment(Int_t iLay, Int_t iSeg) {
  //
  // Return a TGeoVolume* containing a segment of a ladder.
  //

  TGeoMedium *phynoxSDD       = GetMedium("INOX$");
  TGeoMedium *coolerMediumSDD = GetMedium("WATER$");
  TGeoMedium *airSDD          = GetMedium("SDD AIR$");

  Double_t tDY = fgkLadderSegBoxDH/2; //space left on top of the ladder 
  Double_t segmentLength = fgkSegmentLength;
  Double_t spaceBetweenCables = 500*fgkmicron;
  
  //*****************************************
  // Set parameters according to (iLay,iSeg):
  //*****************************************
  Int_t nDetectors          = fgkLay3Ndet;
  Double_t coolPipeSuppH    = fgkLay3CoolPipeSuppH;
  Double_t sensorCenterZPos = fLay3sensorZPos[iSeg]-
                              (fgkSegmentLength*fgkLay3Ndet/2. - 
			       fgkSegmentLength/2-(iSeg)*fgkSegmentLength);
 // sensorCenterZPos = z in segment local coord syst.

  AliITSv11GeomCableFlat *digitCableA = fDigitCableLay3A;
  AliITSv11GeomCableFlat *digitCableB = fDigitCableLay3B;

  if (iLay==3) {
  } else if (iLay==4) {
    nDetectors = fgkLay4Ndet;
    coolPipeSuppH = fgkLay4CoolPipeSuppH;
    sensorCenterZPos = fLay4sensorZPos[iSeg]-
                       (fgkSegmentLength*fgkLay4Ndet/2. -
			fgkSegmentLength/2-(iSeg)*fgkSegmentLength);
    digitCableA = fDigitCableLay4A;
    digitCableB = fDigitCableLay4B;	
  } else
    printf("AliITSv11GeometrySDD::CreateLadderSegment Wrong layer index !");

 
  Double_t cableSideSign = -1;
  if (iSeg<nDetectors/2) cableSideSign = 1;
  Double_t spaceForCables = spaceBetweenCables*
           (nDetectors-TMath::Abs(nDetectors-2*iSeg-1)-1)/2
           +0.1*fgkmicron;
  // gives [0-1-2-2-1-0]*spaceBetweenCables
  // or  [0-1-2-3-3-2-1-0]*spaceBetweenCables
  Int_t iUpdateCableMin;
  Int_t iUpdateCableMax;
  if (cableSideSign==-1) {
    iUpdateCableMin = nDetectors/2;
    iUpdateCableMax = iSeg-1;
  } else {
    iUpdateCableMin = iSeg+1;
    iUpdateCableMax = nDetectors/2-1;
  };

  if(GetDebug(1)){
    cout << "Segment ("<< iLay <<',' << iSeg 
	 << ") : sensor z shift in local segment coord.=" 
	 << sensorCenterZPos << endl;
  };

  //****************************
  // The segment volume
  //****************************

  // Use of TGeoVolumeAssembly increases the calculation time of overlaps and very
  // likely slows down the transport of particles through the geometry
 
  //TGeoVolumeAssembly *virtualSeg = new TGeoVolumeAssembly("ITSsddSegment");

  TGeoBBox *segBox = new TGeoBBox("ITSsddSegBox",
				  fgkLadderWidth/2+fgkPinSuppWidth+fgkLadderSegBoxDW,
				  fgkLadderHeight/2+fgkLadderSegBoxDH/2,
				  segmentLength/2);
  
  TGeoVolume *virtualSeg = new TGeoVolume("ITSsddSegment",
					  segBox, airSDD);
  virtualSeg->SetVisibility(kFALSE);

  //******************************
  // Carbon fiber structure :
  //******************************

   virtualSeg->AddNode(fLaddSegCommonVol[0], 1, fLaddSegCommonTr[0]);
  Int_t volumeIndex = 1;
  for (Int_t i = 1; i<fgkNladdSegCommonVol;i++ ) {
    if (fLaddSegCommonVol[i]==fLaddSegCommonVol[i-1])
      volumeIndex++;
    else
      volumeIndex = 1;
    virtualSeg->AddNode(fLaddSegCommonVol[i], volumeIndex,
		 fLaddSegCommonTr[i]);
  };

  //**********************************
  // Pine support of the sensors :
  //**********************************
  TGeoRotation *rotPS1 = new TGeoRotation("",0,-90,90);
  TGeoRotation *rotPS2 = new TGeoRotation("",0,-90,-90);

  // The use of the following constructor type allow to use rotPS1 and rotPS2
  // (and not copy them) therefore we gain some memory
  TGeoCombiTrans *transPS1 = new TGeoCombiTrans( fgkPinDYOnSensor,
				- fgkLadderHeight/2.-tDY
			        + fgkPinSuppHeight/2.,
				sensorCenterZPos+fgkPinDXminOnSensor,rotPS1);

  TGeoCombiTrans *transPS2 = new TGeoCombiTrans( fgkPinDYOnSensor,
				- fgkLadderHeight/2.-tDY
			        + fgkPinSuppHeight/2.,
				sensorCenterZPos+fgkPinDXminOnSensor,rotPS1);
  AddTranslationToCombiTrans(transPS2, 0, 0, fgkPinPinDDXOnSensor);

  TGeoCombiTrans *transPS3 = new TGeoCombiTrans( fgkPinDYOnSensor,
				- fgkLadderHeight/2.-tDY
			        + fgkPinSuppHeight/2.,
				sensorCenterZPos+fgkPinDXminOnSensor,rotPS1);
  AddTranslationToCombiTrans(transPS3, 0, 0, -2*fgkPinDXminOnSensor);

  TGeoCombiTrans *transPS4 = new TGeoCombiTrans( fgkPinDYOnSensor,
				- fgkLadderHeight/2.-tDY
			        + fgkPinSuppHeight/2.,
				sensorCenterZPos+fgkPinDXminOnSensor,rotPS1);
  AddTranslationToCombiTrans(transPS4, 0, 0, -2*fgkPinDXminOnSensor-fgkPinPinDDXOnSensor);

  TGeoCombiTrans *transPS5 = new TGeoCombiTrans( -fgkPinDYOnSensor,
		       	         - fgkLadderHeight/2. - tDY
			         + fgkPinSuppHeight/2.,
			         sensorCenterZPos+fgkPinDXminOnSensor,rotPS2);

  TGeoCombiTrans *transPS6 = new TGeoCombiTrans( -fgkPinDYOnSensor,
		       	         - fgkLadderHeight/2. - tDY
			         + fgkPinSuppHeight/2.,
			         sensorCenterZPos+fgkPinDXminOnSensor,rotPS2);
  AddTranslationToCombiTrans(transPS6, 0, 0, fgkPinPinDDXOnSensor);

  TGeoCombiTrans *transPS7 = new TGeoCombiTrans( -fgkPinDYOnSensor,
		       	         - fgkLadderHeight/2. - tDY
			         + fgkPinSuppHeight/2.,
			         sensorCenterZPos+fgkPinDXminOnSensor,rotPS2);
  AddTranslationToCombiTrans(transPS7, 0, 0, -2*fgkPinDXminOnSensor);

  TGeoCombiTrans *transPS8 = new TGeoCombiTrans( -fgkPinDYOnSensor,
		       	         - fgkLadderHeight/2. - tDY
			         + fgkPinSuppHeight/2.,
			         sensorCenterZPos+fgkPinDXminOnSensor,rotPS2);
  AddTranslationToCombiTrans(transPS8, 0, 0, -2*fgkPinDXminOnSensor-fgkPinPinDDXOnSensor);
  
  virtualSeg->AddNode(fPinSupport, 1, transPS1);
  virtualSeg->AddNode(fPinSupport, 2, transPS2);
  virtualSeg->AddNode(fPinSupport, 3, transPS3);
  virtualSeg->AddNode(fPinSupport, 4, transPS4);
  virtualSeg->AddNode(fPinSupport, 5, transPS5);
  virtualSeg->AddNode(fPinSupport, 6, transPS6);
  virtualSeg->AddNode(fPinSupport, 7, transPS7);
  virtualSeg->AddNode(fPinSupport, 8, transPS8);

  TGeoMedium *pinMed   = GetMedium("SDDKAPTON (POLYCH2)$");  // medium ???
  Double_t fgkPinHeight = 4.5*fgkmm;
  TGeoTube *pineS = new TGeoTube("ITSsddPin",0,fgkPinR,
				fgkPinHeight/2.);
  TGeoVolume *pineV = new TGeoVolume("ITSsddPinVol", pineS, pinMed);

  TGeoCombiTrans *transPS2b = new TGeoCombiTrans( fgkPinDYOnSensor,
				- fgkLadderHeight/2.-tDY
			        + fgkPinHeight/2.,
				sensorCenterZPos+fgkPinDXminOnSensor,rotPS1);
  AddTranslationToCombiTrans(transPS2b, 0, 0, fgkPinPinDDXOnSensor);
  virtualSeg->AddNode(pineV, 1, transPS2b);

  TGeoCombiTrans *transPS6b = new TGeoCombiTrans( -fgkPinDYOnSensor,
		       	         - fgkLadderHeight/2. - tDY
			         + fgkPinHeight/2.,
			         sensorCenterZPos+fgkPinDXminOnSensor,rotPS2);
  AddTranslationToCombiTrans(transPS6b, 0, 0, fgkPinPinDDXOnSensor);
  virtualSeg->AddNode(pineV, 2, transPS6b);

 
  TGeoCombiTrans *transPS4b = new TGeoCombiTrans( fgkPinDYOnSensor,
				- fgkLadderHeight/2.-tDY
			        + fgkPinHeight/2.,
				sensorCenterZPos+fgkPinDXminOnSensor,rotPS1);
  AddTranslationToCombiTrans(transPS4b, 0, 0, -2*fgkPinDXminOnSensor-fgkPinPinDDXOnSensor);
  virtualSeg->AddNode(pineV, 3, transPS4b);

  TGeoCombiTrans *transPS8b = new TGeoCombiTrans( -fgkPinDYOnSensor,
		       	         - fgkLadderHeight/2. - tDY
			         + fgkPinHeight/2.,
			         sensorCenterZPos+fgkPinDXminOnSensor,rotPS2);
  AddTranslationToCombiTrans(transPS8b, 0, 0, -2*fgkPinDXminOnSensor-fgkPinPinDDXOnSensor);
  virtualSeg->AddNode(pineV, 4, transPS8b);


  //******************************
  // Cooling pipe supports :
  //******************************
  Double_t triangleHeight = fgkLadderHeight - fgkLadderBeamRadius;
  Double_t halfTheta = TMath::ATan( 0.5*fgkLadderWidth/triangleHeight );
  Double_t triangleCPaxeDist = fgkCoolPipeSuppAxeDist-fgkCoolPipeSuppWidthExt-
                               fgkCoolPipeSuppWidthIn+fgkLadderBeamRadius;
  
  Double_t coolPipeSuppL = TMath::Tan(halfTheta)*
                           (triangleHeight+triangleCPaxeDist/
			    TMath::Sin(halfTheta)-coolPipeSuppH);
  if (fAddCoolingSyst) {
  TGeoRotation *rotCPS2 = new TGeoRotation("", -halfTheta*TMath::RadToDeg(), -90,  90);
  TGeoRotation *rotCPS1 = new TGeoRotation("",  halfTheta*TMath::RadToDeg(), -90, -90);
  TGeoCombiTrans *transCPS1 = new TGeoCombiTrans(coolPipeSuppL,
				  -fgkLadderHeight/2. - tDY
				  +coolPipeSuppH+fgkLadderBeamRadius,
			          -segmentLength/2., rotCPS1);

  TGeoCombiTrans *transCPS3 = new TGeoCombiTrans(coolPipeSuppL,
				  -fgkLadderHeight/2. - tDY
				  +coolPipeSuppH+fgkLadderBeamRadius,
			          -segmentLength/2., rotCPS1);
  AddTranslationToCombiTrans(transCPS3, 0, 0, segmentLength);
  
  TGeoCombiTrans *transCPS2 = new TGeoCombiTrans(-coolPipeSuppL,
				  -fgkLadderHeight/2.- tDY
				  +coolPipeSuppH+fgkLadderBeamRadius,
			          segmentLength/2., rotCPS2);

  TGeoCombiTrans *transCPS4 = new TGeoCombiTrans(-coolPipeSuppL,
				  -fgkLadderHeight/2.- tDY
				  +coolPipeSuppH+fgkLadderBeamRadius,
			          segmentLength/2., rotCPS2);
  AddTranslationToCombiTrans(transCPS4, 0, 0, -segmentLength);
  
  virtualSeg->AddNode(fCoolPipeSupportL, 1, transCPS1);
  virtualSeg->AddNode(fCoolPipeSupportL, 2, transCPS2);
  virtualSeg->AddNode(fCoolPipeSupportR, 1, transCPS3);
  virtualSeg->AddNode(fCoolPipeSupportR, 2, transCPS4);
  };
  
  //************************
  // Cooling pipes :
  //************************
  TGeoTranslation *pipeTr1 = new TGeoTranslation(coolPipeSuppL,
				 -fgkLadderHeight/2. - tDY +
				 fgkLadderBeamRadius+coolPipeSuppH, 0);
  TGeoTranslation *pipeTr2 = new TGeoTranslation(-coolPipeSuppL,
			         -fgkLadderHeight/2.- tDY +
			          fgkLadderBeamRadius+coolPipeSuppH, 0);

  if (fAddCoolingSyst) {
    TGeoTube *coolingPipeShape = new TGeoTube( fgkCoolPipeInnerDiam/2,
					       fgkCoolPipeOuterDiam/2,
					       segmentLength/2);
    TGeoTube *coolerShape = new TGeoTube( 0, fgkCoolPipeInnerDiam/2,
					  segmentLength/2);
    
    TGeoVolume *coolingPipe = new TGeoVolume("ITSsddCoolingPipe",
					     coolingPipeShape, phynoxSDD );
    coolingPipe->SetLineColor(fColorPhynox);
    TGeoVolume *cooler = new  TGeoVolume("ITSsddCoolingLiquid",coolerShape,
					 coolerMediumSDD );
    
    
    virtualSeg->AddNode(coolingPipe, 1, pipeTr1);
    virtualSeg->AddNode(coolingPipe, 2, pipeTr2);
    if (fCoolingOn) {
      virtualSeg->AddNode(cooler, 1, pipeTr1);
      virtualSeg->AddNode(cooler, 2, pipeTr2);
    };
  };

  //**********************************
  // Bases of hybrid thermal bridges
  //**********************************
  Double_t shiftHyb = 1.05; // shift between thermal Bridge base and thermal bridge
                           // approx !!! not clear on 0752/14-A
  if (fAddCoolingSyst) {
  TGeoRotation rotHybrid1("", 0,   0, -90 - fgkHybridAngle);
  TGeoRotation rotHybrid2("", 0 ,180,  90 - fgkHybridAngle);
  TGeoCombiTrans *baseTr1 = new TGeoCombiTrans(*pipeTr2, rotHybrid1);
  TGeoCombiTrans *baseTr2 = new TGeoCombiTrans(*pipeTr1, rotHybrid2);
  
  virtualSeg->AddNode(fBaseThermalBridge, 1, baseTr1);
  virtualSeg->AddNode(fBaseThermalBridge, 2, baseTr2);
  };

  //*************************
  // the 2 hybrids :
  //*************************
  Double_t hybDy = ((TGeoBBox*)fHybrid->GetShape())->GetDY();
  Double_t distAxeToHybridCenter = fgkBTBaxisAtoBase+hybDy;
  
  Double_t hybrVolX = ( distAxeToHybridCenter*CosD(fgkHybridAngle) 
			 - shiftHyb*SinD(fgkHybridAngle) );
  Double_t hybrVolY = ( distAxeToHybridCenter*SinD(fgkHybridAngle)
			 + shiftHyb*CosD(fgkHybridAngle) );
  if (fAddHybrids) {
    TGeoRotation rotHybrid3("", 0,   0,  90. - fgkHybridAngle);
    TGeoRotation rotHybrid4("", 0 ,180, -90. - fgkHybridAngle);
    TGeoCombiTrans *hybTr1 = new TGeoCombiTrans(*pipeTr2, rotHybrid3);
    TGeoCombiTrans *hybTr2 = new TGeoCombiTrans(*pipeTr1, rotHybrid4);
    AddTranslationToCombiTrans( hybTr1, -hybrVolX, hybrVolY, 0);
    AddTranslationToCombiTrans( hybTr2,  hybrVolX, hybrVolY, 0);
    
    virtualSeg->AddNode(fHybrid, 1, hybTr1);
    virtualSeg->AddNode(fHybrid, 2, hybTr2);
  };

  //***********
  // cables
  //***********
  if (fAddCables) {
  // Starting from this segment
  Double_t hybDz = ((TGeoBBox*)fHybrid->GetShape())->GetDZ();
  Double_t hybDx = ((TGeoBBox*)fHybrid->GetShape())->GetDX();
  Double_t posDigitCableAlongHyb = shiftHyb+ hybDx 
                                   - digitCableA->GetWidth()/2;
  Double_t distAxeToDigitCableCenter = distAxeToHybridCenter+hybDy
                                       - digitCableA->GetThickness()/2;

  Double_t digitCableX = ( coolPipeSuppL
			   + distAxeToDigitCableCenter*CosD(fgkHybridAngle)
			   - posDigitCableAlongHyb*SinD(fgkHybridAngle) );
  Double_t digitCableY = ( - fgkLadderHeight/2.-TMath::Abs(tDY)
			   + fgkLadderBeamRadius+coolPipeSuppH
			   + distAxeToDigitCableCenter*SinD(fgkHybridAngle)
			   + posDigitCableAlongHyb*CosD(fgkHybridAngle) );


  Double_t digitCableCenterA0[3]={ -cableSideSign*digitCableX,
				   digitCableY, cableSideSign*hybDz };
  Double_t digitCableCenterA1[3] = { 
           -cableSideSign*(digitCableX+spaceForCables*CosD(fgkHybridAngle)),
	   digitCableY+spaceForCables*SinD(fgkHybridAngle),
	   cableSideSign*segmentLength/2 };

  Double_t digitCableCenterB0[3]={ cableSideSign*digitCableX,
				   digitCableY,cableSideSign*hybDz};
  Double_t digitCableCenterB1[3]={ 
           cableSideSign*(digitCableX+spaceForCables*CosD(fgkHybridAngle)),
	   digitCableY+spaceForCables*SinD(fgkHybridAngle),
	   cableSideSign*segmentLength/2 };

  Double_t vZ[3] = {0,0,1};
  digitCableA[iSeg].AddCheckPoint( virtualSeg, 0, digitCableCenterA0, vZ);
  digitCableA[iSeg].AddCheckPoint( virtualSeg, 1, digitCableCenterA1, vZ);
  digitCableB[iSeg].AddCheckPoint( virtualSeg, 0, digitCableCenterB0, vZ);
  digitCableB[iSeg].AddCheckPoint( virtualSeg, 1, digitCableCenterB1, vZ);

  // Updating the other cables
  for (Int_t iCable=iUpdateCableMin; iCable<=iUpdateCableMax; iCable++) {

    Int_t iPoint = TMath::Abs(iCable-iSeg)+1;
    Double_t coord[3];
    digitCableA[iCable].GetPoint( 1, coord);
    digitCableA[iCable].AddCheckPoint( virtualSeg, iPoint, coord, vZ);
    digitCableB[iCable].GetPoint( 1, coord);
    digitCableB[iCable].AddCheckPoint( virtualSeg, iPoint, coord, vZ);
  };
  };

  //**********************************
  if(GetDebug(1)) virtualSeg->CheckOverlaps(0.01);
  return virtualSeg;
}


//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreatePinSupport() {
//
// Create a pine support and its pine
// axis of rotation is the cone axis, center in its middle
//
    TGeoMedium *rytonSDD = GetMedium("SDD C AL (M55J)$"); //medium = ryton ?

    TGeoCone *cone = new TGeoCone("ITSsddPinSuppCone",fgkPinSuppHeight/2.,
                                  0,fgkPinSuppRmax,0,fgkPinSuppRmax-
                                  fgkPinSuppHeight*TanD(fgkPinSuppConeAngle) );
    TGeoBBox *tong = new TGeoBBox("ITSsddPinSuppTong",fgkPinSuppRmax,
                                  fgkPinSuppLength/2.,fgkPinSuppThickness/2.);
    TGeoTube *hole = new TGeoTube("ITSsddPinSuppHole",0,fgkPinR,
                                  fgkPinSuppHeight/2.+0.00001);
    // 0.00001 is for seing the actual hole (avoid viewer artefact)

    if(GetDebug(3)){// Remove compiler warning.
        cone->InspectShape();
        tong->InspectShape();
        hole->InspectShape();
    };

    TGeoTranslation *tongTrans = new TGeoTranslation("ITSsddPinSuppTongTr",0,
                   fgkPinSuppLength/2.,-fgkPinSuppHeight/2.+fgkPinSuppThickness/2.);
    tongTrans->RegisterYourself();
    TGeoCompositeShape *pinSupportShape = new TGeoCompositeShape(
               "ITSsddPinSupportShape","(ITSsddPinSuppCone+"
	       "ITSsddPinSuppTong:ITSsddPinSuppTongTr)-ITSsddPinSuppHole");

    TGeoVolume *pinSupport = new TGeoVolume("ITSsddPinSupport", pinSupportShape,
                                            rytonSDD);
    pinSupport->SetLineColor(fColorRyton);

    return pinSupport;
}


//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateCoolPipeSupportL() {
//
// Create half of the cooling pipe support (ALR-0752/3)
//

  Double_t diffX = fgkCoolPipeSuppHeight*TanD(fgkCoolPipeSuppAngle);
  
  TGeoArb8 *side1 = new TGeoArb8(fgkCoolPipeSuppHeight/2.);
  side1->SetName("ITSsddCPSside1");
  side1->SetVertex( 0, 0, -fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 1, 0,  fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 2, fgkCoolPipeSuppMaxLength/2.-diffX,
		       fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 3, fgkCoolPipeSuppMaxLength/2.-diffX,
		       -fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 4, 0, -fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 5, 0,  fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 6, fgkCoolPipeSuppMaxLength/2.,
		       fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 7, fgkCoolPipeSuppMaxLength/2.,
		       -fgkCoolPipeSuppWidthExt/2.);

  TGeoTranslation *side1Tr = new TGeoTranslation("ITSsddCPStr1",0,
				 - fgkCoolPipeSuppAxeDist
				 + fgkCoolPipeSuppWidthExt/2., 0);
  side1Tr->RegisterYourself();
  TGeoTranslation *side2Tr = new TGeoTranslation("ITSsddCPStr2",0,
				 - fgkCoolPipeSuppAxeDist
                                 + fgkCoolPipeSuppWidthExt*3/2.
		      		 + fgkCoolPipeSuppWidthIn,0);
  side2Tr->RegisterYourself();
  
  TGeoBBox *middle = new TGeoBBox("ITSsddCPSmiddle",
			 (fgkCoolPipeSuppMaxLength/2.-fgkCoolPipeSuppSlitL)/2.,
			 fgkCoolPipeSuppWidthIn/2., fgkCoolPipeSuppHeight/2.);
  TGeoTranslation *middleTr = 
    new TGeoTranslation("ITSsddCPStr3",
			(fgkCoolPipeSuppMaxLength/2.-fgkCoolPipeSuppSlitL)/2.,
			-fgkCoolPipeSuppAxeDist+fgkCoolPipeSuppWidthExt
			+fgkCoolPipeSuppWidthIn/2., 0);
  middleTr->RegisterYourself();
  
  TGeoBBox *axeBox = new TGeoBBox("ITSsddCPSaxeBox",
				  fgkCoolPipeSuppTongW/4.,
				  (fgkCoolPipeSuppFulWidth
				   - 2*fgkCoolPipeSuppWidthExt
				   - fgkCoolPipeSuppWidthIn)/2,
				  fgkCoolPipeSuppHeight/2.);
  
  TGeoTranslation *axeBoxTr = new TGeoTranslation("ITSsddCPSAxBoxTr",
				  fgkCoolPipeSuppTongW/4.,
			       	  - fgkCoolPipeSuppAxeDist
				  + fgkCoolPipeSuppFulWidth
				  - axeBox->GetDY(), 0);
  axeBoxTr->RegisterYourself();

  TGeoTube *axe = new TGeoTube("ITSsddCPSaxe",0,fgkCoolPipeSuppHoleDiam/2.,
			       fgkCoolPipeSuppTongW/4.);

  TGeoRotation *axeRot = new TGeoRotation("ITSsddCPSaxeRot",90,90,0);
  TGeoCombiTrans *axeTrans = new TGeoCombiTrans("ITSsddCPSaxeTr",
				 fgkCoolPipeSuppTongW/4.,0,0,axeRot);
  axeTrans->RegisterYourself();
  //delete axeRot; // make the code crash, no idea of why !!!

  if(GetDebug(3)){
    middle->InspectShape();
    axe->InspectShape();
  };

  TGeoMedium *rytonSDD = GetMedium("SDD C AL (M55J)$"); //medium = ryton ?  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  TGeoCompositeShape *coolPipeSuppShape = new TGeoCompositeShape(
				        "ITSsddCoolPipeSuppShapeL",
				        "ITSsddCPSmiddle:ITSsddCPStr3"
				        "+ITSsddCPSside1:ITSsddCPStr1"
				        "+ITSsddCPSside1:ITSsddCPStr2"
				        "+ITSsddCPSaxeBox:ITSsddCPSAxBoxTr"
			                "-ITSsddCPSaxe:ITSsddCPSaxeTr");
  TGeoVolume *coolPipeSupp = new  TGeoVolume("ITSsddCoolPipeSupportL",
					     coolPipeSuppShape, rytonSDD);

  coolPipeSupp->SetLineColor(fColorRyton);

  return coolPipeSupp;
}


//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateCoolPipeSupportR() {
//
//Create half of the cooling pipe support (ALR-0752/3)
//

  Double_t diffX = fgkCoolPipeSuppHeight*TanD(fgkCoolPipeSuppAngle);
  
  TGeoArb8 *side1 = new TGeoArb8(fgkCoolPipeSuppHeight/2.);
  side1->SetName("ITSsddCPSside1R");
  side1->SetVertex( 0, 0, -fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 1, -(fgkCoolPipeSuppMaxLength/2.-diffX),
		       -fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 2, -(fgkCoolPipeSuppMaxLength/2.-diffX),
		       fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 3, 0,  fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 4, 0, -fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 5, -fgkCoolPipeSuppMaxLength/2.,
		       -fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 6, -fgkCoolPipeSuppMaxLength/2.,
		       fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 7, 0,  fgkCoolPipeSuppWidthExt/2.);

  TGeoTranslation *side1Tr = new TGeoTranslation("ITSsddCPStr1R",0,
				 - fgkCoolPipeSuppAxeDist
				 + fgkCoolPipeSuppWidthExt/2., 0);
  side1Tr->RegisterYourself();
  TGeoTranslation *side2Tr = new TGeoTranslation("ITSsddCPStr2R",0,
				 - fgkCoolPipeSuppAxeDist
				 + fgkCoolPipeSuppWidthExt*3/2.
		      		 + fgkCoolPipeSuppWidthIn, 0);
  side2Tr->RegisterYourself();
  
  TGeoBBox *middle = new TGeoBBox("ITSsddCPSmiddleR",
				  (fgkCoolPipeSuppMaxLength/2.
				   - fgkCoolPipeSuppSlitL)/2.,
				  fgkCoolPipeSuppWidthIn/2., 
				  fgkCoolPipeSuppHeight/2.);
  TGeoTranslation *middleTr = 
    new TGeoTranslation("ITSsddCPStr3R",
			-( fgkCoolPipeSuppMaxLength/2.
			   -fgkCoolPipeSuppSlitL)/2.,
			-fgkCoolPipeSuppAxeDist + fgkCoolPipeSuppWidthExt
			+ fgkCoolPipeSuppWidthIn/2.,0);
  middleTr->RegisterYourself();
  
  TGeoBBox *axeBox = new TGeoBBox("ITSsddCPSaxeBoxR",
				  fgkCoolPipeSuppTongW/4.,
				  (fgkCoolPipeSuppFulWidth
				   - 2*fgkCoolPipeSuppWidthExt
				   - fgkCoolPipeSuppWidthIn)/2,
				  fgkCoolPipeSuppHeight/2.);
  
  TGeoTranslation *axeBoxTr = new TGeoTranslation("ITSsddCPSAxBoxTrR",
				  - fgkCoolPipeSuppTongW/4.,
			       	  - fgkCoolPipeSuppAxeDist
				  + fgkCoolPipeSuppFulWidth
				  - axeBox->GetDY(),0);
  axeBoxTr->RegisterYourself();

  TGeoTube *axe = new TGeoTube("ITSsddCPSaxeR",0,fgkCoolPipeSuppHoleDiam/2.,
			       fgkCoolPipeSuppTongW/4.);

  TGeoRotation *axeRot = new TGeoRotation("ITSsddCPSaxeRotR",90,90,0);
  TGeoCombiTrans *axeTrans = new TGeoCombiTrans("ITSsddCPSaxeTrR",
						-fgkCoolPipeSuppTongW/4.,0,0,axeRot);
  axeTrans->RegisterYourself();
  //delete axeRot;

  if(GetDebug(3)){
    middle->InspectShape();
    axe->InspectShape();
  };
  
  TGeoCompositeShape *coolPipeSuppShape = new TGeoCompositeShape(
			              "ITSsddCoolPipeSuppShapeR",
			              "ITSsddCPSmiddleR:ITSsddCPStr3R"
			              "+ITSsddCPSside1R:ITSsddCPStr1R"
			              "+ITSsddCPSside1R:ITSsddCPStr2R"
				      "+ITSsddCPSaxeBoxR:ITSsddCPSAxBoxTrR"
				      "-ITSsddCPSaxeR:ITSsddCPSaxeTrR");
  
  TGeoMedium *rytonSDD = GetMedium("SDD C AL (M55J)$"); //medium = ryton ? To code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoVolume *coolPipeSupp = new TGeoVolume( "ITSsddCoolPipeSupportR",
					     coolPipeSuppShape, rytonSDD);
  coolPipeSupp->SetLineColor(fColorRyton);

  return coolPipeSupp;
}

//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateBaseThermalBridge() {
  //
  // based on ALR 0752/8
  //

  Double_t dy = fgkBTBaxisAtoBase - fgkRadiusBminBTB - fgkBTBthick;

  Double_t base1width = fgkBTBwidth - fgkBTBaxisAtoBottom - fgkRadiusBminBTB
                        - (fgkRadiusAminBTB+fgkBTBthick);
  TGeoBBox *base1 = new TGeoBBox( "ITSsddBTBbase1", base1width/2.,
				  fgkBTBthick/2., fgkBTBlength/2.);
  TGeoTranslation *base1Tr = new TGeoTranslation("ITSsddBTBtr1",
				 fgkBTBaxisAtoBottom-fgkBTBwidth+base1width/2.,
				 -(fgkBTBaxisAtoBase-fgkBTBthick/2.), 0);
  base1Tr->RegisterYourself();

  Double_t base2width = fgkBTBaxisAtoBottom - fgkRadiusAminBTB - fgkBTBthick
                        - fgkRadiusBminBTB;
  TGeoBBox *base2 = new TGeoBBox( "ITSsddBTBbase2", base2width/2.,
				  fgkBTBthick/2., fgkBTBlength/2.);
  TGeoTranslation *base2Tr = new TGeoTranslation("ITSsddBTBtr2",
				 fgkBTBaxisAtoBottom - base2width/2.,
				 -(fgkBTBaxisAtoBase-fgkBTBthick/2.), 0);
  base2Tr->RegisterYourself();

  TGeoBBox *side = new TGeoBBox( "ITSsddBTBside",
				 fgkBTBthick/2., dy/2., fgkBTBlength/2.);
  TGeoTranslation *sideTr1 = new TGeoTranslation("ITSsddBTBsideTr1",
				 -fgkRadiusAminBTB-fgkBTBthick/2., -dy/2., 0);
  TGeoTranslation *sideTr2 = new TGeoTranslation("ITSsddBTBsideTr2",
				 fgkRadiusAminBTB+fgkBTBthick/2., -dy/2., 0);
  sideTr1->RegisterYourself();
  sideTr2->RegisterYourself();

  TGeoBBox *hole = new TGeoBBox( "ITSsddBTBhole", fgkBTBHolewidth/2.,
				 fgkBTBthick/2., fgkBTBHoleLength/2.);
  TGeoTranslation *holeTr1 = new TGeoTranslation("ITSsddBTBholeTr1",
				 - fgkBTBHoleRefX + fgkBTBHolewidth/2.,
				 - (fgkBTBaxisAtoBase-fgkBTBthick/2.),
				 fgkBTBHoleRefY+(fgkBTBHoleLength-fgkBTBlength)/2.);
  TGeoTranslation *holeTr2 = new TGeoTranslation("ITSsddBTBholeTr2",
				 - fgkBTBHoleRefX + fgkBTBHolewidth/2.,
				 - (fgkBTBaxisAtoBase-fgkBTBthick/2.),
				 - fgkBTBHoleRefY-(fgkBTBHoleLength-fgkBTBlength)/2.);
  holeTr1->RegisterYourself();
  holeTr2->RegisterYourself();

  Double_t radiusAmaxBTB = fgkRadiusAminBTB + fgkBTBthick;
  TGeoTubeSeg *mainAxis = new TGeoTubeSeg( "ITSsddBTBmainAxis",
					   fgkRadiusAminBTB, radiusAmaxBTB,
					   fgkBTBlength/2., 0., 180.);
  TGeoTubeSeg *round1 = new TGeoTubeSeg( "ITSsddBTBround1",
			   fgkRadiusBminBTB, fgkRadiusBminBTB+fgkBTBthick,
			   fgkBTBlength/2., 270., 360.);
  TGeoTranslation *roundTr1 = new TGeoTranslation("ITSsddBTBround1Tr",
				  -(fgkRadiusAminBTB+fgkBTBthick+fgkRadiusBminBTB),
				  -dy, 0);
  roundTr1->RegisterYourself();

  TGeoTubeSeg *round2 = new TGeoTubeSeg( "ITSsddBTBround2",
			   fgkRadiusBminBTB, fgkRadiusBminBTB+fgkBTBthick,
			   fgkBTBlength/2., 180., 270.);
  TGeoTranslation *roundTr2 = new TGeoTranslation("ITSsddBTBround2Tr",
				  (fgkRadiusAminBTB+fgkBTBthick+fgkRadiusBminBTB),
				  -dy, 0);
  roundTr2->RegisterYourself();

  TGeoCompositeShape *sBaseThermalBridge = new TGeoCompositeShape(
				      "ITSsddBaseThermalBridgeShape",
				      "ITSsddBTBbase1:ITSsddBTBtr1"
				      "+ ITSsddBTBbase2:ITSsddBTBtr2"
				      "+ ITSsddBTBround1:ITSsddBTBround1Tr"
				      "+ ITSsddBTBround2:ITSsddBTBround2Tr"
				      "+ ITSsddBTBside:ITSsddBTBsideTr1"
				      "+ ITSsddBTBside:ITSsddBTBsideTr2"
				      "- ITSsddBTBhole:ITSsddBTBholeTr1"
				      "- ITSsddBTBhole:ITSsddBTBholeTr2"
				      "+ ITSsddBTBmainAxis");

    if(GetDebug(3)){// Remove compiler warning.
        base1->InspectShape();
        base2->InspectShape();
        side->InspectShape();
        hole->InspectShape();
        mainAxis->InspectShape();
        round1->InspectShape();
        round2->InspectShape();
    };

  TGeoMedium *carbonFiberLadderStruct = GetMedium("SDD C AL (M55J)$");
  TGeoVolume *vBaseThermalBridge = new TGeoVolume( "ITSsddBaseThermalBridge",
						   sBaseThermalBridge,
						   carbonFiberLadderStruct);

  vBaseThermalBridge->SetLineColor(fColorCarbonFiber);
  return vBaseThermalBridge;
}


//________________________________________________________________________
TGeoVolumeAssembly* AliITSv11GeometrySDD::CreateEndLadder(Int_t iLay) {
  //
  // Return an assembly containing a end of a CF ladder.
  //

  TGeoMedium *carbonFiberLadderStruct = GetMedium("SDD C AL (M55J)$"); // ITSsddCarbonM55J
  TGeoMedium *stesalite       = GetMedium("G10FR4$");
  TGeoMedium *phynoxSDD       = GetMedium("INOX$");
  TGeoMedium *coolerMediumSDD = GetMedium("WATER$");

  Double_t length        = (fgkLay3LadderLength-fgkLay3Ndet*fgkSegmentLength)/2.;
  Double_t coolPipeSuppH = fgkLay3CoolPipeSuppH;
  Double_t underSegDH    = fLay3LadderUnderSegDH;
  Double_t footDZ    = fgkRubyZladd3 - fgkLay3Ndet*fgkSegmentLength/2 - length/2;
  // footDZ is also where to place the ruby's center in local Z
  Double_t coolPipeEndLen = (fgkCoolPipeLay3Len-fgkSegmentLength*fgkLay3Ndet)/2;

  if (iLay==3) {
  } else if (iLay==4) {
    length         = (fgkLay4LadderLength-fgkLay4Ndet*fgkSegmentLength)/2.;
    coolPipeSuppH  = fgkLay4CoolPipeSuppH;
    underSegDH     = fLay4LadderUnderSegDH;
    footDZ         = fgkRubyZladd4 - fgkLay4Ndet*fgkSegmentLength/2 - length/2;
    coolPipeEndLen = (fgkCoolPipeLay4Len-fgkSegmentLength*fgkLay4Ndet)/2;
  } else {
    printf("error in AliITSv11GeometrySDD::CreateEndLadder: Wrong layer");
    return 0;
  };
    
  Double_t tDY = (- fgkLadderSegBoxDH/2       //space left on top of the ladder
		  + underSegDH/2);          //space under ladder segment
        // here tDY is not the same as for the segment because the end ladder
        // does not have a space under it, inside the general ladder volume.
  Double_t segmentLength   = fgkSegmentLength;
  Double_t topCornerLength = fgkSegmentLength/2.-fgkLay4LaddTopCornerEnd;

  TGeoVolumeAssembly *virtualEnd = new TGeoVolumeAssembly("ITSsddEnd");
  
  //**********************************
  // coding real matter :
  //**********************************
  Double_t triangleHeight = fgkLadderHeight - fgkLadderBeamRadius;
  Double_t halfTheta   = TMath::ATan( 0.5*fgkLadderWidth/triangleHeight );
  Double_t beta        = (TMath::Pi()-2.*halfTheta)/4.;
  Double_t alpha       = TMath::Pi()*3./4. - halfTheta/2.;
  
  //--- The 3 V shape corners of the Carbon Fiber Ladder
  //--- the top V
  TGeoArb8 *cfLaddTop1 = CreateLadderSide("CFladdTopCornerV1shape",
					  topCornerLength/2., halfTheta, -1,
					  fgkLadderLa, fgkLadderHa, fgkLadderl);
  TGeoVolume *cfLaddTopVol1 = new TGeoVolume("ITSsddCFladdTopCornerV1",
				  cfLaddTop1,carbonFiberLadderStruct);
  cfLaddTopVol1->SetLineColor(fColorCarbonFiber);
  TGeoArb8 *cfLaddTop2 = CreateLadderSide( "CFladdTopCornerV2shape",
					   topCornerLength/2., halfTheta, 1,
					   fgkLadderLa, fgkLadderHa, fgkLadderl);
  TGeoVolume *cfLaddTopVol2 = new TGeoVolume("ITSsddCFladdTopCornerV2",
			          cfLaddTop2,carbonFiberLadderStruct);
  cfLaddTopVol2->SetLineColor(fColorCarbonFiber);
  TGeoTranslation *trTop1 = new TGeoTranslation(0, fgkLadderHeight/2+tDY,
						-(length-topCornerLength)/2.);
  virtualEnd->AddNode(cfLaddTopVol1, 1, trTop1);
  virtualEnd->AddNode(cfLaddTopVol2, 1, trTop1);

  //--- the 2 side V
  TGeoArb8 *cfLaddSide1 = CreateLadderSide( "CFladdSideCornerV1shape",
					    length/2., beta, -1,
					    fgkLadderLb, fgkLadderHb, fgkLadderl);
  TGeoVolume *cfLaddSideVol1 = new TGeoVolume("ITSsddCFladdSideCornerV1",
				   cfLaddSide1,carbonFiberLadderStruct);
  cfLaddSideVol1->SetLineColor(fColorCarbonFiber);
  TGeoArb8 *cfLaddSide2 = CreateLadderSide( "CFladdSideCornerV2shape",
					    length/2., beta, 1,
					    fgkLadderLb, fgkLadderHb, fgkLadderl);
  TGeoVolume *cfLaddSideVol2 = new TGeoVolume("ITSsddCFladdSideCornerV2",
				   cfLaddSide2,carbonFiberLadderStruct);
  cfLaddSideVol2->SetLineColor(fColorCarbonFiber);
  Double_t dYTranslation = ( fgkLadderHeight/2. - 0.5*fgkLadderWidth*
			     TMath::Tan(beta) - fgkLadderBeamRadius );
  
  // because center of the triangle doesn't correspond to virtual vol. center
  Double_t distCenterSideDown =  0.5*fgkLadderWidth/TMath::Cos(beta);
  TGeoCombiTrans *ctSideR = CreateCombiTrans("", distCenterSideDown, 0,
					     alpha*TMath::RadToDeg());
  AddTranslationToCombiTrans(ctSideR, 0, -dYTranslation+tDY, 0);
  TGeoCombiTrans *ctSideL = CreateCombiTrans("", distCenterSideDown, 0, 
					     -alpha*TMath::RadToDeg());
  AddTranslationToCombiTrans(ctSideL, 0, -dYTranslation+tDY, 0);
  virtualEnd->AddNode(cfLaddSideVol1, 1, ctSideR);
  virtualEnd->AddNode(cfLaddSideVol2, 1, ctSideR);
  virtualEnd->AddNode(cfLaddSideVol1, 2, ctSideL);
  virtualEnd->AddNode(cfLaddSideVol2, 2, ctSideL);
  
  //--- The beams
  // Beams on the sides
  Double_t beamPhiPrime = TMath::ASin(1./TMath::Sqrt( (1+TMath::Sin(2*beta)*
		  TMath::Sin(2*beta)/(TanD(fgkBeamSidePhi)*TanD(fgkBeamSidePhi))) ));

  //Euler rotation : about Z, then new X, then new Z
  TGeoRotation *beamRot1 = new TGeoRotation("", 90-2.*beta*TMath::RadToDeg(),
			-beamPhiPrime*TMath::RadToDeg(), -90);
  TGeoRotation *beamRot2 = new TGeoRotation("", 90-2.*beta*TMath::RadToDeg(),
			beamPhiPrime*TMath::RadToDeg(), -90);
  TGeoRotation *beamRot3 = new TGeoRotation("", 90+2.*beta*TMath::RadToDeg(),
			beamPhiPrime*TMath::RadToDeg(), -90);
  TGeoRotation *beamRot4 = new TGeoRotation("", 90+2.*beta*TMath::RadToDeg(),
			-beamPhiPrime*TMath::RadToDeg(), -90);
  TGeoCombiTrans *beamTransf1 = new TGeoCombiTrans(0.5*triangleHeight*
						   TMath::Tan(halfTheta),
					      fgkLadderBeamRadius/2. + tDY,
			         -length/2 + segmentLength/8, beamRot1);
  TGeoCombiTrans *beamTransf3 = new TGeoCombiTrans( 0.5*triangleHeight*
						    TMath::Tan(halfTheta),
						 fgkLadderBeamRadius/2.+tDY,
				-length/2 + 3*segmentLength/8, beamRot2);
  TGeoCombiTrans *beamTransf5 = new TGeoCombiTrans(-0.5*triangleHeight*
						   TMath::Tan(halfTheta),
						fgkLadderBeamRadius/2.+tDY,
			         -length/2 + segmentLength/8, beamRot3);
  TGeoCombiTrans *beamTransf7 = new TGeoCombiTrans(-0.5*triangleHeight*
						   TMath::Tan(halfTheta),
					      fgkLadderBeamRadius/2. + tDY,
			         -length/2+3*segmentLength/8, beamRot4);

  virtualEnd->AddNode(fLaddSegCommonVol[6], 1, beamTransf1);
  virtualEnd->AddNode(fLaddSegCommonVol[6], 2, beamTransf3);
  virtualEnd->AddNode(fLaddSegCommonVol[6], 3, beamTransf5);
  virtualEnd->AddNode(fLaddSegCommonVol[6], 4, beamTransf7);

  //--- Beams of the bottom
  TGeoRotation *bottomBeamRot1 = new TGeoRotation("",90, 90, 90);

  /* Not there actually
  TGeoTubeSeg *bottomBeam1 = new TGeoTubeSeg(0, fgkLadderBeamRadius,
				 fgkLadderWidth/2.-fgkLadderLb/3, 0, 180);
  TGeoVolume *bottomBeam1Vol = new TGeoVolume("ITSsddBottomBeam1Vol",
			           bottomBeam1, carbonFiberLadderStruct);
  bottomBeam1Vol->SetLineColor(fColorCarbonFiber);

  TGeoCombiTrans *bottomBeamTransf1 = new TGeoCombiTrans(0,
                            -(fgkLadderHeight/2-fgkLadderBeamRadius)+tDY,
			    -length/2+fgkSegmentLength/2, bottomBeamRot1);
  virtualEnd->AddNode(bottomBeam1Vol, 1, bottomBeamTransf1);
*/
  TGeoTubeSeg *bottomBeam2 = new TGeoTubeSeg(0, fgkLadderBeamRadius,
			         fgkLadderWidth/2.-fgkLadderLb/3, 0, 90);
  TGeoVolume *bottomBeam2Vol = new TGeoVolume("ITSsddBottomBeam2Vol",
				   bottomBeam2, carbonFiberLadderStruct);
  bottomBeam2Vol->SetLineColor(fColorCarbonFiber);
  TGeoCombiTrans *bottomBeamTransf2 = new TGeoCombiTrans(0,
     -(fgkLadderHeight/2-fgkLadderBeamRadius)+tDY,-length/2,bottomBeamRot1);
  virtualEnd->AddNode(bottomBeam2Vol, 1, bottomBeamTransf2);

  //**********************************
  //the cooling pipe supports
  Double_t triangleCPaxeDist = fgkCoolPipeSuppAxeDist-fgkCoolPipeSuppWidthExt-
                               fgkCoolPipeSuppWidthIn+fgkLadderBeamRadius;

  Double_t coolPipeSuppL = TMath::Tan(halfTheta)*
                           (triangleHeight+triangleCPaxeDist/
			    TMath::Sin(halfTheta)-coolPipeSuppH);
  
  if (fAddCoolingSyst) {
    TGeoRotation *rotCPS2 = new TGeoRotation("",-halfTheta*TMath::RadToDeg(),-90, 90);
    TGeoRotation *rotCPS1 = new TGeoRotation("", halfTheta*TMath::RadToDeg(),-90,-90);
    TGeoCombiTrans *transCPS1 = new TGeoCombiTrans(coolPipeSuppL,
						   -fgkLadderHeight/2.+ tDY +
						   coolPipeSuppH+fgkLadderBeamRadius,
						   -length/2., rotCPS1);
    TGeoCombiTrans *transCPS4 = new TGeoCombiTrans(-coolPipeSuppL,
						   -fgkLadderHeight/2.+ tDY +
						   coolPipeSuppH+fgkLadderBeamRadius,
						   -length/2., rotCPS2);
    
    virtualEnd->AddNode(fCoolPipeSupportL, 1, transCPS1);
    virtualEnd->AddNode(fCoolPipeSupportR, 1, transCPS4);
  };

  //**********************************
  //--- The stesalite foot of the ladder

  Double_t footDY = -(fgkLadderHeight/2-fgkLadderBeamRadius)+tDY
                    - fgkLadFootY/2+fgkLadFingerPrintY;

  TGeoTranslation *footTr = new TGeoTranslation("SDDfootTr",0,footDY,footDZ);
  virtualEnd->AddNode(fLadderFoot, 1, footTr);

  //=====================================
  //--- cooling pipe

  if (fAddCoolingSyst) {

    TGeoTranslation *pipeTr1 = new TGeoTranslation(coolPipeSuppL,
				       -fgkLadderHeight/2.+ tDY +
				       coolPipeSuppH + fgkLadderBeamRadius,
				       -length/2.+coolPipeEndLen/2.);
    TGeoTranslation *pipeTr2 = new TGeoTranslation(-coolPipeSuppL,
				   -fgkLadderHeight/2. + tDY +
				    fgkLadderBeamRadius + coolPipeSuppH,
				    -length/2.+coolPipeEndLen/2.);

    TGeoTube *coolingPipeShape = new TGeoTube( fgkCoolPipeInnerDiam/2,
					       fgkCoolPipeOuterDiam/2,
					       coolPipeEndLen/2);
    TGeoTube *coolerShape = new TGeoTube( 0, fgkCoolPipeInnerDiam/2,
					  coolPipeEndLen/2);
    
    TGeoVolume *coolingPipe = new TGeoVolume("ITSsddCoolingPipeEnd",
					     coolingPipeShape, phynoxSDD );
    coolingPipe->SetLineColor(fColorPhynox);
    TGeoVolume *cooler = new  TGeoVolume("ITSsddCoolingEndLiquid",coolerShape,
					 coolerMediumSDD );

    virtualEnd->AddNode(coolingPipe, 1, pipeTr1);
    virtualEnd->AddNode(coolingPipe, 2, pipeTr2);
    if (fCoolingOn) {
      virtualEnd->AddNode(cooler, 1, pipeTr1);
      virtualEnd->AddNode(cooler, 2, pipeTr2);
    };
  };

  //=====================================
  //--- HV cable guide


  TGeoBBox* guideHVbox = new TGeoBBox("guideHVbox",fgkHVguideX1/2,
				      fgkHVguideY1/2,fgkHVguideZ1/2);
  TGeoVolume *guideHV = new TGeoVolume("guideHV",guideHVbox,stesalite);

  TGeoTranslation* guideHVtr = new TGeoTranslation(fgkHVguideDX,
     -(fgkLadderHeight/2-fgkLadderBeamRadius)+tDY-fgkHVguideY1/2,
     footDZ+fgkLadFootZ/2+fgkHVguideZ1/2-(fgkHVguideSuppFullZ-fgkHVguideZ2));
  virtualEnd->AddNode(guideHV, 1, guideHVtr);

  //=====================================
  //--- raccordo
  Double_t raccordFullLen = fgkConnectorCoolTubeL1+fgkConnectorCoolTubeL2+fgkConnectorCoolTubeL3;
  TGeoTranslation *trRaccordo1 = new TGeoTranslation("trRaccordo1",-coolPipeSuppL,
						     -fgkLadderHeight/2.+ tDY +
						     coolPipeSuppH+fgkLadderBeamRadius,
					    -length/2.+coolPipeEndLen+raccordFullLen/2);
  TGeoTranslation *trRaccordo2 = new TGeoTranslation("trRaccordo2", coolPipeSuppL,
						     -fgkLadderHeight/2.+ tDY +
						     coolPipeSuppH+fgkLadderBeamRadius,
					    -length/2.+coolPipeEndLen+raccordFullLen/2);

  virtualEnd->AddNode(fRaccordoL, 1, trRaccordo1);
  virtualEnd->AddNode(fRaccordoL, 2, trRaccordo2);

  if(GetDebug(1)) virtualEnd->CheckOverlaps(0.01);

  return virtualEnd;
}

//________________________________________________________________________
TGeoVolumeAssembly* AliITSv11GeometrySDD::CreateLadderFoot() {

  //--- The stesalite foot of the ladder
  // Are missing :
  // The 2 screw holes on the left part
  // the small holes at each corner of the ruby cage (diam 2mm)
  // the really small level difference of 0.3mm on the bottom


  TGeoMedium *stesalite = GetMedium("G10FR4$");

  TGeoVolumeAssembly *virtualFoot = new TGeoVolumeAssembly("ITSsddFoot");

  Double_t epsilon = 2e-10;
  TGeoBBox *ladFootBox1 = new TGeoBBox("ladFootBox1",fgkLadBox1X/2, fgkLadFootY/2,
				       fgkLadFootZ/2);
  TGeoTranslation *ladFootBox1Tr = new TGeoTranslation("ladFootBox1Tr",
						       fgkLadFootX/2-fgkLadBox1X/2,0,0);
  TGeoBBox *ladFingerPrint = new TGeoBBox("ladFingerPrint",fgkLadFingerPrintX/2,
					  fgkLadFingerPrintY/2+epsilon, fgkLadFootZ/2+epsilon);

  TGeoTranslation *ladFingerPrintTr = new TGeoTranslation("ladFingerPrintTr",
			    fgkLadFootX/2-fgkLadFingerPrintBorder-fgkLadFingerPrintX/2,
			    fgkLadFootY/2-fgkLadFingerPrintY/2+epsilon,
			    0);

  TGeoBBox *rubyCageHole = new TGeoBBox("rubyCageHole",fgkRubyCageHoleX/2,
					fgkRubyCageHoleY/2+epsilon, fgkRubyCageHoleZ/2);

  TGeoTranslation *rubyCageHoleTr = new TGeoTranslation("rubyCageHoleTr",
				 fgkLadFootX/2-(fgkLadFootX/2-fgkRubyDX)+fgkRubyCageAxisShift,
				 fgkLadFootY/2-fgkRubyCageHoleY/2,0);

  double rubyScrewHoleLen = fgkLadFootX/2-fgkRubyDX;
  TGeoTube *rubyScrewHole = new TGeoTube("rubyScrewHole", 0,fgkScrewM4diam/2,
			    rubyScrewHoleLen/2);

  TGeoRotation *rot9090 = new TGeoRotation("",90,90,0);
  TGeoCombiTrans *rubyScrewHoleTr = new TGeoCombiTrans("rubyScrewHoleTr",
				    fgkLadFootX/2-rubyScrewHoleLen/2,
				    -fgkRubyScrewShiftToCenterY, 0, rot9090);

  Double_t rubyHoleLen = fgkLadFootY-fgkRubyCageHoleY;
  TGeoTube *rubyHole = new TGeoTube("rubyHole", 0,fgkRubyHoleDiam/2,
				    rubyHoleLen/2);

  TGeoRotation *rot90 = new TGeoRotation("",0,90,0);
  TGeoCombiTrans *rubyHoleTr = new TGeoCombiTrans("rubyHoleTr", fgkRubyDX,
				   -(fgkLadFootY-rubyHoleLen)/2, 0, rot90);

  ladFootBox1Tr->RegisterYourself();
  ladFingerPrintTr->RegisterYourself();
  rubyCageHoleTr->RegisterYourself();
  rubyScrewHoleTr->RegisterYourself();
  rubyHoleTr->RegisterYourself();

  TGeoCompositeShape *footRightPart = new TGeoCompositeShape(
	      "ladFootBox1:ladFootBox1Tr-(ladFingerPrint:ladFingerPrintTr"
	      "+rubyCageHole:rubyCageHoleTr+rubyScrewHole:rubyScrewHoleTr"
	      "+rubyHole:rubyHoleTr)");
  TGeoVolume *vFootRightPart = new TGeoVolume("vFootRightPart",
					      footRightPart,stesalite);
  vFootRightPart->SetLineColor(fColorStesalite);
 
  virtualFoot->AddNode(vFootRightPart, 1, 0);


  //--- This was the right part of the foot, now let's do the middle
  //--- and the right parts

  Double_t middleX = fgkLadFootX-fgkLadBox1X-fgkLadFingerPrintX-fgkLadFingerPrintBorder;
  TGeoBBox *footMiddle = new TGeoBBox("footMiddle", middleX/2, fgkLadFootMiddleY/2,
				      fgkLadFootZ/2);
  TGeoTranslation *middleXTr = new TGeoTranslation("middleXTr",
				   fgkLadFootX/2-fgkLadBox1X-middleX/2,
				   fgkLadFootY/2-fgkLadFootMiddleY/2, 0);

  TGeoVolume *vFootMiddle = new TGeoVolume("vFootMiddle", footMiddle,stesalite);
  vFootMiddle->SetLineColor(fColorStesalite);
  virtualFoot->AddNode(vFootMiddle, 1, middleXTr);
  
  //--
  TGeoBBox *footLeftLadFinger = new TGeoBBox("footLeftLadFinger", fgkLadFingerPrintX/2,
					     (fgkLadFootY-fgkLadFingerPrintY)/2,
					     fgkLadFootZ/2);
  TGeoTranslation *footLeftLadFingerTr = new TGeoTranslation("footLeftLadFingerTr",
				   -fgkLadFootX/2+fgkLadFingerPrintBorder+fgkLadFingerPrintX/2,
				   -fgkLadFingerPrintY/2, 0);
  TGeoVolume *vFootLeftLadFinger = new TGeoVolume("vFootLeftLadFinger",footLeftLadFinger,
						  stesalite);
  vFootLeftLadFinger->SetLineColor(fColorStesalite);
  virtualFoot->AddNode(vFootLeftLadFinger, 1, footLeftLadFingerTr);

  //--
  TGeoBBox *footLeft = new TGeoBBox("footLeft", fgkLadFingerPrintBorder/2,
				    fgkLadFootY/2,
				    fgkLadFootZ/2);
  TGeoTranslation *footLeftTr = new TGeoTranslation("footLeftTr",
				   -fgkLadFootX/2+fgkLadFingerPrintBorder/2,
				    0, 0);
  TGeoVolume *vFootLeft = new TGeoVolume("vFootLeft",footLeft,stesalite);
  vFootLeft->SetLineColor(fColorStesalite);
  virtualFoot->AddNode(vFootLeft, 1, footLeftTr);

  if(GetDebug(3)){ // Remove compiler warning.
    ladFingerPrint->InspectShape();
    ladFootBox1->InspectShape();
    rubyCageHole->InspectShape();
    rubyScrewHole->InspectShape();
    rubyHole->InspectShape();
  }

  return virtualFoot;
}

//________________________________________________________________________
TGeoVolumeAssembly* AliITSv11GeometrySDD::CreateCarlosCard(Int_t iLay) {
  //
  // return an assembly containing the CARLOS end-ladder board
  // and the heat bridge
  //

  (void) iLay;
  TGeoMedium *glassFiber  = GetMedium("SDD SI CHIP$");// glassFiber   TO CODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoMedium *siliconChip = GetMedium("SDD SI CHIP$");// ITSsddSiChip
  TGeoMedium *plastiChip  = GetMedium("SDDKAPTON (POLYCH2)$"); // ITSsddKAPTON_POLYCH2
  TGeoMedium *copper      = GetMedium("COPPER$"); 
  TGeoMedium *alCu12SDD   = GetMedium("INOX$"); // ITSsddAlCu12,  to code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoMedium *stainless   = GetMedium("INOX$"); // for screws, what is the material ???????????

  //=========================================
  // cooling support of the Carlos card (HeatBridge):
  TGeoVolumeAssembly *assemblySupCarlos = new TGeoVolumeAssembly("assemblySupCarlos");

  TGeoBBox *supCarlosBoard1 = new TGeoBBox("",fgkCarlosSuppX1/2,fgkCarlosSuppY1/2,
					   fgkCarlosSuppZ/2);
  TGeoBBox *supCarlosBoard2 = new TGeoBBox("",fgkCarlosSuppX2/2,fgkCarlosSuppY2/2,
					   fgkCarlosSuppZ/2);
  TGeoVolume *vSupCarlosBoard1 = new TGeoVolume("vSupCarlosBoard1",
						supCarlosBoard1, alCu12SDD);
  TGeoVolume *vSupCarlosBoard2 = new TGeoVolume("vSupCarlosBoard2",
						supCarlosBoard2, alCu12SDD);
  vSupCarlosBoard1->SetLineColor(4);
  vSupCarlosBoard2->SetLineColor(4);


  Double_t shiftGlob = -fgkCarlosSuppZ/2+fgkCarlosSuppTopLen;
  // shift of the main planes in the direction of their width 
  // the center is fixed at the center of the 2 small fixing arms on each sides.
  //shiftGlob=0.5;

  shiftGlob+= 0.5*fgkCarlosSuppY3/cos((90-fgkCarlosSuppAngle)*TMath::DegToRad());
  shiftGlob-= 0.5*fgkCarlosSuppY2*tan((90-fgkCarlosSuppAngle)*TMath::DegToRad());
  Double_t shiftGlobY = shiftGlob*sin(fgkCarlosSuppAngle*TMath::DegToRad());
  Double_t shiftGlobZ = shiftGlob*cos(fgkCarlosSuppAngle*TMath::DegToRad());

  TGeoTranslation *carlosSupTr1 = new TGeoTranslation( -fgkCarlosSuppX2/2,
				 (-fgkCarlosSuppY1+fgkCarlosSuppY2)/2+shiftGlobY,
						       +shiftGlobZ);

  TGeoTranslation *carlosSupTr2 = new TGeoTranslation( fgkCarlosSuppX1/2,
						     shiftGlobY,
						     shiftGlobZ);

  assemblySupCarlos->AddNode(vSupCarlosBoard1, 0, carlosSupTr1);
  assemblySupCarlos->AddNode(vSupCarlosBoard2, 0, carlosSupTr2);

  //=========================================
  // fixing arm of the cooling support :
  TGeoBBox *supCarlosBoard3 = new TGeoBBox("",fgkCarlosSuppX3/2,fgkCarlosSuppY3/2,
					   fgkCarlosSuppZ3/2);
  TGeoVolume *vSupCarlosBoard3 = new TGeoVolume("vSupCarlosBoard3",
						supCarlosBoard3, alCu12SDD);
  vSupCarlosBoard3->SetLineColor(4);

  // screw inside :
  TGeoTube *littleScrew = new TGeoTube("littleScrew", 0, fgkLittleScrewR,
				       fgkCarlosSuppY3/2);
  TGeoVolume *vLittleScrew = new TGeoVolume("vLittleScrew",
					    littleScrew, stainless);
  TGeoRotation *rotScrew = new TGeoRotation("",0,90,0);
  TGeoCombiTrans *cbScrew1 = new TGeoCombiTrans(0, 0, fgkCarlosSuppZ3/2 -
						fgkLittleScrewHeadR-0.07, rotScrew);
  TGeoCombiTrans *cbScrew2 = new TGeoCombiTrans(0, 0, -fgkCarlosSuppZ3/2 +
						fgkLittleScrewHeadR+0.07, rotScrew);
  vSupCarlosBoard3->AddNode(vLittleScrew,1, cbScrew1);
  vSupCarlosBoard3->AddNode(vLittleScrew,2, cbScrew2);

  TGeoRotation *carlosSupRot = new TGeoRotation("carlosSuppInvertAngle",
						 0, fgkCarlosSuppAngle, 0);
  TGeoCombiTrans *carlosSupTr3 = new TGeoCombiTrans((fgkCarlosSuppX1+
			         fgkCarlosSuppX2+fgkCarlosSuppX3)/2,0,0, carlosSupRot);
  TGeoCombiTrans *carlosSupTr4 = new TGeoCombiTrans(-(fgkCarlosSuppX1+
				 fgkCarlosSuppX2+fgkCarlosSuppX3)/2,0,0, carlosSupRot);
  assemblySupCarlos->AddNode(vSupCarlosBoard3, 0, carlosSupTr3);
  assemblySupCarlos->AddNode(vSupCarlosBoard3, 1, carlosSupTr4);


  //=========================================
  // screws fixing the board on the U tube
  Double_t aaa = fgkCarlosSuppY3; // ???
  //Double_t aaa = fgkCarlosSuppY3/2 + fgkLittleScrewHeadH/2;
  Double_t bbb = fgkCarlosSuppZ3/2 - fgkLittleScrewHeadR;
  Double_t screw1y = ( aaa*cos(TMath::DegToRad()*fgkCarlosSuppAngle) - 
		       bbb*sin(TMath::DegToRad()*fgkCarlosSuppAngle) );
  Double_t screw1z = ( aaa*sin(TMath::DegToRad()*fgkCarlosSuppAngle) + 
		       bbb*cos(TMath::DegToRad()*fgkCarlosSuppAngle) )-0.07;

  TGeoRotation *CarlosSuppRot = (TGeoRotation *)fCommonTr[0];	

  TGeoCombiTrans* lScrewTr1 = new TGeoCombiTrans((fgkCarlosSuppX1+
			      fgkCarlosSuppX2+fgkCarlosSuppX3)/2,
			      screw1y,screw1z, CarlosSuppRot);

  TGeoCombiTrans* lScrewTr2 = new TGeoCombiTrans((fgkCarlosSuppX1+
			      fgkCarlosSuppX2+fgkCarlosSuppX3)/2,
		       	      screw1z,screw1y, CarlosSuppRot);

  TGeoCombiTrans *lScrewTr3 = new TGeoCombiTrans(-(fgkCarlosSuppX1+
			      fgkCarlosSuppX2+fgkCarlosSuppX3)/2,
			      screw1y,screw1z, CarlosSuppRot);

  TGeoCombiTrans *lScrewTr4 = new TGeoCombiTrans(-(fgkCarlosSuppX1+
			      fgkCarlosSuppX2+fgkCarlosSuppX3)/2,
			      screw1z,screw1y, CarlosSuppRot);

  assemblySupCarlos->AddNode(fCommonVol[0], 1, lScrewTr1);
  assemblySupCarlos->AddNode(fCommonVol[0], 2, lScrewTr2);
  assemblySupCarlos->AddNode(fCommonVol[0], 3, lScrewTr3);
  assemblySupCarlos->AddNode(fCommonVol[0], 4, lScrewTr4);

  //=========================================
  // board
  Double_t p1[3], p2[3], vX[3] = {1,0,0};
  AliITSv11GeomCableFlat card1("cardCarlos1", fgkCarlosCardZ1, fgkCarlosCardY1); // name, width, thickness
  card1.SetNLayers(2);
  card1.SetLayer(0, fgkCarlosCardCuY, copper,     kOrange); // index, thickness, material, color
  card1.SetLayer(1, fgkCarlosCardY1-fgkCarlosCardCuY,   glassFiber, 30);
  card1.SetInitialNode( (TGeoVolume *) assemblySupCarlos);
  p1[0] = -fgkCarlosCardX1/2;
  p1[1] = shiftGlobY - fgkCarlosCard2HeatBridge;
  p1[2] = fgkCarlosCardShift;
  p2[0] = fgkCarlosCardX1/2;
  p2[1] = shiftGlobY - fgkCarlosCard2HeatBridge;
  p2[2] = fgkCarlosCardShift;
  card1.AddCheckPoint( (TGeoVolume *) assemblySupCarlos, 0, p1, vX);
  card1.AddCheckPoint( (TGeoVolume *) assemblySupCarlos, 1, p2, vX);
  card1.CreateAndInsertBoxCableSegment(1,90);

  AliITSv11GeomCableFlat card2("cardCarlos2", fgkCarlosCardZ2, fgkCarlosCardY1); // name, width, thickness
  card2.SetNLayers(2);
  card2.SetLayer(0, fgkCarlosCardCuY, copper,     kOrange); // index, thickness, material, color
  card2.SetLayer(1, fgkCarlosCardY1-fgkCarlosCardCuY,   glassFiber, 30);
  card2.SetInitialNode( (TGeoVolume *) assemblySupCarlos);

  p1[0] = -fgkCarlosCardX1/2;
  p1[1] = shiftGlobY - fgkCarlosCard2HeatBridge;
  p1[2] = fgkCarlosCardShift + fgkCarlosCardZ1/2 + fgkCarlosCardZ2/2;

  p2[0] = -fgkCarlosCardX1/2 + fgkCarlosCardX2;
  p2[1] = shiftGlobY - fgkCarlosCard2HeatBridge;
  p2[2] = fgkCarlosCardShift + fgkCarlosCardZ1/2 + fgkCarlosCardZ2/2;
  card2.AddCheckPoint( (TGeoVolume *) assemblySupCarlos, 0, p1, vX);
  card2.AddCheckPoint( (TGeoVolume *) assemblySupCarlos, 1, p2, vX);
  card2.CreateAndInsertBoxCableSegment(1,90);

  //=========================================
  // some chips on the board 

  AliITSv11GeomCableFlat u1("carlosCardU1", fgkCarlosU1Z, fgkCarlosU1Y); // name, width, thickness
  u1.SetNLayers(2);
  u1.SetLayer(0, fgkCarlosCardChipSiThick, siliconChip, kGreen); // index, thickness, material, color
  u1.SetLayer(1, fgkCarlosU1Y - fgkCarlosCardChipSiThick, plastiChip, kGray+3);
  u1.SetInitialNode( (TGeoVolume *) assemblySupCarlos);

  p1[0] = fgkCarlosU1posX - fgkCarlosU1X/2;
  p1[1] = shiftGlobY - fgkCarlosCard2HeatBridge + fgkCarlosCardY1/2 + fgkCarlosU1Y/2;
  p1[2] = fgkCarlosCardShift + fgkCarlosU1posZ;

  p2[0] = fgkCarlosU1posX + fgkCarlosU1X/2;
  p2[1] = shiftGlobY - fgkCarlosCard2HeatBridge + fgkCarlosCardY1/2 + fgkCarlosU1Y/2;
  p2[2] = fgkCarlosCardShift + fgkCarlosU1posZ;
  u1.AddCheckPoint( (TGeoVolume *) assemblySupCarlos, 0, p1, vX);
  u1.AddCheckPoint( (TGeoVolume *) assemblySupCarlos, 1, p2, vX);
  u1.CreateAndInsertBoxCableSegment(1,90);

  //---
  AliITSv11GeomCableFlat u2("carlosCardU2", fgkCarlosU2Z, fgkCarlosU2Y); // name, width, thickness
  u2.SetNLayers(2);
  u2.SetLayer(0, fgkCarlosCardChipSiThick, siliconChip, kGreen); // index, thickness, material, color
  u2.SetLayer(1, fgkCarlosU2Y - fgkCarlosCardChipSiThick, plastiChip, kGray+3);
  u2.SetInitialNode( (TGeoVolume *) assemblySupCarlos);

  p1[0] = fgkCarlosU2posX - fgkCarlosU2X/2;
  p1[1] = shiftGlobY - fgkCarlosCard2HeatBridge + fgkCarlosCardY1/2 + fgkCarlosU2Y/2;
  p1[2] = fgkCarlosCardShift + fgkCarlosU2posZ;

  p2[0] = fgkCarlosU2posX + fgkCarlosU2X/2;
  p2[1] = shiftGlobY - fgkCarlosCard2HeatBridge + fgkCarlosCardY1/2 + fgkCarlosU2Y/2;
  p2[2] = fgkCarlosCardShift + fgkCarlosU2posZ;
  u2.AddCheckPoint( (TGeoVolume *) assemblySupCarlos, 0, p1, vX);
  u2.AddCheckPoint( (TGeoVolume *) assemblySupCarlos, 1, p2, vX);
  u2.CreateAndInsertBoxCableSegment(1,90);

  //---
  AliITSv11GeomCableFlat u3("carlosCardU3", fgkCarlosU3Z, fgkCarlosU3Y); // name, width, thickness
  u3.SetNLayers(2);
  u3.SetLayer(0, fgkCarlosCardChipSiThick, siliconChip, kGreen); // index, thickness, material, color
  u3.SetLayer(1, fgkCarlosU3Y - fgkCarlosCardChipSiThick, plastiChip, kGray+3);
  u3.SetInitialNode( (TGeoVolume *) assemblySupCarlos);

  Double_t u3Y = shiftGlobY - fgkCarlosCard2HeatBridge + fgkCarlosCardY1/2 + fgkCarlosU3Y/2;
  p1[0] = fgkCarlosU3posX - fgkCarlosU3X/2;
  p1[1] = u3Y;
  p1[2] = fgkCarlosCardShift + fgkCarlosU3posZ;

  p2[0] = fgkCarlosU3posX + fgkCarlosU3X/2;
  p2[1] = u3Y;
  p2[2] = fgkCarlosCardShift + fgkCarlosU3posZ;
  u3.AddCheckPoint( (TGeoVolume *) assemblySupCarlos, 0, p1, vX);
  u3.AddCheckPoint( (TGeoVolume *) assemblySupCarlos, 1, p2, vX);
  TGeoVolume *u3Vol = u3.CreateAndInsertBoxCableSegment(1,90);

  //--- U4 is like U3 (?)
  TGeoCombiTrans *u4Trans = new TGeoCombiTrans;
  u4Trans->RotateX(90);
  u4Trans->SetTranslation(fgkCarlosU4posX, u3Y,
			  fgkCarlosCardShift + fgkCarlosU4posZ);
  assemblySupCarlos->AddNode(u3Vol, 2, u4Trans);
						 
  //---
  AliITSv11GeomCableFlat u17("carlosCardU17", fgkCarlosU17Z, fgkCarlosU17Y); // name, width, thickness
  u17.SetNLayers(2);
  u17.SetLayer(0, fgkCarlosCardChipSiThick, siliconChip, kGreen); // index, thickness, material, color
  u17.SetLayer(1, fgkCarlosU17Y - fgkCarlosCardChipSiThick, plastiChip, kGray+3);
  u17.SetInitialNode( (TGeoVolume *) assemblySupCarlos);

  p1[0] = fgkCarlosU17posX - fgkCarlosU17X/2;
  p1[1] = shiftGlobY - fgkCarlosCard2HeatBridge + fgkCarlosCardY1/2 + fgkCarlosU17Y/2;
  p1[2] = fgkCarlosCardShift + fgkCarlosU17posZ;

  p2[0] = fgkCarlosU17posX + fgkCarlosU17X/2;
  p2[1] = shiftGlobY - fgkCarlosCard2HeatBridge + fgkCarlosCardY1/2 + fgkCarlosU17Y/2;
  p2[2] = fgkCarlosCardShift + fgkCarlosU17posZ;
  u17.AddCheckPoint( (TGeoVolume *) assemblySupCarlos, 0, p1, vX);
  u17.AddCheckPoint( (TGeoVolume *) assemblySupCarlos, 1, p2, vX);
  u17.CreateAndInsertBoxCableSegment(1,90);

  //---
  AliITSv11GeomCableFlat u35("carlosCardU35", fgkCarlosU35Z, fgkCarlosU35Y); // name, width, thickness
  u35.SetNLayers(2);
  u35.SetLayer(0, fgkCarlosCardChipSiThick, siliconChip, kGreen); // index, thickness, material, color
  u35.SetLayer(1, fgkCarlosU35Y - fgkCarlosCardChipSiThick, plastiChip, kGray+3);
  u35.SetInitialNode( (TGeoVolume *) assemblySupCarlos);

  p1[0] = fgkCarlosU35posX - fgkCarlosU35X/2;
  p1[1] = shiftGlobY - fgkCarlosCard2HeatBridge + fgkCarlosCardY1/2 + fgkCarlosU35Y/2;
  p1[2] = fgkCarlosCardShift + fgkCarlosU35posZ;

  p2[0] = fgkCarlosU35posX + fgkCarlosU35X/2;
  p2[1] = shiftGlobY - fgkCarlosCard2HeatBridge + fgkCarlosCardY1/2 + fgkCarlosU35Y/2;
  p2[2] = fgkCarlosCardShift + fgkCarlosU35posZ;
  u35.AddCheckPoint( (TGeoVolume *) assemblySupCarlos, 0, p1, vX);
  u35.AddCheckPoint( (TGeoVolume *) assemblySupCarlos, 1, p2, vX);
  u35.CreateAndInsertBoxCableSegment(1,90);

  //---
  AliITSv11GeomCableFlat u36("carlosCardU36", fgkCarlosU36Z, fgkCarlosU36Y); // name, width, thickness
  u36.SetNLayers(2);
  u36.SetLayer(0, fgkCarlosCardChipSiThick, siliconChip, kGreen); // index, thickness, material, color
  u36.SetLayer(1, fgkCarlosU36Y - fgkCarlosCardChipSiThick, plastiChip, kGray+3);
  u36.SetInitialNode( (TGeoVolume *) assemblySupCarlos);

  p1[0] = fgkCarlosU36posX - fgkCarlosU36X/2;
  p1[1] = shiftGlobY - fgkCarlosCard2HeatBridge + fgkCarlosCardY1/2 + fgkCarlosU36Y/2;
  p1[2] = fgkCarlosCardShift + fgkCarlosU36posZ;

  p2[0] = fgkCarlosU36posX + fgkCarlosU36X/2;
  p2[1] = shiftGlobY - fgkCarlosCard2HeatBridge + fgkCarlosCardY1/2 + fgkCarlosU36Y/2;
  p2[2] = fgkCarlosCardShift + fgkCarlosU36posZ;
  u36.AddCheckPoint( (TGeoVolume *) assemblySupCarlos, 0, p1, vX);
  u36.AddCheckPoint( (TGeoVolume *) assemblySupCarlos, 1, p2, vX);
  u36.CreateAndInsertBoxCableSegment(1,90);

  //--- QZ1
  AliITSv11GeomCableFlat qz1("carlosCardQZ1", fgkCarlosQZ1Z, fgkCarlosQZ1Y); // name, width, thickness
  qz1.SetNLayers(2);
  qz1.SetLayer(0, fgkCarlosCardChipSiThick, siliconChip, kGreen); // index, thickness, material, color
  qz1.SetLayer(1, fgkCarlosQZ1Y - fgkCarlosCardChipSiThick, plastiChip, kGray+3);
  qz1.SetInitialNode( (TGeoVolume *) assemblySupCarlos);

  p1[0] = fgkCarlosQZ1posX - fgkCarlosQZ1X/2;
  p1[1] = shiftGlobY - fgkCarlosCard2HeatBridge + fgkCarlosCardY1/2 + fgkCarlosQZ1Y/2;
  p1[2] = fgkCarlosCardShift + fgkCarlosQZ1posZ;

  p2[0] = fgkCarlosQZ1posX + fgkCarlosQZ1X/2;
  p2[1] = shiftGlobY - fgkCarlosCard2HeatBridge + fgkCarlosCardY1/2 + fgkCarlosQZ1Y/2;
  p2[2] = fgkCarlosCardShift + fgkCarlosQZ1posZ;
  qz1.AddCheckPoint( (TGeoVolume *) assemblySupCarlos, 0, p1, vX);
  qz1.AddCheckPoint( (TGeoVolume *) assemblySupCarlos, 1, p2, vX);
  qz1.CreateAndInsertBoxCableSegment(1,90);

  return assemblySupCarlos;
}

//________________________________________________________________________
Int_t AliITSv11GeometrySDD::CreateLVCard() {
  // 
  // Creates the assemblies containing the LV cards (left and right)
  //

  TGeoMedium *glassFiber  = GetMedium("SDD SI CHIP$");// glassFiber   TO CODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoMedium *siliconChip = GetMedium("SDD SI CHIP$");// ITSsddSiChip
  TGeoMedium *plastiChip  = GetMedium("SDDKAPTON (POLYCH2)$"); // ITSsddKAPTON_POLYCH2
  TGeoMedium *copper      = GetMedium("COPPER$"); 
  TGeoMedium *alCu12SDD   = GetMedium("INOX$"); // ITSsddAlCu12,  to code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoMedium *stainless   = GetMedium("INOX$"); // for screws, what is the material ???????????

  fCardLVL = new TGeoVolumeAssembly("ITSsddLVCardLeft");
  fCardLVR = new TGeoVolumeAssembly("ITSsddLVCardRight");

  // we are going to use flat cable class to create multilayer box,
  // then we can use the pointers to created volumes to place them elsewhere  
  Double_t p1[3], p2[3], vX[3] = {1,0,0};

  Double_t carLVfullThick = fgkLVcardZ+fgkLVcardCuZ;
  AliITSv11GeomCableFlat cardLV("cardLV", fgkLVcardY, carLVfullThick); // name, width, thickness
  cardLV.SetNLayers(2);
  cardLV.SetLayer(0, fgkLVcardCuZ, copper,     30); // index, thickness, material, color
  cardLV.SetLayer(1, fgkLVcardZ,   glassFiber, 30);
  cardLV.SetInitialNode( (TGeoVolume *) fCardLVL);
  p1[0] = 0;
  p1[1] = fgkLVcardY/2;
  p1[2] = 0;
  p2[0] = fgkLVcardX;
  p2[1] = fgkLVcardY/2;
  p2[2] = 0;
  cardLV.AddCheckPoint( (TGeoVolume *) fCardLVL, 0, p1, vX);
  cardLV.AddCheckPoint( (TGeoVolume *) fCardLVL, 1, p2, vX);
  TGeoVolume* boxVol = cardLV.CreateAndInsertBoxCableSegment(1);
  TGeoRotation *rotAdd = new TGeoRotation("",90,0,0);
  TGeoCombiTrans *trCard = new TGeoCombiTrans(-fgkLVcardX/2,fgkLVcardY/2,0,rotAdd);
  fCardLVR->AddNode(boxVol, 1, trCard);

  Double_t chip0fullThick = fgkLVChip0Z + fgkLVChip0SiZ;
  AliITSv11GeomCableFlat chipO("chipO", fgkLVChip0Y, chip0fullThick); // name, width, thickness
  chipO.SetNLayers(2);
  chipO.SetLayer(0, fgkLVChip0SiZ, siliconChip, 8); // index, thickness, material, color
  chipO.SetLayer(1, fgkLVChip0Z,   plastiChip, 12);
  chipO.SetInitialNode( (TGeoVolume *) fCardLVL);
  p1[0] = (fgkLVChip0PosX - fgkLVChip0X/2);
  p1[1] = fgkLVChip0PosY;
  p1[2] = carLVfullThick/2 + chip0fullThick/2;

  p2[0] = (fgkLVChip0PosX + fgkLVChip0X/2);
  p2[1] = fgkLVChip0PosY;
  p2[2] = carLVfullThick/2 + chip0fullThick/2;
  chipO.AddCheckPoint( (TGeoVolume *) fCardLVL, 0, p1, vX);
  chipO.AddCheckPoint( (TGeoVolume *) fCardLVL, 1, p2, vX);
  boxVol = chipO.CreateAndInsertBoxCableSegment(1);
  trCard = new TGeoCombiTrans( -fgkLVChip0PosX,
			       fgkLVChip0PosY,
			       carLVfullThick/2+chip0fullThick/2, rotAdd);
  fCardLVR->AddNode(boxVol, 1, trCard);

  // put also this chip on the other side of the card
  trCard = new TGeoCombiTrans( fgkLVChip0PosX,
			       fgkLVChip0PosY,
			       -carLVfullThick/2-chip0fullThick/2, rotAdd);
  fCardLVL->AddNode(boxVol, 2, trCard);
  trCard = new TGeoCombiTrans( -fgkLVChip0PosX,
			       fgkLVChip0PosY,
			       -carLVfullThick/2-chip0fullThick/2, rotAdd);
  fCardLVR->AddNode(boxVol, 2, trCard);

  Double_t chip1fullThick = fgkLVChip1Z + fgkLVChip1SiZ;
  AliITSv11GeomCableFlat chip1("chip1", fgkLVChip1Y, chip1fullThick);
  chip1.SetNLayers(2);
  chip1.SetLayer(0, fgkLVChip1SiZ, siliconChip, 8);
  chip1.SetLayer(1, fgkLVChip1Z,   plastiChip, 12);
  chip1.SetInitialNode( (TGeoVolume *) fCardLVL);
  p1[0] = (fgkLVChip1PosX-fgkLVChip1X/2);
  p1[1] = fgkLVChip1PosY;
  p1[2] = carLVfullThick/2 + chip1fullThick/2;

  p2[0] = (fgkLVChip1PosX+fgkLVChip1X/2);
  p2[1] = fgkLVChip1PosY;
  p2[2] = carLVfullThick/2 + chip1fullThick/2;
  chip1.AddCheckPoint( (TGeoVolume *) fCardLVL, 0, p1, vX);
  chip1.AddCheckPoint( (TGeoVolume *) fCardLVL, 1, p2, vX);
  boxVol = chip1.CreateAndInsertBoxCableSegment(1);
  trCard = new TGeoCombiTrans( -fgkLVChip1PosX,
			       fgkLVChip1PosY,
			       carLVfullThick/2 + chip1fullThick/2, rotAdd);
  fCardLVR->AddNode(boxVol, 1, trCard);

  Double_t chip2fullThick = fgkLVChip2Z + fgkLVChip2SiZ;
  AliITSv11GeomCableFlat chip2("chip2", fgkLVChip2Y, chip2fullThick);
  chip2.SetNLayers(2);
  chip2.SetLayer(0, fgkLVChip2SiZ, siliconChip, 8);
  chip2.SetLayer(1, fgkLVChip2Z,   plastiChip, 12);
  chip2.SetInitialNode( (TGeoVolume *) fCardLVL);
  p1[0] = (fgkLVChip2PosX-fgkLVChip2X/2);
  p1[1] = fgkLVChip2PosY;
  p1[2] = carLVfullThick/2 + chip2fullThick/2;
  p2[0] = (fgkLVChip2PosX+fgkLVChip2X/2);
  p2[1] = fgkLVChip2PosY;
  p2[2] = carLVfullThick/2 + chip2fullThick/2;
  chip2.AddCheckPoint( (TGeoVolume *) fCardLVL, 0, p1, vX);
  chip2.AddCheckPoint( (TGeoVolume *) fCardLVL, 1, p2, vX);
  boxVol = chip2.CreateAndInsertBoxCableSegment(1);
  trCard = new TGeoCombiTrans( -fgkLVChip2PosX,
			       fgkLVChip2PosY,
			       carLVfullThick/2 + chip2fullThick/2, rotAdd);
  fCardLVR->AddNode(boxVol, 1, trCard);

  Double_t chip3fullThick = fgkLVChip3Z + fgkLVChip3SiZ;
  AliITSv11GeomCableFlat chip3("chip3", fgkLVChip3Y, chip3fullThick);
  chip3.SetNLayers(2);
  chip3.SetLayer(0, fgkLVChip3Z,   plastiChip, 12);
  chip3.SetLayer(1, fgkLVChip3SiZ, siliconChip, 8);
  chip3.SetInitialNode( (TGeoVolume *) fCardLVL);
  p1[0] = (fgkLVChip3PosX-fgkLVChip3X/2);
  p1[1] = fgkLVChip3PosY;
  p1[2] = -carLVfullThick/2 - chip3fullThick/2;
  p2[0] = (fgkLVChip3PosX+fgkLVChip3X/2);
  p2[1] = fgkLVChip3PosY;
  p2[2] = -carLVfullThick/2 - chip3fullThick/2;
  chip3.AddCheckPoint( (TGeoVolume *) fCardLVL, 0, p1, vX);
  chip3.AddCheckPoint( (TGeoVolume *) fCardLVL, 1, p2, vX);
  boxVol = chip3.CreateAndInsertBoxCableSegment(1);
  trCard = new TGeoCombiTrans( -fgkLVChip3PosX,
			       fgkLVChip3PosY,
			       -carLVfullThick/2 - chip3fullThick/2, rotAdd);
  fCardLVR->AddNode(boxVol, 1, trCard);

  // the Al pieces for heat exchange :
  TGeoBBox *alLVcooling1 = new TGeoBBox("alLVcooling1" ,
			   fgkLVcoolX1/2, fgkLVcoolY1/2, fgkLVcoolZ1/2);

  TGeoTranslation *alLVcooling1Tr = new TGeoTranslation("alLVcooling1Tr",
					(fgkLVcoolX1/2+fgkLVcoolX2),
					fgkLVcoolPosY+fgkLVcoolY1/2,
			  carLVfullThick/2+chip0fullThick+fgkLVcoolZ1/2);
  TGeoTranslation *alLVcooling1TrB = new TGeoTranslation("alLVcooling1TrB",
				         (fgkLVcoolX1/2+fgkLVcoolX2),
					 fgkLVcoolPosY+fgkLVcoolY1/2,
			  -(carLVfullThick/2+chip0fullThick+fgkLVcoolZ1/2));

  TGeoVolume *vAlLVcooling1 = new TGeoVolume("vAlLVcooling1",alLVcooling1,
					     alCu12SDD);
  vAlLVcooling1->SetLineColor(2);

  //--
  TGeoBBox * alLVcooling2 = new TGeoBBox("lLVcooling2" ,
				 fgkLVcoolX2/2, fgkLVcoolY2/2, fgkLVcoolZ2/2);
  TGeoTranslation *alLVcooling2Tr = new TGeoTranslation("alLVcooling2Tr",
							(fgkLVcoolX2/2),
					fgkLVcoolPosY+fgkLVcoolY1/2,
		     carLVfullThick/2+chip0fullThick+fgkLVcoolZ1-fgkLVcoolZ2/2);
  TGeoTranslation *alLVcooling2TrB = new TGeoTranslation("alLVcooling2TrB",
							 (fgkLVcoolX2/2),
					fgkLVcoolPosY+fgkLVcoolY1/2,
		   -(carLVfullThick/2+chip0fullThick+fgkLVcoolZ1-fgkLVcoolZ2/2));

  TGeoVolume *vAlLVcooling2 = new TGeoVolume("vAlLVcooling2",alLVcooling2,
					     alCu12SDD);
  vAlLVcooling2->SetLineColor(2);

  //--
  Double_t alLVcoolZ3 = (fgkLVcardCuZ+fgkLVcardZ+2.*(fgkLVChip0SiZ+fgkLVChip0Z)
			  +fgkLVcoolZ1*2.);
  TGeoBBox * alLVcooling3 = new TGeoBBox("lLVcooling3" ,
			   fgkLVcoolX3/2, fgkLVcoolY3/2, alLVcoolZ3/2);
  TGeoTranslation *alLVcooling3Tr = new TGeoTranslation("alLVcooling3Tr",
							(-fgkLVcoolX3/2),
					fgkLVcoolPosY+fgkLVcoolY1-fgkLVcoolY3/2,
					0);
  TGeoVolume *vAlLVcooling3 = new TGeoVolume("vAlLVcooling3",alLVcooling3,alCu12SDD);
  vAlLVcooling3->SetLineColor(2);

  //=== screw fixing th LV card to the U cooling tube :
  TGeoTube *littleScrew = new TGeoTube("littleScrewLV", 0, fgkLittleScrewR,
				       fgkLVcoolY3/2);
  TGeoVolume *vLittleScrew = new TGeoVolume("vLittleScrewLV",
					    littleScrew, stainless);
  TGeoRotation *rotScrew = new TGeoRotation("",0,90,0);

  TGeoCombiTrans *cbScrew = new TGeoCombiTrans(0,0,fgkShiftLittleScrewLV,
					       rotScrew);
  vAlLVcooling3->AddNode(vLittleScrew, 1, cbScrew);

  TGeoTube *littleScrewHead = new TGeoTube("littleScrewLVhead",
					   0, fgkLittleLVScrewHeadR,
					   fgkLittleScrewHeadH/2);
  TGeoVolume *vLittleScrewHead = new TGeoVolume("vLittleScrewLVhead",
						littleScrewHead, stainless);
  vLittleScrewHead->SetLineColor(kGray);
  TGeoCombiTrans *cbScrewHeadL = new TGeoCombiTrans( -fgkLVcoolX3/2,
						     fgkLVcoolPosY+fgkLVcoolY1 + fgkLittleScrewHeadH/2,
						     fgkShiftLittleScrewLV,
						     rotScrew);
  fCardLVL->AddNode(vLittleScrewHead, 1, cbScrewHeadL);

  TGeoCombiTrans *cbScrewHeadR = new TGeoCombiTrans( fgkLVcoolX3/2,
						     fgkLVcoolPosY+fgkLVcoolY1 + fgkLittleScrewHeadH/2,
						     fgkShiftLittleScrewLV,
						     rotScrew);
  fCardLVR->AddNode(vLittleScrewHead, 1, cbScrewHeadR);

  // adding the cooling pieces to the left card
  fCardLVL->AddNode(vAlLVcooling1, 1,alLVcooling1Tr);
  fCardLVL->AddNode(vAlLVcooling1, 2,alLVcooling1TrB);
  fCardLVL->AddNode(vAlLVcooling2, 1,alLVcooling2Tr);
  fCardLVL->AddNode(vAlLVcooling2, 2,alLVcooling2TrB);
  fCardLVL->AddNode(vAlLVcooling3, 1,alLVcooling3Tr);

  TGeoTranslation *alLVcooling1TrR = new TGeoTranslation("alLVcooling1TrR",
					 -(fgkLVcoolX1/2+fgkLVcoolX2),
					 fgkLVcoolPosY+fgkLVcoolY1/2,
				     carLVfullThick/2+chip0fullThick+fgkLVcoolZ1/2);
  TGeoTranslation *alLVcooling1TrBR = new TGeoTranslation("alLVcooling1TrBR",
					  -(fgkLVcoolX1/2+fgkLVcoolX2),
					  fgkLVcoolPosY+fgkLVcoolY1/2,
				     -(carLVfullThick/2+chip0fullThick+fgkLVcoolZ1/2));
  TGeoTranslation *alLVcooling2TrR = new TGeoTranslation("alLVcooling2TrR",
							 -(fgkLVcoolX2/2),
					 fgkLVcoolPosY+fgkLVcoolY1/2,
		     carLVfullThick/2+chip0fullThick+fgkLVcoolZ1-fgkLVcoolZ2/2);
  TGeoTranslation *alLVcooling2TrBR = new TGeoTranslation("alLVcooling2TrBR",
							  -(fgkLVcoolX2/2),
					fgkLVcoolPosY+fgkLVcoolY1/2,
		   -(carLVfullThick/2+chip0fullThick+fgkLVcoolZ1-fgkLVcoolZ2/2));

  TGeoTranslation *alLVcooling3TrR = new TGeoTranslation("alLVcooling3TrR",
					fgkLVcoolX3/2,
					fgkLVcoolPosY+fgkLVcoolY1-fgkLVcoolY3/2,
					0);
  // and to the right card
  fCardLVR->AddNode(vAlLVcooling1, 1,alLVcooling1TrR);
  fCardLVR->AddNode(vAlLVcooling1, 2,alLVcooling1TrBR);
  fCardLVR->AddNode(vAlLVcooling2, 1,alLVcooling2TrR);
  fCardLVR->AddNode(vAlLVcooling2, 2,alLVcooling2TrBR);
  fCardLVR->AddNode(vAlLVcooling3, 1,alLVcooling3TrR);

  return kTRUE;
}

//________________________________________________________________________
TGeoVolumeAssembly*  AliITSv11GeometrySDD::CreateHVCard(Int_t iLay){
  // 
  // return an assembly containing the HV card
  //
  iLay = iLay;

  TGeoMedium *ceramic          = GetMedium("CERAMICS$"); // ceramicHVcard
  TGeoMedium *medSMDcapaMiddle = GetMedium("SDD X7R capacitors$");      //    TO CODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoMedium *medSMDcapaEnd    = GetMedium("SDD X7R capacitors$");      // SDDX7RcapacitorsSDD   TO CODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoMedium *stainless        = GetMedium("INOX$");       // ITSspdStainlesSteal ???????????
  TGeoMedium *plastic          = GetMedium("SDDKAPTON (POLYCH2)$");  // ITS_ITSsddKAPTON_POLYCH2 ???????????
  TGeoMedium *alCu12SDD       = GetMedium("INOX$"); // ITSsddAlCu12  : to code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  TGeoVolumeAssembly *highVCard = new TGeoVolumeAssembly("ITSsddHVCard");

  //====================================
  //--- the card itself
  TGeoBBox *ceramicCard = new TGeoBBox("ceramCard", fgkHVCardCeramX/2,
				       fgkHVCardCeramY/2, fgkHVCardCeramZ/2);
  TGeoVolume *vCeramicCard = new TGeoVolume("vCeramCard", ceramicCard, ceramic);
  vCeramicCard->SetLineColor(38);// or 9 blue slightly dark 

  highVCard->AddNode(vCeramicCard, 1, 0);


  //====================================
  //--- capacitors

  // capa1
  TGeoBBox *capa1Middle = new TGeoBBox("cardHVCapa1Middle", fgkHVCardCapa1X/2,
				       fgkHVCardCapa1Ymid/2, fgkHVCardCapa1Z/2);
  TGeoVolume *vCapa1Middle = new TGeoVolume("vCardHVCapa1Middle",capa1Middle,
					    medSMDcapaMiddle);

  TGeoBBox *capa1End = new TGeoBBox("cardHVCapa1End", fgkHVCardCapa1X/2,
				    fgkHVCardCapa1Yend/2, fgkHVCardCapa1Z/2);
  TGeoVolume *vCapa1End = new TGeoVolume("vCardHVCapa1End",capa1End,
					 medSMDcapaEnd);
  vCapa1End->SetLineColor(18);// grey silver
  TGeoTranslation *capa1EndTr1 = new TGeoTranslation("cardHVcapa1EndTr1", 0,
				     (fgkHVCardCapa1Ymid+fgkHVCardCapa1Yend)/2,0);
  TGeoTranslation *capa1EndTr2 = new TGeoTranslation("cardHVcapa1EndTr2", 0,
				     -(fgkHVCardCapa1Ymid+fgkHVCardCapa1Yend)/2,0);

  TGeoTranslation *capa1PosTr = new TGeoTranslation("cardHVcapa1PosTr",
				    fgkHVCardCapa1PosX, fgkHVCardCapa1PosY,
				    -fgkHVCardCeramZ/2-fgkHVCardCapa1Z/2);

  TGeoVolumeAssembly *capa1  = new TGeoVolumeAssembly("cardHVCapa1");
  capa1->AddNode(vCapa1Middle, 1,0);
  capa1->AddNode(vCapa1End, 1, capa1EndTr1);
  capa1->AddNode(vCapa1End, 2, capa1EndTr2);

  highVCard->AddNode(capa1, 1, capa1PosTr);

  // capa2
  TGeoBBox *capa2Middle = new TGeoBBox("cardHVCapa2Middle", fgkHVCardCapa2X/2,
				       fgkHVCardCapa2Ymid/2, fgkHVCardCapa2Z/2);
  TGeoVolume *vCapa2Middle = new TGeoVolume("vCardHVCapa2Middle",capa2Middle,
					    medSMDcapaMiddle);

  TGeoBBox *capa2End = new TGeoBBox("cardHVCapa2End", fgkHVCardCapa2X/2,
				    fgkHVCardCapa2Yend/2, fgkHVCardCapa2Z/2);
  TGeoVolume *vCapa2End = new TGeoVolume("vCardHVCapa2End",capa2End,
					 medSMDcapaEnd);
  vCapa2End->SetLineColor(18);// grey silver
  TGeoTranslation *capa2EndTr1 = new TGeoTranslation("cardHVcapa2EndTr1", 0,
				     (fgkHVCardCapa2Ymid+fgkHVCardCapa2Yend)/2,0);
  TGeoTranslation *capa2EndTr2 = new TGeoTranslation("cardHVcapa2EndTr2", 0,
				     -(fgkHVCardCapa2Ymid+fgkHVCardCapa2Yend)/2,0);

  TGeoTranslation *capa2PosTr = new TGeoTranslation("cardHVcapa2PosTr",
				    fgkHVCardCapa2PosX, fgkHVCardCapa2PosY,
				    -fgkHVCardCeramZ/2-fgkHVCardCapa2Z/2);

  TGeoVolumeAssembly *capa2  = new TGeoVolumeAssembly("cardHVCapa2");
  capa2->AddNode(vCapa2Middle, 1,0);
  capa2->AddNode(vCapa2End, 1, capa2EndTr1);
  capa2->AddNode(vCapa2End, 2, capa2EndTr2);

  highVCard->AddNode(capa2, 1, capa2PosTr);

  // capa3
  TGeoBBox *capa3Middle = new TGeoBBox("cardHVCapa3Middle", fgkHVCardCapa3Xmid/2,
				       fgkHVCardCapa3Y/2, fgkHVCardCapa3Z/2);
  TGeoVolume *vCapa3Middle = new TGeoVolume("vCardHVCapa3Middle",capa3Middle,
					    medSMDcapaMiddle);

  TGeoBBox *capa3End = new TGeoBBox("cardHVCapa3End", fgkHVCardCapa3Xend/2,
				    fgkHVCardCapa3Y/2, fgkHVCardCapa3Z/2);
  TGeoVolume *vCapa3End = new TGeoVolume("vCardHVCapa3End",capa3End,
					 medSMDcapaEnd);
  vCapa3End->SetLineColor(18);// grey silver

  TGeoTranslation *capa3EndTr1 = new TGeoTranslation("cardHVcapa3EndTr1",
				 (fgkHVCardCapa3Xmid+fgkHVCardCapa3Xend)/2, 0, 0);
  TGeoTranslation *capa3EndTr2 = new TGeoTranslation("cardHVcapa2EndTr2",
				 -(fgkHVCardCapa3Xmid+fgkHVCardCapa3Xend)/2, 0, 0);

  TGeoVolumeAssembly *capa3  = new TGeoVolumeAssembly("cardHVCapa3");
  capa3->AddNode(vCapa3Middle, 1,0);
  capa3->AddNode(vCapa3End, 1, capa3EndTr1);
  capa3->AddNode(vCapa3End, 2, capa3EndTr2);

  TGeoTranslation *capa3PosTr1 = new TGeoTranslation("cardHVcapa3PosTr1",
				    fgkHVCardCapa3PosX1, fgkHVCardCapa3PosY1,
				    -fgkHVCardCeramZ/2-fgkHVCardCapa3Z/2);

  TGeoTranslation *capa3PosTr2 = new TGeoTranslation("cardHVcapa3PosTr2",
				    fgkHVCardCapa3PosX2, fgkHVCardCapa3PosY1,
				    -fgkHVCardCeramZ/2-fgkHVCardCapa3Z/2);

  TGeoTranslation *capa3PosTr3 = new TGeoTranslation("cardHVcapa3PosTr3",
				    fgkHVCardCapa3PosX3, fgkHVCardCapa3PosY2,
				    -fgkHVCardCeramZ/2-fgkHVCardCapa3Z/2);

  TGeoTranslation *capa3PosTr4 = new TGeoTranslation("cardHVcapa3PosTr4",
				    fgkHVCardCapa3PosX4, fgkHVCardCapa3PosY2,
				    -fgkHVCardCeramZ/2-fgkHVCardCapa3Z/2);

  TGeoTranslation *capa3PosTr5 = new TGeoTranslation("cardHVcapa3PosTr5",
				    fgkHVCardCapa3PosX5, fgkHVCardCapa3PosY3,
				    -fgkHVCardCeramZ/2-fgkHVCardCapa3Z/2);

  highVCard->AddNode(capa3, 1, capa3PosTr1);
  highVCard->AddNode(capa3, 2, capa3PosTr2);
  highVCard->AddNode(capa3, 3, capa3PosTr3);
  highVCard->AddNode(capa3, 4, capa3PosTr4);
  highVCard->AddNode(capa3, 5, capa3PosTr5);

  //====================================
  //--- connexions to LV card

  Double_t fgkConnexLVHVdiam1 =  0.8*fgkmm;
  Double_t fgkConnexLVHVdiam2 =  2*fgkmm;
  Double_t fgkConnexLVHVlen   =  6.2*fgkmm;
  Double_t fgkConnexLVHVx     =  3*fgkmm;
  Double_t fgkConnexLVHVy1    =  8*fgkmm;
  Double_t fgkConnexLVHVdy    =  2.5*fgkmm;

  TGeoTube *connexLVHVmetal = new TGeoTube("connexLVHVmetal",0,
				  fgkConnexLVHVdiam1/2,fgkConnexLVHVlen/2);
  TGeoTube *connexLVHVplastic = new TGeoTube("connexLVHVplastic",
					     fgkConnexLVHVdiam1/2,
					     fgkConnexLVHVdiam2/2,
					     fgkConnexLVHVlen/2);
  TGeoVolume *vConnexLVHVmetal = new TGeoVolume("ITSsddConnexLVHVmetal",
						connexLVHVmetal, stainless);
  TGeoVolume *vConnexLVHVplast = new TGeoVolume("ITSsddConnexLVHVplast",
						connexLVHVplastic, plastic);
  vConnexLVHVmetal->SetLineColor(10);// white
  vConnexLVHVplast->SetLineColor(12);  // dark grey

  TGeoVolumeAssembly *connexion = new TGeoVolumeAssembly("ITSsddConnexLVHV");
  connexion->AddNode(vConnexLVHVmetal, 1, 0);
  connexion->AddNode(vConnexLVHVplast, 1, 0);

  TGeoTranslation *trConnexion1 = new TGeoTranslation(-fgkConnexLVHVx,fgkConnexLVHVy1,
				       -fgkHVCardCeramZ/2-fgkConnexLVHVlen/2 );
  TGeoTranslation *trConnexion2 = new TGeoTranslation( fgkConnexLVHVx,fgkConnexLVHVy1,
				       -fgkHVCardCeramZ/2-fgkConnexLVHVlen/2 );

  TGeoTranslation *trConnexion3 = new TGeoTranslation(-fgkConnexLVHVx,
					fgkConnexLVHVy1+fgkConnexLVHVdy,
				       -fgkHVCardCeramZ/2-fgkConnexLVHVlen/2 );
  TGeoTranslation *trConnexion4 = new TGeoTranslation( fgkConnexLVHVx,
					fgkConnexLVHVy1+fgkConnexLVHVdy,
				       -fgkHVCardCeramZ/2-fgkConnexLVHVlen/2 );

  TGeoTranslation *trConnexion5 = new TGeoTranslation(-fgkConnexLVHVx,
					fgkConnexLVHVy1+2*fgkConnexLVHVdy,
				       -fgkHVCardCeramZ/2-fgkConnexLVHVlen/2 );
  TGeoTranslation *trConnexion6 = new TGeoTranslation( fgkConnexLVHVx,
					fgkConnexLVHVy1+2*fgkConnexLVHVdy,
				       -fgkHVCardCeramZ/2-fgkConnexLVHVlen/2 );

  TGeoTranslation *trConnexion7 = new TGeoTranslation(-fgkConnexLVHVx,
					fgkConnexLVHVy1+3*fgkConnexLVHVdy,
				       -fgkHVCardCeramZ/2-fgkConnexLVHVlen/2 );
  TGeoTranslation *trConnexion8 = new TGeoTranslation( fgkConnexLVHVx,
					fgkConnexLVHVy1+3*fgkConnexLVHVdy,
				       -fgkHVCardCeramZ/2-fgkConnexLVHVlen/2 );

  highVCard->AddNode(connexion, 1, trConnexion1);
  highVCard->AddNode(connexion, 2, trConnexion2);
  highVCard->AddNode(connexion, 3, trConnexion3);
  highVCard->AddNode(connexion, 4, trConnexion4);
  highVCard->AddNode(connexion, 5, trConnexion5);
  highVCard->AddNode(connexion, 6, trConnexion6);
  highVCard->AddNode(connexion, 7, trConnexion7);
  highVCard->AddNode(connexion, 8, trConnexion8);

  //====================================
  //--- cooling pieces

  TGeoBBox *cardHVcool1 = new TGeoBBox("cardHVcool1",fgkHVCardCool1X/2,
				       fgkHVCardCool1Y/2, fgkHVCardCool1Z/2);


  TGeoBBox *cardHVcool2 = new TGeoBBox("cardHVcool2",fgkHVCardCool2X/2,
				       fgkHVCardCool2Y/2, fgkHVCardCool2Z/2);

  TGeoBBox *cardHVcool3 = new TGeoBBox("cardHVcool3",fgkHVCardCool3X/2,
				       fgkHVCardCool3Y/2, fgkHVCardCool3Z/2);

  TGeoVolume *vCardHVcool1 = new TGeoVolume("vCardHVcool1",cardHVcool1,
					    alCu12SDD);
  TGeoVolume *vCardHVcool2 = new TGeoVolume("vCardHVcool2",cardHVcool2,
					    alCu12SDD);
  TGeoVolume *vCardHVcool3 = new TGeoVolume("vCardHVcool3",cardHVcool3,
					    alCu12SDD);
  // This last volume contains the screw used for fixing
  // the card to the cooling tube ...
  TGeoTube *littleScrewHV = new TGeoTube("littleScrewHV", 0, fgkLittleScrewR,
					 fgkHVCardCool3Y/2);
  TGeoVolume *vLittleScrewHV = new TGeoVolume("vLittleScrewHV",
					      littleScrewHV, stainless);

  TGeoRotation *rotScrewHead = new TGeoRotation("",0,90,0);
  vCardHVcool3->AddNode(vLittleScrewHV, 1,rotScrewHead);

  vCardHVcool1->SetLineColor(2); //red
  vCardHVcool2->SetLineColor(2); //red
  vCardHVcool3->SetLineColor(2); //red

  TGeoTranslation *cool1Tr1 = new TGeoTranslation("cardHVcool1Tr1",
					     fgkHVCardCeramX/2-fgkHVCardCool1X/2,
					     -fgkHVCardCoolDY+fgkHVCardCool1Y/2,
					     fgkHVCardCeramZ/2+fgkHVCardCool1Z/2);
  TGeoTranslation *cool1Tr2 = new TGeoTranslation("cardHVcool1Tr2",
					    -fgkHVCardCeramX/2+fgkHVCardCool1X/2,
					     -fgkHVCardCoolDY+fgkHVCardCool1Y/2,
					     fgkHVCardCeramZ/2+fgkHVCardCool1Z/2);

  highVCard->AddNode(vCardHVcool1, 1, cool1Tr1);
  highVCard->AddNode(vCardHVcool1, 2, cool1Tr2);

  TGeoTranslation *cool2Tr1 = new TGeoTranslation("cardHVcool2Tr1",
			          fgkHVCardCeramX/2-fgkHVCardCool1X+fgkHVCardCool2X/2,
					     -fgkHVCardCoolDY-fgkHVCardCool2Y/2,
					     fgkHVCardCeramZ/2+fgkHVCardCool2Z/2);

  TGeoTranslation *cool2Tr2 = new TGeoTranslation("cardHVcool2Tr2",
				 -fgkHVCardCeramX/2+fgkHVCardCool1X-fgkHVCardCool2X/2,
					     -fgkHVCardCoolDY-fgkHVCardCool2Y/2,
					     fgkHVCardCeramZ/2+fgkHVCardCool2Z/2);

  highVCard->AddNode(vCardHVcool2, 1, cool2Tr1);
  highVCard->AddNode(vCardHVcool2, 2, cool2Tr2);

  TGeoTranslation *cool3Tr1 = new TGeoTranslation("cardHVcool2Tr1",
	         fgkHVCardCeramX/2-fgkHVCardCool1X+fgkHVCardCool2X+fgkHVCardCool3X/2,
					     -fgkHVCardCoolDY-fgkHVCardCool3Y/2,
				 fgkHVCardCeramZ/2+fgkHVCardCool2Z-fgkHVCardCool3Z/2);

  TGeoTranslation *cool3Tr2 = new TGeoTranslation("cardHVcool2Tr2",
		-fgkHVCardCeramX/2+fgkHVCardCool1X-fgkHVCardCool2X-fgkHVCardCool3X/2,
					     -fgkHVCardCoolDY-fgkHVCardCool3Y/2,
				 fgkHVCardCeramZ/2+fgkHVCardCool2Z-fgkHVCardCool3Z/2);

  highVCard->AddNode(vCardHVcool3, 1, cool3Tr1);
  highVCard->AddNode(vCardHVcool3, 2, cool3Tr2);

  //====================================
  //--- screws
  TGeoCombiTrans *cbScrewHead1 = new TGeoCombiTrans("cardHVscrewHeadTr1",
		 fgkHVCardCeramX/2-fgkHVCardCool1X+fgkHVCardCool2X+fgkHVCardCool3X/2,
				 -fgkHVCardCoolDY+fgkLittleScrewHeadH/2,
				 fgkHVCardCeramZ/2+fgkHVCardCool2Z-fgkHVCardCool3Z/2,
						    rotScrewHead);
  TGeoCombiTrans *cbScrewHead2 = new TGeoCombiTrans("cardHVscrewHeadTr2",
		-fgkHVCardCeramX/2+fgkHVCardCool1X-fgkHVCardCool2X-fgkHVCardCool3X/2,
				 -fgkHVCardCoolDY+fgkLittleScrewHeadH/2,
				 fgkHVCardCeramZ/2+fgkHVCardCool2Z-fgkHVCardCool3Z/2,
						    rotScrewHead);

  highVCard->AddNode(fCommonVol[0], 1, cbScrewHead1);
  highVCard->AddNode(fCommonVol[0], 2, cbScrewHead2);

  return highVCard;
}


//________________________________________________________________________
TGeoVolumeAssembly*  AliITSv11GeometrySDD::CreateEndLadderCards(Int_t iLay) {
// 
// return an assembly containing the LV, HV and Carlos cards of one ladder
// and their cooling system 
//

  TGeoMedium *alCu12SDD       = GetMedium("AL$"); // ITSsddAlCu12 : to code !!!!!!!!!!!!!!
  TGeoMedium *phynoxSDD       = GetMedium("INOX$");
  TGeoMedium *coolerMediumSDD = GetMedium("WATER$");

  TGeoVolumeAssembly *endLadderCards = new TGeoVolumeAssembly("endLadderCards");

  //=*********************************
  //--- The rounded pipe for the end ladder card coooling

  Double_t endLadPipeUlength = fgkEndLadPipeUlengthLay3;
  Double_t endLadPipeArmZ = fgkEndLadPipeArmZLay3;
  Int_t    nCards = 3;

  if (iLay==4) {
    endLadPipeUlength = fgkEndLadPipeUlengthLay4;
    endLadPipeArmZ = fgkEndLadPipeArmZLay4;
    nCards = 4;
  }

  AliITSv11GeomCableRound endLadderPipe("endLadderPipe", fgkEndLadPipeOuterDiam/2);
  endLadderPipe.SetNLayers(2); 
  endLadderPipe.SetLayer(0, fgkEndLadPipeInnerDiam/2, coolerMediumSDD, 4);
  endLadderPipe.SetLayer(1, (fgkEndLadPipeOuterDiam-fgkEndLadPipeInnerDiam)/2, phynoxSDD, fColorPhynox);

  Double_t coolUzPos = fgkEndLadPipeOuterDiam/2+2.*fgkmm; //it is the x coord of the axis
  // of the U colling pipe in its center

  Double_t coordA[3] = { fgkEndLadPipeUwidth/2, 0, endLadPipeUlength+coolUzPos};
  Double_t vectA[3]  = {0,0,1};

  Double_t coordB[3] = { fgkEndLadPipeUwidth/2,0, fgkEndLadPipeRadius+coolUzPos};
  Double_t vectB[3]  = {0,0,1};

  Double_t coordC[3] = { fgkEndLadPipeUwidth/2-fgkEndLadPipeRadius, 0, coolUzPos};
  Double_t vectC[3]  = {1,0,0};

  Double_t coordD[3] = {-fgkEndLadPipeUwidth/2+fgkEndLadPipeRadius, 0, coolUzPos};
  Double_t vectD[3]  = {-1,0,0};

  Double_t coordE[3] = {-fgkEndLadPipeUwidth/2, 0, fgkEndLadPipeRadius+coolUzPos};
  Double_t vectE[3]  = {0,0,-1};

  Double_t coordF[3] = {-fgkEndLadPipeUwidth/2,0, endLadPipeUlength+coolUzPos};
  Double_t vectF[3]  = {0,0,-1};

  endLadderPipe.AddCheckPoint( (TGeoVolume *) endLadderCards, 0, coordA, vectA);
  endLadderPipe.AddCheckPoint( (TGeoVolume *) endLadderCards, 1, coordB, vectB);
  endLadderPipe.AddCheckPoint( (TGeoVolume *) endLadderCards, 2, coordC, vectC);
  endLadderPipe.AddCheckPoint( (TGeoVolume *) endLadderCards, 3, coordD, vectD);
  endLadderPipe.AddCheckPoint( (TGeoVolume *) endLadderCards, 4, coordE, vectE);
  endLadderPipe.AddCheckPoint( (TGeoVolume *) endLadderCards, 5, coordF, vectF);

  endLadderPipe.SetInitialNode((TGeoVolume *) endLadderCards); //Set the root node
  //endLadderPipe.CreateAndInsertCableSegment( 1);
  endLadderPipe.CreateAndInsertTubeSegment( 1);
  //endLadderPipe.CreateAndInsertCableSegment( 2);
  endLadderPipe.CreateAndInsertTorusSegment( 2);
  //endLadderPipe.CreateAndInsertCableSegment( 3);
  endLadderPipe.CreateAndInsertTubeSegment( 3);
  //endLadderPipe.CreateAndInsertCableSegment( 4);
  endLadderPipe.CreateAndInsertTorusSegment( 4);
  //endLadderPipe.CreateAndInsertCableSegment( 5);
  endLadderPipe.CreateAndInsertTubeSegment( 5);

  TGeoBBox *endLadPipeArmBox = new TGeoBBox("endLadPipeArmBox",fgkEndLadPipeArmX/2,
				         fgkEndLadPipeArmY/2, endLadPipeArmZ/2);
  TGeoTube *endLadPipeArmTube = new TGeoTube("endLadPipeArmTube", 0,
				    fgkEndLadPipeOuterDiam/2, endLadPipeArmZ/2);

  TGeoTranslation *endLadPipeArmBoxDY1 = new TGeoTranslation("endLadPipeArmBoxDY1",
							    - fgkEndLadPipeArmBoxDX,
							     fgkEndLadPipeArmBoxDY,0);
  TGeoTranslation *endLadPipeArmBoxDY2 = new TGeoTranslation("endLadPipeArmBoxDY2",
							    fgkEndLadPipeArmBoxDX,
							    fgkEndLadPipeArmBoxDY,0);
  endLadPipeArmBoxDY1->RegisterYourself();
  endLadPipeArmBoxDY2->RegisterYourself();

  if(GetDebug(3)) { // Remove compiler warning.
    endLadPipeArmBox->InspectShape();
    endLadPipeArmTube->InspectShape();
  }

  TGeoCompositeShape *endLadPipeArm1 = new TGeoCompositeShape("ITSsddEndLadPipeArm1",
					   "endLadPipeArmBox:endLadPipeArmBoxDY1"
					   "- endLadPipeArmTube");
  TGeoCompositeShape *endLadPipeArm2 = new TGeoCompositeShape("ITSsddEndLadPipeArm2",
					   "endLadPipeArmBox:endLadPipeArmBoxDY2"
					   "- endLadPipeArmTube");

  TGeoVolume *vEndLadPipeArm1 = new TGeoVolume("ITSsddVolEndLadPipeArm1",
					       endLadPipeArm1, alCu12SDD);
  TGeoVolume *vEndLadPipeArm2 = new TGeoVolume("ITSsddVolEndLadPipeArm2",
					       endLadPipeArm2, alCu12SDD);
  vEndLadPipeArm1->SetLineColor(2);
  vEndLadPipeArm2->SetLineColor(2);

  Double_t armZ = (coolUzPos-fgkEndLadPipeOuterDiam/2+endLadPipeArmZ/2
		   +fgkEndLadPipeArmZpos);

  TGeoTranslation *trEndLadPipeArm1 = new TGeoTranslation("trEndLadPipeArm1",
					  -fgkEndLadPipeUwidth/2,0,armZ);
  TGeoTranslation *trEndLadPipeArm2 = new TGeoTranslation("trEndLadPipeArm2",
					   fgkEndLadPipeUwidth/2,0,armZ);

  endLadderCards->AddNode(vEndLadPipeArm1, 1, trEndLadPipeArm1);
  endLadderCards->AddNode(vEndLadPipeArm2, 1, trEndLadPipeArm2);

  //=*********************************
  //--- LV cards
  TGeoVolumeAssembly *cardLVassemblyR = fCardLVR;
  TGeoVolumeAssembly *cardLVassemblyL = fCardLVL;

  Double_t spaceBetweenCards = 0.2*fgkmm; 

  Double_t cardLVxShift = (fgkEndLadPipeUwidth/2-fgkEndLadPipeArmX/2
			   +fgkEndLadPipeArmBoxDX);
  Double_t cardLVyShift = (-fgkLVcoolPosY-fgkLVcoolY1+fgkLVcoolY3
			   +fgkEndLadPipeArmY/2+fgkEndLadPipeArmBoxDY);

  Double_t alLVcoolZ3 = (fgkLVcardCuZ+fgkLVcardZ+2.*(fgkLVChip0SiZ+fgkLVChip0Z)
			  +fgkLVcoolZ1*2.);

  Double_t firstLVCardZ = fgkEndLadPipeArmZpos-fgkEndLadPipeOuterDiam/2.+alLVcoolZ3/2
                          +coolUzPos+1.25*fgkmm;
  // Position in z of the first LVB with respect to the start of the cooling
  // rectangular arm, coming  (from inside of the ladder)
  // The cards are added one after the other

  for (Int_t iCard=0; iCard<nCards; iCard++) {

    Double_t cardLVzShift = firstLVCardZ + 
      Double_t(iCard)*(alLVcoolZ3 + 2.*spaceBetweenCards+fgkHVCardCool3Z);

    TGeoTranslation *trCardLVassemblyR = new TGeoTranslation(cardLVxShift,
					     cardLVyShift, cardLVzShift);
    TGeoTranslation *trCardLVassemblyL = new TGeoTranslation(-cardLVxShift,
					     cardLVyShift, cardLVzShift);

    endLadderCards->AddNode(cardLVassemblyR, iCard+1, trCardLVassemblyR);
    endLadderCards->AddNode(cardLVassemblyL, iCard+1, trCardLVassemblyL);
  }

  //=*********************************
  //--- HV cards
  TGeoVolumeAssembly *cardHV = fCardHV;

  Double_t coolHVdy = (fgkHVCardCoolDY + fgkHVCardCool3Y 
		       + fgkEndLadPipeArmY/2 + fgkEndLadPipeArmBoxDY);

  Double_t coolHVCenterShift = (fgkHVCardCool3Z/2-fgkHVCardCool2Z
				-(fgkHVCardCeramZ)/2); 

  for (Int_t iCard=0; iCard<nCards; iCard++) {

    Double_t fact = iCard*2.+1.;
    Double_t coolHVdz = (firstLVCardZ + alLVcoolZ3*fact/2 + spaceBetweenCards*fact
			 + fgkHVCardCool3Z*fact/2. + coolHVCenterShift);
    TGeoTranslation *trCardHV = new TGeoTranslation(0,coolHVdy, coolHVdz);
    endLadderCards->AddNode(cardHV, iCard+1, trCardHV);
  }

  //=*********************************
  //--- Carlos card

  TGeoVolumeAssembly *assemblySupCarlos = fCardCarlos;
//   TGeoRotation *carlosSupRot1 = new TGeoRotation("carlosSuppAngle",
// 						 0, -fgkCarlosSuppAngle, 0);

  Double_t spaceBetweenCarlsoCards = 0.1*fgkmm;
  Double_t firstCarlosCardZ = (firstLVCardZ - alLVcoolZ3/2 + alLVcoolZ3*4 +
			       fgkHVCardCool3Z*4 + spaceBetweenCards*7 + 2*fgkmm);
  // position in z of the first Carlos board, coming  from inside of the ladder

  Double_t coolCarlosDy = (fgkCarlosSuppY3/2 + fgkEndLadPipeArmY/2 + 
			   fgkEndLadPipeArmBoxDY);

  for (Int_t iCard=0; iCard<nCards; iCard++) {

    Double_t carloszPos = ( firstCarlosCardZ + fgkCarlosSuppZ3/2 +
			    iCard*(fgkCarlosSuppZ3+spaceBetweenCarlsoCards) );
    TGeoCombiTrans *carlosPos = new TGeoCombiTrans(0,coolCarlosDy,carloszPos,
						   (TGeoRotation*) fCommonTr[0]);

    endLadderCards->AddNode(assemblySupCarlos, iCard, carlosPos);
  }

  return endLadderCards;
}


//________________________________________________________________________
TGeoVolume*  AliITSv11GeometrySDD::CreateEndLadderCardsV(Int_t iLay) {
// 
// return an Pcon containing the LV, HV and Carlos cards of one ladder
// and their cooling system 
// This is the code actually used for the end ladder cards
//

  TGeoMedium *alCu12SDD       = GetMedium("AL$"); // ITSsddAlCu12 : to code !!!!!!!!!!!!!!
  TGeoMedium *phynoxSDD       = GetMedium("INOX$");
  TGeoMedium *coolerMediumSDD = GetMedium("WATER$");
  TGeoMedium *copper          = GetMedium("COPPER$");
  TGeoMedium *plastic         = GetMedium("SDDKAPTON (POLYCH2)$");  // ???
  TGeoMedium *airSDD          = GetMedium("SDD AIR$");
  TGeoMedium *opticalFiber    = GetMedium("SDD OPTICFIB$");
  TGeoMedium *polyurethane    = GetMedium("POLYURETHANE$");

  Double_t endLadPipeUlength = fgkEndLadPipeUlengthLay3;
  Double_t endLadPipeArmZ    = fgkEndLadPipeArmZLay3;
  Int_t    nCards = 3;
  Double_t rREF   = fgkEndLaddCardsShortRadiusLay3;
  Double_t deltaZcables = 0;
  // reference radius corresponding to local y=0

  if (iLay==4) {
    endLadPipeUlength = fgkEndLadPipeUlengthLay4;
    endLadPipeArmZ = fgkEndLadPipeArmZLay4;
    nCards = 4;
    rREF = fgkEndLaddCardsShortRadiusLay4;
    deltaZcables = 2.8*fgkmm;
  }

  Double_t cardLVxShift = (fgkEndLadPipeUwidth/2-fgkEndLadPipeArmX/2
			   +fgkEndLadPipeArmBoxDX);
  Double_t cardLVyShift = (-fgkLVcoolPosY-fgkLVcoolY1+fgkLVcoolY3
			   +fgkEndLadPipeArmY/2+fgkEndLadPipeArmBoxDY);

  Double_t rMin   = rREF + cardLVyShift;
  // (The LV card is defining rMin because it is the lower object)

  Double_t thickTotCable = 0.5;

  //==================================
  //--- The Pcon container

  // minimum angle of the Pcon :
  Double_t tanDPhi = ((fgkEndLadPipeUwidth/2+fgkEndLadPipeArmX/2) /
		     (rREF-fgkEndLadPipeArmY/2) );
  Double_t dphi = 2*TMath::ATan(tanDPhi)*TMath::RadToDeg();
  Double_t phi0 = 90-dphi/2;
  Double_t coolUzPos = fgkEndLadPipeOuterDiam/2 + fgkDistEndLaddCardsLadd; // it is the z coord of the axis
                                                        // of the U colling pipe in its center
  Double_t zMax = endLadPipeUlength+coolUzPos;
  Double_t rMax = rMin + fgkLVcardY;
  rMax = TMath::Sqrt(rMax*rMax + cardLVxShift*cardLVxShift);
  Double_t cablesRadius  = rMax-0.5;

  TGeoPcon *containerShape = new TGeoPcon("EndLadderCcontainerShape", phi0, dphi, 10);
   //DefineSection(Int_t snum, Double_t z, Double_t rmin, Double_t rmax);
  // hard coded numbers are fine tuning to avoid overlaps with other volume in the old geometry
  containerShape->DefineSection(0, fgkDistEndLaddCardsLadd, rREF-fgkEndLadPipeOuterDiam/2-0.2, rMax);
  containerShape->DefineSection(1, fgkDistEndLaddCardsLadd+1.4, rREF-fgkEndLadPipeOuterDiam/2-0.2, rMax);
  containerShape->DefineSection(2,  fgkDistEndLaddCardsLadd+1.4, rMin, rMax);
  containerShape->DefineSection(3,  endLadPipeArmZ+2*fgkEndLadPipeRadius, rMin, rMax);
  containerShape->DefineSection(4,  endLadPipeArmZ+2*fgkEndLadPipeRadius, rREF-1.*fgkmm, rMax);
  containerShape->DefineSection(5,  zMax, rREF-1.*fgkmm, rMax);
  // the following is quite dirty but works for the moment ...
  containerShape->DefineSection(6,  zMax,  rREF+fgkCarlosCardZ1/2, rMax);
  containerShape->DefineSection(7, zMax+1, cablesRadius-thickTotCable/2, rMax);

  // The next parameters define the shape of the Pcon at its end and where cables
  // are escaping...
  Double_t cableSectionR1 = cablesRadius-thickTotCable/2;
  Double_t cableSectionR2 = rMax;
  Double_t cableSectionZ1 = zMax + 23.6*fgkmm + 3.0*fgkcm + deltaZcables;
  Double_t cableSectionZ2 = zMax + 23.6*fgkmm + 4.0*fgkcm + deltaZcables;
  // Those numbers are to be fixed to stick the maximum to the SDD cone
  // (hardcoded numbers are ugly, but it's easier to find where to stop)

  containerShape->DefineSection(8, cableSectionZ1, cableSectionR1, rMax);
  containerShape->DefineSection(9, cableSectionZ2, cableSectionR2, rMax);

  TGeoVolume *endLadderCards = new TGeoVolume("endLadderCards",containerShape,airSDD);
  //endLadderCards->SetVisibility(kFALSE);

  //=*********************************
  //--- The rounded pipe for the end ladder card cooling

  AliITSv11GeomCableRound endLadderPipe("endLadderPipe", fgkEndLadPipeOuterDiam/2);
  endLadderPipe.SetNLayers(2); 
  endLadderPipe.SetLayer(0, fgkEndLadPipeInnerDiam/2, coolerMediumSDD, 4);
  endLadderPipe.SetLayer(1, (fgkEndLadPipeOuterDiam-fgkEndLadPipeInnerDiam)/2, phynoxSDD, fColorPhynox);

  Double_t coordA[3] = { fgkEndLadPipeUwidth/2, rREF, endLadPipeUlength+coolUzPos};
  Double_t vectA[3]  = {0,0,1};

  Double_t coordB[3] = { fgkEndLadPipeUwidth/2,rREF, fgkEndLadPipeRadius+coolUzPos};
  Double_t vectB[3]  = {0,0,1};

  Double_t coordC[3] = { fgkEndLadPipeUwidth/2-fgkEndLadPipeRadius, rREF, coolUzPos};
  Double_t vectC[3]  = {1,0,0};

  Double_t coordD[3] = {-fgkEndLadPipeUwidth/2+fgkEndLadPipeRadius, rREF, coolUzPos};
  Double_t vectD[3]  = {-1,0,0};

  Double_t coordE[3] = {-fgkEndLadPipeUwidth/2, rREF, fgkEndLadPipeRadius+coolUzPos};
  Double_t vectE[3]  = {0,0,-1};

  Double_t coordF[3] = {-fgkEndLadPipeUwidth/2,rREF, endLadPipeUlength+coolUzPos};
  Double_t vectF[3]  = {0,0,-1};

  endLadderPipe.AddCheckPoint( (TGeoVolume *) endLadderCards, 0, coordA, vectA);
  endLadderPipe.AddCheckPoint( (TGeoVolume *) endLadderCards, 1, coordB, vectB);
  endLadderPipe.AddCheckPoint( (TGeoVolume *) endLadderCards, 2, coordC, vectC);
  endLadderPipe.AddCheckPoint( (TGeoVolume *) endLadderCards, 3, coordD, vectD);
  endLadderPipe.AddCheckPoint( (TGeoVolume *) endLadderCards, 4, coordE, vectE);
  endLadderPipe.AddCheckPoint( (TGeoVolume *) endLadderCards, 5, coordF, vectF);

  endLadderPipe.SetInitialNode((TGeoVolume *) endLadderCards); //Set the root node
  //endLadderPipe.CreateAndInsertCableSegment( 1);
  endLadderPipe.CreateAndInsertTubeSegment( 1);
  //endLadderPipe.CreateAndInsertCableSegment( 2);
  endLadderPipe.CreateAndInsertTorusSegment( 2);
  //endLadderPipe.CreateAndInsertCableSegment( 3);
  endLadderPipe.CreateAndInsertTubeSegment( 3);
  //endLadderPipe.CreateAndInsertCableSegment( 4);
  endLadderPipe.CreateAndInsertTorusSegment( 4);
  //endLadderPipe.CreateAndInsertCableSegment( 5);
  endLadderPipe.CreateAndInsertTubeSegment( 5);

  TGeoBBox *endLadPipeArmBox = new TGeoBBox("endLadPipeArmBox",fgkEndLadPipeArmX/2,
				         fgkEndLadPipeArmY/2, endLadPipeArmZ/2);
  TGeoTube *endLadPipeArmTube = new TGeoTube("endLadPipeArmTube", 0,
				    fgkEndLadPipeOuterDiam/2, endLadPipeArmZ/2);

  TGeoTranslation *endLadPipeArmBoxDY1 = new TGeoTranslation("endLadPipeArmBoxDY1",
							    - fgkEndLadPipeArmBoxDX,
							     fgkEndLadPipeArmBoxDY,0);
  TGeoTranslation *endLadPipeArmBoxDY2 = new TGeoTranslation("endLadPipeArmBoxDY2",
							    fgkEndLadPipeArmBoxDX,
							    fgkEndLadPipeArmBoxDY,0);
  endLadPipeArmBoxDY1->RegisterYourself();
  endLadPipeArmBoxDY2->RegisterYourself();

  if(GetDebug(3)) { // Remove compiler warning.
    endLadPipeArmBox->InspectShape();
    endLadPipeArmTube->InspectShape();
  }

  TGeoCompositeShape *endLadPipeArm1 = new TGeoCompositeShape("ITSsddEndLadPipeArm1",
					   "endLadPipeArmBox:endLadPipeArmBoxDY1"
					   "- endLadPipeArmTube");
  TGeoCompositeShape *endLadPipeArm2 = new TGeoCompositeShape("ITSsddEndLadPipeArm2",
					   "endLadPipeArmBox:endLadPipeArmBoxDY2"
					   "- endLadPipeArmTube");

  TGeoVolume *vEndLadPipeArm1 = new TGeoVolume("ITSsddVolEndLadPipeArm1",
					       endLadPipeArm1, alCu12SDD);
  TGeoVolume *vEndLadPipeArm2 = new TGeoVolume("ITSsddVolEndLadPipeArm2",
					       endLadPipeArm2, alCu12SDD);
  vEndLadPipeArm1->SetLineColor(2);
  vEndLadPipeArm2->SetLineColor(2);

  Double_t armZ = (coolUzPos-fgkEndLadPipeOuterDiam/2+endLadPipeArmZ/2
		   +fgkEndLadPipeArmZpos);

  TGeoTranslation *trEndLadPipeArm1 = new TGeoTranslation("trEndLadPipeArm1",
					  -fgkEndLadPipeUwidth/2,rREF,armZ);
  TGeoTranslation *trEndLadPipeArm2 = new TGeoTranslation("trEndLadPipeArm2",
					   fgkEndLadPipeUwidth/2,rREF,armZ);

  endLadderCards->AddNode(vEndLadPipeArm1, 1, trEndLadPipeArm1);
  endLadderCards->AddNode(vEndLadPipeArm2, 1, trEndLadPipeArm2);

  //=*********************************
  //--- LV cards
  TGeoVolumeAssembly *cardLVassemblyR = fCardLVR;
  TGeoVolumeAssembly *cardLVassemblyL = fCardLVL;

  Double_t spaceBetweenCards = 0.2*fgkmm; 


  Double_t alLVcoolZ3 = (fgkLVcardCuZ+fgkLVcardZ+2.*(fgkLVChip0SiZ+fgkLVChip0Z)
			  +fgkLVcoolZ1*2.);

  Double_t firstLVCardZ = fgkEndLadPipeArmZpos-fgkEndLadPipeOuterDiam/2.+alLVcoolZ3/2
                          +coolUzPos+1.25*fgkmm;
  // Position in z of the first LVB with respect to the start of the cooling
  // rectangular arm, coming  (from inside of the ladder)
  // The cards are added one after the other

  for (Int_t iCard=0; iCard<nCards; iCard++) {

    Double_t cardLVzShift = firstLVCardZ + 
      Double_t(iCard)*(alLVcoolZ3 + 2.*spaceBetweenCards+fgkHVCardCool3Z);

    TGeoTranslation *trCardLVassemblyR = new TGeoTranslation(cardLVxShift,
					     cardLVyShift+rREF, cardLVzShift);
    TGeoTranslation *trCardLVassemblyL = new TGeoTranslation(-cardLVxShift,
					     cardLVyShift+rREF, cardLVzShift);

    endLadderCards->AddNode(cardLVassemblyR, iCard+1, trCardLVassemblyR);
    endLadderCards->AddNode(cardLVassemblyL, iCard+1, trCardLVassemblyL);
  }

  //=*********************************
  //--- HV cards
  TGeoVolumeAssembly *cardHV = fCardHV;

  Double_t coolHVdy = (fgkHVCardCoolDY + fgkHVCardCool3Y
		       + fgkEndLadPipeArmY/2 + fgkEndLadPipeArmBoxDY);
  // shift of the HV card in local y w.r.t the local y=0 (center of cooling tube)

  Double_t coolHVCenterShift = (fgkHVCardCool3Z/2-fgkHVCardCool2Z
				-(fgkHVCardCeramZ)/2); 

  for (Int_t iCard=0; iCard<nCards; iCard++) {

    Double_t fact = iCard*2.+1.;
    Double_t coolHVdz = (firstLVCardZ + alLVcoolZ3*fact/2 + spaceBetweenCards*fact
			 + fgkHVCardCool3Z*fact/2. + coolHVCenterShift);
    TGeoTranslation *trCardHV = new TGeoTranslation(0,coolHVdy+rREF, coolHVdz);
    endLadderCards->AddNode(cardHV, iCard+1, trCardHV);
  }

  //=*********************************
  //--- Carlos card

  TGeoVolumeAssembly *assemblySupCarlos = fCardCarlos;
//   TGeoRotation *carlosSupRot1 = new TGeoRotation("carlosSuppAngle",
// 						 0, -fgkCarlosSuppAngle, 0);

  Double_t spaceBetweenCarlsoCards = 0.1*fgkmm;
  Double_t firstCarlosCardZ = (firstLVCardZ - alLVcoolZ3/2 + alLVcoolZ3*4 +
			       fgkHVCardCool3Z*4 + spaceBetweenCards*7 + 2*fgkmm);
  // position in z of the first Carlos board, coming  from inside of the ladder

  Double_t coolCarlosDy = (fgkCarlosSuppY3/2 + fgkEndLadPipeArmY/2 + 
			   fgkEndLadPipeArmBoxDY);

  for (Int_t iCard=0; iCard<nCards; iCard++) {

    Double_t carloszPos = ( firstCarlosCardZ + fgkCarlosSuppZ3/2 +
			    iCard*(fgkCarlosSuppZ3+spaceBetweenCarlsoCards) );
    TGeoCombiTrans *carlosPos = new TGeoCombiTrans(0,coolCarlosDy+rREF,carloszPos,
						   (TGeoRotation*) fCommonTr[0]);

    endLadderCards->AddNode(assemblySupCarlos, iCard, carlosPos);
  }


  //=*********************************
  //--- Cables


  Double_t sectionV   = (fgkSectionCuPerMod+fgkSectionPlastPerMod
			 + fgkSectionGlassPerMod)*nCards;
  // We fix thickness, then width is calculated accordingly
  Double_t width      = sectionV/thickTotCable;
  Double_t thickCu    = thickTotCable*fgkSectionCuPerMod
            / (fgkSectionCuPerMod+fgkSectionPlastPerMod+fgkSectionGlassPerMod);
  Double_t thickPlast = thickTotCable*fgkSectionPlastPerMod
            / (fgkSectionCuPerMod+fgkSectionPlastPerMod+fgkSectionGlassPerMod);
  Double_t thickGlass = thickTotCable - thickCu - thickPlast;

  Double_t thickCoolPolyu = fgkSectionCoolPolyuEL/width;
  Double_t thickCoolWater = fgkSectionCoolWaterEL/width;

  Double_t totalThick = thickTotCable + thickCoolPolyu + thickCoolWater;

  AliITSv11GeomCableFlat cable("SDDcableEndLadder",width,totalThick);
  cable.SetNLayers(5);
  cable.SetLayer(0, thickCu, copper, kRed);
  cable.SetLayer(1, thickPlast, plastic, kYellow);
  cable.SetLayer(2, thickGlass, opticalFiber, kGreen);
  cable.SetLayer(3, thickCoolPolyu, polyurethane, kGray);
  cable.SetLayer(4, thickCoolWater, coolerMediumSDD, kBlue);

  Double_t zVect[3]={0,0,1};
  Double_t xMinCable = firstCarlosCardZ+nCards*(fgkCarlosSuppZ3
		       +spaceBetweenCarlsoCards)/2 + 2.9;
  // the 2.9cm is for taking into account carlos card angle...

  Double_t zEndCable = GetConeZ(cablesRadius-thickTotCable/2, cableSectionR1,
				cableSectionR2,cableSectionZ1,cableSectionZ2);

  Double_t pos1[3] = {0, cablesRadius, xMinCable};
  Double_t pos2[3] = {0, cablesRadius, zEndCable};
  cable.AddCheckPoint( endLadderCards, 0, pos1, zVect );
  cable.AddCheckPoint( endLadderCards, 1, pos2, zVect );
  cable.SetInitialNode(endLadderCards);
  cable.CreateAndInsertCableSegment(1);

  return endLadderCards;
}

//________________________________________________________________________
TGeoVolumeAssembly* AliITSv11GeometrySDD::CreateSupportRing(Int_t iLay) {
//
// return an assembly of the support rings, attaching the ladders to the cone 
//


  iLay = iLay;

  TGeoMedium *stainless = GetMedium("INOX$"); // To code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoVolumeAssembly *supportRing = new TGeoVolumeAssembly("supportRing");
 

  //**********************************
  // ruby cage

  Double_t fgkRubyCageX          = 9*fgkmm;
  Double_t fgkRubyCageY          = 5.5*fgkmm;
  Double_t fgkRubyCageZ          = 8*fgkmm;
  Double_t fgkRubyCageInternSide = 5.*fgkmm; //side of the internal square
  Double_t fgkRubyCageHoleDX     = 2.*fgkmm;
  Double_t fgkRubyCageVIntern    = 5.42*fgkmm;
  Double_t fgkRubyCageScrewHoleR = 4.5/2*fgkmm;
  Double_t fgkRubyCageScrewHoleY = 1.5*fgkmm;

  TGeoBBox *rubyCageBox = new TGeoBBox("rubyCageBox",fgkRubyCageX/2,fgkRubyCageY/2,
				       fgkRubyCageZ/2);

  Double_t epsilon = 1e-10; //dummy epsilon to force the gl viewer to show holes

  // pieces common to both square and V cages
  TGeoBBox *rubyCageInternBox = new TGeoBBox("rubyCageInternBox",fgkRubyCageInternSide/2,
				    fgkRubyCageY/2+epsilon, fgkRubyCageInternSide/2);

  TGeoTube *screwHole = new TGeoTube("screwHole", 0, fgkRubyCageScrewHoleR,
				     fgkRubyCageHoleDX/2+epsilon);

  TGeoRotation *rotV = new TGeoRotation("", 90,90,-90);
  TGeoCombiTrans *trScrewHole = new TGeoCombiTrans("trScrewHole",
				    fgkRubyCageX/2-fgkRubyCageHoleDX/2,
				   -fgkRubyCageY/2+fgkRubyCageScrewHoleY,0,rotV);
  trScrewHole->RegisterYourself();

  TGeoBBox *screwHoleFoot = new TGeoBBox("screwHoleFoot",fgkRubyCageHoleDX/2+epsilon,
				fgkRubyCageScrewHoleY/2+epsilon, fgkRubyCageScrewHoleR);
  TGeoTranslation *trScrewHoleFoot = new TGeoTranslation("trScrewHoleFoot",
					 fgkRubyCageX/2-fgkRubyCageHoleDX/2,
					-fgkRubyCageY/2+fgkRubyCageScrewHoleY/2, 0);
  trScrewHoleFoot->RegisterYourself();


  // pieces which differ
  Double_t rubyCageVInternBoxX = fgkRubyCageVIntern - fgkRubyCageInternSide/2;

  TGeoBBox *rubyCageVInternBox = new TGeoBBox("rubyCageVInternBox",rubyCageVInternBoxX/2,
				     fgkRubyCageY/2+epsilon, fgkRubyCageInternSide/2);

  TGeoTranslation *trRubyCageVInternBox = new TGeoTranslation("trRubyCageVInternB",
			 fgkRubyCageX/2-fgkRubyCageHoleDX-rubyCageVInternBoxX/2,0,0);
  trRubyCageVInternBox->RegisterYourself();

  TGeoTrd1 *rubyCageVInternTriangl = new TGeoTrd1("rubyCageVInternTriangl", 0,
				     fgkRubyCageInternSide/2, fgkRubyCageY/2+epsilon,
						  fgkRubyCageInternSide/4);

  TGeoCombiTrans *trRubyCageVInternTriangl = new TGeoCombiTrans("trRubyCageVInternTriangl",
	fgkRubyCageX/2-fgkRubyCageHoleDX-rubyCageVInternBoxX-fgkRubyCageInternSide/4
								+epsilon,0,0, rotV );
  trRubyCageVInternTriangl->RegisterYourself();

  //---
  TGeoCompositeShape *rubyCageSquare = new TGeoCompositeShape("rubyCageSquare",
				           "rubyCageBox-(rubyCageInternBox"
				 "+screwHole:trScrewHole+screwHoleFoot:trScrewHoleFoot)");

  TGeoVolume *vRubyCageSquare = new TGeoVolume("vRubyCageSquare",
					       rubyCageSquare, stainless);
  vRubyCageSquare->SetLineColor(10);
  
  TGeoCompositeShape *rubyCageV = new TGeoCompositeShape("rubyCageV",
			      	 "rubyCageBox-(rubyCageVInternBox:trRubyCageVInternB"
			       	 "+rubyCageVInternTriangl:trRubyCageVInternTriangl"
				 "+screwHole:trScrewHole+screwHoleFoot:trScrewHoleFoot)");
  TGeoVolume *vRubyCageV = new TGeoVolume("vRubyCageV", rubyCageV, stainless);
  vRubyCageV->SetLineColor(10);

  if(GetDebug(3)) { // Remove compiler warning.
    rubyCageBox->InspectShape();
    rubyCageInternBox->InspectShape();
    screwHole->InspectShape();
    screwHoleFoot->InspectShape();
    rubyCageVInternBox->InspectShape();
    rubyCageVInternTriangl->InspectShape();
  }

  supportRing->AddNode(vRubyCageSquare, 0, 0);
  //supportRing->AddNode(vRubyCageV, 0, 0);
  return supportRing;
}



//________________________________________________________________________
void AliITSv11GeometrySDD::CreateSDDsensor() {
//
// return a box containing the SDD sensor
//

  TGeoMedium *airSDD         = GetMedium("SDD AIR$");
  TGeoMedium *siliconSDD     = GetMedium("SDD SI insensitive$");  // ITSsddSi
  TGeoMedium *siliconSDDsens = GetMedium("SI$");                  // ITSsddSi
  TGeoMedium *alSDD          = GetMedium("AL$");                  // ITSal
  TGeoMedium *polyhamideSDD  = GetMedium("SDDKAPTON (POLYCH2)$"); // ITSsddKAPTON_POLYCH2
  TGeoMedium *glassSDD       = GetMedium("SDD SI insensitive$");  //  To code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Double_t rWraping = fgkWaferThickness/2+fgkWaHVcableAlThick+fgkWaHVcablePolyThick;
  Double_t witdhCableBox = (fgkWaHVcableWitdh - TMath::Pi()*rWraping)/2;
  // width : in the beam direction !

  Double_t sensoxBoxLength = ( fgkWaferLength +
			       2*(rWraping+witdhCableBox-fgkWaHVcableDW) );
  // Makes life easier to include the space for the WA HV cable on both sides 
  Double_t sensoxBoxThick = fgkWaferThickness +
                            2*(fgkWaHVcableAlThick+fgkWaHVcablePolyThick);

//   cout << "fgkWaferLength=" << fgkWaferLength << " sensoxBoxLength="<< sensoxBoxLength <<endl;
//   cout << "fgkWaferThickness=" << fgkWaferThickness << " sensoxBoxThick=" << sensoxBoxThick << endl;

  TGeoBBox *box = new TGeoBBox("ITSsddSensorBox",
		      fgkWaferWidth/2, sensoxBoxThick/2, sensoxBoxLength/2);

  fSDDsensor3 = new TGeoVolume("ITSsddSensor3", box, airSDD);
  fSDDsensor4 = new TGeoVolume("ITSsddSensor4", box, airSDD);


  //****************************
  // silicon wafer
  //****************************
  if (fAddSensors) {
    // we need 2 different sensor objects, because they have to have different names
    // This is required for the step manager

    TGeoBBox *waferShape = new TGeoBBox("ITSsddWaferShape",
			       fgkWaferWidth/2, fgkWaferThickness/2, fgkWaferLength/2);


    TGeoVolume *wafer3 = new TGeoVolume("ITSsddWafer3", waferShape, siliconSDD);
    wafer3->SetLineColor(fColorSilicon);
    TGeoBBox *sensBox3 = new TGeoBBox("ITSsddSensorSensBox3",
		        fgkWaferWidthSens/2, fgkWaferThickSens/2, fgkWaferLengthSens/2);
    TGeoVolume *sensVol3 = new TGeoVolume(fgSDDsensitiveVolName3,sensBox3, siliconSDDsens);
    sensVol3->SetLineColor(fColorSilicon+5);
    wafer3->AddNode(sensVol3, 1, 0);
    fSDDsensor3->AddNode(wafer3, 1, 0);

    TGeoVolume *wafer4 = new TGeoVolume("ITSsddWafer4", waferShape, siliconSDD);
    wafer4->SetLineColor(fColorSilicon);
    TGeoBBox *sensBox4 = new TGeoBBox("ITSsddSensorSensBox4",
		        fgkWaferWidthSens/2, fgkWaferThickSens/2, fgkWaferLengthSens/2);
    TGeoVolume *sensVol4 = new TGeoVolume(fgSDDsensitiveVolName4,sensBox4, siliconSDDsens);
    sensVol4->SetLineColor(fColorSilicon+5);
    wafer4->AddNode(sensVol4, 1, 0);
    fSDDsensor4->AddNode(wafer4, 1, 0);
  };
  
  //****************************
  // glass
  //****************************
  TGeoBBox *glass = new TGeoBBox("ITSsddGlassBox", fgkSensorGlassLX/2,
				 fgkSensorGlassLY/2, fgkSensorGlassLZ/2);
  TGeoVolume *vGlass  = new  TGeoVolume("ITSsddGlass",glass, glassSDD);
  vGlass->SetLineColor(fColorGlass);
  TGeoTranslation *glassTr1 = new TGeoTranslation("",fgkGlassDXOnSensor,
				  fgkWaferThickness/2+fgkSensorGlassLY/2,
				  fgkGlassDZOnSensor);
  TGeoTranslation *glassTr2 = new TGeoTranslation("",-fgkGlassDXOnSensor,
				  fgkWaferThickness/2+fgkSensorGlassLY/2,
				  fgkGlassDZOnSensor);
  TGeoTranslation *glassTr3 = new TGeoTranslation("",fgkGlassDXOnSensor,
				  fgkWaferThickness/2+fgkSensorGlassLY/2,
				  -fgkGlassDZOnSensor);
  TGeoTranslation *glassTr4 = new TGeoTranslation("",-fgkGlassDXOnSensor,
				  fgkWaferThickness/2+fgkSensorGlassLY/2,
				  -fgkGlassDZOnSensor);
  fSDDsensor3->AddNode(vGlass, 1, glassTr1);
  fSDDsensor3->AddNode(vGlass, 2, glassTr2);
  fSDDsensor3->AddNode(vGlass, 3, glassTr3);
  fSDDsensor3->AddNode(vGlass, 4, glassTr4);

  fSDDsensor4->AddNode(vGlass, 1, glassTr1);
  fSDDsensor4->AddNode(vGlass, 2, glassTr2);
  fSDDsensor4->AddNode(vGlass, 3, glassTr3);
  fSDDsensor4->AddNode(vGlass, 4, glassTr4);

  //****************************
  // Wrap-around cable
  //****************************
  if (fAddHVcables) {
  AliITSv11GeomCableFlat waHVCable("ITSsddWaHVCableU",witdhCableBox,
				      fgkWaHVcableAlThick+fgkWaHVcablePolyThick);
  waHVCable.SetNLayers(2);
  waHVCable.SetLayer(0, fgkWaHVcablePolyThick,polyhamideSDD,fColorPolyhamide);
  waHVCable.SetLayer(1, fgkWaHVcableAlThick, alSDD, fColorAl);
  waHVCable.SetInitialNode(fSDDsensor3);

  Double_t x1[3], x2[3], vX[3] = {1,0,0};
  x1[0] = -fgkWaHVcableLength/2;
  x2[0] = -x1[0];
  x1[1] = (fgkWaferThickness + waHVCable.GetThickness())/2;
  x2[1] = x1[1];
  x1[2] = fgkWaferLength/2+waHVCable.GetWidth()/2-fgkWaHVcableDW;
  x2[2] = x1[2];

  waHVCable.AddCheckPoint(fSDDsensor3, 0, x1, vX);
  waHVCable.AddCheckPoint(fSDDsensor3, 1, x2, vX);
  TGeoCombiTrans *ctSegment = 0;
  TGeoVolume* segment = waHVCable.CreateAndInsertBoxCableSegment(1,-90, &ctSegment);
  fSDDsensor4->AddNode(segment, 1, ctSegment);

  x1[1] = -x1[1];
  x2[1] = x1[1];
  waHVCable.SetName("ITSsddWaHVCableD");
  waHVCable.ResetPoints();
  waHVCable.AddCheckPoint(fSDDsensor3, 0, x1, vX);
  waHVCable.AddCheckPoint(fSDDsensor3, 1, x2, vX);
  segment = waHVCable.CreateAndInsertBoxCableSegment(1, 90, &ctSegment);
  fSDDsensor4->AddNode(segment, 1, ctSegment);

  AliITSv11GeomCableRound waHVCableFold("ITSsddWaHVCableFold",
					   rWraping);
  waHVCableFold.SetPhi(180,360);
  waHVCableFold.SetNLayers(2);
  waHVCableFold.SetLayer(0, fgkWaferThickness/2+fgkWaHVcablePolyThick,
			 polyhamideSDD, fColorPolyhamide);
  waHVCableFold.SetLayer(1, fgkWaHVcableAlThick, alSDD, fColorAl);
  waHVCableFold.SetInitialNode(fSDDsensor3);
  x1[1] = 0;
  x2[1] = 0;
  x1[2] = fgkWaferLength/2-fgkWaHVcableDW+witdhCableBox;
  x2[2] = x1[2];
  waHVCableFold.AddCheckPoint(fSDDsensor3, 0, x1, vX);
  waHVCableFold.AddCheckPoint(fSDDsensor3, 1, x2, vX);
  segment = waHVCableFold.CreateAndInsertCableSegment(1, &ctSegment);
  fSDDsensor4->AddNode(segment, 1, ctSegment);


  //****************************
  // transition cable
  //****************************
  Double_t headRadius = (fgkTransitHVHeadLX*fgkTransitHVHeadLX/4.+
			 fgkTransitHVHeadLZ*fgkTransitHVHeadLZ)
                        /(2.*fgkTransitHVHeadLZ);
  Double_t theta = TMath::ATan2(fgkTransitHVHeadLX/2,
				headRadius-fgkTransitHVHeadLZ)
                   *TMath::RadToDeg();

  TGeoTubeSeg *headPoly = new TGeoTubeSeg(0,headRadius,
					  fgkTransitHVPolyThick/2,
					  90-theta,90+theta);
  headPoly->SetName("headPoly");
  TGeoTranslation *headPolyTr = new TGeoTranslation(0,0,
				    -fgkTransitHVPolyThick/2);
  headPolyTr->SetName("headPolyTr");
  headPolyTr->RegisterYourself();

  TGeoTubeSeg *headAl = new TGeoTubeSeg(0,headRadius,
					fgkTransitHVAlThick/2,
					90-theta,90+theta);
  headAl->SetName("headAl");
  TGeoTranslation *headAlTr = new TGeoTranslation(0,0,
				  -fgkTransitHVPolyThick
				  -fgkTransitHVAlThick/2);
  headAlTr->SetName("headAlTr");
  headAlTr->RegisterYourself();

  TGeoBBox *cache = new TGeoBBox(fgkTransitHVHeadLX/2,
				 (headRadius-fgkTransitHVHeadLZ)/2,
		       (fgkTransitHVPolyThick+fgkTransitHVAlThick)/2);
  cache->SetName("cache");

  TGeoTranslation *headCacheTr = new TGeoTranslation(0,
				(headRadius-fgkTransitHVHeadLZ)/2,
				 -(fgkTransitHVPolyThick
				   +fgkTransitHVAlThick)/2);
  headCacheTr->SetName("cacheTr");
  headCacheTr->RegisterYourself();

  TGeoCompositeShape *headPolyComp = new TGeoCompositeShape(
		      "headPoly:headPolyTr-cache:cacheTr");
  TGeoVolume *vHeadPolyComp = new TGeoVolume(
              "ITSsddHVtransitHeadPoly",headPolyComp, polyhamideSDD);
  vHeadPolyComp->SetLineColor(fColorPolyhamide);
  TGeoCompositeShape *headAlComp = new TGeoCompositeShape(
				       "headAl:headAlTr-cache:cacheTr");
  TGeoVolume *vHeadAlComp = new TGeoVolume(
              "ITSsddHVtransitHeadAl",headAlComp, alSDD);
  vHeadAlComp->SetLineColor(fColorAl);


//   TGeoRotation rotHead("",0,90,0);
//   TGeoCombiTrans *rotHeadTr = new TGeoCombiTrans(0,fgkWaferThickness/2,
// 		  -headRadius+fgkTransitHVHeadLZ+fgkTransitHVBondingLZ/2,
// 						 &rotHead);
  TGeoRotation *rotHead = new TGeoRotation("",0,90,0);
  TGeoCombiTrans *rotHeadTr = new TGeoCombiTrans(0,fgkWaferThickness/2,
		  -headRadius+fgkTransitHVHeadLZ+fgkTransitHVBondingLZ/2,
						 rotHead);

  fSDDsensor3->AddNode(vHeadPolyComp,1,rotHeadTr);
  fSDDsensor3->AddNode(vHeadAlComp,1,rotHeadTr);
  fSDDsensor4->AddNode(vHeadPolyComp,1,rotHeadTr);
  fSDDsensor4->AddNode(vHeadAlComp,1,rotHeadTr);

  //---
  AliITSv11GeomCableFlat transitHVCable("ITSsddHVtransitCenter",
					fgkTransitHVBondingLZ,
			    fgkTransitHVPolyThick+fgkTransitHVAlThick);
  transitHVCable.SetNLayers(2);
  transitHVCable.SetLayer(0, fgkTransitHVPolyThick,polyhamideSDD,
			  fColorPolyhamide);
  transitHVCable.SetLayer(1, fgkTransitHVAlThick, alSDD, fColorAl);
  transitHVCable.SetInitialNode(fSDDsensor3);

  x1[0] = -fgkTransitHVHeadLX/2;
  x2[0] = -x1[0];
  x1[1] = (fgkWaferThickness+fgkTransitHVPolyThick+fgkTransitHVAlThick)/2;
  x2[1] = x1[1];
  x1[2] = 0;
  x2[2] = 0;
  transitHVCable.AddCheckPoint(fSDDsensor3, 0, x1, vX);
  transitHVCable.AddCheckPoint(fSDDsensor3, 1, x2, vX);
  segment = transitHVCable.CreateAndInsertBoxCableSegment(1,-90,&ctSegment);
  fSDDsensor4->AddNode(segment, 1, ctSegment);

  transitHVCable.ResetPoints();
  transitHVCable.SetName("ITSsddHVtransitTail");
  transitHVCable.SetWidth(fgkTransitHVtailWidth);
  x1[0] = fgkTransitHVtailXpos;
  x2[0] = fgkTransitHVtailXpos;
  x1[2] = -fgkTransitHVBondingLZ/2;
  x2[2] = -fgkTransitHVBondingLZ/2-fgkTransitHVtailLength;
  Double_t vZ[3] = {0,0,1};
  transitHVCable.AddCheckPoint(fSDDsensor3, 0, x1, vZ);
  transitHVCable.AddCheckPoint(fSDDsensor3, 1, x2, vZ);
  segment = transitHVCable.CreateAndInsertBoxCableSegment(1,0, &ctSegment);
  fSDDsensor4->AddNode(segment, 1, ctSegment);

  //---
  TGeoArb8 *sideLeft = new TGeoArb8( fgkTransitHVPolyThick/2 );
  sideLeft->SetVertex(0, fgkTransitHVtailXpos+fgkTransitHVtailWidth/2,0);
  sideLeft->SetVertex(1, fgkTransitHVtailXpos+fgkTransitHVtailWidth/2,
		      fgkTransitHVsideLZ);
  sideLeft->SetVertex(2, fgkTransitHVHeadLX/2, fgkTransitHVsideLeftZ);
  sideLeft->SetVertex(3, fgkTransitHVHeadLX/2, 0);
  sideLeft->SetVertex(4, fgkTransitHVtailXpos+fgkTransitHVtailWidth/2,0);
  sideLeft->SetVertex(5, fgkTransitHVtailXpos+fgkTransitHVtailWidth/2,
		      fgkTransitHVsideLZ);
  sideLeft->SetVertex(6, fgkTransitHVHeadLX/2, fgkTransitHVsideLeftZ);
  sideLeft->SetVertex(7, fgkTransitHVHeadLX/2, 0);

  TGeoArb8 *sideLeftAl = new TGeoArb8( fgkTransitHVAlThick/2 );
  sideLeftAl->SetVertex(0, fgkTransitHVtailXpos+fgkTransitHVtailWidth/2,0);
  sideLeftAl->SetVertex(1, fgkTransitHVtailXpos+fgkTransitHVtailWidth/2,
			fgkTransitHVsideLZ);
  sideLeftAl->SetVertex(2, fgkTransitHVHeadLX/2, fgkTransitHVsideLeftZ);
  sideLeftAl->SetVertex(3, fgkTransitHVHeadLX/2, 0);
  sideLeftAl->SetVertex(4, fgkTransitHVtailXpos+fgkTransitHVtailWidth/2,0);
  sideLeftAl->SetVertex(5, fgkTransitHVtailXpos+fgkTransitHVtailWidth/2,
			fgkTransitHVsideLZ);
  sideLeftAl->SetVertex(6, fgkTransitHVHeadLX/2, fgkTransitHVsideLeftZ);
  sideLeftAl->SetVertex(7, fgkTransitHVHeadLX/2, 0);

  // sideRight is not there actually
//   TGeoArb8 *sideRight = new TGeoArb8( fgkTransitHVPolyThick/2 );
//   sideRight->SetVertex(0, fgkTransitHVtailXpos-fgkTransitHVtailWidth/2,0);
//   sideRight->SetVertex(1, fgkTransitHVtailXpos-fgkTransitHVtailWidth/2,
// 		       fgkTransitHVsideLZ);
//   sideRight->SetVertex(2, -fgkTransitHVHeadLX/2, fgkTransitHVsideRightZ);
//   sideRight->SetVertex(3, -fgkTransitHVHeadLX/2, 0);
//   sideRight->SetVertex(4, fgkTransitHVtailXpos-fgkTransitHVtailWidth/2,0);
//   sideRight->SetVertex(5, fgkTransitHVtailXpos-fgkTransitHVtailWidth/2,
// 		       fgkTransitHVsideLZ);
//   sideRight->SetVertex(6, -fgkTransitHVHeadLX/2, fgkTransitHVsideRightZ);
//   sideRight->SetVertex(7, -fgkTransitHVHeadLX/2, 0);

//   TGeoRotation rotSide("",0,-90,0);
//   TGeoCombiTrans *sideRightTr = new TGeoCombiTrans(0,
// 				(fgkWaferThickness+fgkTransitHVPolyThick)/2,
// 			        -fgkTransitHVBondingLZ/2,&rotSide);
//   TGeoCombiTrans *sideLeftTr = new TGeoCombiTrans(0,
// 			       (fgkWaferThickness+fgkTransitHVPolyThick)/2,
// 			       -fgkTransitHVBondingLZ/2, &rotSide);
//   TGeoCombiTrans *sideLeftAlTr = new TGeoCombiTrans(0,
// 		  fgkTransitHVPolyThick+(fgkWaferThickness+fgkTransitHVAlThick)/2,
// 	       	  -fgkTransitHVBondingLZ/2, &rotSide);
  TGeoRotation *rotSide = new TGeoRotation("",0,-90,0);
//   TGeoCombiTrans *sideRightTr = new TGeoCombiTrans(0,
// 				(fgkWaferThickness+fgkTransitHVPolyThick)/2,
// 			        -fgkTransitHVBondingLZ/2,rotSide);
  TGeoCombiTrans *sideLeftTr = new TGeoCombiTrans(0,
			       (fgkWaferThickness+fgkTransitHVPolyThick)/2,
			       -fgkTransitHVBondingLZ/2, rotSide);
  TGeoCombiTrans *sideLeftAlTr = new TGeoCombiTrans(0,
		  fgkTransitHVPolyThick+(fgkWaferThickness+fgkTransitHVAlThick)/2,
	       	  -fgkTransitHVBondingLZ/2, rotSide);

  TGeoVolume *vSideLeft = new TGeoVolume("ITSsddHVtransitSideLeft",
					 sideLeft,polyhamideSDD);
  vSideLeft->SetLineColor(fColorPolyhamide);
  TGeoVolume *vSideLeftAl = new TGeoVolume("ITSsddHVtransitSideLeftAl",
					   sideLeftAl,alSDD);
  vSideLeftAl->SetLineColor(fColorAl);

//   TGeoVolume *vSideRight = new TGeoVolume("ITSsddHVtransitSideRight",
// 					  sideRight,polyhamideSDD);
//   vSideRight->SetLineColor(fColorPolyhamide);

  fSDDsensor3->AddNode(vSideLeft,   1, sideLeftTr);
  fSDDsensor3->AddNode(vSideLeftAl, 1, sideLeftAlTr);
//   fSDDsensor3->AddNode(vSideRight,  1, sideRightTr);

  fSDDsensor4->AddNode(vSideLeft,   1, sideLeftTr);
  fSDDsensor4->AddNode(vSideLeftAl, 1, sideLeftAlTr);
//   fSDDsensor4->AddNode(vSideRight,  1, sideRightTr);
  };

  //****************************
  if(GetDebug(1)) {
    fSDDsensor3->CheckOverlaps(0.01);
    fSDDsensor4->CheckOverlaps(0.01);
  }

  fSDDsensor3->SetVisibility(kFALSE);
  fSDDsensor4->SetVisibility(kFALSE);
}

/*
//________________________________________________________________________
TGeoVolume *AliITSv11GeometrySDD::CreateDetectors(Int_t iLay) {
  //
  // return a box volume containing the detectors
  //

  TGeoMedium *airSDD = GetMedium("SDD AIR$");

  Int_t    nDetectors   = fgkLay3Ndet;
  Double_t ladderLength = fgkLay3LadderLength;
  Double_t *sensorZPos  = fLay3sensorZPos;
  
  if (iLay==3) {}
  else if (iLay==4) {
    nDetectors   = fgkLay4Ndet;
    ladderLength = fgkLay4LadderLength;
    sensorZPos   = fLay4sensorZPos;
  } else {
    printf("AliITSv11GeometrySDD::CreateDetectors: Error : Wrong layer");
  };

  char name[30];
  Double_t volThickness = ( fgkLadWaferSep + 2*fgkWaferThickness +
			    2*(fgkWaHVcableAlThick+fgkWaHVcablePolyThick));
  
  sprintf(name,"ITSsddDetBox%i",iLay);
  TGeoBBox *detBox = new TGeoBBox(name, fgkWaferWidth/2, volThickness/2,
			 ladderLength*((nDetectors-0.5)/nDetectors)/2);
  TGeoVolume *virtualDet = new TGeoVolume("ITSsddLadd",detBox, airSDD);
 
    for (Int_t i=0; i<nDetectors; i++) {
        Double_t localZ = sensorZPos[i];
        Double_t localY = fgkLadWaferSep/2+fgkWaferThickness/2;
        if (iLay==3) if (i%2!=0) localY = -localY;
	if (iLay==4) if (i%2==0) localY = -localY;
        sprintf(name, "ITSsddLay%iSensorPos%i",iLay, i);

	if (i >= nDetectors/2) {
	  TGeoTranslation *sensorPos = new TGeoTranslation(0,localY,localZ);
	  sensorPos->SetName(name);
	  virtualDet->AddNode(fSDDsensor, i, sensorPos);
	}
	else {
	  TGeoRotation *rotSensor = new TGeoRotation("",0, 180, 180);
	  TGeoCombiTrans *sensorPos = new TGeoCombiTrans(0,localY,
				      		 localZ, rotSensor);
	  sensorPos->SetName(name);
	  virtualDet->AddNode(fSDDsensor, i, sensorPos);
	};
    }

    if(GetDebug(1)) virtualDet->CheckOverlaps(0.01);
    virtualDet->SetVisibility(kFALSE);
    return virtualDet;
}
*/

//________________________________________________________________________
TGeoVolumeAssembly *AliITSv11GeometrySDD::CreateDetectorsAssembly(Int_t iLay) {
//
// return a box volume containing the detectors
//
  
  Int_t    nDetectors   = fgkLay3Ndet;
  Double_t ladderLength = fgkLay3LadderLength;
  Double_t *sensorZPos  = fLay3sensorZPos;
  TGeoVolume *sensorSDD = fSDDsensor3;

  if (iLay==3) {}
  else if (iLay==4) {
    nDetectors   = fgkLay4Ndet;
    ladderLength = fgkLay4LadderLength;
    sensorZPos   = fLay4sensorZPos;
    sensorSDD    = fSDDsensor4;
  } else {
    printf("AliITSv11GeometrySDD::CreateDetectorsAssembly: Error:Wrong layer");
  };

  char name[30];
  sprintf(name,"ITSsddDetBox%i",iLay);
  
  TGeoVolumeAssembly  *virtualDet = new TGeoVolumeAssembly("ITSsddLadd");

  for (Int_t i=0; i<nDetectors; i++) {
    Double_t localZ = sensorZPos[i];
    Double_t localY = fgkLadWaferSep/2+fgkWaferThickness/2;
    if (iLay==3) if (i%2!=0) localY = -localY;
    if (iLay==4) if (i%2==0) localY = -localY;
    sprintf(name, "ITSsddLay%iSensorPos%i",iLay, i);
 
    if (i >= nDetectors/2) {
      TGeoTranslation *sensorPos = new TGeoTranslation(0,localY,localZ);
      sensorPos->SetName(name);
      virtualDet->AddNode(sensorSDD, i, sensorPos);
    }
    else {
      TGeoRotation *rotSensor = new TGeoRotation("",0, 180, 180);
      TGeoCombiTrans *sensorPos = new TGeoCombiTrans(0,localY,
						     localZ, rotSensor);
      sensorPos->SetName(name);
      virtualDet->AddNode(sensorSDD, i, sensorPos);
    };
  }
  
  if(GetDebug(1)) virtualDet->CheckOverlaps(0.01);
  return virtualDet;
}


//________________________________________________________________________
TGeoVolumeAssembly *AliITSv11GeometrySDD::CreateDetectorsAssemblyLadd2() {
//
// return a box volume containing the detectors
// Special case for Layer 3 Ladder 2 which is rotated (cannot simply
// rotate the standard volume, because the module numbering would be wrong)
// M.Sitta 25 Nov 2009
//
  
  Int_t    nDetectors   = fgkLay3Ndet;
  Double_t *sensorZPos  = fLay3sensorZPos;
  TGeoVolume *sensorSDD = fSDDsensor3;

  char name[30];
  sprintf(name,"ITSsddDetBoxLadd2");
  
  TGeoVolumeAssembly  *virtualDet = new TGeoVolumeAssembly("ITSsddLadd");

  for (Int_t i=0; i<nDetectors; i++) {
    Double_t localZ = (-1.)*sensorZPos[nDetectors-1-i];
    Double_t localY = fgkLadWaferSep/2+fgkWaferThickness/2;
    if (i%2==0) localY = -localY;
    sprintf(name, "ITSsddLayLadd2SensorPos%i", i);
 
    if (i >= nDetectors/2) {
      TGeoTranslation *sensorPos = new TGeoTranslation(0,localY,localZ);
      sensorPos->SetName(name);
      virtualDet->AddNode(sensorSDD, i, sensorPos);
    }
    else {
      TGeoRotation *rotSensor = new TGeoRotation("",0, 180, 180);
      TGeoCombiTrans *sensorPos = new TGeoCombiTrans(0,localY,
						     localZ, rotSensor);
      sensorPos->SetName(name);
      virtualDet->AddNode(sensorSDD, i, sensorPos);
    };
  }
  
  if(GetDebug(1)) virtualDet->CheckOverlaps(0.01);
  return virtualDet;
}


//________________________________________________________________________
Int_t AliITSv11GeometrySDD::ExportSensorGeometry(AliITSgeom *geom, Int_t iLaySDD,
						 Int_t startMod) {
//
// export the geometry in a AliITSgeom object
// Obsolete
//

  if (! geom) {
    printf("error:Try to fill null (AliITSgeom *) object");
    return kFALSE;
  };
  if (! fMotherVol) {
    printf("error:Try to set sensor geometry while geometry is not defined\n");
    return kFALSE;
  };

  const Float_t kDxyz[3] = {fgkWaferWidthSens/2., fgkWaferThickSens/2.,
			    fgkWaferLengthSens/2.};
  if(!(geom->IsShapeDefined(kSDD)))
    geom->ReSetShape(kSDD, new AliITSgeomSDD256(3, kDxyz));

  char layerName[30];
  char ladderName[30];
  char sensorName[30];
  char senstivName[30];
  const Int_t kNLay = 2;
  const Int_t kNLadd[kNLay] = {fgkLay3Nladd, fgkLay4Nladd};
  const Int_t kNDet[kNLay]  = {fgkLay3Ndet,  fgkLay4Ndet};

  if (GetDebug(1))
    printf("AliITSv11GeometrySDD::SetSensorGeometry(), nodes found :\n");

  Int_t firstSDDmod = startMod;
  for (Int_t iLay=0; iLay<kNLay; iLay++) {
    /////////////////////////////////////////
    sprintf(layerName, "ITSsddLayer%i_1",iLay+3);
    TGeoNode *layNode = fMotherVol->GetNode(layerName);
    if (layNode) {
      if (GetDebug(1)) printf("%s\n",layNode->GetName());
      TGeoVolume *layVolume = layNode->GetVolume();
      TGeoHMatrix layMatrix(*layNode->GetMatrix());

      for (Int_t iLadd=0; iLadd<kNLadd[iLay]; iLadd++) {
	/////////////////////////////////////////
	sprintf(ladderName, "ITSsddLadd_%i", iLadd);
	TGeoNode *laddNode = layVolume->GetNode(ladderName);
	if (laddNode) {
	  if (GetDebug(1)) printf("|   %s\n",laddNode->GetName());
	  TGeoVolume *laddVolume = laddNode->GetVolume();
	  TGeoHMatrix laddMatrix(layMatrix);
	  laddMatrix.Multiply(laddNode->GetMatrix());

	  for (Int_t iDet=0; iDet<kNDet[iLay]; iDet++) {
	    /////////////////////////////////////////
	    sprintf(sensorName, "ITSsddSensor_%i",iDet);
	    TGeoNode *detNode = laddVolume->GetNode(sensorName);
	    if (detNode) {
	      if (GetDebug(1)) printf("|   |   %s\n",detNode->GetName());
	      TGeoVolume *detVolume = detNode->GetVolume();
	      TGeoHMatrix detMatrix(laddMatrix);
	      detMatrix.Multiply(detNode->GetMatrix());

	      TGeoNode *wafNode = detVolume->GetNode("ITSsddWafer_1");
	      if (wafNode) {
		TGeoVolume *wafVolume = wafNode->GetVolume();
		TGeoHMatrix wafMatrix(detMatrix);
		detMatrix.Multiply(wafNode->GetMatrix());
		//--------------------------------------------------------
		sprintf(senstivName, "%s%s", fgSDDsensitiveVolName3,"_1");
		TGeoNode *sensitivNode = wafVolume->GetNode(senstivName);
	      if (sensitivNode) {
		TGeoHMatrix sensMatrix(wafMatrix);
		sensMatrix.Multiply(sensitivNode->GetMatrix());

		// Sticking to the convention for local wafer coordinate
		// in AliITSgeom :
		if (iDet >= kNDet[iLay]/2) {
		  //		  TGeoRotation rotY("",0,180,0);
		  TGeoRotation rotY("",-180,-180,0);
		  sensMatrix.Multiply(&rotY);
		};
		// Creating the matrix in AliITSgeom for
		// this sensitive volume :
		Double_t *trans = sensMatrix.GetTranslation();
		Double_t *r     = sensMatrix.GetRotationMatrix();
		Double_t rot[10] = {r[0],r[1],r[2],
				    r[3],r[4],r[5],
				    r[6],r[7],r[8], 1.0};
		//rot[9]!=0.0 => not a unity matrix
		geom->CreateMatrix(startMod,iLay+iLaySDD,iLadd+1,iDet+1,
				  kSDD,trans,rot);
		// iLadd+1, iDet+1 because ladd. and det. start at +1
		// elsewhere
		startMod++;

	      } else
		printf("Error (ExportSensorGeometry) %s not found !\n",
		       senstivName);
	      } else
		printf("Error (ExportSensorGeometry) %s not found !\n",
		       "ITSsddWafer_1");
	    } else
	      printf("Error (ExportSensorGeometry) %s not found !\n",
		     sensorName);
	  };
	} else 
	  printf("Error (ExportSensorGeometry) %s not found !\n",
		 ladderName);
      };  
    } else
      printf("Error (ExportSensorGeometry) %s not found !\n",
	     layerName);
  };

  return (startMod-firstSDDmod);
}


//________________________________________________________________________
Int_t AliITSv11GeometrySDD::
GetCurrentLayLaddDet(Int_t &lay, Int_t &ladd, Int_t&det) const {
//
// Function which gives the layer, ladder and det.
// index of the current volume. To be used in
// AliITS::StepManager()
  //

  if (gGeoManager->GetLevel()<3) return kFALSE;
  // Get the det index :
  TGeoNode *node = gGeoManager->GetMother(2);
  if (!node) return kFALSE;
  det = node->GetNumber()+1;

  // Get the ladder index :
  node = gGeoManager->GetMother(3);
  if (!node) return kFALSE;
  ladd = node->GetNumber()+1;

 // Get the layer index :
  if (node->GetNdaughters()==fgkLay3Ndet)
    lay = 3;            // this has to be equal to the iLaySDD argument given to ExportSensorGeometry() !!!
  else lay = 4;

  return kTRUE;
}


//________________________________________________________________________
TGeoPcon* AliITSv11GeometrySDD::CreateConeConstSection(Double_t r1max, Double_t z1,
						       Double_t r2max, Double_t z2,
						       Double_t section, Int_t nDiv)
{
  // Creates a cone along z where the section is approximately constant
  // with z. This is for simulation of cables, because a cone with a constant
  // radius difference would show a quantity of matter increasing with z...
  // The max radius of the created Pcon is evolving linearly, the min radius
  // is calculated at several steps (nDiv).
  // z2 > z1 (required by the Pcon)

  TGeoPcon *myPcon = new TGeoPcon(0, 360, 1+nDiv);
  
  Double_t dr = (r2max-r1max)/nDiv;
  Double_t dz = (z2-z1)/nDiv;
  Double_t r1minI, r2minI, r1maxI, r2maxI;
  Double_t z1I, z2I;

  Double_t lZ = TMath::Sqrt((r2max-r1max)*(r2max-r1max) + (z2-z1)*(z2-z1));
  Double_t cosAlpha = (z2-z1)/lZ;

  r1minI = TMath::Sqrt(r1max*r1max-section/(TMath::Pi()*cosAlpha));
  myPcon->DefineSection(0, z1, r1minI, r1max);

  for (Int_t i=0; i<nDiv; i++) {
    
    z1I = z1 + i*dz;
    z2I = z1I + dz;
    r1maxI = r1max + i*dr;
    r2maxI = r1maxI + dr;

    r2minI =  TMath::Sqrt(r2maxI*r2maxI-section/(TMath::Pi()*cosAlpha));
    myPcon->DefineSection(i+1, z2I, r2minI, r2maxI);
  }
  return myPcon;
}


//________________________________________________________________________
Double_t AliITSv11GeometrySDD::GetConeZ(Double_t r, Double_t refR1, Double_t refR2, 
					Double_t refZ1, Double_t refZ2) {
  // just a helping function
  return refZ1+(refZ2-refZ1)*(r-refR1)/(refR2-refR1);
}

//________________________________________________________________________
Int_t AliITSv11GeometrySDD::CreateAndInsetConeCablePart(TGeoVolume *mother, Double_t angle,
							Int_t nLay3, Int_t nLay4,
							Double_t r1, Double_t z1,
							Double_t r2, Double_t z2) {
  
  // Create some cables portions from SDD modules grouped
  // and attached at the border of the SSD cone

  TGeoMedium *copper     = GetMedium("COPPER$");
  TGeoMedium *plastic    = GetMedium("SDDKAPTON (POLYCH2)$");
  TGeoMedium *opticalFiber = GetMedium("SDD OPTICFIB$");

  char titleCable[30];
  sprintf(titleCable,"cableSDDport%i",(Int_t)angle);

  //---
  Double_t section = (fgkSectionCuPerMod+fgkSectionPlastPerMod+fgkSectionGlassPerMod)*(nLay3+nLay4);
  Double_t thickness = 1.; // let's fix the thickness, then calculate the width
  Double_t width     = section/thickness;
  Double_t thickCu   = thickness*fgkSectionCuPerMod/(fgkSectionCuPerMod+fgkSectionPlastPerMod
						     +fgkSectionGlassPerMod);

  Double_t thickPlast = thickness*fgkSectionPlastPerMod/(fgkSectionCuPerMod+fgkSectionPlastPerMod
							 +fgkSectionGlassPerMod);

  Double_t thickGlass = thickness*fgkSectionGlassPerMod/(fgkSectionCuPerMod+fgkSectionPlastPerMod
							 +fgkSectionGlassPerMod);

  Double_t hypothenus   = TMath::Sqrt( (r2-r1)*(r2-r1) + (z2-z1)*(z2-z1) );
  Double_t cosAlpha     = (z2-z1)/hypothenus;
  Double_t radius1Cable = TMath::Sqrt(r1*r1 - width*width/4) - 0.5*thickness/cosAlpha;
  Double_t radius2Cable = TMath::Sqrt(r2*r2 - width*width/4) - 0.5*thickness/cosAlpha;
  angle *= TMath::DegToRad();
  Double_t x1 = radius1Cable*TMath::Cos(angle), y1 = radius1Cable*TMath::Sin(angle);
  Double_t x2 = radius2Cable*TMath::Cos(angle), y2 = radius2Cable*TMath::Sin(angle);
  Double_t pos1[3] = {x1,y1,z1};
  Double_t pos2[3] = {x2,y2,z2};
  Double_t zVect[3] = {0,0,1};

  AliITSv11GeomCableFlat cable(titleCable,width,thickness);
  cable.SetNLayers(3);
  cable.SetLayer(0, thickPlast, plastic, kYellow);
  cable.SetLayer(1, thickCu, copper, kRed);
  cable.SetLayer(2, thickGlass, opticalFiber, kGreen);

  cable.AddCheckPoint( mother, 0, pos1, zVect );
  cable.AddCheckPoint( mother, 1, pos2, zVect );
  cable.SetInitialNode(mother);
  cable.CreateAndInsertCableSegment(1);

  return kTRUE;
}



//________________________________________________________________________
void AliITSv11GeometrySDD::SDDCables(TGeoVolume *moth)
{
//
// Creates and inserts the SDD cables running on SDD and SSD cones
//
// Input:
//         moth : the TGeoVolume owing the volume structure
// Output:
//
// Created:         ???       Ludovic Gaudichet
// Updated:      15 Mar 2008  Mario Sitta
// Updated:      14 Apr 2008  Mario Sitta            Overlap fixes
// Updated:      09 May 2008  Mario Sitta            SSD overlap fixes
//

  TGeoMedium *copper       = GetMedium("COPPER$");
  TGeoMedium *plastic      = GetMedium("SDDKAPTON (POLYCH2)$");
  TGeoMedium *opticalFiber = GetMedium("SDD OPTICFIB$");
  TGeoMedium *airSDD       = GetMedium("SDD AIR$");


  //==================================
  //  
  //==================================

  Double_t nModLay3 = fgkLay3Nladd*fgkLay3Ndet;
  Double_t nModLay4 = fgkLay4Nladd*fgkLay4Ndet;

  Double_t sectionLay3Cu      = fgkCableBendRatio*fgkSectionCuPerMod*nModLay3/2;
  Double_t sectionLay3Plastic = fgkCableBendRatio*fgkSectionPlastPerMod*nModLay3/2;
  Double_t sectionLay3Glass   = fgkCableBendRatio*fgkSectionGlassPerMod*nModLay3/2;

  Double_t sectionLay4Cu      = fgkCableBendRatio*fgkSectionCuPerMod*nModLay4/2;
  Double_t sectionLay4Plastic = fgkCableBendRatio*fgkSectionPlastPerMod*nModLay4/2;
  Double_t sectionLay4Glass   = fgkCableBendRatio*fgkSectionGlassPerMod*nModLay4/2;

  // Do not use hardcoded numbers, get them from real shapes - M.S. 15/03/08
  TGeoVolume *sddCone = gGeoManager->GetVolume("SDDCarbonFiberCone");
  TGeoPcon *sddConeShape = (TGeoPcon*)sddCone->GetShape();

  TGeoVolume *sddCylinder = gGeoManager->GetVolume("SDDCarbonFiberCylinder");
  TGeoTube *sddCylinderShape = (TGeoTube*)sddCylinder->GetShape();

  // (were fgkConeSDDr1, fgkConeSDDr2, fgkConeSDDz1, fgkConeSDDz2 hardcoded)
  Double_t coneSDDr1 = sddConeShape->GetRmin(5);
  Double_t coneSDDr2 = sddConeShape->GetRmin(3);

  Double_t coneSDDz1 = sddConeShape->GetZ(9) - sddConeShape->GetZ(5) +
		       sddCylinderShape->GetDz();
  Double_t coneSDDz2 = sddConeShape->GetZ(9) - sddConeShape->GetZ(3) +
		       sddCylinderShape->GetDz();

  // Calculate z1, z2 thanks to R1 and R2
  Double_t sddCableZ1 = GetConeZ(fgkSDDCableR1, coneSDDr1, coneSDDr2,
				                coneSDDz1, coneSDDz2);
  Double_t sddCableZ2 = GetConeZ(fgkSDDCableR2, coneSDDr1, coneSDDr2,
				                coneSDDz1, coneSDDz2);
  Double_t sddCableZ3 = GetConeZ(fgkSDDCableR3, coneSDDr1, coneSDDr2,
				                coneSDDz1, coneSDDz2);

  TGeoRotation *rotCableSDD = new TGeoRotation("rotCableSDD",0,180,0);

  //==================================
  //  first set of cones : cables from layer 3
  //==================================

  TGeoPcon* pcon1all = CreateConeConstSection(fgkSDDCableR1, sddCableZ1,
					      fgkSDDCableR2, sddCableZ2,
		  sectionLay3Plastic+sectionLay3Cu+sectionLay3Glass, 1);

  TGeoPcon* pcon1container = new TGeoPcon(0,360,2);
  pcon1container->DefineSection(0, sddCableZ1, pcon1all->GetRmin(0),
					       pcon1all->GetRmax(0));

  Double_t drMax = pcon1all->GetRmax(0)- pcon1all->GetRmin(0);
  pcon1container->DefineSection(1, sddCableZ2, pcon1all->GetRmax(1)-drMax,
					       pcon1all->GetRmax(1));

  TGeoVolume *vpcon1container = new TGeoVolume("vpcon1container",
					       pcon1container, airSDD);
  vpcon1container->SetVisibility(kFALSE);

  TGeoPcon* pcon1plast = CreateConeConstSection(fgkSDDCableR1, sddCableZ1,
						fgkSDDCableR2, sddCableZ2,
						sectionLay3Plastic, 3);

  TGeoVolume *vpcon1plast = new TGeoVolume("ITScablesSDDpcon1Plast",
					   pcon1plast, plastic);
  vpcon1plast->SetLineColor(kYellow);
  vpcon1container->AddNode(vpcon1plast, 0);

  Double_t dr1a = fgkSDDCableR1 - pcon1plast->GetRmin(0);
  TGeoPcon* pcon1Cu = CreateConeConstSection(fgkSDDCableR1 - dr1a, sddCableZ1,
					     fgkSDDCableR2 - dr1a, sddCableZ2,
					     sectionLay3Cu, 3);

  TGeoVolume *vpcon1Cu = new TGeoVolume("ITScablesSDDpcon1Cu",
					pcon1Cu, copper);
  vpcon1Cu->SetLineColor(kRed);
  vpcon1container->AddNode(vpcon1Cu, 0);

  Double_t dr1b = pcon1Cu->GetRmax(0) - pcon1Cu->GetRmin(0);
  TGeoPcon* pcon1glass = CreateConeConstSection(fgkSDDCableR1-dr1a-dr1b, sddCableZ1,
						fgkSDDCableR2-dr1a-dr1b, sddCableZ2,
						sectionLay3Glass, 3);

  TGeoVolume *vpcon1glass = new TGeoVolume("ITScablesSDDpcon1glass",
					   pcon1glass, opticalFiber);
  vpcon1glass->SetLineColor(kGreen);
  vpcon1container->AddNode(vpcon1glass, 0);

  moth->AddNode(vpcon1container, 1);
  moth->AddNode(vpcon1container, 2, rotCableSDD);

  //==================================
  //  2nd set of cones : cables from layer 3 and layer 4
  //==================================

  TGeoPcon* pcon2all = CreateConeConstSection(fgkSDDCableR2, sddCableZ2,
					      fgkSDDCableR3, sddCableZ3,
				      sectionLay3Plastic+sectionLay4Plastic+
				      sectionLay3Cu+sectionLay4Cu+
				      sectionLay3Glass+sectionLay4Glass, 1);

  TGeoPcon* pcon2container = new TGeoPcon(0,360,2);
  pcon2container->DefineSection(0, sddCableZ2, pcon2all->GetRmin(0),
					       pcon2all->GetRmax(0));

  drMax = pcon2all->GetRmax(0)- pcon2all->GetRmin(0);
  pcon2container->DefineSection(1, sddCableZ3, pcon2all->GetRmax(1)-drMax,
					       pcon2all->GetRmax(1));


  TGeoVolume *vpcon2container = new TGeoVolume("vpcon2container",
					       pcon2container, airSDD);
  vpcon2container->SetVisibility(kFALSE);

  TGeoPcon* pcon2plast = CreateConeConstSection(fgkSDDCableR2, sddCableZ2,
						fgkSDDCableR3, sddCableZ3,
						sectionLay3Plastic+
						sectionLay4Plastic, 3);

  TGeoVolume *vpcon2plast = new TGeoVolume("ITScablesSDDpcon2Plast",
					   pcon2plast, plastic);
  vpcon2plast->SetLineColor(kYellow);
  vpcon2container->AddNode(vpcon2plast, 0);

  Double_t dr2a = fgkSDDCableR2 - pcon2plast->GetRmin(0);
  TGeoPcon* pcon2Cu = CreateConeConstSection(fgkSDDCableR2 - dr2a, sddCableZ2,
					     fgkSDDCableR3 - dr2a, sddCableZ3,
					     sectionLay3Cu+sectionLay4Cu, 3);

  TGeoVolume *vpcon2Cu = new TGeoVolume("ITScablesSDDpcon2Cu",
					pcon2Cu, copper);
  vpcon2Cu->SetLineColor(kRed);
  vpcon2container->AddNode(vpcon2Cu, 0);

  Double_t dr2b = pcon2Cu->GetRmax(0) - pcon2Cu->GetRmin(0);
  TGeoPcon* pcon2glass = CreateConeConstSection(fgkSDDCableR2-dr2a-dr2b, sddCableZ2,
						fgkSDDCableR3-dr2a-dr2b, sddCableZ3,
						sectionLay3Glass+
						sectionLay4Glass, 3);

  TGeoVolume *vpcon2glass = new TGeoVolume("ITScablesSDDpcon2glass",
					   pcon2glass, opticalFiber);
  vpcon2glass->SetLineColor(kGreen);
  vpcon2container->AddNode(vpcon2glass, 0);

  moth->AddNode(vpcon2container, 1);
  moth->AddNode(vpcon2container, 2, rotCableSDD);

  //==================================
  //  intermediate cylinder
  //==================================

  // (was fgkSDDCableDZint hardcoded)
  Double_t sddCableDZint = (sddConeShape->GetZ(9) - sddConeShape->GetZ(0) +
		            sddCylinderShape->GetDz()) - sddCableZ3;

  TGeoTube *interCyl = new TGeoTube("sddCableInterCyl",
				    pcon2container->GetRmin(1),
				    pcon2container->GetRmax(1),
				    sddCableDZint/2);

  TGeoVolume *vInterCyl = new TGeoVolume("vSddCableInterCyl",
					 interCyl, airSDD);
  vInterCyl->SetVisibility(kFALSE);

  Double_t rmaxCylPlast = pcon2container->GetRmax(1);
  Double_t rminCylPlast = TMath::Sqrt(rmaxCylPlast*rmaxCylPlast - 
			(sectionLay3Plastic+sectionLay4Plastic)/TMath::Pi() );

  TGeoTube *interCylPlast = new TGeoTube("sddCableInterCylPlast", rminCylPlast,
					 rmaxCylPlast, sddCableDZint/2);

  TGeoVolume *vInterCylPlast = new TGeoVolume("vSddCableInterCylPlast",
					      interCylPlast, plastic);
  vInterCylPlast->SetLineColor(kYellow);
  vInterCyl->AddNode(vInterCylPlast, 0);

  Double_t rmaxCylCu = pcon2Cu->GetRmax(3);
  Double_t rminCylCu = TMath::Sqrt(rmaxCylCu*rmaxCylCu - 
		       (sectionLay3Cu+sectionLay4Cu)/TMath::Pi() );

  TGeoTube *interCylCu = new TGeoTube("sddCableInterCylCu", rminCylCu,
				      rmaxCylCu, sddCableDZint/2);

  TGeoVolume *vInterCylCu = new TGeoVolume("vSddCableInterCylCu",
					   interCylCu, copper);
  vInterCylCu->SetLineColor(kRed);
  vInterCyl->AddNode(vInterCylCu, 0);

  Double_t rmaxCylGlass = pcon2glass->GetRmax(3);
  Double_t rminCylGlass = TMath::Sqrt(rmaxCylGlass*rmaxCylGlass - 
			  (sectionLay3Glass+sectionLay4Glass)/TMath::Pi() );

  TGeoTube *interCylGlass = new TGeoTube("sddCableInterCylGlass", rminCylGlass,
					 rmaxCylGlass, sddCableDZint/2);

  TGeoVolume *vInterCylGlass = new TGeoVolume("vSddCableInterCylGlass",
					      interCylGlass,opticalFiber);
  vInterCylGlass->SetLineColor(kGreen);
  vInterCyl->AddNode(vInterCylGlass, 0);

  moth->AddNode(vInterCyl, 1, new TGeoTranslation(0, 0,
			      sddCableZ3+sddCableDZint/2));
  moth->AddNode(vInterCyl, 2, new TGeoTranslation(0, 0,
			     -sddCableZ3-sddCableDZint/2));

  //==================================
  // cable cone on the SSD cone
  //==================================

  Double_t sddCableR4 = rmaxCylPlast;
  Double_t sddCableZ4 = sddCableZ3 + sddCableDZint;

  TGeoPcon* pcon3all = CreateConeConstSection(sddCableR4, sddCableZ4,
					      fgkSDDCableR5, fgkSDDCableZ5,
					      sectionLay3Plastic+
					      sectionLay4Plastic+
					      sectionLay3Cu+sectionLay4Cu+
					      sectionLay3Glass+sectionLay4Glass, 1);

  TGeoPcon* pcon3container = new TGeoPcon(0,360,2);
  pcon3container->DefineSection(0, sddCableZ4, pcon3all->GetRmin(0),
					       pcon3all->GetRmax(0));

  drMax = pcon3all->GetRmax(0) - pcon3all->GetRmin(0);
  pcon3container->DefineSection(1, fgkSDDCableZ5, pcon3all->GetRmax(1)-drMax,
					       pcon3all->GetRmax(1));


  TGeoVolume *vpcon3container = new TGeoVolume("vpcon3container",
					       pcon3container, airSDD);
  vpcon3container->SetVisibility(kFALSE);

  TGeoPcon* pcon3plast = CreateConeConstSection(sddCableR4, sddCableZ4,
						fgkSDDCableR5, fgkSDDCableZ5,
						sectionLay3Plastic+
						sectionLay4Plastic, 3);

  TGeoVolume *vpcon3plast = new TGeoVolume("ITScablesSDDpcon3Plast",
					   pcon3plast, plastic);
  vpcon3plast->SetLineColor(kYellow);
  vpcon3container->AddNode(vpcon3plast, 0);

  Double_t dr3a = sddCableR4 - pcon3plast->GetRmin(0);
  TGeoPcon* pcon3Cu = CreateConeConstSection(sddCableR4 - dr3a, sddCableZ4,
					     fgkSDDCableR5 - dr3a, fgkSDDCableZ5,
					     sectionLay3Cu+sectionLay4Cu, 3);

  TGeoVolume *vpcon3Cu = new TGeoVolume("ITScablesSDDpcon3Cu",
					pcon3Cu, copper);
  vpcon3Cu->SetLineColor(kRed);
  vpcon3container->AddNode(vpcon3Cu, 0);

  Double_t dr3b = pcon3Cu->GetRmax(0) - pcon3Cu->GetRmin(0);
  TGeoPcon* pcon3glass = CreateConeConstSection(sddCableR4-dr3a-dr3b, sddCableZ4,
						fgkSDDCableR5-dr3a-dr3b, fgkSDDCableZ5,
						sectionLay3Glass+sectionLay4Glass, 3);

  TGeoVolume *vpcon3glass = new TGeoVolume("ITScablesSDDpcon3glass",
					   pcon3glass,opticalFiber);
  vpcon3glass->SetLineColor(kGreen);
  vpcon3container->AddNode(vpcon3glass, 0);

  moth->AddNode(vpcon3container, 1);
  moth->AddNode(vpcon3container, 2, rotCableSDD);

  //==================================
  // cables that are grouped at the end of SSD cones
  //==================================

//  Double_t fgkSDDCableR6 = fgkSDDCableR5+9;
//  Double_t fgkSDDCableZ6 = fgkSDDCableZ5+8.8;
  Double_t fgkSDDCableR6 = fgkSDDCableR5+8;
  Double_t fgkSDDCableZ6 = fgkSDDCableZ5+8;

  TGeoVolumeAssembly *endConeSDDCable = new TGeoVolumeAssembly("endConeSDDCable");

  // Add some hardcoded shifts to avoid overlaps with SSD pathc panels
  CreateAndInsetConeCablePart(endConeSDDCable, 20, 1*3,2*4, fgkSDDCableR5,
			      fgkSDDCableZ5,fgkSDDCableR6-2.6,fgkSDDCableZ6-2.6);

  CreateAndInsetConeCablePart(endConeSDDCable, 50, 1*3,1*4, fgkSDDCableR5,
			      fgkSDDCableZ5,fgkSDDCableR6,fgkSDDCableZ6);

  CreateAndInsetConeCablePart(endConeSDDCable, 85, 2*3,1*4, fgkSDDCableR5,
			      fgkSDDCableZ5,fgkSDDCableR6,fgkSDDCableZ6);

  CreateAndInsetConeCablePart(endConeSDDCable, 95, 0*3,1*4, fgkSDDCableR5,
			      fgkSDDCableZ5,fgkSDDCableR6,fgkSDDCableZ6);

  CreateAndInsetConeCablePart(endConeSDDCable, 105, 2*3,3*4, fgkSDDCableR5,
			      fgkSDDCableZ5,fgkSDDCableR6-2.6,fgkSDDCableZ6-2.6);

  CreateAndInsetConeCablePart(endConeSDDCable, 129, 0*3,3*4, fgkSDDCableR5,
			      fgkSDDCableZ5,fgkSDDCableR6,fgkSDDCableZ6);

  CreateAndInsetConeCablePart(endConeSDDCable, 176, 0*3,1*4, fgkSDDCableR5,
			      fgkSDDCableZ5,fgkSDDCableR6,fgkSDDCableZ6);

  CreateAndInsetConeCablePart(endConeSDDCable, 190, 2*3,0*4, fgkSDDCableR5,
			      fgkSDDCableZ5,fgkSDDCableR6,fgkSDDCableZ6);

  CreateAndInsetConeCablePart(endConeSDDCable, 210, 1*3,2*4, fgkSDDCableR5,
			      fgkSDDCableZ5,fgkSDDCableR6-2.6,fgkSDDCableZ6-2.6);

  CreateAndInsetConeCablePart(endConeSDDCable, 230, 1*3,2*4, fgkSDDCableR5,
			      fgkSDDCableZ5,fgkSDDCableR6,fgkSDDCableZ6);

  CreateAndInsetConeCablePart(endConeSDDCable, 277, 2*3,2*4, fgkSDDCableR5,
			      fgkSDDCableZ5,fgkSDDCableR6,fgkSDDCableZ6);

  CreateAndInsetConeCablePart(endConeSDDCable, 306, 1*3,1*4, fgkSDDCableR5,
			      fgkSDDCableZ5,fgkSDDCableR6,fgkSDDCableZ6);

  CreateAndInsetConeCablePart(endConeSDDCable, 353, 1*3,3*4, fgkSDDCableR5,
			      fgkSDDCableZ5,fgkSDDCableR6,fgkSDDCableZ6);

  moth->AddNode(endConeSDDCable, 1, 0);

  TGeoRotation* reflect = new TGeoRotation("reflectEndConeSDDCable");
  reflect->ReflectZ(kTRUE);
  moth->AddNode(endConeSDDCable, 2, reflect);


  return;
}
