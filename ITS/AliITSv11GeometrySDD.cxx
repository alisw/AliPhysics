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
// SDD geometry, based on ROOT geometrical modeler
//
// Ludovic Gaudichet                                   gaudichet@to.infn.it
//*************************************************************************



// General Root includes
//#include <Riostream.h>
#include <TMath.h>

// Root Geometry includes
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoCone.h>
#include <TGeoTube.h>
#include <TGeoArb8.h>
#include <TGeoCompositeShape.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>

#include "AliITSgeom.h"
#include "AliITSgeomSDD.h"
#include "AliITSv11GeometrySDD.h"
#include "AliITSv11GeomCableFlat.h"
#include "AliITSv11GeomCableRound.h"


const char*    AliITSv11GeometrySDD::fgSDDsensitiveVolName = "ITSsddSensitiv";
const Double_t AliITSv11GeometrySDD::fgkSegmentLength     = 37.2*2*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadderWidth       = 50.0*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadderHeight      = 30.0*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadderSegBoxDW    =  3.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadderSegBoxDH    =  3.*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkLadderBeamRadius  =  0.6*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadderLa          =  3.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadderHa          =  0.6*fgkmm;     //total ???
const Double_t AliITSv11GeometrySDD::fgkLadderLb          =  3.7*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLadderHb          =  0.6*fgkmm;     //total ???
const Double_t AliITSv11GeometrySDD::fgkLadderl           =  0.25*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkBottomBeamAngle   = 56.5;
const Double_t AliITSv11GeometrySDD::fgkBeamSidePhi       = 65;

const Double_t AliITSv11GeometrySDD::fgkLadWaferSep       = 2*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkPinSuppWidth      = 2.5*fgkmm;       // ???
const Double_t AliITSv11GeometrySDD::fgkPinSuppHeight     = 2.*fgkmm;        // ???
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

const Double_t AliITSv11GeometrySDD::fgkLay3Rmin           = 145.*fgkmm;      // not min! Rmin virtual tube
const Double_t AliITSv11GeometrySDD::fgkLay3Rmax           = 205.*fgkmm;      // not min! Rmax virtual tube
const Double_t AliITSv11GeometrySDD::fgkLay3Length         = (524.+0.)*fgkmm; // ladder+supporting rings (length of the virtual tube)
const Double_t AliITSv11GeometrySDD::fgkLay3LadderLength   = 524.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLay3DetShortRadius = 146.0*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLay3DetLongRadius  = 152.0*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLay3LaddTopCornerEnd = 15.6*fgkmm;
const Int_t    AliITSv11GeometrySDD::fgkLay3Ndet           =  6;
const Int_t    AliITSv11GeometrySDD::fgkLay3Nladd          = 14;
const Double_t AliITSv11GeometrySDD::fgkLay3CoolPipeSuppH  =  7.5*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkLay4Rmin           = 220.*fgkmm;         // not min! Rmin virtual tube
const Double_t AliITSv11GeometrySDD::fgkLay4Rmax           = 290.*fgkmm;         // not min! Rmax virtual tube
const Double_t AliITSv11GeometrySDD::fgkLay4Length         = (671.+0.)*fgkmm;    // ladder+supporting rings (length of the virtual tube)
const Double_t AliITSv11GeometrySDD::fgkLay4LadderLength   = 671.*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLay4DetShortRadius = 235.0*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLay4DetLongRadius  = 240.5*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkLay4LaddTopCornerEnd = 15.6*fgkmm;
const Int_t    AliITSv11GeometrySDD::fgkLay4Ndet           = 8;
const Int_t    AliITSv11GeometrySDD::fgkLay4Nladd          = 22;
const Double_t AliITSv11GeometrySDD::fgkLay4CoolPipeSuppH  = 7.5*fgkmm;

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
const Double_t AliITSv11GeometrySDD::fgkmu = 1*fgkmicron; // can be increase for checking thin objects
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
const Double_t AliITSv11GeometrySDD::fgkWaferLengthSens    =  74.97*fgkmm;

const Double_t AliITSv11GeometrySDD::fgkDigitCablWidth     = 18.4*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkDigitCablAlThick   = (30+30*8./10.)*fgkmicron; // will probably change
const Double_t AliITSv11GeometrySDD::fgkDigitCablPolyThick = (20+12)*fgkmicron;        // will probably change

const Double_t AliITSv11GeometrySDD::fgkWaHVcableAlThick   = 30*2./10.*fgkmu;  // will probably change // Al ratio is random !!!
const Double_t AliITSv11GeometrySDD::fgkWaHVcablePolyThick = 175*fgkmu;        // will probably change
const Double_t AliITSv11GeometrySDD::fgkWaHVcableLength    = 67.08*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkWaHVcableWitdh     = 17.4 *fgkmm;              //  check !!!
const Double_t AliITSv11GeometrySDD::fgkWaHVcableDW        =  5.24*fgkmm;             //  check !!!

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
const Double_t AliITSv11GeometrySDD::fgkTransitHVtailXpos   =   8*fgkmm;              // ???, a mesurer !!!
const Double_t AliITSv11GeometrySDD::fgkTransitHVsideLZ     =  10.34*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkTransitHVsideLeftZ  =   4.11*fgkmm;
const Double_t AliITSv11GeometrySDD::fgkTransitHVsideRightZ =   3.5*fgkmm;           // ???, a mesurer !!!

const Double_t AliITSv11GeometrySDD::fgkLongHVcablePolyThick= (20+30+125+30+20+30+125+30+20)*fgkmu; //  check  // will probably change
const Double_t AliITSv11GeometrySDD::fgkLongHVcableAlThick  = (30+30*2/10+30)*fgkmu;                //  check  // will probably change
const Double_t AliITSv11GeometrySDD::fgkLongHVcableSeparation = 600*fgkmicron;


ClassImp(AliITSv11GeometrySDD)

//________________________________________________________________________
 AliITSv11GeometrySDD::AliITSv11GeometrySDD(): 
   AliITSv11Geometry(), fMotherVol(0), fAddHybrids(kTRUE), fAddSensors(kTRUE),
   fAddHVcables(kTRUE), fAddCables(kTRUE), fAddCoolingSyst(kTRUE),
   fCoolingOn(kTRUE),
   fAddOnlyLadder3min(-1), fAddOnlyLadder3max(-1),
   fAddOnlyLadder4min(-1), fAddOnlyLadder4max(-1)
{
  //
  // Standard constructor
  //

  fDigitCableLay3A = new AliITSv11GeomCableFlat[fgkLay3Ndet];
  fDigitCableLay3B = new AliITSv11GeomCableFlat[fgkLay3Ndet];
  fDigitCableLay4A = new AliITSv11GeomCableFlat[fgkLay4Ndet];
  fDigitCableLay4B = new AliITSv11GeomCableFlat[fgkLay4Ndet];
  SetParameters();
};


//________________________________________________________________________
AliITSv11GeometrySDD::AliITSv11GeometrySDD(Int_t debug) :
  AliITSv11Geometry(debug),fMotherVol(0),fAddHybrids(kTRUE),fAddSensors(kTRUE),
  fAddHVcables(kTRUE), fAddCables(kTRUE), fAddCoolingSyst(kTRUE),
  fCoolingOn(kTRUE),
  fAddOnlyLadder3min(-1), fAddOnlyLadder3max(-1),
  fAddOnlyLadder4min(-1), fAddOnlyLadder4max(-1)
{
  //
  // Constructor setting debugging level
  //

  fDigitCableLay3A = new AliITSv11GeomCableFlat[fgkLay3Ndet];
  fDigitCableLay3B = new AliITSv11GeomCableFlat[fgkLay3Ndet];
  fDigitCableLay4A = new AliITSv11GeomCableFlat[fgkLay4Ndet];
  fDigitCableLay4B = new AliITSv11GeomCableFlat[fgkLay4Ndet];
  SetParameters();
};

//________________________________________________________________________
AliITSv11GeometrySDD::AliITSv11GeometrySDD(const AliITSv11GeometrySDD &s) :
 AliITSv11Geometry(s.GetDebug()),fMotherVol(s.fMotherVol),
 fAddHybrids(s.fAddHybrids),fAddSensors(s.fAddSensors),
 fAddHVcables(s.fAddHVcables), fAddCables(s.fAddCables),
 fAddCoolingSyst(s.fAddCoolingSyst),fCoolingOn(s.fCoolingOn),
 fAddOnlyLadder3min(s.fAddOnlyLadder3min),fAddOnlyLadder3max(s.fAddOnlyLadder3max),
 fAddOnlyLadder4min(s.fAddOnlyLadder4min), fAddOnlyLadder4max(s.fAddOnlyLadder4max)
{
  //     Copy Constructor 
  fDigitCableLay3A = new AliITSv11GeomCableFlat[fgkLay3Ndet];
  fDigitCableLay3B = new AliITSv11GeomCableFlat[fgkLay3Ndet];
  fDigitCableLay4A = new AliITSv11GeomCableFlat[fgkLay4Ndet];
  fDigitCableLay4B = new AliITSv11GeomCableFlat[fgkLay4Ndet];
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
};

//________________________________________________________________________
void AliITSv11GeometrySDD::SetParameters() {
  //
  // Define display colors and the non constant geometry parameters
  //

  fColorCarbonFiber = 4;
  fColorRyton = 5;
  fColorPhynox = 7;
  fColorSilicon = 3;
  fColorAl = 7;
  fColorPolyhamide = 5;
  fColorGlass = 2;
  fColorSMD = 12;
  fColorSMDweld = 17;

  fPinSupport = 0;
  fCoolPipeSupportL = 0;
  fCoolPipeSupportR = 0;
  fSDDsensor = 0;
  fBaseThermalBridge = 0;
  fHybrid = 0;

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
};


//________________________________________________________________________
TGeoMedium* AliITSv11GeometrySDD::GetMedium(const char* mediumName) {
  //
  // Called to get a medium, checks that it exists.
  // If not, prints an error and returns 0
  //

  TGeoMedium* medium =  gGeoManager->GetMedium(mediumName);
  if (! medium)
    printf("Error(AliITSv11GeometrySDD)::medium %s not found !\n", mediumName);

  return medium;
};

//________________________________________________________________________
void AliITSv11GeometrySDD::CreateBasicObjects() {
  //
  // Create basics objets which will be assembled together
  // in Layer3 and Layer4 functions
  //

  fPinSupport = CreatePinSupport();
  fCoolPipeSupportL = CreateCoolPipeSupportL();
  fCoolPipeSupportR = CreateCoolPipeSupportR();
  fSDDsensor = CreateSDDsensor();
  fBaseThermalBridge = CreateBaseThermalBridge();
  fHybrid = CreateHybrid(0);

  TGeoMedium *carbonFiberLadderStruct = GetMedium("ITSsddCarbonM55J");
  TGeoMedium *polyhamideSDD = GetMedium("ITSsddKAPTON_POLYCH2");
  TGeoMedium *alSDD         = GetMedium("ITSal");

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
  TGeoArb8 *cfLaddTop1 = CreateLadderSide( fgkSegmentLength/2., halfTheta, 
			  -1, fgkLadderLa, fgkLadderHa, fgkLadderl);
  TGeoVolume *cfLaddTopVol1 = new TGeoVolume("ITSsddCFladdTopCornerVol1",
				  cfLaddTop1,carbonFiberLadderStruct);
  TGeoArb8 *cfLaddTop2 = CreateLadderSide( fgkSegmentLength/2., halfTheta,
			   1, fgkLadderLa, fgkLadderHa, fgkLadderl);
  TGeoVolume *cfLaddTopVol2 = new TGeoVolume("ITSsddCFladdTopCornerVol2",
				  cfLaddTop2, carbonFiberLadderStruct);
  cfLaddTopVol1->SetLineColor(fColorCarbonFiber);
  cfLaddTopVol2->SetLineColor(fColorCarbonFiber);
  TGeoTranslation *trTop1 = new TGeoTranslation(0, fgkLadderHeight/2-dy, 0);

  //--- the 2 side V
  TGeoArb8 *cfLaddSide1 = CreateLadderSide( fgkSegmentLength/2., beta, -1,
					    fgkLadderLb, fgkLadderHb, fgkLadderl);
  TGeoVolume *cfLaddSideVol1 = new TGeoVolume( "ITSsddCFladdSideCornerVol1",
				   cfLaddSide1,carbonFiberLadderStruct);
  TGeoArb8 *cfLaddSide2 = CreateLadderSide( fgkSegmentLength/2., beta, 1,
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
  TGeoRotation beamRot1("", 90-2.*beta*TMath::RadToDeg(),
			-beamPhiPrime*TMath::RadToDeg(),-90);
  TGeoRotation beamRot2("", 90-2.*beta*TMath::RadToDeg(),
			beamPhiPrime*TMath::RadToDeg(), -90);
  TGeoRotation beamRot3("", 90+2.*beta*TMath::RadToDeg(),
			beamPhiPrime*TMath::RadToDeg(), -90);
  TGeoRotation beamRot4("", 90+2.*beta*TMath::RadToDeg(),
			-beamPhiPrime*TMath::RadToDeg(),-90);

  TGeoCombiTrans *beamTransf[8];
  beamTransf[0] = new TGeoCombiTrans( 0.5*triangleHeight*
				      TMath::Tan(halfTheta),
				      fgkLadderBeamRadius/2. - dy,
				      -3*fgkSegmentLength/8, &beamRot1);
  beamTransf[1] = new TGeoCombiTrans(*beamTransf[0]);
  AddTranslationToCombiTrans(beamTransf[1], 0, 0, fgkSegmentLength/2);

  beamTransf[2] = new TGeoCombiTrans(0.5*triangleHeight*
				     TMath::Tan(halfTheta),
				     fgkLadderBeamRadius/2. - dy,
				     -fgkSegmentLength/8, &beamRot2);
  beamTransf[3] = new TGeoCombiTrans(*beamTransf[2]);
  AddTranslationToCombiTrans(beamTransf[3], 0, 0, fgkSegmentLength/2);

  beamTransf[4] = new TGeoCombiTrans(-0.5*triangleHeight*
				     TMath::Tan(halfTheta),
				     fgkLadderBeamRadius/2. - dy,
				     -3*fgkSegmentLength/8, &beamRot3);
  beamTransf[5] = new TGeoCombiTrans(*beamTransf[4]);
  AddTranslationToCombiTrans(beamTransf[5], 0, 0, fgkSegmentLength/2);

  beamTransf[6] = new TGeoCombiTrans(-0.5*triangleHeight*
  TMath::Tan(halfTheta),fgkLadderBeamRadius/2.-dy, -fgkSegmentLength/8,&beamRot4);
  beamTransf[7] = new TGeoCombiTrans(-0.5*triangleHeight*
  TMath::Tan(halfTheta),fgkLadderBeamRadius/2.-dy,3*fgkSegmentLength/8,&beamRot4);

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

  TGeoRotation bottomBeamRot1("", 90, 90,  90);
  TGeoRotation bottomBeamRot2("",-90, 90, -90);
  TGeoCombiTrans *bottomBeamTransf1 = new TGeoCombiTrans
    (0,-(fgkLadderHeight/2-fgkLadderBeamRadius)-dy,0, &bottomBeamRot1);
  TGeoCombiTrans *bottomBeamTransf2 = new TGeoCombiTrans(0,
				      -(fgkLadderHeight/2-fgkLadderBeamRadius)-dy,
				      -fgkSegmentLength/2, &bottomBeamRot1);
  TGeoCombiTrans *bottomBeamTransf3 = new TGeoCombiTrans(0,
                                      -(fgkLadderHeight/2 - fgkLadderBeamRadius)
				      - dy, fgkSegmentLength/2, &bottomBeamRot2);
  // be careful for beams #3: when "reading" from -z to +z and 
  // from the bottom of the ladder, it should draw a Lambda, and not a V
  TGeoRotation bottomBeamRot4("", -90, fgkBottomBeamAngle, -90);
  TGeoRotation bottomBeamRot5("" ,-90,-fgkBottomBeamAngle, -90);
  TGeoCombiTrans *bottomBeamTransf4 = new TGeoCombiTrans
    (0,-(fgkLadderHeight/2-fgkLadderBeamRadius)-dy,-fgkSegmentLength/4,&bottomBeamRot4);
  TGeoCombiTrans *bottomBeamTransf5 = new TGeoCombiTrans
    (0,-(fgkLadderHeight/2-fgkLadderBeamRadius)-dy,fgkSegmentLength/4, &bottomBeamRot5);

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

};


//________________________________________________________________________
void AliITSv11GeometrySDD::CheckOverlaps(Double_t precision){
  //
  // a debugging function for checking some possible overlaps
  //
  if (fSDDsensor)         fSDDsensor->CheckOverlaps(precision);
  if (fHybrid)            fHybrid->CheckOverlaps(precision);
};


//________________________________________________________________________
TGeoCombiTrans *AliITSv11GeometrySDD::
CreateCombiTrans(const char *name, Double_t dy, Double_t dz, Double_t dphi) {
    //
    // return the TGeoCombiTrans which make a translation in y and z
    // and a rotation in phi in the global coord system
    //

    TGeoTranslation t1(dy*CosD(90.+dphi),dy*SinD(90.+dphi), dz);
    TGeoRotation r1("",0.,0.,dphi);

    TGeoCombiTrans *combiTrans1 = new TGeoCombiTrans(name);
    combiTrans1->SetTranslation(t1);
    combiTrans1->SetRotation(r1);
    return combiTrans1;
};


//________________________________________________________________________
void AliITSv11GeometrySDD::AddTranslationToCombiTrans(TGeoCombiTrans* ct,
                                                       Double_t dx,
                                                       Double_t dy,
                                                       Double_t dz) const{
  // Add a dx,dy,dz translation to the initial TGeoCombiTrans
  const Double_t *vect = ct->GetTranslation();
  Double_t newVect[3] = {vect[0]+dx, vect[1]+dy, vect[2]+dz};
  ct->SetTranslation(newVect);
};


//________________________________________________________________________
void AliITSv11GeometrySDD::ShowOnePiece(TGeoVolume *moth) {
// for code developpment and debugging purposes

 if (! fSDDsensor) CreateBasicObjects();

//     Moth->AddNode(fBaseThermalBridge, 1, 0);
     moth->AddNode(fHybrid,100,0);
//     moth->AddNode(fSDDsensor, 1, 0);

//   TGeoVolume* seg = CreateLadderSegment( 4, 0); //lay 4
//   moth->AddNode(seg, 1, 0);




//  TGeoBBox *box1 = new TGeoBBox("box1", 5,5,5);
//  TGeoMedium *air = GetMedium("ITSair");
//  TGeoVolume *vbox1 = new TGeoVolume("vbox1", box1, air);
//  TGeoBBox *box2 = new TGeoBBox("box2", 6,6,6);
//  TGeoVolume *vbox2 = new TGeoVolume("vbox2", box2, air);
//  TGeoBBox *box3 = new TGeoBBox("box3", 7,7,7);
//  TGeoVolume *vbox3 = new TGeoVolume("vbox3", box3, air);

//  vbox1->AddNode(fHybrid,100,0);
//  vbox2->AddNode(vbox1,1,0);
//  vbox3->AddNode(vbox2,1,0);
//  moth->AddNode(vbox3,1,0);



//   //testing cable
//   TGeoBBox *box1 = new TGeoBBox("box1", 10,10,10);
//   TGeoBBox *box2 = new TGeoBBox("box2", 10,10,10);
//   TGeoBBox *box3 = new TGeoBBox("box3", 20,10,10);
//   TGeoMedium *air = GetMedium("ITSsddAir");
//   TGeoVolume *vbox1 = new TGeoVolume("vbox1", box1, air);
//   TGeoVolume *vbox2 = new TGeoVolume("vbox2", box2, air);
//   TGeoVolume *vbox3 = new TGeoVolume("vbox3", box3, air);

//   TGeoTranslation *tr1 = new TGeoTranslation("merdeneg",-10,0,0);
//   TGeoTranslation *tr2 = new TGeoTranslation("merdepos",10,0,0);

//   AliITSv11GeomCableRound napCable(0.9);
//   //AliITSv11GeomCableFlat napCable(2,0.9);
//   napCable.SetNLayers(3);
//   napCable.SetLayer(0, 0.2, air);
//   napCable.SetLayer(1, 0.2, air);
//   napCable.SetLayer(2, 0.5, air);

//   napCable.SetInitialNode(vbox3);

//   Double_t coord1[3] = {0,-2,-2};
//   Double_t vect1[3]= {1,1,0};
//   napCable.AddCheckPoint( vbox1, 0, coord1, vect1);
//   Double_t coord2[3] = {10,0,0};
//   Double_t vect2[3]= {1,0,0};
//   napCable.AddCheckPoint( vbox1, 1, coord2, vect2);

//   //Double_t coord3[3] = {7,7,7};
//   Double_t coord3[3] = {7,-7,-7};
//   Double_t vect3[3]= {1,0,0};
//   napCable.AddCheckPoint( vbox3, 2, coord3, vect3);

//   Double_t coord4[3] = {19,7,7};
//   Double_t vect4[3]= {-1,0,2};
//   napCable.AddCheckPoint( vbox3, 3, coord4, vect4);

//   Double_t coord5[3] = {1,7,7};
//   Double_t vect5[3]= {1,0,0};
//   napCable.AddCheckPoint( vbox3, 4, coord5, vect5);

 
//   TGeoRotation *rot = new TGeoRotation("",0,0,0);
//   TGeoCombiTrans *combi = new TGeoCombiTrans(*tr1,*rot );
//   //vbox3->AddNode(vbox1,1,tr1);
//   vbox3->AddNode(vbox1,1,combi);
//   moth->AddNode(vbox3,1,0);

// //   napCable.CreateAndInsertCableSegment( 1, 135);
// //   napCable.CreateAndInsertCableSegment( 2, 0);
// //   napCable.CreateAndInsertCableSegment( 3, 0);
// //   napCable.CreateAndInsertCableSegment( 4, 0);
//   napCable.CreateAndInsertCableSegment( 1);
//   napCable.CreateAndInsertCableSegment( 2);
//   napCable.CreateAndInsertCableSegment( 3);
//   napCable.CreateAndInsertCableSegment( 4);
//   napCable.PrintCheckPoints();
};


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

  TGeoMedium *airSDD = GetMedium("ITSair");

  fMotherVol = moth;
  if (! fSDDsensor) CreateBasicObjects();
  
  TGeoVolume *lay3Ladder    = CreateLadder(3);
  TGeoVolume *lay3Detectors = CreateDetectors(3);
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

    sprintf(rotName, "ITSsddLay3Ladd%i",iLadd);
    Double_t minRadiusLadBox = fLay3LaddShortRadius-fLay3LadderUnderSegDH;
    if (iLadd%2 != 0)
      minRadiusLadBox = fLay3LaddLongRadius-fLay3LadderUnderSegDH;
    minRadiusLadBox += ((TGeoBBox*)lay3Ladder->GetShape())->GetDY();
    TGeoCombiTrans *ctLadd = CreateCombiTrans(rotName,minRadiusLadBox,
					      0,-90+iLadd*dPhi);
    virtualLayer3->AddNode(lay3Ladder, iLadd, ctLadd);
    ///////////////////////////////////////////////////
    sprintf(rotName, "ITSsddLay3DetBox%i",iLadd);
    Double_t minRadiusDetBox = fgkLay3DetShortRadius;
    if (iLadd%2 != 0) minRadiusDetBox = fgkLay3DetLongRadius;
    minRadiusDetBox += detectorsThick/2;
    TGeoCombiTrans *ctDet = CreateCombiTrans(rotName, minRadiusDetBox,
					     0,-90+iLadd*dPhi);
    virtualLayer3->AddNode(lay3Detectors, iLadd, ctDet);
    ///////////////////////////////////////////////////
  }

  if(GetDebug(1)) virtualLayer3->CheckOverlaps(0.01);
  virtualLayer3->SetVisibility(kFALSE);
  moth->AddNode(virtualLayer3, 1, 0);
};


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

  if (! fSDDsensor) CreateBasicObjects();

  TGeoTube *virtualLayer4Shape =new TGeoTube("ITSsddLayer4Shape",
				    fgkLay4Rmin,fgkLay4Rmax,fgkLay4Length*0.5);
  TGeoMedium *airSDD = GetMedium("ITSair");
  TGeoVolume *virtualLayer4 = new TGeoVolume("ITSsddLayer4",
					     virtualLayer4Shape, airSDD);
  TGeoVolume *lay4Ladder    = CreateLadder(4);
  TGeoVolume *lay4Detectors = CreateDetectors(4);
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
    sprintf(rotName, "ITSsddLay4Ladd%i",iLadd);
    Double_t minRadiusLadBox = fLay4LaddShortRadius-fLay4LadderUnderSegDH;
    if (iLadd%2 != 0)
      minRadiusLadBox = fLay4LaddLongRadius-fLay4LadderUnderSegDH;
    minRadiusLadBox += ((TGeoBBox*)lay4Ladder->GetShape())->GetDY();
    TGeoCombiTrans *ctLadd = CreateCombiTrans(rotName, minRadiusLadBox,
					      0, -90+iLadd*dPhi);
    virtualLayer4->AddNode(lay4Ladder, iLadd, ctLadd);
    sprintf(rotName, "ITSsddLay4DetBox%i",iLadd);
    Double_t minRadiusDetBox = fgkLay4DetShortRadius;
    if (iLadd%2 != 0)
      minRadiusDetBox = fgkLay4DetLongRadius;
    minRadiusDetBox += detBoxThickness/2;
    TGeoCombiTrans *ctDet = CreateCombiTrans(rotName, minRadiusDetBox,
					     0, -90+iLadd*dPhi);
    virtualLayer4->AddNode(lay4Detectors, iLadd, ctDet);
  }

  if(GetDebug(1)) virtualLayer4->CheckOverlaps(0.01);
  virtualLayer4->SetVisibility(kFALSE);
  moth->AddNode(virtualLayer4,1,0);
};


//________________________________________________________________________
TGeoVolume *AliITSv11GeometrySDD::CreateLadder(Int_t iLay) {
  //
  // return a box volume containing the CF ladder
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
  TGeoBBox *ladBox = new TGeoBBox("ITSsddLadBox",
			 fgkLadderWidth/2+fgkPinSuppWidth+fgkLadderSegBoxDW,
			 ladderBoxDH/2, ladderLength/2);
  TGeoMedium *airSDD = GetMedium("ITSair");
  TGeoVolume *virtualLadder = new TGeoVolume("ITSsddLadder",ladBox, airSDD);
  
  // placing virtual ladder segment following detector ordering convention
  //=======================================================================
  char transName[30];

  // adding segment this way to create cable points in the correct order ...
  for (Int_t iSegment = nDetectors/2-1; iSegment >= 0; iSegment-- ) {

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
  TGeoVolume *endLadder = CreateEndLadder( iLay );
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
      if (iPt>1) rotation = 90-fgkHybridAngle; 
      digitCableA[iSegment].CreateAndInsertCableSegment(iPt, rotation);
    };
    
    for (Int_t iPt=1; iPt<digitCableB[iSegment].GetNCheckPoints(); iPt++ ) {
      Double_t rotation = 0;
      if (iPt>1) rotation = fgkHybridAngle-90; 
      digitCableB[iSegment].CreateAndInsertCableSegment(iPt, rotation);
    };
  };
  
  // HV cable
  //=======================================================================
  TGeoMedium *polyhamideSDD = GetMedium("ITSsddKAPTON_POLYCH2");
  TGeoMedium *alSDD         = GetMedium("ITSal");

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

  x1[0] = fgkTransitHVtailXpos;
  x2[0] = fgkTransitHVtailXpos;
  x3[0] = fgkTransitHVtailXpos;
  for (Int_t iSegment  = nDetectors/2-1; iSegment >= 0; iSegment-- ) {
    Double_t cableSeparation = TMath::Abs(iSegment - (nDetectors/2-1))
                               *fgkLongHVcableSeparation;
    x1[1] = - ladderBoxDH/2;
    x2[1] = - ladderBoxDH/2 + underSegDH - cableSeparation
            - (fgkLongHVcablePolyThick+fgkLongHVcableAlThick)/2;
    x3[1] = x2[1];
    x1[2] = sensorZPos[iSegment]+fgkTransitHVtailLength-5*fgkmm;
    x2[2] =  x1[2]+5*fgkmm;
    x3[2] = ladderLength/2-endLength;
    cableHV[iSegment].AddCheckPoint( virtualLadder, 0, x1, vY );
    cableHV[iSegment].AddCheckPoint( virtualLadder, 1, x2, vYZ );
    cableHV[iSegment].AddCheckPoint( virtualLadder, 2, x3, vZ );

    cableHV[iSegment].CreateAndInsertCableSegment(1,0);
    cableHV[iSegment].CreateAndInsertCableSegment(2,0);
  };

  vYZ[2] = -1;
  for (Int_t iSegment = nDetectors/2; iSegment < nDetectors; iSegment++ ) { 
    Double_t cableSeparation = TMath::Abs(iSegment - (nDetectors/2-1))
                               *fgkLongHVcableSeparation;
    x1[1] = - ladderBoxDH/2;
    x2[1] = - ladderBoxDH/2 + underSegDH - cableSeparation
            - (fgkLongHVcablePolyThick+fgkLongHVcableAlThick)/2;
    x3[1] = x2[1];
    x1[2] = sensorZPos[iSegment]-fgkTransitHVtailLength+5*fgkmm;
    x2[2] =  x1[2]-5*fgkmm;
    x3[2] = -ladderLength/2+endLength;
    cableHV[iSegment].AddCheckPoint( virtualLadder, 0, x1, vY );
    cableHV[iSegment].AddCheckPoint( virtualLadder, 1, x2, vYZ );
    cableHV[iSegment].AddCheckPoint( virtualLadder, 2, x3, vZ );

    cableHV[iSegment].CreateAndInsertCableSegment(1,0);
    cableHV[iSegment].CreateAndInsertCableSegment(2,0);
  };

  //**********************************
  if(GetDebug(1)) virtualLadder->CheckOverlaps(0.01);
  //virtualLadder->SetVisibility(kFALSE);
  return virtualLadder;
};


//________________________________________________________________________
TGeoArb8 *AliITSv11GeometrySDD::CreateLadderSide(Double_t dz, Double_t angle,
                         Double_t xSign, Double_t L, Double_t H, Double_t l) {
    // Create one half of the V shape corner of CF ladder
  
    TGeoArb8 *cfLaddSide = new TGeoArb8(dz);
    cfLaddSide->SetVertex( 0, 0,  0);
    cfLaddSide->SetVertex( 1, 0, -H);
    cfLaddSide->SetVertex( 2, xSign*(L*TMath::Sin(angle)-l*TMath::Cos(angle)),
			   -L*TMath::Cos(angle)-l*TMath::Sin(angle));
    cfLaddSide->SetVertex( 3, xSign*L*TMath::Sin(angle), -L*TMath::Cos(angle));
    cfLaddSide->SetVertex( 4, 0,  0);
    cfLaddSide->SetVertex( 5, 0, -H);
    cfLaddSide->SetVertex( 6, xSign*(L*TMath::Sin(angle)-l*TMath::Cos(angle)),
			   -L*TMath::Cos(angle)-l*TMath::Sin(angle));
    cfLaddSide->SetVertex(7, xSign*L*TMath::Sin(angle), -L*TMath::Cos(angle));
    return cfLaddSide;
};


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
  TGeoMedium *airSDD                  = GetMedium("ITSair");
  TGeoMedium *carbonFiberLadderStruct = GetMedium("ITSsddCarbonM55J");
  TGeoMedium *alSDD                   = GetMedium("ITSal");
  TGeoMedium *alSDD80p100             = GetMedium("ITSal");                 // to code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoMedium *alSDD50p100             = GetMedium("ITSal");                 // to code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoMedium *polyhamideSDD           = GetMedium("ITSsddKAPTON_POLYCH2");
  TGeoMedium *niSDD                   = GetMedium("COPPER");                // to code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoMedium *glueAG                  = GetMedium("ITSsddKAPTON_POLYCH2");  // to code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoMedium *siliconSDD              = GetMedium("ITSsddSiChip");
  TGeoMedium *medSMD                  = GetMedium("SDDX7Rcapacitors");      //     TO CODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoMedium *medSMDweld              = GetMedium("SDDX7Rcapacitors");      //    TO CODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  TGeoRotation rotHole("", 0, 90, 0);
  TGeoCombiTrans *hybHolePos1 = new TGeoCombiTrans(hybHolePos1Tr, rotHole);
  hybHolePos1->SetName("hybHolePos1");
  TGeoCombiTrans *hybHolePos2 = new TGeoCombiTrans(hybHolePos2Tr, rotHole);
  hybHolePos2->SetName("hybHolePos2");

  upGlueScreenTr->RegisterYourself();
  alScreenTr->RegisterYourself();
  hybHolePos1->RegisterYourself();
  hybHolePos2->RegisterYourself();

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

  if(GetDebug(3)){ // Remove compiler warning.
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
  lowFLpiece.CreateAndInsertCableSegment(1);
  lowFLpiece.ResetPoints();

  Double_t piece2width = fgkHybFLlowAmbX-fgkHybFLlowPasX
                         -fgkHybFLlowHolePasDX/2-fgkHybFLlowHoleAmbDX/2;

  lowFLpiece.SetWidth(piece2width);
  lowFLpiece.SetName("lowFLpiece2");
  x1[0] = piece2width/2+fgkHybFLlowPasX+fgkHybFLlowHolePasDX/2-fgkHybridWidth/2;
  x2[0] = x1[0];
  lowFLpiece.AddCheckPoint( hybrid, 0, x2, vZ );
  lowFLpiece.AddCheckPoint( hybrid, 1, x1, vZ );
  lowFLpiece.CreateAndInsertCableSegment(1);
  lowFLpiece.ResetPoints();

  Double_t piece3width = fgkHybridWidth - fgkHybFLlowAmbX 
                         - fgkHybFLlowHoleAmbDX/2;

  lowFLpiece.SetWidth(piece3width);
  lowFLpiece.SetName("lowFLpiece3");
  x1[0] = fgkHybridWidth/2-piece3width/2;
  x2[0] = x1[0];
  lowFLpiece.AddCheckPoint( hybrid, 0, x2, vZ );
  lowFLpiece.AddCheckPoint( hybrid, 1, x1, vZ );
  lowFLpiece.CreateAndInsertCableSegment(1);

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
    lowFLpiece.CreateAndInsertCableSegment(1,90);
    lowFLpiece.ResetPoints();

    sprintf(ch, "lowFLpieceB%i", i+4);
    lowFLpiece.SetName(ch);
    x1[0] = fgkHybridWidth/2 - piece3width;
    x2[0] = x1[0] - fgkHybFLlowHoleAmbDX;
    lowFLpiece.AddCheckPoint( hybrid, 0, x1, vX );
    lowFLpiece.AddCheckPoint( hybrid, 1, x2, vX );
    lowFLpiece.CreateAndInsertCableSegment(1,90);
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
    chip.CreateAndInsertCableSegment(1,-90);
    chip.ResetPoints();

    sprintf(ch, "ambraCC%i", i);
    chip.SetName(ch);
    x1[0] = fgkHybFLlowAmbX - fgkHybridWidth/2 - fgkHybAmbraDX/2;
    x2[0] = x1[0] + fgkHybAmbraDX;
    chip.AddCheckPoint( hybrid, 0, x1, vX );
    chip.AddCheckPoint( hybrid, 1, x2, vX );
    chip.CreateAndInsertCableSegment(1,-90);
    chip.ResetPoints();
  };

  //**************************************************** CC outside chips:
  // je crois qu'il n'y a pas de 2ieme couche d'alu ici ... 
  for (Int_t i = 0; i<4; i++) {
    char ch[20];
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
    ccLayer1.CreateAndInsertCableSegment(1,-90);

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
    ccLayer2.CreateAndInsertCableSegment(1,-90);
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
    ccLayer2.CreateAndInsertCableSegment(1,-90);
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
  TGeoRotation rotAluStrip("rotAluStrip",0, -90, 90);
  Double_t yRotAluStrip = lowLayerYmin+lowFLTotalThick
                          +flUpThick+fgkHybAlThick/2;
  TGeoCombiTrans *aluStripTr1 = new TGeoCombiTrans(
				    fgkHybridWidth/2,yRotAluStrip,
		    fgkHybridLength/2-fgkHybFLlowChipZ1+1*fgkmm, &rotAluStrip);
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
};


//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateLadderSegment(Int_t iLay, Int_t iSeg) {
  //
  // Return a box volume containing a segment of a ladder.
  //

  TGeoMedium *airSDD          = GetMedium("ITSair");
  TGeoMedium *phynoxSDD       = GetMedium("ITSal"); // phynoxSDD To code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoMedium *coolerMediumSDD = GetMedium("WATER");

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
  TGeoBBox *segBox = new TGeoBBox("ITSsddSegBox",
				  fgkLadderWidth/2+fgkPinSuppWidth+fgkLadderSegBoxDW,
				  fgkLadderHeight/2+fgkLadderSegBoxDH/2,
				  segmentLength/2);
  
  TGeoVolume *virtualSeg = new TGeoVolume("ITSsddSegment",
					  segBox, airSDD);

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
  TGeoRotation rotPS1("",0,-90,90);
  TGeoRotation rotPS2("",0,-90,-90);
  TGeoCombiTrans *transPS1 = new TGeoCombiTrans( fgkPinDYOnSensor,
				- fgkLadderHeight/2.-tDY
			        + fgkPinSuppHeight/2.,
				sensorCenterZPos+fgkPinDXminOnSensor,&rotPS1);
  TGeoCombiTrans *transPS2 = new TGeoCombiTrans(*transPS1);
  AddTranslationToCombiTrans(transPS2, 0, 0, fgkPinPinDDXOnSensor);
  TGeoCombiTrans *transPS3 = new TGeoCombiTrans(*transPS1);
  AddTranslationToCombiTrans(transPS3, 0, 0, -2*fgkPinDXminOnSensor);
  TGeoCombiTrans *transPS4 = new TGeoCombiTrans(*transPS3);
  AddTranslationToCombiTrans(transPS4, 0, 0, -fgkPinPinDDXOnSensor);

  TGeoCombiTrans *transPS5 = new TGeoCombiTrans( -fgkPinDYOnSensor,
		       	         - fgkLadderHeight/2. - tDY
			         + fgkPinSuppHeight/2.,
			         sensorCenterZPos+fgkPinDXminOnSensor,&rotPS2);
  TGeoCombiTrans *transPS6 = new TGeoCombiTrans(*transPS5);
  AddTranslationToCombiTrans(transPS6, 0, 0, fgkPinPinDDXOnSensor);
  TGeoCombiTrans *transPS7 = new TGeoCombiTrans(*transPS5);
  AddTranslationToCombiTrans(transPS7, 0, 0, -2*fgkPinDXminOnSensor);
  TGeoCombiTrans *transPS8 = new TGeoCombiTrans(*transPS7);
  AddTranslationToCombiTrans(transPS8, 0, 0, -fgkPinPinDDXOnSensor);
  
  virtualSeg->AddNode(fPinSupport, 1, transPS1);
  virtualSeg->AddNode(fPinSupport, 2, transPS2);
  virtualSeg->AddNode(fPinSupport, 3, transPS3);
  virtualSeg->AddNode(fPinSupport, 4, transPS4);
  virtualSeg->AddNode(fPinSupport, 5, transPS5);
  virtualSeg->AddNode(fPinSupport, 6, transPS6);
  virtualSeg->AddNode(fPinSupport, 7, transPS7);
  virtualSeg->AddNode(fPinSupport, 8, transPS8);
  
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
  TGeoRotation rotCPS2("", -halfTheta*TMath::RadToDeg(), -90,  90);
  TGeoRotation rotCPS1("",  halfTheta*TMath::RadToDeg(), -90, -90);
  TGeoCombiTrans *transCPS1 = new TGeoCombiTrans(coolPipeSuppL,
				  -fgkLadderHeight/2.-TMath::Abs(tDY)
				  +coolPipeSuppH+fgkLadderBeamRadius,
			          -segmentLength/2., &rotCPS1);
  TGeoCombiTrans *transCPS3 = new TGeoCombiTrans(*transCPS1);
  AddTranslationToCombiTrans(transCPS3, 0, 0, segmentLength);
  
  TGeoCombiTrans *transCPS2 = new TGeoCombiTrans(-coolPipeSuppL,
				  -fgkLadderHeight/2.-tDY
				  +coolPipeSuppH+fgkLadderBeamRadius,
			          segmentLength/2., &rotCPS2);
  TGeoCombiTrans *transCPS4 = new TGeoCombiTrans(*transCPS2);
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
  //virtualSeg->SetVisibility(kFALSE);
  return virtualSeg;
};


//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreatePinSupport() {
//
// Create a pine support
// axis of rotation is the cone axis, center in its middle
//
    TGeoCone *cone = new TGeoCone("ITSsddPinSuppCone",fgkPinSuppHeight/2.,
                                  0,fgkPinSuppRmax,0,fgkPinSuppRmax-
                                  fgkPinSuppHeight*TanD(fgkPinSuppConeAngle) );
    TGeoBBox *tong = new TGeoBBox("ITSsddPinSuppTong",fgkPinSuppRmax,
                                  fgkPinSuppLength/2.,fgkPinSuppThickness/2.);
    TGeoTube *hole = new TGeoTube("ITSsddPinSuppHole",0,fgkPinR,
                                  fgkPinSuppHeight/2.);
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

    
    TGeoMedium *rytonSDD = GetMedium("ITSsddCarbonM55J"); //medium = ryton ?  To code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TGeoVolume *pinSupport = new TGeoVolume("ITSsddPinSupport",pinSupportShape,
                                            rytonSDD);
    pinSupport->SetLineColor(fColorRyton);
    return pinSupport;
    // include the pin itself                           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
};


//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateCoolPipeSupportL() {
//
// Create half of the cooling pipe support (ALR-0752/3)
//

  Double_t diffX = fgkCoolPipeSuppHeight*TanD(fgkCoolPipeSuppAngle);
  
  TGeoArb8 *side1 = new TGeoArb8(fgkCoolPipeSuppHeight/2.);
  side1->SetVertex( 0, 0, -fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 1, fgkCoolPipeSuppMaxLength/2.-diffX,
		       -fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 2, fgkCoolPipeSuppMaxLength/2.-diffX,
		       fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 3, 0,  fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 4, 0, -fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 5, fgkCoolPipeSuppMaxLength/2.,
		       -fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 6, fgkCoolPipeSuppMaxLength/2.,
		       fgkCoolPipeSuppWidthExt/2.);
  side1->SetVertex( 7, 0,  fgkCoolPipeSuppWidthExt/2.);
  side1->SetName("ITSsddCPSside1");

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
  TGeoRotation axeRot("ITSsddCPSaxeRot",90,90,0);

  TGeoCombiTrans *axeTrans = new TGeoCombiTrans("ITSsddCPSaxeTr",
				 fgkCoolPipeSuppTongW/4.,0,0,&axeRot);
  axeTrans->RegisterYourself();

  if(GetDebug(3)){
    middle->InspectShape();
    axe->InspectShape();
  };

  TGeoMedium *rytonSDD = GetMedium("ITSsddCarbonM55J"); //medium = ryton ?  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
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
};

//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateCoolPipeSupportR() {
//
//Create half of the cooling pipe support (ALR-0752/3)
//

  Double_t diffX = fgkCoolPipeSuppHeight*TanD(fgkCoolPipeSuppAngle);
  
  TGeoArb8 *side1 = new TGeoArb8(fgkCoolPipeSuppHeight/2.);
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
  side1->SetName("ITSsddCPSside1R");

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
  TGeoRotation axeRot("ITSsddCPSaxeRotR",90,90,0);
  TGeoCombiTrans *axeTrans = new TGeoCombiTrans("ITSsddCPSaxeTrR",
				 -fgkCoolPipeSuppTongW/4.,0,0,&axeRot);
  axeTrans->RegisterYourself();

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
  
  TGeoMedium *rytonSDD = GetMedium("ITSsddCarbonM55J"); //medium = ryton ? To code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoVolume *coolPipeSupp = new TGeoVolume( "ITSsddCoolPipeSupportR",
					     coolPipeSuppShape, rytonSDD);
  coolPipeSupp->SetLineColor(fColorRyton);

  return coolPipeSupp;
};

//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateBaseThermalBridge() {
// ALR 0752/8

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

  TGeoMedium *carbonFiberLadderStruct = GetMedium("ITSsddCarbonM55J");
  TGeoVolume *vBaseThermalBridge = new TGeoVolume( "ITSsddBaseThermalBridge",
						   sBaseThermalBridge,
						   carbonFiberLadderStruct);

  vBaseThermalBridge->SetLineColor(fColorCarbonFiber);
  return vBaseThermalBridge;
};


//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateEndLadder(Int_t iLay) {
  //
  // Return a box volume containing a end of a CF ladder.
  //

  TGeoMedium *airSDD                  = GetMedium("ITSair");
  TGeoMedium *carbonFiberLadderStruct = GetMedium("ITSsddCarbonM55J");

  Double_t length        = (fgkLay3LadderLength-fgkLay3Ndet*fgkSegmentLength)/2.;
  Double_t coolPipeSuppH = fgkLay3CoolPipeSuppH;
  Double_t underSegDH    = fLay3LadderUnderSegDH;
  if (iLay==3) {
  } else if (iLay==4) {
    length        = (fgkLay4LadderLength-fgkLay4Ndet*fgkSegmentLength)/2.;
    coolPipeSuppH = fgkLay4CoolPipeSuppH;
    underSegDH    = fLay4LadderUnderSegDH;
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

  TGeoBBox *endBox = new TGeoBBox("ITSsddEndLaddBox",
			 (fgkLadderWidth)/2,
			 fgkLadderHeight/2+fgkLadderSegBoxDH/2+underSegDH/2,
		         length/2);
  TGeoVolume *virtualEnd = new TGeoVolume("ITSsddEnd",endBox, airSDD);
  
  //**********************************
  // coding real matter :
  //**********************************
  Double_t triangleHeight = fgkLadderHeight - fgkLadderBeamRadius;
  Double_t halfTheta   = TMath::ATan( 0.5*fgkLadderWidth/triangleHeight );
  Double_t beta        = (TMath::Pi()-2.*halfTheta)/4.;
  Double_t alpha       = TMath::Pi()*3./4. - halfTheta/2.;
  
  //--- The 3 V shape corners of the Carbon Fiber Ladder
  //--- the top V
  TGeoArb8 *cfLaddTop1 = CreateLadderSide(topCornerLength/2., halfTheta, -1,
					  fgkLadderLa, fgkLadderHa, fgkLadderl);
  TGeoVolume *cfLaddTopVol1 = new TGeoVolume("ITSsddCFladdTopCornerVol1",
				  cfLaddTop1,carbonFiberLadderStruct);
  cfLaddTopVol1->SetLineColor(fColorCarbonFiber);
  TGeoArb8 *cfLaddTop2 = CreateLadderSide( topCornerLength/2., halfTheta, 1,
					   fgkLadderLa, fgkLadderHa, fgkLadderl);
  TGeoVolume *cfLaddTopVol2 = new TGeoVolume("ITSsddCFladdTopCornerV2",
			          cfLaddTop2,carbonFiberLadderStruct);
  cfLaddTopVol2->SetLineColor(fColorCarbonFiber);
  TGeoTranslation *trTop1 = new TGeoTranslation(0, fgkLadderHeight/2+tDY,
						-(length-topCornerLength)/2.);
  virtualEnd->AddNode(cfLaddTopVol1, 1, trTop1);
  virtualEnd->AddNode(cfLaddTopVol2, 1, trTop1);

  //--- the 2 side V
  TGeoArb8 *cfLaddSide1 = CreateLadderSide( length/2., beta, -1,
					    fgkLadderLb, fgkLadderHb, fgkLadderl);
  TGeoVolume *cfLaddSideVol1 = new TGeoVolume("ITSsddCFladdSideCornerV1",
				   cfLaddSide1,carbonFiberLadderStruct);
  cfLaddSideVol1->SetLineColor(fColorCarbonFiber);
  TGeoArb8 *cfLaddSide2 = CreateLadderSide( length/2., beta, 1,
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
  TGeoRotation beamRot1("", 90-2.*beta*TMath::RadToDeg(),
			-beamPhiPrime*TMath::RadToDeg(), -90);
  TGeoRotation beamRot2("", 90-2.*beta*TMath::RadToDeg(),
			beamPhiPrime*TMath::RadToDeg(), -90);
  TGeoRotation beamRot3("", 90+2.*beta*TMath::RadToDeg(),
			beamPhiPrime*TMath::RadToDeg(), -90);
  TGeoRotation beamRot4("", 90+2.*beta*TMath::RadToDeg(),
			-beamPhiPrime*TMath::RadToDeg(), -90);
  TGeoCombiTrans *beamTransf1 = new TGeoCombiTrans(0.5*triangleHeight*
						   TMath::Tan(halfTheta),
					      fgkLadderBeamRadius/2. + tDY,
			         -length/2 + segmentLength/8, &beamRot1);
  TGeoCombiTrans *beamTransf3 = new TGeoCombiTrans( 0.5*triangleHeight*
						    TMath::Tan(halfTheta),
						 fgkLadderBeamRadius/2.+tDY,
				-length/2 + 3*segmentLength/8, &beamRot2);
  TGeoCombiTrans *beamTransf5 = new TGeoCombiTrans(-0.5*triangleHeight*
						   TMath::Tan(halfTheta),
						fgkLadderBeamRadius/2.+tDY,
			         -length/2 + segmentLength/8, &beamRot3);
  TGeoCombiTrans *beamTransf7 = new TGeoCombiTrans(-0.5*triangleHeight*
						   TMath::Tan(halfTheta),
					      fgkLadderBeamRadius/2. + tDY,
			         -length/2+3*segmentLength/8, &beamRot4);

  virtualEnd->AddNode(fLaddSegCommonVol[6], 1, beamTransf1);
  virtualEnd->AddNode(fLaddSegCommonVol[6], 2, beamTransf3);
  virtualEnd->AddNode(fLaddSegCommonVol[6], 3, beamTransf5);
  virtualEnd->AddNode(fLaddSegCommonVol[6], 4, beamTransf7);

  //--- Beams of the bottom
  TGeoTubeSeg *bottomBeam1 = new TGeoTubeSeg(0, fgkLadderBeamRadius,
				 fgkLadderWidth/2.-fgkLadderLb/3, 0, 180);
  TGeoVolume *bottomBeam1Vol = new TGeoVolume("ITSsddBottomBeam1Vol",
			           bottomBeam1, carbonFiberLadderStruct);
  bottomBeam1Vol->SetLineColor(fColorCarbonFiber);

  TGeoRotation bottomBeamRot1("",90, 90, 90);
  TGeoCombiTrans *bottomBeamTransf1 = new TGeoCombiTrans(0,
                            -(fgkLadderHeight/2-fgkLadderBeamRadius)+tDY,
			    -length/2+fgkSegmentLength/2, &bottomBeamRot1);
  virtualEnd->AddNode(bottomBeam1Vol, 1, bottomBeamTransf1);
  TGeoTubeSeg *bottomBeam2 = new TGeoTubeSeg(0, fgkLadderBeamRadius,
			         fgkLadderWidth/2.-fgkLadderLb/3, 0, 90);
  TGeoVolume *bottomBeam2Vol = new TGeoVolume("ITSsddBottomBeam2Vol",
				   bottomBeam2, carbonFiberLadderStruct);
  bottomBeam2Vol->SetLineColor(fColorCarbonFiber);
  TGeoCombiTrans *bottomBeamTransf2 = new TGeoCombiTrans(0,
     -(fgkLadderHeight/2-fgkLadderBeamRadius)+tDY,-length/2,&bottomBeamRot1);
  virtualEnd->AddNode(bottomBeam2Vol, 1, bottomBeamTransf2);

  //**********************************
  //the cooling pipe supports
  Double_t triangleCPaxeDist = fgkCoolPipeSuppAxeDist-fgkCoolPipeSuppWidthExt-
    fgkCoolPipeSuppWidthIn+fgkLadderBeamRadius;

  Double_t coolPipeSuppL = TMath::Tan(halfTheta)*
    (triangleHeight+triangleCPaxeDist/
     TMath::Sin(halfTheta) - coolPipeSuppH);
  
  if (fAddCoolingSyst) {
  TGeoRotation rotCPS2("",-halfTheta*TMath::RadToDeg(),-90, 90);
  TGeoRotation rotCPS1("", halfTheta*TMath::RadToDeg(),-90,-90);
  TGeoCombiTrans *transCPS1 = new TGeoCombiTrans(coolPipeSuppL,
				  -fgkLadderHeight/2.-TMath::Abs(tDY)+
				  coolPipeSuppH+fgkLadderBeamRadius,
				  -length/2., &rotCPS1);
  TGeoCombiTrans *transCPS4 = new TGeoCombiTrans(-coolPipeSuppL,
				  -fgkLadderHeight/2.-TMath::Abs(tDY)+
				  coolPipeSuppH+fgkLadderBeamRadius,
				  -length/2., &rotCPS2);

  virtualEnd->AddNode(fCoolPipeSupportL, 1, transCPS1);
  virtualEnd->AddNode(fCoolPipeSupportR, 1, transCPS4);
  };

  //**********************************
  if(GetDebug(1)) virtualEnd->CheckOverlaps(0.01);
  //virtualEnd->SetVisibility(kFALSE);
  return virtualEnd;
};


//________________________________________________________________________
TGeoVolume* AliITSv11GeometrySDD::CreateSDDsensor() {
  //
  // return a box containing the SDD sensor
  //

  TGeoMedium *airSDD        = GetMedium("ITSair");
  TGeoMedium *siliconSDD    = GetMedium("ITSsddSi");
  TGeoMedium *alSDD         = GetMedium("ITSal");
  TGeoMedium *polyhamideSDD = GetMedium("ITSsddKAPTON_POLYCH2");
  TGeoMedium *glassSDD      = GetMedium("ITSsddSi");       //  To code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Double_t rWraping = fgkWaferThickness/2+fgkWaHVcableAlThick+fgkWaHVcablePolyThick;
  Double_t witdhCableBox = (fgkWaHVcableWitdh - TMath::Pi()*rWraping)/2;

  Double_t sensoxBoxLength = ( fgkWaferLength +
			       2*(rWraping+witdhCableBox-fgkWaHVcableDW) );
  // Makes life easier to include the space for the WA HV cable on both sides 
  Double_t sensoxBoxThick = fgkWaferThickness +
                            2*(fgkWaHVcableAlThick+fgkWaHVcablePolyThick);

  TGeoBBox *box = new TGeoBBox("ITSsddSensorBox",
		      fgkWaferWidth/2, sensoxBoxThick/2, sensoxBoxLength/2);
  TGeoVolume *virtualSensor = new TGeoVolume("ITSsddSensor",box,airSDD);

  //****************************
  // silicon wafer
  //****************************
  if (fAddSensors) {
    TGeoBBox *waferShape = new TGeoBBox("ITSsddWaferShape",
			   fgkWaferWidth/2, fgkWaferThickness/2, fgkWaferLength/2);
    TGeoVolume *wafer = new TGeoVolume("ITSsddWafer", waferShape, siliconSDD);
    wafer->SetLineColor(fColorSilicon);
    TGeoBBox *sensBox = new TGeoBBox("ITSsddSensorSensBox",
		        fgkWaferWidthSens/2,fgkWaferThickSens/2,fgkWaferLengthSens/2);
    TGeoVolume *sensVol=new TGeoVolume(fgSDDsensitiveVolName,sensBox,siliconSDD);
    sensVol->SetLineColor(fColorSilicon);
    
    wafer->AddNode(sensVol, 1, 0);
    virtualSensor->AddNode(wafer, 1, 0);
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
  virtualSensor->AddNode(vGlass, 1, glassTr1);
  virtualSensor->AddNode(vGlass, 2, glassTr2);
  virtualSensor->AddNode(vGlass, 3, glassTr3);
  virtualSensor->AddNode(vGlass, 4, glassTr4);

  //****************************
  // Wrap-around cable
  //****************************
  if (fAddHVcables) {
  AliITSv11GeomCableFlat waHVCable("ITSsddWaHVCableU",witdhCableBox,
				      fgkWaHVcableAlThick+fgkWaHVcablePolyThick);
  waHVCable.SetNLayers(2);
  waHVCable.SetLayer(0, fgkWaHVcablePolyThick,polyhamideSDD,fColorPolyhamide);
  waHVCable.SetLayer(1, fgkWaHVcableAlThick, alSDD, fColorAl);
  waHVCable.SetInitialNode(virtualSensor);

  Double_t x1[3], x2[3], vX[3] = {1,0,0};
  x1[0] = -fgkWaHVcableLength/2;
  x2[0] = -x1[0];
  x1[1] = (fgkWaferThickness + waHVCable.GetThickness())/2;
  x2[1] = x1[1];
  x1[2] = fgkWaferLength/2+waHVCable.GetWidth()/2-fgkWaHVcableDW;
  x2[2] = x1[2];

  waHVCable.AddCheckPoint(virtualSensor, 0, x1, vX);
  waHVCable.AddCheckPoint(virtualSensor, 1, x2, vX);
  waHVCable.CreateAndInsertCableSegment(1,-90);
  x1[1] = -x1[1];
  x2[1] = x1[1];
  waHVCable.SetName("ITSsddWaHVCableD");
  waHVCable.ResetPoints();
  waHVCable.AddCheckPoint(virtualSensor, 0, x1, vX);
  waHVCable.AddCheckPoint(virtualSensor, 1, x2, vX);
  waHVCable.CreateAndInsertCableSegment(1, 90);

  AliITSv11GeomCableRound waHVCableFold("ITSsddWaHVCableFold",
					   rWraping);
  waHVCableFold.SetPhi(180,360);
  waHVCableFold.SetNLayers(2);
  waHVCableFold.SetLayer(0, fgkWaferThickness/2+fgkWaHVcablePolyThick,
			 polyhamideSDD, fColorPolyhamide);
  waHVCableFold.SetLayer(1, fgkWaHVcableAlThick, alSDD, fColorAl);
  waHVCableFold.SetInitialNode(virtualSensor);
  x1[1] = 0;
  x2[1] = 0;
  x1[2] = fgkWaferLength/2-fgkWaHVcableDW+witdhCableBox;
  x2[2] = x1[2];
  waHVCableFold.AddCheckPoint(virtualSensor, 0, x1, vX);
  waHVCableFold.AddCheckPoint(virtualSensor, 1, x2, vX);
  waHVCableFold.CreateAndInsertCableSegment(1);

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


  TGeoRotation rotHead("",0,90,0);
  TGeoCombiTrans *rotHeadTr = new TGeoCombiTrans(0,fgkWaferThickness/2,
		  -headRadius+fgkTransitHVHeadLZ+fgkTransitHVBondingLZ/2,
						 &rotHead);
  virtualSensor->AddNode(vHeadPolyComp,1,rotHeadTr);
  virtualSensor->AddNode(vHeadAlComp,1,rotHeadTr);

  //---
  AliITSv11GeomCableFlat transitHVCable("ITSsddHVtransitCenter",
					fgkTransitHVBondingLZ,
			    fgkTransitHVPolyThick+fgkTransitHVAlThick);
  transitHVCable.SetNLayers(2);
  transitHVCable.SetLayer(0, fgkTransitHVPolyThick,polyhamideSDD,
			  fColorPolyhamide);
  transitHVCable.SetLayer(1, fgkTransitHVAlThick, alSDD, fColorAl);
  transitHVCable.SetInitialNode(virtualSensor);

  x1[0] = -fgkTransitHVHeadLX/2;
  x2[0] = -x1[0];
  x1[1] = (fgkWaferThickness+fgkTransitHVPolyThick+fgkTransitHVAlThick)/2;
  x2[1] = x1[1];
  x1[2] = 0;
  x2[2] = 0;
  transitHVCable.AddCheckPoint(virtualSensor, 0, x1, vX);
  transitHVCable.AddCheckPoint(virtualSensor, 1, x2, vX);
  transitHVCable.CreateAndInsertCableSegment(1,-90);
  transitHVCable.ResetPoints();
  transitHVCable.SetName("ITSsddHVtransitTail");
  transitHVCable.SetWidth(fgkTransitHVtailWidth);
  x1[0] = fgkTransitHVtailXpos;
  x2[0] = fgkTransitHVtailXpos;
  x1[2] = -fgkTransitHVBondingLZ/2;
  x2[2] = -fgkTransitHVBondingLZ/2-fgkTransitHVtailLength;
  Double_t vZ[3] = {0,0,1};
  transitHVCable.AddCheckPoint(virtualSensor, 0, x1, vZ);
  transitHVCable.AddCheckPoint(virtualSensor, 1, x2, vZ);
  transitHVCable.CreateAndInsertCableSegment(1,0);

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

  TGeoArb8 *sideRight = new TGeoArb8( fgkTransitHVPolyThick/2 );
  sideRight->SetVertex(0, fgkTransitHVtailXpos-fgkTransitHVtailWidth/2,0);
  sideRight->SetVertex(1, fgkTransitHVtailXpos-fgkTransitHVtailWidth/2,
		       fgkTransitHVsideLZ);
  sideRight->SetVertex(2, -fgkTransitHVHeadLX/2, fgkTransitHVsideRightZ);
  sideRight->SetVertex(3, -fgkTransitHVHeadLX/2, 0);
  sideRight->SetVertex(4, fgkTransitHVtailXpos-fgkTransitHVtailWidth/2,0);
  sideRight->SetVertex(5, fgkTransitHVtailXpos-fgkTransitHVtailWidth/2,
		       fgkTransitHVsideLZ);
  sideRight->SetVertex(6, -fgkTransitHVHeadLX/2, fgkTransitHVsideRightZ);
  sideRight->SetVertex(7, -fgkTransitHVHeadLX/2, 0);

  TGeoRotation rotSide("",0,-90,0);
  TGeoCombiTrans *sideRightTr = new TGeoCombiTrans(0,
				(fgkWaferThickness+fgkTransitHVPolyThick)/2,
			        -fgkTransitHVBondingLZ/2,&rotSide);
  TGeoCombiTrans *sideLeftTr = new TGeoCombiTrans(0,
			       (fgkWaferThickness+fgkTransitHVPolyThick)/2,
			       -fgkTransitHVBondingLZ/2, &rotSide);
  TGeoCombiTrans *sideLeftAlTr = new TGeoCombiTrans(0,
		  fgkTransitHVPolyThick+(fgkWaferThickness+fgkTransitHVAlThick)/2,
	       	  -fgkTransitHVBondingLZ/2, &rotSide);

  TGeoVolume *vSideLeft = new TGeoVolume("ITSsddHVtransitSideLeft",
					 sideLeft,polyhamideSDD);
  vSideLeft->SetLineColor(fColorPolyhamide);
  TGeoVolume *vSideLeftAl = new TGeoVolume("ITSsddHVtransitSideLeftAl",
					   sideLeftAl,alSDD);
  vSideLeftAl->SetLineColor(fColorAl);
  TGeoVolume *vSideRight = new TGeoVolume("ITSsddHVtransitSideRight",
					  sideRight,polyhamideSDD);
  vSideRight->SetLineColor(fColorPolyhamide);

  virtualSensor->AddNode(vSideLeft,   1, sideLeftTr);
  virtualSensor->AddNode(vSideLeftAl, 1, sideLeftAlTr);
  virtualSensor->AddNode(vSideRight,  1, sideRightTr);
  };

  //****************************
  if(GetDebug(1)) virtualSensor->CheckOverlaps(0.01);
  virtualSensor->SetVisibility(kFALSE);
  return virtualSensor;
};


//________________________________________________________________________
TGeoVolume *AliITSv11GeometrySDD::CreateDetectors(Int_t iLay) {
  //
  // return a box volume containing the detectors
  //

  TGeoMedium *airSDD = GetMedium("ITSair");

  Int_t    nDetectors   = fgkLay3Ndet;
  Double_t ladderLength = fgkLay3LadderLength;
  Double_t *sensorZPos  = fLay3sensorZPos;
  
  if (iLay==3) {}
  else if (iLay==4) {
    nDetectors   = fgkLay4Ndet;
    ladderLength = fgkLay4LadderLength;
    sensorZPos   = fLay4sensorZPos;
  } else {
    printf("AliITSv11GeometrySDD::CreateLay3Detectors: Error : Wrong layer");
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
	Double_t rotationY = 180;
	if (i >= nDetectors/2) rotationY = 0;
	TGeoRotation rotSensor("",0, rotationY, 0);
	TGeoCombiTrans *sensorPos = new TGeoCombiTrans(0,localY,
						       localZ,&rotSensor);
	sensorPos->SetName(name);
        virtualDet->AddNode(fSDDsensor, i, sensorPos);
    }

    if(GetDebug(1)) virtualDet->CheckOverlaps(0.01);
    //virtualDet->SetVisibility(kFALSE);
    return virtualDet;
};



//________________________________________________________________________
Int_t AliITSv11GeometrySDD::ExportSensorGeometry(AliITSgeom *geom, Int_t iLaySDD,
						 Int_t startMod) {
  //
  // export the geometry in a AliITSgeom object
  //

  if (! geom) {
    printf("error:Try to fill null (AliITSgeom *) object");
    return kFALSE;
  };
  if (! fMotherVol) {
    printf("error:Try to set sensor geometry while geometry is not defined\n");
    return kFALSE;
  };

  Int_t firstSDDmod = startMod;
  const Float_t kDxyz[3] = {fgkWaferWidthSens, fgkWaferThickSens, fgkWaferLengthSens};

  if(!(geom->IsShapeDefined(kSDD)))
    geom->ReSetShape(kSDD, new AliITSgeomSDD256(3, kDxyz));

  char layerName[30];
  char ladderName[30];
  char sensorName[30];
  char senstivName[30];
  const Int_t kNLay = 2;
  const Int_t kNLadd[2] = {fgkLay3Nladd, fgkLay4Nladd};
  const Int_t kNDet[2]  = {fgkLay3Ndet,  fgkLay4Ndet};

  if (GetDebug(1))
    printf("AliITSv11GeometrySDD::SetSensorGeometry(), nodes found :\n");

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
		sprintf(senstivName, "%s%s", fgSDDsensitiveVolName,"_1");
		TGeoNode *sensitivNode = wafVolume->GetNode(senstivName);
	      if (sensitivNode) {
		TGeoHMatrix sensMatrix(wafMatrix);
		sensMatrix.Multiply(sensitivNode->GetMatrix());

		Double_t *trans = sensMatrix.GetTranslation();
		Double_t *r     = sensMatrix.GetRotationMatrix();
		Double_t rot[10] = {r[0],r[1],r[2],
				    r[3],r[4],r[5],
				    r[6],r[7],r[8], 0.0};
		//rot[9]=0.0 => not a unity matrix
		geom->CreatMatrix(startMod,iLay+iLaySDD,iLadd+1,iDet+1,
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
};


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
    lay = 1;
  else lay = 2;

  return kTRUE;
};
