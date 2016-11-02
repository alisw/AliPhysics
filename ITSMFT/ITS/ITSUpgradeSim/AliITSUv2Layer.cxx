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
// This class Defines the Geometry for the ITS Upgrade using TGeo
// This is a work class used to study different configurations
// during the development of the new ITS structure.
//
//  Mario Sitta <sitta@to.infn.it>
//  Chinorat Kobdaj (kobdaj@g.sut.ac.th)
//*************************************************************************


/* $Id: AliITSUv2Layer.cxx  */
// General Root includes
#include <TMath.h>
// Root Geometry includes
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoPcon.h>
#include <TGeoCone.h>
#include <TGeoTube.h> // contaings TGeoTubeSeg
#include <TGeoArb8.h>
#include <TGeoXtru.h>
#include <TGeoCompositeShape.h>
#include <TGeoMatrix.h>
#include "AliITSUv2Layer.h"
#include "AliITSUGeomTGeo.h"
#include <TGeoBBox.h>
#include <TGeoShape.h>
#include <TGeoTrd1.h>

#include <AliLog.h>


using namespace TMath;

// General Parameters
const Double_t AliITSUv2Layer::fgkmicron = 1.0E-4;
const Double_t AliITSUv2Layer::fgkmm = 0.10;
const Double_t AliITSUv2Layer::fgkcm = 1.00;

const Int_t    AliITSUv2Layer::fgkNumberOfInnerLayers =   3;

const Double_t AliITSUv2Layer::fgkDefaultSensorThick  =  18*fgkmicron;
const Double_t AliITSUv2Layer::fgkDefaultChipThick    =  50*fgkmicron;

// Inner Barrel Parameters
const Int_t    AliITSUv2Layer::fgkIBChipsPerRow       =   9;
const Int_t    AliITSUv2Layer::fgkIBNChipRows         =   1;

const Double_t AliITSUv2Layer::fgkIBFlexCableAlThick  =  50.0  *fgkmicron;
const Double_t AliITSUv2Layer::fgkIBFlexCableKapThick = 125.0  *fgkmicron;
const Double_t AliITSUv2Layer::fgkIBGlueThick         = 100.0  *fgkmicron;
const Double_t AliITSUv2Layer::fgkIBCarbonFleeceThick =  20.0  *fgkmicron;
const Double_t AliITSUv2Layer::fgkIBCarbonPaperThick  =  30.0  *fgkmicron;
const Double_t AliITSUv2Layer::fgkIBK13D2UThick       =  70.0  *fgkmicron;
const Double_t AliITSUv2Layer::fgkIBCoolPipeInnerD    =   1.024*fgkmm;
const Double_t AliITSUv2Layer::fgkIBCoolPipeThick     =  25.4  *fgkmicron;
const Double_t AliITSUv2Layer::fgkIBCoolPipeXDist     =   5.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBTopVertexWidth1   =   0.258*fgkmm;
const Double_t AliITSUv2Layer::fgkIBTopVertexWidth2   =   0.072*fgkcm;
const Double_t AliITSUv2Layer::fgkIBTopVertexHeight   =   0.04 *fgkcm;
const Double_t AliITSUv2Layer::fgkIBTopVertexAngle    =  60.0; // Deg
const Double_t AliITSUv2Layer::fgkIBSideVertexWidth   =   0.05 *fgkcm;
const Double_t AliITSUv2Layer::fgkIBSideVertexHeight  =   0.074*fgkcm;
const Double_t AliITSUv2Layer::fgkIBTopFilamentLength =   0.9  *fgkcm;
const Double_t AliITSUv2Layer::fgkIBTopFilamentSide   =   0.02 *fgkcm;
const Double_t AliITSUv2Layer::fgkIBTopFilamentAlpha  =  57.0; // Deg
const Double_t AliITSUv2Layer::fgkIBTopFilamentGamma  =  65.0; // Deg

const Double_t AliITSUv2Layer::fgkIBConnectorXWidth   =  10.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnectorYTot     =   4.7  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnectBlockZLen  =  16.5  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnBodyYHeight   =   2.5  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnTailYShift    =   0.9  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnTailYMid      =   2.5  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnTailZLen      =   2.5  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnTailOpenPhi   = 120.0; // Deg
const Double_t AliITSUv2Layer::fgkIBConnRoundHoleD    =   2.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnRoundHoleZ    =(7.0-2.0)*fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnSquareHoleX   =   3.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnSquareHoleZ   =   3.3  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnSquareHoleZPos=   9.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnInsertHoleD   =   3.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnInsertHoleZPos=   9.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnTubeHole1D    =   1.6  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnTubeHole1ZLen =   3.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnTubeHole2D    =   1.2  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnTubeHole3XPos =   1.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnTubeHole3ZPos =   2.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnTubesXDist    =   5.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnTubesYPos     =   1.25 *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnInsertInnerX  =   2.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnInsertZThick  =   0.7  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnInsertD       =   2.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnInsertHeight  =   2.3  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnectAFitExtD   =   1.65 *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnectAFitIntD   =   1.19 *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnectAFitZLen   =  12.5  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnectAFitZOut   =  10.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnPlugInnerD    =   0.8  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnPlugTotLen    =   1.7  *fgkmm;
const Double_t AliITSUv2Layer::fgkIBConnPlugThick     =   0.5  *fgkmm;

const Double_t AliITSUv2Layer::fgkIBStaveHeight       =   0.5  *fgkcm;

// Outer Barrel Parameters
const Int_t    AliITSUv2Layer::fgkOBChipsPerRow       =   7;
const Int_t    AliITSUv2Layer::fgkOBNChipRows         =   2;

const Double_t AliITSUv2Layer::fgkOBHalfStaveWidth    =   3.01 *fgkcm;
const Double_t AliITSUv2Layer::fgkOBModuleWidth       = fgkOBHalfStaveWidth;
const Double_t AliITSUv2Layer::fgkOBModuleGap         =   0.01 *fgkcm;
const Double_t AliITSUv2Layer::fgkOBChipXGap          =   0.01 *fgkcm;
const Double_t AliITSUv2Layer::fgkOBChipZGap          =   0.01 *fgkcm;
const Double_t AliITSUv2Layer::fgkOBFlexCableAlThick  =   0.005*fgkcm;
const Double_t AliITSUv2Layer::fgkOBFlexCableCuThick  =   0.004*fgkcm;
const Double_t AliITSUv2Layer::fgkOBFlexCableKapThick1=   0.01 *fgkcm;
const Double_t AliITSUv2Layer::fgkOBFlexCableKapThick = 125.0  *fgkmicron;
const Double_t AliITSUv2Layer::fgkOBBusCableAlThick   =   0.02 *fgkcm;
const Double_t AliITSUv2Layer::fgkOBBusCableKapThick  =   0.02 *fgkcm;
const Double_t AliITSUv2Layer::fgkOBColdPlateThick    =   0.012*fgkcm;
const Double_t AliITSUv2Layer::fgkOBCarbonPlateThick  =   0.012*fgkcm;
const Double_t AliITSUv2Layer::fgkOBGlueThickM1       =   0.03 *fgkcm;
const Double_t AliITSUv2Layer::fgkOBGlueThick         =   0.01 *fgkcm;
const Double_t AliITSUv2Layer::fgkOBModuleZLength     =  21.06 *fgkcm;
const Double_t AliITSUv2Layer::fgkOBHalfStaveYPos     =   2.067*fgkcm;
const Double_t AliITSUv2Layer::fgkOBHalfStaveYTrans   =   1.76 *fgkmm;
const Double_t AliITSUv2Layer::fgkOBHalfStaveXOverlap =   4.3  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBGraphiteFoilThick =  30.0  *fgkmicron;
const Double_t AliITSUv2Layer::fgkOBCarbonFleeceThick =  20.0  *fgkmicron;
const Double_t AliITSUv2Layer::fgkOBCoolTubeInnerDM1  =   2.052*fgkmm;
const Double_t AliITSUv2Layer::fgkOBCoolTubeInnerD    =   2.05 *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCoolTubeThick     =  32.0  *fgkmicron;
const Double_t AliITSUv2Layer::fgkOBCoolTubeXDist     =  10.0  *fgkmm;

const Double_t AliITSUv2Layer::fgkOBCPConnectorXWidth =  16.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnBlockZLen   =  15.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnBlockYHei   =   3.6  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnHollowZLen  =   3.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnHollowYHei  =   0.9  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnSquareHoleX =   6.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnSquareHoleZ =   6.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnSqrHoleZPos =   4.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnSqrInsertRZ =   3.5  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnRoundHoleD  =   6.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnRndHoleZPos =   4.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnTubesXDist  =  10.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnTubesYPos   =   1.8  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnTubeHole1D  =   2.6  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnTubeHole1Z  =   3.5  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnTubeHole2D  =   2.2  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnFitHoleD    =   2.8  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnTubeHole3XP =   1.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnTubeHole3ZP =   2.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnInstInnerX  =   4.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnInstInnerR  =   1.5  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnInstZThick  =   1.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnInsertYHei  =   3.4  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnInsertD     =   4.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnAFitExtD    =   2.77 *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnAFitThick   =   0.3  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnAFitZLen    =  25.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnAFitZOut    =  22.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnPlugInnerD  =   0.8  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnPlugTotLen  =   1.7  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBCPConnPlugThick   =   0.5  *fgkmm;

const Double_t AliITSUv2Layer::fgkOBSpaceFrameZLen[2] = { 900.0*fgkmm,
							 1526.0*fgkmm};
const Int_t    AliITSUv2Layer::fgkOBSpaceFrameNUnits[2]= { 23, 39};
const Double_t AliITSUv2Layer::fgkOBSpaceFrameUnitLen =  39.1  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBSpaceFrameWidth   =  42.44 *fgkmm;
const Double_t AliITSUv2Layer::fgkOBSpaceFrameHeight  =  36.45 *fgkmm;
const Double_t AliITSUv2Layer::fgkOBSpaceFrameTopVL   =   4.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBSpaceFrameTopVH   =   0.35 *fgkmm;
const Double_t AliITSUv2Layer::fgkOBSpaceFrameSideVL  =   4.5  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBSpaceFrameSideVH  =   0.35 *fgkmm;
const Double_t AliITSUv2Layer::fgkOBSpaceFrameVAlpha  =  60.0; // deg
const Double_t AliITSUv2Layer::fgkOBSpaceFrameVBeta   =  68.0; // deg
const Double_t AliITSUv2Layer::fgkOBSFrameBaseRibDiam =   1.33 *fgkmm;
const Double_t AliITSUv2Layer::fgkOBSFrameBaseRibPhi  =  54.0; // deg
const Double_t AliITSUv2Layer::fgkOBSFrameSideRibDiam =   1.25 *fgkmm;
const Double_t AliITSUv2Layer::fgkOBSFrameSideRibPhi  =  70.0; // deg
const Double_t AliITSUv2Layer::fgkOBSFrameULegLen     =  14.2  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBSFrameULegWidth   =   1.5  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBSFrameULegHeight1 =   2.7  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBSFrameULegHeight2 =   5.0  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBSFrameULegThick   =   0.3  *fgkmm;
const Double_t AliITSUv2Layer::fgkOBSFrameULegXPos    =  12.9  *fgkmm;

// Barrel Geomatrix Parameters
// #pnamwong

//const Double_t AliITSUv2Layer::wrpZSpan[3] = {32.0, 90.1, 152.7};		// from CreateITSUv2.C
// Inner Barrel
const Double_t AliITSUv2Layer::fgkIBStaveLength       	=   29.0  *fgkcm;
Double_t AliITSUv2Layer::fgkIBConnOffset        		=   0.0; 
//const Double_t AliITSUv2Layer::fgkIBConnOffset        =   15.0  *fgkcm;
//const Double_t AliITSUv2Layer::fgkIBConnOffset        =   fZLength/2; don't work
const Double_t AliITSUv2Layer::fgkIBConnOffsetExt 		= 	5.0	*fgkmm; // -0.1 for avoid the overlaps with stave_struct

// Inner Barrel - End Wheel 
//const TGeoMedium *AliITSUv2Layer::fkIBEWheelMedium = gGeoManager->GetMedium("ITS_KAPTON(POLYCH2)$");
const Double_t AliITSUv2Layer::fkALC_0334_RotateAngle[3] = {30., 22.5, 18.}; // unit=angle
const Double_t AliITSUv2Layer::fkALC_0334_RotateOffset[3] = {-1., -5.5, -7.5}; // unit=angle

const Double_t AliITSUv2Layer::fkALC_0334_PlateRadius[3] = {27.9*fgkmm, 35.9*fgkmm, 43.9*fgkmm};
const Double_t AliITSUv2Layer::fkALC_0334_PlateThick = 0.6*fgkmm;

const Double_t AliITSUv2Layer::fkALC_0334_ContactDX[3] = {-7.*fgkmm, -9.5*fgkmm, -11.5*fgkmm};
//const Double_t AliITSUv2Layer::fkALC_0334_ContactDY[3] = {24.4*fgkmm, 32.1*fgkmm, 39.6*fgkmm}; 
const Double_t AliITSUv2Layer::fkALC_0334_ContactDY[3] = {(24.4*fgkmm)+(0.1*fgkmm), (32.1*fgkmm)+(0.1*fgkmm), (39.6*fgkmm)+(0.1*fgkmm)}; 

const Double_t AliITSUv2Layer::fkALC_0334_HoleDia =   2.3*fgkmm;
const Double_t AliITSUv2Layer::fkALC_0334_HoleFromEdge = 4.0*fgkmm;		// center-of-hole to inside-edge
const Double_t AliITSUv2Layer::fkALC_0334_HoleLength = 10.0*fgkmm;		// approximate for hole object

const Double_t AliITSUv2Layer::fkALC_0334_HndLength = 14.0*fgkmm;
const Double_t AliITSUv2Layer::fkALC_0334_HndWidth = 10.3*fgkmm;
const Double_t AliITSUv2Layer::fkALC_0334_HndThick = 3.5*fgkmm;

// Thermal Connector have to be offset with only Y-axis, not X-axis, of local space due to the extended of DCDC
const Double_t AliITSUv2Layer::fkThermalConnOffsetZ[3] = {+154.8845*fgkmm, +148.8890*fgkmm, +142.8967*fgkmm};
const Double_t AliITSUv2Layer::fkThermalConnOffsetY[3] = {+ 18.8926*fgkmm,  +41.4403*fgkmm,  +64.4276*fgkmm};

const Double_t AliITSUv2Layer::fkDCDCOffsetZ[3] = { +98.9550*fgkmm, +94.5536*fgkmm, +87.7700*fgkmm};
const Double_t AliITSUv2Layer::fkDCDCOffsetY[3] = { +14.4330*fgkmm, (+26.4747+1.0)*fgkmm, +37.6518*fgkmm};

const Double_t AliITSUv2Layer::fkDCDCRotateBeta[3] = {5.1557, 17.0, 28.0};

// ETC
const Double_t AliITSUv2Layer::fgkInnerBarrelOffsetX		= 0.0;
const Double_t AliITSUv2Layer::fgkInnerBarrelOffsetY		= 0.0;
const Double_t AliITSUv2Layer::fgkInnerBarrelOffsetZ		= 3.0 *fgkcm;

const Double_t AliITSUv2Layer::fgkOuterBarrelOffsetX		= 0.0;
const Double_t AliITSUv2Layer::fgkOuterBarrelOffsetY		= 0.0;
const Double_t AliITSUv2Layer::fgkOuterBarrelOffsetZ		= -7.0 *fgkcm;

const Double_t AliITSUv2Layer::fgkBarrelRotateAlpha	= 0.0;
const Double_t AliITSUv2Layer::fgkBarrelRotateBetha	= 180.0; // for A <-> C side
const Double_t AliITSUv2Layer::fgkBarrelRotateGamma	= 0.0;

//#define _ITSSB_SHOW_INDICATOR
//#define _ITSSB_SHOW_COMBISHAPE

ClassImp(AliITSUv2Layer)

#define SQ(A) (A)*(A)

//________________________________________________________________________
AliITSUv2Layer::AliITSUv2Layer(): 
  TObject(),
  fLayerNumber(0),
  fPhi0(0),
  fLayRadius(0),
  fZLength(0),
  fSensorThick(0),
  fChipThick(0),
  fStaveWidth(0),
  fStaveTilt(0),
  fNStaves(0),
  fNModules(0),
  fNChips(0),
  fChipTypeID(0),
  fIsTurbo(0),
  fBuildLevel(0),
  fStaveModel(AliITSUv2::kIBModelDummy),
  fAddGammaConv(kFALSE),
  fGammaConvDiam(0),
  fGammaConvXPos(0),
  fNHandleCreated(0),
  fN_DCDC_Created(0),
  fN_DCCNT_Created(0),
  fN_SamtecCable_Created(0),
  fN_SamtecCables_Created(0)  
{
  //
  // Standard constructor
  for (int i=kNHLevels;i--;) fHierarchy[i] = 0;
  //

}

//________________________________________________________________________
AliITSUv2Layer::AliITSUv2Layer(Int_t debug): 
  TObject(),
  fLayerNumber(0),
  fPhi0(0),
  fLayRadius(0),
  fZLength(0),
  fSensorThick(0),
  fChipThick(0),
  fStaveWidth(0),
  fStaveTilt(0),
  fNStaves(0),
  fNModules(0),
  fNChips(0),
  fChipTypeID(0),
  fIsTurbo(0),
  fBuildLevel(0),
  fStaveModel(AliITSUv2::kIBModelDummy),
  fAddGammaConv(kFALSE),
  fGammaConvDiam(0),
  fGammaConvXPos(0),
  fNHandleCreated(0),
  fN_DCDC_Created(0),
  fN_DCCNT_Created(0),
  fN_SamtecCable_Created(0),
  fN_SamtecCables_Created(0)
{
  //
  // Constructor setting debugging level
  for (int i=kNHLevels;i--;) fHierarchy[i] = 0;
  //

}

//________________________________________________________________________
AliITSUv2Layer::AliITSUv2Layer(Int_t lay, Int_t debug): 
  TObject(),
  fLayerNumber(lay),
  fPhi0(0),
  fLayRadius(0),
  fZLength(0),
  fSensorThick(0),
  fChipThick(0),
  fStaveWidth(0),
  fStaveTilt(0),
  fNStaves(0),
  fNModules(0),
  fNChips(0),
  fChipTypeID(0),
  fIsTurbo(0),
  fBuildLevel(0),
  fStaveModel(AliITSUv2::kIBModelDummy),
  fAddGammaConv(kFALSE),
  fGammaConvDiam(0),
  fGammaConvXPos(0),
  fNHandleCreated(0),
  fN_DCDC_Created(0),
  fN_DCCNT_Created(0),
  fN_SamtecCable_Created(0),
  fN_SamtecCables_Created(0)

{
  //
  // Constructor setting layer number and debugging level
  for (int i=kNHLevels;i--;) fHierarchy[i] = 0;
  //

}

//________________________________________________________________________
AliITSUv2Layer::AliITSUv2Layer(Int_t lay, Bool_t turbo, Int_t debug): 
  TObject(),
  fLayerNumber(lay),
  fPhi0(0),
  fLayRadius(0),
  fZLength(0),
  fSensorThick(0),
  fChipThick(0),
  fStaveWidth(0),
  fStaveTilt(0),
  fNStaves(0),
  fNModules(0),
  fNChips(0),
  fChipTypeID(0),
  fIsTurbo(turbo),
  fBuildLevel(0),
  fStaveModel(AliITSUv2::kIBModelDummy),
  fAddGammaConv(kFALSE),
  fGammaConvDiam(0),
  fGammaConvXPos(0),
  fNHandleCreated(0),
  fN_DCDC_Created(0),
  fN_DCCNT_Created(0),
  fN_SamtecCable_Created(0),
  fN_SamtecCables_Created(0)

{
  //
  // Constructor setting layer number and debugging level
  // for a "turbo" layer (i.e. where staves overlap in phi)
  for (int i=kNHLevels;i--;) fHierarchy[i] = 0;
  //
 
}

//________________________________________________________________________
AliITSUv2Layer::AliITSUv2Layer(const AliITSUv2Layer &s):
  TObject(),
  fLayerNumber(s.fLayerNumber),
  fPhi0(s.fPhi0),
  fLayRadius(s.fLayRadius),
  fZLength(s.fZLength),
  fSensorThick(s.fSensorThick),
  fChipThick(s.fChipThick),
  fStaveWidth(s.fStaveWidth),
  fStaveTilt(s.fStaveTilt),
  fNStaves(s.fNStaves),
  fNModules(s.fNModules),
  fNChips(s.fNChips),
  fChipTypeID(s.fChipTypeID),
  fIsTurbo(s.fIsTurbo),
  fBuildLevel(s.fBuildLevel),
  fStaveModel(s.fStaveModel),
  fAddGammaConv(s.fAddGammaConv),
  fGammaConvDiam(s.fGammaConvDiam),
  fGammaConvXPos(s.fGammaConvXPos),
  fNHandleCreated(s.fNHandleCreated),
  fN_DCDC_Created(s.fN_DCDC_Created),
  fN_DCCNT_Created(s.fN_DCCNT_Created),
  fN_SamtecCable_Created(s.fN_SamtecCable_Created),
  fN_SamtecCables_Created(s.fN_SamtecCables_Created)

{
  //
  // Copy constructor
  for (int i=kNHLevels;i--;) fHierarchy[i] = s.fHierarchy[i];
  //
	  // #pnamwong
  fgkIBConnOffset = s.fZLength/2 + fgkIBConnOffsetExt;
}

//________________________________________________________________________
AliITSUv2Layer& AliITSUv2Layer::operator=(const AliITSUv2Layer &s)
{
  //
  // Assignment operator 
  //
  if(&s == this) return *this;

  fLayerNumber   = s.fLayerNumber;
  fPhi0          = s.fPhi0;
  fLayRadius     = s.fLayRadius;
  fZLength       = s.fZLength;
  fSensorThick   = s.fSensorThick;
  fChipThick     = s.fChipThick;
  fStaveWidth    = s.fStaveWidth;
  fStaveTilt     = s.fStaveTilt;
  fNStaves       = s.fNStaves;
  fNModules      = s.fNModules;
  fNChips        = s.fNChips;
  fIsTurbo       = s.fIsTurbo;
  fChipTypeID    = s.fChipTypeID;
  fBuildLevel    = s.fBuildLevel;
  fStaveModel    = s.fStaveModel;
  fAddGammaConv  = s.fAddGammaConv;
  fGammaConvDiam = s.fGammaConvDiam;
  fGammaConvXPos = s.fGammaConvXPos;

  //
  // #pnamwong
  fNHandleCreated = s.fNHandleCreated;
  fN_DCDC_Created = s.fN_DCDC_Created;
  fN_DCCNT_Created = s.fN_DCCNT_Created;
  fN_SamtecCable_Created = s.fN_SamtecCable_Created;
  fN_SamtecCables_Created = s.fN_SamtecCables_Created;

  for (int i=kNHLevels;i--;) fHierarchy[i] = s.fHierarchy[i];
  //
  // #pnamwong
  fgkIBConnOffset = s.fZLength/2 + fgkIBConnOffsetExt;
  
  return *this;
}

//________________________________________________________________________
AliITSUv2Layer::~AliITSUv2Layer() {
  //
  // Destructor
  //
}

//________________________________________________________________________
void AliITSUv2Layer::AddGammaConversionRods(const Double_t diam,
					    const Double_t xPos){
//
// Adds the gamma conversion rods for the current layer
//
// Input:
//         diam : the rods diameter
//         diam : the rods X position in the Half Stave refernce system
//
// Output:
//
// Return:
//
// Created:      26 Oct 2016  Mario Sitta
//

  if (fLayerNumber < fgkNumberOfInnerLayers) {
    AliError("No Gamma Conversion Rods in the Inner Layers!");
    return;
  }

  fAddGammaConv = kTRUE; // By default set to False in the constructor
  fGammaConvDiam = diam;
  fGammaConvXPos = xPos;

}

//________________________________________________________________________
void AliITSUv2Layer::CreateLayer(TGeoVolume *moth){
//
// Creates the actual Layer and places inside its mother volume
//
// Input:
//         moth : the TGeoVolume owing the volume structure
//
// Output:
//
// Return:
//
// Created:      17 Jun 2011  Mario Sitta
// Updated:      08 Jul 2011  Mario Sitta
// Updated:      20 May 2013  Mario Sitta  Layer is Assembly instead of Tube
//
  // Local variables
  char volname[30];
  Double_t xpos, ypos, zpos;
  Double_t alpha;



  // Check if the user set the proper parameters
  if (fLayRadius<= 0) AliFatal(Form("Wrong layer radius (%f)",fLayRadius));
  if (fZLength  <= 0) AliFatal(Form("Wrong layer length (%f)",fZLength));
  if (fNStaves  <= 0) AliFatal(Form("Wrong number of staves (%d)",fNStaves));
  if (fNChips   <= 0) AliFatal(Form("Wrong number of chips (%d)",fNChips));

  if (fLayerNumber >= fgkNumberOfInnerLayers && fNModules <= 0)
    AliFatal(Form("Wrong number of modules (%d)",fNModules));

  if (fChipThick <= 0) {
    AliInfo(Form("Chip thickness wrong or not set (%f), using default (%f)",
		 fChipThick,fgkDefaultChipThick));
    fChipThick = fgkDefaultChipThick;
  }

  if (fSensorThick <= 0) {
    AliInfo(Form("Sensor thickness wrong or not set (%f), using default (%f)",
		 fSensorThick,fgkDefaultSensorThick));
    fSensorThick = fgkDefaultSensorThick;
  }

  if (fSensorThick > fChipThick) {
    AliWarning(Form("Sensor thickness (%f) is greater than chip thickness (%f), fixing",
		 fSensorThick,fChipThick));
    fSensorThick = fChipThick;
  }

  // #pnamwong
  fgkIBConnOffset = fZLength/2 + fgkIBConnOffsetExt;
  
  // If a Turbo layer is requested, do it and exit
  if (fIsTurbo) {
    CreateLayerTurbo(moth);
    return;
  }


  
  // First create the stave container
  alpha = (360./(2*fNStaves))*DegToRad();

  //  fStaveWidth = fLayRadius*Tan(alpha);

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSLayerPattern(),fLayerNumber);
  TGeoVolume *layVol = new TGeoVolumeAssembly(volname);
  layVol->SetTitle(
	   Form("Model %d Build level %d",(Int_t)fStaveModel,fBuildLevel));
  layVol->SetUniqueID(fChipTypeID);

//  layVol->SetVisibility(kFALSE);
  layVol->SetVisibility(kTRUE);
  layVol->SetLineColor(1);

  TGeoVolume *stavVol = CreateStave();


  // Now build up the layer
  alpha = 360./fNStaves;
  Double_t r = fLayRadius + ((TGeoBBox*)stavVol->GetShape())->GetDY();
  for (Int_t j=0; j<fNStaves; j++) {
    Double_t phi = j*alpha + fPhi0;
    xpos = r*CosD(phi);// r*SinD(-phi);
    ypos = r*SinD(phi);// r*CosD(-phi);
    zpos = 0.;
    phi += 90;
    layVol->AddNode(stavVol, j, new TGeoCombiTrans( xpos, ypos, zpos,
						    new TGeoRotation("",phi,0,0)));
  }

  layVol->GetShape()->ComputeBBox(); //RS: enfore recompting of BBox

  // Finally put everything in the mother volume
  moth->AddNode(layVol, 1, 0);


  // Upgrade geometry is served
  return;
}

//________________________________________________________________________
void AliITSUv2Layer::CreateLayerTurbo(TGeoVolume *moth){
//
// Creates the actual Layer and places inside its mother volume
// A so-called "turbo" layer is a layer where staves overlap in phi
// User can set width and tilt angle, no check is performed here
// to avoid volume overlaps
//
// Input:
//         moth : the TGeoVolume owing the volume structure
//
// Output:
//
// Return:
//
// Created:      08 Jul 2011  Mario Sitta
// Updated:      08 Mar 2012  Mario Sitta  Correct way to compute container R
// Updated:      20 May 2013  Mario Sitta  Layer is Assemgbly instead of Tube
//


  // Local variables
  char volname[30];
  Double_t xpos, ypos, zpos;
  Double_t alpha;


  // Check if the user set the proper (remaining) parameters
  if (fStaveWidth <= 0)
    AliFatal(Form("Wrong stave width (%f)",fStaveWidth));
  if (Abs(fStaveTilt) > 45)
    AliWarning(Form("Stave tilt angle (%f) greater than 45deg",fStaveTilt));


  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSLayerPattern(), fLayerNumber);
  TGeoVolume *layVol = new TGeoVolumeAssembly(volname);
  layVol->SetTitle(
	   Form("Model %d Build level %d",(Int_t)fStaveModel,fBuildLevel));
  layVol->SetUniqueID(fChipTypeID);
  layVol->SetVisibility(kTRUE);
  layVol->SetLineColor(1);
  TGeoVolume *stavVol = CreateStave();


  // Now build up the layer
  alpha = 360./fNStaves;
  Double_t r = fLayRadius; // Chip thick taken into account in ITSUModule
  for (Int_t j=0; j<fNStaves; j++) {
    Double_t phi = j*alpha + fPhi0;
    xpos = r*CosD(phi);// r*SinD(-phi);
    ypos = r*SinD(phi);// r*CosD(-phi);
    zpos = 0.;
    phi += 90;
    layVol->AddNode(stavVol, j, new TGeoCombiTrans( xpos, ypos, zpos,
						    new TGeoRotation("", phi-fStaveTilt,0,0)));
  }

  layVol->GetShape()->ComputeBBox(); //RS: enfore recompting of BBox

  // Finally put everything in the mother volume
  moth->AddNode(layVol, 1, 0);

  return;
}

//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateStave(const TGeoManager * /*mgr*/){
//
// Creates the actual Stave
//
// Input:
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      22 Jun 2011  Mario Sitta
// Updated:      18 Dec 2013  Mario Sitta  Handle IB and OB
// Updated:      12 Jan 2015  Mario Sitta  Fix overlap with new OB space frame
//                            (by moving the latter, not the sensors to avoid
//                             spoiling their position in space)
// Updated:      03 Mar 2015  Mario Sitta  Fix chip position
//

  char volname[30];
 
  Double_t xlen, ylen, zlen;
  Double_t xpos, ypos;
  Double_t alpha;


  // First create all needed shapes
  alpha = (360./(2*fNStaves))*DegToRad();

  // The stave
  xlen = fLayRadius*Tan(alpha);
  if (fIsTurbo) xlen = 0.5*fStaveWidth;
  ylen = 0.5*fChipThick;
  zlen = 0.5*fZLength;

  Double_t yplus = 0.46;
  TGeoXtru *stave = new TGeoXtru(2); //z sections
  Double_t xv[5] = {xlen,xlen,0,-xlen,-xlen};
  Double_t yv[5] = {ylen+0.09,-0.15,-yplus-fSensorThick,-0.15,ylen+0.09};    
  stave->DefinePolygon(5,xv,yv);
  stave->DefineSection(0,-zlen,0,0,1.);
  stave->DefineSection(1,+zlen,0,0,1.);

  // We have all shapes: now create the real volumes

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSStavePattern(), fLayerNumber);
//  TGeoVolume *staveVol = new TGeoVolume(volname, stave, medAir);
  TGeoVolume *staveVol = new TGeoVolumeAssembly(volname);

  //  staveVol->SetVisibility(kFALSE);
  staveVol->SetVisibility(kTRUE);
  staveVol->SetLineColor(2);
  TGeoVolume *mechStaveVol = 0;

  // Now build up the stave
  if (fLayerNumber < fgkNumberOfInnerLayers) {
    TGeoVolume *modVol = CreateStaveInnerB(xlen,ylen,zlen);
    ypos = ((TGeoBBox*)(modVol->GetShape()))->GetDY() - fChipThick; // = 0 if not kIBModel4
    staveVol->AddNode(modVol, 0, new TGeoTranslation(0, ypos, 0));
    fHierarchy[kHalfStave] = 1;
 
  // Mechanical stave structure
    mechStaveVol = CreateStaveStructInnerB(xlen,zlen); 
    if (mechStaveVol) {
      ypos = ((TGeoBBox*)(modVol->GetShape()))->GetDY() - ypos;
      if (fStaveModel != AliITSUv2::kIBModel4)
	ypos += ((TGeoBBox*)(mechStaveVol->GetShape()))->GetDY();
      staveVol->AddNode(mechStaveVol, 1, new TGeoCombiTrans(0, -ypos, 0, new TGeoRotation("",0, 0, 180)));
    }
  }

  else{
    TGeoVolume *hstaveVol = CreateStaveOuterB();
    if (fStaveModel == AliITSUv2::kOBModel0) { // Create simplified stave struct as in v0
      staveVol->AddNode(hstaveVol, 0);
      fHierarchy[kHalfStave] = 1;
    } else { // (if fStaveModel) Create new stave struct as in TDR
      xpos = ((TGeoBBox*)(hstaveVol->GetShape()))->GetDX()
	   - fgkOBHalfStaveXOverlap/2;
      // ypos is now a parameter to avoid HS displacement wrt nominal radii
      ypos = fgkOBHalfStaveYPos;
      staveVol->AddNode(hstaveVol, 0, new TGeoTranslation(-xpos, ypos, 0));
      staveVol->AddNode(hstaveVol, 1, new TGeoTranslation( xpos, ypos+fgkOBHalfStaveYTrans, 0));
      fHierarchy[kHalfStave] = 2; // RS 
      mechStaveVol = CreateSpaceFrameOuterB();
      if (mechStaveVol) {
	if (fBuildLevel < 6)   // Carbon
	  staveVol->AddNode(mechStaveVol, 1,
			    new TGeoCombiTrans(0, -fgkOBSFrameULegHeight1, 0,
					     new TGeoRotation("", 180, 0, 0)));
      }
    } // if (fStaveModel)
  }
  
  staveVol->GetShape()->ComputeBBox(); //RS: enfore recompting of BBox

  // Done, return the stave
  return staveVol;
}

//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateStaveInnerB(const Double_t xsta,
					      const Double_t ysta,
					      const Double_t zsta,
					      const TGeoManager *mgr){
//
// Create the chip stave for the Inner Barrel
// (Here we fake the halfstave volume to have the same
// formal geometry hierarchy as for the Outer Barrel)
//
// Input:
//         xsta, ysta, zsta : X, Y, Z stave half lengths
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      06 Mar 2014  Mario Sitta
//

  // Local variables
  Double_t xmod, ymod, zmod;
  char volname[30];

  // First we create the module (i.e. the HIC with 9 chips)
  TGeoVolume *moduleVol = CreateModuleInnerB(xsta, ysta, zsta);

  // Then we create the fake halfstave and the actual stave
  xmod = ((TGeoBBox*)(moduleVol->GetShape()))->GetDX();
  ymod = ((TGeoBBox*)(moduleVol->GetShape()))->GetDY();
  zmod = ((TGeoBBox*)(moduleVol->GetShape()))->GetDZ();

  TGeoBBox *hstave = new TGeoBBox(xmod, ymod, zmod);

  TGeoMedium *medAir = mgr->GetMedium("ITS_AIR$");

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSHalfStavePattern(), fLayerNumber);
  TGeoVolume *hstaveVol  = new TGeoVolume(volname, hstave, medAir);


  // Finally build it up
  hstaveVol->AddNode(moduleVol, 0);
  fHierarchy[kModule] = 1;

  // Done, return the stave structure
  return hstaveVol;
}

//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateModuleInnerB(Double_t xmod,
					       Double_t ymod,
					       Double_t zmod,
					       const TGeoManager *mgr){
//
// Creates the IB Module: (only the chips for the time being)
//
// Input:
//         xmod, ymod, zmod : X, Y, Z module half lengths
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//         the module as a TGeoVolume
//
// Created:      06 Mar 2014  M. Sitta
// Updated:      03 Mar 2015  Mario Sitta  FPC in right position (beyond chip)
//

  Double_t ytot, zchip;
  Double_t ypos, zpos;
  char volname[30];

  // First create the single chip
  zchip = zmod/fgkIBChipsPerRow;
  TGeoVolume *chipVol = CreateChipInnerB(xmod, ymod, zchip);

  // Then create the module and populate it with the chips
  // (and the FPC Kapton and Aluminum in the most recent IB model)
  ytot = ymod;
  if (fStaveModel == AliITSUv2::kIBModel4)
    ytot += 0.5*(fgkIBFlexCableKapThick + fgkIBFlexCableAlThick);

  TGeoBBox *module = new TGeoBBox(xmod, ytot, zmod);

  TGeoBBox *kapCable = new TGeoBBox(xmod, fgkIBFlexCableKapThick/2, zmod);
  TGeoBBox *aluCable = new TGeoBBox(xmod, fgkIBFlexCableAlThick /2, zmod);

  TGeoMedium *medAir      = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medKapton   = mgr->GetMedium("ITS_KAPTON(POLYCH2)$");
  TGeoMedium *medAluminum = mgr->GetMedium("ITS_ALUMINUM$");

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSModulePattern(), fLayerNumber);
  TGeoVolume *modVol = new TGeoVolume(volname, module, medAir);

  TGeoVolume *kapCableVol = new TGeoVolume("FPCKapton", kapCable, medKapton);
  kapCableVol->SetLineColor(kBlue);
  kapCableVol->SetFillColor(kBlue);

  TGeoVolume *aluCableVol = new TGeoVolume("FPCAluminum",
					   aluCable, medAluminum);
  aluCableVol->SetLineColor(kCyan);
  aluCableVol->SetFillColor(kCyan);

  // build up the module
  ypos = -ytot + ymod; // = 0 if not kIBModel4
  for (Int_t j=0; j<fgkIBChipsPerRow; j++) {
    zpos = -zmod + j*2*zchip + zchip;
    modVol->AddNode(chipVol, j, new TGeoTranslation(0, ypos, zpos));
    fHierarchy[kChip]++;
  }

  if (fStaveModel == AliITSUv2::kIBModel4) {
    ypos += (ymod + aluCable->GetDY());
    if (fBuildLevel < 1)   // Aluminum
      modVol->AddNode(aluCableVol, 1, new TGeoTranslation(0, ypos, 0));

    ypos += (aluCable->GetDY() + kapCable->GetDY());
    if (fBuildLevel < 4)   // Kapton
      modVol->AddNode(kapCableVol, 1, new TGeoTranslation(0, ypos, 0));
  }

  // Done, return the module
  return modVol;
}

//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateStaveStructInnerB(const Double_t xsta,
						    const Double_t zsta,
						    const TGeoManager *mgr){
//
// Create the mechanical stave structure
//
// Input:
//         xsta : X length
//         zsta : Z length
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      22 Mar 2013  Chinorat Kobdaj
// Updated:      26 Apr 2013  Mario Sitta
//

  TGeoVolume *mechStavVol = 0;

  switch (fStaveModel) {
    case AliITSUv2::kIBModelDummy:
      mechStavVol = CreateStaveModelInnerBDummy(xsta,zsta,mgr);
      break;
    case AliITSUv2::kIBModel0:
      mechStavVol = CreateStaveModelInnerB0(xsta,zsta,mgr);
      break;
    case AliITSUv2::kIBModel1:
      mechStavVol = CreateStaveModelInnerB1(xsta,zsta,mgr);
      break;
    case AliITSUv2::kIBModel21:
      mechStavVol = CreateStaveModelInnerB21(xsta,zsta,mgr);
      break;
    case AliITSUv2::kIBModel22:
      mechStavVol = CreateStaveModelInnerB22(xsta,zsta,mgr);
      break;
    case AliITSUv2::kIBModel3:
      mechStavVol = CreateStaveModelInnerB3(xsta,zsta,mgr);
      break;
    case AliITSUv2::kIBModel4:
      mechStavVol = CreateStaveModelInnerB4(xsta,zsta,mgr);
      break;
    default:
      AliFatal(Form("Unknown stave model %d",fStaveModel));
      break;
  }

  return mechStavVol; 
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateStaveModelInnerBDummy(const Double_t ,
							const Double_t ,
							const TGeoManager *) const {
//
// Create dummy stave
//
// Input:
//         xsta : X length
//         zsta : Z length
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      22 Mar 2013  Chinorat Kobdaj
// Updated:      26 Apr 2013  Mario Sitta
//

  // Done, return the stave structur
  return 0;
}

//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateStaveModelInnerB0(const Double_t xsta,
						    const Double_t zsta,
						    const TGeoManager *mgr){
//
// Create the mechanical stave structure for Model 0 of TDR
//
// Input:
//         xsta : X length
//         zsta : Z length
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      22 Mar 2013  Chinorat Kobdaj
// Updated:      26 Apr 2013  Mario Sitta
//
  
  // Materials defined in AliITSUv2
  TGeoMedium *medAir    = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medWater  = mgr->GetMedium("ITS_WATER$");

  TGeoMedium *medM60J3K    = mgr->GetMedium("ITS_M60J3K$"); 
  TGeoMedium *medKapton    = mgr->GetMedium("ITS_KAPTON(POLYCH2)$");
  TGeoMedium *medGlue      = mgr->GetMedium("ITS_GLUE$");
  TGeoMedium *medFlexCable = mgr->GetMedium("ITS_FLEXCABLE$");

  // Local parameters
  Double_t kConeOutRadius = 0.15/2;
  Double_t kConeInRadius = 0.1430/2;
  Double_t kStaveLength = zsta*2;
  Double_t kStaveWidth = xsta*2-kConeOutRadius*2;
  Double_t kWidth = kStaveWidth/4;//1/2 of kWidth
  Double_t kStaveHeight = 0.3;
  Double_t kHeight = kStaveHeight/2;
  Double_t kAlpha = 90-67;//90-33.69;
  Double_t kTheta = kAlpha*TMath::DegToRad();
  Double_t kS1 = kWidth/TMath::Sin(kTheta);
  Double_t kL1 = kWidth/TMath::Tan(kTheta);
  Double_t kS2 = TMath::Sqrt(kHeight*kHeight + kS1*kS1);//TMath::Sin(the2);
  Double_t kThe2 = TMath::ATan(kHeight/kS1);
  Double_t kBeta = kThe2*TMath::RadToDeg();
  // Int_t  loop = kStaveLength/(kL1);
  // Double_t s3 = kWidth/(2*TMath::Sin(kTheta));
  // Double_t s4 = 3*kWidth/(2*TMath::Sin(kTheta));

  AliDebug(1, Form("BuildLevel %d\n",fBuildLevel));

  char volname[30];
  snprintf(volname, 30, "%s%d_StaveStruct", AliITSUGeomTGeo::GetITSStavePattern(), fLayerNumber);

  Double_t z=0, y=-0.011+0.0150, x=0;

   TGeoVolume *mechStavVol = 0;

  if (fBuildLevel < 5) {

    // world (trapezoid)
    TGeoXtru *mechStruct = new TGeoXtru(2); //z sections
    Double_t xv[5] = {kStaveWidth/2+0.1,kStaveWidth/2+0.1,0,-kStaveWidth/2-0.1,-kStaveWidth/2-0.1};
    Double_t yv[5] = {-kConeOutRadius*2-0.07,0,kStaveHeight,0,-kConeOutRadius*2-0.07};    
    mechStruct->DefinePolygon(5,xv,yv);
    mechStruct->DefineSection(0,-kStaveLength-0.1,0,0,1.);
    mechStruct->DefineSection(1,kStaveLength+0.1,0,0,1.);

    mechStavVol = new TGeoVolume(volname, mechStruct, medAir);
    mechStavVol->SetLineColor(12);
    mechStavVol->SetFillColor(12); 
    mechStavVol->SetVisibility(kTRUE);
      
    // detailed structure ++++++++++++++
    //Pipe Kapton grey-35
    TGeoTube *coolTube = new TGeoTube(kConeInRadius,kConeOutRadius,kStaveLength/2);
    TGeoVolume *volCoolTube= new TGeoVolume("pipe", coolTube, medKapton);
    volCoolTube->SetFillColor(35);
    volCoolTube->SetLineColor(35);
    mechStavVol->AddNode(volCoolTube,0,new TGeoTranslation(x+(kStaveWidth/2),y-(kHeight-kConeOutRadius),0));
    mechStavVol->AddNode(volCoolTube,1,new TGeoTranslation(x-(kStaveWidth/2),y-(kHeight-kConeOutRadius),0));
  }

  if (fBuildLevel < 4) {
    TGeoTube *coolTubeW = new TGeoTube(0.,kConeInRadius,kStaveLength/2);
    TGeoVolume *volCoolTubeW= new TGeoVolume("pipeWater", coolTubeW, medWater);
    volCoolTubeW->SetFillColor(4);
    volCoolTubeW->SetLineColor(4);
    mechStavVol->AddNode(volCoolTubeW,0,new TGeoTranslation(x+(kStaveWidth/2),y-(kHeight-kConeOutRadius),0));
    mechStavVol->AddNode(volCoolTubeW,1,new TGeoTranslation(x-(kStaveWidth/2),y-(kHeight-kConeOutRadius),0));
  }

  //frequency of filament
  //n = 4 means very dense(4 filaments per interval)
  //n = 2 means dense(2 filaments per interval)
  Int_t n =4;
  Int_t loop = (Int_t)(kStaveLength/(4*kL1/n) + 2/n)-1;
  if (fBuildLevel < 3) {
    //Top CFRP Filament black-12 Carbon structure TGeoBBox (length,thickness,width)
    TGeoBBox *t2=new TGeoBBox(kS2,0.007/2,0.15/2);//(kS2,0.002,0.02);
    TGeoVolume *volT2=new TGeoVolume("TopFilament", t2, medM60J3K);
    volT2->SetLineColor(12);
    volT2->SetFillColor(12); 

    for(int i=1;i<loop;i++){  //i<60;i++){
      mechStavVol->AddNode(volT2,4*i+0,
				  new TGeoCombiTrans(x+kWidth,y+(2*kConeOutRadius),z-kStaveLength/2+(i*(4/n)*kL1)+kS1/2,//z-14.25+(i*2*kL1),
						     new TGeoRotation("volT2",90,90-kAlpha,90-kBeta)));
      mechStavVol->AddNode(volT2,4*i+1,
				  new TGeoCombiTrans(x-kWidth,y+(2*kConeOutRadius),z-kStaveLength/2+(i*(4/n)*kL1)+kS1/2,//z-14.25+(i*2*kL1),
						     new TGeoRotation("volT2",90,-90+kAlpha,-90+kBeta)));
      mechStavVol->AddNode(volT2,4*i+2,
				  new TGeoCombiTrans(x+kWidth,y+(2*kConeOutRadius),z-kStaveLength/2+(i*(4/n)*kL1)+kS1/2,//z-14.25+(i*2*kL1),
						     new TGeoRotation("volT2",90,-90+kAlpha,90-kBeta)));
      mechStavVol->AddNode(volT2,4*i+3,
				  new TGeoCombiTrans(x-kWidth,y+(2*kConeOutRadius),z-kStaveLength/2+(i*(4/n)*kL1)+kS1/2,//z-14.25+(i*2*kL1),  
						     new TGeoRotation("volT2",90,90-kAlpha,-90+kBeta)));
    }


    //Bottom CFRP Filament black-12 Carbon structure  TGeoBBox (thickness,width,length)
    TGeoBBox *t1=new TGeoBBox(0.007/2,0.15/2,kS1);//(0.002,0.02,kS1);
    TGeoVolume *volT1=new TGeoVolume("CFRPBottom", t1, medM60J3K);
    volT1->SetLineColor(12);
    volT1->SetFillColor(12); 

    for(int i=1;i<loop;i++){
      mechStavVol->AddNode(volT1,4*i+0,
				  new TGeoCombiTrans(x+kWidth,y-kHeight,z-kStaveLength/2+((4/n)*kL1*i)+kS1/2, //z-14.25+(i*2*kL1),  
						     new TGeoRotation("volT1",-90,kAlpha,0)));
      mechStavVol->AddNode(volT1,4*i+1,
				  new TGeoCombiTrans(x-kWidth,y-kHeight,z-kStaveLength/2+((4/n)*kL1*i)+kS1/2,  //z-14.25+(i*2*kL1), 
						     new TGeoRotation("volT1",90,kAlpha,0)));
      mechStavVol->AddNode(volT1,4*i+2,
				  new TGeoCombiTrans(x+kWidth,y-kHeight,z-kStaveLength/2+(i*(4/n)*kL1)+kS1/2, //z-14.25+(i*2*kL1), 
						     new TGeoRotation("volT1",-90,-kAlpha,0)));
      mechStavVol->AddNode(volT1,4*i+3,
				  new TGeoCombiTrans(x-kWidth,y-kHeight,z-kStaveLength/2+(i*(4/n)*kL1)+kS1/2, //z-14.25+(i*2*kL1), 
						     new TGeoRotation("volT1",-90,+kAlpha,0)));
    }
  }
   
  if (fBuildLevel < 2) {
    // Glue CFRP-Silicon layers TGeoBBox(thickness,width,kS1);
    TGeoBBox *tG=new TGeoBBox(0.0075/2,0.18/2,kS1);
    TGeoVolume *volTG=new TGeoVolume("Glue1", tG, medGlue);
    volTG->SetLineColor(5);
    volTG->SetFillColor(5); 

    for(int i=1;i<loop;i++){ //i<60;i++){
      mechStavVol->AddNode(volTG,4*i+0,
				  new TGeoCombiTrans(x+kWidth,y-0.16,z-kStaveLength/2+((4/n)*kL1*i)+kS1/2, //z-14.25+(2*kL1*i), 
						     new TGeoRotation("volTG",-90,kAlpha,0)));
      mechStavVol->AddNode(volTG,4*i+1,
				  new TGeoCombiTrans(x-kWidth,y-0.16,z-kStaveLength/2+((4/n)*kL1*i)+kS1/2, //z-14.25+(2*kL1*i), 
						     new TGeoRotation("volTG",90,kAlpha,0)));
      mechStavVol->AddNode(volTG,4*i+2,
				  new TGeoCombiTrans(x+kWidth,y-0.16,z-kStaveLength/2+((4/n)*i*kL1)+kS1/2, //z-14.25+(i*2*kL1), 
						     new TGeoRotation("volTG",-90,-kAlpha,0)));
      mechStavVol->AddNode(volTG,4*i+3,
				  new TGeoCombiTrans(x-kWidth,y-0.16,z-kStaveLength/2+(i*(4/n)*kL1)+kS1/2, //z-14.25+(i*2*kL1), 
						     new TGeoRotation("volTG",-90,+kAlpha,0)));
    }

    TGeoBBox *glue = new TGeoBBox(xsta, 0.005/2, zsta);
    TGeoVolume *volGlue=new TGeoVolume("Glue2", glue, medGlue);
    volGlue->SetLineColor(5);
    volGlue->SetFillColor(5); 
    //mechStavVol->AddNode(volGlue, 0, new TGeoCombiTrans(x, y-0.16, z, new TGeoRotation("",0, 0, 0)));
    mechStavVol->AddNode(volGlue, 1, new TGeoCombiTrans(x, y-0.165-fSensorThick-0.005, z, new TGeoRotation("",0, 0, 0)));
  }

  if (fBuildLevel < 1) {
    //Flex cable brown-28 TGeoBBox(width,thickness,length); 
    TGeoBBox *kapCable = new TGeoBBox(xsta, 0.01/2, zsta);
    TGeoVolume *volCable=new TGeoVolume("FlexCable", kapCable, medFlexCable);
    volCable->SetLineColor(28);
    volCable->SetFillColor(28); 
    mechStavVol->AddNode(volCable, 0, new TGeoCombiTrans(x, y-0.165-fSensorThick-0.005-0.01, z, new TGeoRotation("",0, 0, 0)));
 }

  // Done, return the stave structur
  return mechStavVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateStaveModelInnerB1(const Double_t xsta,
						    const Double_t zsta,
						    const TGeoManager *mgr){
//
// Create the mechanical stave structure for Model 1 of TDR
//
// Input:
//         xsta : X length
//         zsta : Z length
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      22 Mar 2013  Chinorat Kobdaj
// Updated:      26 Apr 2013  Mario Sitta
//
  
  // Materials defined in AliITSUv2
  TGeoMedium *medAir    = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medWater  = mgr->GetMedium("ITS_WATER$");

  TGeoMedium *medM60J3K    = mgr->GetMedium("ITS_M60J3K$"); 
  TGeoMedium *medKapton    = mgr->GetMedium("ITS_KAPTON(POLYCH2)$");
  TGeoMedium *medGlue      = mgr->GetMedium("ITS_GLUE$");
  TGeoMedium *medFlexCable = mgr->GetMedium("ITS_FLEXCABLE$");

  // Local parameters
  Double_t kConeOutRadius = 0.15/2;
  //    Double_t kConeInRadius = 0.1430/2;
  Double_t kStaveLength = zsta*2;
  //    Double_t kStaveWidth = xsta*2-kConeOutRadius*2;
  Double_t kStaveWidth = xsta*2;
  Double_t kWidth = kStaveWidth/4;//1/2 of kWidth
  Double_t kStaveHeight = 0.3;
  Double_t kHeight = kStaveHeight/2;
  Double_t kAlpha = 90-33.;//90-30;
  Double_t kTheta = kAlpha*TMath::DegToRad();
  Double_t kS1 = kWidth/TMath::Sin(kTheta);
  Double_t kL1 = kWidth/TMath::Tan(kTheta);
  Double_t kS2 = TMath::Sqrt(kHeight*kHeight + kS1*kS1);//TMath::Sin(the2);
  Double_t kThe2 = TMath::ATan(kHeight/kS1);
  Double_t kBeta = kThe2*TMath::RadToDeg();
  Int_t  loop = (Int_t)((kStaveLength/(2*kL1))/2);
  

  TGeoVolume *mechStavVol = 0;

  char volname[30];
  snprintf(volname, 30, "%s%d_StaveStruct", AliITSUGeomTGeo::GetITSStavePattern(), fLayerNumber);
    

  // detailed structure ++++++++++++++
  Double_t z=0, y=-0.011+0.0150, x=0;

  // Polimide micro channels numbers
  Double_t yMC = y-kHeight+0.01;
  Int_t nb = (Int_t)(kStaveWidth/0.1)+1;
  Double_t xstaMC = (nb*0.1-0.08)/2;


  if (fBuildLevel < 5) {
    // world (trapezoid)
    TGeoXtru *mechStruct = new TGeoXtru(2); //z sections
    Double_t xv[5] = {kStaveWidth/2+0.1,kStaveWidth/2+0.1,0,-kStaveWidth/2-0.1,-kStaveWidth/2-0.1};
    Double_t yv[5] = {-kConeOutRadius*2-0.07,0,kStaveHeight,0,-kConeOutRadius*2-0.07};    
    mechStruct->DefinePolygon(5,xv,yv);
    mechStruct->DefineSection(0,-kStaveLength-0.1,0,0,1.);
    mechStruct->DefineSection(1,kStaveLength+0.1,0,0,1.);

    mechStavVol = new TGeoVolume(volname, mechStruct, medAir);
    mechStavVol->SetLineColor(12);
    mechStavVol->SetFillColor(12); 
    mechStavVol->SetVisibility(kTRUE);
      
    // Polimide micro channels numbers
    TGeoBBox *tM0=new TGeoBBox(xstaMC, 0.005/2, zsta);
    TGeoVolume *volTM0=new TGeoVolume("MicroChanCover", tM0, medKapton);
    volTM0->SetLineColor(35);
    volTM0->SetFillColor(35); 
    mechStavVol->AddNode(volTM0, 0, new TGeoCombiTrans(x,-0.0125+yMC, z, new TGeoRotation("",0, 0, 0)));
    mechStavVol->AddNode(volTM0, 1, new TGeoCombiTrans(x,+0.0125+yMC, z, new TGeoRotation("",0, 0, 0)));
      
    TGeoBBox *tM0b=new TGeoBBox(0.02/2, 0.02/2, zsta);
    TGeoVolume *volTM0b=new TGeoVolume("MicroChanWalls", tM0b, medKapton);
    volTM0b->SetLineColor(35);
    volTM0b->SetFillColor(35); 
    for (Int_t ib=0;ib<nb;ib++) {
      mechStavVol->AddNode(volTM0b, ib, new TGeoCombiTrans(x+ib*0.1-xstaMC+0.01,yMC, z, new TGeoRotation("",0, 0, 0)));
    }
      
  }
    
  if (fBuildLevel < 4) {
    // Water in Polimide micro channels
    TGeoBBox *water=new TGeoBBox(0.08/2, 0.02/2, zsta+0.1);
    TGeoVolume *volWater=new TGeoVolume("Water", water, medWater);
    volWater->SetLineColor(4);
    volWater->SetFillColor(4); 
    for (Int_t ib=0;ib<(nb-1);ib++) {
      mechStavVol->AddNode(volWater, ib, new TGeoCombiTrans(x+ib*0.1-xstaMC+0.06,yMC, z, new TGeoRotation("",0, 0, 0)));
    }
  }
    
  if (fBuildLevel < 3) {
    //Bottom filament CFRP black-12 Carbon structure TGeoBBox (thickness,width,length)
    Double_t filWidth = 0.04;
    Double_t filHeight= 0.02;
    TGeoBBox *t1=new TGeoBBox(filHeight/2,filWidth/2,kS1);
    TGeoVolume *volT1=new TGeoVolume("CFRPBottom", t1, medM60J3K);
    volT1->SetLineColor(12);
    volT1->SetFillColor(12); 
    for(int i=0;i<loop;i++){//i<30;i++){
      mechStavVol->AddNode(volT1,4*i+0,
				  new TGeoCombiTrans(x+kWidth,y-kHeight+0.04+filHeight/2,z-kStaveLength/2+(4*kL1)+kS1/2, 
						     new TGeoRotation("volT1",-90,kAlpha,0)));
      mechStavVol->AddNode(volT1,4*i+1,
				  new TGeoCombiTrans(x-kWidth,y-kHeight+0.04+filHeight/2,z-kStaveLength/2+(4*kL1*i)+kS1/2, 
						     new TGeoRotation("volT1",90,kAlpha,0)));
      mechStavVol->AddNode(volT1,4*i+2,
				  new TGeoCombiTrans(x+kWidth,y-kHeight+0.04+filHeight/2,z-kStaveLength/2+2*kL1+(i*4*kL1)+kS1/2, 
						     new TGeoRotation("volT1",-90,-kAlpha,0)));
      mechStavVol->AddNode(volT1,4*i+3,
				  new TGeoCombiTrans(x-kWidth,y-kHeight+0.04+filHeight/2,z-kStaveLength/2+2*kL1+(i*4*kL1)+kS1/2, 
						     new TGeoRotation("volT1",-90,+kAlpha,0)));
    }

      // Top filament CFRP black-12 Carbon structure TGeoBBox (length,thickness,width)
    TGeoBBox *t2=new TGeoBBox(kS2,filHeight/2,filWidth/2);
    TGeoVolume *volT2=new TGeoVolume("CFRPTop", t2, medM60J3K);
    volT2->SetLineColor(12);
    volT2->SetFillColor(12); 
    for(int i=0;i<loop;i++){ //i<30;i++){
      mechStavVol->AddNode(volT2,4*i+0,
				  new TGeoCombiTrans(x+kWidth,y+0.04+filHeight/2,z-kStaveLength/2+(i*4*kL1)+kS1/2,
						     new TGeoRotation("volT2",90,90-kAlpha,90-kBeta)));
      mechStavVol->AddNode(volT2,4*i+1,
				  new TGeoCombiTrans(x-kWidth,y+0.04+filHeight/2,z-kStaveLength/2+(i*4*kL1)+kS1/2,
						     new TGeoRotation("volT2",90,-90+kAlpha,-90+kBeta)));
      mechStavVol->AddNode(volT2,4*i+2,
				  new TGeoCombiTrans(x+kWidth,y+0.04+filHeight/2,z-kStaveLength/2+2*kL1+(i*4*kL1)+kS1/2,
						     new TGeoRotation("volT2",90,-90+kAlpha,90-kBeta)));
      mechStavVol->AddNode(volT2,4*i+3,
				  new TGeoCombiTrans(x-kWidth,y+0.04+filHeight/2,z-kStaveLength/2+2*kL1+(i*4*kL1)+kS1/2, 
						     new TGeoRotation("volT2",90,90-kAlpha,-90+kBeta)));
    }
  }

  if (fBuildLevel < 2) {
    // Glue between filament and polimide micro channel
    TGeoBBox *t3=new TGeoBBox(0.01/2,0.04,kS1);
    TGeoVolume *volT3=new TGeoVolume("FilamentGlue", t3, medGlue);
    volT3->SetLineColor(5);
    volT3->SetFillColor(5); 
    for(int i=0;i<loop;i++){//i<30;i++){
      mechStavVol->AddNode(volT3,4*i+0,
				  new TGeoCombiTrans(x+kWidth,y-kHeight+0.0325,z-kStaveLength/2+(4*kL1*i)+kS1/2, 
						     new TGeoRotation("volT1",-90,kAlpha,0)));
      mechStavVol->AddNode(volT3,4*i+1,
				  new TGeoCombiTrans(x-kWidth,y-kHeight+0.0325,z-kStaveLength/2+(4*kL1*i)+kS1/2, 
						     new TGeoRotation("volT1",90,kAlpha,0)));
      mechStavVol->AddNode(volT3,4*i+2,
				  new TGeoCombiTrans(x+kWidth,y-kHeight+0.0325,z-kStaveLength/2+2*kL1+(i*4*kL1)+kS1/2, 
						     new TGeoRotation("volT1",-90,-kAlpha,0)));
      mechStavVol->AddNode(volT3,4*i+3,
				  new TGeoCombiTrans(x-kWidth,y-kHeight+0.0325,z-kStaveLength/2+2*kL1+(i*4*kL1)+kS1/2, 
						     new TGeoRotation("volT1",-90,+kAlpha,0)));
    }
      
    // Glue microchannel and sensor
    TGeoBBox *glueM = new TGeoBBox(xsta, 0.01/2, zsta);
    TGeoVolume *volGlueM=new TGeoVolume("MicroChanGlue", glueM, medGlue);
    volGlueM->SetLineColor(5);
    volGlueM->SetFillColor(5); 
    mechStavVol->AddNode(volGlueM, 0, new TGeoCombiTrans(x, y-0.16, z, new TGeoRotation("",0, 0, 0)));

    // Glue sensor and kapton
    TGeoBBox *glue = new TGeoBBox(xsta, 0.005/2, zsta);
    TGeoVolume *volGlue=new TGeoVolume("SensorGlue", glue, medGlue);
    volGlue->SetLineColor(5);
    volGlue->SetFillColor(5); 
    mechStavVol->AddNode(volGlue, 1, new TGeoCombiTrans(x, y-0.165-fSensorThick-0.005, z, new TGeoRotation("",0, 0, 0)));
  }

  if (fBuildLevel < 1) {
    TGeoBBox *kapCable = new TGeoBBox(xsta, 0.01/2, zsta);
    TGeoVolume *volCable=new TGeoVolume("FlexCable", kapCable, medFlexCable);
    volCable->SetLineColor(28);
    volCable->SetFillColor(28); 
    mechStavVol->AddNode(volCable, 0, new TGeoCombiTrans(x, y-0.165-fSensorThick-0.005-0.01, z, new TGeoRotation("",0, 0, 0)));
  }
    
  // Done, return the stave structur
  return mechStavVol;

}

//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateStaveModelInnerB21(const Double_t xsta,
						     const Double_t zsta,
						     const TGeoManager *mgr){
//
// Create the mechanical stave structure for Model 2.1 of TDR
//
// Input:
//         xsta : X length
//         zsta : Z length
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      22 Mar 2013  Chinorat Kobdaj
// Updated:      26 Apr 2013  Mario Sitta
//
  
  // Materials defined in AliITSUv2
  TGeoMedium *medAir    = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medWater  = mgr->GetMedium("ITS_WATER$");

  TGeoMedium *medM60J3K    = mgr->GetMedium("ITS_M60J3K$"); 
  TGeoMedium *medKapton    = mgr->GetMedium("ITS_KAPTON(POLYCH2)$");
  TGeoMedium *medGlue      = mgr->GetMedium("ITS_GLUE$");
  TGeoMedium *medFlexCable = mgr->GetMedium("ITS_FLEXCABLE$");
  TGeoMedium *medK13D2U2k  = mgr->GetMedium("ITS_K13D2U2k$");
  TGeoMedium *medFGS003    = mgr->GetMedium("ITS_FGS003$"); 
  TGeoMedium *medCarbonFleece = mgr->GetMedium("ITS_CarbonFleece$"); 

  // Local parameters
  Double_t kConeOutRadius =0.151384/2;
  Double_t kConeInRadius = 0.145034/2;
  Double_t kStaveLength = zsta;
  Double_t kStaveWidth = xsta*2;
  Double_t kWidth = (kStaveWidth+0.005)/4;
  Double_t kStaveHeigth = 0.33;//0.33;
  Double_t kHeight = (kStaveHeigth+0.025)/2;
  Double_t kAlpha = 57; //56.31;
  Double_t kTheta = kAlpha*TMath::DegToRad();
  Double_t kS1 = (kStaveWidth/4)/TMath::Sin(kTheta);
  Double_t kL1 = (kStaveWidth/4)/TMath::Tan(kTheta);
  Double_t kS2 = sqrt(kHeight*kHeight + kS1*kS1);//TMath::Sin(the2);
  Double_t kThe2 = TMath::ATan(kHeight/kS1);
  Double_t kBeta = kThe2*TMath::RadToDeg();
  // Double_t lay1 = 0.003157;
  Double_t kLay1 = 0.003;//Amec carbon
  // Double_t lay2 = 0.0043215;//C Fleece carbon
  Double_t kLay2 = 0.002;//C Fleece carbon
  Double_t kLay3 = 0.007;//K13D2U carbon
  Int_t  loop = (Int_t)(kStaveLength/(2*kL1));


  char volname[30];
  snprintf(volname, 30, "%s%d_StaveStruct", AliITSUGeomTGeo::GetITSStavePattern(), fLayerNumber);

  Double_t z=0, y=-(kConeOutRadius+0.03)+0.0385, x=0;

  TGeoVolume *mechStavVol = 0;

  if (fBuildLevel < 5) {
    // world (trapezoid)
    TGeoXtru *mechStruct = new TGeoXtru(2); //z sections
    Double_t xv[5] = {kStaveWidth/2+0.1,kStaveWidth/2+0.1,0,-kStaveWidth/2-0.1,-kStaveWidth/2-0.1};
    Double_t yv[5] = {-kConeOutRadius*2-0.07,0,kStaveHeigth,0,-kConeOutRadius*2-0.07};    
    mechStruct->DefinePolygon(5,xv,yv);
    mechStruct->DefineSection(0,-kStaveLength-0.1,0,0,1.);
    mechStruct->DefineSection(1,kStaveLength+0.1,0,0,1.);

    mechStavVol = new TGeoVolume(volname, mechStruct, medAir);
    mechStavVol->SetLineColor(12);
    mechStavVol->SetFillColor(12); 
    mechStavVol->SetVisibility(kTRUE);  
      
    //Pipe Kapton grey-35 
    TGeoCone *cone1 = new TGeoCone(kStaveLength,kConeInRadius,kConeOutRadius,kConeInRadius,kConeOutRadius);
    TGeoVolume *volCone1= new TGeoVolume("PolyimidePipe", cone1, medKapton);
    volCone1->SetFillColor(35);
    volCone1->SetLineColor(35);
    mechStavVol->AddNode(volCone1,1,new TGeoTranslation(x+0.25,y,z));
    mechStavVol->AddNode(volCone1,2,new TGeoTranslation(x-0.25,y,z));
  }

  if (fBuildLevel < 4) {
    
    TGeoTube *coolTubeW = new TGeoTube(0.,kConeInRadius,kStaveLength);
    TGeoVolume *volCoolTubeW= new TGeoVolume("Water", coolTubeW, medWater);
    volCoolTubeW->SetFillColor(4);
    volCoolTubeW->SetLineColor(4);
    mechStavVol->AddNode(volCoolTubeW,0,new TGeoTranslation(x-0.25,y,z));
    mechStavVol->AddNode(volCoolTubeW,1,new TGeoTranslation(x+0.25,y,z));
  }

  if (fBuildLevel < 3) {
    //top fillament
    // Top filament M60J black-12 Carbon structure TGeoBBox (length,thickness,width)
    TGeoBBox *t2=new TGeoBBox(kS2,0.02/2,0.04/2); //TGeoBBox *t2=new TGeoBBox(kS2,0.01,0.02);
    TGeoVolume *volT2=new TGeoVolume("TopFilament", t2, medM60J3K);
    volT2->SetLineColor(12);
    volT2->SetFillColor(12); 
    for(int i=0;i<loop;i++){// i<28;i++){
      mechStavVol->AddNode(volT2,i*4+1,new TGeoCombiTrans(x+kWidth,y+kHeight+(0.12/2)-0.014+0.007,z-kStaveLength+(i*4*kL1)+kS1/2, new TGeoRotation("volT2",90,90-kAlpha,90-kBeta)));
      mechStavVol->AddNode(volT2,i*4+2,new TGeoCombiTrans(x-kWidth,y+kHeight+(0.12/2)-0.014+0.007,z-kStaveLength+(i*4*kL1)+kS1/2, new TGeoRotation("volT2",90,-90+kAlpha,-90+kBeta)));
      mechStavVol->AddNode(volT2,i*4+3,new TGeoCombiTrans(x+kWidth,y+kHeight+(0.12/2)-0.014+0.007,z-kStaveLength+2*kL1+(i*4*kL1)+kS1/2, new TGeoRotation("volT2",90,-90+kAlpha,90-kBeta)));
      mechStavVol->AddNode(volT2,i*4+4,new TGeoCombiTrans(x-kWidth,y+kHeight+(0.12/2)-0.014+0.007,z-kStaveLength+2*kL1+(i*4*kL1)+kS1/2, new TGeoRotation("volT2",90,90-kAlpha,-90+kBeta)));
//    mechStavVol->AddNode(volT2,i*4+1,new TGeoCombiTrans(x+kWidth+0.0036,y+kHeight-(0.12/2)+0.072,z+kStaveLength+(i*4*kL1)+kS1/2, new TGeoRotation("volT2",90,90-kAlpha,90-kBeta)));

 }
 
    //wall side structure out
    TGeoBBox *box4 = new TGeoBBox(0.03/2,0.12/2,kStaveLength-0.50);
    TGeoVolume *plate4 = new TGeoVolume("WallOut",box4,medM60J3K);
    plate4->SetFillColor(35);
    plate4->SetLineColor(35);
    mechStavVol->AddNode(plate4,1,new TGeoCombiTrans(x+(2*kStaveWidth/4)-(0.03/2),y-0.0022-kConeOutRadius+0.12/2+0.007,z,new TGeoRotation("plate4",0,0,0)));
    mechStavVol->AddNode(plate4,2,new TGeoCombiTrans(x-(2*kStaveWidth/4)+(0.03/2),y-0.0022-kConeOutRadius+0.12/2+0.007,z,new TGeoRotation("plate4",0,0,0)));
    //wall side in
    TGeoBBox *box5 = new TGeoBBox(0.015/2,0.12/2,kStaveLength-0.50);
    TGeoVolume *plate5 = new TGeoVolume("WallIn",box5,medM60J3K);
    plate5->SetFillColor(12);
    plate5->SetLineColor(12);
    mechStavVol->AddNode(plate5,1,new TGeoCombiTrans(x+(2*kStaveWidth/4)-0.03-0.015/2,y-0.0022-kConeOutRadius+0.12/2+0.007,z,new TGeoRotation("plate5",0,0,0)));
    mechStavVol->AddNode(plate5,2,new TGeoCombiTrans(x-(2*kStaveWidth/4)+0.03+0.015/2,y-0.0022-kConeOutRadius+0.12/2+0.007,z,new TGeoRotation("plate5",0,0,0)));

     //Amec Thermasol red-2 cover tube FGS300
    TGeoConeSeg *cons1 = new TGeoConeSeg(kStaveLength-0.50,kConeOutRadius,kConeOutRadius+kLay1,kConeOutRadius,kConeOutRadius+kLay1,0,180);
    TGeoVolume *cone11 = new TGeoVolume("ThermasolPipeCover",cons1,medFGS003);
    cone11->SetFillColor(2);
    cone11->SetLineColor(2);
    mechStavVol->AddNode(cone11,1,new TGeoCombiTrans(x+0.25,y,z,new TGeoRotation("Cone11",0,0,0)));
    mechStavVol->AddNode(cone11,2,new TGeoCombiTrans(x-0.25,y,z,new TGeoRotation("Cone11",0,0,0)));

    TGeoBBox *box2 = new TGeoBBox((0.50-(2*kConeOutRadius))/2,kLay1/2,kStaveLength-0.50);
    TGeoVolume *plate2 = new TGeoVolume("ThermasolMiddle",box2,medFGS003);
    plate2->SetFillColor(2);
    plate2->SetLineColor(2);
    mechStavVol->AddNode(plate2,1,new TGeoCombiTrans(x,y-kConeOutRadius+(kLay1/2),z,new TGeoRotation("plate2",0,0,0)));

    TGeoBBox *box21 = new TGeoBBox((0.75-0.25-kConeOutRadius-kLay1)/2,kLay1/2,kStaveLength-0.50);
    TGeoVolume *plate21 = new TGeoVolume("ThermasolLeftRight",box21,medFGS003);
    plate21->SetFillColor(2);
    plate21->SetLineColor(2);
    mechStavVol->AddNode(plate21,1,new TGeoCombiTrans(x+0.25+kConeOutRadius+(0.75-0.25-kConeOutRadius)/2-(kLay1/2),y-kConeOutRadius+(kLay1/2),z,new TGeoRotation("plate21",0,0,0)));
    mechStavVol->AddNode(plate21,2,new TGeoCombiTrans(x-0.25-kConeOutRadius-(0.75-0.25-kConeOutRadius)/2+(kLay1/2),y-kConeOutRadius+(kLay1/2),z,new TGeoRotation("plate21",0,0,0)));

    TGeoBBox *box22 = new TGeoBBox((kLay1/2),kConeOutRadius/2,kStaveLength-0.50);
    TGeoVolume *plate22 = new TGeoVolume("ThermasolVertical",box22,medFGS003);
    plate22->SetFillColor(2);
    plate22->SetLineColor(2);
    mechStavVol->AddNode(plate22,1,new TGeoCombiTrans(x+0.25+kConeOutRadius+(kLay1/2),y-kConeOutRadius/2,z,new TGeoRotation("plate22",0,0,0)));
    mechStavVol->AddNode(plate22,2,new TGeoCombiTrans(x+0.25-kConeOutRadius-(kLay1/2),y-kConeOutRadius/2,z,new TGeoRotation("plate22",0,0,0)));
    mechStavVol->AddNode(plate22,3,new TGeoCombiTrans(x-0.25+kConeOutRadius+(kLay1/2),y-kConeOutRadius/2,z,new TGeoRotation("plate22",0,0,0)));
    mechStavVol->AddNode(plate22,4,new TGeoCombiTrans(x-0.25-kConeOutRadius-(kLay1/2),y-kConeOutRadius/2,z,new TGeoRotation("plate22",0,0,0)));

    //C Fleece
    TGeoConeSeg *cons2 = new TGeoConeSeg(kStaveLength-0.50,kConeOutRadius+kLay1,kConeOutRadius+kLay1+kLay2,kConeOutRadius+kLay1,kConeOutRadius+kLay1+kLay2,0,180); 
    TGeoVolume *cone12 = new TGeoVolume("CFleecePipeCover",cons2,medCarbonFleece);
    cone12->SetFillColor(28);
    cone12->SetLineColor(28);
    mechStavVol->AddNode(cone12,1,new TGeoCombiTrans(x+0.25,y,z,new TGeoRotation("Cone12",0,0,0)));
    mechStavVol->AddNode(cone12,2,new TGeoCombiTrans(x-0.25,y,z,new TGeoRotation("Cone12",0,0,0)));

    TGeoBBox *box3 = new TGeoBBox((0.50-(2*(kConeOutRadius+kLay1)))/2,kLay2/2,kStaveLength-0.50);
    TGeoVolume *plate3 = new TGeoVolume("CFleeceMiddle",box3,medCarbonFleece);
    plate3->SetFillColor(28);
    plate3->SetLineColor(28);
    mechStavVol->AddNode(plate3,1,new TGeoCombiTrans(x,y-kConeOutRadius+kLay1+(kLay2/2),z,new TGeoRotation("plate3",0,0,0)));

    TGeoBBox *box31 = new TGeoBBox((0.75-0.25-kConeOutRadius-kLay1)/2,kLay2/2,kStaveLength-0.50);
    TGeoVolume *plate31 = new TGeoVolume("CFleeceLeftRight",box31,medCarbonFleece);
    plate31->SetFillColor(28);
    plate31->SetLineColor(28);
    mechStavVol->AddNode(plate31,1,new TGeoCombiTrans(x+0.25+kConeOutRadius+kLay1+(0.75-0.25-kConeOutRadius-kLay1)/2,y-kConeOutRadius+kLay1+(kLay2/2),z,new TGeoRotation("plate31",0,0,0)));
    mechStavVol->AddNode(plate31,2,new TGeoCombiTrans(x-0.25-kConeOutRadius-kLay1-(0.75-0.25-kConeOutRadius-kLay1)/2,y-kConeOutRadius+kLay1+(kLay2/2),z,new TGeoRotation("plate31",0,0,0)));

    TGeoBBox *box32 = new TGeoBBox((kLay2/2),(kConeOutRadius-kLay1)/2,kStaveLength-0.50);
    TGeoVolume *plate32 = new TGeoVolume("CFleeceVertical",box32,medCarbonFleece);
    plate32->SetFillColor(28);
    plate32->SetLineColor(28);
    mechStavVol->AddNode(plate32,1,new TGeoCombiTrans(x+0.25+kConeOutRadius+kLay1+(kLay2/2),y+(kLay1-kConeOutRadius)/2,z,new TGeoRotation("plate32",0,0,0)));
    mechStavVol->AddNode(plate32,2,new TGeoCombiTrans(x+0.25-kConeOutRadius-kLay1-(kLay2/2),y+(kLay1-kConeOutRadius)/2,z,new TGeoRotation("plate32",0,0,0)));
    mechStavVol->AddNode(plate32,3,new TGeoCombiTrans(x-0.25+kConeOutRadius+kLay1+(kLay2/2),y+(kLay1-kConeOutRadius)/2,z,new TGeoRotation("plate32",0,0,0)));
    mechStavVol->AddNode(plate32,4,new TGeoCombiTrans(x-0.25-kConeOutRadius-kLay1-(kLay2/2),y+(kLay1-kConeOutRadius)/2,z,new TGeoRotation("plate32",0,0,0)));


    //K13D2U carbon plate
    TGeoBBox *box1 = new TGeoBBox(2*kWidth,kLay3/2,kStaveLength-0.50);
    TGeoVolume *plate1 = new TGeoVolume("CarbonPlate",box1,medK13D2U2k);
    plate1->SetFillColor(5);
    plate1->SetLineColor(5);
    mechStavVol->AddNode(plate1,1,new TGeoCombiTrans(x,y-(kConeOutRadius+(kLay3/2)),z,new TGeoRotation("plate1",0,0,0)));

    //C Fleece bottom plate 
    TGeoBBox *box6 = new TGeoBBox(2*kWidth,kLay2/2,kStaveLength-0.50);
    TGeoVolume *plate6 = new TGeoVolume("CFleeceBottom",box6,medCarbonFleece);
    plate6->SetFillColor(2);
    plate6->SetLineColor(2);
    mechStavVol->AddNode(plate6,1,new TGeoCombiTrans(x,y-(kConeOutRadius+kLay3+(kLay2/2)),z,new TGeoRotation("plate1",0,0,0)));
      
      
  }

  if (fBuildLevel < 2) {
    //Glue layers and kapton
    TGeoBBox *glue = new TGeoBBox(kStaveWidth/2, 0.005/2, zsta);
    TGeoVolume *volGlue=new TGeoVolume("Glue", glue, medGlue);
    volGlue->SetLineColor(5);
    volGlue->SetFillColor(5); 
    mechStavVol->AddNode(volGlue, 0, new TGeoCombiTrans(x,y-(kConeOutRadius+kLay3+(kLay2/2)+(0.01/2)), z, new TGeoRotation("",0, 0, 0)));
    mechStavVol->AddNode(volGlue, 1, new TGeoCombiTrans(x,y-(kConeOutRadius+kLay3+(kLay2/2)+0.01+fSensorThick+(0.01/2)), z, new TGeoRotation("",0, 0, 0)));
  }

  if (fBuildLevel < 1) {
    TGeoBBox *kapCable = new TGeoBBox(kStaveWidth/2, 0.01/2, zsta);
    TGeoVolume *volCable=new TGeoVolume("FlexCable", kapCable, medFlexCable);
    volCable->SetLineColor(28);
    volCable->SetFillColor(28); 
    mechStavVol->AddNode(volCable, 0, new TGeoCombiTrans(x, y-(kConeOutRadius+kLay3+(kLay2/2)+0.01+fSensorThick+0.01+(0.01/2)), z, new TGeoRotation("",0, 0, 0)));
  }
    

  // Done, return the stave structure
  return mechStavVol;
  
}

// new model22
//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateStaveModelInnerB22(const Double_t xsta,
						     const Double_t zsta,
						     const TGeoManager *mgr){
//
// Create the mechanical stave structure for Model 2.2 of TDR
//
// Input:
//         xsta : X length
//         zsta : Z length
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      22 Mar 2013  Chinorat Kobdaj
// Updated:      26 Apr 2013  Mario Sitta
// Updated:      30 Apr 2013  Wanchaloem Poonsawat 
//
  
  // Materials defined in AliITSUv2
  TGeoMedium *medAir    = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medWater  = mgr->GetMedium("ITS_WATER$");

  TGeoMedium *medM60J3K    = mgr->GetMedium("ITS_M60J3K$"); 
  TGeoMedium *medKapton    = mgr->GetMedium("ITS_KAPTON(POLYCH2)$");
  TGeoMedium *medGlue      = mgr->GetMedium("ITS_GLUE$");
  TGeoMedium *medFlexCable = mgr->GetMedium("ITS_FLEXCABLE$");
  TGeoMedium *medK13D2U2k  = mgr->GetMedium("ITS_K13D2U2k$");
  TGeoMedium *medFGS003    = mgr->GetMedium("ITS_FGS003$"); 
  TGeoMedium *medCarbonFleece = mgr->GetMedium("ITS_CarbonFleece$"); 

  // Local parameters
  Double_t kConeOutRadius =(0.1024+0.0025)/2;//0.107/2;
  Double_t kConeInRadius = 0.1024/2;//0.10105/2
  Double_t kStaveLength = zsta;
  Double_t kStaveWidth = xsta*2;
  Double_t kWidth = (kStaveWidth)/4;
  Double_t kStaveHeight = 0.283;//0.33;
  Double_t kHeight = (kStaveHeight)/2;
  Double_t kAlpha = 57;//56.31;
  Double_t kTheta = kAlpha*TMath::DegToRad();
  Double_t kS1 = ((kStaveWidth)/4)/TMath::Sin(kTheta);
  Double_t kL1 = (kStaveWidth/4)/TMath::Tan(kTheta);
  Double_t kS2 = sqrt(kHeight*kHeight + kS1*kS1);//TMath::Sin(kThe2);
  Double_t kThe2 = TMath::ATan(kHeight/(0.375-0.036));
  Double_t kBeta = kThe2*TMath::RadToDeg();
  Double_t klay1 = 0.003;//Amec carbon
  Double_t klay2 = 0.002;//C Fleece carbon
  Double_t klay3 = 0.007;//CFplate K13D2U carbon
  Double_t klay4 = 0.007;//GluekStaveLength/2
  Double_t klay5 = 0.01;//Flex cable
  Double_t kTopVertexMaxWidth = 0.072;
  Double_t kTopVertexHeight = 0.04;
  Double_t kSideVertexMWidth = 0.052;
  Double_t kSideVertexHeight = 0.11; 

 
  Int_t  loop = (Int_t)(kStaveLength/(2*kL1));

  char volname[30];
  snprintf(volname, 30, "%s%d_StaveStruct", AliITSUGeomTGeo::GetITSStavePattern(), fLayerNumber);

  Double_t z=0, y=-(2*kConeOutRadius)+klay1+klay2+fSensorThick/2-0.0004, x=0;

  TGeoVolume *mechStavVol = 0;

  if (fBuildLevel < 5) {
    // world (trapezoid)
    TGeoXtru *mechStruct = new TGeoXtru(2); //z sections
    Double_t xv[6] = {kStaveWidth/2,kStaveWidth/2,0.012,-0.012,-kStaveWidth/2,-kStaveWidth/2}; 
    /* Double_t yv[6] = {-2*(kConeOutRadius+klay1+1.5*klay2+klay3+klay4+fSensorThick+klay5),
                        0-0.02,kStaveHeight+0.01,kStaveHeight+0.01,0-0.02,
			-2*(kConeOutRadius+klay1+1.5*klay2+klay3+klay4+fSensorThick+klay5)};  // (kConeOutRadius*2)-0.0635 */
    Double_t yv[6] = {-(kConeOutRadius*2)-0.07295,0-0.02,kStaveHeight+0.01,kStaveHeight+0.01,0-0.02,-(kConeOutRadius*2)-0.07295};  // (kConeOutRadius*2)-0.064
    mechStruct->DefinePolygon(6,xv,yv);
    mechStruct->DefineSection(0,-kStaveLength,0,0,1.);
    mechStruct->DefineSection(1,kStaveLength,0,0,1.);

    mechStavVol = new TGeoVolume(volname, mechStruct, medAir);
    mechStavVol->SetLineColor(12);
    mechStavVol->SetFillColor(12); 
    mechStavVol->SetVisibility(kTRUE);  
      
    //Polyimide Pipe Kapton grey-35 
    TGeoCone *cone1 = new TGeoCone(kStaveLength,kConeInRadius,kConeOutRadius-0.0001,kConeInRadius,kConeOutRadius-0.0001);
    TGeoVolume *volCone1= new TGeoVolume("PolyimidePipe", cone1, medKapton);
    volCone1->SetFillColor(35);
    volCone1->SetLineColor(35);
    mechStavVol->AddNode(volCone1,1,new TGeoTranslation(x+0.25,y,z));
    mechStavVol->AddNode(volCone1,2,new TGeoTranslation(x-0.25,y,z));
    }

  if (fBuildLevel < 4) {
    TGeoTube *coolTubeW = new TGeoTube(0.,kConeInRadius-0.0001,kStaveLength);
    TGeoVolume *volCoolTubeW= new TGeoVolume("Water", coolTubeW, medWater);
    volCoolTubeW->SetFillColor(4);
    volCoolTubeW->SetLineColor(4);
    mechStavVol->AddNode(volCoolTubeW,0,new TGeoTranslation(x-0.25,y,z));
    mechStavVol->AddNode(volCoolTubeW,1,new TGeoTranslation(x+0.25,y,z));
  }

  if (fBuildLevel < 3) {
    //top fillament
    // Top filament M60J black-12 Carbon structure TGeoBBox (length,thickness,width)
    TGeoBBox *t2=new TGeoBBox(kS2-0.028,0.02/2,0.02/2); //0.04/2//TGeoBBox *t2=new TGeoBBox(kS2,0.01,0.02);//kS2-0.03 old Config.C
    TGeoVolume *volT2=new TGeoVolume("TopFilament", t2, medM60J3K);
    volT2->SetLineColor(12);
    volT2->SetFillColor(12); 
    for(int i=0;i<loop;i++){// i<28;i++){
      // 1) Front Left Top Filament
       mechStavVol->AddNode(volT2,i*4+1,new TGeoCombiTrans(x+kWidth+0.0036,y+kHeight+0.01,z-kStaveLength+0.1+(i*4*kL1)+kS1/2, new TGeoRotation("volT2",90,90-kAlpha,90-kBeta)));
      // 2) Front Right Top Filament
      mechStavVol->AddNode(volT2,i*4+2,new TGeoCombiTrans(x-kWidth-0.0036,y+kHeight+0.01,z-kStaveLength+0.1+(i*4*kL1)+kS1/2, new TGeoRotation("volT2",90,-90+kAlpha,-90+kBeta)));
      // 3) Back Left  Top Filament
      mechStavVol->AddNode(volT2,i*4+3,new TGeoCombiTrans(x+kWidth+0.0036,y+kHeight+0.01,z-kStaveLength+0.1+2*kL1+(i*4*kL1)+kS1/2, new TGeoRotation("volT2",90,-90+kAlpha,90-kBeta)));
      // 4) Back Right Top Filament
      mechStavVol->AddNode(volT2,i*4+4,new TGeoCombiTrans(x-kWidth-0.0036,y+kHeight+0.01,z-kStaveLength+0.1+2*kL1+(i*4*kL1)+kS1/2, new TGeoRotation("volT2",90,90-kAlpha,-90+kBeta)));
   }
 
     //Vertex  structure 

      //top ver trd1
      TGeoTrd1 *trd1 = new TGeoTrd1(0,kTopVertexMaxWidth/2,kStaveLength,kTopVertexHeight/2);
      TGeoVolume *ibdv = new TGeoVolume("TopVertex",trd1,medM60J3K);
      ibdv->SetFillColor(12);
      ibdv->SetLineColor(12);
      mechStavVol->AddNode(ibdv,1,new TGeoCombiTrans(x,y+kStaveHeight+0.03,z,new TGeoRotation("ibdv",0.,-90,0)));//y+kStaveHeight+0.056

      //left trd2
      TGeoTrd1 *trd2 = new TGeoTrd1(0,kSideVertexMWidth/2,kStaveLength, kSideVertexHeight/2);
      TGeoVolume *ibdv2 = new TGeoVolume("LeftVertex",trd2,medM60J3K);
      ibdv2->SetFillColor(12);
      ibdv2->SetLineColor(12);
      mechStavVol->AddNode(ibdv2,1,new TGeoCombiTrans(x+kStaveWidth/2-0.06,y-0.0348,z,new TGeoRotation("ibdv2",-103.3,90,0))); //x-kStaveWidth/2-0.09 old Config.C y-0.0355,

      //right trd3
      TGeoTrd1 *trd3 = new TGeoTrd1(0,kSideVertexMWidth/2,kStaveLength, kSideVertexHeight/2);
      TGeoVolume *ibdv3 = new TGeoVolume("RightVertex",trd3,medM60J3K);
      ibdv3->SetFillColor(12);
      ibdv3->SetLineColor(12);
      mechStavVol->AddNode(ibdv3,1,new TGeoCombiTrans(x-kStaveWidth/2+0.06,y-0.0348,z,new TGeoRotation("ibdv3",103.3,90,0))); //x-kStaveWidth/2+0.09 old Config.C
      
     //Carbon Fleece
      TGeoConeSeg *cons2 = new TGeoConeSeg(zsta,kConeOutRadius+klay1,kConeOutRadius+klay1+klay2,kConeOutRadius+klay1,kConeOutRadius+klay1+klay2,0,180); 
      TGeoVolume *cone12 = new TGeoVolume("CarbonFleecePipeCover",cons2,medCarbonFleece);
      cone12->SetFillColor(28);
      cone12->SetLineColor(28);
      mechStavVol->AddNode(cone12,1,new TGeoCombiTrans(x+0.25,y,z,new TGeoRotation("cone12",0,0,0)));
      mechStavVol->AddNode(cone12,2,new TGeoCombiTrans(x-0.25,y,z,new TGeoRotation("cone12",0,0,0)));

      TGeoBBox *box3 = new TGeoBBox((0.50-(2*(kConeOutRadius+klay1)))/2,klay2/2,zsta);//kStaveLength-0.50);
      TGeoVolume *plate3 = new TGeoVolume("CarbonFleeceMiddle",box3,medCarbonFleece);
      plate3->SetFillColor(28);
      plate3->SetLineColor(28);
      mechStavVol->AddNode(plate3,1,new TGeoCombiTrans(x,y-kConeOutRadius+klay1+(klay2/2),z,new TGeoRotation("plate3",0,0,0)));

      TGeoBBox *box31 = new TGeoBBox((0.75-0.25-kConeOutRadius-klay1)/2+0.0025,klay2/2,zsta);
      TGeoVolume *plate31 = new TGeoVolume("CarbonFleeceLeftRight",box31,medCarbonFleece);
      plate31->SetFillColor(28);
      plate31->SetLineColor(28);
      mechStavVol->AddNode(plate31,1,new TGeoCombiTrans(x+0.25+kConeOutRadius+klay1+(0.75-0.25-kConeOutRadius-klay1)/2,y-kConeOutRadius+klay1+(klay2/2),z,new TGeoRotation("plate31",0,0,0)));
      mechStavVol->AddNode(plate31,2,new TGeoCombiTrans(x-0.25-kConeOutRadius-klay1-(0.75-0.25-kConeOutRadius-klay1)/2,y-kConeOutRadius+klay1+(klay2/2),z,new TGeoRotation("plate31",0,0,0)));

      TGeoBBox *box32 = new TGeoBBox((klay2/2),(kConeOutRadius-klay1)/2,zsta);
      TGeoVolume *plate32 = new TGeoVolume("CarbonFleeceVertical",box32,medCarbonFleece);
      plate32->SetFillColor(28);
      plate32->SetLineColor(28);
      mechStavVol->AddNode(plate32,1,new TGeoCombiTrans(x+0.25+kConeOutRadius+klay1+(klay2/2),y+(klay1-kConeOutRadius)/2,z,new TGeoRotation("plate32",0,0,0)));
      mechStavVol->AddNode(plate32,2,new TGeoCombiTrans(x+0.25-kConeOutRadius-klay1-(klay2/2),y+(klay1-kConeOutRadius)/2,z,new TGeoRotation("plate32",0,0,0)));
      mechStavVol->AddNode(plate32,3,new TGeoCombiTrans(x-0.25+kConeOutRadius+klay1+(klay2/2),y+(klay1-kConeOutRadius)/2,z,new TGeoRotation("plate32",0,0,0)));
      mechStavVol->AddNode(plate32,4,new TGeoCombiTrans(x-0.25-kConeOutRadius-klay1-(klay2/2),y+(klay1-kConeOutRadius)/2,z,new TGeoRotation("plate32",0,0,0)));

     //Amec Thermasol red-2 cover tube FGS300 or Carbon Paper
      TGeoConeSeg *cons1 = new TGeoConeSeg(zsta,kConeOutRadius,kConeOutRadius+klay1-0.0001,kConeOutRadius,kConeOutRadius+klay1-0.0001,0,180);//kConeOutRadius+klay1-0.0001
      TGeoVolume *cone11 = new TGeoVolume("ThermasolPipeCover",cons1,medFGS003);
      cone11->SetFillColor(2);
      cone11->SetLineColor(2);
      mechStavVol->AddNode(cone11,1,new TGeoCombiTrans(x+0.25,y,z,new TGeoRotation("cone11",0,0,0)));
      mechStavVol->AddNode(cone11,2,new TGeoCombiTrans(x-0.25,y,z,new TGeoRotation("cone11",0,0,0)));

      TGeoBBox *box2 = new TGeoBBox((0.50-(2*kConeOutRadius))/2,(klay1/2),zsta);//kStaveLength-0.50);
      TGeoVolume *plate2 = new TGeoVolume("ThermasolMiddle",box2,medFGS003);
      plate2->SetFillColor(2);
      plate2->SetLineColor(2);
      mechStavVol->AddNode(plate2,1,new TGeoCombiTrans(x,y-kConeOutRadius+(klay1/2),z,new TGeoRotation("plate2",0,0,0)));

      TGeoBBox *box21 = new TGeoBBox((0.75-0.25-kConeOutRadius-klay1)/2+0.0025,(klay1/2),zsta);
      TGeoVolume *plate21 = new TGeoVolume("ThermasolLeftRight",box21,medFGS003);
      plate21->SetFillColor(2);
      plate21->SetLineColor(2);
      mechStavVol->AddNode(plate21,1,new TGeoCombiTrans(x+0.25+kConeOutRadius+(0.75-0.25-kConeOutRadius)/2-(klay1/2)+0.0025,y-kConeOutRadius+(klay1/2),z,new TGeoRotation("plate21",0,0,0)));
      mechStavVol->AddNode(plate21,2,new TGeoCombiTrans(x-0.25-kConeOutRadius-(0.75-0.25-kConeOutRadius)/2+(klay1/2)-0.0025,y-kConeOutRadius+(klay1/2),z,new TGeoRotation("plate21",0,0,0)));

      TGeoBBox *box22 = new TGeoBBox((klay1/2),kConeOutRadius/2,zsta);
      TGeoVolume *plate22 = new TGeoVolume("ThermasolVertical",box22,medFGS003);
      plate22->SetFillColor(2);
      plate22->SetLineColor(2);
      mechStavVol->AddNode(plate22,1,new TGeoCombiTrans(x+0.25+kConeOutRadius+(klay1/2),y-kConeOutRadius/2,z,new TGeoRotation("plate22",0,0,0)));
      mechStavVol->AddNode(plate22,2,new TGeoCombiTrans(x+0.25-kConeOutRadius-(klay1/2),y-kConeOutRadius/2,z,new TGeoRotation("plate22",0,0,0)));
      mechStavVol->AddNode(plate22,3,new TGeoCombiTrans(x-0.25+kConeOutRadius+(klay1/2),y-kConeOutRadius/2,z,new TGeoRotation("plate22",0,0,0)));
      mechStavVol->AddNode(plate22,4,new TGeoCombiTrans(x-0.25-kConeOutRadius-(klay1/2),y-kConeOutRadius/2,z,new TGeoRotation("plate22",0,0,0)));

     //K13D2U CF plate
      TGeoBBox *box1 = new TGeoBBox(2*kWidth,(klay3)/2,zsta);
      TGeoVolume *plate1 = new TGeoVolume("CFPlate",box1,medK13D2U2k);
      plate1->SetFillColor(5);
      plate1->SetLineColor(5);
      mechStavVol->AddNode(plate1,1,new TGeoCombiTrans(x,y-(kConeOutRadius+(klay3/2)),z,new TGeoRotation("plate1",0,0,0)));

     //C Fleece bottom plate 
      TGeoBBox *box6 = new TGeoBBox(2*kWidth,(klay2)/2,zsta);
      TGeoVolume *plate6 = new TGeoVolume("CarbonFleeceBottom",box6,medCarbonFleece);
      plate6->SetFillColor(2);
      plate6->SetLineColor(2);
      mechStavVol->AddNode(plate6,1,new TGeoCombiTrans(x,y-(kConeOutRadius+klay3+(klay2/2)),z,new TGeoRotation("plate6",0,0,0)));

    }
   if (fBuildLevel < 2) {
      //Glue klayers and kapton
     TGeoBBox *glue = new TGeoBBox(kStaveWidth/2, (klay4)/2, zsta);
      TGeoVolume *volGlue=new TGeoVolume("Glue", glue, medGlue);
      volGlue->SetLineColor(5);
      volGlue->SetFillColor(5); 
      // mechStavVol->AddNode(volGlue, 0, new TGeoCombiTrans(x,y-(kConeOutRadius+klay3+klay2+(klay4/2)), z, new TGeoRotation("",0, 0, 0)));
      mechStavVol->AddNode(volGlue, 0, new TGeoCombiTrans(x,y-(kConeOutRadius+klay3+klay2+(klay4)/2)+0.00005, z, new TGeoRotation("",0, 0, 0)));
    }

     if (fBuildLevel < 1) {
     //Flex Cable or Bus
      TGeoBBox *kapCable = new TGeoBBox(kStaveWidth/2, klay5/2, zsta);//klay5/2
      TGeoVolume *volCable=new TGeoVolume("FlexCable", kapCable, medFlexCable);
      volCable->SetLineColor(28);
      volCable->SetFillColor(28); 
      //      mechStavVol->AddNode(volCable, 0, new TGeoCombiTrans(x, y-(kConeOutRadius+klay3+klay2+klay4+fSensorThick+(klay5)/2)+0.0002, z, new TGeoRotation("",0, 0, 0)));
      mechStavVol->AddNode(volCable, 0, new TGeoCombiTrans(x, y-(kConeOutRadius+klay3+klay2+klay4+(klay5)/2), z, new TGeoRotation("",0, 0, 0)));
      }
    // Done, return the stave structe
    return mechStavVol;
}

// model3
//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateStaveModelInnerB3(const Double_t xsta,
						    const Double_t zsta,
						    const TGeoManager *mgr){
//
// Create the mechanical stave structure for Model 3 of TDR
//
// Input:
//         xsta : X length
//         zsta : Z length
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      28 May 2013  Chinorat Kobdaj
// Updated:                   Mario Sitta
// Updated:                   Wanchaloem Poonsawat 
//
  
  // Materials defined in AliITSUv2
  TGeoMedium *medAir    = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medWater  = mgr->GetMedium("ITS_WATER$");

  TGeoMedium *medM60J3K    = mgr->GetMedium("ITS_M60J3K$"); 
  TGeoMedium *medKapton    = mgr->GetMedium("ITS_KAPTON(POLYCH2)$");
  TGeoMedium *medGlue      = mgr->GetMedium("ITS_GLUE$");
  TGeoMedium *medFlexCable = mgr->GetMedium("ITS_FLEXCABLE$");
  //TGeoMedium *medK13D2U2k  = mgr->GetMedium("ITS_K13D2U2k$");
  //TGeoMedium *medFGS003    = mgr->GetMedium("ITS_FGS003$"); 
  //TGeoMedium *medCarbonFleece = mgr->GetMedium("ITS_CarbonFleece$"); 

  // Local parameters
    Double_t kConeOutRadius = 0.15/2;
    Double_t kStaveLength = zsta*2;
    Double_t kStaveWidth = xsta*2;
    Double_t w = kStaveWidth/4;//1/2 of W
    Double_t staveHeight = 0.3;
    Double_t h = staveHeight/2;
    Double_t alpha = 90-33.;//90-30;
    Double_t the1 = alpha*TMath::DegToRad();
    Double_t s1 = w/TMath::Sin(the1);
    Double_t l = w/TMath::Tan(the1);
    Double_t s2 = TMath::Sqrt(h*h + s1*s1);//TMath::Sin(the2);
    Double_t the2 = TMath::ATan(h/s1);
    Double_t beta = the2*TMath::RadToDeg();
    Double_t klay4 = 0.007; //Glue
    Double_t klay5 = 0.01; //Flexcable
    Int_t  loop = (Int_t)((kStaveLength/(2*l))/2);
    Double_t hh = 0.01;
       Double_t ang1 = 0*TMath::DegToRad();
       Double_t ang2 = 0*TMath::DegToRad();
       Double_t ang3 = 0*TMath::DegToRad();
       Int_t chips = 4;
       Double_t headWidth=0.25;
       Double_t smcLength=kStaveLength/chips-2*headWidth;//6.25;
       Double_t smcWidth=kStaveWidth;
       Double_t smcSide1Thick=0.03;
       Double_t vaporThick=0.032;
       Double_t liquidThick=0.028;
       Double_t smcSide2Thick=0.01;
       Double_t smcSide3Thick=0.0055;
       Double_t smcSide4Thick=0.0095;
       Double_t smcSide5Thick=0.0075;
       Double_t smcSpace=0.01;


    char volname[30];
    snprintf(volname, 30, "%s%d_StaveStruct", AliITSUGeomTGeo::GetITSStavePattern(), fLayerNumber);
    
    // detailed structure ++++++++++++++
    Double_t z=0, y=0-0.007, x=0;

    // Polimide micro channels numbers
    Double_t yMC = y-h+0.01;
    Int_t nb = (Int_t)(kStaveWidth/0.1)+1;
    Double_t xstaMC = (nb*0.1-0.08)/2;


    TGeoVolume *mechStavVol = 0;
    if (fBuildLevel < 5) {
      // world (trapezoid)
      TGeoXtru *mechStruct = new TGeoXtru(2); //z sections
      Double_t xv[5] = {kStaveWidth/2+0.1,kStaveWidth/2+0.1,0,-kStaveWidth/2-0.1,-kStaveWidth/2-0.1};
      Double_t yv[5] = {-kConeOutRadius*2-0.07,0,staveHeight,0,-kConeOutRadius*2-0.07};    
      mechStruct->DefinePolygon(5,xv,yv);
      mechStruct->DefineSection(0,-kStaveLength-0.1,0,0,1.);
      mechStruct->DefineSection(1,kStaveLength+0.1,0,0,1.);
      mechStavVol = new TGeoVolume(volname, mechStruct, medAir);
      mechStavVol->SetLineColor(12);
      mechStavVol->SetFillColor(12); 
      mechStavVol->SetVisibility(kTRUE);

       // Silicon micro channels numbers
      
      TGeoBBox *tM0a=new TGeoBBox(smcWidth/2, 0.003/2, headWidth/2);
      TGeoVolume *volTM0a=new TGeoVolume("microChanTop1", tM0a, medKapton);
      volTM0a->SetLineColor(35);
      volTM0a->SetFillColor(35); 

      for(Int_t  mo=1; mo<=chips; mo++) {
      mechStavVol->AddNode(volTM0a, 0, new TGeoCombiTrans(x,yMC+0.03, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth+smcLength/2+(headWidth/2), new TGeoRotation("",ang1, ang2, ang3)));
      mechStavVol->AddNode(volTM0a, 1, new TGeoCombiTrans(x,yMC+0.03, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth-smcLength/2-(headWidth/2), new TGeoRotation("",ang1, ang2, ang3)));
      }
      TGeoBBox *tM0c=new TGeoBBox(0.3/2, 0.003/2,smcLength/2);
      TGeoVolume *volTM0c=new TGeoVolume("microChanTop2", tM0c, medKapton);
      volTM0c->SetLineColor(35);
      volTM0c->SetFillColor(35); 
      for(Int_t  mo=1; mo<=chips; mo++) {
      mechStavVol->AddNode(volTM0c, 0, new TGeoCombiTrans(x+(smcWidth/2)-(0.3/2),yMC+0.03, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));
      mechStavVol->AddNode(volTM0c, 1, new TGeoCombiTrans(x-(smcWidth/2)+(0.3/2),yMC+0.03, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0c1=new TGeoBBox(0.2225/2, 0.003/2,smcLength/2);
      TGeoVolume *volTM0c1=new TGeoVolume("microChanBot1", tM0c1, medKapton);
      volTM0c1->SetLineColor(6);
      volTM0c1->SetFillColor(6); 
      for(Int_t  mo=1; mo<=chips; mo++) {
      mechStavVol->AddNode(volTM0c1, 0, new TGeoCombiTrans(x+smcWidth/2-(smcSide1Thick)-(vaporThick)-(smcSide2Thick)-(smcSide3Thick)-(0.2225/2),yMC+0.03-hh-(0.003), z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      mechStavVol->AddNode(volTM0c1, 1, new TGeoCombiTrans(x-smcWidth/2+(smcSide1Thick)+(liquidThick)+(smcSide2Thick)+(smcSide4Thick)+(0.2225/2),yMC+0.03-hh-(0.003), z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0c2=new TGeoBBox(0.072/2, 0.003/2,smcLength/2);
      TGeoVolume *volTM0c2=new TGeoVolume("microChanBot2", tM0c2, medKapton);
      volTM0c2->SetLineColor(35);
      volTM0c2->SetFillColor(35); 
      for(Int_t  mo=1; mo<=chips; mo++) {
      mechStavVol->AddNode(volTM0c2, 0, new TGeoCombiTrans(x+smcWidth/2-(0.072/2),yMC+0.03-(0.035+0.0015)-(0.003)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0c2r=new TGeoBBox(0.068/2, 0.003/2,smcLength/2);
      TGeoVolume *volTM0c2r=new TGeoVolume("microChanBot3", tM0c2r, medKapton);
      volTM0c2r->SetLineColor(35);
      volTM0c2r->SetFillColor(35); 
      for(Int_t  mo=1; mo<=chips; mo++) {      
      mechStavVol->AddNode(volTM0c2r, 0, new TGeoCombiTrans(x-smcWidth/2+(0.068/2),yMC+0.03-(0.035+0.0015)-(0.003)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0d=new TGeoBBox(smcSide1Thick/2, 0.035/2,smcLength/2);
      TGeoVolume *volTM0d=new TGeoVolume("microChanSide1", tM0d, medKapton);
      volTM0d->SetLineColor(12);
      volTM0d->SetFillColor(12); 
      for(Int_t  mo=1; mo<=chips; mo++) {
      mechStavVol->AddNode(volTM0d, 0, new TGeoCombiTrans(x+smcWidth/2-(smcSide1Thick/2),yMC+0.03-0.0015-(0.035)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      mechStavVol->AddNode(volTM0d, 1, new TGeoCombiTrans(x-smcWidth/2+(smcSide1Thick/2),yMC+0.03-0.0015-(0.035)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }

      TGeoBBox *tM0d1=new TGeoBBox(smcSide2Thick/2, 0.035/2,smcLength/2);
      TGeoVolume *volTM0d1=new TGeoVolume("microChanSide2", tM0d1, medKapton);
      volTM0d1->SetLineColor(12);
      volTM0d1->SetFillColor(12); 
      for(Int_t  mo=1; mo<=chips; mo++) {
      mechStavVol->AddNode(volTM0d1, 0, new TGeoCombiTrans(x+smcWidth/2-(smcSide1Thick)-(vaporThick)-(smcSide2Thick/2),yMC+0.03-(0.003+0.035)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      mechStavVol->AddNode(volTM0d1, 1, new TGeoCombiTrans(x-smcWidth/2+(smcSide1Thick)+(liquidThick)+(smcSide2Thick/2),yMC+0.03-(0.003+0.035)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0d2=new TGeoBBox(smcSide3Thick/2, (hh+0.003)/2, smcLength/2);
      TGeoVolume *volTM0d2=new TGeoVolume("microChanSide3", tM0d2, medKapton);
      volTM0d2->SetLineColor(12);
      volTM0d2->SetFillColor(12); 
      for(Int_t  mo=1; mo<=chips; mo++) {
      mechStavVol->AddNode(volTM0d2, 0, new TGeoCombiTrans(x+smcWidth/2-(smcSide1Thick)-(vaporThick)-(smcSide2Thick)-(smcSide3Thick/2),yMC+0.03-(0.003+hh+0.003)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0d2r=new TGeoBBox(smcSide4Thick/2, (hh+0.003)/2, smcLength/2);
      TGeoVolume *volTM0d2r=new TGeoVolume("microChanSide4", tM0d2r, medKapton);
      volTM0d2r->SetLineColor(12);
      volTM0d2r->SetFillColor(12); 
      for(Int_t  mo=1; mo<=chips; mo++) {
      mechStavVol->AddNode(volTM0d2r, 0, new TGeoCombiTrans(x-smcWidth/2+(smcSide1Thick)+(liquidThick)+(smcSide2Thick)+(smcSide4Thick/2),yMC+0.03-(0.003+hh+0.003)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0e=new TGeoBBox(smcSide5Thick/2, hh/2,smcLength/2);
      TGeoVolume *volTM0e=new TGeoVolume("microChanSide5", tM0e, medKapton);    
      volTM0e->SetLineColor(12);
      volTM0e->SetFillColor(12); 
      for(Int_t  mo=1; mo<=chips; mo++) {
      for (Int_t ie=0;ie<11;ie++) {
	mechStavVol->AddNode(volTM0e, 0, new TGeoCombiTrans(x-(ie*(smcSpace+smcSide5Thick))+smcWidth/2-(smcSide1Thick)-(vaporThick)-(smcSide2Thick)-(smcSide3Thick)-smcSpace-(smcSide5Thick/2),yMC+0.03-(0.003+hh)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
	mechStavVol->AddNode(volTM0e, 1, new TGeoCombiTrans(x+(ie*(smcSpace+smcSide5Thick))-smcWidth/2+(smcSide1Thick)+(liquidThick)+(smcSide2Thick)+(smcSide4Thick)+smcSpace+(smcSide5Thick/2),yMC+0.03-(0.003+hh)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
         }
      }
      
      TGeoBBox *tM0f=new TGeoBBox(0.02/2, hh/2, smcLength/2);
      TGeoVolume *volTM0f=new TGeoVolume("microChanTop3", tM0f, medKapton);
      //Double_t smcChannels=12;
      Double_t smcCloseWallvapor=smcWidth/2-smcSide1Thick-vaporThick-smcSide2Thick-smcSide3Thick-12*smcSpace-11*smcSide5Thick;
      Double_t smcCloseWallliquid=smcWidth/2-smcSide1Thick-liquidThick-smcSide2Thick-smcSide4Thick-12*smcSpace-11*smcSide5Thick;
      volTM0f->SetLineColor(12);
      volTM0f->SetFillColor(12); 
      for(Int_t  mo=1; mo<=chips; mo++) {
       mechStavVol->AddNode(volTM0f, 0, new TGeoCombiTrans(x+smcCloseWallvapor-(0.02)/2,yMC+0.03-(0.003+hh)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
       mechStavVol->AddNode(volTM0f, 1, new TGeoCombiTrans(x-smcCloseWallliquid+(0.02)/2,yMC+0.03-(0.003+hh)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      //Head(back) microchannel

      TGeoBBox *tM0hb=new TGeoBBox(smcWidth/2, 0.025/2, headWidth/2);
      TGeoVolume *volTM0hb=new TGeoVolume("microChanHeadBackBottom1", tM0hb, medKapton);
      volTM0hb->SetLineColor(4);
      volTM0hb->SetFillColor(4); 
      for(Int_t  mo=1; mo<=chips; mo++) {
      mechStavVol->AddNode(volTM0hb, 0, new TGeoCombiTrans(x,yMC+0.03-0.0145-(0.025/2), z+(mo-3)*kStaveLength/4+smcLength/2+headWidth+smcLength/2+(headWidth/2), new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      mechStavVol->AddNode(volTM0hb, 1, new TGeoCombiTrans(x,yMC+0.03-0.0145-(0.025)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth-smcLength/2-(headWidth/2), new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0h1=new TGeoBBox(smcWidth/2, 0.013/2, 0.05/2);
      TGeoVolume *volTM0h1=new TGeoVolume("microChanHeadBackBottom2", tM0h1, medKapton);
      volTM0h1->SetLineColor(5);
      volTM0h1->SetFillColor(5); 
      for(Int_t  mo=1; mo<=chips; mo++) {
      mechStavVol->AddNode(volTM0h1, 0, new TGeoCombiTrans(x,yMC+0.03-0.0015-(0.013/2), z+(mo-3)*kStaveLength/4+smcLength/2+headWidth-smcLength/2-headWidth+(0.05/2), new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0h2=new TGeoBBox(smcWidth/2, 0.003/2, 0.18/2);
      TGeoVolume *volTM0h2=new TGeoVolume("microChanHeadBackBottom7", tM0h2, medKapton);
      volTM0h2->SetLineColor(6);
      volTM0h2->SetFillColor(6);
      for(Int_t  mo=1; mo<=chips; mo++) {
      mechStavVol->AddNode(volTM0h2, 0, new TGeoCombiTrans(x,yMC+0.03-0.0015-0.01-(0.003/2), z+(mo-3)*kStaveLength/4+smcLength/2+headWidth-smcLength/2-0.02-(0.18/2), new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0h3=new TGeoBBox(smcWidth/2, 0.013/2, 0.02/2);
      TGeoVolume *volTM0h3=new TGeoVolume("microChanHeadBackBottom3", tM0h3, medKapton);
      volTM0h3->SetLineColor(5);
      volTM0h3->SetFillColor(5); 
      for(Int_t  mo=1; mo<=chips; mo++) {
      mechStavVol->AddNode(volTM0h3, 0, new TGeoCombiTrans(x,yMC+0.03-0.0015-(0.013/2), z+(mo-3)*kStaveLength/4+smcLength/2+headWidth-smcLength/2-(0.02/2), new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0b1=new TGeoBBox(smcWidth/2, 0.013/2, 0.03/2);
      TGeoVolume *volTM0b1=new TGeoVolume("microChanHeadBackBottom4", tM0b1, medKapton);
      volTM0b1->SetLineColor(5);
      volTM0b1->SetFillColor(5); 
      for(Int_t  mo=1; mo<=chips; mo++) {
      mechStavVol->AddNode(volTM0b1, 0, new TGeoCombiTrans(x,yMC+0.03-0.0015-(0.013/2), z+(mo-3)*kStaveLength/4+smcLength/2+headWidth+smcLength/2+headWidth-(0.03/2), new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0b2=new TGeoBBox(smcWidth/2, 0.003/2, 0.2/2);
      TGeoVolume *volTM0b2=new TGeoVolume("microChanHeadBackBottom5", tM0b2, medKapton);
      volTM0b2->SetLineColor(6);
      volTM0b2->SetFillColor(6); 
      for(Int_t  mo=1; mo<=chips; mo++) {
      mechStavVol->AddNode(volTM0b2, 0, new TGeoCombiTrans(x,yMC+0.03-0.0015-0.01-(0.003/2), z+(mo-3)*kStaveLength/4+smcLength/2+headWidth+smcLength/2+0.02+(0.2/2), new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0b3=new TGeoBBox(smcWidth/2, 0.013/2, 0.02/2);
      TGeoVolume *volTM0b3=new TGeoVolume("microChanHeadBackBottom6", tM0b3, medKapton);
      volTM0b3->SetLineColor(5);
      volTM0b3->SetFillColor(5); 
      for(Int_t  mo=1; mo<=chips; mo++) {
      mechStavVol->AddNode(volTM0b3, 0, new TGeoCombiTrans(x,yMC+0.03-0.0015-(0.013/2), z+(mo-3)*kStaveLength/4+smcLength/2+headWidth+smcLength/2+(0.02/2), new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
     
      TGeoBBox *tM0b=new TGeoBBox(0.02/2, 0.02/2, zsta);
      TGeoVolume *volTM0b=new TGeoVolume("microChanWalls", tM0b, medKapton);
      volTM0b->SetLineColor(35);
      volTM0b->SetFillColor(35); 
      for (Int_t ib=0;ib<nb;ib++) {
	//mechStavVol->AddNode(volTM0b, ib, new TGeoCombiTrans(x+ib*0.1-xstaMC+0.01,yMC, z, new TGeoRotation("",0, 0, 0)));
      }
      
      } 
    
    if (fBuildLevel < 4) {

      //**********cooling  inlet outlet

      TGeoBBox *tM0dv=new TGeoBBox(vaporThick/2, 0.035/2,smcLength/2);
      TGeoVolume *volTM0dv=new TGeoVolume("microChanVapor", tM0dv, medWater);
      volTM0dv->SetLineColor(2);
      volTM0dv->SetFillColor(2);
      for(Int_t  mo=1; mo<=chips; mo++) {
      mechStavVol->AddNode(volTM0dv, 0, new TGeoCombiTrans(x+smcWidth/2-(smcSide1Thick)-(vaporThick/2),yMC+0.03-0.0015-(0.035)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0dl=new TGeoBBox(liquidThick/2, 0.035/2,smcLength/2);
      TGeoVolume *volTM0dl=new TGeoVolume("microChanLiquid", tM0dl, medWater);
      volTM0dl->SetLineColor(3);
      volTM0dl->SetFillColor(3); 
      for(Int_t  mo=1; mo<=chips; mo++) {
      mechStavVol->AddNode(volTM0dl, 0, new TGeoCombiTrans(x-smcWidth/2+(smcSide1Thick)+(liquidThick/2),yMC+0.03-0.0015-(0.035)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      // small cooling fluid now using water wait for freeon value  
      TGeoBBox *tM0dlq=new TGeoBBox(smcSpace/2, hh/2,smcLength/2);
      TGeoVolume *volTM0dlq=new TGeoVolume("smallLiquid", tM0dlq, medWater);
      volTM0dlq->SetLineColor(3);
      volTM0dlq->SetFillColor(3); 
      TGeoBBox *tM0dvp=new TGeoBBox(smcSpace/2, hh/2,smcLength/2);
      TGeoVolume *volTM0dvp=new TGeoVolume("microChanVapor", tM0dvp, medWater);
      volTM0dvp->SetLineColor(2);
      volTM0dvp->SetFillColor(2); 
      for(Int_t  mo=1; mo<=chips; mo++) {
      for (Int_t is=0;is<12;is++) {
	mechStavVol->AddNode(volTM0dlq, 0, new TGeoCombiTrans(x+(is*(smcSpace+smcSide5Thick))-smcWidth/2+(smcSide1Thick)+(vaporThick)+(smcSide2Thick)+(smcSide3Thick)+smcSpace/2,yMC+0.03-(0.003+hh)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
	mechStavVol->AddNode(volTM0dvp, 1, new TGeoCombiTrans(x-(is*(smcSpace+smcSide5Thick))+smcWidth/2-(smcSide1Thick)-(vaporThick)-(smcSide2Thick)-(smcSide3Thick)-smcSpace/2,yMC+0.03-(0.003+hh)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      }

      //*************

    }
    
    if (fBuildLevel < 3) {

      //Bottom filament CFRP black-12 Carbon structure TGeoBBox (thickness,width,length)
 
      Double_t filWidth = 0.04;
      Double_t filHeight= 0.02;
      TGeoBBox *t1=new TGeoBBox(filHeight/2,filWidth/2,s1);
      TGeoVolume *volT1=new TGeoVolume("bottomFilament", t1, medM60J3K);
      volT1->SetLineColor(12);
      volT1->SetFillColor(12); 
      for(int i=0;i<loop;i++){//i<30;i++){
       	mechStavVol->AddNode(volT1,4*i+0,
				    new TGeoCombiTrans(x+w,y-h+0.04+filHeight/2,z-kStaveLength/2+(4*l*i)+s1/2, 
						       new TGeoRotation("volT1",-90,alpha,0)));
	mechStavVol->AddNode(volT1,4*i+1,
				    new TGeoCombiTrans(x-w,y-h+0.04+filHeight/2,z-kStaveLength/2+(4*l*i)+s1/2, 
						       new TGeoRotation("volT1",90,alpha,0)));
	mechStavVol->AddNode(volT1,4*i+2,
				    new TGeoCombiTrans(x+w,y-h+0.04+filHeight/2,z-kStaveLength/2+2*l+(i*4*l)+s1/2, 
						       new TGeoRotation("volT1",-90,-alpha,0)));
	mechStavVol->AddNode(volT1,4*i+3,
				    new TGeoCombiTrans(x-w,y-h+0.04+filHeight/2,z-kStaveLength/2+2*l+(i*4*l)+s1/2, 
						       new TGeoRotation("volT1",-90,+alpha,0)));
	}
 
     // Top filament CERP black-12 Carbon structure TGeoBBox (length,thickness,width)

      TGeoBBox *t2=new TGeoBBox(s2,filHeight/2,filWidth/2);
      TGeoVolume *volT2=new TGeoVolume("topFilament", t2, medM60J3K);
      volT2->SetLineColor(12);
      volT2->SetFillColor(12); 
      for(int i=0;i<loop;i++){ //i<30;i++){
       	mechStavVol->AddNode(volT2,4*i+0,
				    new TGeoCombiTrans(x+w,y+0.04+filHeight/2,z-kStaveLength/2+(i*4*l)+s1/2,
						       new TGeoRotation("volT2",90,90-alpha,90-beta)));
	mechStavVol->AddNode(volT2,4*i+1,
				    new TGeoCombiTrans(x-w,y+0.04+filHeight/2,z-kStaveLength/2+(i*4*l)+s1/2,
						       new TGeoRotation("volT2",90,-90+alpha,-90+beta)));
	mechStavVol->AddNode(volT2,4*i+2,
				    new TGeoCombiTrans(x+w,y+0.04+filHeight/2,z-kStaveLength/2+2*l+(i*4*l)+s1/2,
						       new TGeoRotation("volT2",90,-90+alpha,90-beta)));
	mechStavVol->AddNode(volT2,4*i+3,
				    new TGeoCombiTrans(x-w,y+0.04+filHeight/2,z-kStaveLength/2+2*l+(i*4*l)+s1/2, 
						       new TGeoRotation("volT2",90,90-alpha,-90+beta)));
	}
    }

    if (fBuildLevel < 2) {

      // Glue Filament and Silicon MicroChannel
      TGeoBBox *tM0=new TGeoBBox(xstaMC/5, klay4/2, zsta);
      TGeoVolume *volTM0=new TGeoVolume("glueFM", tM0,medGlue );
      volTM0->SetLineColor(5);
      volTM0->SetFillColor(5); 
      mechStavVol->AddNode(volTM0, 0, new TGeoCombiTrans(x-xsta/2-0.25,0.03+yMC, z, new TGeoRotation("",0, 0, 0)));
      mechStavVol->AddNode(volTM0, 1, new TGeoCombiTrans(x+xsta/2+0.25,0.03+yMC, z, new TGeoRotation("",0, 0, 0)));

            
      // Glue microchannel and sensor
      TGeoBBox *glueM = new TGeoBBox(xstaMC/5, klay4/2, zsta);
      TGeoVolume *volGlueM=new TGeoVolume("glueMS", glueM, medGlue);
      volGlueM->SetLineColor(5);
      volGlueM->SetFillColor(5); 
      mechStavVol->AddNode(volGlueM, 0, new TGeoCombiTrans(x-xsta/2-0.25,yMC-0.01, z, new TGeoRotation("",0, 0, 0)));
      mechStavVol->AddNode(volGlueM, 1, new TGeoCombiTrans(x+xsta/2+0.25,yMC-0.01, z, new TGeoRotation("",0, 0, 0)));
     
       // Glue sensor and kapton
      TGeoBBox *glue = new TGeoBBox(xsta, klay4/2, zsta);
      TGeoVolume *volGlue=new TGeoVolume("glueSensorBus", glue, medGlue);
      volGlue->SetLineColor(5);
      volGlue->SetFillColor(5); 
       mechStavVol->AddNode(volGlue, 1, new TGeoCombiTrans(x, y-0.154-fSensorThick-klay4/2, z, new TGeoRotation("",0, 0, 0)));
    }      

    if (fBuildLevel < 1) {
      TGeoBBox *kapCable = new TGeoBBox(xsta, klay5/2, zsta);
      TGeoVolume *volCable=new TGeoVolume("Flexcable", kapCable, medFlexCable);
      volCable->SetLineColor(28);
      volCable->SetFillColor(28); 
      mechStavVol->AddNode(volCable, 0, new TGeoCombiTrans(x, y-0.154-fSensorThick-klay4-klay5/2, z, new TGeoRotation("",0, 0, 0)));
    }

  // Done, return the stave structur
    return mechStavVol;
 }

// model4
//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateStaveModelInnerB4(const Double_t xstave,
						    const Double_t zstave,
						    const TGeoManager *mgr){
//
// Create the mechanical stave structure for Model 2.2 of TDR
// Logic is similar to method CreateStaveModelInnerB22
// but completely rewritten:
// - code completely revised, made systematic and more similar to OB
// - fix some inconsistencies (stave element sequence, empty space)
// - use static const as parameters
// - comply with latest (nov '14) C.Gargiulo data
//
// Input:
//         xstave : stave X half length
//         zstave : stave Z half length
//         mgr    : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      04 Dec 2014  Mario Sitta
// Updated:      03 Mar 2015  Mario Sitta  FPC in right position (beyond chip)
// Updated:      06 Mar 2015  Mario Sitta  Space Frame corrected (C.G. data)
// Updated:      30 Apr 2015  Mario Sitta  End-stave connectors added
//

  
  // Local parameters
  Double_t layerHeight = 0.;

  Double_t rPipeMin = fgkIBCoolPipeInnerD/2;
  Double_t rPipeMax = rPipeMin + fgkIBCoolPipeThick;

  Double_t topFilTheta = fgkIBTopFilamentAlpha*TMath::DegToRad();
  Double_t topFilLProj = xstave/TMath::Sin(topFilTheta); // Top filament length projected on stave XZ plane
  Double_t topFilYLen = xstave/TMath::Tan(topFilTheta); // Filament length on Y
  Int_t  nFilaments = (Int_t)(zstave/topFilYLen);
  // Question: would it be better to fix the number of filaments and
  // compute the angle alpha from it, or leave as it is now, i.e. fix the
  // filament inclination angle alpha and compute their number ?

  const Int_t nv = 6;
  Double_t xv[nv], yv[nv]; // The stave container Xtru
  Double_t xlen, ylen, zlen;
  Double_t xpos, ypos, zpos, ylay;
  Double_t beta, gamma, theta;


  // First create all needed shapes
  TGeoBBox *glue     = new TGeoBBox(xstave, fgkIBGlueThick/2, zstave);

  TGeoBBox *fleecbot = new TGeoBBox(xstave, fgkIBCarbonFleeceThick/2, zstave);

  TGeoBBox *cfplate  = new TGeoBBox(xstave, fgkIBK13D2UThick/2, zstave);

  TGeoTube *pipe     = new TGeoTube(rPipeMin, rPipeMax, zstave);

  TGeoTube *water    = new TGeoTube(0., rPipeMin, zstave);

  TGeoTubeSeg *cpaptub  = new TGeoTubeSeg(rPipeMax,
					  rPipeMax + fgkIBCarbonPaperThick,
					  zstave, 0, 180);

  TGeoBBox *cpapvert = new TGeoBBox(fgkIBCarbonPaperThick/2,
				    pipe->GetRmax()/2, zstave);

  xlen = fgkIBCoolPipeXDist/2 - pipe->GetRmax() - fgkIBCarbonPaperThick;
  TGeoBBox *cpapmid  = new TGeoBBox(xlen, fgkIBCarbonPaperThick/2, zstave);

  xlen = xstave -fgkIBCoolPipeXDist/2 -pipe->GetRmax() -fgkIBCarbonPaperThick;
  TGeoBBox *cpaplr   = new TGeoBBox(xlen/2, fgkIBCarbonPaperThick/2, zstave);

  TGeoTubeSeg *fleecpipe = new TGeoTubeSeg(cpaptub->GetRmax(),
			       cpaptub->GetRmax() + fgkIBCarbonFleeceThick,
					   zstave, 0, 180); 

  TGeoBBox *fleecvert = new TGeoBBox(fgkIBCarbonFleeceThick/2,
			 	     (pipe->GetRmax()-fgkIBCarbonPaperThick)/2,
				     zstave);

  xlen = fgkIBCoolPipeXDist/2 - pipe->GetRmax() - fgkIBCarbonPaperThick
       - fgkIBCarbonFleeceThick;
  TGeoBBox *fleecmid  = new TGeoBBox(xlen, fgkIBCarbonFleeceThick/2, zstave);

  xlen = xstave - fgkIBCoolPipeXDist/2 - pipe->GetRmax()
       - fgkIBCarbonPaperThick - fgkIBCarbonFleeceThick;
  TGeoBBox *fleeclr   = new TGeoBBox(xlen/2, fgkIBCarbonFleeceThick/2, zstave);

  // The spaceframe structure
  TGeoTrd1 *topv  = new TGeoTrd1(fgkIBTopVertexWidth1/2,
				 fgkIBTopVertexWidth2/2, zstave,
				 fgkIBTopVertexHeight/2);

  xv[0] = 0;
  yv[0] = 0;
  xv[1] = fgkIBSideVertexWidth;
  yv[1] = yv[0];
  xv[2] = xv[0];
  yv[2] = fgkIBSideVertexHeight;

  TGeoXtru *sidev = new TGeoXtru(2);
  sidev->DefinePolygon(3, xv, yv);
  sidev->DefineSection(0,-zstave);
  sidev->DefineSection(1, zstave);

  TGeoBBox *topfil = new TGeoBBox(fgkIBTopFilamentLength/2,
				  fgkIBTopFilamentSide/2,
				  fgkIBTopFilamentSide/2);


  // The half stave container (an XTru to avoid overlaps between neighbours)
  layerHeight = 2*(    glue->GetDY() + fleecbot->GetDY() + cfplate->GetDY()
                   + cpaplr->GetDY() +  fleeclr->GetDY() );

  xv[0] = xstave;
  yv[0] = 0;
  xv[1] = xv[0];
  yv[1] = layerHeight + fgkIBSideVertexHeight + topfil->GetDZ();;
  xv[2] = fgkIBTopVertexWidth2/2;
  yv[2] = fgkIBStaveHeight;
  for (Int_t i = 0; i<nv/2; i++) {
    xv[3+i] = -xv[2-i];
    yv[3+i] =  yv[2-i];
  }

  TGeoXtru *mechStruct = new TGeoXtru(2);
  mechStruct->DefinePolygon(nv, xv, yv);
  mechStruct->SetName("mechStruct");
  mechStruct->DefineSection(0,-zstave);
  mechStruct->DefineSection(1, zstave);

  // The connectors' containers
  zlen = fgkIBConnectBlockZLen - fgkIBConnTailZLen + fgkIBConnectAFitZOut;
  TGeoBBox *connAside = new TGeoBBox("connAsideIB",fgkIBConnectorXWidth/2,
				     fgkIBConnectorYTot/2, zlen/2);

  zlen = fgkIBConnectBlockZLen - fgkIBConnTailZLen;
  TGeoBBox *connCside = new TGeoBBox("connCsideIB",fgkIBConnectorXWidth/2,
				     fgkIBConnectorYTot/2, zlen/2);

  // The StaveStruct container, a Composite Shape
  ypos = connAside->GetDY() - fgkIBConnTailYShift + layerHeight;
  zpos = zstave + connAside->GetDZ();
  TGeoTranslation *transAside = new TGeoTranslation("transAsideIB",
						    0, ypos, zpos);
  transAside->RegisterYourself();

  ypos = connCside->GetDY() - fgkIBConnTailYShift + layerHeight;
  zpos = zstave + connCside->GetDZ();
  TGeoTranslation *transCside = new TGeoTranslation("transCsideIB",
						    0, ypos,-zpos);
  transCside->RegisterYourself();

  TGeoCompositeShape *mechStavSh = new TGeoCompositeShape(
	  "mechStruct+connAsideIB:transAsideIB+connCsideIB:transCsideIB");


  // We have all shapes: now create the real volumes

  TGeoMedium *medAir          = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medWater        = mgr->GetMedium("ITS_WATER$");
  TGeoMedium *medM55J6K       = mgr->GetMedium("ITS_M55J6K$"); 
  TGeoMedium *medM60J3K       = mgr->GetMedium("ITS_M60J3K$"); 
  TGeoMedium *medKapton       = mgr->GetMedium("ITS_KAPTON(POLYCH2)$");
  TGeoMedium *medGlue         = mgr->GetMedium("ITS_GLUE$");
  TGeoMedium *medK13D2U2k     = mgr->GetMedium("ITS_K13D2U2k$");
  TGeoMedium *medFGS003       = mgr->GetMedium("ITS_FGS003$"); 
  TGeoMedium *medCarbonFleece = mgr->GetMedium("ITS_CarbonFleece$"); 


  char volname[30];
  snprintf(volname, 30, "%s%d_StaveStruct",
	   AliITSUGeomTGeo::GetITSStavePattern(), fLayerNumber);
  TGeoVolume *mechStavVol = new TGeoVolume(volname, mechStavSh, medAir);
  mechStavVol->SetLineColor(12);
  mechStavVol->SetFillColor(12); 
  mechStavVol->SetVisibility(kFALSE);

  TGeoVolume *glueVol = new TGeoVolume("Glue", glue, medGlue);
  glueVol->SetLineColor(kBlack);
  glueVol->SetFillColor(kBlack);

  TGeoVolume *fleecbotVol = new TGeoVolume("CarbonFleeceBottom",
					   fleecbot, medCarbonFleece);
  fleecbotVol->SetFillColor(kViolet);
  fleecbotVol->SetLineColor(kViolet);

  TGeoVolume *cfplateVol = new TGeoVolume("CFPlate", cfplate, medK13D2U2k);
  cfplateVol->SetFillColor(5);  // Yellow
  cfplateVol->SetLineColor(5);

  TGeoVolume *pipeVol = new TGeoVolume("PolyimidePipe", pipe, medKapton);
  pipeVol->SetFillColor(35);  // Blue shade
  pipeVol->SetLineColor(35);

  TGeoVolume *waterVol= new TGeoVolume("Water", water, medWater);
  waterVol->SetFillColor(4);  // Bright blue
  waterVol->SetLineColor(4);

  TGeoVolume *cpaptubVol = new TGeoVolume("ThermasolPipeCover",
					  cpaptub, medFGS003);
  cpaptubVol->SetFillColor(2);  // Red
  cpaptubVol->SetLineColor(2);

  TGeoVolume *cpapvertVol = new TGeoVolume("ThermasolVertical",
					   cpapvert, medFGS003);
  cpapvertVol->SetFillColor(2);  // Red
  cpapvertVol->SetLineColor(2);

  TGeoVolume *cpapmidVol = new TGeoVolume("ThermasolMiddle",
					  cpapmid, medFGS003);
  cpapmidVol->SetFillColor(2);  // Red
  cpapmidVol->SetLineColor(2);

  TGeoVolume *cpaplrVol = new TGeoVolume("ThermasolLeftRight",
					 cpaplr, medFGS003);
  cpaplrVol->SetFillColor(2);  // Red
  cpaplrVol->SetLineColor(2);

  TGeoVolume *fleecpipeVol = new TGeoVolume("CarbonFleecePipeCover",
					    fleecpipe, medCarbonFleece);
  fleecpipeVol->SetFillColor(28);  // Brown shade
  fleecpipeVol->SetLineColor(28);

  TGeoVolume *fleecvertVol = new TGeoVolume("CarbonFleeceVertical",
					    fleecvert, medCarbonFleece);
  fleecvertVol->SetFillColor(28);  // Brown shade
  fleecvertVol->SetLineColor(28);

  TGeoVolume *fleecmidVol = new TGeoVolume("CarbonFleeceMiddle",
					   fleecmid, medCarbonFleece);
  fleecmidVol->SetFillColor(28);  // Brown shade
  fleecmidVol->SetLineColor(28);

  TGeoVolume *fleeclrVol = new TGeoVolume("CarbonFleeceLeftRight",
					  fleeclr, medCarbonFleece);
  fleeclrVol->SetFillColor(28);  // Brown shade
  fleeclrVol->SetLineColor(28);

  TGeoVolume *topvVol = new TGeoVolume("TopVertex", topv, medM55J6K);
  topvVol->SetFillColor(12);  // Gray shade
  topvVol->SetLineColor(12);
  
  TGeoVolume *sidevVol = new TGeoVolume("SideVertex", sidev, medM55J6K);
  sidevVol->SetFillColor(12);  // Gray shade
  sidevVol->SetLineColor(12);
  
  TGeoVolume *topfilVol = new TGeoVolume("TopFilament", topfil, medM60J3K);
  topfilVol->SetFillColor(12);  // Gray shade
  topfilVol->SetLineColor(12);
  

  // Now build up the half stave
  ypos = glue->GetDY();
  if (fBuildLevel < 2)   // Glue
    mechStavVol->AddNode(glueVol, 1, new TGeoTranslation(0, ypos, 0));

  ypos += (glue->GetDY() + fleecbot->GetDY());
  if (fBuildLevel < 5)   // Carbon
    mechStavVol->AddNode(fleecbotVol, 1, new TGeoTranslation(0, ypos, 0));

  ypos += (fleecbot->GetDY() + cfplate->GetDY());
  if (fBuildLevel < 5)   // Carbon
    mechStavVol->AddNode(cfplateVol, 1, new TGeoTranslation(0, ypos, 0));

  ylay = ypos + cfplate->GetDY(); // The level where tubes etc. lay

  xpos = fgkIBCoolPipeXDist/2;
  ypos = ylay + pipe->GetRmax();
  if (fBuildLevel < 4) { // Kapton
    mechStavVol->AddNode(pipeVol, 1, new TGeoTranslation(-xpos, ypos, 0));
    mechStavVol->AddNode(pipeVol, 2, new TGeoTranslation( xpos, ypos, 0));
  }

  if (fBuildLevel < 3) { // Water
    mechStavVol->AddNode(waterVol, 1, new TGeoTranslation(-xpos, ypos, 0));
    mechStavVol->AddNode(waterVol, 2, new TGeoTranslation( xpos, ypos, 0));
  }

  if (fBuildLevel < 5) { // Carbon (stave components)
    mechStavVol->AddNode(cpaptubVol, 1, new TGeoTranslation(-xpos, ypos, 0));
    mechStavVol->AddNode(cpaptubVol, 2, new TGeoTranslation( xpos, ypos, 0));

    mechStavVol->AddNode(fleecpipeVol,1, new TGeoTranslation(-xpos, ypos, 0));
    mechStavVol->AddNode(fleecpipeVol,2, new TGeoTranslation( xpos, ypos, 0));

    xpos = fgkIBCoolPipeXDist/2 - pipe->GetRmax() - cpapvert->GetDX();
    ypos = ylay + cpapvert->GetDY();
    mechStavVol->AddNode(cpapvertVol, 1, new TGeoTranslation(-xpos, ypos, 0));
    mechStavVol->AddNode(cpapvertVol, 2, new TGeoTranslation( xpos, ypos, 0));

    xpos = fgkIBCoolPipeXDist/2 + pipe->GetRmax() + cpapvert->GetDX();
    mechStavVol->AddNode(cpapvertVol, 3, new TGeoTranslation(-xpos, ypos, 0));
    mechStavVol->AddNode(cpapvertVol, 4, new TGeoTranslation( xpos, ypos, 0));

    ypos = ylay + fgkIBCarbonPaperThick/2;
    mechStavVol->AddNode(cpapmidVol, 1, new TGeoTranslation(0, ypos, 0));

    xpos = xstave - cpaplr->GetDX();
    mechStavVol->AddNode(cpaplrVol, 1, new TGeoTranslation(-xpos, ypos, 0));
    mechStavVol->AddNode(cpaplrVol, 2, new TGeoTranslation( xpos, ypos, 0));

    xpos = fgkIBCoolPipeXDist/2 - pipe->GetRmax() - 2*cpapvert->GetDX()
         - fleecvert->GetDX();
    ypos = ylay + fgkIBCarbonPaperThick + fleecvert->GetDY();
    mechStavVol->AddNode(fleecvertVol, 1, new TGeoTranslation(-xpos, ypos, 0));
    mechStavVol->AddNode(fleecvertVol, 2, new TGeoTranslation( xpos, ypos, 0));

    xpos = fgkIBCoolPipeXDist/2 + pipe->GetRmax() + 2*cpapvert->GetDX()
         + fleecvert->GetDX();
    mechStavVol->AddNode(fleecvertVol, 3, new TGeoTranslation(-xpos, ypos, 0));
    mechStavVol->AddNode(fleecvertVol, 4, new TGeoTranslation( xpos, ypos, 0));

    ypos = ylay + fgkIBCarbonPaperThick + fgkIBCarbonFleeceThick/2;
    mechStavVol->AddNode(fleecmidVol, 1, new TGeoTranslation(0, ypos, 0));

    xpos = xstave - fleeclr->GetDX();
    mechStavVol->AddNode(fleeclrVol, 1, new TGeoTranslation(-xpos, ypos, 0));
    mechStavVol->AddNode(fleeclrVol, 2, new TGeoTranslation( xpos, ypos, 0));
  }

  ylay += (fgkIBCarbonPaperThick + fgkIBCarbonFleeceThick);

  if (fBuildLevel < 5) { // Carbon (spaceframe)
    ypos = fgkIBStaveHeight - fgkIBTopFilamentSide - topv->GetDz(); // Due to rotation, z is on Y
    mechStavVol->AddNode(topvVol, 1,
			 new TGeoCombiTrans(0, ypos, 0,
					    new TGeoRotation("",0,-90,0)));

    xpos = xstave - sidev->GetX(1);
    ypos = ylay;
    mechStavVol->AddNode(sidevVol, 1, new TGeoTranslation( xpos, ypos, 0));
    mechStavVol->AddNode(sidevVol, 2, new TGeoCombiTrans(-xpos, ypos, 0,
					    new TGeoRotation("",90,180,-90)));

    gamma = fgkIBTopFilamentGamma;
    theta = 90. - fgkIBTopFilamentAlpha;
    xpos = xstave/2 + topfil->GetDZ();
    ypos = ( layerHeight + fgkIBStaveHeight )/2 +
	   fgkIBSideVertexWidth/TMath::Sin(gamma*TMath::DegToRad())/2 ;
    for(int i=0; i<nFilaments; i++){ // i<28 (?)
      // 1) Front Left Top Filament
//      zpos = -zstave + (i*2*topFilYLen) + topFilLProj/4; // ?????
//      zpos = -zstave + (i*2*topFilYLen) + topFilLProj/2;
      zpos = -zstave + (i*2*topFilYLen) + topFilLProj/4 + topfil->GetDY();
      mechStavVol->AddNode(topfilVol, i*4+1,
			 new TGeoCombiTrans( xpos, ypos, zpos,
			      new TGeoRotation("", 90, theta, gamma)));
      // 2) Front Right Top Filament
      mechStavVol->AddNode(topfilVol, i*4+2,
			 new TGeoCombiTrans(-xpos, ypos, zpos,
			      new TGeoRotation("", 90,-theta,-gamma)));
      // 3) Back Left  Top Filament
      zpos += topFilYLen;
      mechStavVol->AddNode(topfilVol, i*4+3,
			 new TGeoCombiTrans( xpos, ypos, zpos,
			      new TGeoRotation("", 90,-theta, gamma)));
      // 4) Back Right Top Filament
      mechStavVol->AddNode(topfilVol, i*4+4,
			 new TGeoCombiTrans(-xpos, ypos, zpos,
			      new TGeoRotation("", 90, theta,-gamma)));
    }
  }


  // Add the end-stave connectors
  TGeoVolume *connectorASide, *connectorCSide;

  // Check whether we have already all pieces
  // Otherwise create them
  connectorASide = mgr->GetVolume("IBConnectorASide");

  if (!connectorASide) {
    CreateIBConnectors(mgr);
    connectorASide = mgr->GetVolume("IBConnectorASide");
  }
  connectorCSide = mgr->GetVolume("IBConnectorCSide");

  ypos = ((TGeoBBox*)connectorASide->GetShape())->GetDY()
       - fgkIBConnTailYShift + ylay;
  zpos = zstave +
        (fgkIBConnectBlockZLen - fgkIBConnTailZLen + fgkIBConnectAFitZOut)/2;
  mechStavVol->AddNode(connectorASide, 1, new TGeoTranslation(0, ypos, zpos));

  zpos = zstave + (fgkIBConnectBlockZLen - fgkIBConnTailZLen)/2;
  mechStavVol->AddNode(connectorCSide, 1, new TGeoCombiTrans(0, ypos,-zpos,
					     new TGeoRotation("",90,180,-90)));


  // Done, return the stave structure
  return mechStavVol;
}

//________________________________________________________________________
void AliITSUv2Layer::CreateIBConnectors(const TGeoManager *mgr){
//
// Create the end-stave connectors for IB staves
// (simply call the actual creator methods)
//
// Input:
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      20 Apr 2015  Mario Sitta
//

  CreateIBConnectorsASide(mgr);
  CreateIBConnectorsCSide(mgr);

}

//________________________________________________________________________
void AliITSUv2Layer::CreateIBConnectorsASide(const TGeoManager *mgr){
//
// Create the A-Side end-stave connectors for IB staves
//
// Input:
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      22 Apr 2015  Mario Sitta
//

  // Local variables
  const Int_t nv = 8;
  Double_t xv[nv], yv[nv];
  Double_t xlen, ylen, zlen;
  Double_t xpos, ypos, zpos;
  char volname[30];


  // Gather all material pointers
  TGeoMedium *medAir      = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medPEEK     = mgr->GetMedium("ITS_PEEKCF30$");
  TGeoMedium *medWC       = mgr->GetMedium("ITS_TUNGCARB$");
  TGeoMedium *medInox304  = mgr->GetMedium("ITS_INOX304$");


  // First create all elements

  // The connector block, two Composite Shapes:
  // the body...
  xlen = fgkIBConnectorXWidth;
  ylen = fgkIBConnBodyYHeight;
  zlen = fgkIBConnectBlockZLen - fgkIBConnTailZLen;
  TGeoBBox *connBody = new TGeoBBox("connBodyA", xlen/2, ylen/2, zlen/2);


  TGeoTube *connRoundHole = new TGeoTube("connRoundHoleA", 0.,
					 fgkIBConnRoundHoleD /2,
					 fgkIBConnBodyYHeight/1.5);

  zpos = -connBody->GetDZ() + fgkIBConnRoundHoleZ;
  TGeoCombiTrans *connRoundHoleTrans = new TGeoCombiTrans ("roundHoleTransA",
					 0, 0, zpos,
					   new TGeoRotation("",0,90,0));
  connRoundHoleTrans->RegisterYourself();


  xlen = fgkIBConnSquareHoleX/2;
  ylen = fgkIBConnBodyYHeight/1.5;
  zlen = fgkIBConnSquareHoleZ/2;
  TGeoBBox *connSquareHole = new TGeoBBox("connSquareHoleA", xlen, ylen, zlen);

  zpos = -connBody->GetDZ() + fgkIBConnSquareHoleZPos;
  TGeoTranslation *connSquareHoleTrans = new TGeoTranslation(
					 "squareHoleTransA", 0, 0, zpos);
  connSquareHoleTrans->RegisterYourself();


  TGeoTube *connTubeHole2 = new TGeoTube("tube2HoleA", 0, fgkIBConnTubeHole2D/2, connBody->GetDZ());

  xpos = fgkIBConnTubesXDist/2;
  ypos = -connBody->GetDY() + fgkIBConnTubesYPos;
  
  TGeoTranslation *connTubes2Trans1 = new TGeoTranslation("tubes2Trans1A",-xpos, ypos, 0);
  connTubes2Trans1->RegisterYourself();
  
  TGeoTranslation *connTubes2Trans2 = new TGeoTranslation("tubes2Trans2A", xpos, ypos, 0);
  connTubes2Trans2->RegisterYourself();


  zlen = fgkIBConnTubeHole1ZLen - fgkIBConnTailZLen;
  TGeoTube *connTubeHole3 = new TGeoTube("tube3HoleA", 0, fgkIBConnTubeHole1D/2, zlen);

  zpos = connBody->GetDZ();
  TGeoTranslation *connTubes3Trans1 = new TGeoTranslation("tubes3Trans1A", -xpos, ypos,-zpos);
  connTubes3Trans1->RegisterYourself();
  TGeoTranslation *connTubes3Trans2 = new TGeoTranslation("tubes3Trans2A", xpos, ypos,-zpos);
  connTubes3Trans2->RegisterYourself();

  zlen = fgkIBConnectAFitZLen - fgkIBConnectAFitZOut;
  TGeoTube *connFitHole = new TGeoTube("fitHoleA", 0, fgkIBConnectAFitExtD/2, zlen);

  TGeoTranslation *connFitHoleTrans1 = new TGeoTranslation("fitTrans1A", -xpos, ypos, zpos);
  connFitHoleTrans1->RegisterYourself();
  TGeoTranslation *connFitHoleTrans2 = new TGeoTranslation("fitTrans2A", xpos, ypos, zpos);
  connFitHoleTrans2->RegisterYourself();


  TGeoCompositeShape *connBodySh = new TGeoCompositeShape(
   "connBodyA-connRoundHoleA:roundHoleTransA-connSquareHoleA:squareHoleTransA-tube2HoleA:tubes2Trans1A-tube2HoleA:tubes2Trans2A-fitHoleA:fitTrans1A-fitHoleA:fitTrans2A-tube3HoleA:tubes3Trans1A-tube3HoleA:tubes3Trans2A");


  TGeoVolume *connBlockBody = new TGeoVolume("IBConnectorBlockBodyASide",
					     connBodySh,medPEEK);
  connBlockBody->SetFillColor(42);  // Brownish shade
  connBlockBody->SetLineColor(42);

  // ...and the tail
  xv[0] = fgkIBConnectorXWidth/2;
  yv[0] = fgkIBConnTailYShift;
  xv[1] = xv[0];
  yv[1] = fgkIBConnTailYMid;
  xv[2] = xv[1] -
      (fgkIBConnectorYTot - fgkIBConnTailYMid)/TanD(90-fgkIBConnTailOpenPhi/2);
  yv[2] = fgkIBConnectorYTot;

  for (Int_t i = 0; i<3; i++) {
    xv[3+i] = -xv[2-i];
    yv[3+i] =  yv[2-i];
  }

  TGeoXtru *connTail = new TGeoXtru(2);
  connTail->SetName("connTailA");
  connTail->DefinePolygon(6, xv, yv);
  connTail->DefineSection(0, 0);
  connTail->DefineSection(1, fgkIBConnTailZLen);


  TGeoTube *connTubeHole1 = new TGeoTube("tube1HoleA", 0,
					 fgkIBConnTubeHole1D/2,
					 fgkIBConnTubeHole1ZLen/1.5);

  xpos = fgkIBConnTubesXDist/2;
  ypos = fgkIBConnTubesYPos;
  zpos = connTail->GetZ(1)/2;
  TGeoTranslation *connTubes1Trans1 = new TGeoTranslation("tubes1Trans1A",
							  -xpos, ypos, zpos);
  connTubes1Trans1->RegisterYourself();
  TGeoTranslation *connTubes1Trans2 = new TGeoTranslation("tubes1Trans2A",
							   xpos, ypos, zpos);
  connTubes1Trans2->RegisterYourself();


  TGeoCompositeShape *connTailSh = new TGeoCompositeShape(
		"connTailA-tube1HoleA:tubes1Trans1A-tube1HoleA:tubes1Trans2A");


  TGeoVolume *connBlockTail = new TGeoVolume("IBConnectorBlockTailASide",
					     connTailSh,medPEEK);
  connBlockTail->SetFillColor(42);  // Brownish shade
  connBlockTail->SetLineColor(42);


  // The steel insert, an Xtru
  xv[0] = (fgkIBConnSquareHoleX - fgkIBConnInsertInnerX)/2;
  yv[0] =  fgkIBConnSquareHoleZ/2 - fgkIBConnInsertZThick;
  xv[1] = xv[0];
  yv[1] = -fgkIBConnSquareHoleZ/2;
  xv[2] =  fgkIBConnSquareHoleX/2;
  yv[2] = yv[1];
  xv[3] = xv[2];
  yv[3] =  fgkIBConnSquareHoleZ/2;

  for (Int_t i = 0; i<nv/2; i++) {
    xv[4+i] = -xv[3-i];
    yv[4+i] =  yv[3-i];
  }

  TGeoXtru *connInsertSh = new TGeoXtru(2);
  connInsertSh->DefinePolygon(nv, xv, yv);
  connInsertSh->DefineSection(0,-fgkIBConnInsertHeight/2);
  connInsertSh->DefineSection(1, fgkIBConnInsertHeight/2);

  TGeoVolume *connInsert = new TGeoVolume("IBConnectorInsertASide",
					  connInsertSh, medWC);
  connInsert->SetFillColor(kGray);
  connInsert->SetLineColor(kGray);


  // The fitting tubes, a Tube
  TGeoTube *connFitSh = new TGeoTube(fgkIBConnectAFitIntD/2,
				     fgkIBConnectAFitExtD/2,
				     fgkIBConnectAFitZLen/2);

  TGeoVolume *connFit = new TGeoVolume("IBConnectorFitting",
				       connFitSh, medInox304);
  connFit->SetFillColor(kGray);
  connFit->SetLineColor(kGray);


  // Now create the container: cannot be a simple box
  // to avoid fake overlaps with stave elements
  xlen = fgkIBConnectorXWidth;
  ylen = fgkIBConnectorYTot;
  zlen = fgkIBConnectBlockZLen - fgkIBConnTailZLen + fgkIBConnectAFitZOut;

  TGeoBBox *connBox = new TGeoBBox("connBoxA", xlen/2, ylen/2, zlen/2);

  ypos = -connBox->GetDY();
  zpos = -connBox->GetDZ() - connTail->GetZ(1);
  TGeoTranslation *transTailA = new TGeoTranslation("transTailA",
						     0, ypos, zpos);
  transTailA->RegisterYourself();

  TGeoTube *connTubeHollow = new TGeoTube("tubeHollowA", 0,
					 fgkIBConnTubeHole1D/2,
					 fgkIBConnTubeHole1ZLen/2);

  xpos = fgkIBConnTubesXDist/2;
  ypos = -connBox->GetDY() + fgkIBConnTubesYPos;
  zpos = -connBox->GetDZ() - connTail->GetZ(1) + fgkIBConnTubeHole1ZLen/2;
  TGeoTranslation *connTubeHollTrans1 = new TGeoTranslation("tubeHollTrans1A",
							    -xpos, ypos, zpos);
  connTubeHollTrans1->RegisterYourself();
  TGeoTranslation *connTubeHollTrans2 = new TGeoTranslation("tubeHollTrans2A",
							     xpos, ypos, zpos);
  connTubeHollTrans2->RegisterYourself();

  TGeoCompositeShape *connBoxSh = new TGeoCompositeShape(
      "connBoxA+connTailA:transTailA-tubeHollowA:tubeHollTrans1A-tubeHollowA:tubeHollTrans2A");

  TGeoVolume *connBoxASide = new TGeoVolume("IBConnectorASide",
					    connBoxSh, medAir);


  // Finally build up the connector
  // (NB: the origin is in the connBox, i.e. w/o the tail in Z)
  ypos = -connBox->GetDY();
  zpos = -connBox->GetDZ() - connTail->GetZ(1);
  connBoxASide->AddNode(connBlockTail, 1, new TGeoTranslation(0, ypos, zpos));

  ypos = -connBox->GetDY() + connBody->GetDY();
  zpos = -connBox->GetDZ() + connBody->GetDZ();
  connBoxASide->AddNode(connBlockBody, 1, new TGeoTranslation(0, ypos, zpos));

  zpos = -connBox->GetDZ() + fgkIBConnSquareHoleZPos;
  connBoxASide->AddNode(connInsert, 1, new TGeoCombiTrans(0, ypos, zpos,
					   new TGeoRotation("",0,-90,0)));

  xpos = fgkIBConnTubesXDist/2;
  ypos = -connBox->GetDY() + fgkIBConnTubesYPos;
  zpos =  connBox->GetDZ() - connFitSh->GetDz();
  connBoxASide->AddNode(connFit, 1, new TGeoTranslation( xpos, ypos, zpos));
  connBoxASide->AddNode(connFit, 2, new TGeoTranslation(-xpos, ypos, zpos));

}

//________________________________________________________________________
void AliITSUv2Layer::CreateIBConnectorsCSide(const TGeoManager *mgr){
//
// Create the C-Side end-stave connectors for IB staves
//
// Input:
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      05 May 2015  Mario Sitta
//

  // Local variables
  const Int_t nv = 8;
  Double_t xv[nv], yv[nv];
  Double_t xlen, ylen, zlen;
  Double_t xpos, ypos, zpos;
  char volname[30];


  // Gather all material pointers
  TGeoMedium *medAir      = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medPEEK     = mgr->GetMedium("ITS_PEEKCF30$");
  TGeoMedium *medWC       = mgr->GetMedium("ITS_TUNGCARB$");


  // First create all elements

  // The connector block, two Composite Shapes:
  // the body...
  xlen = fgkIBConnectorXWidth;
  ylen = fgkIBConnBodyYHeight;
  zlen = fgkIBConnectBlockZLen - fgkIBConnTailZLen;
  TGeoBBox *connBody = new TGeoBBox("connBodyC", xlen/2, ylen/2, zlen/2);


  TGeoTube *connRoundHole = new TGeoTube("connRoundHoleC", 0.,
					 fgkIBConnRoundHoleD /2,
					 fgkIBConnBodyYHeight/1.5);

  zpos = -connBody->GetDZ() + fgkIBConnRoundHoleZ;
  TGeoCombiTrans *connRoundHoleTrans = new TGeoCombiTrans ("roundHoleTransC",
					 0, 0, zpos,
					   new TGeoRotation("",0,90,0));
  connRoundHoleTrans->RegisterYourself();


  TGeoTube *connInsertHole = new TGeoTube("connInsertHoleC", 0,
					  fgkIBConnInsertHoleD/2,
					  fgkIBConnBodyYHeight/1.5);

  zpos = -connBody->GetDZ() + fgkIBConnInsertHoleZPos;
  TGeoCombiTrans *connInsertHoleTrans = new TGeoCombiTrans(
					   "insertHoleTransC", 0, 0, zpos,
					      new TGeoRotation("",0,90,0));
  connInsertHoleTrans->RegisterYourself();


  TGeoTube *connTubeHole2 = new TGeoTube("tube2HoleC", 0,
					 fgkIBConnTubeHole2D/2,
					 connBody->GetDZ());

  xpos = fgkIBConnTubesXDist/2;
  ypos = -connBody->GetDY() + fgkIBConnTubesYPos;
  zpos = fgkIBConnTubeHole3ZPos;
  TGeoTranslation *connTubes2Trans1 = new TGeoTranslation("tubes2Trans1C",
							  -xpos, ypos,-zpos);
  connTubes2Trans1->RegisterYourself();
  TGeoTranslation *connTubes2Trans2 = new TGeoTranslation("tubes2Trans2C",
							   xpos, ypos,-zpos);
  connTubes2Trans2->RegisterYourself();


  zlen = fgkIBConnectorXWidth;
  TGeoTube *connTubeHole3 = new TGeoTube("tube3HoleC", 0,
					 fgkIBConnTubeHole2D/2,
					 zlen/2);

  xpos = fgkIBConnTubeHole3XPos;
  zpos = connBody->GetDZ() - fgkIBConnTubeHole3ZPos;
  TGeoCombiTrans *connTubes3Trans = new TGeoCombiTrans("tubes3TransC",
						        xpos, ypos, zpos,
					       new TGeoRotation("",90,-90,90));
  connTubes3Trans->RegisterYourself();


  zlen = fgkIBConnTubeHole1ZLen - fgkIBConnTailZLen;
  TGeoTube *connTubeHole4 = new TGeoTube("tube4HoleC", 0,
					 fgkIBConnTubeHole1D/2,
					 zlen);

  xpos = fgkIBConnTubesXDist/2;
  zpos = connBody->GetDZ();
  TGeoTranslation *connTubes4Trans1 = new TGeoTranslation("tubes4Trans1C",
							  -xpos, ypos,-zpos);
  connTubes4Trans1->RegisterYourself();
  TGeoTranslation *connTubes4Trans2 = new TGeoTranslation("tubes4Trans2C",
							   xpos, ypos,-zpos);
  connTubes4Trans2->RegisterYourself();


  TGeoCompositeShape *connBodySh = new TGeoCompositeShape(
   "connBodyC-connRoundHoleC:roundHoleTransC-connInsertHoleC:insertHoleTransC-tube2HoleC:tubes2Trans1C-tube2HoleC:tubes2Trans2C-tube3HoleC:tubes3TransC-tube4HoleC:tubes4Trans1C-tube4HoleC:tubes4Trans2C");


  TGeoVolume *connBlockBody = new TGeoVolume("IBConnectorBlockBodyCSide",
					     connBodySh,medPEEK);
  connBlockBody->SetFillColor(42);  // Brownish shade
  connBlockBody->SetLineColor(42);

  // ...and the tail
  xv[0] = fgkIBConnectorXWidth/2;
  yv[0] = fgkIBConnTailYShift;
  xv[1] = xv[0];
  yv[1] = fgkIBConnTailYMid;
  xv[2] = xv[1] -
      (fgkIBConnectorYTot - fgkIBConnTailYMid)/TanD(90-fgkIBConnTailOpenPhi/2);
  yv[2] = fgkIBConnectorYTot;

  for (Int_t i = 0; i<3; i++) {
    xv[3+i] = -xv[2-i];
    yv[3+i] =  yv[2-i];
  }

  TGeoXtru *connTail = new TGeoXtru(2);
  connTail->SetName("connTailC");
  connTail->DefinePolygon(6, xv, yv);
  connTail->DefineSection(0, 0);
  connTail->DefineSection(1, fgkIBConnTailZLen);


  TGeoTube *connTubeHole1 = new TGeoTube("tube1HoleC", 0,
					 fgkIBConnTubeHole1D/2,
					 fgkIBConnTubeHole1ZLen/1.5);

  xpos = fgkIBConnTubesXDist/2;
  ypos = fgkIBConnTubesYPos;
  zpos = connTail->GetZ(1)/2;
  TGeoTranslation *connTubes1Trans1 = new TGeoTranslation("tubes1Trans1C",
							  -xpos, ypos, zpos);
  connTubes1Trans1->RegisterYourself();
  TGeoTranslation *connTubes1Trans2 = new TGeoTranslation("tubes1Trans2C",
							   xpos, ypos, zpos);
  connTubes1Trans2->RegisterYourself();


  TGeoCompositeShape *connTailSh = new TGeoCompositeShape(
		"connTailC-tube1HoleC:tubes1Trans1C-tube1HoleC:tubes1Trans2C");


  TGeoVolume *connBlockTail = new TGeoVolume("IBConnectorBlockTailCSide",
					     connTailSh,medPEEK);
  connBlockTail->SetFillColor(42);  // Brownish shade
  connBlockTail->SetLineColor(42);


  // The steel insert, an Tube
  TGeoTube *connInsertSh = new TGeoTube(fgkIBConnInsertD/2,
					fgkIBConnInsertHoleD/2,
					fgkIBConnInsertHeight/2);

  TGeoVolume *connInsert = new TGeoVolume("IBConnectorInsertCSide",
					  connInsertSh, medWC);
  connInsert->SetFillColor(kGray);
  connInsert->SetLineColor(kGray);


  // The plug, a Pcon
  TGeoPcon *connPlugSh = new TGeoPcon(0,360,4);
  connPlugSh->DefineSection(0,                 0., 0., fgkIBConnTubeHole2D/2);
  connPlugSh->DefineSection(1, fgkIBConnPlugThick, 0., fgkIBConnTubeHole2D/2);
  connPlugSh->DefineSection(2, fgkIBConnPlugThick,
			        fgkIBConnPlugInnerD/2, fgkIBConnTubeHole2D/2);
  connPlugSh->DefineSection(3, fgkIBConnPlugTotLen,
			        fgkIBConnPlugInnerD/2, fgkIBConnTubeHole2D/2);

  TGeoVolume *connPlug = new TGeoVolume("IBConnectorPlugC",
					connPlugSh,medPEEK);
  connPlug->SetFillColor(44);  // Brownish shade (a bit darker to spot it)
  connPlug->SetLineColor(44);


  // Now create the container: cannot be a simple box
  // to avoid fake overlaps with stave elements
  xlen = fgkIBConnectorXWidth;
  ylen = fgkIBConnectorYTot;
  zlen = fgkIBConnectBlockZLen - fgkIBConnTailZLen;

  TGeoBBox *connBox = new TGeoBBox("connBoxC", xlen/2, ylen/2, zlen/2);

  ypos = -connBox->GetDY();
  zpos = -connBox->GetDZ() - connTail->GetZ(1);
  TGeoTranslation *transTailC = new TGeoTranslation("transTailC",
						     0, ypos, zpos);
  transTailC->RegisterYourself();

  TGeoTube *connTubeHollow = new TGeoTube("tubeHollowC", 0,
					 fgkIBConnTubeHole1D/2,
					 fgkIBConnTubeHole1ZLen/2);

  xpos = fgkIBConnTubesXDist/2;
  ypos = -connBox->GetDY() + fgkIBConnTubesYPos;
  zpos = -connBox->GetDZ() - connTail->GetZ(1) + fgkIBConnTubeHole1ZLen/2;
  TGeoTranslation *connTubeHollTrans1 = new TGeoTranslation("tubeHollTrans1C",
							    -xpos, ypos, zpos);
  connTubeHollTrans1->RegisterYourself();
  TGeoTranslation *connTubeHollTrans2 = new TGeoTranslation("tubeHollTrans2C",
							     xpos, ypos, zpos);
  connTubeHollTrans2->RegisterYourself();

  TGeoCompositeShape *connBoxSh = new TGeoCompositeShape(
      "connBoxC+connTailC:transTailC-tubeHollowC:tubeHollTrans1C-tubeHollowC:tubeHollTrans2C");

  TGeoVolume *connBoxCSide = new TGeoVolume("IBConnectorCSide",
					    connBoxSh, medAir);


  // Finally build up the connector
  // (NB: the origin is in the connBox, i.e. w/o the tail in Z)
  ypos = -connBoxSh->GetDY();
  zpos = -connBodySh->GetDZ() - connTail->GetZ(1);
  connBoxCSide->AddNode(connBlockTail, 1, new TGeoTranslation(0, ypos, zpos));

  ypos = -connBoxSh->GetDY() + connBodySh->GetDY();
  connBoxCSide->AddNode(connBlockBody, 1, new TGeoTranslation(0, ypos, 0));

  zpos = -connBox->GetDZ() + fgkIBConnInsertHoleZPos;
  connBoxCSide->AddNode(connInsert, 1, new TGeoCombiTrans(0, ypos, zpos,
					   new TGeoRotation("",0,90,0)));

  xpos =  connBox->GetDX();
  ypos = -connBox->GetDY() + fgkIBConnTubesYPos;
  zpos =  connBox->GetDZ() - fgkIBConnTubeHole3ZPos;;
  connBoxCSide->AddNode(connPlug, 1, new TGeoCombiTrans(xpos, ypos, zpos,
					   new TGeoRotation("",90,-90,90)));

}

//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateStaveOuterB(const TGeoManager *mgr){
//
// Create the chip stave for the Outer Barrel
//
// Input:
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      20 Dec 2013  Mario Sitta
// Updated:      12 Mar 2014  Mario Sitta
//

  TGeoVolume *mechStavVol = 0;

  switch (fStaveModel) {
    case AliITSUv2::kOBModelDummy:
      mechStavVol = CreateStaveModelOuterBDummy(mgr);
      break;
    case AliITSUv2::kOBModel0:
      mechStavVol = CreateStaveModelOuterB0(mgr);
      break;
    case AliITSUv2::kOBModel1:
    case AliITSUv2::kOBModel2:
      mechStavVol = CreateStaveModelOuterB12(mgr);
      break;
    default:
      AliFatal(Form("Unknown stave model %d",fStaveModel));
      break;
  }

  return mechStavVol; 
}

//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateStaveModelOuterBDummy(const TGeoManager *) const {
//
// Create dummy stave
//
// Input:
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      20 Dec 2013  Mario Sitta
//


  // Done, return the stave structure
  return 0;
}

//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateStaveModelOuterB0(const TGeoManager *mgr){
//
// Creation of the mechanical stave structure for the Outer Barrel as in v0
// (we fake the module and halfstave volumes to have always
// the same formal geometry hierarchy)
//
// Input:
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      20 Dec 2013  Mario Sitta
// Updated:      12 Mar 2014  Mario Sitta
//

  // Local variables
  Double_t xmod, ymod, zmod;
  Double_t xlen, ylen, zlen;
  Double_t ypos, zpos;
  char volname[30];
  
  // First create all needed shapes

  // The chip
  xlen = fgkOBHalfStaveWidth;
  ylen = 0.5*fChipThick;
  zlen = fgkOBModuleZLength/2;

  TGeoVolume *chipVol = CreateChipInnerB(xlen, ylen, zlen);

  xmod = ((TGeoBBox*)chipVol->GetShape())->GetDX();
  ymod = ((TGeoBBox*)chipVol->GetShape())->GetDY();
  zmod = ((TGeoBBox*)chipVol->GetShape())->GetDZ();

  TGeoBBox *module = new TGeoBBox(xmod, ymod, zmod);

  zlen = fgkOBModuleZLength*fNModules;
  TGeoBBox *hstave = new TGeoBBox(xlen, ylen, zlen/2);


  // We have all shapes: now create the real volumes

  TGeoMedium *medAir = mgr->GetMedium("ITS_AIR$");

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSModulePattern(), fLayerNumber);
  TGeoVolume *modVol = new TGeoVolume(volname, module, medAir);
  modVol->SetVisibility(kTRUE);

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSHalfStavePattern(), fLayerNumber);
  TGeoVolume *hstaveVol  = new TGeoVolume(volname, hstave, medAir);


  // Finally build it up
  modVol->AddNode(chipVol, 0);
  fHierarchy[kChip]=1;

  for (Int_t j=0; j<fNModules; j++) {
    ypos = 0.021;  // Remove small overlap - M.S: 21may13
    zpos = -hstave->GetDZ() + j*2*zmod + zmod;
    hstaveVol->AddNode(modVol, j, new TGeoTranslation(0, ypos, zpos));
    fHierarchy[kModule]++;
  }


  // Done, return the stave structure
  return hstaveVol;
}

//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateStaveModelOuterB12(const TGeoManager *mgr){
//
// Create the mechanical half stave structure
// for the Outer Barrel as in TDR
//
// Input:
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      20 Nov 2013  Anastasia Barbano
// Updated:      16 Jan 2014  Mario Sitta
// Updated:      24 Feb 2014  Mario Sitta
// Updated:      11 Nov 2014  Mario Sitta  Model2
// Updated:      03 Dec 2014  Mario Sitta  Revised with C.Gargiulo latest infos
//


  // Local parameters
  Double_t yFlex1      = fgkOBFlexCableAlThick;
  Double_t yFlex2      = fgkOBFlexCableKapThick;
  Double_t flexOverlap = 5; // to be checked - unused for the time being
  Double_t xHalfSt     = fgkOBHalfStaveWidth/2;
  Double_t yCFleece    = fgkOBCarbonFleeceThick;
  Double_t yGraph      = fgkOBGraphiteFoilThick;
  Double_t yHalfSt     = 0; // will be computed

  Double_t ymod, zmod;
  Double_t xtru[12], ytru[12];
  Double_t xpos, ypos, ypos1, zpos/*, zpos5cm*/;
  Double_t xlen, ylen, zlen;
  char volname[30];

  Double_t rCoolMin, rCoolMax;
  if (fStaveModel == AliITSUv2::kOBModel1)
    rCoolMin = fgkOBCoolTubeInnerDM1/2;
  else
    rCoolMin = fgkOBCoolTubeInnerD/2;

  rCoolMax = rCoolMin + fgkOBCoolTubeThick;

  zlen = (fNModules*fgkOBModuleZLength + (fNModules-1)*fgkOBModuleGap)/2;


  // First create all needed shapes

  TGeoVolume *moduleVol = CreateModuleOuterB();
  moduleVol->SetVisibility(kTRUE);
  ymod = ((TGeoBBox*)(moduleVol->GetShape()))->GetDY();
  zmod = ((TGeoBBox*)(moduleVol->GetShape()))->GetDZ();

  TGeoBBox *busAl  = new TGeoBBox("BusAl",  xHalfSt, fgkOBBusCableAlThick/2,
					    zlen);
  TGeoBBox *busKap = new TGeoBBox("BusKap", xHalfSt, fgkOBBusCableKapThick/2,
					    zlen);

  TGeoBBox *glue = new TGeoBBox("Glue", xHalfSt, fgkOBGlueThick/2, zlen);

  TGeoBBox *coldPlate  = new TGeoBBox("ColdPlate", fgkOBHalfStaveWidth/2,
				      fgkOBColdPlateThick/2, zlen);

  TGeoBBox *fleeccent  = new TGeoBBox("FleeceCent",  xHalfSt,
				      yCFleece/2, zlen);

  TGeoTube *coolTube   = new TGeoTube("CoolingTube", rCoolMin, rCoolMax, zlen);
  TGeoTube *coolWater  = new TGeoTube("CoolingWater", 0., rCoolMin, zlen);

  xlen = xHalfSt - fgkOBCoolTubeXDist/2 - coolTube->GetRmax();
  TGeoBBox *graphlat   = new TGeoBBox("GraphLateral", xlen/2, yGraph/2, zlen);

  xlen = fgkOBCoolTubeXDist/2 - coolTube->GetRmax();
  TGeoBBox *graphmid   = new TGeoBBox("GraphMiddle", xlen, yGraph/2, zlen);

  ylen = coolTube->GetRmax() - yGraph;
  TGeoBBox *graphvert  = new TGeoBBox("GraphVertical", yGraph/2, ylen/2, zlen);

  TGeoTubeSeg *graphtub = new TGeoTubeSeg("GraphTube", rCoolMax,
					  rCoolMax+yGraph, zlen,
					  180., 360.);

  xlen = xHalfSt - fgkOBCoolTubeXDist/2 - coolTube->GetRmax() - yGraph;
  TGeoBBox *fleeclat  = new TGeoBBox("FleecLateral", xlen/2, yCFleece/2, zlen);

  xlen = fgkOBCoolTubeXDist/2 - coolTube->GetRmax() - yGraph;
  TGeoBBox *fleecmid  = new TGeoBBox("FleecMiddle", xlen, yCFleece/2, zlen);

  ylen = coolTube->GetRmax() - yGraph - yCFleece;
  TGeoBBox *fleecvert = new TGeoBBox("FleecVertical", yCFleece/2, ylen/2,
						      zlen);

  TGeoTubeSeg *fleectub = new TGeoTubeSeg("FleecTube", rCoolMax+yGraph,
					  rCoolMax+yCFleece+yGraph,
					  zlen, 180., 360.);

  TGeoTube *gammaConvRod;
  if (fAddGammaConv)
    gammaConvRod = new TGeoTube("GammaConver", 0, 0.5*fGammaConvDiam, zlen);

  TGeoBBox *flex1_5cm  = new TGeoBBox("Flex1MV_5cm",xHalfSt,yFlex1/2,flexOverlap/2);
  TGeoBBox *flex2_5cm  = new TGeoBBox("Flex2MV_5cm",xHalfSt,yFlex2/2,flexOverlap/2);

  // The half stave container (an XTru to avoid overlaps between neightbours)
  yHalfSt = ymod + busAl->GetDY() + busKap->GetDY() + coldPlate->GetDY()
	  + fleeccent->GetDY() + graphlat->GetDY() + fleeclat->GetDY();
  if (fStaveModel == AliITSUv2::kOBModel2)
    yHalfSt += 2*glue->GetDY();
//IB  if (fAddGammaConv)
//IB   yHalfSt += fGammaConvDiam;

  xtru[0] = xHalfSt;
  ytru[0] = 0;
  xtru[1] = xtru[0];
  ytru[1] = -2*yHalfSt;
  xtru[2] = fgkOBCoolTubeXDist/2 + fleectub->GetRmax();
  ytru[2] = ytru[1];
  xtru[3] = xtru[2];
  ytru[3] = ytru[2] - (coolTube->GetRmax() + fleectub->GetRmax());
  xtru[4] = fgkOBCoolTubeXDist/2 - fleectub->GetRmax();
  ytru[4] = ytru[3];
  xtru[5] = xtru[4];
  ytru[5] = ytru[2];
  for (Int_t i = 0; i<6; i++) {
    xtru[6+i] = -xtru[5-i];
    ytru[6+i] =  ytru[5-i];
  }
  TGeoXtru *halfStaveCent = new TGeoXtru(2);
  halfStaveCent->DefinePolygon(12, xtru, ytru);
  halfStaveCent->DefineSection(0,-fZLength/2);
  halfStaveCent->DefineSection(1, fZLength/2);
  snprintf(volname, 30, "staveCentral%d", fLayerNumber);
  halfStaveCent->SetName(volname);

  // The connectors' containers
  TGeoBBox *connAside = new TGeoBBox("connAsideOB",fgkOBCPConnectorXWidth/2,
				     fgkOBCPConnBlockYHei/2,
			     (fgkOBCPConnBlockZLen + fgkOBCPConnAFitZOut)/2);

  TGeoBBox *connCside = new TGeoBBox("connCsideOB",fgkOBCPConnectorXWidth/2,
				     fgkOBCPConnBlockYHei/2,
				     fgkOBCPConnBlockZLen/2);

  // The StaveStruct container, a Composite Shape
  ypos = 2*yHalfSt + connAside->GetDY() - fgkOBCPConnHollowYHei;
  zpos = zlen + connAside->GetDZ() - fgkOBCPConnHollowZLen;
  snprintf(volname, 30, "transAsideOB%d", fLayerNumber);
  TGeoTranslation *transAside = new TGeoTranslation(volname, 0,-ypos, zpos);
  transAside->RegisterYourself();

  zpos = zlen + connCside->GetDZ() - fgkOBCPConnHollowZLen;
  snprintf(volname, 30, "transCsideOB%d", fLayerNumber);
  TGeoTranslation *transCside = new TGeoTranslation(volname, 0,-ypos,-zpos);
  transCside->RegisterYourself();

  char componame[70];
  snprintf(componame, 70,
     "staveCentral%d+connAsideOB:transAsideOB%d+connCsideOB:transCsideOB%d",
	   fLayerNumber,fLayerNumber,fLayerNumber);

  TGeoCompositeShape *halfStave = new TGeoCompositeShape(componame);


  // We have all shapes: now create the real volumes

  TGeoMedium *medAluminum     = mgr->GetMedium("ITS_ALUMINUM$");
  TGeoMedium *medK13D2U120    = mgr->GetMedium("ITS_K13D2U120$");
  TGeoMedium *medKapton       = mgr->GetMedium("ITS_KAPTON(POLYCH2)$");
  TGeoMedium *medWater        = mgr->GetMedium("ITS_WATER$");
  TGeoMedium *medCarbonFleece = mgr->GetMedium("ITS_CarbonFleece$"); 
  TGeoMedium *medFGS003       = mgr->GetMedium("ITS_FGS003$"); //amec thermasol
  TGeoMedium *medAir          = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medGlue         = mgr->GetMedium("ITS_GLUE$");
  TGeoMedium *medTungsten     = mgr->GetMedium("ITS_TUNGSTEN$");


  TGeoVolume *busAlVol = new TGeoVolume("PowerBusAlVol", busAl , medAluminum);
  busAlVol->SetLineColor(kCyan);
  busAlVol->SetFillColor(busAlVol->GetLineColor());
  busAlVol->SetFillStyle(4000); // 0% transparent

  TGeoVolume *busKapVol = new TGeoVolume("PowerBusKapVol", busKap, medKapton);
  busKapVol->SetLineColor(kBlue);
  busKapVol->SetFillColor(busKapVol->GetLineColor());
  busKapVol->SetFillStyle(4000); // 0% transparent

  TGeoVolume *coldPlateVol = new TGeoVolume("ColdPlateVol",
					    coldPlate, medK13D2U120);
  coldPlateVol->SetLineColor(kYellow-3);
  coldPlateVol->SetFillColor(coldPlateVol->GetLineColor());
  coldPlateVol->SetFillStyle(4000); // 0% transparent

  TGeoVolume *fleeccentVol = new TGeoVolume("CarbonFleeceCentral",
					    fleeccent, medCarbonFleece);
  fleeccentVol->SetLineColor(kViolet);
  fleeccentVol->SetFillColor(fleeccentVol->GetLineColor());
  fleeccentVol->SetFillStyle(4000); // 0% transparent

  TGeoVolume *glueVol = new TGeoVolume("GlueVol", glue, medGlue);
  glueVol->SetLineColor(kBlack);
  glueVol->SetFillColor(glueVol->GetLineColor());
  glueVol->SetFillStyle(4000); // 0% transparent

  TGeoVolume *coolTubeVol  = new TGeoVolume("CoolingTubeVol",
					    coolTube, medKapton);
  coolTubeVol->SetLineColor(kGray);
  coolTubeVol->SetFillColor(coolTubeVol->GetLineColor());
  coolTubeVol->SetFillStyle(4000); // 0% transparent

  TGeoVolume *coolWaterVol;
  coolWaterVol = new TGeoVolume("CoolingWaterVol", coolWater, medWater);
  coolWaterVol->SetLineColor(kBlue);
  coolWaterVol->SetFillColor(coolWaterVol->GetLineColor());
  coolWaterVol->SetFillStyle(4000); // 0% transparent

  TGeoVolume *graphlatVol = new TGeoVolume("GraphiteFoilLateral",
					   graphlat, medFGS003);
  graphlatVol->SetLineColor(kGreen);
  graphlatVol->SetFillColor(graphlatVol->GetLineColor());
  graphlatVol->SetFillStyle(4000); // 0% transparent

  TGeoVolume *graphmidVol = new TGeoVolume("GraphiteFoilMiddle",
					   graphmid, medFGS003);
  graphmidVol->SetLineColor(kGreen);
  graphmidVol->SetFillColor(graphmidVol->GetLineColor());
  graphmidVol->SetFillStyle(4000); // 0% transparent

  TGeoVolume *graphvertVol = new TGeoVolume("GraphiteFoilVertical",
					    graphvert, medFGS003);
  graphvertVol->SetLineColor(kGreen);
  graphvertVol->SetFillColor(graphvertVol->GetLineColor());
  graphvertVol->SetFillStyle(4000); // 0% transparent

  TGeoVolume *graphtubVol = new TGeoVolume("GraphiteFoilPipeCover",
					   graphtub, medFGS003);
  graphtubVol->SetLineColor(kGreen);
  graphtubVol->SetFillColor(graphtubVol->GetLineColor());
  graphtubVol->SetFillStyle(4000); // 0% transparent

  TGeoVolume *fleeclatVol = new TGeoVolume("CarbonFleeceLateral",
					   fleeclat, medCarbonFleece);
  fleeclatVol->SetLineColor(kViolet);
  fleeclatVol->SetFillColor(fleeclatVol->GetLineColor());
  fleeclatVol->SetFillStyle(4000); // 0% transparent

  TGeoVolume *fleecmidVol = new TGeoVolume("CarbonFleeceMiddle",
					   fleecmid, medCarbonFleece);
  fleecmidVol->SetLineColor(kViolet);
  fleecmidVol->SetFillColor(fleecmidVol->GetLineColor());
  fleecmidVol->SetFillStyle(4000); // 0% transparent

  TGeoVolume *fleecvertVol = new TGeoVolume("CarbonFleeceVertical",
					    fleecvert, medCarbonFleece);
  fleecvertVol->SetLineColor(kViolet);
  fleecvertVol->SetFillColor(fleecvertVol->GetLineColor());
  fleecvertVol->SetFillStyle(4000); // 0% transparent

  TGeoVolume *fleectubVol = new TGeoVolume("CarbonFleecePipeCover",
					   fleectub, medCarbonFleece);
  fleectubVol->SetLineColor(kViolet);
  fleectubVol->SetFillColor(fleectubVol->GetLineColor());
  fleectubVol->SetFillStyle(4000); // 0% transparent

  TGeoVolume *gammaConvRodVol;
  if (fAddGammaConv) {
    gammaConvRodVol = new TGeoVolume("GammaConversionRod",
				     gammaConvRod, medTungsten);
    gammaConvRodVol->SetLineColor(kBlack);
    gammaConvRodVol->SetFillColor(gammaConvRodVol->GetLineColor());
    gammaConvRodVol->SetFillStyle(4000); // 0% transparent
  }

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSHalfStavePattern(), fLayerNumber);
  TGeoVolume *halfStaveVol  = new TGeoVolume(volname, halfStave, medAir);
//   halfStaveVol->SetLineColor(12);
//   halfStaveVol->SetFillColor(12); 
//   halfStaveVol->SetVisibility(kTRUE);


  TGeoVolume *flex1_5cmVol = new TGeoVolume("Flex1Vol5cm",flex1_5cm,medAluminum);
  TGeoVolume *flex2_5cmVol = new TGeoVolume("Flex2Vol5cm",flex2_5cm,medKapton);


  flex1_5cmVol->SetLineColor(kRed);
  flex2_5cmVol->SetLineColor(kGreen);
  

  // Now build up the half stave
  ypos = - busKap->GetDY();
  if ( (fStaveModel == AliITSUv2::kOBModel1 && fBuildLevel < 4) || // Kapton
       (fStaveModel == AliITSUv2::kOBModel2 && fBuildLevel < 5) )
    halfStaveVol->AddNode(busKapVol, 1, new TGeoTranslation(0, ypos, 0));

  ypos -= (busKap->GetDY() + busAl->GetDY());
  if ( (fStaveModel == AliITSUv2::kOBModel1 && fBuildLevel < 1) || // Aluminum
       (fStaveModel == AliITSUv2::kOBModel2 && fBuildLevel < 2) )
    halfStaveVol->AddNode(busAlVol, 1, new TGeoTranslation(0, ypos, 0));

  ypos -= busAl->GetDY();

  if (fStaveModel == AliITSUv2::kOBModel2) {
    ypos -= glue->GetDY();
    if (fBuildLevel < 3) // Glue
      halfStaveVol->AddNode(glueVol, 1, new TGeoTranslation(0, ypos, 0));
    ypos -= glue->GetDY();
  }

  ypos -= ymod;
  for (Int_t j=0; j<fNModules; j++) {
    zpos = -zlen + j*(2*zmod + fgkOBModuleGap) + zmod;
    halfStaveVol->AddNode(moduleVol, j, new TGeoTranslation(0, ypos, zpos));
    fHierarchy[kModule]++;
  }

  ypos -= ymod;

  if (fStaveModel == AliITSUv2::kOBModel2) {
    ypos -= glue->GetDY();
    if (fBuildLevel < 3) // Glue
      halfStaveVol->AddNode(glueVol, 2, new TGeoTranslation(0, ypos, 0));
    ypos -= glue->GetDY();
  }

  ypos -= fleeccent->GetDY();
  if ( (fStaveModel == AliITSUv2::kOBModel1 && fBuildLevel < 5) || // Carbon
       (fStaveModel == AliITSUv2::kOBModel2 && fBuildLevel < 6) )
    halfStaveVol->AddNode(fleeccentVol, 1, new TGeoTranslation(0, ypos, 0));
  ypos -= fleeccent->GetDY();

  ypos -= coldPlate->GetDY();
  if ( (fStaveModel == AliITSUv2::kOBModel1 && fBuildLevel < 5) || // Carbon
       (fStaveModel == AliITSUv2::kOBModel2 && fBuildLevel < 6) )
    halfStaveVol->AddNode(coldPlateVol, 1, new TGeoTranslation(0, ypos, 0));

  xpos = fgkOBCoolTubeXDist/2;
  ypos1 = ypos - (coldPlate->GetDY() + coolTube->GetRmax());
  if ( (fStaveModel == AliITSUv2::kOBModel1 && fBuildLevel < 3) || // Water
       (fStaveModel == AliITSUv2::kOBModel2 && fBuildLevel < 4) ) {
    halfStaveVol->AddNode(coolWaterVol, 1, new TGeoTranslation(-xpos, ypos1, 0));
    halfStaveVol->AddNode(coolWaterVol, 2, new TGeoTranslation( xpos, ypos1, 0));
  }

  if ( (fStaveModel == AliITSUv2::kOBModel1 && fBuildLevel < 4) || // Kapton
       (fStaveModel == AliITSUv2::kOBModel2 && fBuildLevel < 5) ) {
    halfStaveVol->AddNode(coolTubeVol,1, new TGeoTranslation(-xpos, ypos1, 0));
    halfStaveVol->AddNode(coolTubeVol,2, new TGeoTranslation( xpos, ypos1, 0));
  }

  if ( (fStaveModel == AliITSUv2::kOBModel1 && fBuildLevel < 5) || // Carbon
       (fStaveModel == AliITSUv2::kOBModel2 && fBuildLevel < 6) ) {
    halfStaveVol->AddNode(graphtubVol,1, new TGeoTranslation(-xpos, ypos1, 0));
    halfStaveVol->AddNode(graphtubVol,2, new TGeoTranslation( xpos, ypos1, 0));

    halfStaveVol->AddNode(fleectubVol,1, new TGeoTranslation(-xpos, ypos1, 0));
    halfStaveVol->AddNode(fleectubVol,2, new TGeoTranslation( xpos, ypos1, 0));
  }

  xpos = xHalfSt - graphlat->GetDX();
  ypos1 = ypos - (coldPlate->GetDY() + graphlat->GetDY());
  if ( (fStaveModel == AliITSUv2::kOBModel1 && fBuildLevel < 5) || // Carbon
       (fStaveModel == AliITSUv2::kOBModel2 && fBuildLevel < 6) ) {
    halfStaveVol->AddNode(graphlatVol,1, new TGeoTranslation(-xpos, ypos1, 0));
    halfStaveVol->AddNode(graphlatVol,2, new TGeoTranslation( xpos, ypos1, 0));

    halfStaveVol->AddNode(graphmidVol,1, new TGeoTranslation(0, ypos1, 0));

    xpos = xHalfSt - 2*graphlat->GetDX() + graphvert->GetDX();
    ypos1 = ypos - (coldPlate->GetDY()+2*graphlat->GetDY()+graphvert->GetDY());
    halfStaveVol->AddNode(graphvertVol,1,new TGeoTranslation(-xpos, ypos1, 0));
    halfStaveVol->AddNode(graphvertVol,2,new TGeoTranslation( xpos, ypos1, 0));
    xpos = graphmid->GetDX() - graphvert->GetDX();
    halfStaveVol->AddNode(graphvertVol,3,new TGeoTranslation(-xpos, ypos1, 0));
    halfStaveVol->AddNode(graphvertVol,4,new TGeoTranslation( xpos, ypos1, 0));
  }

  xpos = xHalfSt - fleeclat->GetDX();
  ypos1 = ypos - (coldPlate->GetDY() +2*graphlat->GetDY() +fleeclat->GetDY());
  if ( (fStaveModel == AliITSUv2::kOBModel1 && fBuildLevel < 5) || // Carbon
       (fStaveModel == AliITSUv2::kOBModel2 && fBuildLevel < 6) ) {
    halfStaveVol->AddNode(fleeclatVol,1, new TGeoTranslation(-xpos, ypos1, 0));
    halfStaveVol->AddNode(fleeclatVol,2, new TGeoTranslation( xpos, ypos1, 0));

    halfStaveVol->AddNode(fleecmidVol, 1, new TGeoTranslation(0, ypos1, 0));

    xpos = xHalfSt - 2*fleeclat->GetDX() + fleecvert->GetDX();
    ypos1 = ypos - (coldPlate->GetDY() +2*graphlat->GetDY()
                 + 2*fleeclat->GetDY() + fleecvert->GetDY());
    halfStaveVol->AddNode(fleecvertVol,1,new TGeoTranslation(-xpos, ypos1, 0));
    halfStaveVol->AddNode(fleecvertVol,2,new TGeoTranslation( xpos, ypos1, 0));
    xpos = fleecmid->GetDX() - fleecvert->GetDX();
    halfStaveVol->AddNode(fleecvertVol,3,new TGeoTranslation(-xpos, ypos1, 0));
    halfStaveVol->AddNode(fleecvertVol,4,new TGeoTranslation( xpos, ypos1, 0));
  }

  // Add the Gamma Converter Rod (only on Layer 3) - M.S. 17 Oct 2016
    if (fAddGammaConv) {
    xpos = fGammaConvXPos;
    ypos1 = ypos - (coldPlate->GetDY() +2*graphlat->GetDY()
                  +2*fleeclat->GetDY() +gammaConvRod->GetRmax());
    halfStaveVol->AddNode(gammaConvRodVol,1,new TGeoTranslation(xpos, ypos1, 0));
  }

  // Add the end-stave connectors
  TGeoVolume *connectorASide, *connectorCSide;

  // Check whether we have already all pieces
  // Otherwise create them
  connectorASide = mgr->GetVolume("OBColdPlateConnectorASide");

  if (!connectorASide) {
    CreateOBColdPlateConnectors();
    connectorASide = mgr->GetVolume("OBColdPlateConnectorASide");
  }
  connectorCSide = mgr->GetVolume("OBColdPlateConnectorCSide");

  ypos = 2*yHalfSt + ((TGeoBBox*)connectorASide->GetShape())->GetDY()
       - fgkOBCPConnHollowYHei;
  zpos = zlen + ((TGeoBBox*)connectorASide->GetShape())->GetDZ()
       - fgkOBCPConnHollowZLen;
  halfStaveVol->AddNode(connectorASide, 1, new TGeoCombiTrans(0,-ypos, zpos,
					      new TGeoRotation("",180,0,0)));

  zpos = zlen + ((TGeoBBox*)connectorCSide->GetShape())->GetDZ()
       - fgkOBCPConnHollowZLen;
  halfStaveVol->AddNode(connectorCSide, 1, new TGeoCombiTrans(0,-ypos,-zpos,
					      new TGeoRotation("",180,0,0)));


  // Done, return the half stave structure
  return halfStaveVol;
}

//________________________________________________________________________
void AliITSUv2Layer::CreateOBColdPlateConnectors(){
//
// Create the Cold Plate connectors for OB half staves
// (simply call the actual creator methods)
//
// Input:
//
// Output:
//
// Return:
//
// Created:      26 May 2015  Mario Sitta
//

  CreateOBColdPlateConnectorsASide();
  CreateOBColdPlateConnectorsCSide();

}


//________________________________________________________________________
void AliITSUv2Layer::CreateOBColdPlateConnectorsASide(){
//
// Create the A-Side end-stave connectors for IB staves
//
// Input:
//
// Output:
//
// Return:
//
// Created:      26 May 2015  Mario Sitta
//

  // The geoManager
  const TGeoManager *mgr = gGeoManager;

  // Local variables
  const Int_t nv = 16;
  Double_t xv[nv], yv[nv];
  Double_t xlen, ylen, zlen;
  Double_t xpos, ypos, zpos;
  char volname[30];


  // Gather all material pointers
  TGeoMedium *medAir      = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medPEEK     = mgr->GetMedium("ITS_PEEKCF30$");
  TGeoMedium *medWC       = mgr->GetMedium("ITS_TUNGCARB$");
  TGeoMedium *medInox304  = mgr->GetMedium("ITS_INOX304$");


  // First create all elements

  // The connector block, a Composite Shape
  xlen = fgkOBCPConnectorXWidth;
  ylen = fgkOBCPConnBlockYHei;
  zlen = fgkOBCPConnBlockZLen;
  TGeoBBox *connBlock = new TGeoBBox("connBlockA", xlen/2, ylen/2, zlen/2);


  xv[0] =  fgkOBCPConnectorXWidth*0.6;
  yv[0] = -fgkOBCPConnHollowYHei;
  xv[1] = xv[0];
  yv[1] =  fgkOBCPConnHollowYHei;
  xv[2] =  fgkOBCPConnTubesXDist/2 + fgkOBCPConnTubeHole1D/2;
  yv[2] = yv[1];
  xv[3] = xv[2];
  yv[3] =  fgkOBCPConnTubesYPos;
  xv[4] =  fgkOBCPConnTubesXDist/2 - fgkOBCPConnTubeHole1D/2;
  yv[4] = yv[3];
  xv[5] = xv[4];
  yv[5] = yv[2];

  for (Int_t i = 0; i<6; i++) {
    xv[6+i] = -xv[5-i];
    yv[6+i] =  yv[5-i];
  }

  TGeoXtru *connBlockHoll = new TGeoXtru(2);
  connBlockHoll->SetName("connBlockHollA");
  connBlockHoll->DefinePolygon(12, xv, yv);
  connBlockHoll->DefineSection(0,-fgkOBCPConnHollowZLen);
  connBlockHoll->DefineSection(1, fgkOBCPConnHollowZLen);

  ypos = -connBlock->GetDY();
  zpos = -connBlock->GetDZ();
  TGeoTranslation *transBlockHoll = new TGeoTranslation("transBlockHollA",
							0, ypos, zpos);
  transBlockHoll->RegisterYourself();


  xlen = fgkOBCPConnSquareHoleX/2;
  ylen = fgkOBCPConnBlockYHei/1.5;
  zlen = fgkOBCPConnSquareHoleZ/2;
  TGeoBBox *connSquareHole = new TGeoBBox("connASquareHole", xlen, ylen, zlen);

  zpos = -connBlock->GetDZ() + ( fgkOBCPConnHollowZLen + fgkOBCPConnSqrHoleZPos
         + fgkOBCPConnSqrInsertRZ - connSquareHole->GetDZ() );
  TGeoTranslation *transSquareHole = new TGeoTranslation(
					 "transASquareHole", 0, 0, zpos);
  transSquareHole->RegisterYourself();


  zlen = fgkOBCPConnTubeHole1Z;
  TGeoTube *connTubeHole1 = new TGeoTube("tube1AHole", 0,
					 fgkOBCPConnTubeHole1D/2,
					 zlen);

  xpos = fgkOBCPConnTubesXDist/2;
  ypos = -connBlock->GetDY() + fgkOBCPConnTubesYPos;
  zpos =  connBlock->GetDZ();
  TGeoTranslation *trans1Tube1AHole = new TGeoTranslation("trans1Tube1AHole",
							  -xpos, ypos,-zpos);
  trans1Tube1AHole->RegisterYourself();
  TGeoTranslation *trans2Tube1AHole = new TGeoTranslation("trans2Tube1AHole",
							   xpos, ypos,-zpos);
  trans2Tube1AHole->RegisterYourself();


  zlen = fgkOBCPConnBlockZLen;
  TGeoTube *connTubeHole2 = new TGeoTube("tube2AHole", 0,
					 fgkOBCPConnTubeHole2D/2,
					 zlen);

  TGeoTranslation *trans1Tube2AHole = new TGeoTranslation("trans1Tube2AHole",
							  -xpos, ypos, 0);
  trans1Tube2AHole->RegisterYourself();
  TGeoTranslation *trans2Tube2AHole = new TGeoTranslation("trans2Tube2AHole",
							   xpos, ypos, 0);
  trans2Tube2AHole->RegisterYourself();


  zlen = fgkOBCPConnAFitZLen - fgkOBCPConnAFitZOut;
  TGeoTube *connFitHole = new TGeoTube("fitAHole", 0,
				       fgkOBCPConnFitHoleD/2,
				       zlen);

  TGeoTranslation *trans1FitAHole = new TGeoTranslation("trans1FitAHole",
							-xpos, ypos, zpos);
  trans1FitAHole->RegisterYourself();
  TGeoTranslation *trans2FitAHole = new TGeoTranslation("trans2FitAHole",
							 xpos, ypos, zpos);
  trans2FitAHole->RegisterYourself();


  TGeoCompositeShape *connBlockSh = new TGeoCompositeShape(
    "connBlockA-connBlockHollA:transBlockHollA-connASquareHole:transASquareHole-tube1AHole:trans1Tube1AHole-tube1AHole:trans2Tube1AHole-tube2AHole:trans1Tube2AHole-tube2AHole:trans2Tube2AHole-fitAHole:trans1FitAHole-fitAHole:trans2FitAHole");


  TGeoVolume *connBlockA = new TGeoVolume("OBColdPlateConnectorBlockASide",
					  connBlockSh,medPEEK);
  connBlockA->SetFillColor(42);  // Brownish shade
  connBlockA->SetLineColor(42);


  // The steel insert, an Xtru
  xv[0] =  fgkOBCPConnSquareHoleX/2;
  yv[0] = -fgkOBCPConnSqrInsertRZ;
  xv[1] = xv[0];
  yv[1] = yv[0] + fgkOBCPConnSquareHoleZ;
  xv[2] = fgkOBCPConnInstInnerX/2;
  yv[2] = yv[1];
  xv[3] = xv[2];
  yv[3] = yv[0] + fgkOBCPConnInstZThick + fgkOBCPConnInstInnerR;
  for (Int_t i = 1; i<5 ; i++) {
    Double_t alpha = TMath::PiOver2()*i/4.;
    xv[3+i] = xv[3] - (1 - TMath::Cos(alpha))*fgkOBCPConnInstInnerR;
    yv[3+i] = yv[3]      - TMath::Sin(alpha) *fgkOBCPConnInstInnerR;
  }

  for (Int_t i = 0; i<8; i++) {
    xv[8+i] = -xv[7-i];
    yv[8+i] =  yv[7-i];
  }

  TGeoXtru *connInsertSh = new TGeoXtru(2);
  connInsertSh->DefinePolygon(16, xv, yv);
  connInsertSh->DefineSection(0,-fgkOBCPConnInsertYHei/2);
  connInsertSh->DefineSection(1, fgkOBCPConnInsertYHei/2);

  TGeoVolume *connInsert = new TGeoVolume("OBColdPlateConnectorInsertSquare",
					  connInsertSh, medWC);
  connInsert->SetFillColor(kGray);
  connInsert->SetLineColor(kGray);


  // The fitting tubes, a Tube
  Double_t rmin = fgkOBCPConnAFitExtD/2 - fgkOBCPConnAFitThick;
  TGeoTube *connFitSh = new TGeoTube(rmin, fgkOBCPConnAFitExtD/2, 
				     fgkOBCPConnAFitZLen/2);

  TGeoVolume *connFit = new TGeoVolume("OBColdPlateConnectorFitting",
				       connFitSh, medInox304);
  connFit->SetFillColor(kGray);
  connFit->SetLineColor(kGray);


  // Now create the container: cannot be a simple box
  // to avoid fake overlaps with stave elements
  xlen = fgkOBCPConnectorXWidth;
  ylen = fgkOBCPConnBlockYHei;
  zlen = fgkOBCPConnBlockZLen + fgkOBCPConnAFitZOut;
  TGeoBBox *connBox = new TGeoBBox("connectorOBCPA", xlen/2, ylen/2, zlen/2);

  ypos = -connBox->GetDY();
  zpos = -connBox->GetDZ();
  TGeoTranslation *transBoxHoll = new TGeoTranslation("transBoxHollA",
						      0, ypos, zpos);
  transBoxHoll->RegisterYourself();

  xpos = fgkOBCPConnTubesXDist/2;
  ypos = -connBox->GetDY() + fgkOBCPConnTubesYPos;
  zpos =  connBox->GetDZ();
  TGeoTranslation *trans1BoxHole = new TGeoTranslation("trans1BoxAHole",
						      -xpos, ypos,-zpos);
  trans1BoxHole->RegisterYourself();
  TGeoTranslation *trans2BoxHole = new TGeoTranslation("trans2BoxAHole",
						       xpos, ypos,-zpos);
  trans2BoxHole->RegisterYourself();


  TGeoCompositeShape *connectSh = new TGeoCompositeShape(
     "connectorOBCPA-connBlockHollA:transBoxHollA-tube1AHole:trans1BoxAHole-tube1AHole:trans2BoxAHole");


  TGeoVolume *connectorASide = new TGeoVolume("OBColdPlateConnectorASide",
					      connectSh, medAir);


  // Finally build up the connector
  zpos = -connectSh->GetDZ() + connBlock->GetDZ();
  connectorASide->AddNode(connBlockA, 1, new TGeoTranslation(0, 0, zpos));

  zpos = -connectSh->GetDZ() + fgkOBCPConnHollowZLen + fgkOBCPConnSqrHoleZPos;
  connectorASide->AddNode(connInsert, 1, new TGeoCombiTrans(0, 0, zpos,
					     new TGeoRotation("",0,-90,0)));

  xpos = fgkOBCPConnTubesXDist/2;
  ypos = -connBlock->GetDY() + fgkOBCPConnTubesYPos;
  zpos = connectSh->GetDZ() - connFitSh->GetDz();
  connectorASide->AddNode(connFit, 1, new TGeoTranslation(-xpos, ypos, zpos));
  connectorASide->AddNode(connFit, 2, new TGeoTranslation( xpos, ypos, zpos));


}

//________________________________________________________________________
void AliITSUv2Layer::CreateOBColdPlateConnectorsCSide(){
//
// Create the C-Side end-stave connectors for IB staves
//
// Input:
//
// Output:
//
// Return:
//
// Created:      29 May 2015  Mario Sitta
//

  // The geoManager
  const TGeoManager *mgr = gGeoManager;

  // Local variables
  const Int_t nv = 16;
  Double_t xv[nv], yv[nv];
  Double_t xlen, ylen, zlen;
  Double_t xpos, ypos, zpos;
  char volname[30];


  // Gather all material pointers
  TGeoMedium *medAir      = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medPEEK     = mgr->GetMedium("ITS_PEEKCF30$");
  TGeoMedium *medWC       = mgr->GetMedium("ITS_TUNGCARB$");


  // First create all elements

  // The connector block, a Composite Shape
  xlen = fgkOBCPConnectorXWidth;
  ylen = fgkOBCPConnBlockYHei;
  zlen = fgkOBCPConnBlockZLen;
  TGeoBBox *connBlock = new TGeoBBox("connBlockC", xlen/2, ylen/2, zlen/2);


  xv[0] =  fgkOBCPConnectorXWidth*0.6;
  yv[0] = -fgkOBCPConnHollowYHei;
  xv[1] = xv[0];
  yv[1] =  fgkOBCPConnHollowYHei;
  xv[2] =  fgkOBCPConnTubesXDist/2 + fgkOBCPConnTubeHole1D/2;
  yv[2] = yv[1];
  xv[3] = xv[2];
  yv[3] =  fgkOBCPConnTubesYPos;
  xv[4] =  fgkOBCPConnTubesXDist/2 - fgkOBCPConnTubeHole1D/2;
  yv[4] = yv[3];
  xv[5] = xv[4];
  yv[5] = yv[2];

  for (Int_t i = 0; i<6; i++) {
    xv[6+i] = -xv[5-i];
    yv[6+i] =  yv[5-i];
  }

  TGeoXtru *connBlockHoll = new TGeoXtru(2);
  connBlockHoll->SetName("connBlockHollC");
  connBlockHoll->DefinePolygon(12, xv, yv);
  connBlockHoll->DefineSection(0,-fgkOBCPConnHollowZLen);
  connBlockHoll->DefineSection(1, fgkOBCPConnHollowZLen);

  ypos = -connBlock->GetDY();
  zpos =  connBlock->GetDZ();
  TGeoTranslation *transBlockHoll = new TGeoTranslation("transBlockHollC",
							0, ypos, zpos);
  transBlockHoll->RegisterYourself();


  TGeoTube *connRoundHole = new TGeoTube("connCRoundHole", 0,
					 fgkOBCPConnRoundHoleD/2,
					 fgkOBCPConnBlockYHei/1.5);

  zpos = connBlock->GetDZ() - fgkOBCPConnHollowZLen - fgkOBCPConnRndHoleZPos;
  TGeoCombiTrans *transRoundHole = new TGeoCombiTrans("transCRoundHole",
						      0, 0, zpos,
						new TGeoRotation("",0,90,0));
  transRoundHole->RegisterYourself();


  zlen = fgkOBCPConnTubeHole1Z;
  TGeoTube *connTubeHole1 = new TGeoTube("tube1CHole", 0,
					 fgkOBCPConnTubeHole1D/2,
					 zlen);

  xpos = fgkOBCPConnTubesXDist/2;
  ypos = -connBlock->GetDY() + fgkOBCPConnTubesYPos;
  zpos =  connBlock->GetDZ();
  TGeoTranslation *trans1Tube1AHole = new TGeoTranslation("trans1Tube1CHole",
							  -xpos, ypos, zpos);
  trans1Tube1AHole->RegisterYourself();
  TGeoTranslation *trans2Tube1AHole = new TGeoTranslation("trans2Tube1CHole",
							   xpos, ypos, zpos);
  trans2Tube1AHole->RegisterYourself();


  TGeoTube *connTubeHole2 = new TGeoTube("tube2CHole", 0,
					 fgkOBCPConnTubeHole2D/2,
					 connBlock->GetDZ());

  zpos = fgkOBCPConnTubeHole3ZP;
  TGeoTranslation *connTubes2Trans1 = new TGeoTranslation("trans1Tube2CHole",
							  -xpos, ypos, zpos);
  connTubes2Trans1->RegisterYourself();
  TGeoTranslation *connTubes2Trans2 = new TGeoTranslation("trans2Tube2CHole",
							   xpos, ypos, zpos);
  connTubes2Trans2->RegisterYourself();


  TGeoTube *connTubeHole3 = new TGeoTube("tube3CHole", 0,
					 fgkOBCPConnTubeHole2D/2,
					 connBlock->GetDX());

  xpos = -fgkOBCPConnTubeHole3XP;
  zpos = -connBlock->GetDZ() + fgkOBCPConnTubeHole3ZP;
  TGeoCombiTrans *connTubes3Trans = new TGeoCombiTrans("transTube3CHole",
						        xpos, ypos, zpos,
					       new TGeoRotation("",90,-90,90));
  connTubes3Trans->RegisterYourself();


  TGeoCompositeShape *connBlockSh = new TGeoCompositeShape(
    "connBlockC-connBlockHollC:transBlockHollC-connCRoundHole:transCRoundHole-tube1CHole:trans1Tube1CHole-tube1CHole:trans2Tube1CHole-tube2CHole:trans1Tube2CHole-tube2CHole:trans2Tube2CHole-tube3CHole:transTube3CHole");


  TGeoVolume *connBlockC = new TGeoVolume("OBColdPlateConnectorBlockCSide",
					  connBlockSh,medPEEK);
  connBlockC->SetFillColor(42);  // Brownish shade
  connBlockC->SetLineColor(42);


  // The steel insert, a Tube
  TGeoTube *connInsertSh = new TGeoTube(fgkOBCPConnInsertD/2,
					fgkOBCPConnRoundHoleD/2,
					fgkOBCPConnInsertYHei/2);

  TGeoVolume *connInsert = new TGeoVolume("OBColdPlateConnectorInsertRound",
					  connInsertSh, medWC);
  connInsert->SetFillColor(kGray);
  connInsert->SetLineColor(kGray);


  // The plug, a Pcon
  TGeoPcon *connPlugSh = new TGeoPcon(0,360,4);
  connPlugSh->DefineSection(0,                   0., 0.,
			       fgkOBCPConnTubeHole2D/2);
  connPlugSh->DefineSection(1, fgkOBCPConnPlugThick, 0.,
			       fgkOBCPConnTubeHole2D/2);
  connPlugSh->DefineSection(2, fgkOBCPConnPlugThick , fgkOBCPConnPlugInnerD/2,
			       fgkOBCPConnTubeHole2D/2);
  connPlugSh->DefineSection(3, fgkOBCPConnPlugTotLen, fgkOBCPConnPlugInnerD/2,
			       fgkOBCPConnTubeHole2D/2);

  TGeoVolume *connPlug = new TGeoVolume("OBCPConnectorPlugC",
					connPlugSh,medPEEK);
  connPlug->SetFillColor(44);  // Brownish shade (a bit darker to spot it)
  connPlug->SetLineColor(44);


  // Now create the container: cannot be a simple box
  // to avoid fake overlaps with stave elements
  xlen = fgkOBCPConnectorXWidth;
  ylen = fgkOBCPConnBlockYHei;
  zlen = fgkOBCPConnBlockZLen;
  TGeoBBox *connBox = new TGeoBBox("connectorOBCPC", xlen/2, ylen/2, zlen/2);

  ypos = -connBox->GetDY();
  zpos =  connBox->GetDZ();
  TGeoTranslation *transBoxHoll = new TGeoTranslation("transBoxHollC",
						      0, ypos, zpos);
  transBoxHoll->RegisterYourself();

  xpos = fgkOBCPConnTubesXDist/2;
  ypos = -connBox->GetDY() + fgkOBCPConnTubesYPos;
  zpos =  connBox->GetDZ();
  TGeoTranslation *trans1BoxHole = new TGeoTranslation("trans1BoxCHole",
						      -xpos, ypos, zpos);
  trans1BoxHole->RegisterYourself();
  TGeoTranslation *trans2BoxHole = new TGeoTranslation("trans2BoxCHole",
						       xpos, ypos, zpos);
  trans2BoxHole->RegisterYourself();


  TGeoCompositeShape *connectSh = new TGeoCompositeShape(
    "connectorOBCPC-connBlockHollC:transBoxHollC-tube1CHole:trans1BoxCHole-tube1CHole:trans2BoxCHole");


  TGeoVolume *connectorCSide = new TGeoVolume("OBColdPlateConnectorCSide",
					      connectSh, medAir);


  // Finally build up the connector
  connectorCSide->AddNode(connBlockC, 1);

  zpos = connBlock->GetDZ() - fgkOBCPConnHollowZLen - fgkOBCPConnRndHoleZPos;
  connectorCSide->AddNode(connInsert, 1, new TGeoCombiTrans(0, 0, zpos,
					     new TGeoRotation("", 0,90, 0)));

  xpos = -connBlock->GetDX();
  ypos = -connBlock->GetDY() + fgkOBCPConnTubesYPos;
  zpos = -connBlock->GetDZ() + fgkOBCPConnTubeHole3ZP;
  connectorCSide->AddNode(connPlug, 1, new TGeoCombiTrans(xpos, ypos, zpos,
					     new TGeoRotation("",90,90,90)));


}

//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateSpaceFrameOuterB(const TGeoManager *mgr){
//
// Create the space frame for the Outer Barrel
//
// Input:
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
//

  TGeoVolume *mechStavVol = 0;

  switch (fStaveModel) {
    case AliITSUv2::kOBModelDummy:
    case AliITSUv2::kOBModel0:
      mechStavVol = CreateSpaceFrameOuterBDummy(mgr);
      break;
    case AliITSUv2::kOBModel1:
    case AliITSUv2::kOBModel2:
      mechStavVol = CreateSpaceFrameOuterB1(mgr);
      break;
    default:
      AliFatal(Form("Unknown stave model %d",fStaveModel));
      break;
  }

  return mechStavVol; 
}

//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateSpaceFrameOuterBDummy(const TGeoManager *) const {
//
// Create dummy stave
//
// Input:
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//


  // Done, return the stave structur
  return 0;
}

//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateSpaceFrameOuterB1(const TGeoManager *mgr){
//
// Create the space frame for the Outer Barrel (Model 1)
// The building blocks are created in another method to avoid
// replicating the same volumes for all OB staves
//
// Input:
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//         a TGeoVolume with the Space Frame of a stave
//
// Created:      03 Feb 2015  Mario Sitta
// Updated:      04 Jun 2015  Mario Sitta  Change container to avoid overlaps
//

  TGeoMedium *medAir = mgr->GetMedium("ITS_AIR$");

  TGeoVolume *unitVol[2], *next2EndVol[2], *endVol[2];
  Double_t *xtru, *ytru;
  Double_t zlen, zpos;
  Int_t nPoints;
  char volname[30];


  // Check whether we have already all pieces
  // Otherwise create them
  unitVol[0] = mgr->GetVolume("SpaceFrameUnit0");

  if (!unitVol[0]) {
    CreateOBSpaceFrameObjects(mgr);
    unitVol[0] = mgr->GetVolume("SpaceFrameUnit0");
  }

  unitVol[1] = mgr->GetVolume("SpaceFrameUnit1");

  next2EndVol[0]  = mgr->GetVolume("SpaceFrameNext2EndUnit0");
  next2EndVol[1]  = mgr->GetVolume("SpaceFrameNext2EndUnit1");

  endVol[0]  = mgr->GetVolume("SpaceFrameEndUnit0");
  endVol[1]  = mgr->GetVolume("SpaceFrameEndUnit1");

  // Get the shape of the units
  // and create a similar shape for the Space Frame container
  TGeoXtru *volShape = (TGeoXtru*)(unitVol[0]->GetShape());

  nPoints = volShape->GetNvert();
  xtru = new Double_t[nPoints];
  ytru = new Double_t[nPoints];

  for (Int_t i=0; i<nPoints; i++) {
    xtru[i] = volShape->GetX(i);
    ytru[i] = volShape->GetY(i);
  }

  Int_t nUnits = fgkOBSpaceFrameNUnits[fLayerNumber/5]; // 3,4 -> 0 - 5,6 -> 1
  zlen = (nUnits-2)*fgkOBSpaceFrameUnitLen; // Take end units out

  TGeoXtru *spaceFrameCentral = new TGeoXtru(2);
  spaceFrameCentral->DefinePolygon(nPoints, xtru, ytru);
  spaceFrameCentral->DefineSection(0,-zlen/2);
  spaceFrameCentral->DefineSection(1, zlen/2);
  snprintf(volname, 30, "sframecentral%d", fLayerNumber);
  spaceFrameCentral->SetName(volname);

  zpos = zlen/2 + fgkOBSpaceFrameUnitLen/2;
  snprintf(volname, 30, "endUnit0Trans%d", fLayerNumber);
  TGeoTranslation *endUnit0Trans = new TGeoTranslation(volname,
						       0, 0,-zpos);
  endUnit0Trans->RegisterYourself();
  snprintf(volname, 30, "endUnit1Trans%d", fLayerNumber);
  TGeoTranslation *endUnit1Trans = new TGeoTranslation(volname,
						       0, 0, zpos);
  endUnit1Trans->RegisterYourself();

  // The Space Frame container: a Composite Shape to avoid overlaps
  // between the U-legs space and the end-stave connectors
  // "endunitcontainer" is defined in CreateOBSpaceFrameObjects)
  char componame[100];
  snprintf(componame, 100,
	   "sframecentral%d+endunitcontainer:endUnit0Trans%d+endunitcontainer:endUnit1Trans%d",
	   fLayerNumber,fLayerNumber,fLayerNumber);

  TGeoCompositeShape *spaceFrame = new TGeoCompositeShape(componame);


  snprintf(volname, 30, "SpaceFrameVolumeLay%d", fLayerNumber);
  TGeoVolume *spaceFrameVol = new TGeoVolume(volname, spaceFrame, medAir);
  spaceFrameVol->SetVisibility(kFALSE);


  // Finally build up the space frame
  TGeoXtru *frameUnit = (TGeoXtru*)(unitVol[0]->GetShape());
  TGeoXtru *endUnit   = (TGeoXtru*)( endVol[0]->GetShape());

  zpos = -spaceFrame->GetDZ() + endUnit->GetDZ();
//  zpos = -fgkOBSpaceFrameZLen[fLayerNumber/5]/2 + endUnit->GetDZ();
  spaceFrameVol->AddNode(endVol[0], 1, new TGeoTranslation(0, 0, zpos));

  zpos += (endUnit->GetDZ() + frameUnit->GetDZ());
  spaceFrameVol->AddNode(next2EndVol[0], 1, new TGeoTranslation(0, 0, zpos));

  for(Int_t i=2; i<nUnits-2; i++){
    zpos = -spaceFrame->GetDZ() + (1 + 2*i)*frameUnit->GetDZ();
//    zpos = -spaceFrame->GetDZ() + 2*i*frameUnit->GetDZ();
    Int_t j = i/2;
    Int_t k = i - j*2;  // alternatively 0 or 1
    spaceFrameVol->AddNode(unitVol[k], j, new TGeoTranslation(0, 0, zpos));
  }

  zpos = -spaceFrame->GetDZ() + (2*nUnits - 3)*frameUnit->GetDZ();
  spaceFrameVol->AddNode(next2EndVol[1], 1, new TGeoTranslation(0, 0, zpos));

  zpos = -spaceFrame->GetDZ() + (2*nUnits - 1)*endUnit->GetDZ();
//  zpos = -fgkOBSpaceFrameZLen[fLayerNumber/5]/2 +
//         (2*nUnits - 1)*endUnit->GetDZ();
  spaceFrameVol->AddNode(endVol[1], 1, new TGeoTranslation(0, 0, zpos));


  // Done, clean up and return the space frame structure
  delete [] xtru;
  delete [] ytru;

  return spaceFrameVol;
}

//________________________________________________________________________
void AliITSUv2Layer::CreateOBSpaceFrameObjects(const TGeoManager *mgr){
//
// Create the space frame building blocks for the Outer Barrel
// This method is practically identical to previous versions of
// CreateSpaceFrameOuterB1
//
// Input:
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//         a TGeoVolume with the Space Frame of a stave
//
// Created:      03 Feb 2015  Mario Sitta
// Updated:      03 Jun 2015  Mario Sitta  End units w/o U-legs
//


  // Materials defined in AliITSUv2
  TGeoMedium *medCarbon       = mgr->GetMedium("ITS_M55J6K$");
  TGeoMedium *medF6151B05M    = mgr->GetMedium("ITS_F6151B05M$");
  TGeoMedium *medAir          = mgr->GetMedium("ITS_AIR$");


  // Local parameters
  Double_t halfFrameWidth  = fgkOBSpaceFrameWidth/2;
  Double_t triangleHeight  = fgkOBSpaceFrameHeight;
  Double_t sframeHeight    = triangleHeight + fgkOBSFrameBaseRibDiam
                                            + fgkOBSFrameULegHeight2*2;
  Double_t staveLa         = fgkOBSpaceFrameTopVL;
  Double_t staveHa         = fgkOBSpaceFrameTopVH;
  Double_t staveLb         = fgkOBSpaceFrameSideVL;
  Double_t staveHb         = fgkOBSpaceFrameSideVH;
  Double_t alphaDeg        = fgkOBSpaceFrameVAlpha;
  Double_t alphaRad        = alphaDeg*TMath::DegToRad()/2;
  Double_t beta            = fgkOBSpaceFrameVBeta*TMath::DegToRad()/2;
  Double_t sideRibRadius   = fgkOBSFrameSideRibDiam/2;
  Double_t sidePhiDeg      = fgkOBSFrameSideRibPhi;
  Double_t sidePhiRad      = sidePhiDeg*TMath::DegToRad();
  Double_t baseRibRadius   = fgkOBSFrameBaseRibDiam/2;
  Double_t basePhiDeg      = fgkOBSFrameBaseRibPhi;
  Double_t basePhiRad      = basePhiDeg*TMath::DegToRad();
  Double_t ulegHalfLen     = fgkOBSFrameULegLen/2;
  Double_t ulegHalfWidth   = fgkOBSFrameULegWidth/2;
  Double_t ulegHigh1       = fgkOBSFrameULegHeight1;
  Double_t ulegHigh2       = fgkOBSFrameULegHeight2;
  Double_t ulegThick       = fgkOBSFrameULegThick;

  Double_t xlen, zlen;
  Double_t xpos, ypos, zpos;
  Double_t unitlen;
  Double_t xtru[22], ytru[22];
  char volname[30];


  unitlen = fgkOBSpaceFrameUnitLen;

  xlen = halfFrameWidth + sideRibRadius;

  // We need a properly shaped Xtru to accomodate the ribs avoiding
  // overlaps with the HalfStave cooling tubes
  xtru[ 0] = fgkOBSFrameULegXPos - ulegHalfLen;
  ytru[ 0] = -(triangleHeight/2 + baseRibRadius);
  xtru[ 1] = xtru[0];
  ytru[ 1] = ytru[0] - ulegHigh1;
  xtru[ 2] = xtru[1] + ulegThick;
  ytru[ 2] = ytru[1];
  xtru[ 3] = xtru[2];
  ytru[ 3] = ytru[0] - ulegThick;
  xtru[ 7] = fgkOBSFrameULegXPos + ulegHalfLen;
  ytru[ 7] = ytru[0];
  xtru[ 6] = xtru[7];
  ytru[ 6] = ytru[1];
  xtru[ 5] = xtru[6] - ulegThick;
  ytru[ 5] = ytru[6];
  xtru[ 4] = xtru[5];
  ytru[ 4] = ytru[3];
  xtru[ 8] = xlen;
  ytru[ 8] = ytru[7];
  xtru[ 9] = xtru[8];
  ytru[ 9] = 0.9*ytru[8];
  xtru[10] = 0.3*xtru[8];
  ytru[10] = triangleHeight/2;
  for (Int_t i=0; i<11; i++) { // Reflect on the X negative side
    xtru[i+11] = -xtru[10-i];
    ytru[i+11] =  ytru[10-i];
  }
  ytru[15] = ytru[0] - ulegHigh2; // U-legs on negative X are longer
  ytru[16] = ytru[15];
  ytru[19] = ytru[15];
  ytru[20] = ytru[15];


  // The space frame single units
  // We need two units because the base ribs are alternately oriented
  // The next-to-end units are slightly different
  TGeoXtru *frameUnit = new TGeoXtru(2);
  frameUnit->DefinePolygon(22, xtru, ytru);
  frameUnit->DefineSection(0,-unitlen/2);
  frameUnit->DefineSection(1, unitlen/2);

  TGeoXtru *next2EndUnit = new TGeoXtru(2);
  next2EndUnit->DefinePolygon(22, xtru, ytru);
  next2EndUnit->DefineSection(0,-unitlen/2);
  next2EndUnit->DefineSection(1, unitlen/2);


  // The end units have no U-legs, so a simpler Xtru
  // to avoid overlaps with the end-stave connectors
  xtru[0] = xlen;
  ytru[0] = -(triangleHeight/2 + baseRibRadius);
  xtru[1] = xtru[0];
  ytru[1] = 0.9*ytru[0];
  xtru[2] = 0.3*xtru[0];
  ytru[2] = triangleHeight/2;
  for (Int_t i=0; i<3; i++) { // Reflect on the X negative side
    xtru[i+3] = -xtru[2-i];
    ytru[i+3] =  ytru[2-i];
  }

  TGeoXtru *endUnit = new TGeoXtru(2);
  endUnit->SetName("endunitcontainer");  // Will be used when create spaceframe
  endUnit->DefinePolygon(6, xtru, ytru);
  endUnit->DefineSection(0,-unitlen/2);
  endUnit->DefineSection(1, unitlen/2);


  // The air containers
  TGeoVolume *unitVol[2];
  unitVol[0] = new TGeoVolume("SpaceFrameUnit0", frameUnit, medAir);
  unitVol[1] = new TGeoVolume("SpaceFrameUnit1", frameUnit, medAir);

  TGeoVolume *next2EndVol[2];
  next2EndVol[0]  = new TGeoVolume("SpaceFrameNext2EndUnit0",
				   next2EndUnit, medAir);
  next2EndVol[1]  = new TGeoVolume("SpaceFrameNext2EndUnit1",
				   next2EndUnit, medAir);

  TGeoVolume *endVol[2];
  endVol[0]  = new TGeoVolume("SpaceFrameEndUnit0", endUnit, medAir);
  endVol[1]  = new TGeoVolume("SpaceFrameEndUnit1", endUnit, medAir);

  // The actual volumes

  //--- The top V of the Carbon Fiber Stave (segment)
  TGeoXtru *cfStavTop = CreateStaveSide("CFstavTopCornerVolshape",
			unitlen/2., alphaRad, beta, staveLa, staveHa, kTRUE);

  TGeoVolume *cfStavTopVol = new TGeoVolume("CFstavTopCornerVol",
					    cfStavTop, medCarbon);
  cfStavTopVol->SetLineColor(35);

  unitVol[0]->AddNode(cfStavTopVol, 1,
		      new TGeoTranslation(0, triangleHeight/2, 0));

  unitVol[1]->AddNode(cfStavTopVol, 1,
		      new TGeoTranslation(0, triangleHeight/2, 0));
  
  next2EndVol[0]->AddNode(cfStavTopVol, 1,
		      new TGeoTranslation(0, triangleHeight/2, 0));

  next2EndVol[1]->AddNode(cfStavTopVol, 1,
		      new TGeoTranslation(0, triangleHeight/2, 0));
  
  endVol[0]->AddNode(cfStavTopVol, 1,
		      new TGeoTranslation(0, triangleHeight/2, 0));

  endVol[1]->AddNode(cfStavTopVol, 1,
		      new TGeoTranslation(0, triangleHeight/2, 0));
  
  //--- The two side V's
  TGeoXtru *cfStavSide = CreateStaveSide("CFstavSideCornerVolshape",
			 unitlen/2., alphaRad, beta, staveLb, staveHb, kFALSE);

  TGeoVolume *cfStavSideVol = new TGeoVolume("CFstavSideCornerVol",
					     cfStavSide, medCarbon);
  cfStavSideVol->SetLineColor(35);

  unitVol[0]->AddNode(cfStavSideVol, 1,
		   new TGeoTranslation( halfFrameWidth, -triangleHeight/2, 0));
  unitVol[0]->AddNode(cfStavSideVol, 2,
		   new TGeoCombiTrans( -halfFrameWidth, -triangleHeight/2, 0,
					 new TGeoRotation("",90,180,-90)));

  unitVol[1]->AddNode(cfStavSideVol, 1,
		   new TGeoTranslation( halfFrameWidth, -triangleHeight/2, 0));
  unitVol[1]->AddNode(cfStavSideVol, 2,
		   new TGeoCombiTrans( -halfFrameWidth, -triangleHeight/2, 0,
					 new TGeoRotation("",90,180,-90)));

  next2EndVol[0]->AddNode(cfStavSideVol, 1,
		   new TGeoTranslation( halfFrameWidth, -triangleHeight/2, 0));
  next2EndVol[0]->AddNode(cfStavSideVol, 2,
		   new TGeoCombiTrans( -halfFrameWidth, -triangleHeight/2, 0,
					 new TGeoRotation("",90,180,-90)));

  next2EndVol[1]->AddNode(cfStavSideVol, 1,
		   new TGeoTranslation( halfFrameWidth, -triangleHeight/2, 0));
  next2EndVol[1]->AddNode(cfStavSideVol, 2,
		   new TGeoCombiTrans( -halfFrameWidth, -triangleHeight/2, 0,
					 new TGeoRotation("",90,180,-90)));

  endVol[0]->AddNode(cfStavSideVol, 1,
		   new TGeoTranslation( halfFrameWidth, -triangleHeight/2, 0));
  endVol[0]->AddNode(cfStavSideVol, 2,
		   new TGeoCombiTrans( -halfFrameWidth, -triangleHeight/2, 0,
					 new TGeoRotation("",90,180,-90)));

  endVol[1]->AddNode(cfStavSideVol, 1,
		   new TGeoTranslation( halfFrameWidth, -triangleHeight/2, 0));
  endVol[1]->AddNode(cfStavSideVol, 2,
		   new TGeoCombiTrans( -halfFrameWidth, -triangleHeight/2, 0,
					 new TGeoRotation("",90,180,-90)));


  //--- The beams
  // Ribs on the sides
  Double_t ribZProj = triangleHeight/TMath::Tan(sidePhiRad);
  Double_t sideRibLen = TMath::Sqrt( ribZProj*ribZProj             +
				     triangleHeight*triangleHeight +
				     halfFrameWidth*halfFrameWidth );

  TGeoTubeSeg *sideRib = new TGeoTubeSeg(0, sideRibRadius,
					 sideRibLen/2, 0, 180);
  TGeoVolume *sideRibVol = new TGeoVolume("CFstavSideBeamVol",
					  sideRib, medCarbon);
  sideRibVol->SetLineColor(35);

  TGeoCombiTrans *sideTransf[4];
  xpos = halfFrameWidth/2 + 0.8*staveHa*TMath::Cos(alphaRad/2);
  ypos = -sideRibRadius/2;
  zpos = unitlen/4;

  sideTransf[0] = new TGeoCombiTrans( xpos, ypos,-zpos,
				      new TGeoRotation("", 90-alphaDeg,
						       -sidePhiDeg, -90));
  sideTransf[1] = new TGeoCombiTrans( xpos, ypos, zpos,
				      new TGeoRotation("", 90-alphaDeg,
						        sidePhiDeg, -90));
  sideTransf[2] = new TGeoCombiTrans(-xpos, ypos,-zpos,
				      new TGeoRotation("", 90+alphaDeg,
						        sidePhiDeg, -90));
  sideTransf[3] = new TGeoCombiTrans(-xpos, ypos, zpos,
				      new TGeoRotation("", 90+alphaDeg,
						       -sidePhiDeg, -90));

  unitVol[0]->AddNode(sideRibVol, 1, sideTransf[0]);
  unitVol[0]->AddNode(sideRibVol, 2, sideTransf[1]);
  unitVol[0]->AddNode(sideRibVol, 3, sideTransf[2]);
  unitVol[0]->AddNode(sideRibVol, 4, sideTransf[3]);

  unitVol[1]->AddNode(sideRibVol, 1, sideTransf[0]);
  unitVol[1]->AddNode(sideRibVol, 2, sideTransf[1]);
  unitVol[1]->AddNode(sideRibVol, 3, sideTransf[2]);
  unitVol[1]->AddNode(sideRibVol, 4, sideTransf[3]);

  next2EndVol[0]->AddNode(sideRibVol, 1, sideTransf[0]);
  next2EndVol[0]->AddNode(sideRibVol, 2, sideTransf[1]);
  next2EndVol[0]->AddNode(sideRibVol, 3, sideTransf[2]);
  next2EndVol[0]->AddNode(sideRibVol, 4, sideTransf[3]);

  next2EndVol[1]->AddNode(sideRibVol, 1, sideTransf[0]);
  next2EndVol[1]->AddNode(sideRibVol, 2, sideTransf[1]);
  next2EndVol[1]->AddNode(sideRibVol, 3, sideTransf[2]);
  next2EndVol[1]->AddNode(sideRibVol, 4, sideTransf[3]);

  endVol[0]->AddNode(sideRibVol, 1, sideTransf[0]);
  endVol[0]->AddNode(sideRibVol, 2, sideTransf[1]);
  endVol[0]->AddNode(sideRibVol, 3, sideTransf[2]);
  endVol[0]->AddNode(sideRibVol, 4, sideTransf[3]);

  endVol[1]->AddNode(sideRibVol, 1, sideTransf[0]);
  endVol[1]->AddNode(sideRibVol, 2, sideTransf[1]);
  endVol[1]->AddNode(sideRibVol, 3, sideTransf[2]);
  endVol[1]->AddNode(sideRibVol, 4, sideTransf[3]);


  // Ribs on the bottom
  // Rib1 are the inclined ones, Rib2 the straight ones
  Double_t baseRibLen = 0.98*2*halfFrameWidth/TMath::Sin(basePhiRad);

  TGeoTubeSeg *baseRib1 = new TGeoTubeSeg(0, baseRibRadius,
					  baseRibLen/2, 0, 180);
  TGeoVolume *baseRib1Vol = new TGeoVolume("CFstavBaseBeam1Vol",
					   baseRib1, medCarbon);
  baseRib1Vol->SetLineColor(35);

  TGeoTubeSeg *baseRib2 = new TGeoTubeSeg(0, baseRibRadius,
					  halfFrameWidth, 0, 90);
  TGeoVolume *baseRib2Vol = new TGeoVolume("CFstavBaseBeam2Vol",
					   baseRib2, medCarbon);
  baseRib2Vol->SetLineColor(35);

  TGeoTubeSeg *baseEndRib = new TGeoTubeSeg(0, baseRibRadius,
					    halfFrameWidth, 0, 180);
  TGeoVolume *baseEndRibVol = new TGeoVolume("CFstavBaseEndBeamVol",
					     baseEndRib, medCarbon);
  baseEndRibVol->SetLineColor(35);

  TGeoCombiTrans *baseTransf[6];
  ypos = triangleHeight/2;
  zpos = unitlen/2;

  baseTransf[0] = new TGeoCombiTrans("", 0, -ypos, -zpos,
				     new TGeoRotation("", 90, 90,  90));
  baseTransf[1] = new TGeoCombiTrans("", 0, -ypos,  zpos,
				     new TGeoRotation("",-90, 90, -90));
  baseTransf[2] = new TGeoCombiTrans(0, -ypos, 0,
				     new TGeoRotation("",-90, basePhiDeg,-90));
  baseTransf[3] = new TGeoCombiTrans(0, -ypos, 0,
				     new TGeoRotation("",-90,-basePhiDeg,-90));
  zpos -= baseEndRib->GetRmax();
  baseTransf[4] = new TGeoCombiTrans("", 0, -ypos, -zpos,
				     new TGeoRotation("", 90, 90,  90));
  baseTransf[5] = new TGeoCombiTrans("", 0, -ypos,  zpos,
				     new TGeoRotation("", 90, 90,  90));

  unitVol[0]->AddNode(baseRib2Vol, 1, baseTransf[0]);
  unitVol[0]->AddNode(baseRib2Vol, 2, baseTransf[1]);
  unitVol[0]->AddNode(baseRib1Vol, 1, baseTransf[2]);

  unitVol[1]->AddNode(baseRib2Vol, 1, baseTransf[0]);
  unitVol[1]->AddNode(baseRib2Vol, 2, baseTransf[1]);
  unitVol[1]->AddNode(baseRib1Vol, 1, baseTransf[3]);

  next2EndVol[0]->AddNode(baseRib2Vol, 1, baseTransf[0]);
  next2EndVol[0]->AddNode(baseRib2Vol, 2, baseTransf[1]);
  next2EndVol[0]->AddNode(baseRib1Vol, 1, baseTransf[3]);

  next2EndVol[1]->AddNode(baseRib2Vol, 1, baseTransf[0]);
  next2EndVol[1]->AddNode(baseRib2Vol, 2, baseTransf[1]);
  next2EndVol[1]->AddNode(baseRib1Vol, 1, baseTransf[3]);

  endVol[0]->AddNode(baseEndRibVol, 1, baseTransf[4]);
  endVol[0]->AddNode(baseRib2Vol,   1, baseTransf[1]);
  endVol[0]->AddNode(baseRib1Vol,   1, baseTransf[2]);

  endVol[1]->AddNode(baseEndRibVol, 1, baseTransf[5]);
  endVol[1]->AddNode(baseRib2Vol,   1, baseTransf[0]);
  endVol[1]->AddNode(baseRib1Vol,   1, baseTransf[2]);


  // U-Legs
  // The shorter
  xtru[0] = ulegHalfLen;
  ytru[0] = 0;
  xtru[1] = xtru[0];
  ytru[1] = -ulegHigh1;
  xtru[2] = xtru[1] - ulegThick;
  ytru[2] = ytru[1];
  xtru[3] = xtru[2];
  ytru[3] = ytru[0] - ulegThick;
  for (Int_t i=0; i<4; i++) { // Reflect on the X negative side
    xtru[i+4] = -xtru[3-i];
    ytru[i+4] =  ytru[3-i];
  }

  TGeoXtru *uleg1full = new TGeoXtru(2);  // This will go in the next end units
  uleg1full->DefinePolygon(8, xtru, ytru);
  uleg1full->DefineSection(0,-ulegHalfWidth);
  uleg1full->DefineSection(1, ulegHalfWidth);

  TGeoXtru *uleg1half = new TGeoXtru(2);  // This will go in the middle unitys
  uleg1half->DefinePolygon(8, xtru, ytru);
  uleg1half->DefineSection(0,-ulegHalfWidth/2);
  uleg1half->DefineSection(1, ulegHalfWidth/2);

  TGeoVolume *uleg1fullVol = new TGeoVolume("CFstavULeg1FullVol",
					    uleg1full, medF6151B05M);
  uleg1fullVol->SetLineColor(35);

  TGeoVolume *uleg1halfVol = new TGeoVolume("CFstavULeg1HalfVol",
					    uleg1half, medF6151B05M);
  uleg1halfVol->SetLineColor(35);

  // The longer
  ytru[1] = -ulegHigh2;
  ytru[2] = -ulegHigh2;
  ytru[5] = -ulegHigh2;
  ytru[6] = -ulegHigh2;

  TGeoXtru *uleg2full = new TGeoXtru(2);  // This will go in the next end units
  uleg2full->DefinePolygon(8, xtru, ytru);
  uleg2full->DefineSection(0,-ulegHalfWidth);
  uleg2full->DefineSection(1, ulegHalfWidth);

  TGeoXtru *uleg2half = new TGeoXtru(2);  // This will go in the middle unitys
  uleg2half->DefinePolygon(8, xtru, ytru);
  uleg2half->DefineSection(0,-ulegHalfWidth/2);
  uleg2half->DefineSection(1, ulegHalfWidth/2);

  TGeoVolume *uleg2fullVol = new TGeoVolume("CFstavULeg2FullVol",
					    uleg2full, medF6151B05M);
  uleg2fullVol->SetLineColor(35);

  TGeoVolume *uleg2halfVol = new TGeoVolume("CFstavULeg2HalfVol",
					    uleg2half, medF6151B05M);
  uleg2halfVol->SetLineColor(35);


  xpos = fgkOBSFrameULegXPos;
  ypos = triangleHeight/2 + baseRibRadius;
  zpos = unitlen/2 - uleg1half->GetZ(1);

  unitVol[0]->AddNode(uleg1halfVol, 1,  // Shorter on +X
		      new TGeoTranslation(  xpos, -ypos, -zpos));
  unitVol[0]->AddNode(uleg1halfVol, 2,
		      new TGeoTranslation(  xpos, -ypos,  zpos));

  unitVol[1]->AddNode(uleg1halfVol, 1,
		      new TGeoTranslation(  xpos, -ypos, -zpos));
  unitVol[1]->AddNode(uleg1halfVol, 2,
		      new TGeoTranslation(  xpos, -ypos,  zpos));

  unitVol[0]->AddNode(uleg2halfVol, 1,  // Longer on -X
		      new TGeoTranslation( -xpos, -ypos, -zpos));
  unitVol[0]->AddNode(uleg2halfVol, 2,
		      new TGeoTranslation( -xpos, -ypos,  zpos));

  unitVol[1]->AddNode(uleg2halfVol, 1,
		      new TGeoTranslation( -xpos, -ypos, -zpos));
  unitVol[1]->AddNode(uleg2halfVol, 2,
		      new TGeoTranslation( -xpos, -ypos,  zpos));

  next2EndVol[0]->AddNode(uleg1halfVol, 1,
		  new TGeoTranslation(  xpos, -ypos,  zpos));
  next2EndVol[0]->AddNode(uleg2halfVol, 1,
		  new TGeoTranslation( -xpos, -ypos,  zpos));

  next2EndVol[1]->AddNode(uleg1halfVol, 1,
		  new TGeoTranslation(  xpos, -ypos, -zpos));
  next2EndVol[1]->AddNode(uleg2halfVol, 1,
		  new TGeoTranslation( -xpos, -ypos, -zpos));

  zpos = unitlen/2 - uleg1full->GetZ(1);
  next2EndVol[0]->AddNode(uleg1fullVol, 1,
		  new TGeoTranslation(  xpos, -ypos, -zpos));
  next2EndVol[0]->AddNode(uleg2fullVol, 1,
 		  new TGeoTranslation( -xpos, -ypos, -zpos));

  next2EndVol[1]->AddNode(uleg1fullVol, 1,
		  new TGeoTranslation(  xpos, -ypos,  zpos));
  next2EndVol[1]->AddNode(uleg2fullVol, 1,
		  new TGeoTranslation( -xpos, -ypos,  zpos));


  // Done
  return;
}

//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateChipInnerB(const Double_t xchip,
					     const Double_t ychip,   
					     const Double_t zchip,
					     const TGeoManager *mgr){
//
// Creates the actual Chip
//
// Input:
//         xchip,ychip,zchip : the chip dimensions
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      22 Jun 2011  Mario Sitta
//

  char volname[30];
  Double_t xlen, ylen, zlen;
  Double_t xpos, ypos, zpos;


  // First create all needed shapes

  // The chip
  TGeoBBox *chip = new TGeoBBox(xchip,  ychip, zchip);

  // The sensor
  xlen = chip->GetDX();
  ylen = 0.5*fSensorThick;
  zlen = chip->GetDZ();
  TGeoBBox *sensor = new TGeoBBox(xlen, ylen, zlen);


  // We have all shapes: now create the real volumes
  TGeoMedium *medSi  = mgr->GetMedium("ITS_SI$");
  TGeoMedium *medAir = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medChip;

  if ( (fLayerNumber <  fgkNumberOfInnerLayers & fBuildLevel < 6) ||
       (fLayerNumber >= fgkNumberOfInnerLayers & fBuildLevel < 7) )
    medChip = medSi;
  else
    medChip = medAir;

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSChipPattern(), fLayerNumber);
  TGeoVolume *chipVol = new TGeoVolume(volname, chip, medChip);
  chipVol->SetVisibility(kTRUE);
  chipVol->SetLineColor(1);

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSSensorPattern(), fLayerNumber);
  TGeoVolume *sensVol = new TGeoVolume(volname, sensor, medChip);
  sensVol->SetVisibility(kTRUE);
  sensVol->SetLineColor(8);
  sensVol->SetLineWidth(1);
  sensVol->SetFillColor(sensVol->GetLineColor());
  sensVol->SetFillStyle(4000); // 0% transparent


  // Now build up the chip
  xpos = 0.;
  ypos = -chip->GetDY() + sensor->GetDY();
  zpos = 0.;

  chipVol->AddNode(sensVol, 1, new TGeoTranslation(xpos, ypos, zpos));

  // Done, return the chip
  return chipVol;
}

//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateModuleOuterB(const TGeoManager *mgr){
//
// Creates the OB Module: HIC + FPC + Carbon plate
//
// Input:
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//         the module as a TGeoVolume
//
// Created:      18 Dec 2013  M. Sitta, A. Barbano
// Updated:      26 Feb 2014  M. Sitta
// Updated:      12 Nov 2014  M. Sitta  Model2 is w/o Carbon Plate and Glue
//                                      and Cu instead of Al
//


  char volname[30];

  Double_t xGap  = fgkOBChipXGap;
  Double_t zGap  = fgkOBChipZGap;

  Double_t xchip, ychip, zchip;
  Double_t xlen, ylen, zlen;
  Double_t xpos, ypos, zpos;
  
  // First create all needed shapes

  // The chip (the same as for IB)
  xlen = (fgkOBHalfStaveWidth/2-xGap/2)/fgkOBNChipRows;
  ylen = 0.5*fChipThick;
  zlen = (fgkOBModuleZLength - (fgkOBChipsPerRow-1)*zGap)/(2*fgkOBChipsPerRow);

  TGeoVolume *chipVol = CreateChipInnerB(xlen, ylen, zlen);

  xchip = ((TGeoBBox*)chipVol->GetShape())->GetDX();
  ychip = ((TGeoBBox*)chipVol->GetShape())->GetDY();
  zchip = ((TGeoBBox*)chipVol->GetShape())->GetDZ();

  // The module carbon plate
  xlen = fgkOBHalfStaveWidth/2;
  ylen = fgkOBCarbonPlateThick/2;
  zlen = fgkOBModuleZLength/2;
  TGeoBBox *modPlate = new TGeoBBox("CarbonPlate", xlen, ylen, zlen);

  // The glue
  ylen = fgkOBGlueThickM1/2;
  TGeoBBox *glue = new TGeoBBox("Glue", xlen, ylen, zlen);

  // The FPC cables
  TGeoBBox *flexMetal;
  if (fStaveModel == AliITSUv2::kOBModel1) {
    ylen = fgkOBFlexCableAlThick/2;
    flexMetal = new TGeoBBox("FlexAl", xlen, ylen, zlen);
  } else {
    ylen = fgkOBFlexCableCuThick/2;
    flexMetal = new TGeoBBox("FlexCu", xlen, ylen, zlen);
  }

  if (fStaveModel == AliITSUv2::kOBModel1)
    ylen = fgkOBFlexCableKapThick1/2;
  else
    ylen = fgkOBFlexCableKapThick/2;
  TGeoBBox *flexKap = new TGeoBBox("FlexKap", xlen, ylen, zlen);

  // The module
  xlen = fgkOBHalfStaveWidth/2;
  ylen = ychip + flexMetal->GetDY() + flexKap->GetDY();
  if (fStaveModel == AliITSUv2::kOBModel1)
    ylen += (modPlate->GetDY() + glue->GetDY());
  zlen = fgkOBModuleZLength/2;
  TGeoBBox *module = new TGeoBBox("OBModule", xlen,  ylen, zlen);


  // We have all shapes: now create the real volumes
 
  TGeoMedium *medAir      = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medCarbon   = mgr->GetMedium("ITS_CARBON$");
  TGeoMedium *medGlue     = mgr->GetMedium("ITS_GLUE$");
  TGeoMedium *medAluminum = mgr->GetMedium("ITS_ALUMINUM$");
  TGeoMedium *medCopper   = mgr->GetMedium("ITS_COPPER$");
  TGeoMedium *medKapton   = mgr->GetMedium("ITS_KAPTON(POLYCH2)$");

  TGeoVolume *modPlateVol = new TGeoVolume("CarbonPlateVol",
					    modPlate, medCarbon);
  modPlateVol->SetLineColor(kMagenta-8);
  modPlateVol->SetFillColor(modPlateVol->GetLineColor());
  modPlateVol->SetFillStyle(4000); // 0% transparent

  TGeoVolume *glueVol = new TGeoVolume("GlueVol", glue, medGlue);
  glueVol->SetLineColor(kBlack);
  glueVol->SetFillColor(glueVol->GetLineColor());
  glueVol->SetFillStyle(4000); // 0% transparent

  TGeoVolume *flexMetalVol;
  if (fStaveModel == AliITSUv2::kOBModel1)
    flexMetalVol = new TGeoVolume("FPCAlVol", flexMetal, medAluminum);
  else
    flexMetalVol = new TGeoVolume("FPCCuVol", flexMetal, medCopper);
  flexMetalVol->SetLineColor(kRed);
  flexMetalVol->SetFillColor(flexMetalVol->GetLineColor());
  flexMetalVol->SetFillStyle(4000); // 0% transparent

  TGeoVolume *flexKapVol = new TGeoVolume("FPCKapVol", flexKap, medKapton);
  flexKapVol->SetLineColor(kGreen);
  flexKapVol->SetFillColor(flexKapVol->GetLineColor());
  flexKapVol->SetFillStyle(4000); // 0% transparent

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSModulePattern(), fLayerNumber);
  TGeoVolume *modVol = new TGeoVolume(volname, module, medAir);
  modVol->SetVisibility(kTRUE);
  

  // Now build up the module
  // Model1 : CarbonPlate-Glue-Chips-AluminumFPC-KaptonFPC
  // Model2 : Chips-CopperFPC-KaptonFPCChips
  ypos = -module->GetDY();

  if (fStaveModel == AliITSUv2::kOBModel1) {
    ypos += modPlate->GetDY();
    if (fBuildLevel < 5) // Carbon
      modVol->AddNode(modPlateVol, 1, new TGeoTranslation(0, ypos, 0));

    ypos += (modPlate->GetDY() + glue->GetDY());
    if (fBuildLevel < 2) // Glue
      modVol->AddNode(glueVol, 1, new TGeoTranslation(0, ypos, 0));

    ypos += glue->GetDY();
  }

  xpos = -module->GetDX() + xchip;
  ypos += ychip;
  for(Int_t k=0; k<fgkOBChipsPerRow; k++)   //put 7x2 chip into one module
    {
      zpos = -module->GetDZ() + zchip + k*(2*zchip + zGap);
      modVol->AddNode(chipVol, 2*k  , new TGeoTranslation( xpos, ypos, zpos));
      modVol->AddNode(chipVol, 2*k+1, 
      new TGeoCombiTrans(-xpos, ypos, zpos, new TGeoRotation("",0,180,180)));
      fHierarchy[kChip]+=2;
    }

  ypos += (ychip + flexMetal->GetDY());
  if (fBuildLevel < 1) // Model1: Aluminum  Model2: Copper
    modVol->AddNode(flexMetalVol, 1, new TGeoTranslation(0, ypos, 0));

  ypos += (flexMetal->GetDY() + flexKap->GetDY());
  if ( (fStaveModel == AliITSUv2::kOBModel1 && fBuildLevel < 4) || // Kapton
       (fStaveModel == AliITSUv2::kOBModel2 && fBuildLevel < 5) )
    modVol->AddNode(flexKapVol, 1, new TGeoTranslation(0, ypos, 0));

  // Done, return the module
  return modVol;
}

//________________________________________________________________________
Double_t AliITSUv2Layer::GetGammaConversionRodDiam()
{
//
// Gets the diameter of the gamma conversion rods, if defined
//
//
// Input:
//
// Output:
//
// Return:
//         the diameter of the gamma conversion rods for this layer
//
// Created:      26 Oct 2016  Mario Sitta
//

  if (fAddGammaConv)
    return fGammaConvDiam;
  else
    AliWarning("Gamma Conversion rods not defined for this layer");  
}

//________________________________________________________________________
Double_t AliITSUv2Layer::GetGammaConversionRodXPos()
{
//
// Gets the X position of the gamma conversion rods, if defined
//
//
// Input:
//
// Output:
//
// Return:
//         the X position of the gamma conversion rods for this layer
//         in the Half Stave reference system
//
// Created:      26 Oct 2016  Mario Sitta
//

  if (fAddGammaConv)
    return fGammaConvXPos;
  else
    AliWarning("Gamma Conversion rods not defined for this layer");  
}

//________________________________________________________________________
void AliITSUv2Layer::SetNUnits(Int_t u)
{
//
// Sets the number of units in a stave:
//      for the Inner Barrel: the number of chips per stave
//      for the Outer Barrel: the number of modules per half stave
//
//
// Input:
//         u :  the number of units
//
// Output:
//
// Return:
//
// Created:      18 Feb 2013  Mario Sitta (was already SetNChips)
//

  if (fLayerNumber < fgkNumberOfInnerLayers)
    fNChips = u;
  else {
    fNModules = u;
    fNChips = fgkOBChipsPerRow;
  }

}

//________________________________________________________________________
void AliITSUv2Layer::SetStaveTilt(const Double_t t)
{
//
// Sets the Stave tilt angle (for turbo layers only)
//
// Input:
//         t :  the stave tilt angle
//
// Output:
//
// Return:
//
// Created:      08 Jul 2011  Mario Sitta
//

  if (fIsTurbo)
    fStaveTilt = t;
  else
    AliError("Not a Turbo layer");

}

//________________________________________________________________________
void AliITSUv2Layer::SetStaveWidth(const Double_t w){
//
// Sets the Stave width (for turbo layers only)
//
// Input:
//         w :  the stave width
//
// Output:
//
// Return:
//
// Created:      08 Jul 2011  Mario Sitta
//

  if (fIsTurbo)
    fStaveWidth = w;
  else
    AliError("Not a Turbo layer");

}

//________________________________________________________________________
TGeoXtru *AliITSUv2Layer::CreateStaveSide(const char *name,
                         Double_t dz, Double_t alpha, Double_t beta,
                         Double_t L, Double_t H, Bool_t top) {
//
// Creates the V-shaped sides of the OB space frame
// (from a similar method with same name and function
// in AliITSv11GeometrySDD class by L.Gaudichet)
//
// Updated:      15 Dec 2014  Mario Sitta  Rewritten using Xtru
// Updated:      09 Jan 2015  Mario Sitta  Rewritten again using different
//                                         aperture angles (info by C.Gargiulo)
//

    // Create the V shape corner of CF stave

    const Int_t nv = 6;
    Double_t xv[nv], yv[nv];
  
    TGeoXtru *cfStavSide = new TGeoXtru(2);
    cfStavSide->SetName(name);

    Double_t theta = TMath::PiOver2() - beta;
    Double_t gamma = beta - alpha;
    // Points must be in clockwise order
    if (top) { // TOP - vertices not in order
      xv[3] = 0;
      yv[3] = 0;
      xv[2] =  L*TMath::Sin(alpha);
      yv[2] = -L*TMath::Cos(alpha);
      xv[1] = xv[2] - H*TMath::Cos(alpha);
      yv[1] = yv[2] - H*TMath::Sin(alpha);
      xv[0] = 0;
      yv[0] = yv[1] + TMath::Tan(theta)*xv[1];
      xv[4] = -xv[2];  // Reflect
      yv[4] =  yv[2];
      xv[5] = -xv[1];
      yv[5] =  yv[1];
    } else { // SIDE
      Double_t m = -TMath::Tan(alpha), n = TMath::Tan(gamma);
      xv[0] = 0;
      yv[0] = 0;
      xv[1] = -L*TMath::Cos(2*alpha);
      yv[1] =  L*TMath::Sin(2*alpha);
      xv[2] = xv[1] - H*TMath::Sin(2*alpha);
      yv[2] = yv[1] - H*TMath::Cos(2*alpha);
      xv[4] = -L;
      yv[4] =  H;
      xv[5] = xv[4];
      yv[5] = 0;
      xv[3] = (yv[4] - n*xv[4])/(m - n);
      yv[3] = m*xv[3];
    }

    cfStavSide->DefinePolygon(nv, xv, yv);
    cfStavSide->DefineSection(0,-dz);
    cfStavSide->DefineSection(1, dz);

    return cfStavSide;
}

//________________________________________________________________________
TGeoCombiTrans *AliITSUv2Layer::CreateCombiTrans(const char *name,
					 Double_t dy, Double_t dz,
					 Double_t dphi, Bool_t planeSym) {
//
// Help method to create a TGeoCombiTrans matrix
// (from a similar method with same name and function
// in AliITSv11GeometrySDD class by L.Gaudichet)
//

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
void AliITSUv2Layer::AddTranslationToCombiTrans(TGeoCombiTrans* ct,
                                                       Double_t dx,
                                                       Double_t dy,
                                                       Double_t dz) const{
//
// Help method to add a translation to a TGeoCombiTrans matrix
// (from a similar method with same name and function
// in AliITSv11GeometrySDD class by L.Gaudichet)
//

  // Add a dx,dy,dz translation to the initial TGeoCombiTrans
  const Double_t *vect = ct->GetTranslation();
  Double_t newVect[3] = {vect[0]+dx, vect[1]+dy, vect[2]+dz};
  ct->SetTranslation(newVect);
}

// #pnamwong


//________________________________________________________________________
void AliITSUv2Layer::CreateBarrelLayer(TGeoVolume* moth){
	
	const Double_t kHandleOffsetASide = 400*fgkmm/2 -85.0*fgkmm; // Offset from the center of InnerBCShell but measure from the edge of A-Side

	const Double_t kDCPSOffsetZ0 =  347.5*fgkmm/2 - (400.0-347.5)*fgkmm/2;
	const Double_t kDCPSOffsetZ[7] = {49.3*fgkmm, 45.23*fgkmm, 39.67*fgkmm, 0, 0, 0, 0};
	const Color_t kEWheelConeColor[7] = {kMagenta, kBlue, kYellow, 0, 0, 0, 0};

	TGeoVolume *layVol = new TGeoVolumeAssembly(TString("Barrel") += fLayerNumber);
	
	  
	if (fLayerNumber < fgkNumberOfInnerLayers)
	{
		
		// half down barrel
		//layVol->AddNode(CreateInnerBEWheelA2(kEWheelConeColor[(int)fLayerNumber],gGeoManager->GetMedium("ITS_KAPTON(POLYCH2)$")), 1, 
			//new TGeoCombiTrans(fgkInnerBarrelOffsetX, fgkInnerBarrelOffsetY, fgkInnerBarrelOffsetZ +kHandleOffsetASide, 
			//new TGeoRotation("", fgkBarrelRotateAlpha + 180.0, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));			
		

		// half barrel
		layVol->AddNode(CreateInnerBEWheelA3(kEWheelConeColor[(int)fLayerNumber], gGeoManager->GetMedium("ITS_PEEKCF30$")), 1,
			new TGeoCombiTrans(0, 0, fgkIBConnOffset, 
			new TGeoRotation("", 180, 0, 0)));			
		layVol->AddNode(CreateInnerBEWheelC3(kEWheelConeColor[(int)fLayerNumber], gGeoManager->GetMedium("ITS_PEEKCF30$")), 1,
			new TGeoCombiTrans(0, 0, -fgkIBConnOffset, 
			new TGeoRotation("", 180, 0, 0)));			
		layVol->AddNode(CreateInnerDCDCLayer3(), 1,
			new TGeoCombiTrans(0, 0, fgkIBConnOffset, 
			new TGeoRotation("", 180, 0, 0)));			

		// another half barrel
		layVol->AddNode(CreateInnerBEWheelA3(kEWheelConeColor[(int)fLayerNumber], gGeoManager->GetMedium("ITS_PEEKCF30$")), 1,
			new TGeoCombiTrans(0, 0, fgkIBConnOffset, 
			new TGeoRotation("", 0, 0, 0)));			
		layVol->AddNode(CreateInnerBEWheelC3(kEWheelConeColor[(int)fLayerNumber], gGeoManager->GetMedium("ITS_PEEKCF30$")), 1,
			new TGeoCombiTrans(0, 0, -fgkIBConnOffset, 
			new TGeoRotation("", 0, 0, 0)));			
		layVol->AddNode(CreateInnerDCDCLayer3(), 1,
			new TGeoCombiTrans(0, 0, fgkIBConnOffset, 
			new TGeoRotation("", 0, 0, 0)));			



		// still not modify
		//layVol->AddNode(CreateInnerBEWheelC2(kEWheelConeColor[(int)fLayerNumber],gGeoManager->GetMedium("ITS_KAPTON(POLYCH2)$")), 1, 
			//new TGeoCombiTrans(fgkInnerBarrelOffsetX, fgkInnerBarrelOffsetY, fgkInnerBarrelOffsetZ +kHandleOffsetASide, 
			//new TGeoRotation("", fgkBarrelRotateAlpha + 180.0, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));			
		//layVol->AddNode(CreateInnerBEWheelC(gGeoManager->GetMedium("ITS_KAPTON(POLYCH2)$")), 1, 
			//new TGeoCombiTrans(fgkInnerBarrelOffsetX, fgkInnerBarrelOffsetY, fgkInnerBarrelOffsetZ, 
			//new TGeoRotation("", fgkBarrelRotateAlpha + 0.0, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));
		//layVol->AddNode(CreateDCPSLayer(), 1, 
			//new TGeoCombiTrans(fgkInnerBarrelOffsetX, fgkInnerBarrelOffsetY, fgkInnerBarrelOffsetZ +kDCPSOffsetZ0 +kDCPSOffsetZ[(int)fLayerNumber], 
			//new TGeoRotation("", fgkBarrelRotateAlpha +90, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));

		// half up barrel
		//layVol->AddNode(CreateInnerBEWheelA2(kEWheelConeColor[(int)fLayerNumber],gGeoManager->GetMedium("ITS_KAPTON(POLYCH2)$")), 1, 
			//new TGeoCombiTrans(fgkInnerBarrelOffsetX, fgkInnerBarrelOffsetY, fgkInnerBarrelOffsetZ +kHandleOffsetASide, 
			//new TGeoRotation("", fgkBarrelRotateAlpha + 0.0, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));			
		//layVol->AddNode(CreateInnerBEWheelC(gGeoManager->GetMedium("ITS_KAPTON(POLYCH2)$")), 1, 
			//new TGeoCombiTrans(fgkInnerBarrelOffsetX, fgkInnerBarrelOffsetY, fgkInnerBarrelOffsetZ, 
			//new TGeoRotation("", fgkBarrelRotateAlpha + 180.0, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));
		//layVol->AddNode(CreateDCPSLayer(), 1, 
			//new TGeoCombiTrans(fgkInnerBarrelOffsetX, fgkInnerBarrelOffsetY, fgkInnerBarrelOffsetZ +kDCPSOffsetZ0 +kDCPSOffsetZ[(int)fLayerNumber], 
			//new TGeoRotation("", fgkBarrelRotateAlpha -90, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));

		if (fLayerNumber==0)
		{
			//layVol->AddNode(CreateInnerBCShell(), 1, 
				//new TGeoCombiTrans(fgkInnerBarrelOffsetX, fgkInnerBarrelOffsetY, fgkInnerBarrelOffsetZ, 
				//new TGeoRotation("", fgkBarrelRotateAlpha, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));
			//layVol->AddNode(CreateInnerBServiceB2(), 1,
				//new TGeoCombiTrans(fgkInnerBarrelOffsetX, fgkInnerBarrelOffsetY, fgkInnerBarrelOffsetZ +kHandleOffsetASide, 
				//new TGeoRotation("", fgkBarrelRotateAlpha + 180.0, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));
			/*
			TGeoVolume *volzTube;		
			
			// create XY axis at interaction point.
			volzTube= new TGeoVolume("IP0-axis", new TGeoTube(0,0.5*fgkmm,50*fgkcm), 0);
			volzTube->SetFillColor(2);
			volzTube->SetLineColor(2);
			moth->AddNode(volzTube,0, new TGeoCombiTrans(0, 0, 0, new TGeoRotation("", 0, 90, 0)));
			moth->AddNode(volzTube,0, new TGeoCombiTrans(0, 0, 0, new TGeoRotation("", 90, 90, 0)));
			moth->AddNode(volzTube,0, new TGeoCombiTrans(0, 0, 0, new TGeoRotation("", 0, 0, 0)));


			volzTube= new TGeoVolume("IP0", new TGeoTube(0, 100*fgkmm, 0.05*fgkmm), 0);
			volzTube->SetFillColor(2);
			volzTube->SetLineColor(2);
			moth->AddNode(volzTube,0, new TGeoCombiTrans(0, 0, 0, new TGeoRotation("", 0, 0, 0)));

			volzTube= new TGeoVolume("CalMatZone", new TGeoTube(0, 50*fgkmm, 0.1*fgkmm), 0);
			volzTube->SetFillColor(3);
			volzTube->SetLineColor(3);
			moth->AddNode(volzTube,0, new TGeoCombiTrans(0, 0, -fgkIBConnOffset, new TGeoRotation("", 0, 0, 0)));
			*/
			//volzTube= new TGeoVolume(TString("IP+")+=(fgkIBConnOffset*10), new TGeoTube(0, 50*fgkmm, 0.1*fgkmm), 0);
			//volzTube->SetFillColor(4);
			//volzTube->SetLineColor(4);
			//moth->AddNode(volzTube,0, new TGeoCombiTrans(0, 0, fgkIBConnOffset, new TGeoRotation("", 0, 0, 0)));

			//volzTube= new TGeoVolume(TString("IP-")+=(fgkIBConnOffset*10), new TGeoTube(0, 50*fgkmm, 0.1*fgkmm), 0);
			//volzTube->SetFillColor(4);
			//volzTube->SetLineColor(4);
			//moth->AddNode(volzTube,0, new TGeoCombiTrans(0, 0, -fgkIBConnOffset, new TGeoRotation("", 0, 0, 0)));

			// approximation test
			//const Double_t approxOff = 9.5*fgkmm;
			//volzTube= new TGeoVolume((TString("IP+")+=(fgkIBConnOffset*10)) + (TString("-")+=(approxOff*10)), new TGeoTube(0, 50*fgkmm, 0.1*fgkmm), 0);
			//volzTube->SetFillColor(5);
			//volzTube->SetLineColor(5);
			//moth->AddNode(volzTube,0, new TGeoCombiTrans(0, 0, fgkIBConnOffset-approxOff, new TGeoRotation("", 0, 0, 0)));

			//volzTube= new TGeoVolume((TString("IP-")+=(fgkIBConnOffset*10)) + (TString("+")+=(approxOff*10)), new TGeoTube(0, 50*fgkmm, 0.1*fgkmm), 0);
			//volzTube->SetFillColor(5);
			//volzTube->SetLineColor(5);
			//moth->AddNode(volzTube,0, new TGeoCombiTrans(0, 0, -fgkIBConnOffset+approxOff, new TGeoRotation("", 0, 0, 0)));


		}
		

	}
	else
	{
	  /*
		layVol->AddNode(CreateOuterBEWheelA(), 1, 
			new TGeoCombiTrans(fgkOuterBarrelOffsetX, fgkOuterBarrelOffsetY, fgkOuterBarrelOffsetZ, 
			new TGeoRotation("", fgkBarrelRotateAlpha, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));
		
		layVol->AddNode(CreateOuterBEWheelC(), 1, 
			new TGeoCombiTrans(fgkOuterBarrelOffsetX, fgkOuterBarrelOffsetY, fgkOuterBarrelOffsetZ, 
			new TGeoRotation("", fgkBarrelRotateAlpha, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));
	  */

	}
	
	moth->AddNode(layVol, 1, 0);
}


//______________________________________________________________________
void AliITSUv2Layer::CreateBarrel(TGeoVolume* moth){	

	// -----------------------------------------------------------------
	const Double_t kHandleOffsetASide = -400*fgkmm/2 +85.0*fgkmm; // Offset from the center of InnerBCShell but measure from the edge of A-Side
	const Double_t kDCPSOffsetZ0 = -347.5*fgkmm/2 + (400.0-347.5)*fgkmm/2;

	
	// Supporter of inner barrel
	moth->AddNode(CreateInnerBSupporterC(), 1, 
		new TGeoCombiTrans(fgkInnerBarrelOffsetX, fgkInnerBarrelOffsetY, fgkInnerBarrelOffsetZ + (40.0*fgkcm)/2, 
		new TGeoRotation("", fgkBarrelRotateAlpha, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));
	moth->AddNode(CreateInnerBSupporterC(), 1, 
		new TGeoCombiTrans(fgkInnerBarrelOffsetX, fgkInnerBarrelOffsetY, fgkInnerBarrelOffsetZ + (40.0*fgkcm)/2, 
		new TGeoRotation("", fgkBarrelRotateAlpha + 180, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));

	moth->AddNode(CreateInnerBSupporterA(), 1, 
		new TGeoCombiTrans(fgkInnerBarrelOffsetX, fgkInnerBarrelOffsetY, fgkInnerBarrelOffsetZ - (13.0*fgkcm), 
		new TGeoRotation("", fgkBarrelRotateAlpha, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));
	moth->AddNode(CreateInnerBSupporterA(), 1, 
		new TGeoCombiTrans(fgkInnerBarrelOffsetX, fgkInnerBarrelOffsetY, fgkInnerBarrelOffsetZ - (13.0*fgkcm), 
		new TGeoRotation("", fgkBarrelRotateAlpha + 180, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));
	
	// Service Barrel
	moth->AddNode(CreateInnerBServiceB2(), 1,
		new TGeoCombiTrans(fgkInnerBarrelOffsetX, fgkInnerBarrelOffsetY, fgkInnerBarrelOffsetZ +kHandleOffsetASide, 
		new TGeoRotation("", fgkBarrelRotateAlpha + 180.0, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));
	moth->AddNode(CreateInnerBServiceB2(), 1,
		new TGeoCombiTrans(fgkInnerBarrelOffsetX, fgkInnerBarrelOffsetY, fgkInnerBarrelOffsetZ +kHandleOffsetASide, 
		new TGeoRotation("", fgkBarrelRotateAlpha, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));
	// moth-> CreateOutterServiceB() NOT YET DONE.

	// -----------------------------------------------------------------
	// C-Shell of the detector barrel
	moth->AddNode(CreateInnerBCShell(), 1, 
		new TGeoCombiTrans(fgkInnerBarrelOffsetX, fgkInnerBarrelOffsetY, fgkInnerBarrelOffsetZ, 
		new TGeoRotation("", fgkBarrelRotateAlpha, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));
	moth->AddNode(CreateInnerBCShell(), 1, 
		new TGeoCombiTrans(fgkInnerBarrelOffsetX, fgkInnerBarrelOffsetY, fgkInnerBarrelOffsetZ, 
		new TGeoRotation("", fgkBarrelRotateAlpha + 180, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));
	
	moth->AddNode(CreateOuterBCShell(), 1, 
		new TGeoCombiTrans(fgkInnerBarrelOffsetX, fgkInnerBarrelOffsetY, fgkInnerBarrelOffsetZ, 
		new TGeoRotation("", fgkBarrelRotateAlpha, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));
	moth->AddNode(CreateOuterBCShell(), 1, 
		new TGeoCombiTrans(fgkInnerBarrelOffsetX, fgkInnerBarrelOffsetY, fgkInnerBarrelOffsetZ, 
		new TGeoRotation("", fgkBarrelRotateAlpha + 180, fgkBarrelRotateBetha, fgkBarrelRotateGamma)));

	// -----------------------------------------------------------------
	// ----- NOT USE ---- THIS IS A PROTOTYPE OF DCDC FOR TESTING ------
	// DC set
	// moth->AddNode(CreateDCSet(), 1, new TGeoCombiTrans(fgkInnerBarrelOffsetX*fgkmm, fgkInnerBarrelOffsetY-100*fgkmm, fgkInnerBarrelOffsetZ-690*fgkmm, new TGeoRotation("", fgkBarrelRotateAlpha+180, fgkBarrelRotateBetha+180+15, fgkBarrelRotateGamma)));

	//const Double_t xfRadius = 100*fgkmm;
	//const Int_t fNumOfStave = 12;
	//const Int_t NDegPerStave = 360/fNumOfStave;
	//const Double_t Deg2RadFactor = TMath::Pi()/180;

	
	//for (int i=0; i<fNumOfStave; i++)
	//{
		//Double_t CosDegIter = TMath::Cos(i*NDegPerStave*Deg2RadFactor);
		//Double_t SinDegIter = TMath::Sin(i*NDegPerStave*Deg2RadFactor);

		//moth->AddNode(CreateDCSet(), 1, 
			//new TGeoCombiTrans(fgkInnerBarrelOffsetX + (xfRadius * SinDegIter), fgkInnerBarrelOffsetY - (xfRadius * CosDegIter), fgkInnerBarrelOffsetZ-690*fgkmm, 
			//new TGeoRotation("", fgkBarrelRotateAlpha+180+i*NDegPerStave , fgkBarrelRotateBetha+180+15, fgkBarrelRotateGamma)));
	//}
	
	//, new TGeoRotation("", i*NDegPerStave + xfEachRotate, 90.0, 0.0)); // perpendicular to z axis

}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateInnerBEWheelA2(const Color_t color, const TGeoMedium* med){

	TGeoVolume* ewhVol = new TGeoVolumeAssembly(TString("EndWheelA") += fLayerNumber);

	const Double_t kOffsetZ = -18.0*fgkmm;
	const Double_t fkPlateThick = 0.4*fgkmm;
	const Double_t fkHndLength = 13.0*fgkmm; // also use in CreateInnerBHandleA.
		
	const Double_t kHndPlateRadius[7]={27.9*fgkmm, 35.9*fgkmm, 43.4*fgkmm,0,0,0,0};
	const Double_t kContactRadius[7]={25.9*fgkmm, 33.8*fgkmm, 41.1*fgkmm,0,0,0,0}; // approx : -2mm from Radius
		
	const Double_t fkHndWidthAngle = 360/fNStaves;

	static const Double_t kConeArray[(int)fgkNumberOfInnerLayers][6][3] =
		{		// Layer 0
			{{ 	  -0.00*fgkmm,  26.90*fgkmm,  28.90*fgkmm },	// Section 0
			{ 	 -20.35*fgkmm,  26.90*fgkmm,  28.90*fgkmm },	// Section 1
			{ 	 -23.50*fgkmm,  27.85*fgkmm,  28.90*fgkmm },	// Section 2
			{ 	-126.85*fgkmm,  58.95*fgkmm,  60.50*fgkmm },	// Section 3
			{ 	-127.00*fgkmm,  59.50*fgkmm,  60.50*fgkmm },	// Section 4
			{ 	-167.00*fgkmm,  59.50*fgkmm,  60.50*fgkmm }}	// Section 5
			,	// Layer 1
			{{ 	  -0.00*fgkmm,  34.90*fgkmm,  36.90*fgkmm },	// Section 0
			{ 	 -21.95*fgkmm,  34.90*fgkmm,  36.90*fgkmm },	// Section 1
			{ 	 -23.50*fgkmm,  35.75*fgkmm,  36.90*fgkmm },	// Section 2
			{ 	-121.56*fgkmm,  89.96*fgkmm,  91.10*fgkmm },	// Section 3
			{ 	-121.84*fgkmm,  90.10*fgkmm,  91.10*fgkmm },	// Section 4
			{ 	-161.80*fgkmm,  90.10*fgkmm,  91.10*fgkmm }}	// Section 5
			,	// Layer 2
			{{ 	  -0.00*fgkmm,  42.90*fgkmm,  44.90*fgkmm },	// Section 0
			{ 	 -22.72*fgkmm,  42.90*fgkmm,  44.90*fgkmm },	// Section 1
			{ 	 -23.50*fgkmm,  43.57*fgkmm,  44.90*fgkmm },	// Section 2
			{ 	-110.30*fgkmm, 118.77*fgkmm, 120.10*fgkmm },	// Section 3
			{ 	-110.73*fgkmm, 119.10*fgkmm, 120.10*fgkmm },	// Section 4
			{ 	-150.70*fgkmm, 119.10*fgkmm, 120.10*fgkmm }}	// Section 5
		};
		
	static const Double_t kConeTubeArray[(int)fgkNumberOfInnerLayers][2][3] =
		{		// Layer 0
			{{ 	-127.00*fgkmm,  28.60*fgkmm,  29.60*fgkmm },	// Section 0
			{ 	-167.00*fgkmm,  28.60*fgkmm,  29.60*fgkmm }}	// Section 1
			,	// Layer 1
			{{ 	-127.00*fgkmm,  60.50*fgkmm,  61.50*fgkmm },	// Section 0
			{ 	-161.80*fgkmm,  60.50*fgkmm,  61.50*fgkmm }}	// Section 1
			,	// Layer 2
			{{ 	-121.84*fgkmm,  91.10*fgkmm,  92.10*fgkmm },	// Section 0
			{ 	-150.70*fgkmm,  91.10*fgkmm,  92.10*fgkmm }}	// Section 1
		};
	
	
	if (fLayerNumber < fgkNumberOfInnerLayers)
	{
		TGeoVolume* vol;
		TGeoShape* shp;
		
		// create handles
		for (Int_t i=0; i<fNStaves/2; i++)
		{
			vol = CreateInnerBHandleA((TString("_")+=(int)i).Data(), kHndPlateRadius[(int)fLayerNumber], kContactRadius[(int)fLayerNumber], med);
			vol->SetFillColor(color);
			vol->SetLineColor(color);
			ewhVol->AddNode(vol, 1, new TGeoCombiTrans(0, 0, 0, new TGeoRotation("", 0, 0, +fkHndWidthAngle/2 +i*fkHndWidthAngle)));
		}
				
		shp = new TGeoPcon(TString("WheelConeA")+=fLayerNumber, 0, 180, 6);
		for (int i=0; i<6; i++)
			((TGeoPcon*)shp)->DefineSection(i, kConeArray[(int)fLayerNumber][i][0], kConeArray[(int)fLayerNumber][i][1], kConeArray[(int)fLayerNumber][i][2]);	
		vol = new TGeoVolume(TString("WheelConeA")+=fLayerNumber, shp, med);
		vol->SetFillColor(color);
		vol->SetLineColor(color);
		ewhVol->AddNode(vol, 1, new TGeoCombiTrans(0, 0, +kOffsetZ, new TGeoRotation("", 0, 0, 0)));
		
		if (fLayerNumber==0)
		{
			shp = new TGeoPcon(TString("WheelConeTubeA")+=fLayerNumber, 0, 180, 2);
			for (int i=0; i<2; i++)
				((TGeoPcon*)shp)->DefineSection(i, kConeTubeArray[(int)fLayerNumber][i][0], kConeArray[(int)fLayerNumber][i][1], kConeArray[(int)fLayerNumber][i][2]);	
			vol = new TGeoVolume(TString("WheelConeTubeA")+=fLayerNumber, shp, med);
			vol->SetFillColor(color);
			vol->SetLineColor(color);
			ewhVol->AddNode(vol, 1, new TGeoCombiTrans(0, 0, +kOffsetZ, new TGeoRotation("", 0, 0, 0)));
		}
	}

	return ewhVol;
}


TGeoVolume* AliITSUv2Layer::CreateInnerBEWheelA3(const Color_t color, const TGeoMedium* med){
	// Endwheel for half stave
	TGeoVolume *vALIC = gGeoManager->GetVolume("ALIC");
	
	TGeoVolume* ewhVol = new TGeoVolumeAssembly(TString("EndWheelA") += fLayerNumber);

	TGeoVolume* vol;
	TGeoShape* shp;

	vol = Create_ALC_0336_xxA_ins(med);					// Endwheel, inside the wrap volume
	vol->SetFillColor(color);
	vol->SetLineColor(color);
	ewhVol->AddNode(vol, 1,	new TGeoCombiTrans(0, 0, -fkALC_0334_HoleFromEdge, new TGeoRotation("", 0, 0, 0)));			

	
	for (int i=0; i<fNStaves/2; i++)					// Handles
	{
		vol = Create_ALC_0334_xxA(med);
		vol->SetFillColor(color);
		vol->SetLineColor(color);
		ewhVol->AddNode(vol, 1, 
			new TGeoCombiTrans(0, 0, 0, 
			new TGeoRotation("", 0, 0, -90. + (fkALC_0334_RotateAngle[fLayerNumber]*i) + fkALC_0334_RotateOffset[fLayerNumber])));						
	}		

	ewhVol->AddNode(Create_ALC_0336_xxA_ous(color, med), 1, // Endwheel, outside the wrap volume
		new TGeoCombiTrans(0,0,0,
		new TGeoRotation("", 0,0,0)));

	//ewhVol->GetShape()->ComputeBBox(); //RS: enfore recompting of BBox
	
	return ewhVol;
}

TGeoVolume* AliITSUv2Layer::CreateInnerDCDCLayer3(const TGeoManager *mgr){
	// DCDC Set of half barrel
	//TGeoVolume *vALIC = gGeoManager->GetVolume("ALIC");
	
	TGeoVolume* ewhVol = new TGeoVolumeAssembly(TString("DCDCSET") += fLayerNumber);

	TGeoVolume* vol;
	TGeoShape* shp;

	for (int i=0; i<fNStaves/2; i++)					// Handles
	{
		ewhVol->AddNode(CreateInnerDCCNT3(), 1,
			new TGeoCombiTrans(0, 0, +fkThermalConnOffsetZ[fLayerNumber], 
			new TGeoRotation("", 0, 0, -90. + (fkALC_0334_RotateAngle[fLayerNumber]*i) + fkALC_0334_RotateOffset[fLayerNumber])));			

		ewhVol->AddNode(CreateInnerDCDC3(), 1,
			new TGeoCombiTrans(0, 0, +fkDCDCOffsetZ[fLayerNumber], 
			new TGeoRotation("", 0, 0, -90. + (fkALC_0334_RotateAngle[fLayerNumber]*i) + fkALC_0334_RotateOffset[fLayerNumber])));			

	}		



	
	return ewhVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateInnerBEWheelC(const TGeoMedium* med){
	//
	//
	//
	Double_t fHoleRadius =   5.*fgkmm;
	Double_t fHoleFromEdge = 4.*fgkmm;		// the center of the hole
	
	Double_t fOffsetZ = 0.;
	Double_t fRadius = 	0.;
	Double_t fLength1 = 0.;
	Double_t fLength2 = 0.;
	Double_t fThick1 = 	0.;
	Double_t fThick2 = 	0.;
	Double_t fColor = 	0.;


	Double_t fHoleFromEdge1 = fHoleFromEdge - fHoleRadius/2;	// the border of the hole
	Double_t fHoleFromEdge2 = fHoleFromEdge + fHoleRadius/2;	// the border of the hole

	TGeoVolume* ewhVol = new TGeoVolumeAssembly(TString("EndWheelC") += fLayerNumber);

	switch (fLayerNumber)
	{
		case 0:
			fOffsetZ = (400.0/2-4.0)*fgkmm;
			fRadius = 23.5*fgkmm;
			fLength1 = 0.6*fgkmm;
			fLength2 = 18.4*fgkmm;
			fThick1 = 4.4*fgkmm;
			fThick2 = 0.6*fgkmm;
			fColor = 6;

			// Supporter (C-Side)
			ewhVol->AddNode(CreateInnerBSupporterC(), 1, new TGeoTranslation(0.,0.,40.*fgkcm/2));
		break;
		
		case 1:
			fOffsetZ = (400.0/2-6.0)*fgkmm;
			fRadius = 31.5*fgkmm;
			fLength1 = 0.6*fgkmm;
			fLength2 = 16.4*fgkmm;
			fThick1 = 4.4*fgkmm;
			fThick2 = 0.6*fgkmm;
			fColor = 5;

		break;

		case 2:
			fOffsetZ = (400.0/2-8.0)*fgkmm;
			fRadius = 38.0*fgkmm;
			fLength1 = 0.6*fgkmm;
			fLength2 = 14.4*fgkmm;
			fThick1 = 5.4*fgkmm;
			fThick2 = 0.6*fgkmm;	
			fColor = 7;		
		break;
		
		default:;
		
	}
	

	ewhVol->AddNode(CreateLTube("EndWheelObj[1]", fRadius, fLength1, fLength2-fHoleFromEdge2, fThick1, fThick2, kTRUE, kFALSE, fColor, med), 1, 
		new TGeoCombiTrans(0., 0., fOffsetZ, new TGeoRotation("",0.,0.,180.)));

	ewhVol->AddNode(CreateTube("EndWheelObj[2]", fRadius+fThick1, fHoleFromEdge1, 0.6*fgkmm, kFALSE, fColor, 180, med), 1, 
		new TGeoCombiTrans(0., 0., fOffsetZ-fLength1-fLength2+fHoleFromEdge1, new TGeoRotation("",0.,0.,180.)));



	return ewhVol;
}

TGeoVolume* AliITSUv2Layer::CreateInnerBEWheelC3(const Color_t color, const TGeoMedium* med){
	// Endwheel for half stave
	TGeoVolume *vALIC = gGeoManager->GetVolume("ALIC");
	
	TGeoVolume* ewhVol = new TGeoVolumeAssembly(TString("EndWheelC") += fLayerNumber);

	TGeoVolume* vol;
	TGeoShape* shp;

	vol = Create_ALC_0337_xxC_ins(med);					// Plate inside warp volume
	vol->SetFillColor(color);
	vol->SetLineColor(color);
	ewhVol->AddNode(vol, 1,	new TGeoCombiTrans(0, 0, +fkALC_0334_HoleFromEdge, new TGeoRotation("", 0, 0, 0)));			

	
	for (int i=0; i<fNStaves/2; i++)					
	{
		vol = Create_ALC_0334_xxC(med);					// Handles
		vol->SetFillColor(color);
		vol->SetLineColor(color);
		ewhVol->AddNode(vol, 1, new TGeoCombiTrans(0, 0, 0, 
			new TGeoRotation("", 0, 0, -90. + (fkALC_0334_RotateAngle[fLayerNumber]*i) + fkALC_0334_RotateOffset[fLayerNumber])));			
	}		

	ewhVol->AddNode(Create_ALC_0337_xxC_ous(color, med), 1, // Endwheel, outside wrap volume
		new TGeoCombiTrans(0, 0, 0,
		new TGeoRotation("", 0,0,0)));

	//ewhVol->GetShape()->ComputeBBox(); //RS: enfore recompting of BBox

	return ewhVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateOuterBEWheelC(const TGeoManager *mgr){
	//
	//
	//
	

	TGeoVolume* ewhVol = new TGeoVolumeAssembly(TString("EndWheelC") += fLayerNumber);

	/*
	const Double_t ITSSBOB01SWSkinThick = 0.60 *fgkmm;
	const Double_t ITSSBOB01SWCoreThick = 8.80 *fgkmm;
	const Double_t ITSSBOB01SWThick = ITSSBOB01SWCoreThick + 2*ITSSBOB01SWSkinThick;

	// GLOBAL VALUE FOR ALL SUB OBJECTS
	const Double_t ITSSBOB01CutOffset = 0.20;
	const Double_t ITSSBOB01Length = 1740.00 *fgkmm;
	const Double_t ITSSBOB01Radius =  449.00 *fgkmm; 
	const Double_t ITSSBOB01MinRadius = 2.20;
	const Double_t ITSSBOB01MaxRadius = ITSSBOB01Radius + ITSSBOB01SWThick;

	// box to div
	TGeoBBox *xbboxdiv = new TGeoBBox("xbboxdiv2", ITSSBOB01MaxRadius+1.0, ITSSBOB01MaxRadius+1.0, (ITSSBOB01Length/2)+1.0); // extend edge about 1.0
	TGeoTranslation *xtrandiv = new TGeoTranslation("xtrandiv2", 0, 0-xbboxdiv->GetDY()-ITSSBOB01CutOffset, 0);
	xtrandiv->RegisterYourself();			


	if (fLayerNumber==3)
	{ // 0113x : Layer 3 handle

		// 01131 are the base ( L tube )
		// 01132 : NOT USE
		// 01133x are the handles
		// 01134 is the tube for div to make a hole in each base 1113x
		// 01135 is a tube for div to cut unwanted outside the tube
		// to make sure we can welded 01133x with (01131&01132), set 01135Radius to have more a bit value than 01132Radius 
		
		Double_t ITSSBOB01131Length1 =  3.0 *fgkmm; // fixed parameter
		Double_t ITSSBOB01131Length2 = 24.0 *fgkmm; // fixed parameter

		Double_t ITSSBOB01131Radius = 212.0 *fgkmm;
		
		Double_t ITSSBOB01131Thick1 =  35.0 *fgkmm;
		Double_t ITSSBOB01131Thick2 =  14.0 *fgkmm;
		
		// offset is compared with volume's shape of section1
		TGeoTranslation *ITSSBOB01131Offset = new TGeoTranslation(0.0, 0.0, 0.0 + ITSSBOB01Length/2 - 315.0*fgkmm);

		Double_t ITSSBOB01131Z0 = 0.00; // The center of object is not the same as the center of space.
		Double_t ITSSBOB01131Z1 = -(ITSSBOB01131Length1);
		Double_t ITSSBOB01131Z2 = -(ITSSBOB01131Length1);
		Double_t ITSSBOB01131Z3 = -(ITSSBOB01131Length1 + ITSSBOB01131Length2);

		Double_t ITSSBOB01131RMin0 = ITSSBOB01131Radius, ITSSBOB01131RMax0 = ITSSBOB01131Radius + ITSSBOB01131Thick1;
		Double_t ITSSBOB01131RMin1 = ITSSBOB01131Radius, ITSSBOB01131RMax1 = ITSSBOB01131Radius + ITSSBOB01131Thick1;
		Double_t ITSSBOB01131RMin2 = ITSSBOB01131Radius, ITSSBOB01131RMax2 = ITSSBOB01131Radius + ITSSBOB01131Thick2;
		Double_t ITSSBOB01131RMin3 = ITSSBOB01131Radius, ITSSBOB01131RMax3 = ITSSBOB01131Radius + ITSSBOB01131Thick2;										

		TGeoPcon *xpconeITSSBOB01131 = new TGeoPcon("xpconeITSSBOB01131", 0, 360, 4);
		xpconeITSSBOB01131->DefineSection(0, ITSSBOB01131Z0, ITSSBOB01131RMin0, ITSSBOB01131RMax0);
		xpconeITSSBOB01131->DefineSection(1, ITSSBOB01131Z1, ITSSBOB01131RMin1, ITSSBOB01131RMax1);
		xpconeITSSBOB01131->DefineSection(2, ITSSBOB01131Z2, ITSSBOB01131RMin2, ITSSBOB01131RMax2);
		xpconeITSSBOB01131->DefineSection(3, ITSSBOB01131Z3, ITSSBOB01131RMin3, ITSSBOB01131RMax3);

		// --------------- CREATE HOLES ---------------------

		
		TGeoRotation *xrotITSSBOB0113Circle = new TGeoRotation("xrotITSSBOB0113Circle", 0.0, 0.0, 0.0); // this rotate BEFORE div
		xrotITSSBOB0113Circle->RegisterYourself();

		vector<TString> xsholes = CreateDBOuterHoles("ITSSBOB01131"
			, 210.5 *fgkmm	// faces radius
			,  (10.0/2) *fgkmm  // holes radius
			,  (15.0*2.5) *fgkmm	// over approximate, for div
			, 0.0 // not rotate
			, - 7.0*fgkmm - ((ITSSBOB01131Length1+ITSSBOB01131Length2)/2));


		// --------------- CREATE HANDLES -------------------
		
		vector<TString> xshandles = CreateDBOuterHandles("ITSSBOB01131"
			//, 210.5 *fgkmm	// faces radius
			, 207.5 *fgkmm	// faces radius, JUST FOR SHOW
			,  27.0 *fgkmm  // length, the same as base length
			,  50.5 *fgkmm  // width, overvalue, FOR SHOW : JUST REAL VALUE
			,   3.0 *fgkmm  // thick, overvalue
			, 0.0 // not rotate
			, + (27.0*fgkmm)/2 - (ITSSBOB01131Length1+ITSSBOB01131Length2));		

		// ----------------------- combine all 0111x objects
		TString xstrITSSBOB0113cs("(");
		xstrITSSBOB0113cs += "xpconeITSSBOB01131";
	
		
		for (Int_t i=0; i<fNStaves; i++) // ITSSBOB01133
		{ // *** IF PLACE CODE 114 BELOW 113, ITS POSITION MIGHT BE UNDETERMINED ***
			xstrITSSBOB0113cs += " + (";
			xstrITSSBOB0113cs += xshandles[i];
			
			#ifdef _ITSSB_SHOW_COMBISHAPE
				xstrITSSBOB0113cs += " - ";		// make a hole
			#else
				xstrITSSBOB0113cs += " + ";		// just show
			#endif
			xstrITSSBOB0113cs += xsholes[i];												
			xstrITSSBOB0113cs += ")"; 
		}

														
		xstrITSSBOB0113cs += "):xrotITSSBOB0113Circle";
		
		#ifdef _ITSSB_SHOW_WITHBOXDIV
		xstrITSSBOB0113cs += " * (xbboxdiv2:xtrandiv2)";
		#endif 
		

		TGeoVolume *xvolITSSBOB0113 = new TGeoVolume("ITSSBOB0113x", 
			new TGeoCompositeShape("xcompITSSB0B0113cs", xstrITSSBOB0113cs), 
			mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 
		xvolITSSBOB0113->SetFillColor(9);
		xvolITSSBOB0113->SetLineColor(9);
		ewhVol->AddNode(xvolITSSBOB0113, 1, ITSSBOB01131Offset);
		
	}
	
	if (fLayerNumber==4)
	{ 


		// 01141 is the base ( L tube )
		// 01142 is the base ( L tube )
		// 01143x are the handles
		// 01144 is the tube for div to make a hole in each base 1113x
		// 01145 is a tube for div to cut unwanted outside the tube
		// to make sure we can welded 01143x with (01141&01142), set 01145Radius to have more a bit value than 01142Radius 
		
		
		Double_t ITSSBOB01141Length1 =  8.0 *fgkmm; // fixed parameter
		Double_t ITSSBOB01141Length2 = 16.0 *fgkmm; // fixed parameter

		Double_t ITSSBOB01141Radius = 226.0 *fgkmm;
		
		Double_t ITSSBOB01141Thick1 =  48.0 *fgkmm;
		Double_t ITSSBOB01141Thick2 =   9.5 *fgkmm;

		
		// offset is compared with volume's shape of section1
		TGeoTranslation *ITSSBOB01141Offset = new TGeoTranslation(0.0, 0.0, 0.0 + ITSSBOB01Length/2 - 318.0*fgkmm);

		Double_t ITSSBOB01141Z0 = 0.00; // The center of object is not the same as the center of space.
		Double_t ITSSBOB01141Z1 = -(ITSSBOB01141Length1);
		Double_t ITSSBOB01141Z2 = -(ITSSBOB01141Length1);
		Double_t ITSSBOB01141Z3 = -(ITSSBOB01141Length1 + ITSSBOB01141Length2);

		Double_t ITSSBOB01141RMin0 = ITSSBOB01141Radius, ITSSBOB01141RMax0 = ITSSBOB01141Radius + ITSSBOB01141Thick1;
		Double_t ITSSBOB01141RMin1 = ITSSBOB01141Radius, ITSSBOB01141RMax1 = ITSSBOB01141Radius + ITSSBOB01141Thick1;
		Double_t ITSSBOB01141RMin2 = ITSSBOB01141Radius + ITSSBOB01141Thick1 - ITSSBOB01141Thick2, ITSSBOB01141RMax2 = ITSSBOB01141Radius + ITSSBOB01141Thick1;
		Double_t ITSSBOB01141RMin3 = ITSSBOB01141Radius + ITSSBOB01141Thick1 - ITSSBOB01141Thick2, ITSSBOB01141RMax3 = ITSSBOB01141Radius + ITSSBOB01141Thick1;										

		TGeoPcon *xpconeITSSBOB01141 = new TGeoPcon("xpconeITSSBOB01141", 0, 360, 4);
		xpconeITSSBOB01141->DefineSection(0, ITSSBOB01141Z0, ITSSBOB01141RMin0, ITSSBOB01141RMax0);
		xpconeITSSBOB01141->DefineSection(1, ITSSBOB01141Z1, ITSSBOB01141RMin1, ITSSBOB01141RMax1);
		xpconeITSSBOB01141->DefineSection(2, ITSSBOB01141Z2, ITSSBOB01141RMin2, ITSSBOB01141RMax2);
		xpconeITSSBOB01141->DefineSection(3, ITSSBOB01141Z3, ITSSBOB01141RMin3, ITSSBOB01141RMax3);

		// --------------- CREATE HOLES ---------------------

		
		TGeoRotation *xrotITSSBOB0114Circle = new TGeoRotation("xrotITSSBOB0114Circle", 0.0, 0.0, 0.0); // this rotate BEFORE div
		xrotITSSBOB0114Circle->RegisterYourself();

		vector<TString> xsholes = CreateDBOuterHoles("ITSSBIB01141"
			, 263.5 *fgkmm	// faces radius
			,  (10.0/2) *fgkmm  // holes radius
			,  (15.0*2.5) *fgkmm	// over approximate, for div
			, 0.0 // not rotate
			, - 7.0*fgkmm - ((ITSSBOB01141Length1+ITSSBOB01141Length2)/2));

		// --------------- CREATE HANDLES -------------------
		
		vector<TString> xshandles = CreateDBOuterHandles("ITSSBOB01141"
			//, 263.5 *fgkmm	// faces radius
			, 260.5 *fgkmm	// faces radius, JUST FOR SHOW
			,  14.0 *fgkmm  // length, the same as base length
			,  46.2 *fgkmm  // width, overvalue, FOR SHOW : JUST REAL VALUE
			,   3.0 *fgkmm  // thick, overvalue
			, 0.0 // not rotate
			, + (16.0*fgkmm)/2 - (ITSSBOB01141Length1+ITSSBOB01141Length2));		

		// ------------------------- combine all 0111x objects
		TString xstrITSSBOB0114cs("(");
		xstrITSSBOB0114cs += "xpconeITSSBOB01141";
	
		
		for (Int_t i=0; i<fNStaves; i++) // ITSSBOB01143
		{ // *** IF PLACE CODE 114 BELOW 113, ITS POSITION MIGHT BE UNDETERMINED ***
			xstrITSSBOB0114cs += " + (";
			xstrITSSBOB0114cs += xshandles[i];
			
			#ifdef _ITSSB_SHOW_COMBISHAPE
				xstrITSSBOB0114cs += " - ";		// make a hole
			#else
				xstrITSSBOB0114cs += " + ";		// just show
			#endif
			xstrITSSBOB0114cs += xsholes[i];												
			xstrITSSBOB0114cs += ")"; 
		}

														
		xstrITSSBOB0114cs += "):xrotITSSBOB0114Circle";
		
		#ifdef _ITSSB_SHOW_WITHBOXDIV
		xstrITSSBOB0114cs += " * (xbboxdiv2:xtrandiv2)";
		#endif 
		

		TGeoVolume *xvolITSSBOB01141 = new TGeoVolume("ITSSBOB0114x", new TGeoCompositeShape("xcompITSSB0B0114cs", xstrITSSBOB0114cs), mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 
		xvolITSSBOB01141->SetFillColor(7);
		xvolITSSBOB01141->SetLineColor(7);
		ewhVol->AddNode(xvolITSSBOB01141, 1, ITSSBOB01141Offset);
		
		
		// ----------- 01142
		Double_t ITSSBOB01142Length1 =  3.0 *fgkmm; // fixed parameter
		Double_t ITSSBOB01142Length2 = 15.0 *fgkmm; // fixed parameter

		Double_t ITSSBOB01142Radius = 321.0 *fgkmm;
		
		Double_t ITSSBOB01142Thick1 =  30.0 *fgkmm;
		Double_t ITSSBOB01142Thick2 =   3.0 *fgkmm;

		// offset is compared with volume's shape of section1
		TGeoTranslation *ITSSBOB01142Offset = new TGeoTranslation(0.0, 0.0, 0.0 + ITSSBOB01Length/2 - 2.0*fgkmm);

		Double_t ITSSBOB01142Z0 = 0.00; // The center of object is not the same as the center of space.
		Double_t ITSSBOB01142Z1 = -(ITSSBOB01142Length1);
		Double_t ITSSBOB01142Z2 = -(ITSSBOB01142Length1);
		Double_t ITSSBOB01142Z3 = -(ITSSBOB01142Length1 + ITSSBOB01142Length2);

		Double_t ITSSBOB01142RMin0 = ITSSBOB01142Radius, ITSSBOB01142RMax0 = ITSSBOB01142Radius + ITSSBOB01142Thick1;
		Double_t ITSSBOB01142RMin1 = ITSSBOB01142Radius, ITSSBOB01142RMax1 = ITSSBOB01142Radius + ITSSBOB01142Thick1;
		Double_t ITSSBOB01142RMin2 = ITSSBOB01142Radius, ITSSBOB01142RMax2 = ITSSBOB01142Radius + ITSSBOB01142Thick2;
		Double_t ITSSBOB01142RMin3 = ITSSBOB01142Radius, ITSSBOB01142RMax3 = ITSSBOB01142Radius + ITSSBOB01142Thick2;										

		TGeoPcon *xpconeITSSBOB01142 = new TGeoPcon("xpconeITSSBOB01142", 0, 360, 4);
		xpconeITSSBOB01142->DefineSection(0, ITSSBOB01142Z0, ITSSBOB01142RMin0, ITSSBOB01142RMax0);
		xpconeITSSBOB01142->DefineSection(1, ITSSBOB01142Z1, ITSSBOB01142RMin1, ITSSBOB01142RMax1);
		xpconeITSSBOB01142->DefineSection(2, ITSSBOB01142Z2, ITSSBOB01142RMin2, ITSSBOB01142RMax2);
		xpconeITSSBOB01142->DefineSection(3, ITSSBOB01142Z3, ITSSBOB01142RMin3, ITSSBOB01142RMax3);


		//TGeoCompositeShape *xcompITSSBOB0113 = new TGeoCompositeShape("xcompITSSBOB0113cs", xstrITSSBOB0113cs);
		//TGeoVolume *xvolITSSBOB0113 = new TGeoVolume("ITSSBOB0113x", xcompITSSBOB0113, med); 			
		#ifdef _ITSSB_SHOW_WITHBOXDIV
		TGeoVolume *xvolITSSBOB01142 = new TGeoVolume("ITSSBOB01142", new TGeoCompositeShape("xpconeITSSBOB01142*(xbboxdiv2:xtrandiv2)"), mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 					
		#else		
		TGeoVolume *xvolITSSBOB01142 = new TGeoVolume("ITSSBOB01142", xpconeITSSBOB01142, mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 
		#endif
		xvolITSSBOB01142->SetFillColor(7);
		xvolITSSBOB01142->SetLineColor(7);
		ewhVol->AddNode(xvolITSSBOB01142, 1, ITSSBOB01142Offset);
		
		//ewhVol->AddNode(CreateSTube45("OuterLTube", 274.0*fgkmm, 324.0*fgkmm, 337.0*fgkmm, 151.5 *fgkmm, 5.0*fgkmm + 2*0.5*fgkmm, 0, 7), 1, new TGeoCombiTrans(0.0, 0.0, 0.0 + (174.0*fgkcm)/2 - (5.0*fgkmm), new TGeoRotation("",0,0,180)), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));


		//ewhVol->AddNode(CreateLTube("Test", 100, 2, 10, 10, 2, 0, 1, 7), 1, new TGeoCombiTrans(0.0, 0.0, 0.0, new TGeoRotation("",0,0,180)), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	}


	if (fLayerNumber==5)
	{ // 0115x : Layer 5 handle

		// 01151 is the base ( tube )
		// 01152 is the base
		// 01153x are the handles
		// 01154 is the tube for div to make a hole in each base 1113x
		// 01155 is a tube for div to cut unwanted outside the tube
		// to make sure we can welded 01153x with (01151&01152), set 01155Radius to have more a bit value than 01152Radius 

		// ------------------------ 01151 ---------------------
							
		Double_t ITSSBOB01151Length1 =  3.0 *fgkmm; // fixed parameter
		Double_t ITSSBOB01151Length2 =  8.0 *fgkmm; // fixed parameter
		Double_t ITSSBOB01151Length3 = 16.0 *fgkmm; // fixed parameter

		Double_t ITSSBOB01151Radius = 330.0 *fgkmm;
		
		Double_t ITSSBOB01151Thick1 =  33.0 *fgkmm;
		Double_t ITSSBOB01151Thick2 =  77.0 *fgkmm;
		Double_t ITSSBOB01151Thick3 =  53.0 *fgkmm;
		Double_t ITSSBOB01151Thick4 =  41.0 *fgkmm;
		
		
		// offset is compared with volume's shape of section1
		TGeoTranslation *ITSSBOB01151Offset = new TGeoTranslation(0.0, 0.0, 0.0 + ITSSBOB01Length/2 - 2.0*fgkmm);

		Double_t ITSSBOB01151Z0 = 0.00; // The center of object is not the same as the center of space.
		Double_t ITSSBOB01151Z1 = -(ITSSBOB01151Length1);
		Double_t ITSSBOB01151Z2 = -(ITSSBOB01151Length1);
		Double_t ITSSBOB01151Z3 = -(ITSSBOB01151Length1 + ITSSBOB01151Length2);
		Double_t ITSSBOB01151Z4 = -(ITSSBOB01151Length1 + ITSSBOB01151Length2);
		Double_t ITSSBOB01151Z5 = -(ITSSBOB01151Length1 + ITSSBOB01151Length2 + ITSSBOB01151Length3);

		Double_t ITSSBOB01151RMin0 = ITSSBOB01151Radius + ITSSBOB01151Thick1, ITSSBOB01151RMax0 = ITSSBOB01151Radius + ITSSBOB01151Thick2;
		Double_t ITSSBOB01151RMin1 = ITSSBOB01151Radius + ITSSBOB01151Thick1, ITSSBOB01151RMax1 = ITSSBOB01151Radius + ITSSBOB01151Thick2;
		Double_t ITSSBOB01151RMin2 = ITSSBOB01151Radius, ITSSBOB01151RMax2 = ITSSBOB01151Radius + ITSSBOB01151Thick3;
		Double_t ITSSBOB01151RMin3 = ITSSBOB01151Radius, ITSSBOB01151RMax3 = ITSSBOB01151Radius + ITSSBOB01151Thick3;
		Double_t ITSSBOB01151RMin4 = ITSSBOB01151Radius + ITSSBOB01151Thick4, ITSSBOB01151RMax4 = ITSSBOB01151Radius + ITSSBOB01151Thick3;
		Double_t ITSSBOB01151RMin5 = ITSSBOB01151Radius + ITSSBOB01151Thick4, ITSSBOB01151RMax5 = ITSSBOB01151Radius + ITSSBOB01151Thick3;
												
		TGeoPcon *xpconeITSSBOB01151 = new TGeoPcon("xpconeITSSBOB01151", 0, 360, 6);
		xpconeITSSBOB01151->DefineSection(0, ITSSBOB01151Z0, ITSSBOB01151RMin0, ITSSBOB01151RMax0);
		xpconeITSSBOB01151->DefineSection(1, ITSSBOB01151Z1, ITSSBOB01151RMin1, ITSSBOB01151RMax1);
		xpconeITSSBOB01151->DefineSection(2, ITSSBOB01151Z2, ITSSBOB01151RMin2, ITSSBOB01151RMax2);
		xpconeITSSBOB01151->DefineSection(3, ITSSBOB01151Z3, ITSSBOB01151RMin3, ITSSBOB01151RMax3);
		xpconeITSSBOB01151->DefineSection(4, ITSSBOB01151Z4, ITSSBOB01151RMin4, ITSSBOB01151RMax4);
		xpconeITSSBOB01151->DefineSection(5, ITSSBOB01151Z5, ITSSBOB01151RMin5, ITSSBOB01151RMax5);

		// --------------- CREATE HOLES ---------------------

		
		TGeoRotation *xrotITSSBOB01151Circle = new TGeoRotation("xrotITSSBOB01151Circle", 0.0, 0.0, 0.0); // this rotate BEFORE div
		xrotITSSBOB01151Circle->RegisterYourself();

		vector<TString> xsholes = CreateDBOuterHoles("ITSSBIB01151"
			, 370.0 *fgkmm	// faces radius
			,  (10.0/2) *fgkmm  // holes radius
			,  (15.0*2.5) *fgkmm	// over approximate, for div
			, 0.0 // not rotate
			, - 7.0*fgkmm - ((ITSSBOB01151Length1+ITSSBOB01151Length2+ITSSBOB01151Length3)/2));

		// --------------- CREATE HANDLES -------------------
		
		vector<TString> xshandles = CreateDBOuterHandles("ITSSBOB01151"
			//, 370.0 *fgkmm	// faces radius
			, 367.0 *fgkmm	// faces radius, JUST FOR SHOW
			,  14.0 *fgkmm  // length, the same as base length
			,  54.4 *fgkmm  // width, overvalue, FOR SHOW : JUST REAL VALUE
			,   3.0 *fgkmm  // thick, overvalue
			, 0.0 // not rotate
			, + (14.0*fgkmm)/2 - (ITSSBOB01151Length1+ITSSBOB01151Length2+ITSSBOB01151Length3));		

		// ------------------------ combine all 0111x objects
		TString xstrITSSBOB01151cs("(");
		xstrITSSBOB01151cs += "xpconeITSSBOB01151";
		
		for (Int_t i=0; i<fNStaves; i++) // ITSSBOB011513
		{ // *** IF PLACE CODE 114 BELOW 113, ITS POSITION MIGHT BE UNDETERMINED ***
			xstrITSSBOB01151cs += " + (";
			xstrITSSBOB01151cs += xshandles[i];
			
			#ifdef _ITSSB_SHOW_COMBISHAPE
				xstrITSSBOB01151cs += " - ";		// make a hole
			#else
				xstrITSSBOB01151cs += " + ";		// just show
			#endif
			xstrITSSBOB01151cs += xsholes[i];												
			xstrITSSBOB01151cs += ")"; 
		}
														
		xstrITSSBOB01151cs += "):xrotITSSBOB01151Circle";
		
		#ifdef _ITSSB_SHOW_WITHBOXDIV
		xstrITSSBOB01151cs += " * (xbboxdiv2:xtrandiv2)";
		#endif 

		TGeoVolume *xvolITSSBOB01151 = new TGeoVolume("ITSSBOB01151x", 
			new TGeoCompositeShape("xcompITSSB0B01151cs", xstrITSSBOB01151cs), 
			mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 
		xvolITSSBOB01151->SetFillColor(8);
		xvolITSSBOB01151->SetLineColor(8);
		ewhVol->AddNode(xvolITSSBOB01151, 1, ITSSBOB01151Offset);
	}
	
	if (fLayerNumber==6)
	{ // 0116x : Layer 6 handle

		// 01161 is the base ( tube )
		// 01162 is the base
		// 01163x are the handles
		// 01164 is the tube for div to make a hole in each base 1113x
		// 01165 is a tube for div to cut unwanted outside the tube
		// to make sure we can welded 01163x with (01161&01162), set 01165Radius to have more a bit value than 01162Radius 



		// ------------------------ 01161 ---------------------
							
		Double_t ITSSBOB01161Length1 =  5.0 *fgkmm; // fixed parameter
		Double_t ITSSBOB01161Length2 =  8.0 *fgkmm; // fixed parameter
		Double_t ITSSBOB01161Length3 = 16.0 *fgkmm; // fixed parameter

		Double_t ITSSBOB01161Radius = 383.0 *fgkmm;
		
		Double_t ITSSBOB01161Thick1 =  64.0 *fgkmm;
		Double_t ITSSBOB01161Thick2 =  66.0 *fgkmm;
		Double_t ITSSBOB01161Thick3 =  39.3 *fgkmm;
		Double_t ITSSBOB01161Thick4 =  53.0 *fgkmm;
		
		// offset is compared with volume's shape of section1
		TGeoTranslation *ITSSBOB01161Offset = new TGeoTranslation(0.0, 0.0, 0.0 + ITSSBOB01Length/2 - 0.0*fgkmm);

		Double_t ITSSBOB01161Z0 = 0.00; // The center of object is not the same as the center of space.
		Double_t ITSSBOB01161Z1 = -(ITSSBOB01161Length1);
		Double_t ITSSBOB01161Z2 = -(ITSSBOB01161Length1);
		Double_t ITSSBOB01161Z3 = -(ITSSBOB01161Length1 + ITSSBOB01161Length2);
		Double_t ITSSBOB01161Z4 = -(ITSSBOB01161Length1 + ITSSBOB01161Length2);
		Double_t ITSSBOB01161Z5 = -(ITSSBOB01161Length1 + ITSSBOB01161Length2 + ITSSBOB01161Length3);

		Double_t ITSSBOB01161RMin0 = ITSSBOB01161Radius + ITSSBOB01161Thick1, ITSSBOB01161RMax0 = ITSSBOB01161Radius + ITSSBOB01161Thick2;
		Double_t ITSSBOB01161RMin1 = ITSSBOB01161Radius + ITSSBOB01161Thick1, ITSSBOB01161RMax1 = ITSSBOB01161Radius + ITSSBOB01161Thick2;
		Double_t ITSSBOB01161RMin2 = ITSSBOB01161Radius, ITSSBOB01161RMax2 = ITSSBOB01161Radius + ITSSBOB01161Thick2;
		Double_t ITSSBOB01161RMin3 = ITSSBOB01161Radius, ITSSBOB01161RMax3 = ITSSBOB01161Radius + ITSSBOB01161Thick2;
		Double_t ITSSBOB01161RMin4 = ITSSBOB01161Radius + ITSSBOB01161Thick3, ITSSBOB01161RMax4 = ITSSBOB01161Radius + ITSSBOB01161Thick4;
		Double_t ITSSBOB01161RMin5 = ITSSBOB01161Radius + ITSSBOB01161Thick3, ITSSBOB01161RMax5 = ITSSBOB01161Radius + ITSSBOB01161Thick4;
												
		TGeoPcon *xpconeITSSBOB01161 = new TGeoPcon("xpconeITSSBOB01161", 0, 360, 6);
		xpconeITSSBOB01161->DefineSection(0, ITSSBOB01161Z0, ITSSBOB01161RMin0, ITSSBOB01161RMax0);
		xpconeITSSBOB01161->DefineSection(1, ITSSBOB01161Z1, ITSSBOB01161RMin1, ITSSBOB01161RMax1);
		xpconeITSSBOB01161->DefineSection(2, ITSSBOB01161Z2, ITSSBOB01161RMin2, ITSSBOB01161RMax2);
		xpconeITSSBOB01161->DefineSection(3, ITSSBOB01161Z3, ITSSBOB01161RMin3, ITSSBOB01161RMax3);
		xpconeITSSBOB01161->DefineSection(4, ITSSBOB01161Z4, ITSSBOB01161RMin4, ITSSBOB01161RMax4);
		xpconeITSSBOB01161->DefineSection(5, ITSSBOB01161Z5, ITSSBOB01161RMin5, ITSSBOB01161RMax5);

		// --------------- CREATE HOLES ---------------------

		
		TGeoRotation *xrotITSSBOB01161Circle = new TGeoRotation("xrotITSSBOB01161Circle", 0.0, 0.0, 0.0); // this rotate BEFORE div
		xrotITSSBOB01161Circle->RegisterYourself();

		vector<TString> xsholes = CreateDBOuterHoles("ITSSBIB01161"
			, 421.48 *fgkmm	// faces radius
			,  (10.0/2) *fgkmm  // holes radius
			,  (15.0*2.5) *fgkmm	// over approximate, for div
			, 0.0 // not rotate
			, - 7.0*fgkmm - ((ITSSBOB01161Length1+ITSSBOB01161Length2+ITSSBOB01161Length3)/2));

		// --------------- CREATE HANDLES -------------------
		
		vector<TString> xshandles = CreateDBOuterHandles("ITSSBOB01161"
			//, 421.48 *fgkmm	// faces radius
			, 418.48 *fgkmm	// faces radius, JUST FOR SHOW
			,  14.00 *fgkmm  // length, the same as base length
			,  52.00 *fgkmm  // width, overvalue, FOR SHOW : JUST REAL VALUE
			,   3.00 *fgkmm  // thick, overvalue
			, 0.0 // not rotate
			, + (14.0*fgkmm)/2 - (ITSSBOB01161Length1+ITSSBOB01161Length2+ITSSBOB01161Length3));		

		// ----------------------  combine all 0111x objects
		TString xstrITSSBOB01161cs("(");
		xstrITSSBOB01161cs += "xpconeITSSBOB01161";
	
		for (Int_t i=0; i<fNStaves; i++) // ITSSBOB011613
		{ // *** IF PLACE CODE 114 BELOW 113, ITS POSITION MIGHT BE UNDETERMINED ***
			xstrITSSBOB01161cs += " + (";
			xstrITSSBOB01161cs += xshandles[i];
			
			#ifdef _ITSSB_SHOW_COMBISHAPE
				xstrITSSBOB01161cs += " - ";		// make a hole
			#else
				xstrITSSBOB01161cs += " + ";		// just show
			#endif
			xstrITSSBOB01161cs += xsholes[i];												
			xstrITSSBOB01161cs += ")"; 
		}
					
		xstrITSSBOB01161cs += "):xrotITSSBOB01161Circle";
		
		#ifdef _ITSSB_SHOW_WITHBOXDIV
		xstrITSSBOB01161cs += " * (xbboxdiv2:xtrandiv2)";
		#endif 


		TGeoVolume *xvolITSSBOB01161 = new TGeoVolume("ITSSBOB01161x", 
			new TGeoCompositeShape("xcompITSSB0B01161cs", xstrITSSBOB01161cs), 
			mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 
		xvolITSSBOB01161->SetFillColor(3);
		xvolITSSBOB01161->SetLineColor(3);
		ewhVol->AddNode(xvolITSSBOB01161, 1, ITSSBOB01161Offset);

		

	}
	

	*/
	return ewhVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateOuterBEWheelA(const TGeoManager *mgr){
	//
	//
	//

	
	TGeoVolume* ewhVol = new TGeoVolumeAssembly(TString("EndWheelA") += fLayerNumber);

	/*
	const Double_t ITSSBOB01SWSkinThick = 0.60 *fgkmm;
	const Double_t ITSSBOB01SWCoreThick = 8.80 *fgkmm;
	const Double_t ITSSBOB01SWThick = ITSSBOB01SWCoreThick + 2*ITSSBOB01SWSkinThick;

	// GLOBAL VALUE FOR ALL SUB OBJECTS
	const Double_t ITSSBOB01CutOffset = 0.20;
	const Double_t ITSSBOB01Length = 1740.00 *fgkmm;
	const Double_t ITSSBOB01Radius =  449.00 *fgkmm; 
	const Double_t ITSSBOB01MinRadius = 2.20;
	const Double_t ITSSBOB01MaxRadius = ITSSBOB01Radius + ITSSBOB01SWThick;

	// box to div
	TGeoBBox *xbboxdiv = new TGeoBBox("xbboxdiv2", ITSSBOB01MaxRadius+1.0, ITSSBOB01MaxRadius+1.0, (ITSSBOB01Length/2)+1.0); // extend edge about 1.0
	TGeoTranslation *xtrandiv = new TGeoTranslation("xtrandiv2", 0, 0-xbboxdiv->GetDY()-ITSSBOB01CutOffset, 0);
	xtrandiv->RegisterYourself();			
	
	

	if (fLayerNumber==3)
	{ // 0123x : Layer 3 handle

		// 01231 is the base ( tube )
		// 01232 is the base
		// 01233x are the handles
		// 01234 is the tube for div to make a hole in each base 1113x
		// 01235 is a tube for div to cut unwanted outside the tube
		// to make sure we can welded 01233x with (01231&01232), set 01235Radius to have more a bit value than 01232Radius 
		
		
		Double_t ITSSBOB01231Length1 = 14.0 *fgkmm; // fixed parameter

		Double_t ITSSBOB01231Radius = 212.0 *fgkmm;
		
		Double_t ITSSBOB01231Thick1 =   7.5 *fgkmm;

		
		// offset is compared with volume's shape of section1
		TGeoTranslation *ITSSBOB01231Offset = new TGeoTranslation(0.0, 0.0, 0.0 - 372.0*fgkmm);

		Double_t ITSSBOB01231Z0 = 0.00; // The center of object is not the same as the center of space.
		Double_t ITSSBOB01231Z1 = -(ITSSBOB01231Length1);

		Double_t ITSSBOB01231RMin0 = ITSSBOB01231Radius, ITSSBOB01231RMax0 = ITSSBOB01231Radius + ITSSBOB01231Thick1;
		Double_t ITSSBOB01231RMin1 = ITSSBOB01231Radius, ITSSBOB01231RMax1 = ITSSBOB01231Radius + ITSSBOB01231Thick1;

		TGeoPcon *xpconeITSSBOB01231 = new TGeoPcon("xpconeITSSBOB01231", 0, 360, 2);
		xpconeITSSBOB01231->DefineSection(0, ITSSBOB01231Z0, ITSSBOB01231RMin0, ITSSBOB01231RMax0);
		xpconeITSSBOB01231->DefineSection(1, ITSSBOB01231Z1, ITSSBOB01231RMin1, ITSSBOB01231RMax1);


		// --------------- CREATE HOLES ---------------------

		
		TGeoRotation *xrotITSSBOB01231Circle = new TGeoRotation("xrotITSSBOB01231Circle", 0.0, 0.0, 0.0); // this rotate BEFORE div
		xrotITSSBOB01231Circle->RegisterYourself();

		vector<TString> xsholes = CreateDBOuterHoles("ITSSBIB01231"
			, 210.5 *fgkmm	// faces radius
			,  (10.0/2) *fgkmm  // holes radius
			,  (15.0*2.5) *fgkmm	// over approximate, for div
			, 0.0 // not rotate
			, - 7.0*fgkmm);

		// --------------- CREATE HANDLES -------------------
		
		vector<TString> xshandles = CreateDBOuterHandles("ITSSBOB01231"
			//, 210.5 *fgkmm	// faces radius
			, 207.5 *fgkmm	// faces radius, JUST FOR SHOW
			,  14.0 *fgkmm  // length, the same as base length
			,  50.5 *fgkmm  // width, overvalue, FOR SHOW : JUST REAL VALUE
			,   3.0 *fgkmm  // thick, overvalue
			, 0.0 // not rotate
			, - (14.0*fgkmm)/2);		

		// -----------------------  combine all 0111x objects
		TString xstrITSSBOB01231cs("(");
		xstrITSSBOB01231cs += "xpconeITSSBOB01231";
		
		for (Int_t i=0; i<fNStaves; i++) // ITSSBOB012313
		{ // *** IF PLACE CODE 114 BELOW 113, ITS POSITION MIGHT BE UNDETERMINED ***
			xstrITSSBOB01231cs += " + (";
			xstrITSSBOB01231cs += xshandles[i];
			
			#ifdef _ITSSB_SHOW_COMBISHAPE
				xstrITSSBOB01231cs += " - ";		// make a hole
			#else
				xstrITSSBOB01231cs += " + ";		// just show
			#endif
			xstrITSSBOB01231cs += xsholes[i];												
			xstrITSSBOB01231cs += ")"; 
		}
														
		xstrITSSBOB01231cs += "):xrotITSSBOB01231Circle";
		
		#ifdef _ITSSB_SHOW_WITHBOXDIV
		xstrITSSBOB01231cs += " * (xbboxdiv2:xtrandiv2)";
		#endif 
		
		TGeoVolume *xvolITSSBOB01231 = new TGeoVolume("ITSSBOB01231x", new TGeoCompositeShape("xcompITSSB0B01231cs", xstrITSSBOB01231cs), mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 
		xvolITSSBOB01231->SetFillColor(9);
		xvolITSSBOB01231->SetLineColor(9);
		ewhVol->AddNode(xvolITSSBOB01231, 1, ITSSBOB01231Offset);
		
		
		// ------------------------ 01232-1 ---------------------
							
		Double_t ITSSBOB012321Length1 = 15.0 *fgkmm; // fixed parameter
		Double_t ITSSBOB012321Length2 = 10.0 *fgkmm; // fixed parameter

		Double_t ITSSBOB012321Radius = 180.0 *fgkmm; // old one is 180.5, fixing by overlaps to be 180.0
		
		Double_t ITSSBOB012321Thick1 =   2.0 *fgkmm;
		Double_t ITSSBOB012321Thick2 =  40.0 *fgkmm;

		// offset is compared with volume's shape of section1
		TGeoTranslation *ITSSBOB012321Offset = new TGeoTranslation(0.0, 0.0, 0.0 - 482.0*fgkmm);

		Double_t ITSSBOB012321Z0 = 0.00; // The center of object is not the same as the center of space.
		Double_t ITSSBOB012321Z1 = -(ITSSBOB012321Length1);
		Double_t ITSSBOB012321Z2 = -(ITSSBOB012321Length1);
		Double_t ITSSBOB012321Z3 = -(ITSSBOB012321Length1 + ITSSBOB012321Length2);

		Double_t ITSSBOB012321RMin0 = ITSSBOB012321Radius + ITSSBOB012321Thick2 - ITSSBOB012321Thick1, ITSSBOB012321RMax0 = ITSSBOB012321Radius + ITSSBOB012321Thick2;
		Double_t ITSSBOB012321RMin1 = ITSSBOB012321Radius + ITSSBOB012321Thick2 - ITSSBOB012321Thick1, ITSSBOB012321RMax1 = ITSSBOB012321Radius + ITSSBOB012321Thick2;
		Double_t ITSSBOB012321RMin2 = ITSSBOB012321Radius, ITSSBOB012321RMax2 = ITSSBOB012321Radius + ITSSBOB012321Thick2;
		Double_t ITSSBOB012321RMin3 = ITSSBOB012321Radius, ITSSBOB012321RMax3 = ITSSBOB012321Radius + ITSSBOB012321Thick2;										

		TGeoPcon *xpconeITSSBOB012321 = new TGeoPcon("xpconeITSSBOB012321", 0, 360, 4);
		xpconeITSSBOB012321->DefineSection(0, ITSSBOB012321Z0, ITSSBOB012321RMin0, ITSSBOB012321RMax0);
		xpconeITSSBOB012321->DefineSection(1, ITSSBOB012321Z1, ITSSBOB012321RMin1, ITSSBOB012321RMax1);
		xpconeITSSBOB012321->DefineSection(2, ITSSBOB012321Z2, ITSSBOB012321RMin2, ITSSBOB012321RMax2);
		xpconeITSSBOB012321->DefineSection(3, ITSSBOB012321Z3, ITSSBOB012321RMin3, ITSSBOB012321RMax3);


		//TGeoCompositeShape *xcompITSSBOB0114 = new TGeoCompositeShape("xcompITSSBOB0114cs", xstrITSSBOB0114cs);
		//TGeoVolume *xvolITSSBOB0114 = new TGeoVolume("ITSSBOB0114x", xcompITSSBOB0114, med); 	
		#ifdef _ITSSB_SHOW_WITHBOXDIV
		TGeoVolume *xvolITSSBOB012321 = new TGeoVolume("ITSSBOB012321", new TGeoCompositeShape("xpconeITSSBOB012321*(xbboxdiv2:xtrandiv2)"), mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 					
		#else		
		TGeoVolume *xvolITSSBOB012321 = new TGeoVolume("ITSSBOB012321", xpconeITSSBOB012321, mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 
		#endif
		xvolITSSBOB012321->SetFillColor(9);
		xvolITSSBOB012321->SetLineColor(9);
		ewhVol->AddNode(xvolITSSBOB012321, 1, ITSSBOB012321Offset);

		// ------------------------ 01232-2 ---------------------
							
		Double_t ITSSBOB012322Length1 = 15.0 *fgkmm; // fixed parameter
		Double_t ITSSBOB012322Length2 =  3.0 *fgkmm; // fixed parameter

		Double_t ITSSBOB012322Radius = 216.5 *fgkmm;
		
		Double_t ITSSBOB012322Thick1 =   3.0 *fgkmm;
		Double_t ITSSBOB012322Thick2 =  25.0 *fgkmm;
		
		// offset is compared with volume's shape of section1
		TGeoTranslation *ITSSBOB012322Offset = new TGeoTranslation(0.0, 0.0, 0.0 - 509.0*fgkmm);

		Double_t ITSSBOB012322Z0 = 0.00; // The center of object is not the same as the center of space.
		Double_t ITSSBOB012322Z1 = -(ITSSBOB012322Length1);
		Double_t ITSSBOB012322Z2 = -(ITSSBOB012322Length1);
		Double_t ITSSBOB012322Z3 = -(ITSSBOB012322Length1 + ITSSBOB012322Length2);

		Double_t ITSSBOB012322RMin0 = ITSSBOB012322Radius, ITSSBOB012322RMax0 = ITSSBOB012322Radius + ITSSBOB012322Thick1;
		Double_t ITSSBOB012322RMin1 = ITSSBOB012322Radius, ITSSBOB012322RMax1 = ITSSBOB012322Radius + ITSSBOB012322Thick1;
		Double_t ITSSBOB012322RMin2 = ITSSBOB012322Radius, ITSSBOB012322RMax2 = ITSSBOB012322Radius + ITSSBOB012322Thick2;
		Double_t ITSSBOB012322RMin3 = ITSSBOB012322Radius, ITSSBOB012322RMax3 = ITSSBOB012322Radius + ITSSBOB012322Thick2;										

		TGeoPcon *xpconeITSSBOB012322 = new TGeoPcon("xpconeITSSBOB012322", 0, 360, 4);
		xpconeITSSBOB012322->DefineSection(0, ITSSBOB012322Z0, ITSSBOB012322RMin0, ITSSBOB012322RMax0);
		xpconeITSSBOB012322->DefineSection(1, ITSSBOB012322Z1, ITSSBOB012322RMin1, ITSSBOB012322RMax1);
		xpconeITSSBOB012322->DefineSection(2, ITSSBOB012322Z2, ITSSBOB012322RMin2, ITSSBOB012322RMax2);
		xpconeITSSBOB012322->DefineSection(3, ITSSBOB012322Z3, ITSSBOB012322RMin3, ITSSBOB012322RMax3);


		//TGeoCompositeShape *xcompITSSBOB0114 = new TGeoCompositeShape("xcompITSSBOB0114cs", xstrITSSBOB0114cs);
		//TGeoVolume *xvolITSSBOB0114 = new TGeoVolume("ITSSBOB0114x", xcompITSSBOB0114, med); 	
		#ifdef _ITSSB_SHOW_WITHBOXDIV
		TGeoVolume *xvolITSSBOB012322 = new TGeoVolume("ITSSBOB012322", new TGeoCompositeShape("xpconeITSSBOB012322*(xbboxdiv2:xtrandiv2)"), mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 					
		#else		
		TGeoVolume *xvolITSSBOB012322 = new TGeoVolume("ITSSBOB012322", xpconeITSSBOB012322, mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 
		#endif
		xvolITSSBOB012322->SetFillColor(9);
		xvolITSSBOB012322->SetLineColor(9);
		ewhVol->AddNode(xvolITSSBOB012322, 1, ITSSBOB012322Offset);
		
		//------------------------------------ Base Plate ---------------------
		
		const Double_t ITSSBOB01130SWSkinThick = 0.5 *fgkmm;
		const Double_t ITSSBOB01130SWCoreThick = 5.0 *fgkmm;
		const Double_t ITSSBOB01130SWThick = ITSSBOB01130SWCoreThick + 2*ITSSBOB01130SWSkinThick;
		
		const Double_t ITSSBOB01130Length = 152.0 *fgkmm;
		const Double_t ITSSBOB01130Radius = 219.5 *fgkmm;
		const Double_t ITSSBOB01130MaxRadius = ITSSBOB01130Radius + ITSSBOB01130SWThick;
		
		// offset is compared with volume's shape of section1
		TGeoTranslation *ITSSBOB01130Offset = new TGeoTranslation(0.0, 0.0, 0.0 - ITSSBOB01130Length/2 - 372.0*fgkmm);
		
		TGeoMedium *medITSSBOB01130SWSkin = mgr->GetMedium("ITS_KAPTON(POLYCH2)$");
		TGeoMedium *medITSSBOB01130SWCore = mgr->GetMedium("ITS_KAPTON(POLYCH2)$");

		
		// 011301 : Inside
		new TGeoTube("xtubeITSSBOB011301", ITSSBOB01130Radius, ITSSBOB01130Radius + ITSSBOB01130SWSkinThick, ITSSBOB01130Length/2);
		TGeoVolume *xvolITSSBOB011301 = new TGeoVolume("ITSSBOB011301", new TGeoCompositeShape("xtubeITSSBOB011301*(xbboxdiv2:xtrandiv2)"), medITSSBOB01130SWSkin); 			
		xvolITSSBOB011301->SetFillColor(9);
		xvolITSSBOB011301->SetLineColor(9);
		ewhVol->AddNode(xvolITSSBOB011301, 1, ITSSBOB01130Offset);
		
		// 011302 : Foam
		new TGeoTube("xtubeITSSBOB011302", ITSSBOB01130Radius + ITSSBOB01130SWSkinThick, ITSSBOB01130Radius + ITSSBOB01130SWSkinThick + ITSSBOB01130SWCoreThick, ITSSBOB01130Length/2);
		TGeoVolume *xvolITSSBOB011302 = new TGeoVolume("ITSSBOB011302", new TGeoCompositeShape("xtubeITSSBOB011302*(xbboxdiv2:xtrandiv2)"), medITSSBOB01130SWCore); 			
		xvolITSSBOB011302->SetFillColor(5);
		xvolITSSBOB011302->SetLineColor(5);
		ewhVol->AddNode(xvolITSSBOB011302, 1, ITSSBOB01130Offset);
		
		// 011303 : Outside
		new TGeoTube("xtubeITSSBOB011303", ITSSBOB01130Radius + ITSSBOB01130SWSkinThick + ITSSBOB01130SWCoreThick, ITSSBOB01130Radius + 2*ITSSBOB01130SWSkinThick + ITSSBOB01130SWCoreThick, ITSSBOB01130Length/2);
		TGeoVolume *xvolITSSBOB011303 = new TGeoVolume("ITSSBOB011303", new TGeoCompositeShape("xtubeITSSBOB011303*(xbboxdiv2:xtrandiv2)"), medITSSBOB01130SWSkin); 			
		xvolITSSBOB011303->SetFillColor(9);
		xvolITSSBOB011303->SetLineColor(9);
		ewhVol->AddNode(xvolITSSBOB011303, 1, ITSSBOB01130Offset);
		
	}
	
	if (fLayerNumber==4)
	{ // 0124x : Layer 4 handle

		// 01241 is the base ( tube )
		// 01242 is the base
		// 01243x are the handles
		// 01244 is the tube for div to make a hole in each base 1113x
		// 01245 is a tube for div to cut unwanted outside the tube
		// to make sure we can welded 01243x with (01241&01242), set 01245Radius to have more a bit value than 01242Radius 
		
		
		Double_t ITSSBOB01241Length1 = 14.0 *fgkmm; // fixed parameter

		Double_t ITSSBOB01241Radius = 265.0 *fgkmm;
		
		Double_t ITSSBOB01241Thick1 =   9.0 *fgkmm;
		
		// offset is compared with volume's shape of section1
		TGeoTranslation *ITSSBOB01241Offset = new TGeoTranslation(0.0, 0.0, 0.0 - 372.0*fgkmm);

		Double_t ITSSBOB01241Z0 = 0.00; // The center of object is not the same as the center of space.
		Double_t ITSSBOB01241Z1 = -(ITSSBOB01241Length1);

		Double_t ITSSBOB01241RMin0 = ITSSBOB01241Radius, ITSSBOB01241RMax0 = ITSSBOB01241Radius + ITSSBOB01241Thick1;
		Double_t ITSSBOB01241RMin1 = ITSSBOB01241Radius, ITSSBOB01241RMax1 = ITSSBOB01241Radius + ITSSBOB01241Thick1;

		TGeoPcon *xpconeITSSBOB01241 = new TGeoPcon("xpconeITSSBOB01241", 0, 360, 2);
		xpconeITSSBOB01241->DefineSection(0, ITSSBOB01241Z0, ITSSBOB01241RMin0, ITSSBOB01241RMax0);
		xpconeITSSBOB01241->DefineSection(1, ITSSBOB01241Z1, ITSSBOB01241RMin1, ITSSBOB01241RMax1);

		// --------------- CREATE HOLES ---------------------

		
		TGeoRotation *xrotITSSBOB01241Circle = new TGeoRotation("xrotITSSBOB01241Circle", 0.0, 0.0, 0.0); // this rotate BEFORE div
		xrotITSSBOB01241Circle->RegisterYourself();

		vector<TString> xsholes = CreateDBOuterHoles("ITSSBIB01241"
			, 263.5 *fgkmm	// faces radius
			,  (10.0/2) *fgkmm  // holes radius
			,  (15.0*2.5) *fgkmm	// over approximate, for div
			, 0.0 // not rotate
			, - 7.0*fgkmm);

		// --------------- CREATE HANDLES -------------------
		
		vector<TString> xshandles = CreateDBOuterHandles("ITSSBOB01241"
			//, 263.5 *fgkmm	// faces radius
			, 260.5 *fgkmm	// faces radius, JUST FOR SHOW
			,  14.0 *fgkmm  // length, the same as base length
			,  56.5 *fgkmm  // width, overvalue, FOR SHOW : JUST REAL VALUE
			,   3.0 *fgkmm  // thick, overvalue
			, 0.0 // not rotate
			, - (14.0*fgkmm)/2);		

		// ---------------------- combine all 0111x objects
		TString xstrITSSBOB01241cs("(");
		xstrITSSBOB01241cs += "xpconeITSSBOB01241";
		
		for (Int_t i=0; i<fNStaves; i++) // ITSSBOB012413
		{ // *** IF PLACE CODE 114 BELOW 113, ITS POSITION MIGHT BE UNDETERMINED ***
			xstrITSSBOB01241cs += " + (";
			xstrITSSBOB01241cs += xshandles[i];
			
			#ifdef _ITSSB_SHOW_COMBISHAPE
				xstrITSSBOB01241cs += " - ";		// make a hole
			#else
				xstrITSSBOB01241cs += " + ";		// just show
			#endif
			xstrITSSBOB01241cs += xsholes[i];												
			xstrITSSBOB01241cs += ")"; 
		}
														
		xstrITSSBOB01241cs += "):xrotITSSBOB01241Circle";
		
		#ifdef _ITSSB_SHOW_WITHBOXDIV
		xstrITSSBOB01241cs += " * (xbboxdiv2:xtrandiv2)";
		#endif 
		
		TGeoVolume *xvolITSSBOB01241 = new TGeoVolume("ITSSBOB01241x", new TGeoCompositeShape("xcompITSSB0B01241cs", xstrITSSBOB01241cs), mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 
		xvolITSSBOB01241->SetFillColor(7);
		xvolITSSBOB01241->SetLineColor(7);
		ewhVol->AddNode(xvolITSSBOB01241, 1, ITSSBOB01241Offset);
		
		
		// ------------------------ 01242-1 ---------------------
							
		Double_t ITSSBOB012421Length1 = 15.0 *fgkmm; // fixed parameter
		Double_t ITSSBOB012421Length2 = 10.0 *fgkmm; // fixed parameter

		Double_t ITSSBOB012421Radius = 225.5 *fgkmm;
		
		Double_t ITSSBOB012421Thick1 =   2.0 *fgkmm;
		Double_t ITSSBOB012421Thick2 =  48.5 *fgkmm;

		// offset is compared with volume's shape of section1
		TGeoTranslation *ITSSBOB012421Offset = new TGeoTranslation(0.0, 0.0, 0.0 - 498.0*fgkmm);

		Double_t ITSSBOB012421Z0 = 0.00; // The center of object is not the same as the center of space.
		Double_t ITSSBOB012421Z1 = -(ITSSBOB012421Length1);
		Double_t ITSSBOB012421Z2 = -(ITSSBOB012421Length1);
		Double_t ITSSBOB012421Z3 = -(ITSSBOB012421Length1 + ITSSBOB012421Length2);

		Double_t ITSSBOB012421RMin0 = ITSSBOB012421Radius + ITSSBOB012421Thick2 - ITSSBOB012421Thick1, ITSSBOB012421RMax0 = ITSSBOB012421Radius + ITSSBOB012421Thick2;
		Double_t ITSSBOB012421RMin1 = ITSSBOB012421Radius + ITSSBOB012421Thick2 - ITSSBOB012421Thick1, ITSSBOB012421RMax1 = ITSSBOB012421Radius + ITSSBOB012421Thick2;
		Double_t ITSSBOB012421RMin2 = ITSSBOB012421Radius, ITSSBOB012421RMax2 = ITSSBOB012421Radius + ITSSBOB012421Thick2;
		Double_t ITSSBOB012421RMin3 = ITSSBOB012421Radius, ITSSBOB012421RMax3 = ITSSBOB012421Radius + ITSSBOB012421Thick2;										

		TGeoPcon *xpconeITSSBOB012421 = new TGeoPcon("xpconeITSSBOB012421", 0, 360, 4);
		xpconeITSSBOB012421->DefineSection(0, ITSSBOB012421Z0, ITSSBOB012421RMin0, ITSSBOB012421RMax0);
		xpconeITSSBOB012421->DefineSection(1, ITSSBOB012421Z1, ITSSBOB012421RMin1, ITSSBOB012421RMax1);
		xpconeITSSBOB012421->DefineSection(2, ITSSBOB012421Z2, ITSSBOB012421RMin2, ITSSBOB012421RMax2);
		xpconeITSSBOB012421->DefineSection(3, ITSSBOB012421Z3, ITSSBOB012421RMin3, ITSSBOB012421RMax3);


		//TGeoCompositeShape *xcompITSSBOB0114 = new TGeoCompositeShape("xcompITSSBOB0114cs", xstrITSSBOB0114cs);
		//TGeoVolume *xvolITSSBOB0114 = new TGeoVolume("ITSSBOB0114x", xcompITSSBOB0114, med); 	
		#ifdef _ITSSB_SHOW_WITHBOXDIV
		TGeoVolume *xvolITSSBOB012421 = new TGeoVolume("ITSSBOB012421", new TGeoCompositeShape("xpconeITSSBOB012421*(xbboxdiv2:xtrandiv2)"), mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 					
		#else		
		TGeoVolume *xvolITSSBOB012421 = new TGeoVolume("ITSSBOB012421", xpconeITSSBOB012421, mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 
		#endif
		xvolITSSBOB012421->SetFillColor(7);
		xvolITSSBOB012421->SetLineColor(7);
		ewhVol->AddNode(xvolITSSBOB012421, 1, ITSSBOB012421Offset);

		// ------------------------ 01242-2 ---------------------
							
		Double_t ITSSBOB012422Length1 = 15.0 *fgkmm; // fixed parameter
		Double_t ITSSBOB012422Length2 =  3.0 *fgkmm; // fixed parameter

		Double_t ITSSBOB012422Radius = 320.0 *fgkmm;
		
		Double_t ITSSBOB012422Thick1 =   3.0 *fgkmm;
		Double_t ITSSBOB012422Thick2 =  25.0 *fgkmm;
		
		// offset is compared with volume's shape of section1
		TGeoTranslation *ITSSBOB012422Offset = new TGeoTranslation(0.0, 0.0, 0.0 - 783.0*fgkmm);

		Double_t ITSSBOB012422Z0 = 0.00; // The center of object is not the same as the center of space.
		Double_t ITSSBOB012422Z1 = -(ITSSBOB012422Length1);
		Double_t ITSSBOB012422Z2 = -(ITSSBOB012422Length1);
		Double_t ITSSBOB012422Z3 = -(ITSSBOB012422Length1 + ITSSBOB012422Length2);

		Double_t ITSSBOB012422RMin0 = ITSSBOB012422Radius, ITSSBOB012422RMax0 = ITSSBOB012422Radius + ITSSBOB012422Thick1;
		Double_t ITSSBOB012422RMin1 = ITSSBOB012422Radius, ITSSBOB012422RMax1 = ITSSBOB012422Radius + ITSSBOB012422Thick1;
		Double_t ITSSBOB012422RMin2 = ITSSBOB012422Radius, ITSSBOB012422RMax2 = ITSSBOB012422Radius + ITSSBOB012422Thick2;
		Double_t ITSSBOB012422RMin3 = ITSSBOB012422Radius, ITSSBOB012422RMax3 = ITSSBOB012422Radius + ITSSBOB012422Thick2;										

		TGeoPcon *xpconeITSSBOB012422 = new TGeoPcon("xpconeITSSBOB012422", 0, 360, 4);
		xpconeITSSBOB012422->DefineSection(0, ITSSBOB012422Z0, ITSSBOB012422RMin0, ITSSBOB012422RMax0);
		xpconeITSSBOB012422->DefineSection(1, ITSSBOB012422Z1, ITSSBOB012422RMin1, ITSSBOB012422RMax1);
		xpconeITSSBOB012422->DefineSection(2, ITSSBOB012422Z2, ITSSBOB012422RMin2, ITSSBOB012422RMax2);
		xpconeITSSBOB012422->DefineSection(3, ITSSBOB012422Z3, ITSSBOB012422RMin3, ITSSBOB012422RMax3);


		//TGeoCompositeShape *xcompITSSBOB0114 = new TGeoCompositeShape("xcompITSSBOB0114cs", xstrITSSBOB0114cs);
		//TGeoVolume *xvolITSSBOB0114 = new TGeoVolume("ITSSBOB0114x", xcompITSSBOB0114, med); 	
		#ifdef _ITSSB_SHOW_WITHBOXDIV
		TGeoVolume *xvolITSSBOB012422 = new TGeoVolume("ITSSBOB012422", new TGeoCompositeShape("xpconeITSSBOB012422*(xbboxdiv2:xtrandiv2)"), mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 					
		#else		
		TGeoVolume *xvolITSSBOB012422 = new TGeoVolume("ITSSBOB012422", xpconeITSSBOB012422, mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 
		#endif
		xvolITSSBOB012422->SetFillColor(7);
		xvolITSSBOB012422->SetLineColor(7);
		ewhVol->AddNode(xvolITSSBOB012422, 1, ITSSBOB012422Offset);
		
		//ewhVol->AddNode(CreateSTube45("OuterLTube", 274.0*fgkmm, 323.0*fgkmm, 426.0*fgkmm, 170*fgkmm, 5.0*fgkmm + 2*0.5*fgkmm, 0, 7), 1, new TGeoCombiTrans(0.0, 0.0, 0.0 - (174.0*fgkcm)/2 + (72.0*fgkmm), new TGeoRotation("",0,180,0)), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));

	}
	
	if (fLayerNumber==5)
	{ // 0125x : Layer 5 handle

		// 01251 is the base ( tube )
		// 01252 is the base
		// 01253x are the handles
		// 01254 is the tube for div to make a hole in each base 1113x
		// 01255 is a tube for div to cut unwanted outside the tube
		// to make sure we can welded 01253x with (01251&01252), set 01255Radius to have more a bit value than 01252Radius 
		
		
		Double_t ITSSBOB01251Length1 = 14.0 *fgkmm; // fixed parameter

		Double_t ITSSBOB01251Radius = 371.0 *fgkmm;
		
		Double_t ITSSBOB01251Thick1 =   6.0 *fgkmm;

		
		// offset is compared with volume's shape of section1
		TGeoTranslation *ITSSBOB01251Offset = new TGeoTranslation(0.0, 0.0, 0.0 - 685.0*fgkmm);

		Double_t ITSSBOB01251Z0 = 0.00; // The center of object is not the same as the center of space.
		Double_t ITSSBOB01251Z1 = -(ITSSBOB01251Length1);

		Double_t ITSSBOB01251RMin0 = ITSSBOB01251Radius, ITSSBOB01251RMax0 = ITSSBOB01251Radius + ITSSBOB01251Thick1;
		Double_t ITSSBOB01251RMin1 = ITSSBOB01251Radius, ITSSBOB01251RMax1 = ITSSBOB01251Radius + ITSSBOB01251Thick1;

		TGeoPcon *xpconeITSSBOB01251 = new TGeoPcon("xpconeITSSBOB01251", 0, 360, 2);
		xpconeITSSBOB01251->DefineSection(0, ITSSBOB01251Z0, ITSSBOB01251RMin0, ITSSBOB01251RMax0);
		xpconeITSSBOB01251->DefineSection(1, ITSSBOB01251Z1, ITSSBOB01251RMin1, ITSSBOB01251RMax1);

		// --------------- CREATE HOLES ---------------------

		
		TGeoRotation *xrotITSSBOB01251Circle = new TGeoRotation("xrotITSSBOB01251Circle", 0.0, 0.0, 0.0); // this rotate BEFORE div
		xrotITSSBOB01251Circle->RegisterYourself();

		vector<TString> xsholes = CreateDBOuterHoles("ITSSBIB01251"
			, 370.0 *fgkmm	// faces radius
			,  (10.0/2) *fgkmm  // holes radius
			,  (15.0*2.5) *fgkmm	// over approximate, for div
			, 0.0 // not rotate
			, - 7.0*fgkmm);

		// --------------- CREATE HANDLES -------------------
		
		vector<TString> xshandles = CreateDBOuterHandles("ITSSBOB01251"
			//, 370.0 *fgkmm	// faces radius
			, 367.0 *fgkmm	// faces radius, JUST FOR SHOW
			,  14.0 *fgkmm  // length, the same as base length
			,  54.8 *fgkmm  // width, overvalue, FOR SHOW : JUST REAL VALUE
			,   3.0 *fgkmm  // thick, overvalue
			, 0.0 // not rotate
			, - (14.0*fgkmm)/2 );		


		// combine all 0111x objects
		TString xstrITSSBOB01251cs("(");
		xstrITSSBOB01251cs += "xpconeITSSBOB01251";
		
		for (Int_t i=0; i<fNStaves; i++) // ITSSBOB012513
		{ // *** IF PLACE CODE 114 BELOW 113, ITS POSITION MIGHT BE UNDETERMINED ***
			xstrITSSBOB01251cs += " + (";
			xstrITSSBOB01251cs += xshandles[i];
			
			#ifdef _ITSSB_SHOW_COMBISHAPE
				xstrITSSBOB01251cs += " - ";		// make a hole
			#else
				xstrITSSBOB01251cs += " + ";		// just show
			#endif
			xstrITSSBOB01251cs += xsholes[i];												
			xstrITSSBOB01251cs += ")"; 
		}
														
		xstrITSSBOB01251cs += "):xrotITSSBOB01251Circle";
		
		#ifdef _ITSSB_SHOW_WITHBOXDIV
		xstrITSSBOB01251cs += " * (xbboxdiv2:xtrandiv2)";
		#endif 
		
		TGeoVolume *xvolITSSBOB01251 = new TGeoVolume("ITSSBOB01251x", new TGeoCompositeShape("xcompITSSB0B01251cs", xstrITSSBOB01251cs), mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 
		xvolITSSBOB01251->SetFillColor(8);
		xvolITSSBOB01251->SetLineColor(8);
		ewhVol->AddNode(xvolITSSBOB01251, 1, ITSSBOB01251Offset);
		
		
		// ------------------------ 01252-1 ---------------------
							
		Double_t ITSSBOB012521Length1 = 15.0 *fgkmm; // fixed parameter
		Double_t ITSSBOB012521Length2 = 10.0 *fgkmm; // fixed parameter

		Double_t ITSSBOB012521Radius = 329.0 *fgkmm;
		
		Double_t ITSSBOB012521Thick1 =   2.0 *fgkmm;
		Double_t ITSSBOB012521Thick2 =  48.0 *fgkmm;

		// offset is compared with volume's shape of section1
		TGeoTranslation *ITSSBOB012521Offset = new TGeoTranslation(0.0, 0.0, 0.0 - 772.0*fgkmm);

		Double_t ITSSBOB012521Z0 = 0.00; // The center of object is not the same as the center of space.
		Double_t ITSSBOB012521Z1 = -(ITSSBOB012521Length1);
		Double_t ITSSBOB012521Z2 = -(ITSSBOB012521Length1);
		Double_t ITSSBOB012521Z3 = -(ITSSBOB012521Length1 + ITSSBOB012521Length2);

		Double_t ITSSBOB012521RMin0 = ITSSBOB012521Radius + ITSSBOB012521Thick2 - ITSSBOB012521Thick1, ITSSBOB012521RMax0 = ITSSBOB012521Radius + ITSSBOB012521Thick2;
		Double_t ITSSBOB012521RMin1 = ITSSBOB012521Radius + ITSSBOB012521Thick2 - ITSSBOB012521Thick1, ITSSBOB012521RMax1 = ITSSBOB012521Radius + ITSSBOB012521Thick2;
		Double_t ITSSBOB012521RMin2 = ITSSBOB012521Radius, ITSSBOB012521RMax2 = ITSSBOB012521Radius + ITSSBOB012521Thick2;
		Double_t ITSSBOB012521RMin3 = ITSSBOB012521Radius, ITSSBOB012521RMax3 = ITSSBOB012521Radius + ITSSBOB012521Thick2;										

		TGeoPcon *xpconeITSSBOB012521 = new TGeoPcon("xpconeITSSBOB012521", 0, 360, 4);
		xpconeITSSBOB012521->DefineSection(0, ITSSBOB012521Z0, ITSSBOB012521RMin0, ITSSBOB012521RMax0);
		xpconeITSSBOB012521->DefineSection(1, ITSSBOB012521Z1, ITSSBOB012521RMin1, ITSSBOB012521RMax1);
		xpconeITSSBOB012521->DefineSection(2, ITSSBOB012521Z2, ITSSBOB012521RMin2, ITSSBOB012521RMax2);
		xpconeITSSBOB012521->DefineSection(3, ITSSBOB012521Z3, ITSSBOB012521RMin3, ITSSBOB012521RMax3);


		//TGeoCompositeShape *xcompITSSBOB0114 = new TGeoCompositeShape("xcompITSSBOB0114cs", xstrITSSBOB0114cs);
		//TGeoVolume *xvolITSSBOB0114 = new TGeoVolume("ITSSBOB0114x", xcompITSSBOB0114, med); 	
		#ifdef _ITSSB_SHOW_WITHBOXDIV
		TGeoVolume *xvolITSSBOB012521 = new TGeoVolume("ITSSBOB012521", new TGeoCompositeShape("xpconeITSSBOB012521*(xbboxdiv2:xtrandiv2)"), mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 					
		#else		
		TGeoVolume *xvolITSSBOB012521 = new TGeoVolume("ITSSBOB012521", xpconeITSSBOB012521, mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 
		#endif
		xvolITSSBOB012521->SetFillColor(8);
		xvolITSSBOB012521->SetLineColor(8);
		ewhVol->AddNode(xvolITSSBOB012521, 1, ITSSBOB012521Offset);

		// ------------------------ 01252-2 ---------------------
							
		Double_t ITSSBOB012522Length1 = 15.0 *fgkmm; // fixed parameter
		Double_t ITSSBOB012522Length2 =  3.0 *fgkmm; // fixed parameter

		Double_t ITSSBOB012522Radius = 374.0 *fgkmm;
		
		Double_t ITSSBOB012522Thick1 =   3.0 *fgkmm;
		Double_t ITSSBOB012522Thick2 =  25.0 *fgkmm;
		
		// offset is compared with volume's shape of section1
		TGeoTranslation *ITSSBOB012522Offset = new TGeoTranslation(0.0, 0.0, 0.0 - 800.0*fgkmm);

		Double_t ITSSBOB012522Z0 = 0.00; // The center of object is not the same as the center of space.
		Double_t ITSSBOB012522Z1 = -(ITSSBOB012522Length1);
		Double_t ITSSBOB012522Z2 = -(ITSSBOB012522Length1);
		Double_t ITSSBOB012522Z3 = -(ITSSBOB012522Length1 + ITSSBOB012522Length2);

		Double_t ITSSBOB012522RMin0 = ITSSBOB012522Radius, ITSSBOB012522RMax0 = ITSSBOB012522Radius + ITSSBOB012522Thick1;
		Double_t ITSSBOB012522RMin1 = ITSSBOB012522Radius, ITSSBOB012522RMax1 = ITSSBOB012522Radius + ITSSBOB012522Thick1;
		Double_t ITSSBOB012522RMin2 = ITSSBOB012522Radius, ITSSBOB012522RMax2 = ITSSBOB012522Radius + ITSSBOB012522Thick2;
		Double_t ITSSBOB012522RMin3 = ITSSBOB012522Radius, ITSSBOB012522RMax3 = ITSSBOB012522Radius + ITSSBOB012522Thick2;										

		TGeoPcon *xpconeITSSBOB012522 = new TGeoPcon("xpconeITSSBOB012522", 0, 360, 4);
		xpconeITSSBOB012522->DefineSection(0, ITSSBOB012522Z0, ITSSBOB012522RMin0, ITSSBOB012522RMax0);
		xpconeITSSBOB012522->DefineSection(1, ITSSBOB012522Z1, ITSSBOB012522RMin1, ITSSBOB012522RMax1);
		xpconeITSSBOB012522->DefineSection(2, ITSSBOB012522Z2, ITSSBOB012522RMin2, ITSSBOB012522RMax2);
		xpconeITSSBOB012522->DefineSection(3, ITSSBOB012522Z3, ITSSBOB012522RMin3, ITSSBOB012522RMax3);


		//TGeoCompositeShape *xcompITSSBOB0114 = new TGeoCompositeShape("xcompITSSBOB0114cs", xstrITSSBOB0114cs);
		//TGeoVolume *xvolITSSBOB0114 = new TGeoVolume("ITSSBOB0114x", xcompITSSBOB0114, med); 	
		#ifdef _ITSSB_SHOW_WITHBOXDIV
		TGeoVolume *xvolITSSBOB012522 = new TGeoVolume("ITSSBOB012522", new TGeoCompositeShape("xpconeITSSBOB012522*(xbboxdiv2:xtrandiv2)"), mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 					
		#else		
		TGeoVolume *xvolITSSBOB012522 = new TGeoVolume("ITSSBOB012522", xpconeITSSBOB012522, mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 
		#endif
		xvolITSSBOB012522->SetFillColor(8);
		xvolITSSBOB012522->SetLineColor(8);
		ewhVol->AddNode(xvolITSSBOB012522, 1, ITSSBOB012522Offset);
		
		//------------------------------------ Base Plate ---------------------
		
		const Double_t ITSSBOB01250SWSkinThick = 0.5 *fgkmm;
		const Double_t ITSSBOB01250SWCoreThick = 5.0 *fgkmm;
		const Double_t ITSSBOB01250SWThick = ITSSBOB01250SWCoreThick + 2*ITSSBOB01250SWSkinThick;
		
		const Double_t ITSSBOB01250Length = 130.0 *fgkmm;
		const Double_t ITSSBOB01250Radius = 377.0 *fgkmm;
		const Double_t ITSSBOB01250MaxRadius = ITSSBOB01250Radius + ITSSBOB01250SWThick;
		
		// offset is compared with volume's shape of section1
		TGeoTranslation *ITSSBOB01250Offset = new TGeoTranslation(0.0, 0.0, 0.0 - ITSSBOB01250Length/2 - 685.0*fgkmm);
		
		TGeoMedium *medITSSBOB01250SWSkin = mgr->GetMedium("ITS_KAPTON(POLYCH2)$");
		TGeoMedium *medITSSBOB01250SWCore = mgr->GetMedium("ITS_KAPTON(POLYCH2)$");

		
		// 012501 : Inside
		new TGeoTube("xtubeITSSBOB012501", ITSSBOB01250Radius, ITSSBOB01250Radius + ITSSBOB01250SWSkinThick, ITSSBOB01250Length/2);
		TGeoVolume *xvolITSSBOB012501 = new TGeoVolume("ITSSBOB012501", new TGeoCompositeShape("xtubeITSSBOB012501*(xbboxdiv2:xtrandiv2)"), medITSSBOB01250SWSkin); 			
		xvolITSSBOB012501->SetFillColor(8);
		xvolITSSBOB012501->SetLineColor(8);
		ewhVol->AddNode(xvolITSSBOB012501, 1, ITSSBOB01250Offset);
		
		// 012502 : Foam
		new TGeoTube("xtubeITSSBOB012502", ITSSBOB01250Radius + ITSSBOB01250SWSkinThick, ITSSBOB01250Radius + ITSSBOB01250SWSkinThick + ITSSBOB01250SWCoreThick, ITSSBOB01250Length/2);
		TGeoVolume *xvolITSSBOB012502 = new TGeoVolume("ITSSBOB012502", new TGeoCompositeShape("xtubeITSSBOB012502*(xbboxdiv2:xtrandiv2)"), medITSSBOB01250SWCore); 			
		xvolITSSBOB012502->SetFillColor(5);
		xvolITSSBOB012502->SetLineColor(5);
		ewhVol->AddNode(xvolITSSBOB012502, 1, ITSSBOB01250Offset);
		
		// 012503 : Outside
		new TGeoTube("xtubeITSSBOB012503", ITSSBOB01250Radius + ITSSBOB01250SWSkinThick + ITSSBOB01250SWCoreThick, ITSSBOB01250Radius + 2*ITSSBOB01250SWSkinThick + ITSSBOB01250SWCoreThick, ITSSBOB01250Length/2);
		TGeoVolume *xvolITSSBOB012503 = new TGeoVolume("ITSSBOB012503", new TGeoCompositeShape("xtubeITSSBOB012503*(xbboxdiv2:xtrandiv2)"), medITSSBOB01250SWSkin); 			
		xvolITSSBOB012503->SetFillColor(8);
		xvolITSSBOB012503->SetLineColor(8);
		ewhVol->AddNode(xvolITSSBOB012503, 1, ITSSBOB01250Offset);
		

	}
	
	if (fLayerNumber==6)
	{ // 0126x : Layer 6 handle

		// 01261 is the base ( tube )
		// 01262 is the base
		// 01263x are the handles
		// 01264 is the tube for div to make a hole in each base 1113x
		// 01265 is a tube for div to cut unwanted outside the tube
		// to make sure we can welded 01263x with (01261&01262), set 01265Radius to have more a bit value than 01262Radius 

		// ------------------------ 01261 ---------------------
							
		Double_t ITSSBOB01261Length1 = 14.0 *fgkmm; // fixed parameter
		Double_t ITSSBOB01261Length2 =  8.0 *fgkmm; // fixed parameter
		Double_t ITSSBOB01261Length3 =  5.0 *fgkmm; // fixed parameter

		Double_t ITSSBOB01261Radius = 422.3 *fgkmm;
		
		Double_t ITSSBOB01261Thick1 =  13.7 *fgkmm;
		Double_t ITSSBOB01261Thick2 =  26.7 *fgkmm;
		Double_t ITSSBOB01261Thick3 =   2.0 *fgkmm;
		
		
		// offset is compared with volume's shape of section1
		TGeoTranslation *ITSSBOB01261Offset = new TGeoTranslation(0.0, 0.0, 0.0 - 685.0*fgkmm);

		Double_t ITSSBOB01261Z0 = 0.00; // The center of object is not the same as the center of space.
		Double_t ITSSBOB01261Z1 = -(ITSSBOB01261Length1);
		Double_t ITSSBOB01261Z2 = -(ITSSBOB01261Length1);
		Double_t ITSSBOB01261Z3 = -(ITSSBOB01261Length1 + ITSSBOB01261Length2);
		Double_t ITSSBOB01261Z4 = -(ITSSBOB01261Length1 + ITSSBOB01261Length2);
		Double_t ITSSBOB01261Z5 = -(ITSSBOB01261Length1 + ITSSBOB01261Length2 + ITSSBOB01261Length3);

		Double_t ITSSBOB01261RMin0 = ITSSBOB01261Radius, ITSSBOB01261RMax0 = ITSSBOB01261Radius + ITSSBOB01261Thick1;
		Double_t ITSSBOB01261RMin1 = ITSSBOB01261Radius, ITSSBOB01261RMax1 = ITSSBOB01261Radius + ITSSBOB01261Thick1;
		Double_t ITSSBOB01261RMin2 = ITSSBOB01261Radius, ITSSBOB01261RMax2 = ITSSBOB01261Radius + ITSSBOB01261Thick2;
		Double_t ITSSBOB01261RMin3 = ITSSBOB01261Radius, ITSSBOB01261RMax3 = ITSSBOB01261Radius + ITSSBOB01261Thick2;
		Double_t ITSSBOB01261RMin4 = ITSSBOB01261Radius + ITSSBOB01261Thick2 - ITSSBOB01261Thick3, ITSSBOB01261RMax4 = ITSSBOB01261Radius + ITSSBOB01261Thick2;
		Double_t ITSSBOB01261RMin5 = ITSSBOB01261Radius + ITSSBOB01261Thick2 - ITSSBOB01261Thick3, ITSSBOB01261RMax5 = ITSSBOB01261Radius + ITSSBOB01261Thick2;
												
		TGeoPcon *xpconeITSSBOB01261 = new TGeoPcon("xpconeITSSBOB01261", 0, 360, 6);
		xpconeITSSBOB01261->DefineSection(0, ITSSBOB01261Z0, ITSSBOB01261RMin0, ITSSBOB01261RMax0);
		xpconeITSSBOB01261->DefineSection(1, ITSSBOB01261Z1, ITSSBOB01261RMin1, ITSSBOB01261RMax1);
		xpconeITSSBOB01261->DefineSection(2, ITSSBOB01261Z2, ITSSBOB01261RMin2, ITSSBOB01261RMax2);
		xpconeITSSBOB01261->DefineSection(3, ITSSBOB01261Z3, ITSSBOB01261RMin3, ITSSBOB01261RMax3);
		xpconeITSSBOB01261->DefineSection(4, ITSSBOB01261Z4, ITSSBOB01261RMin4, ITSSBOB01261RMax4);
		xpconeITSSBOB01261->DefineSection(5, ITSSBOB01261Z5, ITSSBOB01261RMin5, ITSSBOB01261RMax5);

		// --------------- CREATE HOLES ---------------------

		
		TGeoRotation *xrotITSSBOB01261Circle = new TGeoRotation("xrotITSSBOB01261Circle", 0.0, 0.0, 0.0); // this rotate BEFORE div
		xrotITSSBOB01261Circle->RegisterYourself();

		vector<TString> xsholes = CreateDBOuterHoles("ITSSBIB01261"
			, 421.48 *fgkmm	// faces radius
			,  (10.0/2) *fgkmm  // holes radius
			,  (15.0*2.5) *fgkmm	// over approximate, for div
			, 0.0 // not rotate
			, - 7.0*fgkmm);

		// --------------- CREATE HANDLES -------------------
		
		vector<TString> xshandles = CreateDBOuterHandles("ITSSBOB01261"
			//, 421.48 *fgkmm	// faces radius
			, 418.48 *fgkmm	// faces radius, JUST FOR SHOW
			,  22.00 *fgkmm  // length, the same as base length
			,  52.00 *fgkmm  // width, overvalue, FOR SHOW : JUST REAL VALUE
			,   3.00 *fgkmm  // thick, overvalue
			, 0.0 // not rotate
			, - (22.0*fgkmm)/2 );		

		// ------------------------ combine all 0111x objects
		TString xstrITSSBOB01261cs("(");
		xstrITSSBOB01261cs += "xpconeITSSBOB01261";
		
		for (Int_t i=0; i<fNStaves; i++) // ITSSBOB012613
		{ // *** IF PLACE CODE 114 BELOW 113, ITS POSITION MIGHT BE UNDETERMINED ***
			xstrITSSBOB01261cs += " + (";
			xstrITSSBOB01261cs += xshandles[i];
			
			#ifdef _ITSSB_SHOW_COMBISHAPE
				xstrITSSBOB01261cs += " - ";		// make a hole
			#else
				xstrITSSBOB01261cs += " + ";		// just show
			#endif
			xstrITSSBOB01261cs += xsholes[i];												
			xstrITSSBOB01261cs += ")"; 
		}
														
		xstrITSSBOB01261cs += "):xrotITSSBOB01261Circle";
		
		#ifdef _ITSSB_SHOW_WITHBOXDIV
		xstrITSSBOB01261cs += " * (xbboxdiv2:xtrandiv2)";
		#endif 
		
		TGeoVolume *xvolITSSBOB01261 = new TGeoVolume("ITSSBOB01261x", new TGeoCompositeShape("xcompITSSB0B01261cs", xstrITSSBOB01261cs), mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 
		xvolITSSBOB01261->SetFillColor(3);
		xvolITSSBOB01261->SetLineColor(3);
		ewhVol->AddNode(xvolITSSBOB01261, 1, ITSSBOB01261Offset);

		// ------------------------ 01262 ---------------------
							
		Double_t ITSSBOB01262Length1 = 15.0 *fgkmm; // fixed parameter
		Double_t ITSSBOB01262Length2 = 10.0 *fgkmm; // fixed parameter

		Double_t ITSSBOB01262Radius = 383.0 *fgkmm;
		
		Double_t ITSSBOB01262Thick1 =   2.0 *fgkmm;
		Double_t ITSSBOB01262Thick2 =  66.0 *fgkmm;

		// offset is compared with volume's shape of section1
		TGeoTranslation *ITSSBOB01262Offset = new TGeoTranslation(0.0, 0.0, 0.0 - 789.0*fgkmm);

		Double_t ITSSBOB01262Z0 = 0.00; // The center of object is not the same as the center of space.
		Double_t ITSSBOB01262Z1 = -(ITSSBOB01262Length1);
		Double_t ITSSBOB01262Z2 = -(ITSSBOB01262Length1);
		Double_t ITSSBOB01262Z3 = -(ITSSBOB01262Length1 + ITSSBOB01262Length2);

		Double_t ITSSBOB01262RMin0 = ITSSBOB01262Radius + ITSSBOB01262Thick2 - ITSSBOB01262Thick1, ITSSBOB01262RMax0 = ITSSBOB01262Radius + ITSSBOB01262Thick2;
		Double_t ITSSBOB01262RMin1 = ITSSBOB01262Radius + ITSSBOB01262Thick2 - ITSSBOB01262Thick1, ITSSBOB01262RMax1 = ITSSBOB01262Radius + ITSSBOB01262Thick2;
		Double_t ITSSBOB01262RMin2 = ITSSBOB01262Radius, ITSSBOB01262RMax2 = ITSSBOB01262Radius + ITSSBOB01262Thick2;
		Double_t ITSSBOB01262RMin3 = ITSSBOB01262Radius, ITSSBOB01262RMax3 = ITSSBOB01262Radius + ITSSBOB01262Thick2;										

		TGeoPcon *xpconeITSSBOB01262 = new TGeoPcon("xpconeITSSBOB01262", 0, 360, 4);
		xpconeITSSBOB01262->DefineSection(0, ITSSBOB01262Z0, ITSSBOB01262RMin0, ITSSBOB01262RMax0);
		xpconeITSSBOB01262->DefineSection(1, ITSSBOB01262Z1, ITSSBOB01262RMin1, ITSSBOB01262RMax1);
		xpconeITSSBOB01262->DefineSection(2, ITSSBOB01262Z2, ITSSBOB01262RMin2, ITSSBOB01262RMax2);
		xpconeITSSBOB01262->DefineSection(3, ITSSBOB01262Z3, ITSSBOB01262RMin3, ITSSBOB01262RMax3);


		//TGeoCompositeShape *xcompITSSBOB0114 = new TGeoCompositeShape("xcompITSSBOB0114cs", xstrITSSBOB0114cs);
		//TGeoVolume *xvolITSSBOB0114 = new TGeoVolume("ITSSBOB0114x", xcompITSSBOB0114, med); 	
		#ifdef _ITSSB_SHOW_WITHBOXDIV
		TGeoVolume *xvolITSSBOB01262 = new TGeoVolume("ITSSBOB01262", new TGeoCompositeShape("xpconeITSSBOB01262*(xbboxdiv2:xtrandiv2)"), mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 					
		#else		
		TGeoVolume *xvolITSSBOB01262 = new TGeoVolume("ITSSBOB01262", xpconeITSSBOB01262, mgr->GetMedium("ITS_KAPTON(POLYCH2)$")); 
		#endif
		xvolITSSBOB01262->SetFillColor(3);
		xvolITSSBOB01262->SetLineColor(3);
		ewhVol->AddNode(xvolITSSBOB01262, 1, ITSSBOB01262Offset);

	}
	
	*/
	return ewhVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateInnerBHandleA(const char* fkName, 
												const Double_t fkHndPlateRadius, 
												const Double_t fkContactRadius, 
												const TGeoMedium* med){
	//
	//
	//
	

	const Double_t fkHoleDia =   5.0*fgkmm;
	const Double_t fkHoleRadius = fkHoleDia/2;
	const Double_t fkHoleFromEdge = 4.0*fgkmm;		// center-of-hole to inside-edge
	const Double_t fkHoleLength = 10.0*fgkmm;		// approximate

	const Double_t fkHndLength = 13.0*fgkmm;
	const Double_t fkHndWidth = 10.0*fgkmm;
	const Double_t fkHndThick = 5.0*fgkmm;
	const Double_t fkHndPlateThick = 0.6*fgkmm;

	const Double_t fkHoleFromEdge1 = fkHoleFromEdge - fkHoleRadius;
	const Double_t fkHoleFromEdge2 = fkHoleFromEdge + fkHoleRadius;	

	const Double_t fkHndWidthAngle = 360/fNStaves;
	
	const Double_t fkEachRotateAngle = 15.0;

	TString xStr, xTmp1, xTmp2;
	TGeoShape* shp;
	
	// ----------------- Create the plate -----------------
	xTmp1 = TString("shpHndPlate")+=fLayerNumber;
	shp = new TGeoPcon(xTmp1, 0, fkHndWidthAngle, 2);
	((TGeoPcon*)shp)->DefineSection(0, 0.0					, fkHndPlateRadius	, fkHndPlateRadius + fkHndPlateThick);
	((TGeoPcon*)shp)->DefineSection(1, (-1)*(fkHndLength)	, fkHndPlateRadius	, fkHndPlateRadius + fkHndPlateThick);
	xTmp2 = TString("xtrHndPlate")+=fLayerNumber;
	(new TGeoCombiTrans(xTmp2, 0, 0, 0, new TGeoRotation("", 0, 0, -fkHndWidthAngle/2)))->RegisterYourself();
	xStr = xTmp1 + ":" + xTmp2;
	
	// ----------------- Create the base -----------------
	xTmp1 = TString("shpHndBase")+=fLayerNumber;
	shp = new TGeoBBox(xTmp1, fkHndWidth/2, fkHndThick/2, fkHndLength/2);
	xTmp2 = TString("xtrHndBase")+=fLayerNumber;
	(new TGeoCombiTrans(xTmp2, fkContactRadius+fkHndThick/2, 0, -fkHndLength/2, new TGeoRotation("", 90, 0, fkEachRotateAngle)))->RegisterYourself();

	#ifdef _ITSSB_SHOW_COMBISHAPE
	{
		// ----------------- Cut the base -----------------
		TString xStrAA, xTmp1AA, xTmp2AA;
		TGeoShape* shpAA;
		
		xTmp1AA = TString("shpHndCutPlate")+=fLayerNumber;
		shpAA = new TGeoPcon(xTmp1AA, 0, fkHndWidthAngle, 2);
		((TGeoPcon*)shpAA)->DefineSection(0, 0.0					, fkHndPlateRadius + fkHndPlateThick/2	, fkHndPlateRadius + fkHndWidth);
		((TGeoPcon*)shpAA)->DefineSection(1, -(fkHndLength+1*fgkmm)	, fkHndPlateRadius + fkHndPlateThick/2	, fkHndPlateRadius + fkHndWidth);
		xTmp2AA = TString("xtrHndCutPlate") += fLayerNumber;
		(new TGeoCombiTrans(xTmp2AA, 0, 0, 0.5*fgkmm, new TGeoRotation("", 0, 0, -fkHndWidthAngle/2)))->RegisterYourself();		
		xStrAA = TString("(") + xTmp1 + ":" + xTmp2 + " - " + xTmp1AA + ":" + xTmp2AA + ")";
		xStr += TString(" + ") + xStrAA;
	}
	#else
		xStr += TString(" + ") + xTmp1 + ":" + xTmp2;
	#endif
	
	// ----------------- Make a hole -----------------
	xTmp1 = TString("shpHndHole")+=fLayerNumber;
	shp = new TGeoTube(xTmp1, 0, fkHoleRadius, fkHoleLength);
	xTmp2 = TString("xtrHndHole")+=fLayerNumber;
	(new TGeoCombiTrans(xTmp2, fkContactRadius, 0, -fkHoleFromEdge, new TGeoRotation("", 90+fkEachRotateAngle, 90, 0)))->RegisterYourself();
	#ifdef _ITSSB_SHOW_COMBISHAPE
		xStr = TString("(")+ xStr + TString(")") + TString(" - ") + xTmp1 + ":" + xTmp2;
	#else
		xStr += TString(" + ") + xTmp1 + ":" + xTmp2;
	#endif
	
	// ----------------- Show indicator -----------------	
	#ifdef _ITSSB_SHOW_INDICATOR
		xTmp1 = TString("shpHndInd")+=fLayerNumber;
		shp = new TGeoTube(xTmp1, 0, 0.1*fgkmm, fkHoleLength);
		xTmp2 = TString("xtrHndInd")+=fLayerNumber;
		(new TGeoCombiTrans(xTmp2, fkContactRadius, 0, -fkHoleFromEdge, new TGeoRotation("", 0, 0, 0)))->RegisterYourself();
		xStr += TString(" + ") + xTmp1 + ":" + xTmp2;
	#endif
	
	return new TGeoVolume((TString("HandleAEach")+=fLayerNumber)+=TString(fkName), new TGeoCompositeShape("", xStr), med); 			
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateInnerBHandleC(const TGeoMedium* med){
	//
	//
	//

	TGeoVolume* ewhVol = new TGeoVolumeAssembly(TString("HandleA") += fLayerNumber);
	
	Double_t fHoleRadius =   5.*fgkmm;
	Double_t fHoleFromEdge = 4.*fgkmm;		// the center of the hole

	Double_t fHoleFromEdge1 = fHoleFromEdge - fHoleRadius/2;	// the border of the hole
	Double_t fHoleFromEdge2 = fHoleFromEdge + fHoleRadius/2;	// the border of the hole	

	if (fLayerNumber==0)
	{
		//ewhVol->AddNode(CreateTube("EndWheelObj[1]", 28.1*fgkmm, 5.7, 0.4*fgkmm, kFALSE, 6, 180, med), 1, 
			//new TGeoCombiTrans(0., 0., 0.-(400./2-57.-28.)*fgkmm-fHoleFromEdge2, new TGeoRotation("",0.,0.,180.)));

		//ewhVol->AddNode(CreateTube("EndWheelObj[2]", 28.1*fgkmm, 1.3, 0.4*fgkmm, kFALSE, 6, 360./fNStaves, med), 1, 
			//new TGeoCombiTrans(0., 0., 0.-(400./2-57.-28.)*fgkmm, new TGeoRotation("",0.,0.,180.)));
		/*
			TString xStr;
			xStr += CreateInnerBHandle("EndWheel", 2.576, 1.30, 1.10, 0.45, 15.0, 0);
			TGeoCompositeShape *xComp = new TGeoCompositeShape("xcom0", xStr);
			TGeoVolume *xVol = new TGeoVolume("XCOM0", xComp, med); 			
			xVol->SetFillColor(6);
			xVol->SetLineColor(6);
		ewhVol->AddNode(xVol, 1, new TGeoCombiTrans(0., 0., 0.-(400./2-57.-28.)*fgkmm, new TGeoRotation("",0.,0.,0.)));
*/
	}
	
	if (fLayerNumber==1)
	{
		//ewhVol->AddNode(CreateSTube45("EndWheelObj[1]", 28.5*fgkmm, (36.5-0.4)*fgkmm, 52.*fgkmm-fHoleFromEdge2, 3.4*fgkmm, 0.4*fgkmm, kFALSE, 5, med), 1, 
			//new TGeoCombiTrans(0., 0., 0.-(400./2-52.-33.)*fgkmm-fHoleFromEdge2, new TGeoRotation("",0.,0.,180.)));

		//ewhVol->AddNode(CreateTube("EndWheelObj[2]", (36.5-0.4)*fgkmm, fHoleFromEdge1, 0.4*fgkmm, kFALSE, 5, 180, med), 1, 
			//new TGeoCombiTrans(0., 0., 0.-(400./2-52.-33.)*fgkmm, new TGeoRotation("",0.,0.,180.)));
		
	}

	
	if (fLayerNumber==2)
	{
		//ewhVol->AddNode(CreateSTube45("EndWheelObj[1]", 36.5*fgkmm, (44.-0.4)*fgkmm, 36.*fgkmm-fHoleFromEdge2, 4.*fgkmm, 0.4*fgkmm, kFALSE, 7, med), 1, 
			//new TGeoCombiTrans(0., 0., 0.-(400./2-36.-49.)*fgkmm-fHoleFromEdge2, new TGeoRotation("",0.,0.,180.)));

		//ewhVol->AddNode(CreateTube("EndWheelObj[2]", (44.-0.4)*fgkmm, fHoleFromEdge1, 0.4*fgkmm, kFALSE, 7, 180, med), 1, 
			//new TGeoCombiTrans(0., 0., 0.-(400./2-36.-49.)*fgkmm, new TGeoRotation("",0.,0.,180.)));

		
		// Supporter (A-Side)
		//ewhVol->AddNode(CreateInnerBSupporterA(), 1, new TGeoTranslation(0.,0.,0.-(400./2-5.-65.)*fgkmm));
	}
	

	return ewhVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::Create_ALC_0334_xxA(const TGeoMedium* med){
	//
	//
	//
	const Double_t fkHoleRadius = fkALC_0334_HoleDia/2;

	const Double_t fkHoleFromEdge1 = fkALC_0334_HoleFromEdge - fkHoleRadius;
	const Double_t fkHoleFromEdge2 = fkALC_0334_HoleFromEdge + fkHoleRadius;		

	TString xStr, xTmp1, xTmp2;
	TGeoShape* shp;
	
	// ----------------- Create the plate -----------------
	//xTmp1 = TString("shpHndPlate")+=fLayerNumber;
	//shp = new TGeoPcon(xTmp1, 0, 360., 2);
		//((TGeoPcon*)shp)->DefineSection(0, 0.0						, fkALC_0334_PlateRadius[fLayerNumber]	, fkALC_0334_PlateRadius[fLayerNumber] + fkALC_0334_PlateThick);
		//((TGeoPcon*)shp)->DefineSection(1, -fkALC_0334_HndLength	, fkALC_0334_PlateRadius[fLayerNumber]	, fkALC_0334_PlateRadius[fLayerNumber] + fkALC_0334_PlateThick);
	//xTmp2 = TString("xtrHndPlate")+=fLayerNumber;
		//(new TGeoCombiTrans(xTmp2, 0, 0, fkALC_0334_HoleFromEdge, new TGeoRotation("", 0, 0, 0)))->RegisterYourself();
	//xStr = xTmp1 + ":" + xTmp2;

	// ----------------- Create the base -----------------
	xTmp1 = (TString("shpHndBaseA")+=fLayerNumber) + (TString("_")+=(int)fNHandleCreated);
		shp = new TGeoBBox(xTmp1, fkALC_0334_HndWidth/2, fkALC_0334_HndThick/2, fkALC_0334_HndLength/2);
	xTmp2 = (TString("xtrHndBaseA")+=fLayerNumber) + (TString("_")+=(int)fNHandleCreated);
		(new TGeoCombiTrans(xTmp2, -fkALC_0334_HndWidth/2+fkALC_0334_HoleFromEdge+fkALC_0334_ContactDX[fLayerNumber], +fkALC_0334_HndThick/2+fkALC_0334_ContactDY[fLayerNumber], +fkALC_0334_HndLength/2-fkALC_0334_HoleFromEdge, new TGeoRotation("", 0, 0, 0)))->RegisterYourself();

	#ifdef _ITSSB_SHOW_COMBISHAPE
		xStr = xTmp1 + ":" + xTmp2;							// without plate
		//xStr += TString(" + ") + xTmp1 + ":" + xTmp2;		// with plate
	#else
	{
		TString xStrAA, xTmp1AA, xTmp2AA;
		TGeoShape* shpAA;
		
		xTmp1AA = (TString("shpHndCutPlateA")+=fLayerNumber) + (TString("_")+=(int)fNHandleCreated);
		shpAA = new TGeoPcon(xTmp1AA, 0, 360., 2);
			((TGeoPcon*)shpAA)->DefineSection(0, 0.0							, fkALC_0334_PlateRadius[fLayerNumber] 	, fkALC_0334_PlateRadius[fLayerNumber] + fkALC_0334_HndThick);
			((TGeoPcon*)shpAA)->DefineSection(1, +fkALC_0334_HndLength+2*fgkmm	, fkALC_0334_PlateRadius[fLayerNumber] 	, fkALC_0334_PlateRadius[fLayerNumber] + fkALC_0334_HndThick);
		
		xTmp2AA = (TString("xtrHndCutPlateA") += fLayerNumber) + (TString("_")+=(int)fNHandleCreated);
			(new TGeoCombiTrans(xTmp2AA, 0, 0, -(fkALC_0334_HoleFromEdge+1*fgkmm), new TGeoRotation("", 0, 0, 0)))->RegisterYourself();		
		xStrAA = TString("(") + xTmp1 + ":" + xTmp2 + " - " + xTmp1AA + ":" + xTmp2AA + ")";
		xStr = xStrAA;										// without plate
		//xStr += TString(" + ") + xStrAA;					// with plate

	}
	#endif


	// ----------------- Make a hole -----------------
	xTmp1 = (TString("shpHndHoleA")+=fLayerNumber) + (TString("_")+=(int)fNHandleCreated);
		shp = new TGeoTube(xTmp1, 0, fkHoleRadius, fkALC_0334_HoleLength);
	xTmp2 = (TString("xtrHndHoleA")+=fLayerNumber) + (TString("_")+=(int)fNHandleCreated);
		(new TGeoCombiTrans(xTmp2, +fkALC_0334_ContactDX[fLayerNumber], +fkALC_0334_ContactDY[fLayerNumber], 0, new TGeoRotation("", 0, 90, 0)))->RegisterYourself();
	#ifdef _ITSSB_SHOW_COMBISHAPE
		xStr += TString(" + ") + xTmp1 + ":" + xTmp2;
	#else
		xStr = TString("(")+ xStr + TString(")") + TString(" - ") + xTmp1 + ":" + xTmp2;
	#endif
	
	return new TGeoVolume((TString("ALC_0334_0") += (int)fLayerNumber) + (TString("A_") += (int)fNHandleCreated++), new TGeoCompositeShape("", xStr), med); 			
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::Create_ALC_0336_xxA_ins(const TGeoMedium* med){

	TString xStr, xTmp1, xTmp2;
	TGeoShape* shp;
	
	// ----------------- Create the plate -----------------
	xTmp1 = (TString("shpHndPlateA")+=fLayerNumber) + (TString("_")+=(int)fNHandleCreated);
	shp = new TGeoPcon(xTmp1, 0, 180., 2);
		((TGeoPcon*)shp)->DefineSection(0, 0.0						, fkALC_0334_PlateRadius[fLayerNumber]	, fkALC_0334_PlateRadius[fLayerNumber] + fkALC_0334_PlateThick);
		((TGeoPcon*)shp)->DefineSection(1, +fkALC_0334_HndLength	, fkALC_0334_PlateRadius[fLayerNumber]	, fkALC_0334_PlateRadius[fLayerNumber] + fkALC_0334_PlateThick);
	//xTmp2 = (TString("xtrHndPlateA")+=fLayerNumber)  + (TString("_")+=(int)fNHandleCreated);
		//(new TGeoCombiTrans(xTmp2, 0, 0, fkALC_0334_HoleFromEdge, new TGeoRotation("", 0, 0, 0)))->RegisterYourself();
	//xStr = xTmp1 + ":" + xTmp2;
	
		
	//return new TGeoVolume((TString("ALC_0337_0") += (int)fLayerNumber) + TString("C"), new TGeoCompositeShape("", xStr), med); 			
	return new TGeoVolume((TString("ALC_0336_0") += (int)fLayerNumber) + TString("A"), shp, med);
}

//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::Create_ALC_0336_xxA_ous(const Int_t color, const TGeoMedium* med){

	TGeoVolume* ewhVol = new TGeoVolumeAssembly(TString("EndWheelA") += fLayerNumber);

	//const Double_t kOffsetZ = -18.0*fgkmm;
	//const Double_t fkPlateThick = 0.4*fgkmm;
	//const Double_t fkHndLength = 13.0*fgkmm; // also use in CreateInnerBHandleA.
		
	//const Double_t kHndPlateRadius[7]={27.9*fgkmm, 35.9*fgkmm, 43.4*fgkmm,0,0,0,0};
	//const Double_t kContactRadius[7]={25.9*fgkmm, 33.8*fgkmm, 41.1*fgkmm,0,0,0,0}; // approx : -2mm from Radius
		
	//const Double_t fkHndWidthAngle = 360/fNStaves;
	static const Int_t kNumOfConeSection = 12;
	static const Double_t kConeArray[(int)fgkNumberOfInnerLayers][kNumOfConeSection][3] =
	{		// Layer 0
		{{ 	  10.00*fgkmm,  27.90*fgkmm,  28.50*fgkmm },	// Section 0
		{ 	  38.00*fgkmm,  27.90*fgkmm,  28.50*fgkmm },	// Section 1
		{ 	  38.00*fgkmm,  27.90*fgkmm,  29.30*fgkmm },	// Section 2
		{ 	  43.00*fgkmm,  27.90*fgkmm,  29.30*fgkmm },	// Section 3
		{ 	  43.00*fgkmm,  28.50*fgkmm,  29.30*fgkmm },	// Section 4
		{ 	 118.96*fgkmm,  56.15*fgkmm,  57.00*fgkmm },	// Section 5
		{ 	 119.11*fgkmm,  56.20*fgkmm,  57.00*fgkmm },	// Section 6
		{ 	 121.00*fgkmm,  56.20*fgkmm,  57.00*fgkmm },	// Section 7
		{ 	 121.00*fgkmm,  56.20*fgkmm,  58.00*fgkmm },	// Section 8
		{ 	 127.00*fgkmm,  56.20*fgkmm,  58.00*fgkmm },	// Section 9
		{ 	 127.00*fgkmm,  57.00*fgkmm,  58.00*fgkmm },	// Section 10
		{ 	 187.00*fgkmm,  57.00*fgkmm,  58.00*fgkmm }}	// Section 11

		,	// Layer 1
		{{ 	  10.00*fgkmm,  35.90*fgkmm,  36.50*fgkmm },	// Section 0
		{ 	  35.00*fgkmm,  35.90*fgkmm,  36.50*fgkmm },	// Section 1
		{ 	  35.00*fgkmm,  35.90*fgkmm,  37.30*fgkmm },	// Section 2
		{ 	  40.00*fgkmm,  35.90*fgkmm,  37.30*fgkmm },	// Section 3
		{ 	  40.00*fgkmm,  36.50*fgkmm,  37.30*fgkmm },	// Section 4
		{ 	 125.87*fgkmm,  86.08*fgkmm,  87.00*fgkmm },	// Section 5
		{ 	 126.08*fgkmm,  86.20*fgkmm,  87.00*fgkmm },	// Section 6
		{ 	 128.00*fgkmm,  86.20*fgkmm,  87.00*fgkmm },	// Section 7
		{ 	 128.00*fgkmm,  86.20*fgkmm,  88.00*fgkmm },	// Section 8
		{ 	 133.00*fgkmm,  86.20*fgkmm,  88.00*fgkmm },	// Section 9
		{ 	 133.00*fgkmm,  87.00*fgkmm,  88.00*fgkmm },	// Section 10
		{ 	 180.00*fgkmm,  87.00*fgkmm,  88.00*fgkmm }}	// Section 11

		,	// Layer 2
		{{ 	  10.00*fgkmm,  43.90*fgkmm,  44.50*fgkmm },	// Section 0
		{ 	  38.00*fgkmm,  43.90*fgkmm,  44.50*fgkmm },	// Section 1
		{ 	  38.00*fgkmm,  43.90*fgkmm,  45.30*fgkmm },	// Section 2
		{ 	  43.00*fgkmm,  43.90*fgkmm,  45.30*fgkmm },	// Section 3
		{ 	  43.00*fgkmm,  44.50*fgkmm,  45.30*fgkmm },	// Section 4
		{ 	 118.96*fgkmm, 115.96*fgkmm, 117.00*fgkmm },	// Section 5
		{ 	 119.11*fgkmm, 116.20*fgkmm, 117.00*fgkmm },	// Section 6
		{ 	 121.00*fgkmm, 116.20*fgkmm, 117.00*fgkmm },	// Section 7
		{ 	 121.00*fgkmm, 116.20*fgkmm, 118.00*fgkmm },	// Section 8
		{ 	 127.00*fgkmm, 116.20*fgkmm, 118.00*fgkmm },	// Section 9
		{ 	 127.00*fgkmm, 117.00*fgkmm, 118.00*fgkmm },	// Section 10
		{ 	 173.00*fgkmm, 117.00*fgkmm, 118.00*fgkmm }}	// Section 11
	};
		
	static const Int_t kNumOfConeTubeSection = 2;
	static const Double_t kConeTubeArray[(int)fgkNumberOfInnerLayers][kNumOfConeTubeSection][3] =
	{		// Layer 0
		{{ 	 140.00*fgkmm,  28.00*fgkmm,  29.00*fgkmm },	// Section 0
		{ 	 182.00*fgkmm,  28.00*fgkmm,  29.00*fgkmm }}	// Section 1
		,	// Layer 1
		{{ 	 138.00*fgkmm,  58.00*fgkmm,  59.00*fgkmm },	// Section 0
		{ 	 180.00*fgkmm,  58.00*fgkmm,  59.00*fgkmm }}	// Section 1
		,	// Layer 2
		{{ 	 131.00*fgkmm,  88.00*fgkmm,  89.00*fgkmm },	// Section 0
		{ 	 173.00*fgkmm,  88.00*fgkmm,  89.00*fgkmm }}	// Section 1
	};
	
	static const Int_t kNumOfConeEdgeSection = 5;
	static const Double_t kConeEdgeArray[(int)fgkNumberOfInnerLayers][kNumOfConeEdgeSection][3] =
	{		// Layer 0
		{{ 	 176.00*fgkmm,  59.20*fgkmm,  59.50*fgkmm },	// Section 0
		{ 	 176.20*fgkmm,  59.00*fgkmm,  59.50*fgkmm },	// Section 1
		{ 	 181.00*fgkmm,  59.00*fgkmm,  59.50*fgkmm },	// Section 2
		{ 	 181.00*fgkmm,  58.00*fgkmm,  59.50*fgkmm },	// Section 3
		{ 	 187.00*fgkmm,  58.00*fgkmm,  59.50*fgkmm }}	// Section 4
		,	// Layer 1
		{{ 	 169.00*fgkmm,  89.20*fgkmm,  89.50*fgkmm },	// Section 0
		{ 	 169.20*fgkmm,  89.00*fgkmm,  89.50*fgkmm },	// Section 1
		{ 	 174.00*fgkmm,  89.00*fgkmm,  89.50*fgkmm },	// Section 2
		{ 	 174.00*fgkmm,  88.00*fgkmm,  89.50*fgkmm },	// Section 3
		{ 	 180.00*fgkmm,  88.00*fgkmm,  89.50*fgkmm }}	// Section 4
		,	// Layer 2
		{{ 	 162.00*fgkmm, 119.20*fgkmm, 119.50*fgkmm },	// Section 0
		{ 	 162.20*fgkmm, 119.00*fgkmm, 119.50*fgkmm },	// Section 1
		{ 	 167.00*fgkmm, 119.00*fgkmm, 119.50*fgkmm },	// Section 2
		{ 	 167.00*fgkmm, 118.00*fgkmm, 119.50*fgkmm },	// Section 3
		{ 	 173.00*fgkmm, 118.00*fgkmm, 119.50*fgkmm }}	// Section 4
	};
	
	
	static const Int_t kNumOfConeLSection = 4;
	static const Double_t kConeLArray[(int)fgkNumberOfInnerLayers][kNumOfConeLSection][3] =
	{		// Layer 0
		{{ 	  10.20*fgkmm,  26.50*fgkmm,  27.90*fgkmm },	// Section 0
		{ 	  10.60*fgkmm,  26.50*fgkmm,  27.90*fgkmm },	// Section 1
		{ 	  10.60*fgkmm,  27.50*fgkmm,  27.90*fgkmm },	// Section 2
		{ 	  13.20*fgkmm,  27.50*fgkmm,  27.90*fgkmm }}	// Section 3
		,	// Layer 1
		{{ 	  10.20*fgkmm,  34.50*fgkmm,  35.90*fgkmm },	// Section 0
		{ 	  10.60*fgkmm,  34.50*fgkmm,  35.90*fgkmm },	// Section 1
		{ 	  10.60*fgkmm,  35.50*fgkmm,  35.90*fgkmm },	// Section 2
		{ 	  13.20*fgkmm,  35.50*fgkmm,  35.90*fgkmm }}	// Section 3
		,	// Layer 2
		{{ 	  10.20*fgkmm,  42.50*fgkmm,  43.90*fgkmm },	// Section 0
		{ 	  10.60*fgkmm,  42.50*fgkmm,  43.90*fgkmm },	// Section 1
		{ 	  10.60*fgkmm,  43.50*fgkmm,  43.90*fgkmm },	// Section 2
		{ 	  13.20*fgkmm,  43.50*fgkmm,  43.90*fgkmm }}	// Section 3
	};
	
		
	static const Int_t kNumOfConePPanelSection = 2;
	static const Double_t kConePPanelArray[(int)fgkNumberOfInnerLayers][kNumOfConePPanelSection][3] =
	{		// Layer 0
		{{ 	 170.00*fgkmm,  29.50*fgkmm,  56.50*fgkmm },	// Section 0
		{ 	 182.00*fgkmm,  29.50*fgkmm,  56.50*fgkmm }}	// Section 1
		,	// Layer 1
		{{ 	 164.00*fgkmm,  60.00*fgkmm,  89.50*fgkmm },	// Section 0
		{ 	 176.00*fgkmm,  60.00*fgkmm,  89.50*fgkmm }}	// Section 1
		,	// Layer 2
		{{ 	 158.00*fgkmm,  90.00*fgkmm, 119.50*fgkmm },	// Section 0
		{ 	 170.00*fgkmm,  90.00*fgkmm, 119.50*fgkmm }}	// Section 1
	};


	TGeoVolume* vol;
	TGeoShape* shp;
	
			
	shp = new TGeoPcon(TString("WheelConeA")+=fLayerNumber, 0, 180, kNumOfConeSection);
	for (int i=0; i<kNumOfConeSection; i++)
		((TGeoPcon*)shp)->DefineSection(i, kConeArray[(int)fLayerNumber][i][0], kConeArray[(int)fLayerNumber][i][1], kConeArray[(int)fLayerNumber][i][2]);	
	vol = new TGeoVolume(TString("WheelConeA")+=fLayerNumber, shp, med);
		vol->SetFillColor(color);
		vol->SetLineColor(color);
	ewhVol->AddNode(vol, 1, new TGeoCombiTrans(0, 0, 0, new TGeoRotation("", 0, 0, 0)));

	shp = new TGeoPcon(TString("WheelConeTubeA")+=fLayerNumber, 0, 180, kNumOfConeTubeSection);
	for (int i=0; i<kNumOfConeTubeSection; i++)
		((TGeoPcon*)shp)->DefineSection(i, kConeTubeArray[(int)fLayerNumber][i][0], kConeTubeArray[(int)fLayerNumber][i][1], kConeTubeArray[(int)fLayerNumber][i][2]);	
	vol = new TGeoVolume(TString("WheelConeTubeA")+=fLayerNumber, shp, med);
		vol->SetFillColor(color);
		vol->SetLineColor(color);
	ewhVol->AddNode(vol, 1, new TGeoCombiTrans(0, 0, 0, new TGeoRotation("", 0, 0, 0)));

	shp = new TGeoPcon(TString("WheelConeEdgeA")+=fLayerNumber, 0, 180, kNumOfConeEdgeSection);
	for (int i=0; i<kNumOfConeEdgeSection; i++)
		((TGeoPcon*)shp)->DefineSection(i, kConeEdgeArray[(int)fLayerNumber][i][0], kConeEdgeArray[(int)fLayerNumber][i][1], kConeEdgeArray[(int)fLayerNumber][i][2]);	
	vol = new TGeoVolume(TString("WheelConeEdgeA")+=fLayerNumber, shp, med);
		vol->SetFillColor(color);
		vol->SetLineColor(color);
	ewhVol->AddNode(vol, 1, new TGeoCombiTrans(0, 0, 0, new TGeoRotation("", 0, 0, 0)));
	
	// -- OVERLAPS WITH STAVE_STRUCT --
	shp = new TGeoPcon(TString("WheelConeLA")+=fLayerNumber, 0, 180, kNumOfConeLSection);
	for (int i=0; i<kNumOfConeLSection; i++)
		((TGeoPcon*)shp)->DefineSection(i, kConeLArray[(int)fLayerNumber][i][0], kConeLArray[(int)fLayerNumber][i][1], kConeLArray[(int)fLayerNumber][i][2]);	
	vol = new TGeoVolume(TString("WheelConeLA")+=fLayerNumber, shp, med);
	vol->SetFillColor(color);	vol->SetLineColor(color);
	ewhVol->AddNode(vol, 1, new TGeoCombiTrans(0, 0, 0, new TGeoRotation("", 0, 0, 0)));

	// Patch Panel
	//shp = new TGeoPcon(TString("WheelConePPanelA")+=fLayerNumber, 0, 180, kNumOfConePPanelSection);
	//for (int i=0; i<kNumOfConePPanelSection; i++)
		//((TGeoPcon*)shp)->DefineSection(i, kConePPanelArray[(int)fLayerNumber][i][0], kConePPanelArray[(int)fLayerNumber][i][1], kConePPanelArray[(int)fLayerNumber][i][2]);	
	//vol = new TGeoVolume(TString("WheelConePPanelA")+=fLayerNumber, shp, med);
		//vol->SetFillColor(color);
		//vol->SetLineColor(color);
	//ewhVol->AddNode(vol, 1, new TGeoCombiTrans(0, 0, 0, new TGeoRotation("", 0, 0, 0)));
	
	//if (fLayerNumber==0)
	//{
		//shp = new TGeoPcon(TString("WheelConeTubeA")+=fLayerNumber, 0, 180, 2);
		//for (int i=0; i<2; i++)
			//((TGeoPcon*)shp)->DefineSection(i, kConeTubeArray[(int)fLayerNumber][i][0], kConeArray[(int)fLayerNumber][i][1], kConeArray[(int)fLayerNumber][i][2]);	
		//vol = new TGeoVolume(TString("WheelConeTubeA")+=fLayerNumber, shp, med);
		//vol->SetFillColor(color);
		//vol->SetLineColor(color);
		//ewhVol->AddNode(vol, 1, new TGeoCombiTrans(0, 0, +kOffsetZ, new TGeoRotation("", 0, 0, 0)));
	//}

	return ewhVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::Create_ALC_0337_xxC_ous(const Int_t color, const TGeoMedium* med){

	TGeoVolume* ewhVol = new TGeoVolumeAssembly(TString("EndWheelC") += fLayerNumber);

	//const Double_t kOffsetZ = -18.0*fgkmm;
	//const Double_t fkPlateThick = 0.4*fgkmm;
	//const Double_t fkHndLength = 13.0*fgkmm; // also use in CreateInnerBHandleA.
		
	//const Double_t kHndPlateRadius[7]={27.9*fgkmm, 35.9*fgkmm, 43.4*fgkmm,0,0,0,0};
	//const Double_t kContactRadius[7]={25.9*fgkmm, 33.8*fgkmm, 41.1*fgkmm,0,0,0,0}; // approx : -2mm from Radius
		
	//const Double_t fkHndWidthAngle = 360/fNStaves;
	static const Int_t kNumOfConeTubeSection = 4;
	static const Double_t kConeTubeArray[(int)fgkNumberOfInnerLayers][kNumOfConeTubeSection][3] =
	{		// Layer 0
		{{ 	 -10.00*fgkmm,  27.90*fgkmm,  28.50*fgkmm },	// Section 0
		{ 	 -20.20*fgkmm,  27.90*fgkmm,  28.50*fgkmm },	// Section 1
		{ 	 -20.20*fgkmm,  22.25*fgkmm,  28.50*fgkmm },	// Section 2
		{ 	 -21.00*fgkmm,  22.25*fgkmm,  28.50*fgkmm }}	// Section 3
		,	// Layer 1
		{{ 	 -10.00*fgkmm,  35.90*fgkmm,  36.50*fgkmm },	// Section 0
		{ 	 -17.00*fgkmm,  35.90*fgkmm,  36.50*fgkmm },	// Section 1
		{ 	 -17.00*fgkmm,  29.00*fgkmm,  36.50*fgkmm },	// Section 2
		{ 	 -18.50*fgkmm,  29.00*fgkmm,  36.50*fgkmm }}	// Section 3
		,	// Layer 2
		{{ 	 -10.00*fgkmm,  43.90*fgkmm,  44.50*fgkmm },	// Section 0
		{ 	 -15.20*fgkmm,  43.90*fgkmm,  44.50*fgkmm },	// Section 1
		{ 	 -15.20*fgkmm,  37.00*fgkmm,  44.50*fgkmm },	// Section 2
		{ 	 -16.00*fgkmm,  37.00*fgkmm,  44.50*fgkmm }}	// Section 3
	};
	
	static const Int_t kNumOfConeLSection = 4;
	static const Double_t kConeLArray[(int)fgkNumberOfInnerLayers][kNumOfConeLSection][3] =
	{		// Layer 0
		{{ 	 -10.20*fgkmm,  22.00*fgkmm,  27.90*fgkmm },	// Section 0
		{ 	 -10.60*fgkmm,  22.00*fgkmm,  27.90*fgkmm },	// Section 1
		{ 	 -10.60*fgkmm,  27.50*fgkmm,  27.90*fgkmm },	// Section 2
		{ 	 -14.20*fgkmm,  27.50*fgkmm,  27.90*fgkmm }}	// Section 3
		,	// Layer 1
		{{ 	 -10.20*fgkmm,  29.00*fgkmm,  35.90*fgkmm },	// Section 0
		{ 	 -10.60*fgkmm,  29.00*fgkmm,  35.90*fgkmm },	// Section 1
		{ 	 -10.60*fgkmm,  35.50*fgkmm,  35.90*fgkmm },	// Section 2
		{ 	 -14.20*fgkmm,  35.50*fgkmm,  35.90*fgkmm }}	// Section 3
		,	// Layer 2
		{{ 	 -10.20*fgkmm,  37.00*fgkmm,  43.90*fgkmm },	// Section 0
		{ 	 -10.60*fgkmm,  37.00*fgkmm,  43.90*fgkmm },	// Section 1
		{ 	 -10.60*fgkmm,  43.50*fgkmm,  43.90*fgkmm },	// Section 2
		{ 	 -14.20*fgkmm,  43.50*fgkmm,  43.90*fgkmm }}	// Section 3
	};
	

	TGeoVolume* vol;
	TGeoShape* shp;
	
			

	shp = new TGeoPcon(TString("WheelConeTubeC")+=fLayerNumber, 0, 180, kNumOfConeTubeSection);
	for (int i=0; i<kNumOfConeTubeSection; i++)
		((TGeoPcon*)shp)->DefineSection(i, kConeTubeArray[(int)fLayerNumber][i][0], kConeTubeArray[(int)fLayerNumber][i][1], kConeTubeArray[(int)fLayerNumber][i][2]);	
	vol = new TGeoVolume(TString("WheelConeTubeC")+=fLayerNumber, shp, med);
		vol->SetFillColor(color);
		vol->SetLineColor(color);
	ewhVol->AddNode(vol, 1, new TGeoCombiTrans(0, 0, 0, new TGeoRotation("", 0, 0, 0)));

	
	shp = new TGeoPcon(TString("WheelConeLC")+=fLayerNumber, 0, 180, kNumOfConeLSection);
	for (int i=0; i<kNumOfConeLSection; i++)
		((TGeoPcon*)shp)->DefineSection(i, kConeLArray[(int)fLayerNumber][i][0], kConeLArray[(int)fLayerNumber][i][1], kConeLArray[(int)fLayerNumber][i][2]);	
	vol = new TGeoVolume(TString("WheelConeLC")+=fLayerNumber, shp, med);
		vol->SetFillColor(color);
		vol->SetLineColor(color);
	ewhVol->AddNode(vol, 1, new TGeoCombiTrans(0, 0, 0, new TGeoRotation("", 0, 0, 0)));


	return ewhVol;
}



//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::Create_ALC_0337_xxC_ins(const TGeoMedium* med){

	TString xStr, xTmp1, xTmp2;
	TGeoShape* shp;
	
	// ----------------- Create the plate -----------------
	xTmp1 = (TString("shpHndPlateC")+=fLayerNumber) + (TString("_")+=(int)fNHandleCreated);
	shp = new TGeoPcon(xTmp1, 0, 180., 2);
		((TGeoPcon*)shp)->DefineSection(0, 0.0						, fkALC_0334_PlateRadius[fLayerNumber]	, fkALC_0334_PlateRadius[fLayerNumber] + fkALC_0334_PlateThick);
		((TGeoPcon*)shp)->DefineSection(1, -fkALC_0334_HndLength	, fkALC_0334_PlateRadius[fLayerNumber]	, fkALC_0334_PlateRadius[fLayerNumber] + fkALC_0334_PlateThick);
	//xTmp2 = (TString("xtrHndPlateC")+=fLayerNumber)  + (TString("_")+=(int)fNHandleCreated);
		//(new TGeoCombiTrans(xTmp2, 0, 0, fkALC_0334_HoleFromEdge, new TGeoRotation("", 0, 0, 0)))->RegisterYourself();
	//xStr = xTmp1 + ":" + xTmp2;
	
		
	//return new TGeoVolume((TString("ALC_0337_0") += (int)fLayerNumber) + TString("C"), new TGeoCompositeShape("", xStr), med); 			
	return new TGeoVolume((TString("ALC_0337_0") += (int)fLayerNumber) + TString("C"), shp, med);
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::Create_ALC_0334_xxC(const TGeoMedium* med){
	//
	//
	//
	const Double_t fkHoleRadius = fkALC_0334_HoleDia/2;

	const Double_t fkHoleFromEdge1 = fkALC_0334_HoleFromEdge - fkHoleRadius;
	const Double_t fkHoleFromEdge2 = fkALC_0334_HoleFromEdge + fkHoleRadius;		

	TString xStr, xTmp1, xTmp2;
	TGeoShape* shp;
	
	// ----------------- Create the plate -----------------
	//xTmp1 = TString("shpHndPlate")+=fLayerNumber;
	//shp = new TGeoPcon(xTmp1, 0, 360., 2);
		//((TGeoPcon*)shp)->DefineSection(0, 0.0						, fkALC_0334_PlateRadius[fLayerNumber]	, fkALC_0334_PlateRadius[fLayerNumber] + fkALC_0334_PlateThick);
		//((TGeoPcon*)shp)->DefineSection(1, -fkALC_0334_HndLength	, fkALC_0334_PlateRadius[fLayerNumber]	, fkALC_0334_PlateRadius[fLayerNumber] + fkALC_0334_PlateThick);
	//xTmp2 = TString("xtrHndPlate")+=fLayerNumber;
		//(new TGeoCombiTrans(xTmp2, 0, 0, fkALC_0334_HoleFromEdge, new TGeoRotation("", 0, 0, 0)))->RegisterYourself();
	//xStr = xTmp1 + ":" + xTmp2;

	// ----------------- Create the base -----------------
	xTmp1 = (TString("shpHndBaseC")+=fLayerNumber) + (TString("_")+=(int)fNHandleCreated);
		shp = new TGeoBBox(xTmp1, fkALC_0334_HndWidth/2, fkALC_0334_HndThick/2, fkALC_0334_HndLength/2);
	xTmp2 = (TString("xtrHndBaseC")+=fLayerNumber) + (TString("_")+=(int)fNHandleCreated);
		(new TGeoCombiTrans(xTmp2, -fkALC_0334_HndWidth/2+fkALC_0334_HoleFromEdge+fkALC_0334_ContactDX[fLayerNumber], +fkALC_0334_HndThick/2+fkALC_0334_ContactDY[fLayerNumber], -fkALC_0334_HndLength/2+fkALC_0334_HoleFromEdge, new TGeoRotation("", 0, 0, 0)))->RegisterYourself();

	#ifdef _ITSSB_SHOW_COMBISHAPE
		xStr = xTmp1 + ":" + xTmp2;							// without plate
		//xStr += TString(" + ") + xTmp1 + ":" + xTmp2;		// with plate
	#else
	{
		TString xStrAA, xTmp1AA, xTmp2AA;
		TGeoShape* shpAA;
		
		xTmp1AA = (TString("shpHndCutPlateC")+=fLayerNumber) + (TString("_")+=(int)fNHandleCreated);
		shpAA = new TGeoPcon(xTmp1AA, 0, 360., 2);
			((TGeoPcon*)shpAA)->DefineSection(0, 0.0							, fkALC_0334_PlateRadius[fLayerNumber] 	, fkALC_0334_PlateRadius[fLayerNumber] + fkALC_0334_HndThick);
			((TGeoPcon*)shpAA)->DefineSection(1, -fkALC_0334_HndLength-2*fgkmm	, fkALC_0334_PlateRadius[fLayerNumber] 	, fkALC_0334_PlateRadius[fLayerNumber] + fkALC_0334_HndThick);
		
		xTmp2AA = (TString("xtrHndCutPlateC") += fLayerNumber) + (TString("_")+=(int)fNHandleCreated);
			(new TGeoCombiTrans(xTmp2AA, 0, 0, fkALC_0334_HoleFromEdge+1*fgkmm, new TGeoRotation("", 0, 0, 0)))->RegisterYourself();		
		xStrAA = TString("(") + xTmp1 + ":" + xTmp2 + " - " + xTmp1AA + ":" + xTmp2AA + ")";
		xStr = xStrAA;										// without plate
		//xStr += TString(" + ") + xStrAA;					// with plate

	}
	#endif


	// ----------------- Make a hole -----------------
	xTmp1 = (TString("shpHndHoleC")+=fLayerNumber) + (TString("_")+=(int)fNHandleCreated);
		shp = new TGeoTube(xTmp1, 0, fkHoleRadius, fkALC_0334_HoleLength);
	xTmp2 = (TString("xtrHndHoleC")+=fLayerNumber) + (TString("_")+=(int)fNHandleCreated);
		(new TGeoCombiTrans(xTmp2, +fkALC_0334_ContactDX[fLayerNumber], +fkALC_0334_ContactDY[fLayerNumber], 0, new TGeoRotation("", 0, 90, 0)))->RegisterYourself();
	#ifdef _ITSSB_SHOW_COMBISHAPE
		xStr += TString(" + ") + xTmp1 + ":" + xTmp2;
	#else
		xStr = TString("(")+ xStr + TString(")") + TString(" - ") + xTmp1 + ":" + xTmp2;
	#endif
	
	return new TGeoVolume((TString("ALC_0334_0") += (int)fLayerNumber) + (TString("C_") += (int)fNHandleCreated++), new TGeoCompositeShape("", xStr), med); 			
}



//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateInnerBEWheelA0(const TGeoMedium* med){
	//
	//
	//


	TGeoVolume* ewhVol = new TGeoVolumeAssembly(TString("EndWheelA") += fLayerNumber);


	const Double_t ITSSBIB01SWSkinThick = 0.01 *fgkcm;
	const Double_t ITSSBIB01SWCoreThick = 0.18 *fgkcm;
	const Double_t ITSSBIB01SWThick = ITSSBIB01SWCoreThick + 2*ITSSBIB01SWSkinThick;

	const Double_t ITSSBIB01CutOffset = 0.20;
	const Double_t ITSSBIB01Length = 40.00;
	const Double_t ITSSBIB01Radius = 4.78; 
	const Double_t ITSSBIB01MinRadius = 2.20;
	const Double_t ITSSBIB01MaxRadius = ITSSBIB01Radius + ITSSBIB01SWThick;

	// box to div
	TGeoBBox *xbboxdiv = new TGeoBBox("xbboxdiv", ITSSBIB01MaxRadius+1.0, ITSSBIB01MaxRadius+1.0, (ITSSBIB01Length/2)+1.0); // extend edge about 1.0
	TGeoTranslation *xtrandiv = new TGeoTranslation("xtrandiv", 0, 0-xbboxdiv->GetDY()-ITSSBIB01CutOffset, 0);
	xtrandiv->RegisterYourself();			


	if (fLayerNumber==0)
	{ // 0121x : Layer 0
		// 01211 = the object
		// 01212 = -------------
		// 01213 = the base
		// 01214 = the tube for holes
		// 01215 = the tube for cut
		
		Double_t ITSSBIB01211Thick = 0.04;
		Double_t ITSSBIB01213Thick = 0.45; // over approximate, vary this to offset for contact with stave

		Double_t ITSSBIB01211Length = 5.70;
		Double_t ITSSBIB01212Length = 1.90;
		Double_t ITSSBIB01213Length = 1.30;
		Double_t ITSSBIB01214Length = ITSSBIB01213Thick*2.5; // 2.5 time is easy to make a hole without offset this tube

		Double_t ITSSBIB0121Length = ITSSBIB01211Length;
		Double_t ITSSBIB01215Length = ITSSBIB0121Length*1.5; // 1.5 time is easy to make tube for div without offset its

		Double_t ITSSBIB01213Width = 1.10; // over approximate, vary this to offset for contact with stave

		Double_t ITSSBIB01211Radius = 2.85-ITSSBIB01211Thick; // = ITSSBIB01103Radius-ITSSBIB01211Thick; 
		Double_t ITSSBIB01212Radius = ITSSBIB01211Radius + ITSSBIB01211Length - ITSSBIB01211Thick;
		Double_t ITSSBIB01213Radius = 2.576; // 2.57673598, L1 radius for circle, offset at top surface of each object
		Double_t ITSSBIB01214Radius = 0.25; // ***IMPORTANT*** this is radius each hole, for circle we must use ITSSBIB01213Radius ***IMPORTANT***
		Double_t ITSSBIB01215Radius = ITSSBIB01211Radius + 0.01; // TO MAKE SURE WE CAN WELDED TO 1111&1112 (A bit should < ITSSBIB01211Radius+ITSSBIB01211Thick

		Double_t ITSSBIB0121Radius = ITSSBIB01211Radius;
		
		TGeoRotation *xrotITSSBIB0121Circle = new TGeoRotation("xrotITSSBIB0121Circle", 0.0, 0.0,  - 90.0 + 16.37); // this rotate BEFORE div
		xrotITSSBIB0121Circle->RegisterYourself();
		
		Double_t ITSSBIB01213EachRotate = 15.0; // degree of rotate each object of 01213
		Double_t ITSSBIB01214EachRotate = 15.0;
		
		Double_t ITSSBIB0121MaxRadius = ITSSBIB01211Radius + ITSSBIB01211Length;
		
		// offset is compared with volume's shape of section1
		Double_t ITSSBIB0121OffsetX = 0.0;
		Double_t ITSSBIB0121OffsetY = 0.0;
		Double_t ITSSBIB0121OffsetZ = 0.0 - 11.5 - (ITSSBIB0121Length/2); // +OffsetAtTheHead
		TGeoRotation *ITSSBIB0121Rotate = new TGeoRotation("", 0.0, 0.0, 0.0); // this rotate AFTER div
		TGeoCombiTrans *ITSSBIB0121CombT = new TGeoCombiTrans(ITSSBIB0121OffsetX, ITSSBIB0121OffsetY, ITSSBIB0121OffsetZ, ITSSBIB0121Rotate);

		TGeoTube *xtubeITSSBIB01211 = new TGeoTube("xtubeITSSBIB01211"
			, ITSSBIB01211Radius, ITSSBIB01211Radius + ITSSBIB01211Thick, ITSSBIB01211Length/2);	
		TGeoTranslation *xtranITSSBIB01211 = new TGeoTranslation("xtranITSSBIB01211", 0.0, 0.0, 0.0);
		xtranITSSBIB01211->RegisterYourself();			


		Int_t NDegPerStave = 360/fNStaves;
		Double_t OffNegSideX01213 = ITSSBIB01213Width/2 - 0.5;
		Double_t Deg2RadFactor = TMath::Pi()/180;
		Double_t CosEachRotate = TMath::Cos(ITSSBIB01213EachRotate*Deg2RadFactor);
		Double_t SinEachRotate = TMath::Sin(ITSSBIB01213EachRotate*Deg2RadFactor);
		Double_t FaceCirRadi01213 = ITSSBIB01213Thick/2 + ITSSBIB01213Radius;
		
		for (Int_t i=0; i<fNStaves; i++)
		{
			Double_t CosDegIter = TMath::Cos(i*NDegPerStave*Deg2RadFactor);
			Double_t SinDegIter = TMath::Sin(i*NDegPerStave*Deg2RadFactor);					

			{ // 01213 : base
				TString str("xbboxITSSBIB01213");
				str += i;
				TGeoBBox *xbboxITSSBIB001213x = new TGeoBBox(str, ITSSBIB01213Width/2, ITSSBIB01213Thick/2, ITSSBIB01213Length/2);
				
				TString str2("xcombtITSSBIB01213");
				str2 += i;
				TGeoCombiTrans *xcombtITSSBIB01213x = new TGeoCombiTrans(str2
					, 0 + (FaceCirRadi01213 * SinDegIter)	// circle with top face offset
						+ (OffNegSideX01213 * CosEachRotate) * CosDegIter // x-offset
						- (OffNegSideX01213 * SinEachRotate) * SinDegIter // y-offset
					, 0 - (FaceCirRadi01213 * CosDegIter)	// circle	with top face offset
						+ (OffNegSideX01213 * CosEachRotate) * SinDegIter // x-offset
						+ (OffNegSideX01213 * SinEachRotate) * CosDegIter // y-offset  (*** ??? why positive, not negative ??? ***)
					, 0 - (ITSSBIB01213Length/2) + (ITSSBIB0121Length/2) 	// set z to the head of parent object
					
					, new TGeoRotation("", 0.0, 0.0, i*NDegPerStave + ITSSBIB01213EachRotate));
				
				xcombtITSSBIB01213x->RegisterYourself();
			}
			
			{ // 01214 : tube for holes
				TString str("xtubeITSSBIB01214");
				str += i;
				TGeoTube *xtubeITSSBIB001214x = new TGeoTube(str, 0.0, ITSSBIB01214Radius, ITSSBIB01214Length/2);
				
				TString str2("xcombtITSSBIB01214");
				str2 += i;
				TGeoCombiTrans *xcombtITSSBIB01214x = new TGeoCombiTrans(str2
					, 0 + (ITSSBIB01213Radius * SinDegIter)	// circle *** USE ITSSBIB01213Radius ***
					, 0 - (ITSSBIB01213Radius * CosDegIter)	// circle	*** USE ITSSBIB01213Radius ***
					, 0 - .40 + (ITSSBIB0121Length/2) 		// set z to the head of parent object
					// this have been set the position to circle around tube, also offset to the top face(Y-POS) *** BUT NOT X-POS ***
					, new TGeoRotation("", i*NDegPerStave + ITSSBIB01214EachRotate, 90.0, 0.0)); // perpendicular to 01213
				
				xcombtITSSBIB01214x->RegisterYourself();
			}					
			
			#ifdef _ITSSB_SHOW_INDICATOR
			{ // #indicator
				TString str("xtubeITSSBIB01211xind");
				str += i;
				TGeoTube *xtubeITSSBIB001211xind = new TGeoTube(str, 0.0, 0.01, ITSSBIB01213Length);
				
				TString str2("xcombtITSSBIB01211xind");
				str2 += i;
				TGeoCombiTrans *xcombtITSSBIB01211xind = new TGeoCombiTrans(str2
					, 0 + (ITSSBIB01213Radius * SinDegIter)	// circle *** USE ITSSBIB01213Radius ***
					, 0 - (ITSSBIB01213Radius * CosDegIter)	// circle	*** USE ITSSBIB01213Radius ***
					, 0 + (ITSSBIB0121Length/2)		// set z to the head of parent object
					// this have been set the position to circle around tube, also offset to the top face(Y-POS) *** BUT NOT X-POS ***
					, new TGeoRotation("", 0.0, 0.0, 0.0)); 
				
				xcombtITSSBIB01211xind->RegisterYourself();
			}
			#endif					
			
		}

		TGeoTube *xtubeITSSBIB01215 = new TGeoTube("xtubeITSSBIB01215", 0.0, ITSSBIB01215Radius, ITSSBIB01215Length/2);
		TGeoTranslation *xtranITSSBIB01215 = new TGeoTranslation("xtranITSSBIB01215", 0.0, 0.0, 0.0);
		xtranITSSBIB01215->RegisterYourself();			


		// combine all 0121x objects
		TString xstrITSSBIB0121cs("(");				
		TString xstrITSSBIB01214trunks("("); // collect 14, all holes to div
		for (Int_t i=0; i<fNStaves; i++) // ITSSBIB01213
		{ // *** IF PLACE CODE 114 BELOW 113, ITS POSITION MIGHT BE UNDETERMINED ***
			
			TString xstrITSSBIB01214trunk("xtubeITSSBIB01214");
			xstrITSSBIB01214trunk += i;
			xstrITSSBIB01214trunk += ":xcombtITSSBIB01214";
			xstrITSSBIB01214trunk += i;
			
			if (i>0) xstrITSSBIB01214trunks += " + ";
			xstrITSSBIB01214trunks += xstrITSSBIB01214trunk;

			if (i>0) xstrITSSBIB0121cs += " + ";
			xstrITSSBIB0121cs += "(xbboxITSSBIB01213";
			xstrITSSBIB0121cs += i;
			xstrITSSBIB0121cs += ":xcombtITSSBIB01213";
			xstrITSSBIB0121cs += i; 
			
			#ifdef _ITSSB_SHOW_COMBISHAPE
			xstrITSSBIB0121cs += "*xtubeITSSBIB01215:xtranITSSBIB01215 - ";	// cut & make a hole
			//xstrITSSBIB0121cs += " - ";	// cut & make a hole
			#else
			xstrITSSBIB0121cs += " + ";		// just show
			#endif
			xstrITSSBIB0121cs += xstrITSSBIB01214trunk;												
			
			xstrITSSBIB0121cs += ")"; 
			
		}
		xstrITSSBIB01214trunks += ")";
		//xstrITSSBIB0121cs += "+ xtubeITSSBIB01215:xtranITSSBIB01215";	// cut & make a hole
		
		// 01211
		xstrITSSBIB0121cs += " + (xtubeITSSBIB01211:xtranITSSBIB01211";
		#ifdef _ITSSB_SHOW_COMBISHAPE
		xstrITSSBIB0121cs += " - ";
		xstrITSSBIB0121cs += xstrITSSBIB01214trunks;	// make holes
		#endif
		xstrITSSBIB0121cs += ")";
		
														
		xstrITSSBIB0121cs += "):xrotITSSBIB0121Circle";
		#ifdef _ITSSB_SHOW_WITHBOXDIV
		xstrITSSBIB0121cs += " * (xbboxdiv:xtrandiv)";
		#endif 
		
		TGeoCompositeShape *xcompITSSBIB0121 = new TGeoCompositeShape("xcompITSSBIB0121cs", xstrITSSBIB0121cs);
		TGeoVolume *xvolITSSBIB0121 = new TGeoVolume("ITSSBIB0121x", xcompITSSBIB0121, med); 			
		xvolITSSBIB0121->SetFillColor(7);
		xvolITSSBIB0121->SetLineColor(7);
		ewhVol->AddNode(xvolITSSBIB0121, 1, ITSSBIB0121CombT);
		
		#ifdef _ITSSB_SHOW_INDICATOR
		{ // #indicator
			TString xstrITSSBIB01211xind("(");			
			for (Int_t i=0; i<fNStaves; i++) // ITSSBIB01211xind
			{
				if (i>0) xstrITSSBIB01211xind += " + ";
				xstrITSSBIB01211xind += "xtubeITSSBIB01211xind";
				xstrITSSBIB01211xind += i;
				xstrITSSBIB01211xind += ":xcombtITSSBIB01211xind";
				xstrITSSBIB01211xind += i;
			}	
			xstrITSSBIB01211xind += ")";

			TGeoCompositeShape *xcompITSSBIB01211xind = new TGeoCompositeShape("xcompITSSBIB01211xind", xstrITSSBIB01211xind);
			TGeoVolume *xvolITSSBIB01211xind = new TGeoVolume("ITSSBIB01211xind", xcompITSSBIB01211xind, med); 			
			xvolITSSBIB01211xind->SetFillColor(2);
			xvolITSSBIB01211xind->SetLineColor(2);
			ewhVol->AddNode(xvolITSSBIB01211xind, 1, ITSSBIB0121CombT);
		}
		#endif
	}

	if (fLayerNumber==1)
	{ // 0122x : Layer 1
		// 01221 = the object
		// 01222 = -------------
		// 01223 = the base
		// 01224 = the tube for holes
		// 01224 = the tube for cut

		Double_t ITSSBIB01221Thick = 0.04;
		Double_t ITSSBIB01223Thick = 0.45; // over approximate, vary this to offset for contact with stave

		Double_t ITSSBIB01221Radius1 = 2.85; // = ITSSBIB01103Radius //
		Double_t ITSSBIB01221Radius2 = 3.65 - ITSSBIB01221Thick; // = ITSSBIB01105Radius-ITSSBIB01221Thick; //

		Double_t ITSSBIB01221Length = 5.20;
		Double_t ITSSBIB01222Length = 1.90;
		Double_t ITSSBIB01223Length = 1.30;
		Double_t ITSSBIB01224Length = ITSSBIB01223Thick*2.5; // 2.5 time is easy to make a hole without offset this tube

		Double_t ITSSBIB0122Length = ITSSBIB01221Length;
		Double_t ITSSBIB01225Length = ITSSBIB0122Length*1.5; // 1.5 time is easy to make tube for div without offset its

		Double_t ITSSBIB01223Width = 1.10; // over approximate, vary this to offset for contact with stave
		Double_t ITSSBIB01221CenterWidth = (ITSSBIB01221Radius2+(ITSSBIB01221Thick/2))-(ITSSBIB01221Radius1-(ITSSBIB01221Thick/2));

		Double_t ITSSBIB01221Radius = 2.85-ITSSBIB01221Thick; // = ITSSBIB01103Radius-ITSSBIB01221Thick; 
		Double_t ITSSBIB01222Radius = ITSSBIB01221Radius + ITSSBIB01221Length - ITSSBIB01221Thick;
		Double_t ITSSBIB01223Radius = 3.38; // 3.38059, radius for circle, offset at top surface of each object
		Double_t ITSSBIB01224Radius = 0.25; // ***IMPORTANT*** this is radius each hole, for circle we must use ITSSBIB01223Radius ***IMPORTANT***
		Double_t ITSSBIB01225Radius = ITSSBIB01221Radius2 + 0.01; // TO MAKE SURE WE CAN WELDED TO 1111&1112 (A bit should < ITSSBIB01221Radius+ITSSBIB01221Thick

		Double_t ITSSBIB0122Radius = ITSSBIB01221Radius1;

		// calculate at 45 degree
		Double_t ITSSBIB01221Length1 = 0.34; // fixed parameter
		Double_t ITSSBIB01221Length2 = ITSSBIB01221CenterWidth-(ITSSBIB01221Thick/2)-(ITSSBIB01221Thick/2); // followed parameter
		Double_t ITSSBIB01221Length3 = ITSSBIB01221Length-ITSSBIB01221Length1-ITSSBIB01221Length2; // followed parameter

		Double_t ITSSBIB01221Z0 = 0.00; // The center of object is not the same as the center of space.
		Double_t ITSSBIB01221Z1 = -(ITSSBIB01221Length3 - (ITSSBIB01221Thick/2));
		Double_t ITSSBIB01221Z2 = -(ITSSBIB01221Length3 + (ITSSBIB01221Thick/2));
		Double_t ITSSBIB01221Z3 = -(ITSSBIB01221Length3 + ITSSBIB01221Length2 - (ITSSBIB01221Thick/2));
		Double_t ITSSBIB01221Z4 = -(ITSSBIB01221Length3 + ITSSBIB01221Length2 + (ITSSBIB01221Thick/2));
		Double_t ITSSBIB01221Z5 = -(ITSSBIB01221Length);

		Double_t ITSSBIB01221RMin0 = ITSSBIB01221Radius2, ITSSBIB01221RMax0 = ITSSBIB01221Radius2+ITSSBIB01221Thick;
		Double_t ITSSBIB01221RMin1 = ITSSBIB01221Radius2, ITSSBIB01221RMax1 = ITSSBIB01221Radius2+ITSSBIB01221Thick;
		Double_t ITSSBIB01221RMin2 = ITSSBIB01221Radius2 - ITSSBIB01221Thick, ITSSBIB01221RMax2 = ITSSBIB01221Radius2+ITSSBIB01221Thick;
		Double_t ITSSBIB01221RMin3 = ITSSBIB01221Radius1, ITSSBIB01221RMax3 = ITSSBIB01221Radius1+(2*ITSSBIB01221Thick);
		Double_t ITSSBIB01221RMin4 = ITSSBIB01221Radius1, ITSSBIB01221RMax4 = ITSSBIB01221Radius1+ITSSBIB01221Thick;
		Double_t ITSSBIB01221RMin5 = ITSSBIB01221Radius1, ITSSBIB01221RMax5 = ITSSBIB01221Radius1+ITSSBIB01221Thick;

		
		TGeoRotation *xrotITSSBIB0122Circle = new TGeoRotation("xrotITSSBIB0122Circle", 0.0, 0.0,  - 90.0 + 12.03); // this rotate BEFORE div
		xrotITSSBIB0122Circle->RegisterYourself();
		
		Double_t ITSSBIB01223EachRotate = 15.0; // degree of rotate each object of 01223
		Double_t ITSSBIB01224EachRotate = 15.0;
		
		Double_t ITSSBIB0122MaxRadius = ITSSBIB01221Radius + ITSSBIB01221Length;
		
		// offset is compared with volume's shape of section1
		Double_t ITSSBIB0122OffsetX = 0.0;
		Double_t ITSSBIB0122OffsetY = 0.0;
		Double_t ITSSBIB0122OffsetZ = 0.0 - 11.5; // original at z=0, not the center of the object
		TGeoRotation *ITSSBIB0122Rotate = new TGeoRotation("", 0.0, 0.0, 0.0); // this rotate AFTER div
		TGeoCombiTrans *ITSSBIB0122CombT = new TGeoCombiTrans(ITSSBIB0122OffsetX, ITSSBIB0122OffsetY, ITSSBIB0122OffsetZ, ITSSBIB0122Rotate);


		TGeoPcon *xpconeITSSBIB01221 = new TGeoPcon("xpconeITSSBIB01221", 0, 360, 6);
		xpconeITSSBIB01221->DefineSection(0, ITSSBIB01221Z0, ITSSBIB01221RMin0, ITSSBIB01221RMax0);
		xpconeITSSBIB01221->DefineSection(1, ITSSBIB01221Z1, ITSSBIB01221RMin1, ITSSBIB01221RMax1);
		xpconeITSSBIB01221->DefineSection(2, ITSSBIB01221Z2, ITSSBIB01221RMin2, ITSSBIB01221RMax2);
		xpconeITSSBIB01221->DefineSection(3, ITSSBIB01221Z3, ITSSBIB01221RMin3, ITSSBIB01221RMax3);
		xpconeITSSBIB01221->DefineSection(4, ITSSBIB01221Z4, ITSSBIB01221RMin4, ITSSBIB01221RMax4);
		xpconeITSSBIB01221->DefineSection(5, ITSSBIB01221Z5, ITSSBIB01221RMin5, ITSSBIB01221RMax5);			


		TGeoTranslation *xtranITSSBIB01221 = new TGeoTranslation("xtranITSSBIB01221", 0.0, 0.0, 0.0);
		xtranITSSBIB01221->RegisterYourself();			

		Int_t NDegPerStave = 360/fNStaves;
		Double_t OffNegSideX01223 = ITSSBIB01223Width/2 - 0.5;
		Double_t Deg2RadFactor = TMath::Pi()/180;
		Double_t CosEachRotate = TMath::Cos(ITSSBIB01223EachRotate*Deg2RadFactor);
		Double_t SinEachRotate = TMath::Sin(ITSSBIB01223EachRotate*Deg2RadFactor);
		Double_t FaceCirRadi01223 = ITSSBIB01223Thick/2 + ITSSBIB01223Radius;
		
		for (Int_t i=0; i<fNStaves; i++)
		{
			Double_t CosDegIter = TMath::Cos(i*NDegPerStave*Deg2RadFactor);
			Double_t SinDegIter = TMath::Sin(i*NDegPerStave*Deg2RadFactor);					

			{ // 012213
				TString str("xbboxITSSBIB01223");
				str += i;
				TGeoBBox *xbboxITSSBIB001223x = new TGeoBBox(str, ITSSBIB01223Width/2, ITSSBIB01223Thick/2, ITSSBIB01223Length/2);
				
				TString str2("xcombtITSSBIB01223");
				str2 += i;
				TGeoCombiTrans *xcombtITSSBIB01223x = new TGeoCombiTrans(str2
					, 0 + (FaceCirRadi01223 * SinDegIter)	// circle with top face offset
						+ (OffNegSideX01223 * CosEachRotate) * CosDegIter // x-offset
						- (OffNegSideX01223 * SinEachRotate) * SinDegIter // y-offset
					, 0 - (FaceCirRadi01223 * CosDegIter)	// circle	with top face offset
						+ (OffNegSideX01223 * CosEachRotate) * SinDegIter // x-offset
						+ (OffNegSideX01223 * SinEachRotate) * CosDegIter // y-offset  (*** ??? why positive, not negative ??? ***)
					, 0 - (ITSSBIB01223Length/2) 	// set z to the head of parent object
					
					, new TGeoRotation("", 0.0, 0.0, i*NDegPerStave + ITSSBIB01223EachRotate));
				
				xcombtITSSBIB01223x->RegisterYourself();
			}
			
			{ // 012214
				TString str("xtubeITSSBIB01224");
				str += i;
				TGeoTube *xtubeITSSBIB001224x = new TGeoTube(str, 0.0, ITSSBIB01224Radius, ITSSBIB01224Length/2);
				
				TString str2("xcombtITSSBIB01224");
				str2 += i;
				TGeoCombiTrans *xcombtITSSBIB01224x = new TGeoCombiTrans(str2
					, 0 + (ITSSBIB01223Radius * SinDegIter)	// circle *** USE ITSSBIB01223Radius ***
					, 0 - (ITSSBIB01223Radius * CosDegIter)	// circle	*** USE ITSSBIB01223Radius ***
					, 0 - .40 		// set z to the head of parent object
					// this have been set the position to circle around tube, also offset to the top face(Y-POS) *** BUT NOT X-POS ***
					, new TGeoRotation("", i*NDegPerStave + ITSSBIB01224EachRotate, 90.0, 0.0)); // perpendicular to 01223
				
				xcombtITSSBIB01224x->RegisterYourself();
			}					
			
			#ifdef _ITSSB_SHOW_INDICATOR
			{ // #indicator
				TString str("xtubeITSSBIB01221xind");
				str += i;
				TGeoTube *xtubeITSSBIB001221xind = new TGeoTube(str, 0.0, 0.01, ITSSBIB01223Length);
				
				TString str2("xcombtITSSBIB01221xind");
				str2 += i;
				TGeoCombiTrans *xcombtITSSBIB01221xind = new TGeoCombiTrans(str2
					, 0 + (ITSSBIB01223Radius * SinDegIter)	// circle *** USE ITSSBIB01223Radius ***
					, 0 - (ITSSBIB01223Radius * CosDegIter)	// circle	*** USE ITSSBIB01223Radius ***
					, 0		// set z to the head of parent object
					// this have been set the position to circle around tube, also offset to the top face(Y-POS) *** BUT NOT X-POS ***
					, new TGeoRotation("", 0.0, 0.0, 0.0)); 
				
				xcombtITSSBIB01221xind->RegisterYourself();
			}
			#endif					
			
		}

		TGeoTube *xtubeITSSBIB01225 = new TGeoTube("xtubeITSSBIB01225", 0.0, ITSSBIB01225Radius, ITSSBIB01225Length/2);
		TGeoTranslation *xtranITSSBIB01225 = new TGeoTranslation("xtranITSSBIB01225", 0.0, 0.0, 0.0);
		xtranITSSBIB01225->RegisterYourself();			


		// combine all 0122x objects
		TString xstrITSSBIB0122cs("(");				
		TString xstrITSSBIB01224trunks("("); // collect 24, all holes to div
		for (Int_t i=0; i<fNStaves; i++) // ITSSBIB01223
		{ // *** IF PLACE CODE 1224 BELOW 1223, ITS POSITION MIGHT BE UNDETERMINED ***
			
			TString xstrITSSBIB01224trunk("xtubeITSSBIB01224");
			xstrITSSBIB01224trunk += i;
			xstrITSSBIB01224trunk += ":xcombtITSSBIB01224";
			xstrITSSBIB01224trunk += i;
			
			if (i>0) xstrITSSBIB01224trunks += " + ";
			xstrITSSBIB01224trunks += xstrITSSBIB01224trunk;

			if (i>0) xstrITSSBIB0122cs += " + ";
			xstrITSSBIB0122cs += "(xbboxITSSBIB01223";
			xstrITSSBIB0122cs += i;
			xstrITSSBIB0122cs += ":xcombtITSSBIB01223";
			xstrITSSBIB0122cs += i; 
			
			#ifdef _ITSSB_SHOW_COMBISHAPE
			xstrITSSBIB0122cs += "*xtubeITSSBIB01225:xtranITSSBIB01225 - ";	// cut & make a hole
			//xstrITSSBIB0122cs += " - ";	// cut & make a hole
			#else
			xstrITSSBIB0122cs += " + ";		// just show
			#endif
			xstrITSSBIB0122cs += xstrITSSBIB01224trunk;												
			
			xstrITSSBIB0122cs += ")"; 
			
		}
		xstrITSSBIB01224trunks += ")";
		//xstrITSSBIB0122cs += " + xtubeITSSBIB01225:xtranITSSBIB01225";	// cut & make a hole

		// 01221
		xstrITSSBIB0122cs += " + (xpconeITSSBIB01221:xtranITSSBIB01221";
		#ifdef _ITSSB_SHOW_COMBISHAPE
		xstrITSSBIB0122cs += " - ";
		xstrITSSBIB0122cs += xstrITSSBIB01224trunks;	// make holes
		#endif
		xstrITSSBIB0122cs += ")";
		
														
		xstrITSSBIB0122cs += "):xrotITSSBIB0122Circle";
		#ifdef _ITSSB_SHOW_WITHBOXDIV
		xstrITSSBIB0122cs += " * (xbboxdiv:xtrandiv)";
		#endif 
		
		TGeoCompositeShape *xcompITSSBIB0122 = new TGeoCompositeShape("xcompITSSBIB0122cs", xstrITSSBIB0122cs);
		TGeoVolume *xvolITSSBIB0122 = new TGeoVolume("ITSSBIB0122x", xcompITSSBIB0122, med); 			
		xvolITSSBIB0122->SetFillColor(8);
		xvolITSSBIB0122->SetLineColor(8);
		ewhVol->AddNode(xvolITSSBIB0122, 1, ITSSBIB0122CombT);
		
		#ifdef _ITSSB_SHOW_INDICATOR
		{ // #indicator
			TString xstrITSSBIB01221xind("(");			
			for (Int_t i=0; i<fNStaves; i++) // ITSSBIB01221xind
			{
				if (i>0) xstrITSSBIB01221xind += " + ";
				xstrITSSBIB01221xind += "xtubeITSSBIB01221xind";
				xstrITSSBIB01221xind += i;
				xstrITSSBIB01221xind += ":xcombtITSSBIB01221xind";
				xstrITSSBIB01221xind += i;
			}	
			xstrITSSBIB01221xind += ")";

			TGeoCompositeShape *xcompITSSBIB01221xind = new TGeoCompositeShape("xcompITSSBIB01221xind", xstrITSSBIB01221xind);
			TGeoVolume *xvolITSSBIB01221xind = new TGeoVolume("ITSSBIB01221xind", xcompITSSBIB01221xind, med); 			
			xvolITSSBIB01221xind->SetFillColor(2);
			xvolITSSBIB01221xind->SetLineColor(2);
			ewhVol->AddNode(xvolITSSBIB01221xind, 1, ITSSBIB0122CombT);
		}
		#endif
	}

	if (fLayerNumber==2)
	{ // 0123x : Layer 2
		// 01231 = the object
		// 01232 = -------------
		// 01233 = the base
		// 01234 = the tube for holes
		// 01234 = the tube for cut

		Double_t ITSSBIB01231Thick = 0.04;
		Double_t ITSSBIB01233Thick = 0.45; // over approximate, vary this to offset for contact with stave

		Double_t ITSSBIB01231Radius1 = 3.65; // = ITSSBIB01105Radius //
		Double_t ITSSBIB01231Radius2 = 4.40 - ITSSBIB01231Thick; // = ITSSBIB01107Radius-ITSSBIB0123Thick; //

		Double_t ITSSBIB01231Length = 3.60;
		Double_t ITSSBIB01232Length = 1.90;
		Double_t ITSSBIB01233Length = 1.30;
		Double_t ITSSBIB01234Length = ITSSBIB01233Thick*2.5; // 2.5 time is easy to make a hole without offset this tube

		Double_t ITSSBIB0123Length = ITSSBIB01231Length;
		Double_t ITSSBIB01235Length = ITSSBIB0123Length*1.5; // 1.5 time is easy to make tube for div without offset its

		Double_t ITSSBIB01233Width = 1.10; // over approximate, vary this to offset for contact with stave
		Double_t ITSSBIB01231CenterWidth = (ITSSBIB01231Radius2+(ITSSBIB01231Thick/2))-(ITSSBIB01231Radius1-(ITSSBIB01231Thick/2));

		Double_t ITSSBIB01231Radius = 2.85-ITSSBIB01231Thick; // = ITSSBIB01103Radius-ITSSBIB01231Thick; 
		Double_t ITSSBIB01232Radius = ITSSBIB01231Radius + ITSSBIB01231Length - ITSSBIB01231Thick;
		Double_t ITSSBIB01233Radius = 4.11; // 4.107681742, radius for circle, offset at top surface of each object
		Double_t ITSSBIB01234Radius = 0.25; // ***IMPORTANT*** this is radius each hole, for circle we must use ITSSBIB01233Radius ***IMPORTANT***
		Double_t ITSSBIB01235Radius = ITSSBIB01231Radius2 + 0.01; // TO MAKE SURE WE CAN WELDED TO 1111&1112 (A bit should < ITSSBIB01231Radius+ITSSBIB01231Thick

		Double_t ITSSBIB0123Radius = ITSSBIB01231Radius1;

		// calculate at 45 degree
		Double_t ITSSBIB01231Length1 = 0.40; // fixed parameter
		Double_t ITSSBIB01231Length2 = ITSSBIB01231CenterWidth-(ITSSBIB01231Thick/2)-(ITSSBIB01231Thick/2); // followed parameter
		Double_t ITSSBIB01231Length3 = ITSSBIB01231Length-ITSSBIB01231Length1-ITSSBIB01231Length2; // followed parameter

		Double_t ITSSBIB01231Z0 = 0.00; // The center of object is not the same as the center of space.
		Double_t ITSSBIB01231Z1 = -(ITSSBIB01231Length3 - (ITSSBIB01231Thick/2));
		Double_t ITSSBIB01231Z2 = -(ITSSBIB01231Length3 + (ITSSBIB01231Thick/2));
		Double_t ITSSBIB01231Z3 = -(ITSSBIB01231Length3 + ITSSBIB01231Length2 - (ITSSBIB01231Thick/2));
		Double_t ITSSBIB01231Z4 = -(ITSSBIB01231Length3 + ITSSBIB01231Length2 + (ITSSBIB01231Thick/2));
		Double_t ITSSBIB01231Z5 = -(ITSSBIB01231Length);

		Double_t ITSSBIB01231RMin0 = ITSSBIB01231Radius2, ITSSBIB01231RMax0 = ITSSBIB01231Radius2+ITSSBIB01231Thick;
		Double_t ITSSBIB01231RMin1 = ITSSBIB01231Radius2, ITSSBIB01231RMax1 = ITSSBIB01231Radius2+ITSSBIB01231Thick;
		Double_t ITSSBIB01231RMin2 = ITSSBIB01231Radius2 - ITSSBIB01231Thick, ITSSBIB01231RMax2 = ITSSBIB01231Radius2+ITSSBIB01231Thick;
		Double_t ITSSBIB01231RMin3 = ITSSBIB01231Radius1, ITSSBIB01231RMax3 = ITSSBIB01231Radius1+(2*ITSSBIB01231Thick);
		Double_t ITSSBIB01231RMin4 = ITSSBIB01231Radius1, ITSSBIB01231RMax4 = ITSSBIB01231Radius1+ITSSBIB01231Thick;
		Double_t ITSSBIB01231RMin5 = ITSSBIB01231Radius1, ITSSBIB01231RMax5 = ITSSBIB01231Radius1+ITSSBIB01231Thick;

		
		TGeoRotation *xrotITSSBIB0123Circle = new TGeoRotation("xrotITSSBIB0123Circle", 0.0, 0.0,  - 90.0 + 10.02); // this rotate BEFORE div
		xrotITSSBIB0123Circle->RegisterYourself();
		
		Double_t ITSSBIB01233EachRotate = 15.0; // degree of rotate each object of 01233
		Double_t ITSSBIB01234EachRotate = 15.0;
		
		Double_t ITSSBIB0123MaxRadius = ITSSBIB01231Radius + ITSSBIB01231Length;
		
		// offset is compared with volume's shape of section1
		Double_t ITSSBIB0123OffsetX = 0.0;
		Double_t ITSSBIB0123OffsetY = 0.0;
		Double_t ITSSBIB0123OffsetZ = 0.0 - 11.5; // original at z=0, not the center of the object
		TGeoRotation *ITSSBIB0123Rotate = new TGeoRotation("", 0.0, 0.0, 0.0); // this rotate AFTER div
		TGeoCombiTrans *ITSSBIB0123CombT = new TGeoCombiTrans(ITSSBIB0123OffsetX, ITSSBIB0123OffsetY, ITSSBIB0123OffsetZ, ITSSBIB0123Rotate);

		TGeoTube *xtubeITSSBIB01231 = new TGeoTube("xtubeITSSBIB01231"
			, ITSSBIB01231Radius, ITSSBIB01231Radius + ITSSBIB01231Thick, ITSSBIB01231Length/2);

		TGeoPcon *xpconeITSSBIB01231 = new TGeoPcon("xpconeITSSBIB01231", 0, 360, 6);
		xpconeITSSBIB01231->DefineSection(0, ITSSBIB01231Z0, ITSSBIB01231RMin0, ITSSBIB01231RMax0);
		xpconeITSSBIB01231->DefineSection(1, ITSSBIB01231Z1, ITSSBIB01231RMin1, ITSSBIB01231RMax1);
		xpconeITSSBIB01231->DefineSection(2, ITSSBIB01231Z2, ITSSBIB01231RMin2, ITSSBIB01231RMax2);
		xpconeITSSBIB01231->DefineSection(3, ITSSBIB01231Z3, ITSSBIB01231RMin3, ITSSBIB01231RMax3);
		xpconeITSSBIB01231->DefineSection(4, ITSSBIB01231Z4, ITSSBIB01231RMin4, ITSSBIB01231RMax4);
		xpconeITSSBIB01231->DefineSection(5, ITSSBIB01231Z5, ITSSBIB01231RMin5, ITSSBIB01231RMax5);			

		TGeoTranslation *xtranITSSBIB01231 = new TGeoTranslation("xtranITSSBIB01231", 0.0, 0.0, 0.0);
		xtranITSSBIB01231->RegisterYourself();			

		Int_t NDegPerStave = 360/fNStaves;
		Double_t OffNegSideX01233 = ITSSBIB01233Width/2 - 0.5;
		Double_t Deg2RadFactor = TMath::Pi()/180;
		Double_t CosEachRotate = TMath::Cos(ITSSBIB01233EachRotate*Deg2RadFactor);
		Double_t SinEachRotate = TMath::Sin(ITSSBIB01233EachRotate*Deg2RadFactor);
		Double_t FaceCirRadi01233 = ITSSBIB01233Thick/2 + ITSSBIB01233Radius;
		
		for (Int_t i=0; i<fNStaves; i++)
		{
			Double_t CosDegIter = TMath::Cos(i*NDegPerStave*Deg2RadFactor);
			Double_t SinDegIter = TMath::Sin(i*NDegPerStave*Deg2RadFactor);					

			{ // 012313
				TString str("xbboxITSSBIB01233");
				str += i;
				TGeoBBox *xbboxITSSBIB001233x = new TGeoBBox(str, ITSSBIB01233Width/2, ITSSBIB01233Thick/2, ITSSBIB01233Length/2);
				
				TString str2("xcombtITSSBIB01233");
				str2 += i;
				TGeoCombiTrans *xcombtITSSBIB01233x = new TGeoCombiTrans(str2
					, 0 + (FaceCirRadi01233 * SinDegIter)	// circle with top face offset
						+ (OffNegSideX01233 * CosEachRotate) * CosDegIter // x-offset
						- (OffNegSideX01233 * SinEachRotate) * SinDegIter // y-offset
					, 0 - (FaceCirRadi01233 * CosDegIter)	// circle	with top face offset
						+ (OffNegSideX01233 * CosEachRotate) * SinDegIter // x-offset
						+ (OffNegSideX01233 * SinEachRotate) * CosDegIter // y-offset  (*** ??? why positive, not negative ??? ***)
					, 0 - (ITSSBIB01233Length/2) 	// set z to the head of parent object
					
					, new TGeoRotation("", 0.0, 0.0, i*NDegPerStave + ITSSBIB01233EachRotate));
				
				xcombtITSSBIB01233x->RegisterYourself();
			}
			
			{ // 012314
				TString str("xtubeITSSBIB01234");
				str += i;
				TGeoTube *xtubeITSSBIB001234x = new TGeoTube(str, 0.0, ITSSBIB01234Radius, ITSSBIB01234Length/2);
				
				TString str2("xcombtITSSBIB01234");
				str2 += i;
				TGeoCombiTrans *xcombtITSSBIB01234x = new TGeoCombiTrans(str2
					, 0 + (ITSSBIB01233Radius * SinDegIter)	// circle *** USE ITSSBIB01233Radius ***
					, 0 - (ITSSBIB01233Radius * CosDegIter)	// circle	*** USE ITSSBIB01233Radius ***
					, 0 - .40 		// set z to the head of parent object
					// this have been set the position to circle around tube, also offset to the top face(Y-POS) *** BUT NOT X-POS ***
					, new TGeoRotation("", i*NDegPerStave + ITSSBIB01234EachRotate, 90.0, 0.0)); // perpendicular to 01233
				
				xcombtITSSBIB01234x->RegisterYourself();
			}					
			
			#ifdef _ITSSB_SHOW_INDICATOR
			{ // #indicator
				TString str("xtubeITSSBIB01231xind");
				str += i;
				TGeoTube *xtubeITSSBIB001231xind = new TGeoTube(str, 0.0, 0.01, ITSSBIB01233Length);
				
				TString str2("xcombtITSSBIB01231xind");
				str2 += i;
				TGeoCombiTrans *xcombtITSSBIB01231xind = new TGeoCombiTrans(str2
					, 0 + (ITSSBIB01233Radius * SinDegIter)	// circle *** USE ITSSBIB01233Radius ***
					, 0 - (ITSSBIB01233Radius * CosDegIter)	// circle	*** USE ITSSBIB01233Radius ***
					, 0		// set z to the head of parent object
					// this have been set the position to circle around tube, also offset to the top face(Y-POS) *** BUT NOT X-POS ***
					, new TGeoRotation("", 0.0, 0.0, 0.0)); 
				
				xcombtITSSBIB01231xind->RegisterYourself();
			}
			#endif					
			
		}

		TGeoTube *xtubeITSSBIB01235 = new TGeoTube("xtubeITSSBIB01235", 0.0, ITSSBIB01235Radius, ITSSBIB01235Length/2);
		TGeoTranslation *xtranITSSBIB01235 = new TGeoTranslation("xtranITSSBIB01235", 0.0, 0.0, 0.0);
		xtranITSSBIB01235->RegisterYourself();			


		// combine all 0123x objects
		TString xstrITSSBIB0123cs("(");				
		TString xstrITSSBIB01234trunks("("); // collect 14, all holes to div
		for (Int_t i=0; i<fNStaves; i++) // ITSSBIB01233
		{ // *** IF PLACE CODE 114 BELOW 113, ITS POSITION MIGHT BE UNDETERMINED ***
			
			TString xstrITSSBIB01234trunk("xtubeITSSBIB01234");
			xstrITSSBIB01234trunk += i;
			xstrITSSBIB01234trunk += ":xcombtITSSBIB01234";
			xstrITSSBIB01234trunk += i;
			
			if (i>0) xstrITSSBIB01234trunks += " + ";
			xstrITSSBIB01234trunks += xstrITSSBIB01234trunk;

			if (i>0) xstrITSSBIB0123cs += " + ";
			xstrITSSBIB0123cs += "(xbboxITSSBIB01233";
			xstrITSSBIB0123cs += i;
			xstrITSSBIB0123cs += ":xcombtITSSBIB01233";
			xstrITSSBIB0123cs += i; 
			
			#ifdef _ITSSB_SHOW_COMBISHAPE
			xstrITSSBIB0123cs += "*xtubeITSSBIB01235:xtranITSSBIB01235 - ";	// cut & make a hole
			//xstrITSSBIB0123cs += " - ";	// cut & make a hole
			#else
			xstrITSSBIB0123cs += " + ";		// just show
			#endif
			xstrITSSBIB0123cs += xstrITSSBIB01234trunk;												
			
			xstrITSSBIB0123cs += ")"; 
			
		}
		xstrITSSBIB01234trunks += ")";
		//xstrITSSBIB0123cs += " + xtubeITSSBIB01235:xtranITSSBIB01235";	// cut & make a hole

		// 01231
		xstrITSSBIB0123cs += " + (xpconeITSSBIB01231:xtranITSSBIB01231";
		#ifdef _ITSSB_SHOW_COMBISHAPE
		xstrITSSBIB0123cs += " - ";
		xstrITSSBIB0123cs += xstrITSSBIB01234trunks;	// make holes
		#endif
		xstrITSSBIB0123cs += ")";
		
														
		xstrITSSBIB0123cs += "):xrotITSSBIB0123Circle";
		#ifdef _ITSSB_SHOW_WITHBOXDIV
		xstrITSSBIB0123cs += " * (xbboxdiv:xtrandiv)";
		#endif 
		
		TGeoCompositeShape *xcompITSSBIB0123 = new TGeoCompositeShape("xcompITSSBIB0123cs", xstrITSSBIB0123cs);
		TGeoVolume *xvolITSSBIB0123 = new TGeoVolume("ITSSBIB0123x", xcompITSSBIB0123, med); 			
		xvolITSSBIB0123->SetFillColor(9);
		xvolITSSBIB0123->SetLineColor(9);
		ewhVol->AddNode(xvolITSSBIB0123, 1, ITSSBIB0123CombT);
		
		#ifdef _ITSSB_SHOW_INDICATOR
		{ // #indicator
			TString xstrITSSBIB01231xind("(");			
			for (Int_t i=0; i<fNStaves; i++) // ITSSBIB01231xind
			{
				if (i>0) xstrITSSBIB01231xind += " + ";
				xstrITSSBIB01231xind += "xtubeITSSBIB01231xind";
				xstrITSSBIB01231xind += i;
				xstrITSSBIB01231xind += ":xcombtITSSBIB01231xind";
				xstrITSSBIB01231xind += i;
			}	
			xstrITSSBIB01231xind += ")";

			TGeoCompositeShape *xcompITSSBIB01231xind = new TGeoCompositeShape("xcompITSSBIB01231xind", xstrITSSBIB01231xind);
			TGeoVolume *xvolITSSBIB01231xind = new TGeoVolume("ITSSBIB01231xind", xcompITSSBIB01231xind, med); 			
			xvolITSSBIB01231xind->SetFillColor(2);
			xvolITSSBIB01231xind->SetLineColor(2);
			ewhVol->AddNode(xvolITSSBIB01231xind, 1, ITSSBIB0123CombT);
		}
		#endif
	}



	return ewhVol;
}



//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateInnerBEWheelA1(const TGeoMedium* med){
	//
	//
	//

	TGeoVolume* ewhVol = new TGeoVolumeAssembly(TString("EndWheelA") += fLayerNumber);

	const Double_t fkPlateThick = 0.4*fgkmm;
	
	const Double_t fkHndLength = 13.0*fgkmm; // also use in CreateInnerBHandleA.
	
	// Layer 0
	const Double_t fkPlateL0Radius = 28.1*fgkmm;
	const Double_t fkPlateL0Length = 57.0*fgkmm;
	const Double_t fkHndPlateL0Radius = 27.9*fgkmm;
	const Double_t fkContactL0Radius = 25.9*fgkmm; // approx : -2mm from Radius
	
	// Layer 1
	const Double_t fkPlateL1Radius1 = fkPlateL0Radius + fkPlateThick;
	const Double_t fkPlateL1Radius2 = 36.1*fgkmm;
	const Double_t fkPlateL1Length1 = 52.0*fgkmm - fkHndLength;
	const Double_t fkPlateL1Length2 = 3.4*fgkmm;
	const Double_t fkHndPlateL1Radius = 35.9*fgkmm;
	const Double_t fkContactL1Radius = 33.8*fgkmm;
	
	// Layer 2
	const Double_t fkPlateL2Radius1 = fkPlateL1Radius2 + fkPlateThick;
	const Double_t fkPlateL2Radius2 = 43.6*fgkmm;
	const Double_t fkPlateL2Length1 = 36.0*fgkmm - fkHndLength;
	const Double_t fkPlateL2Length2 = 4.0*fgkmm;
	const Double_t fkHndPlateL2Radius = 43.4*fgkmm;
	const Double_t fkContactL2Radius = 41.1*fgkmm;
	
	if (fLayerNumber==0)
	{
		const Double_t fkHndWidthAngle = 360/fNStaves;

		
		ewhVol->AddNode(CreateTube("EndWheelObj[1]", fkPlateL0Radius, fkPlateL0Length-fkHndLength, fkPlateThick, kFALSE, 6, 180, med), 1, 
			new TGeoCombiTrans(0, 0, 0-fkHndLength, new TGeoRotation("", 0, 0, 180)));

		for (Int_t i=0; i<fNStaves/2; i++)
		ewhVol->AddNode(CreateInnerBHandleA("", fkHndPlateL0Radius, fkContactL0Radius, med), 1, 
			new TGeoCombiTrans(0, 0, 0, new TGeoRotation("", 0, 0, 180 +fkHndWidthAngle/2 +i*fkHndWidthAngle)));
	}
	
	if (fLayerNumber==1)
	{
		const Double_t fkHndWidthAngle = 360/fNStaves;

		ewhVol->AddNode(CreateSTube45("EndWheelObj[1]", fkPlateL1Radius1, fkPlateL1Radius2, fkPlateL1Length1, fkPlateL1Length2, fkPlateThick, kFALSE, 5, med), 1, 
			new TGeoCombiTrans(0, 0, 0-fkHndLength, new TGeoRotation("", 0, 0, 180)));

		for (Int_t i=0; i<fNStaves/2; i++)
		ewhVol->AddNode(CreateInnerBHandleA("", fkHndPlateL1Radius, fkContactL1Radius, med), 1, 
			new TGeoCombiTrans(0, 0, 0, new TGeoRotation("", 0, 0, 180+fkHndWidthAngle/2+i*fkHndWidthAngle)));
	}

	
	if (fLayerNumber==2)
	{
		const Double_t fkHndWidthAngle = 360/fNStaves;

		ewhVol->AddNode(CreateSTube45("EndWheelObj[1]", fkPlateL2Radius1, fkPlateL2Radius2, fkPlateL2Length1, fkPlateL2Length2, fkPlateThick, kFALSE, 7, med), 1, 
			new TGeoCombiTrans(0, 0, 0-fkHndLength, new TGeoRotation("", 0, 0, 180)));

		for (Int_t i=0; i<fNStaves/2; i++)
		ewhVol->AddNode(CreateInnerBHandleA("", fkHndPlateL2Radius, fkContactL2Radius, med), 1, 
			new TGeoCombiTrans(0, 0, 0, new TGeoRotation("", 0, 0, 180+fkHndWidthAngle/2+i*fkHndWidthAngle)));

		// Supporter (A-Side)
		ewhVol->AddNode(CreateInnerBSupporterA(), 1, new TGeoTranslation(0.,0.,0.-(400./2-5.-65.)*fgkmm));
	}
	

	return ewhVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateInnerBSupporterA(const TGeoManager *mgr){
//
//
//
	const Double_t kRadius = 44.0*fgkmm;
	const Double_t kThick1 = 3.8*fgkmm;
	const Double_t kLength1 = 2.0*fgkmm;
	const Double_t kThick2 = 0.4*fgkmm;
	const Double_t kLength2 = 3.0*fgkmm;

	TGeoMedium *medium = mgr->GetMedium("ITS_PEEKCF30$");

	TGeoVolume* suppVol = new TGeoVolumeAssembly("InnerBSupporterA");
		
	TGeoPcon *xpCone = new TGeoPcon("", 0, 180.0, 4);
		xpCone->DefineSection(0, 0.00, kRadius, kRadius + kThick1);
		xpCone->DefineSection(1, -kLength1, kRadius, kRadius + kThick1);
		xpCone->DefineSection(2, -kLength1, kRadius, kRadius + kThick2);
		xpCone->DefineSection(3, -(kLength1 + kLength2), kRadius, kRadius + kThick2);

	TGeoVolume *medVol = new TGeoVolume("InnerBSupporterA", xpCone, medium); 			
	medVol->SetFillColor(7);
	medVol->SetLineColor(7);
	suppVol->AddNode(medVol, 0, new TGeoCombiTrans( 0, 0, 0, new TGeoRotation("",0,0,180.0)));

	return suppVol;

}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateInnerBSupporterC(const TGeoManager *mgr){
//
//
//
	TGeoMedium *medium = mgr->GetMedium("ITS_PEEKCF30$");

	const Double_t kThick  = 1.00 *fgkmm;
	const Double_t kLength = 10.0 *fgkmm;
	const Double_t kRadius = 22.0 *fgkmm;
	const Double_t kMaxRadius = 47.8 *fgkmm;
											
	const Double_t kPosL1R = 28.5 *fgkmm;
	const Double_t kPosL1Z = 4.00 *fgkmm;
	const Double_t kPosL2R = 36.5 *fgkmm;
	const Double_t kPosL2Z = 6.00 *fgkmm;
	const Double_t kPosL3R = 44.0 *fgkmm;
	const Double_t kPosL3Z = 8.00 *fgkmm;
	
	TGeoVolume* suppVol = new TGeoVolumeAssembly("InnerBSupporterC");

	TGeoPcon *xpCone = new TGeoPcon("", 0, 180, 16);
		xpCone->DefineSection(0, 	0.00				, kRadius	, kRadius + kThick);
		xpCone->DefineSection(1, 	-(kPosL1Z - kThick)	, kRadius	, kRadius + kThick);
		xpCone->DefineSection(2, 	-(kPosL1Z - kThick)	, kRadius	, kPosL1R + kThick);
		xpCone->DefineSection(3, 	-kPosL1Z			, kRadius	, kPosL1R + kThick);
		xpCone->DefineSection(4, 	-kPosL1Z			, kPosL1R	, kPosL1R + kThick);
		xpCone->DefineSection(5, 	-(kPosL2Z - kThick)	, kPosL1R	, kPosL1R + kThick);
		xpCone->DefineSection(6, 	-(kPosL2Z - kThick)	, kPosL1R	, kPosL2R + kThick);
		xpCone->DefineSection(7, 	-kPosL2Z			, kPosL1R	, kPosL2R + kThick);
		xpCone->DefineSection(8, 	-kPosL2Z			, kPosL2R	, kPosL2R + kThick);
		xpCone->DefineSection(9, 	-(kPosL3Z - kThick)	, kPosL2R	, kPosL2R + kThick);
		xpCone->DefineSection(10, 	-(kPosL3Z - kThick)	, kPosL2R	, kPosL3R + kThick);
		xpCone->DefineSection(11, 	-kPosL3Z			, kPosL2R	, kPosL3R + kThick);
		xpCone->DefineSection(12, 	-kPosL3Z			, kPosL3R	, kPosL3R + kThick);
		xpCone->DefineSection(13, 	-(kLength - kThick)	, kPosL3R	, kPosL3R + kThick);
		xpCone->DefineSection(14, 	-(kLength - kThick)	, kPosL3R	, kMaxRadius);
		xpCone->DefineSection(15, 	-kLength			, kPosL3R	, kMaxRadius);

	TGeoVolume *medVol = new TGeoVolume("InnerBSupporterCLadder", xpCone, medium); 		
	medVol->SetFillColor(6);
	medVol->SetLineColor(6);
	suppVol->AddNode(medVol, 0, new TGeoCombiTrans( 0, 0, 0, new TGeoRotation("",0,0,180.0)));
	
	xpCone = new TGeoPcon("", 0, 180, 2);
		xpCone->DefineSection(0, 0.0					, kMaxRadius - kThick	, kMaxRadius);
		xpCone->DefineSection(1, -(kLength - kThick)	, kMaxRadius - kThick	, kMaxRadius);
	medVol = new TGeoVolume("InnerBSupporterCBase", xpCone, medium); 			
	medVol->SetFillColor(6);
	medVol->SetLineColor(6);
	suppVol->AddNode(medVol, 1, new TGeoCombiTrans( 0, 0, 0, new TGeoRotation("",0,0,180.0)));


	return suppVol;

}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateInnerBCShell(const TGeoManager *mgr){

	const Double_t kRadius =   47.8*fgkmm;
	//const Double_t kLength =  400.0*fgkmm;
	const Double_t kLength =  347.5*fgkmm;
	const Double_t kSkinThick = 0.1*fgkmm;
	const Double_t kCoreThick = 1.8*fgkmm;

	TGeoMedium* skinMedium = mgr->GetMedium("ITS_PEEKCF30$");
	TGeoMedium* coreMedium = mgr->GetMedium("ITS_CarbonFoam$");

	TGeoVolume* shellVol = new TGeoVolumeAssembly("InnerBCShell");
	TGeoCombiTrans* combTrans = new TGeoCombiTrans( 0, 0, kLength/2 + (400.0*fgkmm - kLength)/2, new TGeoRotation("",0,0,180.0));

	TGeoPcon* xpCone = new TGeoPcon("", 0, 180.0, 2);
		xpCone->DefineSection(0, 0.00, kRadius, kRadius + kSkinThick);
		xpCone->DefineSection(1, -kLength, kRadius, kRadius + kSkinThick);
	TGeoVolume* volume = new TGeoVolume("InnerBCShellInside", xpCone, skinMedium);
	volume->SetFillColor(5);
	volume->SetLineColor(5);
	shellVol->AddNode(volume, 0, combTrans);
	
	xpCone = new TGeoPcon("", 0, 180.0, 2);
		xpCone->DefineSection(0, 0.00, kRadius + kSkinThick, kRadius + kSkinThick + kCoreThick);
		xpCone->DefineSection(1, -kLength, kRadius + kSkinThick, kRadius + kSkinThick + kCoreThick);
	volume = new TGeoVolume("InnerBCShellFoam", xpCone, coreMedium); 			
	volume->SetFillColor(41);
	volume->SetLineColor(41);
	shellVol->AddNode(volume, 1, combTrans);
	
	xpCone = new TGeoPcon("", 0, 180.0, 2);
		xpCone->DefineSection(0, 0.00, kRadius + kSkinThick + kCoreThick, kRadius + 2*kSkinThick + kCoreThick);
		xpCone->DefineSection(1, -kLength, kRadius + kSkinThick + kCoreThick, kRadius + 2*kSkinThick + kCoreThick);
	volume = new TGeoVolume("InnerBCShellOutside", xpCone, skinMedium); 			
	volume->SetFillColor(5);
	volume->SetLineColor(5);
	shellVol->AddNode(volume, 2, combTrans);

	return shellVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateOuterBCShell(const TGeoManager *mgr){

	const Double_t kRadius = 44.9*fgkcm;
	const Double_t kLength = 174.0*fgkcm;
	const Double_t kSkinThick = 0.6*fgkmm;
	const Double_t kCoreThick = 8.8*fgkmm;

	TGeoMedium* skinMedium = mgr->GetMedium("ITS_PEEKCF30$");
	TGeoMedium* coreMedium = mgr->GetMedium("ITS_CarbonFoam$");

	TGeoVolume* shellVol = new TGeoVolumeAssembly("OuterBCShell");
	TGeoCombiTrans* combTrans = new TGeoCombiTrans( 0, 0, kLength/2, new TGeoRotation("",0,0,180.0));

	TGeoPcon* xpCone = new TGeoPcon("", 0, 180.0, 2);
		xpCone->DefineSection(0, 0.00, kRadius, kRadius + kSkinThick);
		xpCone->DefineSection(1, -kLength, kRadius, kRadius + kSkinThick);
	TGeoVolume* volume = new TGeoVolume("OuterBCShellInside", xpCone, skinMedium);
	volume->SetFillColor(3);
	volume->SetLineColor(3);
	shellVol->AddNode(volume, 0, combTrans);
	
	xpCone = new TGeoPcon("", 0, 180.0, 2);
		xpCone->DefineSection(0, 0.00, kRadius + kSkinThick, kRadius + kSkinThick + kCoreThick);
		xpCone->DefineSection(1, -kLength, kRadius + kSkinThick, kRadius + kSkinThick + kCoreThick);
	volume = new TGeoVolume("OuterBCShellFoam", xpCone, coreMedium); 			
	volume->SetFillColor(5);
	volume->SetLineColor(5);
	shellVol->AddNode(volume, 1, combTrans);
	
	xpCone = new TGeoPcon("", 0, 180.0, 2);
		xpCone->DefineSection(0, 0.00, kRadius + kSkinThick + kCoreThick, kRadius + 2*kSkinThick + kCoreThick);
		xpCone->DefineSection(1, -kLength, kRadius + kSkinThick + kCoreThick, kRadius + 2*kSkinThick + kCoreThick);
	volume = new TGeoVolume("OuterBCShellOutside", xpCone, skinMedium); 			
	volume->SetFillColor(3);
	volume->SetLineColor(3);
	shellVol->AddNode(volume, 2, combTrans);

	return shellVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateInnerBServiceB(const TGeoManager *mgr){


	TGeoVolume* svbVol = new TGeoVolumeAssembly("InnerBServiceB");


	const Double_t ITSSBIB02OffsetX = 0.0;
	const Double_t ITSSBIB02OffsetY = 0.0;
	const Double_t ITSSBIB02OffsetZ = 0.0; 
	
	const Double_t ITSSBIB02RotateAlpha = fgkBarrelRotateAlpha;
	const Double_t ITSSBIB02RotateBetha = fgkBarrelRotateBetha;
	const Double_t ITSSBIB02RotateGamma = fgkBarrelRotateGamma;

	const Double_t ITSSBIB02CutOffset1 = 10.0 *fgkmm;
	const Double_t ITSSBIB02CutOffset2 = 0.00;
	const Double_t ITSSBIB02CutOffset3 = 16.00 *fgkmm;	
	const Double_t ITSSBIB02CutOffset4 = 15.50 *fgkmm;
	const Double_t ITSSBIB02CutOffset5 = 10.50 *fgkmm;		
	const Double_t ITSSBIB02CutOffset6 = 10.00 *fgkmm;		
	const Double_t ITSSBIB02Length = 270.0 *fgkcm;
	const Double_t ITSSBIB02Radius = 4.78; //ITSSBIB01Radius
	
	// Section 01
	const Double_t ITSSBIB02SWSkinThick1 = 0.10 *fgkmm;
	const Double_t ITSSBIB02SWCoreThick1 = 0.00 *fgkmm;
	const Double_t ITSSBIB02SWSkinThick2 = 0.50 *fgkmm;
	const Double_t ITSSBIB02SWCoreThick2 = 0.00 *fgkmm;
	
	TGeoMedium *MediumWall = mgr->GetMedium("ITS_KAPTON(POLYCH2)$");	
	TGeoMedium *MediumFoam = mgr->GetMedium("ITS_KAPTON(POLYCH2)$");	
		
	// 0200x
	if (true)
	{ 	
		//---------------------------- box div --------------------------
		//---------------------- GLOBAL PARAMETER------------------------
		// box to div 0200xa
		TGeoBBox *xdiv0200xa = new TGeoBBox("xdiv0200xa", 110.0*fgkmm, 110.0*fgkmm, (374.3252*fgkmm)/2); 
		(new TGeoTranslation("xtran0200xa", 0.0, 0.0 + xdiv0200xa->GetDY() - ITSSBIB02CutOffset1, ITSSBIB02OffsetZ - xdiv0200xa->GetDZ()))->RegisterYourself();

		// box to div 0200xb : TEMPORARY
		TGeoBBox *xdiv0200xb = new TGeoBBox("xdiv0200xb", 240.0*fgkmm, 240.0*fgkmm, ((574.7274 - 374.3252 + 2.0)*fgkmm)/2); 
		(new TGeoTranslation("xtran0200xb", 0.0, 0.0 + xdiv0200xb->GetDY() - ITSSBIB02CutOffset1, ITSSBIB02OffsetZ - xdiv0200xb->GetDZ() - 372.3252*fgkmm - 1.0*fgkmm))->RegisterYourself();

		// box to div 0200xc
		TGeoBBox *xdiv0200xc = new TGeoBBox("xdiv0200xc", 500.0*fgkmm, 500.0*fgkmm, (152.2762*fgkmm)/2);
		(new TGeoTranslation("xtran0200xc", 0.0, 0.0 + xdiv0200xc->GetDY() - ITSSBIB02CutOffset3, ITSSBIB02OffsetZ - xdiv0200xc->GetDZ() - 574.7274*fgkmm))->RegisterYourself(); // OffsetZ = --02001:Pos5

		// box to div 0200xd : depend on obj
		// box to div 0200xe : depend on obj
		// box to div 0200xf : depend on obj

		// box to div 0200xt : cut a little bit of the tail (RAPID MODEL)
		TGeoBBox *xdiv0200xt = new TGeoBBox("xdiv0200xt", 500.0*fgkmm, 500.0*fgkmm, (10.0*fgkmm)/2);
		(new TGeoTranslation("xtran0200xt", 0.0, 0.0, ITSSBIB02OffsetZ - xdiv0200xt->GetDZ() - 2700.0*fgkmm + 0.5*fgkmm))->RegisterYourself();
		

		//----------------------------------------------------------------

		// for div with the wing of 02001
		TGeoPcon *xpcone02000s = new TGeoPcon("xpconeITSSBIB02000s", 0, 360, 2);
		xpcone02000s->DefineSection(0,  -576.8138*fgkmm,  0.0, 227.40995840*fgkmm);
		xpcone02000s->DefineSection(1, -2700.0000*fgkmm,  0.0, 473.00390000*fgkmm);
		(new TGeoTranslation("xtranITSSBIB02000s", ITSSBIB02OffsetX, ITSSBIB02OffsetY, ITSSBIB02OffsetZ))->RegisterYourself();
		
		// box to div 02000e : crop the tail
		TGeoBBox *xdiv02000e = new TGeoBBox("xdiv02000e", ((879.0 - 0.5*2 -5.0*2)*fgkmm)/2, 500.0*fgkmm, (2700.0+1.0*fgkmm)/2);
		(new TGeoTranslation("xtran02000e", 0.0, 0.0, ITSSBIB02OffsetZ - xdiv02000e->GetDZ() - 1.0*fgkmm))->RegisterYourself(); // OffsetZ = to the tail

		//-------------------------------------------------------------------------------------------------
		// 02001 - The inner
		
		TGeoPcon *xpcone02001 = new TGeoPcon("xpconeITSSBIB02001", 0, 360, 8);
		xpcone02001->DefineSection(0,    -0.0000*fgkmm,  47.80000000*fgkmm,  47.90000000*fgkmm);
		xpcone02001->DefineSection(1,   -74.5463*fgkmm,  47.80000000*fgkmm,  47.90000000*fgkmm);
		xpcone02001->DefineSection(2,   -77.0018*fgkmm,  47.80000000*fgkmm,  48.30815396*fgkmm);
		xpcone02001->DefineSection(3,  -376.3568*fgkmm,  97.35912424*fgkmm,  97.86584406*fgkmm);
		xpcone02001->DefineSection(4,  -376.5415*fgkmm,  97.38916457*fgkmm,  97.98531620*fgkmm);
		xpcone02001->DefineSection(5,  -576.6399*fgkmm, 227.29703910*fgkmm, 227.89316060*fgkmm);
		xpcone02001->DefineSection(6,  -576.8138*fgkmm, 227.40995840*fgkmm, 227.91578880*fgkmm);
		xpcone02001->DefineSection(7, -2700.0000*fgkmm, 472.50000000*fgkmm, 473.00390000*fgkmm);
		(new TGeoTranslation("xtranITSSBIB02001", ITSSBIB02OffsetX, ITSSBIB02OffsetY, ITSSBIB02OffsetZ))->RegisterYourself();

		// for div with the wing of 02002
		TGeoPcon *xpcone02001s = new TGeoPcon("xpconeITSSBIB02001s", 0, 360, 2);
		xpcone02001s->DefineSection(0,  -576.8138*fgkmm,  0.0, 227.91578880*fgkmm);
		xpcone02001s->DefineSection(1, -2700.0000*fgkmm,  0.0, 473.00390000*fgkmm);
		(new TGeoTranslation("xtranITSSBIB02001s", ITSSBIB02OffsetX, ITSSBIB02OffsetY, ITSSBIB02OffsetZ))->RegisterYourself();

		//---------------------------- box div --------------------------

		// box to div 02001d
		TGeoBBox *xdiv02001d = new TGeoBBox("xdiv02001d", 500.0*fgkmm, 500.0*fgkmm, (2000.0*fgkmm)/2);
		(new TGeoTranslation("xtran02001d", 0.0, 0.0 + xdiv02001d->GetDY() - ITSSBIB02CutOffset6, ITSSBIB02OffsetZ + xdiv02001d->GetDZ() - (2700+1.0)*fgkmm))->RegisterYourself(); // OffsetZ = to the tail
		
		// box to div 02001e : crop the tail
		TGeoBBox *xdiv02001e = new TGeoBBox("xdiv02001e", ((879.0 - 0.5*2 - 5.0*2)*fgkmm)/2, 500.0*fgkmm, ((2700.0+1.0)*fgkmm)/2);
		(new TGeoTranslation("xtran02001e", 0.0, 0.0, ITSSBIB02OffsetZ - xdiv02001e->GetDZ() - 1.0*fgkmm))->RegisterYourself(); // OffsetZ = to the tail
		
		// box to div 02001f : crop the paste of the tail
		TGeoPcon *xdiv02001f = new TGeoPcon("xdiv02001f", 0, 360, 2);
		xdiv02001f->DefineSection(0,  -574.9012*fgkmm,  226.76805760*fgkmm,  (226.76805760+50.0)*fgkmm);
		xdiv02001f->DefineSection(1, -2700.0000*fgkmm,  473.00390000*fgkmm,  (473.00390000+50.0)*fgkmm);
		(new TGeoTranslation("xtran02001f", 0.0, 0.0, ITSSBIB02OffsetZ))->RegisterYourself(); // OffsetZ = to the tail
		
					
		//---------------------------- box add --------------------------								
		// box to add 02001a : paste the tail
		TGeoBBox *xadd02001a = new TGeoBBox("xadd02001a", (0.5*fgkmm)/2, 100*fgkmm, (350*fgkmm)/2);
		(new TGeoTranslation("xtranadd02001aR", 0.0 - xadd02001a->GetDX() + ((879.0 - 0.5*2 - 5.0*2)*fgkmm)/2
			, 0.0 - xadd02001a->GetDY(), ITSSBIB02OffsetZ + xadd02001a->GetDZ() - 2700*fgkmm))->RegisterYourself();
		(new TGeoTranslation("xtranadd02001aL", 0.0 + xadd02001a->GetDX() - ((879.0 - 0.5*2 - 5.0*2)*fgkmm)/2
			, 0.0 - xadd02001a->GetDY(), ITSSBIB02OffsetZ + xadd02001a->GetDZ() - 2700*fgkmm))->RegisterYourself();

		// box to add 02001b : roller handle
		TGeoBBox *xadd02001b = new TGeoBBox("xadd02001b", (0.5*fgkmm)/2, 53*fgkmm, ((2700.0-(574.7274+152.2762))*fgkmm)/2); // -(wing length in Z)
		(new TGeoTranslation("xtranadd02001bR", 0.0 - xadd02001b->GetDX() + ((917.0 + 9.0*2 + 0.5*2)*fgkmm)/2
			, 0.0 - xadd02001b->GetDY(), ITSSBIB02OffsetZ + xadd02001b->GetDZ() - 2700*fgkmm))->RegisterYourself();
		(new TGeoTranslation("xtranadd02001bL", 0.0 + xadd02001b->GetDX() - ((917.0 + 9.0*2 + 0.5*2)*fgkmm)/2
			, 0.0 - xadd02001b->GetDY(), ITSSBIB02OffsetZ + xadd02001b->GetDZ() - 2700*fgkmm))->RegisterYourself();
		
		// box to add 02001c : the wing
		TGeoBBox *xadd02001c = new TGeoBBox("xadd02001c", ((917.0 + 9.0*2 + 0.5*2)*fgkmm)/2, ((5.0+0.5)*fgkmm)/2
			, ((2700.000-(574.7274+152.2762))*fgkmm)/2); // -(wing length in Z)
		(new TGeoTranslation("xtranadd02001c", 0.0, 0.0 - xadd02001c->GetDY() - ITSSBIB02CutOffset6
			, ITSSBIB02OffsetZ + xadd02001c->GetDZ() - 2700*fgkmm))->RegisterYourself(); // OffsetZ = to the tail

		
		//#ifdef _ITSSB_SHOW_WITHBOXDIV
		{
			TString sobj("(");
				sobj += "(";
					sobj += "(";
						sobj += "(xpconeITSSBIB02001:xtranITSSBIB02001 * xdiv02001e:xtran02001e)";
						sobj += " + ((xadd02001a:xtranadd02001aL + xadd02001a:xtranadd02001aR) - xdiv02001f:xtran02001f)";
						sobj += " + (xadd02001b:xtranadd02001bL + xadd02001b:xtranadd02001bR)";
					sobj += ")";
					sobj += " - (xdiv0200xa:xtran0200xa + xdiv0200xb:xtran0200xb + xdiv0200xc:xtran0200xc + xdiv02001d:xtran02001d)";
				sobj += ")";						
				sobj += " + (xadd02001c:xtranadd02001c - (xpconeITSSBIB02000s:xtranITSSBIB02000s * xdiv02000e:xtran02000e) )";
			sobj += ")";
			//sobj += " - (xdiv0200xt:xtran0200xt)";

		
			TGeoCompositeShape *xcomp = new TGeoCompositeShape("", sobj);
			TGeoVolume *xvol = new TGeoVolume("ITSSBIB02001", xcomp, MediumWall);
			xvol->SetFillColor(6);
			xvol->SetLineColor(6);				
			svbVol->AddNode(xvol, 1, 0);

			// show wing
			//xvol = new TGeoVolume("ITSSBIB02001xxx", xadd02001c, MediumWall);
			//xvol->SetFillColor(6);
			//xvol->SetLineColor(6);				
			//svbVol->AddNode(xvol, 1, new TGeoTranslation(0.0, 0.0 - xadd02001c->GetDY() - ITSSBIB02CutOffset6, ITSSBIB02OffsetZ + xadd02001c->GetDZ() - 2700*fgkmm));				

		}

		//#else
		//#endif 

		
		//-------------------------------------------------------------------------------------------------
		// 02002 - The foam
		
		TGeoPcon *xpcone02002 = new TGeoPcon("xpconeITSSBIB02002", 0, 360, 8);
		xpcone02002->DefineSection(0,    -0.0000*fgkmm,  47.90000000*fgkmm,  49.70000000*fgkmm);
		xpcone02002->DefineSection(1,   -54.8037*fgkmm,  47.90000000*fgkmm,  49.70000000*fgkmm);
		xpcone02002->DefineSection(2,   -74.5463*fgkmm,  47.90000000*fgkmm,  52.96979799*fgkmm);
		xpcone02002->DefineSection(3,  -374.5100*fgkmm,  97.56015756*fgkmm, 102.62775960*fgkmm);
		xpcone02002->DefineSection(4,  -376.3568*fgkmm,  97.86584406*fgkmm, 103.82729100*fgkmm);
		xpcone02002->DefineSection(5,  -574.9012*fgkmm, 226.76805760*fgkmm, 232.72565530*fgkmm);
		xpcone02002->DefineSection(6,  -576.6399*fgkmm, 227.89316060*fgkmm, 232.92637960*fgkmm);
		xpcone02002->DefineSection(7, -2700.0000*fgkmm, 473.00390000*fgkmm, 478.03710000*fgkmm);
		(new TGeoTranslation("xtranITSSBIB02002", ITSSBIB02OffsetX, ITSSBIB02OffsetY, ITSSBIB02OffsetZ))->RegisterYourself();

		// for div with the wing of 02003
		TGeoPcon *xpcone02002s = new TGeoPcon("xpconeITSSBIB02002s", 0, 360, 2);
		xpcone02002s->DefineSection(0,  -576.6399*fgkmm,  0.0, 232.92637960*fgkmm);
		xpcone02002s->DefineSection(1, -2700.0000*fgkmm,  0.0, 478.03710000*fgkmm);
		(new TGeoTranslation("xtranITSSBIB02002s", ITSSBIB02OffsetX, ITSSBIB02OffsetY, ITSSBIB02OffsetZ))->RegisterYourself();

		//---------------------------- box div --------------------------

		// box to div 02002d
		TGeoBBox *xdiv02002d = new TGeoBBox("xdiv02002d", 500.0*fgkmm, 500.0*fgkmm, (2000.0*fgkmm)/2);
		(new TGeoTranslation("xtran02002d", 0.0, 0.0 + xdiv02002d->GetDY() - ITSSBIB02CutOffset5, ITSSBIB02OffsetZ + xdiv02002d->GetDZ() - (2700+1.0)*fgkmm))->RegisterYourself(); // OffsetZ = to the tail
	
		// box to div 02002e : crop the tail
		TGeoBBox *xdiv02002e = new TGeoBBox("xdiv02002e", ((879.0 - 0.5*2)*fgkmm)/2, 500.0*fgkmm, ((2700.0+1.0)*fgkmm)/2);
		(new TGeoTranslation("xtran02002e", 0.0, 0.0, ITSSBIB02OffsetZ - xdiv02002e->GetDZ() - 1.0*fgkmm))->RegisterYourself(); // OffsetZ = to the tail
		
		// box to div 02002f : crop the paste of the tail
		TGeoPcon *xdiv02002f = new TGeoPcon("xdiv02002f", 0, 360, 2);
		xdiv02002f->DefineSection(0,  -574.7274*fgkmm,  232.72565530*fgkmm,  (232.72565530+50.0)*fgkmm);
		xdiv02002f->DefineSection(1, -2700.0000*fgkmm,  478.03710000*fgkmm,  (478.03710000+50.0)*fgkmm);
		(new TGeoTranslation("xtran02002f", 0.0, 0.0, ITSSBIB02OffsetZ))->RegisterYourself(); // OffsetZ = to the tail
		
					
		//---------------------------- box add --------------------------								
		// box to add 02002a : paste the tail
		TGeoBBox *xadd02002a = new TGeoBBox("xadd02002a", (5.0*fgkmm)/2, 100*fgkmm, (350*fgkmm)/2);
		(new TGeoTranslation("xtranadd02002aR", 0.0 - xadd02002a->GetDX() + ((879.0 - 0.5*2)*fgkmm)/2, 0.0 - xadd02002a->GetDY(), ITSSBIB02OffsetZ + xadd02002a->GetDZ() - 2700*fgkmm))->RegisterYourself(); // OffsetZ = to the tail
		(new TGeoTranslation("xtranadd02002aL", 0.0 + xadd02002a->GetDX() - ((879.0 - 0.5*2)*fgkmm)/2, 0.0 - xadd02002a->GetDY(), ITSSBIB02OffsetZ + xadd02002a->GetDZ() - 2700*fgkmm))->RegisterYourself(); // OffsetZ = to the tail

		// box to add 02002b : roller handle
		TGeoBBox *xadd02002b = new TGeoBBox("xadd02002b", (9.0*fgkmm)/2, 53*fgkmm, ((2700.0-(574.7274+152.2762))*fgkmm)/2); // -(wing length in Z)
		(new TGeoTranslation("xtranadd02002bR", 0.0 - xadd02002b->GetDX() + ((917.0 + 9.0*2)*fgkmm)/2, 0.0 - xadd02002b->GetDY(), ITSSBIB02OffsetZ + xadd02002b->GetDZ() - 2700*fgkmm))->RegisterYourself(); // OffsetZ = to the tail
		(new TGeoTranslation("xtranadd02002bL", 0.0 + xadd02002b->GetDX() - ((917.0 + 9.0*2)*fgkmm)/2, 0.0 - xadd02002b->GetDY(), ITSSBIB02OffsetZ + xadd02002b->GetDZ() - 2700*fgkmm))->RegisterYourself(); // OffsetZ = to the tail
		
		// box to add 02002c : the wing
		TGeoBBox *xadd02002c = new TGeoBBox("xadd02002c", ((917.0 + 9.0*2)*fgkmm)/2, (5.0*fgkmm)/2, ((2700.000-(574.7274+152.2762))*fgkmm)/2); // -(wing length in Z)
		(new TGeoTranslation("xtranadd02002c", 0.0, 0.0 - xadd02002c->GetDY() - ITSSBIB02CutOffset5, ITSSBIB02OffsetZ + xadd02002c->GetDZ() - 2700*fgkmm))->RegisterYourself(); // OffsetZ = to the tail

		
		//#ifdef _ITSSB_SHOW_WITHBOXDIV

		{		
			TString sobj("(");
				sobj += "(";
					sobj += "(";
						sobj += "(xpconeITSSBIB02002:xtranITSSBIB02002 * xdiv02002e:xtran02002e)";
						sobj += " + ((xadd02002a:xtranadd02002aL + xadd02002a:xtranadd02002aR) - xdiv02002f:xtran02002f)";
						sobj += " + (xadd02002b:xtranadd02002bL + xadd02002b:xtranadd02002bR)";
					sobj += ")";
					sobj += " - (xdiv0200xa:xtran0200xa + xdiv0200xb:xtran0200xb + xdiv0200xc:xtran0200xc + xdiv02002d:xtran02002d)";
				sobj += ")";						
				sobj += " + (xadd02002c:xtranadd02002c - (xpconeITSSBIB02001s:xtranITSSBIB02001s * xdiv02001e:xtran02001e) )";
			sobj += ")";
			//sobj += " - (xdiv0200xt:xtran0200xt)";

		
			TGeoCompositeShape *xcomp = new TGeoCompositeShape("", sobj);
			TGeoVolume *xvol = new TGeoVolume("ITSSBIB02002", xcomp, MediumFoam);
			xvol->SetFillColor(5);
			xvol->SetLineColor(5);				
			svbVol->AddNode(xvol, 1, 0);

			// show wing
			//xvol = new TGeoVolume("ITSSBIB02002xxx", xadd02002c, MediumFoam);
			//xvol->SetFillColor(5);
			//xvol->SetLineColor(5);				
			//svbVol->AddNode(xvol, 1, new TGeoTranslation(0.0, 0.0 - xadd02002c->GetDY() - ITSSBIB02CutOffset5, ITSSBIB02OffsetZ + xadd02002c->GetDZ() - 2700*fgkmm));				

		}

		//#else
		//#endif 
		
		//-------------------------------------------------------------------------------------------------
		// 02003 - The outer																											
		
		TGeoPcon *xpcone02003 = new TGeoPcon("xpconeITSSBIB02003", 0, 360, 8);
		xpcone02003->DefineSection(0,    -0.0000*fgkmm,  49.70000000*fgkmm,  49.80000000*fgkmm);
		xpcone02003->DefineSection(1,   -52.3480*fgkmm,  49.70000000*fgkmm,  49.80000000*fgkmm);
		xpcone02003->DefineSection(2,   -54.8037*fgkmm,  49.70000000*fgkmm,  50.20813166*fgkmm);
		xpcone02003->DefineSection(3,  -374.3252*fgkmm, 102.59770330*fgkmm, 103.10449080*fgkmm);
		xpcone02003->DefineSection(4,  -374.5100*fgkmm, 102.62775960*fgkmm, 103.22392570*fgkmm);
		xpcone02003->DefineSection(5,  -574.7274*fgkmm, 232.61228230*fgkmm, 233.20891150*fgkmm);
		xpcone02003->DefineSection(6,  -574.9012*fgkmm, 232.72565530*fgkmm, 233.22906390*fgkmm);
		xpcone02003->DefineSection(7, -2700.0000*fgkmm, 478.03710000*fgkmm, 478.54040000*fgkmm);

		(new TGeoTranslation("xtranITSSBIB02003", ITSSBIB02OffsetX, ITSSBIB02OffsetY, ITSSBIB02OffsetZ))->RegisterYourself();

		//---------------------------- box div --------------------------

		// box to div 02003d
		TGeoBBox *xdiv02003d = new TGeoBBox("xdiv02003d", 500.0*fgkmm, 500.0*fgkmm, (2000.0*fgkmm)/2);
		(new TGeoTranslation("xtran02003d", 0.0, 0.0 + xdiv02003d->GetDY() - ITSSBIB02CutOffset4, ITSSBIB02OffsetZ + xdiv02003d->GetDZ() - (2700+1.0)*fgkmm))->RegisterYourself(); // OffsetZ = to the tail
		
		// box to div 02003e : crop the tail
		TGeoBBox *xdiv02003e = new TGeoBBox("xdiv02003e", (879.0*fgkmm)/2, 500.0*fgkmm, ((2700.0+1.0)*fgkmm)/2);
		(new TGeoTranslation("xtran02003e", 0.0, 0.0, ITSSBIB02OffsetZ - xdiv02003e->GetDZ()))->RegisterYourself(); // OffsetZ = to the tail
		
		// box to div 02003f : crop the paste of the tail
		TGeoPcon *xdiv02003f = new TGeoPcon("xdiv02003f", 0, 360, 2);
		xdiv02003f->DefineSection(0,  -574.7274*fgkmm,  233.20891150*fgkmm,  (233.20891150+50.0)*fgkmm);
		xdiv02003f->DefineSection(1, -2700.0000*fgkmm,  478.54040000*fgkmm,  (478.54040000+50.0)*fgkmm);
		(new TGeoTranslation("xtran02003f", 0.0, 0.0, ITSSBIB02OffsetZ))->RegisterYourself(); // OffsetZ = to the tail
		
					
		//---------------------------- box add --------------------------								
		// box to add 02003a : paste the tail
		TGeoBBox *xadd02003a = new TGeoBBox("xadd02003a", (0.5*fgkmm)/2, 200*fgkmm, (350*fgkmm)/2); 
		(new TGeoTranslation("xtranadd02003aR", 0.0 - xadd02003a->GetDX() + (879.0*fgkmm)/2, 0.0, ITSSBIB02OffsetZ + xadd02003a->GetDZ() - 2700*fgkmm))->RegisterYourself(); // OffsetZ = to the tail
		(new TGeoTranslation("xtranadd02003aL", 0.0 + xadd02003a->GetDX() - (879.0*fgkmm)/2, 0.0, ITSSBIB02OffsetZ + xadd02003a->GetDZ() - 2700*fgkmm))->RegisterYourself(); // OffsetZ = to the tail

		// box to add 02003b : roller handle
		TGeoBBox *xadd02003b = new TGeoBBox("xadd02003b", (0.5*fgkmm)/2, 53*2*fgkmm, ((2700.0-(574.7274+152.2762))*fgkmm)/2); // -(wing length in Z)
		(new TGeoTranslation("xtranadd02003bR", 0.0 - xadd02003b->GetDX() + (917.0*fgkmm)/2, 0.0, ITSSBIB02OffsetZ + xadd02003b->GetDZ() - 2700*fgkmm))->RegisterYourself(); // OffsetZ = to the tail
		(new TGeoTranslation("xtranadd02003bL", 0.0 + xadd02003b->GetDX() - (917.0*fgkmm)/2, 0.0, ITSSBIB02OffsetZ + xadd02003b->GetDZ() - 2700*fgkmm))->RegisterYourself(); // OffsetZ = to the tail
		
		// box to add 02003c : the wing
		TGeoBBox *xadd02003c = new TGeoBBox("xadd02003c", (917.0*fgkmm)/2, (0.5*fgkmm)/2, ((2700.000-(574.7274+152.2762))*fgkmm)/2); // -(wing length in Z)
		(new TGeoTranslation("xtranadd02003c", 0.0, 0.0 - xadd02003c->GetDY() - ITSSBIB02CutOffset4, ITSSBIB02OffsetZ + xadd02003c->GetDZ() - 2700*fgkmm))->RegisterYourself(); // OffsetZ = to the tail

		
		//#ifdef _ITSSB_SHOW_WITHBOXDIV
		{
			/*TString sobj("(");
				sobj += "(";
					sobj += "(";
						sobj += "(xpconeITSSBIB02003:xtranITSSBIB02003 * xdiv02003e:xtran02003e)";
						sobj += " + ((xadd02003a:xtranadd02003aL + xadd02003a:xtranadd02003aR) - xdiv02003f:xtran02003f)";
						sobj += " + (xadd02003b:xtranadd02003bL + xadd02003b:xtranadd02003bR)";
					sobj += ")";
					sobj += " - (xdiv0200xa:xtran0200xa + xdiv0200xb:xtran0200xb + xdiv0200xc:xtran0200xc + xdiv02003d:xtran02003d)";
				sobj += ")";						
				sobj += " + (xadd02003c:xtranadd02003c - (xpconeITSSBIB02002s:xtranITSSBIB02002s * xdiv02002e:xtran02002e) )";
			sobj += ")";
			sobj += " - (xdiv0200xt:xtran0200xt)";
			*/
			
			TString sobj("(");
				sobj += "(";
					sobj += "(";
						sobj += "(xpconeITSSBIB02003:xtranITSSBIB02003 * xdiv02003e:xtran02003e)";
						sobj += " + ((xadd02003a:xtranadd02003aL + xadd02003a:xtranadd02003aR) - xdiv02003f:xtran02003f)";
						sobj += " + (xadd02003b:xtranadd02003bL + xadd02003b:xtranadd02003bR)";
					sobj += ")";
					sobj += " - (xdiv0200xa:xtran0200xa + xdiv0200xb:xtran0200xb + xdiv0200xc:xtran0200xc + xdiv02003d:xtran02003d)";
				sobj += ")";						
				sobj += " + (xadd02003c:xtranadd02003c - (xpconeITSSBIB02002s:xtranITSSBIB02002s * xdiv02002e:xtran02002e) )";
			sobj += ")";
			//sobj += " - (xdiv0200xt:xtran0200xt)";
									
			TGeoCompositeShape *xcomp = new TGeoCompositeShape("", sobj);
			TGeoVolume *xvol = new TGeoVolume("ITSSBIB020031", xcomp, MediumWall);
			xvol->SetFillColor(6);
			xvol->SetLineColor(6);
			svbVol->AddNode(xvol, 1, 0);

/*						TString sobj2("(");
				sobj2 += "(";
					sobj2 += "(";
						//sobj2 += "(xpconeITSSBIB02003:xtranITSSBIB02003 * xdiv02003e:xtran02003e)";
						//sobj2 += " + ((xadd02003a:xtranadd02003aL + xadd02003a:xtranadd02003aR) - xdiv02003f:xtran02003f)";
						sobj2 += "(xadd02003b:xtranadd02003bL + xadd02003b:xtranadd02003bR)";
					sobj2 += ")";
					sobj2 += " - (xdiv0200xa:xtran0200xa + xdiv0200xb:xtran0200xb + xdiv0200xc:xtran0200xc + xdiv02003d:xtran02003d)";
				sobj2 += ")";						
				sobj2 += " + (xadd02003c:xtranadd02003c - (xpconeITSSBIB02002s:xtranITSSBIB02002s * xdiv02002e:xtran02002e) )";
			sobj2 += ")";
			sobj2 += " - (xdiv0200xt:xtran0200xt)";
			

			TGeoCompositeShape *xcomp2 = new TGeoCompositeShape("", sobj2);
			TGeoVolume *xvol2 = new TGeoVolume("ITSSBIB020032", xcomp2, MediumWall);
			xvol2->SetFillColor(5);
			xvol2->SetLineColor(5);				
			svbVol->AddNode(xvol2, 1, 0);
*/
			// show cone
			//xvol = new TGeoVolume("ITSSBIB02003xx1", xpcone02003, MediumWall);
			//xvol->SetFillColor(6);
			//xvol->SetLineColor(6);				
			//svbVol->AddNode(xvol, 1, new TGeoTranslation("", ITSSBIB02OffsetX, ITSSBIB02OffsetY, ITSSBIB02OffsetZ));				
			
			// show box
			//xvol = new TGeoVolume("ITSSBIB02003xx2", xdiv02003e, MediumWall);
			//xvol->SetFillColor(5);
			//xvol->SetLineColor(5);				
			//svbVol->AddNode(xvol, 1, new TGeoTranslation("", 0.0, 0.0, ITSSBIB02OffsetZ - xdiv02003e->GetDZ()));				

			// show wing
			//xvol = new TGeoVolume("ITSSBIB02003xxx", xadd02003c, Medium);
			//xvol->SetFillColor(5);
			//xvol->SetLineColor(5);				
			//svbVol->AddNode(xvol, 1, new TGeoTranslation(0.0, 0.0 - xadd02003c->GetDY() - ITSSBIB02CutOffset4, ITSSBIB02OffsetZ + xadd02003c->GetDZ() - 2700*fgkmm));				
			
			// show s
			//TGeoCompositeShape *xcomp2 = new TGeoCompositeShape("", "(xpconeITSSBIB02002s:xtranITSSBIB02002s * xdiv02002e:xtran02002e)");
			//xvol = new TGeoVolume("ITSSBIB02003xxx", xcomp2, MediumWall);
			//xvol->SetFillColor(5);
			//xvol->SetLineColor(5);				
			//svbVol->AddNode(xvol, 1, new TGeoTranslation(ITSSBIB02OffsetX, ITSSBIB02OffsetY, ITSSBIB02OffsetZ));				

		//#else
		//#endif 
		}
		

	}
	
	// 0201x
	if (true)
	{ 	
		//---------------------------- box div --------------------------
		//------------------ GLOBAL VALUE -------------------------------

		// box to div 0201xd
		TGeoBBox *xdiv0201xd = new TGeoBBox("xdiv0201xd", 500.0*fgkmm, 500.0*fgkmm, (2500.0*fgkmm)/2);
		(new TGeoTranslation("xtran0201xd", 0.0, 0.0 + xdiv0201xd->GetDY() - ITSSBIB02CutOffset3, ITSSBIB02OffsetZ + xdiv0201xd->GetDZ() - (2700+1.0)*fgkmm))->RegisterYourself(); // OffsetZ = to the tail
		
		// box to div 0201xe : crop the tail
		TGeoBBox *xdiv0201xe = new TGeoBBox("xdiv0201xe", ((879.0 - 0.5*2 - 5.0*2)*fgkmm)/2, 500.0*fgkmm, ((2700.0+1.0)*fgkmm)/2);
		(new TGeoTranslation("xtran0201xe", 0.0, 0.0, ITSSBIB02OffsetZ - xdiv0201xe->GetDZ()))->RegisterYourself(); // OffsetZ = to the tail
	

		//-------------------------------------------------------------------------------------------------
		// 02011					
		
		TGeoPcon *xpcone02011 = new TGeoPcon("xpconeITSSBIB02011", 0, 360, 3);
		xpcone02011->DefineSection(0,  -542.1782*fgkmm, 204.92400000*fgkmm, 204.92400000*fgkmm);
		xpcone02011->DefineSection(1,  -543.1212*fgkmm, 205.03300000*fgkmm, 205.53600000*fgkmm);
		xpcone02011->DefineSection(2, -2700.0000*fgkmm, 454.01280000*fgkmm, 454.51610000*fgkmm);
		(new TGeoTranslation("xtranITSSBIB02011", ITSSBIB02OffsetX, ITSSBIB02OffsetY, ITSSBIB02OffsetZ))->RegisterYourself();


									
		//---------------------------- box add --------------------------								

		
		//#ifdef _ITSSB_SHOW_WITHBOXDIV
		{
			TString sobj("");
			sobj += "((xpconeITSSBIB02011:xtranITSSBIB02011)";
			sobj += " - (xdiv0201xd:xtran0201xd)";
			//sobj += " - ((xdiv02011a:xtran02011a)+(xdiv02011d:xtran02011d))";
			sobj += " * (xdiv0201xe:xtran0201xe))";
		
			TGeoCompositeShape *xcomp = new TGeoCompositeShape("", sobj);
			TGeoVolume *xvol = new TGeoVolume("ITSSBIB02011", xcomp, MediumWall);
			xvol->SetFillColor(6);
			xvol->SetLineColor(6);				
			svbVol->AddNode(xvol, 1, 0);

			// show wing
			//xvol = new TGeoVolume("ITSSBIB02011xxx", xadd02011c, MediumWall);
			//xvol->SetFillColor(6);
			//xvol->SetLineColor(6);				
			//svbVol->AddNode(xvol, 1, new TGeoTranslation(0.0, 0.0 - xadd02011c->GetDY() - ITSSBIB02CutOffset6, ITSSBIB02OffsetZ + xadd02011c->GetDZ() - 2700*fgkmm));				

		}

		//#else
		//#endif 

		//-------------------------------------------------------------------------------------------------
		// 02012					
		
		TGeoPcon *xpcone02012 = new TGeoPcon("xpconeITSSBIB02012", 0, 360, 3);
		xpcone02012->DefineSection(0,  -532.7489*fgkmm, 198.80200000*fgkmm, 198.80200000*fgkmm);
		xpcone02012->DefineSection(1,  -542.1782*fgkmm, 199.89100000*fgkmm, 205.03300000*fgkmm);
		xpcone02012->DefineSection(2, -2700.0000*fgkmm, 448.97960000*fgkmm, 454.01280000*fgkmm);
		(new TGeoTranslation("xtranITSSBIB02012", ITSSBIB02OffsetX, ITSSBIB02OffsetY, ITSSBIB02OffsetZ))->RegisterYourself();


									
		//---------------------------- box add --------------------------								

		
		//#ifdef _ITSSB_SHOW_WITHBOXDIV
		{
			TString sobj("");
			sobj += "((xpconeITSSBIB02012:xtranITSSBIB02012)";
			sobj += " - (xdiv0201xd:xtran0201xd)";
			//sobj += " - ((xdiv02012a:xtran02012a)+(xdiv02012d:xtran02012d))";
			sobj += " * (xdiv0201xe:xtran0201xe))";
		
			TGeoCompositeShape *xcomp = new TGeoCompositeShape("", sobj);
			TGeoVolume *xvol = new TGeoVolume("ITSSBIB02012", xcomp, MediumFoam);
			xvol->SetFillColor(5);
			xvol->SetLineColor(5);				
			svbVol->AddNode(xvol, 1, 0);

			// show wing
			//xvol = new TGeoVolume("ITSSBIB02012xxx", xadd02012c, MediumWall);
			//xvol->SetFillColor(6);
			//xvol->SetLineColor(6);				
			//svbVol->AddNode(xvol, 1, new TGeoTranslation(0.0, 0.0 - xadd02012c->GetDY() - ITSSBIB02CutOffset6, ITSSBIB02OffsetZ + xadd02012c->GetDZ() - 2700*fgkmm));				

		}

		//#else
		//#endif 

		//-------------------------------------------------------------------------------------------------
		// 02013					
		
		TGeoPcon *xpcone02013 = new TGeoPcon("xpconeITSSBIB02013", 0, 360, 3);
		xpcone02013->DefineSection(0,  -531.8061*fgkmm, 198.19000000*fgkmm, 198.19000000*fgkmm);
		xpcone02013->DefineSection(1,  -532.7489*fgkmm, 198.29900000*fgkmm, 198.80200000*fgkmm);
		xpcone02013->DefineSection(2, -2700.0000*fgkmm, 448.47630000*fgkmm, 448.97960000*fgkmm);
		(new TGeoTranslation("xtranITSSBIB02013", ITSSBIB02OffsetX, ITSSBIB02OffsetY, ITSSBIB02OffsetZ))->RegisterYourself();


									
		//---------------------------- box add --------------------------								

		
		//#ifdef _ITSSB_SHOW_WITHBOXDIV
		{
			TString sobj("");
			sobj += "((xpconeITSSBIB02013:xtranITSSBIB02013)";
			sobj += " - (xdiv0201xd:xtran0201xd)";
			//sobj += " - ((xdiv02013a:xtran02013a)+(xdiv02013d:xtran02013d))";
			sobj += " * (xdiv0201xe:xtran0201xe))";
		
			TGeoCompositeShape *xcomp = new TGeoCompositeShape("", sobj);
			TGeoVolume *xvol = new TGeoVolume("ITSSBIB02013", xcomp, MediumWall);
			xvol->SetFillColor(6);
			xvol->SetLineColor(6);				
			svbVol->AddNode(xvol, 1, 0);

			// show wing
			//xvol = new TGeoVolume("ITSSBIB02013xxx", xadd02013c, MediumWall);
			//xvol->SetFillColor(6);
			//xvol->SetLineColor(6);				
			//svbVol->AddNode(xvol, 1, new TGeoTranslation(0.0, 0.0 - xadd02013c->GetDY() - ITSSBIB02CutOffset6, ITSSBIB02OffsetZ + xadd02013c->GetDZ() - 2700*fgkmm));				

		}

		//#else
		//#endif 


	}
	

	return svbVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateInnerBServiceB2(const TGeoManager *mgr){

	const Double_t kOffsetZ = -32.1*fgkmm;

	TGeoVolume* ewhVol = new TGeoVolumeAssembly("InnerBServiceB");
	
	TGeoMedium *MediumWall = mgr->GetMedium("ITS_PEEKCF30$");	
	TGeoMedium *MediumFoam = mgr->GetMedium("ITS_CarbonFoam$");	
		
	TGeoShape* shp;
	TGeoVolume* vol;
	
	// -- Wall Inside
	shp = new TGeoPcon("", 0, 180, 6);
	((TGeoPcon*)shp)->DefineSection(0,    -2.00*fgkmm,  47.30*fgkmm,  47.80*fgkmm);
	((TGeoPcon*)shp)->DefineSection(1,   -11.13*fgkmm,  47.30*fgkmm,  47.80*fgkmm);
	((TGeoPcon*)shp)->DefineSection(2,   -11.34*fgkmm,  47.30*fgkmm,  48.00*fgkmm);
	((TGeoPcon*)shp)->DefineSection(3,   -86.78*fgkmm, 123.63*fgkmm, 124.34*fgkmm);
	((TGeoPcon*)shp)->DefineSection(4,   -87.00*fgkmm, 123.00*fgkmm, 124.34*fgkmm);
	((TGeoPcon*)shp)->DefineSection(5,  -371.35*fgkmm, 123.00*fgkmm, 124.34*fgkmm);
	vol = new TGeoVolume("WallInside", shp, MediumWall);
	vol->SetFillColor(56);
	vol->SetLineColor(56);
	ewhVol->AddNode(vol, 1, new TGeoCombiTrans(0, 0, +kOffsetZ, new TGeoRotation("", 0, 0, 0)));


	// -- Wall Middle
	shp = new TGeoPcon("", 0, 180, 8);
	((TGeoPcon*)shp)->DefineSection(0,    -2.00*fgkmm,  47.80*fgkmm,  49.00*fgkmm);
	((TGeoPcon*)shp)->DefineSection(1,    -6.50*fgkmm,  47.80*fgkmm,  49.00*fgkmm);
	((TGeoPcon*)shp)->DefineSection(2,    -6.50*fgkmm,  47.80*fgkmm,  52.80*fgkmm);
	((TGeoPcon*)shp)->DefineSection(3,    -9.06*fgkmm,  47.80*fgkmm,  52.80*fgkmm);
	((TGeoPcon*)shp)->DefineSection(4,   -11.13*fgkmm,  47.80*fgkmm,  54.00*fgkmm); // rmax fix ovl 55.0 -> 54.0
	((TGeoPcon*)shp)->DefineSection(5,   -84.71*fgkmm, 124.34*fgkmm, 128.50*fgkmm);
	((TGeoPcon*)shp)->DefineSection(6,   -86.78*fgkmm, 124.34*fgkmm, 128.50*fgkmm);
	((TGeoPcon*)shp)->DefineSection(7,  -371.35*fgkmm, 124.34*fgkmm, 128.50*fgkmm);
	vol = new TGeoVolume("WallFoam", shp, MediumFoam);
	vol->SetFillColor(5);
	vol->SetLineColor(5);
	ewhVol->AddNode(vol, 1, new TGeoCombiTrans(0, 0, +kOffsetZ, new TGeoRotation("", 0, 0, 0)));


	// -- Wall Outside
	shp = new TGeoPcon("", 0, 180, 6);
	((TGeoPcon*)shp)->DefineSection(0,    -6.50*fgkmm,  52.80*fgkmm,  53.30*fgkmm);
	((TGeoPcon*)shp)->DefineSection(1,    -8.85*fgkmm,  52.80*fgkmm,  53.30*fgkmm);
	((TGeoPcon*)shp)->DefineSection(2,   -09.06*fgkmm,  52.80*fgkmm,  53.52*fgkmm);
	((TGeoPcon*)shp)->DefineSection(3,   -84.50*fgkmm, 128.50*fgkmm, 129.15*fgkmm); // rmin fix ovl 128.44 -> 128.50
	((TGeoPcon*)shp)->DefineSection(4,   -84.71*fgkmm, 128.50*fgkmm, 129.15*fgkmm);
	((TGeoPcon*)shp)->DefineSection(5,  -371.35*fgkmm, 128.50*fgkmm, 129.15*fgkmm);
	vol = new TGeoVolume("WallOutside", shp, MediumWall);
	vol->SetFillColor(56);
	vol->SetLineColor(56);
	ewhVol->AddNode(vol, 1, new TGeoCombiTrans(0, 0, +kOffsetZ, new TGeoRotation("", 0, 0, 0)));



	return ewhVol;
}


TGeoVolume* AliITSUv2Layer::CreateOuterBServiceB(const TGeoManager *mgr){


	TGeoVolume* svbVol = new TGeoVolumeAssembly("OuterBServiceB");

	return svbVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateDCSet(const char *nid, const TGeoManager *mgr){

	TGeoVolume* setVol = new TGeoVolumeAssembly(TString("DCSet")+=nid);
	
	TGeoVolume* unitVol;
	TGeoVolume* baseVol;
	TGeoVolume* rndpVol;
	
	baseVol = CreateDCBase();
	setVol->AddNode(baseVol, 1, 0);

	unitVol = CreateDCUnit("1");
	setVol->AddNode(unitVol, 1, new TGeoTranslation(
		0, 
		+ ((TGeoBBox*)(baseVol->GetShape()))->GetDY() + 00.74*fgkmm/2, 
		- ((TGeoBBox*)(baseVol->GetShape()))->GetDZ() + ((TGeoBBox*)(unitVol->GetShape()))->GetDZ() + 13.3*fgkmm/2 ));

	unitVol = CreateDCUnit("2");
	setVol->AddNode(unitVol, 1, new TGeoCombiTrans(
		0, 
		- ((TGeoBBox*)(baseVol->GetShape()))->GetDY() - 00.74*fgkmm/2, 
		- ((TGeoBBox*)(baseVol->GetShape()))->GetDZ() + ((TGeoBBox*)(unitVol->GetShape()))->GetDZ() + 13.3*fgkmm/2 + 29.6*fgkmm
		, new TGeoRotation("", 0, 0, 180) ));
	
	rndpVol = CreateRoundP();
	setVol->AddNode(rndpVol, 1, new TGeoTranslation(
		0, 
		+ ((TGeoBBox*)(baseVol->GetShape()))->GetDY() 
			+ ((TGeoBBox*)(unitVol->GetShape()))->GetDY()*2
			+ 1.5*fgkmm, 
		0));
	
	return setVol;
}



//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateDCUnit(const char *nid, const TGeoManager *mgr){

	const Double_t kPlateDX = 13.30*2*fgkmm;
	const Double_t kPlateDY = 00.74*fgkmm;
	const Double_t kPlateDZ = 68.00*fgkmm;
	
	const Double_t kBoxDX = 10.35*2*fgkmm;
	const Double_t kBoxDY = 10.35*fgkmm;
	const Double_t kBoxDZ = 29.60*fgkmm;

	const Double_t kTabDX = 10.35*2*fgkmm;
	const Double_t kTabDY = 03.00*fgkmm;
	const Double_t kTabDZ = 04.44*fgkmm;
	
	const Double_t kHoleR = 02.00*fgkmm;

	TGeoVolume* unitVol = new TGeoVolumeAssembly(TString("DCUnit")+=nid);
	TGeoVolume* volume;
	
	// -- Plate --
	volume = new TGeoVolume("DCPlate", new TGeoBBox(kPlateDX/2, kPlateDY/2, kPlateDZ/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	volume->SetFillColor(5);
	volume->SetLineColor(5);
	unitVol->AddNode(volume, 0, 0);
	
	// -- Box --
	volume = new TGeoVolume("DCBox", new TGeoBBox(kBoxDX/2, kBoxDY/2, kBoxDZ/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	volume->SetFillColor(6);
	volume->SetLineColor(6);
	unitVol->AddNode(volume, 0, new TGeoTranslation(0, + kBoxDY/2 + kPlateDY/2, - kPlateDZ/2 + kBoxDZ/2 + 28.10*fgkmm));
	
	// -- Tab --
	volume = new TGeoVolume("DCTab", new TGeoBBox(kTabDX/2, kTabDY/2, kTabDZ/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	volume->SetFillColor(7);
	volume->SetLineColor(7);
	unitVol->AddNode(volume, 0, new TGeoTranslation(0, + kTabDY/2 + kPlateDY/2, - kPlateDZ/2 + kTabDZ/2 + 10.35*fgkmm));
	
	
	return unitVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateDCBase(const char *nid, const TGeoManager *mgr){
	
	const Double_t kColdPipeR = 02.00*fgkmm;
	const Double_t kColdPipeDZ = 121.30*fgkmm;
	
	const Double_t kWidthBtwPipes = 13.30*2*fgkmm; // kPlateDX
	
	const Double_t kColdPlateDX = kWidthBtwPipes - kColdPipeR*2;
	const Double_t kColdPlateDY = kColdPipeR*2;
	const Double_t kColdPlateDZ = 108.00*fgkmm;
	
	TGeoVolume* baseVol = new TGeoVolumeAssembly(TString("DCBase")+=nid);
	TGeoVolume* volume;

	// -- Pipes --
	volume = new TGeoVolume("DCColdPipe1", new TGeoTube(kColdPipeR - 0.5*fgkmm, kColdPipeR, kColdPipeDZ/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	volume->SetFillColor(8);
	volume->SetLineColor(8);
	baseVol->AddNode(volume, 0, new TGeoTranslation(- kWidthBtwPipes/2, 0, 0));

	volume = new TGeoVolume("DCColdPipe2", new TGeoTube(kColdPipeR - 0.5*fgkmm, kColdPipeR, kColdPipeDZ/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	volume->SetFillColor(8);
	volume->SetLineColor(8);
	baseVol->AddNode(volume, 0, new TGeoTranslation(+ kWidthBtwPipes/2, 0, 0));
	
	// -- Plate --
	volume = new TGeoVolume("DCColdPlate", new TGeoBBox(kColdPlateDX/2, kColdPlateDY/2, kColdPlateDZ/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	volume->SetFillColor(9);
	volume->SetLineColor(9);
	baseVol->AddNode(volume, 0, 0);
	
	
	return baseVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateRoundP(const char *nid, const TGeoManager *mgr){
	
	const Double_t kRPipesR = 1.5*fgkmm;
	const Double_t kRPipesDZ = 121.30*fgkmm;
	
	TGeoVolume* baseVol = new TGeoVolumeAssembly(TString("RoundPipes")+=nid);
	TGeoVolume* volume;

	for (int i=-2; i<=2; i++) // first stack
	{
		volume = new TGeoVolume("DCRoundPipe", new TGeoTube(0, kRPipesR, kRPipesDZ/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
		volume->SetFillColor(3);
		volume->SetLineColor(3);
		baseVol->AddNode(volume, 0, new TGeoTranslation(i*kRPipesR*2, 0, 0));
	}

	for (int i=-5; i<=5; i+=2) // second stack
	{
		volume = new TGeoVolume("DCRoundPipe", new TGeoTube(0, kRPipesR, kRPipesDZ/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
		volume->SetFillColor(3);
		volume->SetLineColor(3);
		baseVol->AddNode(volume, 0, new TGeoTranslation(i*kRPipesR, kRPipesR*2, 0));
	}


	return baseVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateDCPSLayer(const char *nid, const TGeoManager *mgr){

	TGeoVolume* setVol = new TGeoVolumeAssembly(TString("DCPSSet")+TString(nid)+(TString("L")+=fLayerNumber));

	const Double_t kTwistDegree = 15.0;
	const Double_t kDCOffsetRadius[7] = {36.5*fgkmm, 52.3*fgkmm, 66.51*fgkmm, 0, 0, 0, 0};
	const Double_t kSlopeDegree[7] = {5.2, 17.2, 29.2, 0, 0, 0, 0};

	if (fLayerNumber < fgkNumberOfInnerLayers)
	{		
		const Double_t fkDegPerStave = 360/fNStaves;
		const Double_t fkDeg2RadFactor = TMath::Pi()/180;

		for (int i=0; i<fNStaves/2; i++)
		{
			Double_t fCosDegIter = TMath::Cos(i*fkDegPerStave*fkDeg2RadFactor + (fkDegPerStave*fkDeg2RadFactor)/2);
			Double_t fSinDegIter = TMath::Sin(i*fkDegPerStave*fkDeg2RadFactor + (fkDegPerStave*fkDeg2RadFactor)/2);

			setVol->AddNode(CreateDCPS2((TString("L")+=fLayerNumber)+=(TString("_")+=i)), 1, 
				new TGeoCombiTrans((kDCOffsetRadius[(int)fLayerNumber] * fSinDegIter), (kDCOffsetRadius[(int)fLayerNumber] * fCosDegIter), 0, 
				new TGeoRotation("", -i*fkDegPerStave -fkDegPerStave/2 -kTwistDegree +180.0, -kSlopeDegree[(int)fLayerNumber], 0)));
		}
		
	}
	
	return setVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateDCPS1(const char *nid, const TGeoManager *mgr){

	TGeoVolume* setVol = new TGeoVolumeAssembly(TString("DCPS")+=nid);
	TGeoVolume* volume;

	// -- Pipes --	
	const Double_t kPipeHoleRadius = 0.6*fgkmm;
	const Double_t kPipeThick = 0.2*fgkmm;
	const Double_t kPipeLength = 77.0*fgkmm;
	const Double_t kWidthBtwnPipes = 12.4*fgkmm;
	volume = new TGeoVolume("CoolingPipe", new TGeoTube(kPipeHoleRadius, kPipeHoleRadius+kPipeThick, kPipeLength/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	volume->SetFillColor(8);
	volume->SetLineColor(8);
	setVol->AddNode(volume, 0, new TGeoTranslation(- kWidthBtwnPipes/2, 0, 0));
	setVol->AddNode(volume, 0, new TGeoTranslation(+ kWidthBtwnPipes/2, 0, 0));

	// -- Plate --
	const Double_t kPlateLength = 54.0*fgkmm;
	const Double_t kPlateThick = 2*(kPipeHoleRadius + kPipeThick);
	const Double_t kPlateWidth = kWidthBtwnPipes - 2*(kPipeHoleRadius + kPipeThick); // instead of using the composite object
	const Double_t kPlateOffset = 2.5*fgkmm;
	volume = new TGeoVolume("CoolingPlate", new TGeoBBox(kPlateWidth/2, kPlateThick/2, kPlateLength/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	volume->SetFillColor(9);
	volume->SetLineColor(9);
	setVol->AddNode(volume, 0, new TGeoTranslation(0, 0, -kPlateOffset));

	// -- Up and Down Plates
	const Double_t kWidthBtwnUDPlates = 1.7*fgkmm;

	// -- Up Plate --
	const Double_t kUpPlateLength = 76.8*fgkmm;
	const Double_t kUpPlateThick = 0.1*fgkmm;
	const Double_t kUpPlateWidth = 16.0*fgkmm;
	const Double_t kUpPlateOffset = 4.5*fgkmm;
	volume = new TGeoVolume("UpPlate", new TGeoBBox(kUpPlateWidth/2, kUpPlateThick/2, kUpPlateLength/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	volume->SetFillColor(9);
	volume->SetLineColor(9);
	setVol->AddNode(volume, 0, new TGeoTranslation(0, +kWidthBtwnUDPlates/2+kUpPlateThick/2, -(kPipeLength/2-kUpPlateLength/2)-kUpPlateOffset));
	
	// -- Up Extend Plate 0 --
	const Double_t kUpEPlateLength = 32.3*fgkmm;
	const Double_t kUpEPlateThick = 0.1*fgkmm;
	const Double_t kUpEPlateWidth = 16.0*fgkmm;
	volume = new TGeoVolume("UpEPlate", new TGeoBBox(kUpEPlateWidth/2, kUpEPlateThick/2, kUpEPlateLength/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	volume->SetFillColor(8);
	volume->SetLineColor(8);
	TGeoVolume* asmVol = new TGeoVolumeAssembly("");
	asmVol->AddNode(volume, 0, new TGeoTranslation(0, -kUpEPlateThick/2, +kUpEPlateLength/2));
	setVol->AddNode(asmVol, 0, new TGeoCombiTrans(0, +kWidthBtwnUDPlates/2+kUpPlateThick/2 +kUpEPlateThick/2, -kPipeLength/2 -kUpPlateOffset +kUpPlateLength, new TGeoRotation("", 0, -13, 0)));
	
	// -- Dn Plate --
	const Double_t kDnPlateLength = 69.2*fgkmm;
	const Double_t kDnPlateThick = 0.1*fgkmm;
	const Double_t kDnPlateWidth = 16.0*fgkmm;
	const Double_t kDnPlateOffset = 4.5*fgkmm;
	volume = new TGeoVolume("DnPlate", new TGeoBBox(kDnPlateWidth/2, kDnPlateThick/2, kDnPlateLength/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	volume->SetFillColor(9);
	volume->SetLineColor(9);
	setVol->AddNode(volume, 0, new TGeoTranslation(0, -kWidthBtwnUDPlates/2-kDnPlateThick/2, -(kPipeLength/2-kDnPlateLength/2)-kDnPlateOffset));
	
	// -- Power Supply Box
	const Double_t kUpBoxOffset = 27.3*fgkmm;
	const Double_t kDnBoxOffset = 09.1*fgkmm;
	setVol->AddNode(CreateDCBox1("Up"), 1, new TGeoTranslation(0, +kWidthBtwnUDPlates/2 +kUpPlateThick, -kPipeLength/2 +kUpBoxOffset));
	setVol->AddNode(CreateDCBox1("Dn"), 1, new TGeoCombiTrans(0, -kWidthBtwnUDPlates/2 -kDnPlateThick, -kPipeLength/2 +kDnBoxOffset,new TGeoRotation("", 0, 0, 180) ));
		
	// -- Pipes
	const Double_t kPipesOffsetZ = 0*fgkmm;
	const Double_t kPipesOffsetY = 8.5*fgkmm +kWidthBtwnUDPlates/2 +kUpPlateThick; // 8.5 = DCBox = 8.0 + 0.4 + empty space
	setVol->AddNode(CreateDCPipes(), 1, new TGeoTranslation(0, +kPipesOffsetY, +kPipesOffsetZ));
	
	// -- FPGA Connector
	const Double_t kCntOffsetZ = +kPipeLength/2 -7.0*fgkmm;
	const Double_t kCntOffsetY = 8.4*fgkmm +kWidthBtwnUDPlates/2 +kUpPlateThick -0.5*fgkmm; // 8.4 =DCBox, 0.5 = obj thick
	const Double_t kCntOffsetZ2 = 60.0*fgkmm;
	setVol->AddNode(CreateFPGAConnector("1"), 1, new TGeoTranslation(0, +kCntOffsetY, +kCntOffsetZ));
	setVol->AddNode(CreateFPGAConnector("2"), 1, new TGeoCombiTrans(0, +kCntOffsetY, -kCntOffsetZ2 +2*fgkmm, new TGeoRotation("",0,12,0)));
	setVol->AddNode(CreateFPGAConnector("3"), 1, new TGeoCombiTrans(0, +kCntOffsetY/2, -kCntOffsetZ2, new TGeoRotation("",0,12,0)));
	
	TGeoVolume* newOffsetVol = new TGeoVolumeAssembly(TString("DCPS")+=nid);
	newOffsetVol->AddNode(setVol, 1, new TGeoTranslation(0, 0, -kPipeLength/2));
	
	
	return newOffsetVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateDCBox1(const char *nid, const TGeoManager *mgr){
	
	TGeoVolume* unitVol = new TGeoVolumeAssembly(TString("DCBox")+=nid);
	TGeoVolume* volume;
	
	// -- Base --
	const Double_t kBaseDX = 16.9*fgkmm/2;
	const Double_t kBaseDY = 00.4*fgkmm/2;
	const Double_t kBaseDZ = 43.5*fgkmm/2; // add all translation +Z for change the center point to the edge
	volume = new TGeoVolume("BoxBase", new TGeoBBox(kBaseDX, kBaseDY, kBaseDZ), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	volume->SetFillColor(5);
	volume->SetLineColor(5);
	unitVol->AddNode(volume, 0, new TGeoTranslation(0, kBaseDY, +kBaseDZ));
	
	// -- Box --
	const Double_t kBoxDX = 14.0*fgkmm/2;
	const Double_t kBoxDY = 08.0*fgkmm/2;
	const Double_t kBoxDZ = 18.5*fgkmm/2;
	const Double_t kBoxOffset = 6.7*fgkmm; // offset from the A-Side of the edge of the base
	volume = new TGeoVolume("BoxBox", new TGeoBBox(kBoxDX, kBoxDY, kBoxDZ), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	volume->SetFillColor(6);
	volume->SetLineColor(6);
	unitVol->AddNode(volume, 0, new TGeoTranslation(0, +kBoxDY +kBaseDY, -kBaseDZ +kBoxDZ +kBoxOffset +kBaseDZ));
	
	// -- Bar --
	const Double_t kBarDX = 15.0*fgkmm/2;
	const Double_t kBarDY = 02.3*fgkmm/2;
	const Double_t kBarDZ = 03.0*fgkmm/2;	
	const Double_t kBarOffset = 6.5*fgkmm; // offset from the C-Side of the edge of the base
	volume = new TGeoVolume("BoxBar", new TGeoBBox(kBarDX, kBarDY, kBarDZ), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	volume->SetFillColor(7);
	volume->SetLineColor(7);
	unitVol->AddNode(volume, 0, new TGeoTranslation(0, +kBarDY +kBaseDY, +kBaseDZ -kBarDZ -kBarOffset +kBaseDZ));
	
	// -- Holes --
	//const Double_t kHoleRadius = 2.0*fgkmm;
	//const Double_t kHoleAOffset = 3.6*fgkmm;
	//const Double_t kHoleCOffset = 2.2*fgkmm;
	//const Double_t kWidthBtwnAHoles = 9.4*fgkmm;
	//const Double_t kHoleLength = 2.0*fgkmm;
	//volume = new TGeoVolume("BoxBaseHole", new TGeoTube(0, kHoleRadius, kHoleLength/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	//volume->SetFillColor(8);
	//volume->SetLineColor(8);
	//setVol->AddNode(volume, 0, new TGeoCombiTrans(0, 0, 0, new TGeoRotations("", 0, 0, 0))); // not finished


	
	return unitVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateInnerDCCNT3(const TGeoManager *mgr){
		
	TString strName = (TString("DCCNT")+=fLayerNumber) + (TString("_")+=(int)fN_DCCNT_Created);
	TString strObjName, strObjTran, strObjTran2;
	
	TGeoCombiTrans *cbtObj;	TGeoVolume* setVol = new TGeoVolumeAssembly(strName);
	TGeoVolume* vol;

	TString strVol = TString("(");

	// -- Cooling Plate --
	const Double_t kCPlateWidth = 15.0*fgkmm;
	const Double_t kCPlateThick = 3.0*fgkmm;
	const Double_t kCPlateLength = 21.22*fgkmm;

	vol = new TGeoVolume(strName+TString("_CPlate0Cnt"), new TGeoBBox(kCPlateWidth/2, kCPlateThick/2, kCPlateLength/2), gGeoManager->GetMedium("ITS_ALUMINUM$"));
	vol->SetFillColor(5);	vol->SetLineColor(5);
	//setVol->AddNode(vol, 0, new TGeoCombiTrans(0, 0, +kCPlateLength/2, new TGeoRotation("", 0, 0, 0)));
	setVol->AddNode(vol, 0, new TGeoCombiTrans(0, +8.15*fgkmm, +kCPlateLength/2 /*+ 0.26*fgkmm*/, new TGeoRotation("", 0, 0, 0)));

	strObjName = strName+TString("_CPlateCnt");
	strObjTran = strObjName+TString("_tr");
	new TGeoBBox(strObjName, kCPlateWidth/2, kCPlateThick/2, kCPlateLength/2);	
	(new TGeoCombiTrans(strObjTran, 0, 0, +kCPlateLength/2, new TGeoRotation("", 0, 0, 0)))->RegisterYourself();
	strVol += strObjName;	strVol += ":";	strVol += strObjTran;	
	
	// -- Cooling Pipe Hole  --	
	const Double_t kCPipeHoleRadius = 1.0*fgkmm;
	const Double_t kCPipeHoleExt = 2.0*fgkmm;
	const Double_t kCPipeHoleLength = kCPlateLength + kCPipeHoleExt;
	strObjName = strName+TString("_CPipeHoleCnt");
	strObjTran = strObjName+TString("_tr");
	strObjTran2 = strObjName+TString("_tr2");
	new TGeoTube(strObjName, 0, kCPipeHoleRadius, kCPipeHoleLength/2);
	(new TGeoCombiTrans(strObjTran, +5.0*fgkmm, 0, +kCPipeHoleLength/2-kCPipeHoleExt/2, new TGeoRotation("", 0, 0, 0)))->RegisterYourself();
	strVol += " - ";	strVol += strObjName;	strVol += ":";	strVol += strObjTran;
	(new TGeoCombiTrans(strObjTran2,-5.0*fgkmm, 0, +kCPipeHoleLength/2-kCPipeHoleExt/2, new TGeoRotation("", 0, 0, 0)))->RegisterYourself();
	strVol += " - ";	strVol += strObjName;	strVol += ":";	strVol += strObjTran2;
	
	// -- Cooling Pipe Nook --
	const Double_t kCPipeNookWidth = 2.5*fgkmm;
	const Double_t kCPipeNookThick = 1.5*fgkmm;
	const Double_t kCPipeNookExt = 2.0*fgkmm;
	const Double_t kCPipeNookLength = kCPlateLength + kCPipeNookExt;
	strObjName = strName+TString("_CPipeNookCnt");
	strObjTran = strObjName+TString("_tr");
	new TGeoBBox(strObjName, kCPipeNookWidth/2, kCPipeNookThick/2, kCPipeNookLength/2);
	(new TGeoCombiTrans(strObjTran, +kCPlateWidth/2 -kCPipeNookWidth/2 +0.5*fgkmm, 0, +kCPipeNookLength/2-kCPipeNookExt/2, new TGeoRotation("", 0, 0, 0)))->RegisterYourself();
	strVol += " - ";	strVol += strObjName;	strVol += ":";	strVol += strObjTran;
	strObjTran2 = strObjName+TString("_tr2");
	(new TGeoCombiTrans(strObjTran2, -kCPlateWidth/2 +kCPipeNookWidth/2 -0.5*fgkmm, 0, +kCPipeNookLength/2-kCPipeNookExt/2, new TGeoRotation("", 0, 0, 0)))->RegisterYourself();
	strVol += " - ";	strVol += strObjName;	strVol += ":";	strVol += strObjTran2;

	strVol += ")";
	vol = new TGeoVolume(strName+TString("_CoolingPlateCnt"), new TGeoCompositeShape("_", strVol), gGeoManager->GetMedium("ITS_ALUMINUM$"));
	vol->SetFillColor(5);	vol->SetLineColor(5);
	setVol->AddNode(vol, 0, 0);

	// -- Cooling Pipes
	const Double_t kCPipeThick = 0.4*fgkmm;
	const Double_t kCPipeInnerRadius = 1.2*fgkmm/2;
	const Double_t kCPipeCntLength = 51.24*fgkmm;
	vol = new TGeoVolume(strName+TString("_CPipeCntRSide"), new TGeoTube(kCPipeInnerRadius, kCPipeInnerRadius+kCPipeThick, kCPipeCntLength/2), gGeoManager->GetMedium("ITS_PUR$"));
	vol->SetFillColor(4);	vol->SetLineColor(4);
	setVol->AddNode(vol, 0, new TGeoCombiTrans(+5.0*fgkmm, 0, +kCPipeCntLength/2 -(2.5*fgkmm), new TGeoRotation("", 0, 0, 0)));
	vol = new TGeoVolume(strName+TString("_CPipeCntLSide"), new TGeoTube(kCPipeInnerRadius, kCPipeInnerRadius+kCPipeThick, kCPipeCntLength/2), gGeoManager->GetMedium("ITS_PUR$"));
	vol->SetFillColor(2);	vol->SetLineColor(2);
	setVol->AddNode(vol, 0, new TGeoCombiTrans(-5.0*fgkmm, 0, +kCPipeCntLength/2 -(2.5*fgkmm), new TGeoRotation("", 0, 0, 0)));

	
	// -- Cooling Plate Tail --
	const Double_t kCPlateTailWidth = 7.0*fgkmm;
	const Double_t kCPlateTailThick = 3.0*fgkmm;
	const Double_t kCPlateTailLength = 5.0*fgkmm;
	vol = new TGeoVolume(strName+TString("_CPlateCntTail"), new TGeoBBox(kCPlateTailWidth/2, kCPlateTailThick/2, kCPlateTailLength/2), gGeoManager->GetMedium("ITS_ALUMINUM$"));
	vol->SetFillColor(5);	vol->SetLineColor(5);
	setVol->AddNode(vol, 0, new TGeoCombiTrans(0, 0, +kCPlateLength+kCPlateTailLength/2, new TGeoRotation("", 0, 0, 0)));
	setVol->AddNode(vol, 0, new TGeoCombiTrans(0, +8.15*fgkmm, +kCPlateLength + kCPlateTailLength/2, new TGeoRotation("", 0, 0, 0)));

	// -- FMC cable head --
	//const Double_t kFmcCableHeadWidth = 8.0*fgkmm;
	//const Double_t kFmcCableHeadThick = 0.44*fgkmm;
	//const Double_t kFmcCableHeadLength = 3.5*fgkmm;
	//vol = new TGeoVolume(strName+TString("_FmcCableHead"), new TGeoBBox(kFmcCableHeadWidth/2, kFmcCableHeadThick/2, kFmcCableHeadLength/2), gGeoManager->GetMedium("ITS_KAPTON(POLYCH2)$"));
	//vol->SetFillColor(8);	vol->SetLineColor(8);
	//setVol->AddNode(vol, 0, new TGeoCombiTrans(0,  +3.18*fgkmm, +kFmcCableHeadLength/2+5.88*fgkmm, new TGeoRotation("", 0, 0, 0)));
	//setVol->AddNode(vol, 0, new TGeoCombiTrans(0, +11.42*fgkmm, +kFmcCableHeadLength/2+5.88*fgkmm, new TGeoRotation("", 0, 0, 0)));
	//setVol->AddNode(vol, 0, new TGeoCombiTrans(0,  -3.18*fgkmm, +kFmcCableHeadLength/2+5.88*fgkmm, new TGeoRotation("", 0, 0, 0)));

	// -- FMC cable --
	//const Double_t kFmcCableWidth = 15.0*fgkmm;
	//const Double_t kFmcCableThick = 0.44*fgkmm;
	//const Double_t kFmcCableLength = 17.83*fgkmm;
	//vol = new TGeoVolume(strName+TString("_FmcCable"), new TGeoBBox(kFmcCableWidth/2, kFmcCableThick/2, kFmcCableLength/2), gGeoManager->GetMedium("ITS_KAPTON(POLYCH2)$"));
	//vol->SetFillColor(8);	vol->SetLineColor(8);
	//setVol->AddNode(vol, 0, new TGeoCombiTrans(0,  +3.18*fgkmm, +kFmcCableLength/2+5.88*fgkmm+kFmcCableHeadLength, new TGeoRotation("", 0, 0, 0)));
	//setVol->AddNode(vol, 0, new TGeoCombiTrans(0, +11.42*fgkmm, +kFmcCableLength/2+5.88*fgkmm+kFmcCableHeadLength, new TGeoRotation("", 0, 0, 0)));
	//setVol->AddNode(vol, 0, new TGeoCombiTrans(0,  -3.18*fgkmm, +kFmcCableLength/2+5.88*fgkmm+kFmcCableHeadLength, new TGeoRotation("", 0, 0, 0)));
	
	// -- FPGA --
	//const Double_t kFPGAWidth = 16.0*fgkmm;
	//const Double_t kFPGAThick = 0.1*fgkmm;
	//const Double_t kFPGALength = 10.00*fgkmm;
	//vol = new TGeoVolume(strName+TString("_FPGA"), new TGeoBBox(kFPGAWidth/2, kFPGAThick/2, kFPGALength/2), gGeoManager->GetMedium("ITS_KAPTON(POLYCH2)$"));
	//vol->SetFillColor(43);	vol->SetLineColor(43);
	//setVol->AddNode(vol, 0, new TGeoCombiTrans(0,  +1.55*fgkmm, +kFPGALength/2, new TGeoRotation("", 0, 0, 0)));
	//setVol->AddNode(vol, 0, new TGeoCombiTrans(0,  +9.70*fgkmm, +kFPGALength/2, new TGeoRotation("", 0, 0, 0)));
	//setVol->AddNode(vol, 0, new TGeoCombiTrans(0,  -1.55*fgkmm, +kFPGALength/2, new TGeoRotation("", 0, 0, 0)));
	
	
	// -- TBC Wire  --	
	//const Double_t kTBCWireRadius = 1.0*fgkmm;
	//const Double_t kTBCWireLength = 25.0*fgkmm;
	//vol = new TGeoVolume(strName+TString("_TBCWire"), new TGeoTube(0, kTBCWireRadius, kTBCWireLength/2), gGeoManager->GetMedium("ITS_KAPTON(POLYCH2)$"));
	//vol->SetFillColor(6);	vol->SetLineColor(6);
	//setVol->AddNode(vol, 0, new TGeoCombiTrans(0, +4.5*fgkmm, +kTBCWireLength/2+24.38*fgkmm, new TGeoRotation("", 0, 0, 0)));
	//setVol->AddNode(vol, 0, new TGeoCombiTrans(0, -4.5*fgkmm, +kTBCWireLength/2+24.38*fgkmm, new TGeoRotation("", 0, 0, 0)));

	// -- Samtec Cables --
	setVol->AddNode(CreateSamtecCables(), 0, new TGeoCombiTrans( 0,  0,  0, new TGeoRotation("",0,0,0)));

	
	// This will offset only X-axis and Y-axis. the Z-axis is to be offset outside the function.
	TGeoVolume* ofsVol = new TGeoVolumeAssembly(strName);
	ofsVol->AddNode(setVol, 1, 
		new TGeoCombiTrans(fkALC_0334_ContactDX[fLayerNumber], fkALC_0334_ContactDY[fLayerNumber]+fkThermalConnOffsetY[fLayerNumber], 0, 
		new TGeoRotation("", 180.0, 0, 0)));
	
	fN_DCCNT_Created++;
	
	return ofsVol;
	//return newOffsetVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateInnerDCDC3(const TGeoManager *mgr){
	
	TString strName = (TString("DCDC")+=fLayerNumber) + (TString("_")+=(int)fN_DCDC_Created);
	TString strObjName, strObjTran;
	
	TGeoCombiTrans *cbtObj;
	TGeoVolume* setVol = new TGeoVolumeAssembly(strName);
	TGeoVolume* vol;


	TString strVol = TString("(");

	// -- Cooling Plate --
	const Double_t kCPlateWidth = 16.0*fgkmm;
	const Double_t kCPlateThick = 3.0*fgkmm;
	const Double_t kCPlateLength = 36.75*fgkmm;
	strObjName = strName+TString("_CPlate");
	strObjTran = strObjName+TString("_tr");
	new TGeoBBox(strObjName ,kCPlateWidth/2, kCPlateThick/2, kCPlateLength/2);	
	(new TGeoCombiTrans(strObjTran, 0, 0, +kCPlateLength/2, new TGeoRotation("", 0, 0, 0)))->RegisterYourself();
	strVol += strObjName;	strVol += ":";	strVol += strObjTran;

	// -- Cooling Pipe Hole R-Side (Cold Pipe) --		
	const Double_t kCPipeHoleRRadius = 1.0*fgkmm;
	const Double_t kCPipeHoleRExt = 2.0*fgkmm;
	const Double_t kCPipeHoleRLength = kCPlateLength + kCPipeHoleRExt;
	strObjName = strName+TString("_CPipeHoleR");
	strObjTran = strObjName+TString("_tr");
	new TGeoTube(strObjName, 0, kCPipeHoleRRadius, kCPipeHoleRLength/2);
	(new TGeoCombiTrans(strObjTran, +6.5*fgkmm, 0, +kCPipeHoleRLength/2-kCPipeHoleRExt/2, new TGeoRotation("", 0, 0, 0)))->RegisterYourself();
	strVol += " - ";	strVol += strObjName;	strVol += ":";	strVol += strObjTran;

	// -- Cooling Pipe Nook R-Side (Cold Pipe) --
	const Double_t kCPipeNookRWidth = 3.0*fgkmm;
	const Double_t kCPipeNookRThick = 1.8*fgkmm;
	const Double_t kCPipeNookRExt = 2.0*fgkmm;
	const Double_t kCPipeNookRLength = kCPlateLength + kCPipeNookRExt;
	strObjName = strName+TString("_CPipeNook");
	strObjTran = strObjName+TString("_tr");
	new TGeoBBox(strObjName, kCPipeNookRWidth/2, kCPipeNookRThick/2, kCPipeNookRLength/2);
	(new TGeoCombiTrans(strObjTran, +kCPlateWidth/2, 0, +kCPipeNookRLength/2-kCPipeNookRExt/2, new TGeoRotation("", 0, 0, 0)))->RegisterYourself();
	strVol += " - ";	strVol += strObjName;	strVol += ":";	strVol += strObjTran;

	// -- Combine: Cooling Plate - ( Pipe Hole + Pipe Nook ) -- (Cold Pipe)
	strVol += ")";
	vol = new TGeoVolume(strName+TString("_CoolingPlate"), new TGeoCompositeShape("_", strVol), gGeoManager->GetMedium("ITS_ALUMINUM$"));
	vol->SetFillColor(11);	vol->SetLineColor(11);
	setVol->AddNode(vol, 0, 0);
	
	// -- Cooling Pipe R-Side (Cold Pipe) --
	const Double_t kCPipeThick = 0.4*fgkmm;
	const Double_t kCPipeInnerRadius = 1.2*fgkmm/2;
	const Double_t kCPipePlateRSideLength = 50.85*fgkmm + (11.0*fgkmm); // The same length as kPCBPlateUpLength - 21.0, 11.0 is just an extend.
	vol = new TGeoVolume(strName+TString("_CPipeRSide"), new TGeoTube(kCPipeInnerRadius, kCPipeInnerRadius+kCPipeThick, kCPipePlateRSideLength/2), gGeoManager->GetMedium("ITS_PUR$"));
	vol->SetFillColor(4);	vol->SetLineColor(4);
	setVol->AddNode(vol, 0, new TGeoCombiTrans(+6.5*fgkmm, 0, +kCPipePlateRSideLength/2 - (11.0*fgkmm), new TGeoRotation("", 0, 0, 0)));

	// -- Cooling Pipe L-Side (Hot Pipe) -- Head
	const Double_t kCPipePlateLSideHeadLength = 10.0*fgkmm; 
	vol = new TGeoVolume(strName+TString("_CPipeLSideHead"), new TGeoTube(kCPipeInnerRadius, kCPipeInnerRadius+kCPipeThick, kCPipePlateLSideHeadLength/2), gGeoManager->GetMedium("ITS_PUR$"));
	vol->SetFillColor(2);	vol->SetLineColor(2);
	setVol->AddNode(vol, 0, new TGeoCombiTrans(-6.5*fgkmm, 0, +kCPipePlateLSideHeadLength/2 - (11.0*fgkmm), new TGeoRotation("", 0, 0, 0)));

	// -- Cooling Pipe L-Side (Hot Pipe) -- Tail
	const Double_t kCPipePlateLSideTailLength = 13.10*fgkmm; // 50.85 - kCPlateLength - 1.0(offset)
	vol = new TGeoVolume(strName+TString("_CPipeLSideTail"), new TGeoTube(kCPipeInnerRadius, kCPipeInnerRadius+kCPipeThick, kCPipePlateLSideTailLength/2), gGeoManager->GetMedium("ITS_PUR$"));
	vol->SetFillColor(2);	vol->SetLineColor(2);
	setVol->AddNode(vol, 0, new TGeoCombiTrans(-6.5*fgkmm, 0, +kCPipePlateLSideTailLength/2 + kCPlateLength + (1.0*fgkmm), new TGeoRotation("", 0, 0, 0)));
		
		
	// -- PCB Plate --
	const Double_t kPCBPlateWidth = 16.0*fgkmm;
	const Double_t kPCBPlateThick = 0.1*fgkmm;
	const Double_t kPCBPlateUpLength = 71.85*fgkmm;
	vol = new TGeoVolume(strName+TString("_PCBPlateUp"), new TGeoBBox(kPCBPlateWidth/2, kPCBPlateThick/2, kPCBPlateUpLength/2), gGeoManager->GetMedium("ITS_G10FR4CU$"));
	vol->SetFillColor(5);	vol->SetLineColor(5);
	setVol->AddNode(vol, 0, new TGeoCombiTrans(0, +kCPlateThick/2+kPCBPlateThick/2, +kPCBPlateUpLength/2-21.0*fgkmm, new TGeoRotation("", 0, 0, 0)));	
	const Double_t kPCBPlateDnLength = 58.20*fgkmm;
	vol = new TGeoVolume(strName+TString("_PCBPlateDn"), new TGeoBBox(kPCBPlateWidth/2, kPCBPlateThick/2, kPCBPlateDnLength/2), gGeoManager->GetMedium("ITS_G10FR4CU$"));
	vol->SetFillColor(5);	vol->SetLineColor(5);
	setVol->AddNode(vol, 0, new TGeoCombiTrans(0, -kCPlateThick/2-kPCBPlateThick/2, +kPCBPlateDnLength/2-7.45*fgkmm, new TGeoRotation("", 0, 0, 0)));
	
	// -- DCDC --
	const Double_t kDCDCWidth = 14.0*fgkmm;
	const Double_t kDCDCThick = 8.0*fgkmm;
	const Double_t kDCDCLength = 18.50*fgkmm;
	vol = new TGeoVolume(strName+TString("_DCDC"), new TGeoBBox(kDCDCWidth/2, kDCDCThick/2, kDCDCLength/2), gGeoManager->GetMedium("ITS_DCDCSHIELD$"));
	vol->SetFillColor(6); vol->SetLineColor(6);
	setVol->AddNode(vol, 0, new TGeoCombiTrans(0, +kCPlateThick/2+kPCBPlateThick+kDCDCThick/2, +kDCDCLength/2, new TGeoRotation("", 0, 0, 0)));
	setVol->AddNode(vol, 0, new TGeoCombiTrans(0, -kCPlateThick/2-kPCBPlateThick-kDCDCThick/2, +kDCDCLength/2+18.25*fgkmm, new TGeoRotation("", 0, 0, 0)));

	//vol = new TGeoVolume(strName+TString("_DCDCDn"), new TGeoBBox(kDCDCWidth/2, kDCDCThick/2, kDCDCLength/2), gGeoManager->GetMedium("ITS_DCDCSHIELD$"));
	//vol->SetFillColor(6);	vol->SetLineColor(6);
	//setVol->AddNode(vol, 0, new TGeoCombiTrans(0, -kCPlateThick/2-kPCBPlateThick-kDCDCThick/2, +kDCDCLength/2+18.25*fgkmm, new TGeoRotation("", 0, 0, 0)));

	// -- DCDC Connectors / DCDCx --
	const Double_t kDCDCxWidth = 16.0*fgkmm;
	const Double_t kDCDCxThick = 2.3*fgkmm;
	const Double_t kDCDCxLength = 9.0*fgkmm;
	vol = new TGeoVolume(strName+TString("_DCDCx"), new TGeoBBox(kDCDCxWidth/2, kDCDCxThick/2, kDCDCxLength/2), gGeoManager->GetMedium("ITS_DCDCPASSCNT$"));
	vol->SetFillColor(4);	vol->SetLineColor(4);
	setVol->AddNode(vol, 0, new TGeoCombiTrans(0, +kCPlateThick/2+kPCBPlateThick+kDCDCxThick/2, +kDCDCxLength/2-9.0*fgkmm, new TGeoRotation("", 0, 0, 0)));
	setVol->AddNode(vol, 0, new TGeoCombiTrans(0, -kCPlateThick/2-kPCBPlateThick-kDCDCxThick/2, +kDCDCxLength/2+9.25*fgkmm, new TGeoRotation("", 0, 0, 0)));

	//vol = new TGeoVolume(strName+TString("_DCDCxDn"), new TGeoBBox(kDCDCxWidth/2, kDCDCxThick/2, kDCDCxLength/2), gGeoManager->GetMedium("ITS_DCDCPASSCNT$"));
	//vol->SetFillColor(6);	vol->SetLineColor(6);
	//setVol->AddNode(vol, 0, new TGeoCombiTrans(0, -kCPlateThick/2-kPCBPlateThick-kDCDCxThick/2, +kDCDCxLength/2+9.25*fgkmm, new TGeoRotation("", 0, 0, 0)));


	// -- FPGA Connector -- DO NOT USE IN THE LASTEST DESIGN
	//const Double_t kCntOffsetZ = +kPipeLength/2 -7.0*fgkmm;
	//const Double_t kCntOffsetY = 8.4*fgkmm +kWidthBtwnUDPlates/2 +kUpPlateThick -0.5*fgkmm; // 8.4 =DCBox, 0.5 = obj thick
	//const Double_t kCntOffsetZ2 = 61.0*fgkmm;
	//setVol->AddNode(CreateFPGAConnector("1"), 1, new TGeoTranslation(0, +kCntOffsetY, +kCntOffsetZ));
	//setVol->AddNode(CreateFPGAConnector("2"), 1, new TGeoCombiTrans(0, +kCntOffsetY, -kCntOffsetZ2 +2*fgkmm, new TGeoRotation("",0,12,0)));
	//setVol->AddNode(CreateFPGAConnector("3"), 1, new TGeoCombiTrans(0, +kCntOffsetY/2, -kCntOffsetZ2, new TGeoRotation("",0,12,0)));
	
	//TGeoVolume* newOffsetVol = new TGeoVolumeAssembly(TString("DCPS")+=(int)fLayerNumber);
	//newOffsetVol->AddNode(setVol, 1, new TGeoTranslation(0, 0, 0));
	//vALIC->AddNode(setVol, 1, new TGeoCombiTrans(0,0,0,new TGeoRotation("",0,0,0)));
	
	// This will offset only X-axis and Y-axis. 
	// The Z-axis is to be offset outside the function.
	TGeoVolume* ofsVol = new TGeoVolumeAssembly(strName);
	ofsVol->AddNode(setVol, 1, 
		new TGeoCombiTrans(fkALC_0334_ContactDX[fLayerNumber], fkALC_0334_ContactDY[fLayerNumber]+fkDCDCOffsetY[fLayerNumber], 0, 
		new TGeoRotation("", 180.0, fkDCDCRotateBeta[fLayerNumber], 0)));
	
	fN_DCDC_Created++;
	
	return ofsVol;
	//return newOffsetVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateDCPS2(const char *nid, const TGeoManager *mgr){

	TGeoVolume* setVol = new TGeoVolumeAssembly(TString("DCPS")+=nid);
	TGeoVolume* vol;

	// -- Pipes --	
	const Double_t kPipeHoleRadius = 0.6*fgkmm;
	const Double_t kPipeThick = 0.2*fgkmm;
	const Double_t kPipeLength = 77.0*fgkmm;
	const Double_t kWidthBtwnPipes = 12.4*fgkmm;
	vol = new TGeoVolume("CoolingPipe", new TGeoTube(kPipeHoleRadius, kPipeHoleRadius+kPipeThick, kPipeLength/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	vol->SetFillColor(5);	vol->SetLineColor(5);
	setVol->AddNode(vol, 0, new TGeoTranslation(- kWidthBtwnPipes/2, 0, 0));
	setVol->AddNode(vol, 0, new TGeoTranslation(+ kWidthBtwnPipes/2, 0, 0));

	// -- Plate --
	const Double_t kPlateLength = 54.0*fgkmm;
	const Double_t kPlateThick = 2*(kPipeHoleRadius + kPipeThick);
	const Double_t kPlateWidth = kWidthBtwnPipes - 2*(kPipeHoleRadius + kPipeThick); // instead of using the composite object
	const Double_t kPlateOffset = 2.5*fgkmm;
	vol = new TGeoVolume("CoolingPlate", new TGeoBBox(kPlateWidth/2, kPlateThick/2, kPlateLength/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	vol->SetFillColor(5);	vol->SetLineColor(5);
	setVol->AddNode(vol, 0, new TGeoTranslation(0, 0, -kPlateOffset));

	// -- Up and Down Plates
	const Double_t kWidthBtwnUDPlates = 1.7*fgkmm;

	// -- Up Plate --
	const Double_t kUpPlateLength = 76.8*fgkmm;
	const Double_t kUpPlateThick = 0.1*fgkmm;
	const Double_t kUpPlateWidth = 16.0*fgkmm;
	const Double_t kUpPlateOffset = 4.5*fgkmm;
	vol = new TGeoVolume("UpPlate", new TGeoBBox(kUpPlateWidth/2, kUpPlateThick/2, kUpPlateLength/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	vol->SetFillColor(2);	vol->SetLineColor(2);
	setVol->AddNode(vol, 0, new TGeoTranslation(0, +kWidthBtwnUDPlates/2+kUpPlateThick/2, -(kPipeLength/2-kUpPlateLength/2)-kUpPlateOffset));

	// -- Up Extend Plate 0 --
	const Double_t kUpEPlateLength = 32.3*fgkmm;
	const Double_t kUpEPlateThick = 0.1*fgkmm;
	const Double_t kUpEPlateWidth = 16.0*fgkmm;
	vol = new TGeoVolume("UpEPlate", new TGeoBBox(kUpEPlateWidth/2, kUpEPlateThick/2, kUpEPlateLength/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	vol->SetFillColor(2);	vol->SetLineColor(2);
	TGeoVolume* asmVol = new TGeoVolumeAssembly("");
	asmVol->AddNode(vol, 0, new TGeoTranslation(0, -kUpEPlateThick/2, +kUpEPlateLength/2));
	setVol->AddNode(asmVol, 0, new TGeoCombiTrans(0, +kWidthBtwnUDPlates/2+kUpPlateThick/2 +kUpEPlateThick/2, -kPipeLength/2 -kUpPlateOffset +kUpPlateLength, new TGeoRotation("", 0, -13, 0)));
	
	// -- Dn Plate --
	const Double_t kDnPlateLength = 69.2*fgkmm;
	const Double_t kDnPlateThick = 0.1*fgkmm;
	const Double_t kDnPlateWidth = 16.0*fgkmm;
	const Double_t kDnPlateOffset = 4.5*fgkmm;
	vol = new TGeoVolume("DnPlate", new TGeoBBox(kDnPlateWidth/2, kDnPlateThick/2, kDnPlateLength/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	vol->SetFillColor(2);	vol->SetLineColor(2);
	setVol->AddNode(vol, 0, new TGeoTranslation(0, -kWidthBtwnUDPlates/2-kDnPlateThick/2, -(kPipeLength/2-kDnPlateLength/2)-kDnPlateOffset));
	
	// -- Power Supply Box
	const Double_t kUpBoxOffset = +kUpPlateThick +kWidthBtwnUDPlates/2;
	const Double_t kDnBoxOffset = +kDnPlateThick +kWidthBtwnUDPlates/2;
	//	setVol->AddNode(CreateDCBox2("Up"), 1, new TGeoTranslation(0, +kWidthBtwnUDPlates/2 +kUpPlateThick, -kPipeLength/2 +kUpBoxOffset));
	//	setVol->AddNode(CreateDCBox2("Dn"), 1, new TGeoCombiTrans(0, -kWidthBtwnUDPlates/2 -kDnPlateThick, -kPipeLength/2 +kDnBoxOffset,new TGeoRotation("", 0, 0, 180) ));
	
	// -- Box --
	const Double_t kBoxDX = 14.0*fgkmm/2;
	const Double_t kBoxDY = 07.0*fgkmm/2;
	const Double_t kBoxDZ = 18.5*fgkmm/2;
	const Double_t kBoxUpOffset = 13.5*fgkmm; // offset from the A-Side of the edge of the pipe
	const Double_t kBoxDnOffset = 22.5*fgkmm; // offset from the A-Side of the edge of the pipe
	vol = new TGeoVolume("Box", new TGeoBBox(kBoxDX, kBoxDY, kBoxDZ), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	vol->SetFillColor(4);	vol->SetLineColor(4);
	setVol->AddNode(vol, 0, new TGeoTranslation(0, +kBoxDY +kUpBoxOffset, +kBoxDZ +kBoxUpOffset -kPipeLength/2));
	setVol->AddNode(vol, 0, new TGeoTranslation(0, -kBoxDY -kDnBoxOffset, +kBoxDZ +kBoxDnOffset -kPipeLength/2));
	
	// -- Bar --
	const Double_t kBarDX = 16.0*fgkmm/2;
	const Double_t kBarDY = 03.0*fgkmm/2;
	const Double_t kBarDZ = 09.0*fgkmm/2;	
	const Double_t kBarUpOffset = kBoxUpOffset + kBoxDZ*2; // offset from the A-Side of the edge of the pipe
	const Double_t kBarDnOffset = kBoxDnOffset + kBoxDZ*2; // offset from the A-Side of the edge of the pipe
	vol = new TGeoVolume("Bar", new TGeoBBox(kBarDX, kBarDY, kBarDZ), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	vol->SetFillColor(9);	vol->SetLineColor(9);
	setVol->AddNode(vol, 0, new TGeoTranslation(0, +kBarDY +kUpBoxOffset, +kBarDZ +kBarUpOffset -kPipeLength/2));
	setVol->AddNode(vol, 0, new TGeoTranslation(0, -kBarDY -kDnBoxOffset, +kBarDZ +kBarDnOffset -kPipeLength/2));

	// -- E Plate 1 --
	const Double_t kEPlate1Length = 49.0*fgkmm;
	const Double_t kEPlate1Thick = 0.1*fgkmm;
	const Double_t kEPlate1Width = 16.0*fgkmm;
	const Double_t kEPlate1Offset = 4.5*fgkmm;
	vol = new TGeoVolume("EPlate1", new TGeoBBox(kEPlate1Width/2, kEPlate1Thick/2, kEPlate1Length/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	vol->SetFillColor(9);
	vol->SetLineColor(9);
	setVol->AddNode(vol, 0, new TGeoTranslation(0, +kWidthBtwnUDPlates/2+kEPlate1Thick+(kEPlate1Thick/2)+kBoxDY*2, -(kPipeLength/2-kEPlate1Length/2)-kEPlate1Offset));

	// -- E Plate 2 --
	const Double_t kEPlate2Length = 20.0*fgkmm;
	const Double_t kEPlate2Thick = 0.1*fgkmm;
	const Double_t kEPlate2Width = 16.0*fgkmm;
	const Double_t kEPlate2Offset = 4.5*fgkmm;
	vol = new TGeoVolume("EPlate2", new TGeoBBox(kEPlate2Width/2, kEPlate2Thick/2, kEPlate2Length/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	vol->SetFillColor(9);
	vol->SetLineColor(9);
	setVol->AddNode(vol, 0, new TGeoTranslation(0, +kWidthBtwnUDPlates/2+kEPlate2Thick/2+kEPlate2Thick/2+kBoxDY*2, -(kPipeLength/2-kEPlate2Length/2)-kEPlate2Offset+kEPlate1Length));
	
	// -- Pipes
	//const Double_t kPipesOffsetZ = 0*fgkmm;
	//const Double_t kPipesOffsetY = 8.5*fgkmm +kWidthBtwnUDPlates/2 +kUpPlateThick; // 8.5 = DCBox = 8.0 + 0.4 + empty space
	//setVol->AddNode(CreateDCPipes(), 1, new TGeoTranslation(0, +kPipesOffsetY, +kPipesOffsetZ));
	
	// -- FPGA Connector
	const Double_t kCntOffsetZ = +kPipeLength/2 -7.0*fgkmm;
	const Double_t kCntOffsetY = 8.4*fgkmm +kWidthBtwnUDPlates/2 +kUpPlateThick -0.5*fgkmm; // 8.4 =DCBox, 0.5 = obj thick
	const Double_t kCntOffsetZ2 = 61.0*fgkmm;
	setVol->AddNode(CreateFPGAConnector("1"), 1, new TGeoTranslation(0, +kCntOffsetY, +kCntOffsetZ));
	setVol->AddNode(CreateFPGAConnector("2"), 1, new TGeoCombiTrans(0, +kCntOffsetY, -kCntOffsetZ2 +2*fgkmm, new TGeoRotation("",0,12,0)));
	setVol->AddNode(CreateFPGAConnector("3"), 1, new TGeoCombiTrans(0, +kCntOffsetY/2, -kCntOffsetZ2, new TGeoRotation("",0,12,0)));
	
	TGeoVolume* newOffsetVol = new TGeoVolumeAssembly(TString("DCPS")+=nid);
	newOffsetVol->AddNode(setVol, 1, new TGeoTranslation(0, 0, -kPipeLength/2));
	
	
	
	return newOffsetVol;
}

//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateSamtecCables(const TGeoManager *mgr){
	
	TString strName = (TString("SamtecCables")+=fLayerNumber) + (TString("_")+=(int)fN_SamtecCables_Created);
	TGeoVolume* setVol = new TGeoVolumeAssembly(strName);
	TGeoVolume* vol;

	setVol->AddNode(CreateSamtecCable(), 0, new TGeoCombiTrans( 4.12*fgkmm,  12.63*fgkmm,  11.60*fgkmm, new TGeoRotation("",0,0,0)));
	setVol->AddNode(CreateSamtecCable(), 0, new TGeoCombiTrans( 2.47*fgkmm,  12.63*fgkmm,  11.60*fgkmm, new TGeoRotation("",0,0,0)));
	setVol->AddNode(CreateSamtecCable(), 0, new TGeoCombiTrans( 0.82*fgkmm,  12.63*fgkmm,  11.60*fgkmm, new TGeoRotation("",0,0,0)));
	setVol->AddNode(CreateSamtecCable(), 0, new TGeoCombiTrans(-0.83*fgkmm,  12.63*fgkmm,  11.60*fgkmm, new TGeoRotation("",0,0,0)));
	setVol->AddNode(CreateSamtecCable(), 0, new TGeoCombiTrans(-2.48*fgkmm,  12.63*fgkmm,  11.60*fgkmm, new TGeoRotation("",0,0,0)));
	setVol->AddNode(CreateSamtecCable(), 0, new TGeoCombiTrans(-4.13*fgkmm,  12.63*fgkmm,  11.60*fgkmm, new TGeoRotation("",0,0,0)));

	setVol->AddNode(CreateSamtecCable(), 0, new TGeoCombiTrans( 4.12*fgkmm,  13.53*fgkmm,  11.60*fgkmm, new TGeoRotation("",0,0,0)));
	setVol->AddNode(CreateSamtecCable(), 0, new TGeoCombiTrans( 2.47*fgkmm,  13.53*fgkmm,  11.60*fgkmm, new TGeoRotation("",0,0,0)));
	setVol->AddNode(CreateSamtecCable(), 0, new TGeoCombiTrans( 0.82*fgkmm,  13.53*fgkmm,  11.60*fgkmm, new TGeoRotation("",0,0,0)));
	setVol->AddNode(CreateSamtecCable(), 0, new TGeoCombiTrans(-0.83*fgkmm,  13.53*fgkmm,  11.60*fgkmm, new TGeoRotation("",0,0,0)));
	setVol->AddNode(CreateSamtecCable(), 0, new TGeoCombiTrans(-2.48*fgkmm,  13.53*fgkmm,  11.60*fgkmm, new TGeoRotation("",0,0,0)));
	setVol->AddNode(CreateSamtecCable(), 0, new TGeoCombiTrans(-4.13*fgkmm,  13.53*fgkmm,  11.60*fgkmm, new TGeoRotation("",0,0,0)));

	fN_SamtecCables_Created++;
	
	return setVol;
}



//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateSamtecCable(const TGeoManager *mgr){
	
	TString strName = (TString("Samtec")+=fLayerNumber) + (TString("_")+=(int)fN_SamtecCable_Created);
	TGeoVolume* setVol = new TGeoVolumeAssembly(strName);
	TGeoVolume* vol;
	
	TGeoMedium *medInner      = mgr->GetMedium("ITS_COPPER$");
	TGeoMedium *medOuter	  = mgr->GetMedium("ITS_PUR$");

	const Double_t kHoleRadius = 0.125*fgkmm; 
	
	const Double_t kPipeThick = 0.325*fgkmm;
	const Double_t kPipeLength = 25.0*fgkmm;
	const Double_t kWidthBtwnHoles = 0.65*fgkmm;
	
	// -- Middle Cable Wall
	vol = new TGeoVolume(strName+TString("_MidleCableWall"), new TGeoBBox(kWidthBtwnHoles/2-kHoleRadius, kHoleRadius, kPipeLength/2), medOuter);
	vol->SetFillColor(41);	vol->SetLineColor(41);
	setVol->AddNode(vol, 0, new TGeoCombiTrans(0, 0, +kPipeLength/2, new TGeoRotation("",0,0,0)));
	
	// -- Outer Cable Wall --
	vol = new TGeoVolume(strName+TString("_OuterCableWall"), new TGeoBBox(kWidthBtwnHoles/2, kPipeThick/2, kPipeLength/2), medOuter);
	vol->SetFillColor(41);	vol->SetLineColor(41);
	setVol->AddNode(vol, 0, new TGeoTranslation(0, +kHoleRadius +kPipeThick/2, +kPipeLength/2));
	setVol->AddNode(vol, 0, new TGeoTranslation(0, -kHoleRadius -kPipeThick/2, +kPipeLength/2));
	
	// -- Outer cable --
	TGeoPcon *volObj = new TGeoPcon(strName+TString("_OuterCableCurveTube"), 0, 180, 2);
	volObj->DefineSection(0, 0.0		, kHoleRadius, kHoleRadius + kPipeThick);
	volObj->DefineSection(1, kPipeLength, kHoleRadius, kHoleRadius + kPipeThick);
	vol = new TGeoVolume(strName+TString("_OuterCableCurve"), volObj, medOuter);
	vol->SetFillColor(41);	vol->SetLineColor(41);
	setVol->AddNode(vol, 0, new TGeoCombiTrans(-kWidthBtwnHoles/2, 0, 0, new TGeoRotation("",0,0,90)));
	setVol->AddNode(vol, 0, new TGeoCombiTrans(+kWidthBtwnHoles/2, 0, 0, new TGeoRotation("",0,0,-90)));

	// -- Inner cable --
	volObj = new TGeoPcon(strName+TString("_InnerCableTube"), 0, 360, 2);
	volObj->DefineSection(0, 0.0		, 0, kHoleRadius);
	volObj->DefineSection(1, kPipeLength, 0, kHoleRadius);
	vol = new TGeoVolume(strName+TString("_InnerCable"), volObj, medInner);
	vol->SetFillColor(46);	vol->SetLineColor(46);
	setVol->AddNode(vol, 0, new TGeoCombiTrans(-kWidthBtwnHoles/2, 0, 0, new TGeoRotation("",0,0,0)));
	setVol->AddNode(vol, 0, new TGeoCombiTrans(+kWidthBtwnHoles/2, 0, 0, new TGeoRotation("",0,0,0)));

	
	fN_SamtecCable_Created++;
	
	return setVol;
}



//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateDCPipe(const char *nid, const TGeoManager *mgr){
	
	TGeoVolume* unitVol = new TGeoVolumeAssembly(TString("DCPipe")+=nid);
	TGeoVolume* volume;
	
	const Double_t kHoleRadius = 0.1*fgkmm; 
	
	const Double_t kPipeThick = 0.35*fgkmm;
	const Double_t kPipeLength = 75.0*fgkmm;
	const Double_t kWidthBtwnHoles = 0.65*fgkmm;
	const Double_t kOffsetCenterY = kHoleRadius + kPipeThick;
	// -- Pipe Wall --
	volume = new TGeoVolume("PipeWall", new TGeoBBox(kWidthBtwnHoles/2, kPipeThick/2, kPipeLength/2), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	volume->SetFillColor(5);
	volume->SetLineColor(5);
	unitVol->AddNode(volume, 0, new TGeoTranslation(0, +kHoleRadius +kPipeThick/2 +kOffsetCenterY, 0));
	unitVol->AddNode(volume, 0, new TGeoTranslation(0, -kHoleRadius -kPipeThick/2 +kOffsetCenterY, 0));
	// -- Pipe Curve --
	TGeoPcon *volObj = new TGeoPcon("PipeCurveTube", 0, 180, 2);
	volObj->DefineSection(0, 0.0		, kHoleRadius, kHoleRadius + kPipeThick);
	volObj->DefineSection(1, kPipeLength, kHoleRadius, kHoleRadius + kPipeThick);
	volume = new TGeoVolume("PipeCurve", volObj, mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	volume->SetFillColor(5);
	volume->SetLineColor(5);
	unitVol->AddNode(volume, 0, new TGeoCombiTrans(-kWidthBtwnHoles/2, +kOffsetCenterY, -kPipeLength/2, new TGeoRotation("",0,0,90)));
	unitVol->AddNode(volume, 0, new TGeoCombiTrans(+kWidthBtwnHoles/2, +kOffsetCenterY, -kPipeLength/2, new TGeoRotation("",0,0,-90)));
	
	return unitVol;
}


TGeoVolume* AliITSUv2Layer::CreateDCPipes(const char *nid, const TGeoManager *mgr){
	
	TGeoVolume* unitVol = new TGeoVolumeAssembly(TString("DCPipes")+=nid);
	TGeoVolume* volume;
	
	// -- Pipes --
	const Double_t kPipeLX = 1.55*fgkmm;
	const Double_t kPipeLY = 0.90*fgkmm;
	const Double_t kPipeDX = kPipeLX/2;
	const Double_t kPipeDY = kPipeLY/2;
	const Double_t kHorGap = 0.05*fgkmm;
	const Double_t kVerGap = 0.02*fgkmm;
	const Double_t kHorWidthIter = kPipeLX + kHorGap;
	const Double_t kVerWidthIter = kPipeLY + kVerGap;
	// -- Bottom Row --
	unitVol->AddNode(volume = CreateDCPipe("3"), 1, new TGeoTranslation(0, +kPipeDY, 0));
	unitVol->AddNode(volume = CreateDCPipe("2"), 1, new TGeoTranslation(-kHorWidthIter, +kPipeDY, 0));
	unitVol->AddNode(volume = CreateDCPipe("4"), 1, new TGeoTranslation(+kHorWidthIter, +kPipeDY, 0));
	unitVol->AddNode(volume = CreateDCPipe("1"), 1, new TGeoTranslation(-2*kHorWidthIter, +kPipeDY, 0));
	unitVol->AddNode(volume = CreateDCPipe("5"), 1, new TGeoTranslation(+2*kHorWidthIter, +kPipeDY, 0));
	// -- Upper Row
	unitVol->AddNode(volume = CreateDCPipe("9"), 1, new TGeoTranslation(-kHorWidthIter/2, +kVerWidthIter +kPipeDY, 0));
	unitVol->AddNode(volume = CreateDCPipe("8"), 1, new TGeoTranslation(+kHorWidthIter/2, +kVerWidthIter +kPipeDY, 0));
	unitVol->AddNode(volume = CreateDCPipe("10"), 1, new TGeoTranslation(-3*kHorWidthIter/2 , +kVerWidthIter +kPipeDY, 0));
	unitVol->AddNode(volume = CreateDCPipe("7"), 1, new TGeoTranslation(+3*kHorWidthIter/2, +kVerWidthIter +kPipeDY, 0));
	unitVol->AddNode(volume = CreateDCPipe("11"), 1, new TGeoTranslation(-5*kHorWidthIter/2 , +kVerWidthIter +kPipeDY, 0));
	unitVol->AddNode(volume = CreateDCPipe("6"), 1, new TGeoTranslation(+5*kHorWidthIter/2, +kVerWidthIter +kPipeDY, 0));

	return unitVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateFPGAConnector(const char *nid, const TGeoManager *mgr){
	
	TGeoVolume* unitVol = new TGeoVolumeAssembly(TString("FPGACnt")+=nid);
	TGeoVolume* volume;
	
	const Double_t kDY = 1.0*fgkmm/2;
	
	const Double_t kB1DX = 11.0*fgkmm/2;
	const Double_t kB1DZ = 7.0*fgkmm/2;
	volume = new TGeoVolume("B1", new TGeoBBox(kB1DX, kDY, kB1DZ), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	unitVol->AddNode(volume, 0, new TGeoTranslation(0, 0, kB1DZ));
	
	const Double_t kB2DX = 8.0*fgkmm/2;
	const Double_t kB2DZ = 5.0*fgkmm/2;
	volume = new TGeoVolume("B2", new TGeoBBox(kB2DX, kDY, kB2DZ), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	unitVol->AddNode(volume, 0, new TGeoTranslation(0, 0, kB1DZ*2 +kB2DZ));

	const Double_t kB3DX = 11.0*fgkmm/2;
	const Double_t kB3DZ = 5.0*fgkmm/2;
	volume = new TGeoVolume("B3", new TGeoBBox(kB3DX, kDY, kB3DZ), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	unitVol->AddNode(volume, 0, new TGeoTranslation(0, 0, kB1DZ*2 +kB2DZ*2 + kB3DZ));

	const Double_t kB4DX = 10.0*fgkmm/2;
	const Double_t kB4DZ = 1*fgkmm/2;
	volume = new TGeoVolume("B4", new TGeoBBox(kB4DX, kDY, kB4DZ), mgr->GetMedium("ITS_KAPTON(POLYCH2)$"));
	unitVol->AddNode(volume, 0, new TGeoTranslation(0, 0, kB1DZ*2 +kB2DZ*2 +kB3DZ*2 +kB4DZ));
	
	return unitVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateSTube45(
	TString volName, 
	Double_t radius1, 
	Double_t radius2, 	// radius2 > radius1
	Double_t length, 
	Double_t length1, 
	Double_t thick,
	Bool_t isPosExt, 	// Is Positive Extrude ? 
	Int_t color, 
	const TGeoMedium* med)
{
	TGeoPcon *volObj = new TGeoPcon("volObj", 0, 180, 6);
	volObj->DefineSection(0, 0.0																		, radius2			, radius2+thick);
	volObj->DefineSection(1, (isPosExt?(1):(-1))*(length - length1 - (radius2 - radius1) - (thick/2))	, radius2			, radius2+thick);
	volObj->DefineSection(2, (isPosExt?(1):(-1))*(length - length1 - (radius2 - radius1) + (thick/2))	, radius2 - thick	, radius2+thick);
	volObj->DefineSection(3, (isPosExt?(1):(-1))*(length - length1 - (thick/2))							, radius1			, radius1+(2*thick));
	volObj->DefineSection(4, (isPosExt?(1):(-1))*(length - length1 + (thick/2))							, radius1			, radius1+thick);
	volObj->DefineSection(5, (isPosExt?(1):(-1))*(length)												, radius1			, radius1+thick);			

	TGeoVolume *volMed = new TGeoVolume(volName, volObj, med);
	volMed->SetFillColor(color);
	volMed->SetLineColor(color);

	TGeoVolume *volAsm = new TGeoVolumeAssembly(volName);
	volAsm->AddNode(volMed, 1, 0);
	return volAsm;
}


//______________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateTube(
	TString volName, 
	Double_t radius, 
	Double_t length, 
	Double_t thick, 
	Bool_t isPosExt, 	// Is Positive Extrude ? 
	Int_t color, 
	Double_t r,
	const TGeoMedium* med)
{
	
	TGeoPcon *volObj = new TGeoPcon("volObj", 0, r, 2);
	volObj->DefineSection(0, 0.0						, radius	, radius + thick);
	volObj->DefineSection(1, (isPosExt?1:-1)*(length)	, radius	, radius + thick);

	TGeoVolume *volMed = new TGeoVolume(volName, volObj, med);
	volMed->SetFillColor(color);
	volMed->SetLineColor(color);

	TGeoVolume *volAsm = new TGeoVolumeAssembly(volName);
	volAsm->AddNode(volMed, 1, 0);
	return volAsm;

}


//______________________________________________________________________
TGeoVolume* AliITSUv2Layer::CreateLTube(
	TString volName, 
	Double_t radius, 
	Double_t length1, 
	Double_t length2, 
	Double_t thick1, 
	Double_t thick2,
	Bool_t isHeadIn,
	Bool_t isPosExt, 	// Is Positive Extrude ? 
	Int_t color, 
	const TGeoMedium* med)
{
	
	TGeoPcon *volObj = new TGeoPcon("volObj", 0, 180, 4);
	volObj->DefineSection(0, 0.0								, radius							, radius + thick1 + (isHeadIn?1:0)*thick2);
	volObj->DefineSection(1, (isPosExt?1:-1)*(length1)			, radius							, radius + thick1 + (isHeadIn?1:0)*thick2);
	volObj->DefineSection(2, (isPosExt?1:-1)*(length1)			, radius + (isHeadIn?1:0)*thick1	, radius + (isHeadIn?1:0)*thick1 + thick2);
	volObj->DefineSection(3, (isPosExt?1:-1)*(length1+length2)	, radius + (isHeadIn?1:0)*thick1	, radius + (isHeadIn?1:0)*thick1 + thick2);

	TGeoVolume *volMed = new TGeoVolume(volName, volObj, med);
	volMed->SetFillColor(color);
	volMed->SetLineColor(color);

	TGeoVolume *volAsm = new TGeoVolumeAssembly(volName);
	volAsm->AddNode(volMed, 1, 0);
	return volAsm;

}

/*				
//________________________________________________________________________
vector<TString> AliITSUv2Layer::CreateDBOuterHoles(
	TString xsName, 
	Double_t xfRadius,
	Double_t xfHoleRadius, 
	Double_t xfLength,
	Double_t xfEachRotate,
	Double_t xfOffsetZ
	)
{ 					
	vector<TString> retStrings;

	Int_t NDegPerStave = 360/fNStaves;
	Double_t Deg2RadFactor = TMath::Pi()/180;
	Double_t HalfWidthBtwnHoles = (24.0 *fgkmm)/2;
	
	for (Int_t i=0; i<fNStaves; i++)
	{
		Double_t CosDegIter = TMath::Cos(i*NDegPerStave*Deg2RadFactor);
		Double_t SinDegIter = TMath::Sin(i*NDegPerStave*Deg2RadFactor);

		// --------------------------- Left ---------
		TString str_l("xtube_l");
		str_l += xsName;
		str_l += i;
		new TGeoTube(str_l, 0.0, xfHoleRadius, xfLength/2);
		
		TString str_l2("xcombt_l");
		str_l2 += xsName;
		str_l2 += i;
		TGeoCombiTrans *xcombt_l = new TGeoCombiTrans(str_l2
			, 0 + (xfRadius * SinDegIter) + (HalfWidthBtwnHoles * CosDegIter)	// circle *** USE xfRadius ***
			, 0 - (xfRadius * CosDegIter) + (HalfWidthBtwnHoles * SinDegIter)	// circle	*** USE xfRadius ***
			, 0 + xfOffsetZ		// set z to the tail of parent object
			// this have been set the position to circle around tube, also offset to the top face(Y-POS) *** BUT NOT X-POS ***
			, new TGeoRotation("", i*NDegPerStave + xfEachRotate, 90.0, 0.0)); // perpendicular to z axis
		
		xcombt_l->RegisterYourself();
		
		// ------------------------- Right -------------
		TString str_r("xtube_r");
		str_r += xsName;
		str_r += i;
		new TGeoTube(str_r, 0.0, xfHoleRadius, xfLength/2);
		
		TString str_r2("xcombt_r");
		str_r2 += xsName;
		str_r2 += i;
		TGeoCombiTrans *xcombt_r = new TGeoCombiTrans(str_r2
			, 0 + (xfRadius * SinDegIter) - (HalfWidthBtwnHoles * CosDegIter)	// circle *** USE xfRadius ***
			, 0 - (xfRadius * CosDegIter) - (HalfWidthBtwnHoles * SinDegIter)	// circle	*** USE xfRadius ***
			, 0 + xfOffsetZ		// set z to the tail of parent object
			// this have been set the position to circle around tube, also offset to the top face(Y-POS) *** BUT NOT X-POS ***
			, new TGeoRotation("", i*NDegPerStave + xfEachRotate, 90.0, 0.0)); // perpendicular to z axis
		
		xcombt_r->RegisterYourself();
		
		// ------------------------- Result ------------
		TString str3("(");
		str3 += str_l;
		str3 += ":";
		str3 += str_l2;
		str3 += " + ";
		str3 += str_r;
		str3 += ":";
		str3 += str_r2;
		str3 += ")";
		
		retStrings.push_back(str3);		
	}
	
	return retStrings;
}	
*/

/*
//________________________________________________________________________
vector<TString> AliITSUv2Layer::CreateDBOuterHandles(
	TString xsName, 
	Double_t xfRadius, 
	Double_t xfLength,
	Double_t xfWidth,
	Double_t xfThick,
	Double_t xfEachRotate,
	Double_t xfOffsetZ
	)
{ 					
	vector<TString> retStrings;

	Int_t NDegPerStave = 360/fNStaves;
	Double_t OffNegSideX = 0.0;
	Double_t Deg2RadFactor = TMath::Pi()/180;
	Double_t CosEachRotate = TMath::Cos(xfEachRotate*Deg2RadFactor);
	Double_t SinEachRotate = TMath::Sin(xfEachRotate*Deg2RadFactor);
	Double_t FaceCirRadius = xfThick/2 + xfRadius;
	
	for (Int_t i=0; i<fNStaves; i++)
	{
		Double_t CosDegIter = TMath::Cos(i*NDegPerStave*Deg2RadFactor);
		Double_t SinDegIter = TMath::Sin(i*NDegPerStave*Deg2RadFactor);					

		TString str("xbbox");
		str += xsName;
		str += i;
		new TGeoBBox(str, xfWidth/2, xfThick/2, xfLength/2);
		
		TString str2("xcombt");
		str2 += xsName;
		str2 += i;
		TGeoCombiTrans *xcombt = new TGeoCombiTrans(str2
			, 0 + (FaceCirRadius * SinDegIter)	// circle with top face offset
				+ (OffNegSideX * CosEachRotate) * CosDegIter // x-offset
				- (OffNegSideX * SinEachRotate) * SinDegIter // y-offset
			, 0 - (FaceCirRadius * CosDegIter)	// circle	with top face offset
				+ (OffNegSideX * CosEachRotate) * SinDegIter // x-offset
				+ (OffNegSideX * SinEachRotate) * CosDegIter // y-offset  (*** ??? why positive, not negative ??? ***)
			, 0 + xfOffsetZ		// set z to the tail of parent object
			, new TGeoRotation("", 0.0, 0.0, i*NDegPerStave + xfEachRotate));
		
		xcombt->RegisterYourself();

		TString str3;
		str3 += str;
		str3 += ":";
		str3 += str2;
		
		retStrings.push_back(str3);		
	}
	
	return retStrings;
}	
*/
