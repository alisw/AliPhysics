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
  fStaveModel(AliITSUv2::kIBModelDummy)
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
  fStaveModel(AliITSUv2::kIBModelDummy)
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
  fStaveModel(AliITSUv2::kIBModelDummy)
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
  fStaveModel(AliITSUv2::kIBModelDummy)
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
  fStaveModel(s.fStaveModel)
{
  //
  // Copy constructor
  for (int i=kNHLevels;i--;) fHierarchy[i] = s.fHierarchy[i];
  //
}

//________________________________________________________________________
AliITSUv2Layer& AliITSUv2Layer::operator=(const AliITSUv2Layer &s)
{
  //
  // Assignment operator 
  //
  if(&s == this) return *this;

  fLayerNumber = s.fLayerNumber;
  fPhi0        = s.fPhi0;
  fLayRadius   = s.fLayRadius;
  fZLength     = s.fZLength;
  fSensorThick = s.fSensorThick;
  fChipThick   = s.fChipThick;
  fStaveWidth  = s.fStaveWidth;
  fStaveTilt   = s.fStaveTilt;
  fNStaves     = s.fNStaves;
  fNModules    = s.fNModules;
  fNChips      = s.fNChips;
  fIsTurbo     = s.fIsTurbo;
  fChipTypeID  = s.fChipTypeID;
  fBuildLevel  = s.fBuildLevel;
  fStaveModel  = s.fStaveModel;
  for (int i=kNHLevels;i--;) fHierarchy[i] = s.fHierarchy[i];
  //
  return *this;
}

//________________________________________________________________________
AliITSUv2Layer::~AliITSUv2Layer() {
  //
  // Destructor
  //
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


  TGeoTube *connTubeHole2 = new TGeoTube("tube2HoleA", 0,
					 fgkIBConnTubeHole2D/2,
					 connBody->GetDZ());

  xpos = fgkIBConnTubesXDist/2;
  ypos = -connBody->GetDY() + fgkIBConnTubesYPos;
  TGeoTranslation *connTubes2Trans1 = new TGeoTranslation("tubes2Trans1A",
							  -xpos, ypos, 0);
  connTubes2Trans1->RegisterYourself();
  TGeoTranslation *connTubes2Trans2 = new TGeoTranslation("tubes2Trans2A",
							   xpos, ypos, 0);
  connTubes2Trans2->RegisterYourself();


  zlen = fgkIBConnTubeHole1ZLen - fgkIBConnTailZLen;
  TGeoTube *connTubeHole3 = new TGeoTube("tube3HoleA", 0,
					 fgkIBConnTubeHole1D/2,
					 zlen);

  zpos = connBody->GetDZ();
  TGeoTranslation *connTubes3Trans1 = new TGeoTranslation("tubes3Trans1A",
							  -xpos, ypos,-zpos);
  connTubes3Trans1->RegisterYourself();
  TGeoTranslation *connTubes3Trans2 = new TGeoTranslation("tubes3Trans2A",
							   xpos, ypos,-zpos);
  connTubes3Trans2->RegisterYourself();


  zlen = fgkIBConnectAFitZLen - fgkIBConnectAFitZOut;
  TGeoTube *connFitHole = new TGeoTube("fitHoleA", 0, fgkIBConnectAFitExtD/2,
				       zlen);

  TGeoTranslation *connFitHoleTrans1 = new TGeoTranslation("fitTrans1A",
							  -xpos, ypos, zpos);
  connFitHoleTrans1->RegisterYourself();
  TGeoTranslation *connFitHoleTrans2 = new TGeoTranslation("fitTrans2A",
							   xpos, ypos, zpos);
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

  TGeoBBox *flex1_5cm  = new TGeoBBox("Flex1MV_5cm",xHalfSt,yFlex1/2,flexOverlap/2);
  TGeoBBox *flex2_5cm  = new TGeoBBox("Flex2MV_5cm",xHalfSt,yFlex2/2,flexOverlap/2);

  // The half stave container (an XTru to avoid overlaps between neightbours)
  yHalfSt = ymod + busAl->GetDY() + busKap->GetDY() + coldPlate->GetDY()
	  + fleeccent->GetDY() + graphlat->GetDY() + fleeclat->GetDY();
  if (fStaveModel == AliITSUv2::kOBModel2)
    yHalfSt += 2*glue->GetDY();

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
