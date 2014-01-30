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


/* $Id: AliITSUv1Layer.cxx  */
// General Root includes
#include <TMath.h>
// Root Geometry includes
//#include <AliLog.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoPcon.h>
#include <TGeoCone.h>
#include <TGeoTube.h> // contaings TGeoTubeSeg
#include <TGeoArb8.h>
#include <TGeoXtru.h>
#include <TGeoCompositeShape.h>
#include <TGeoMatrix.h>
#include "AliITSUv1Layer.h"
#include "AliITSUGeomTGeo.h"
#include <TGeoBBox.h>
#include <TGeoShape.h>
#include <TGeoTrd1.h>
using namespace TMath;

const Double_t AliITSUv1Layer::fgkDefaultSensorThick = 300*fgkmicron;
const Double_t AliITSUv1Layer::fgkDefaultLadderThick =   1*fgkcm;

const Double_t AliITSUv1Layer::fgkOBHalfStaveWidth   =   3.01 *fgkcm;
const Double_t AliITSUv1Layer::fgkOBModuleGap        =   0.01 *fgkcm;
const Double_t AliITSUv1Layer::fgkOBFlexCable1Thick  =   0.005*fgkcm;
const Double_t AliITSUv1Layer::fgkOBFlexCable2Thick  =   0.01 *fgkcm;
const Double_t AliITSUv1Layer::fgkOBBusCable1Thick   =   0.02 *fgkcm;
const Double_t AliITSUv1Layer::fgkOBBusCable2Thick   =   0.02 *fgkcm;
const Double_t AliITSUv1Layer::fgkOBColdPlateThick   =   0.012*fgkcm;
const Double_t AliITSUv1Layer::fgkOBCarbonPlateThick =   0.012*fgkcm;
const Double_t AliITSUv1Layer::fgkOBGlueThick        =   0.03 *fgkcm;
const Double_t AliITSUv1Layer::fgkOBModuleZLength    =  21.06 *fgkcm;


ClassImp(AliITSUv1Layer)

#define SQ(A) (A)*(A)

//________________________________________________________________________
AliITSUv1Layer::AliITSUv1Layer(): 
  AliITSv11Geometry(),
  fLayerNumber(0),
  fPhi0(0),
  fLayRadius(0),
  fZLength(0),
  fSensorThick(0),
  fLadderThick(0),
  fLadderWidth(0),
  fLadderTilt(0),
  fNLadders(0),
  fNModules(0),
  fDetTypeID(0),
  fIsTurbo(0),
  fBuildLevel(0),
  fStaveModel(AliITSUv1::kIBModelDummy)
{
  //
  // Standard constructor
  //
}

//________________________________________________________________________
AliITSUv1Layer::AliITSUv1Layer(Int_t debug): 
  AliITSv11Geometry(debug),
  fLayerNumber(0),
  fPhi0(0),
  fLayRadius(0),
  fZLength(0),
  fSensorThick(0),
  fLadderThick(0),
  fLadderWidth(0),
  fLadderTilt(0),
  fNLadders(0),
  fNModules(0),
  fDetTypeID(0),
  fIsTurbo(0),
  fBuildLevel(0),
  fStaveModel(AliITSUv1::kIBModelDummy)
{
  //
  // Constructor setting debugging level
  //
}

//________________________________________________________________________
AliITSUv1Layer::AliITSUv1Layer(Int_t lay, Int_t debug): 
  AliITSv11Geometry(debug),
  fLayerNumber(lay),
  fPhi0(0),
  fLayRadius(0),
  fZLength(0),
  fSensorThick(0),
  fLadderThick(0),
  fLadderWidth(0),
  fLadderTilt(0),
  fNLadders(0),
  fNModules(0),
  fDetTypeID(0),
  fIsTurbo(0),
  fBuildLevel(0),
  fStaveModel(AliITSUv1::kIBModelDummy)
{
  //
  // Constructor setting layer number and debugging level
  //
}

//________________________________________________________________________
AliITSUv1Layer::AliITSUv1Layer(Int_t lay, Bool_t turbo, Int_t debug): 
  AliITSv11Geometry(debug),
  fLayerNumber(lay),
  fPhi0(0),
  fLayRadius(0),
  fZLength(0),
  fSensorThick(0),
  fLadderThick(0),
  fLadderWidth(0),
  fLadderTilt(0),
  fNLadders(0),
  fNModules(0),
  fDetTypeID(0),
  fIsTurbo(turbo),
  fBuildLevel(0),
  fStaveModel(AliITSUv1::kIBModelDummy)
{
  //
  // Constructor setting layer number and debugging level
  // for a "turbo" layer (i.e. where ladders overlap in phi)
  //
}

//________________________________________________________________________
AliITSUv1Layer::AliITSUv1Layer(const AliITSUv1Layer &s):
  AliITSv11Geometry(s.GetDebug()),
  fLayerNumber(s.fLayerNumber),
  fPhi0(s.fPhi0),
  fLayRadius(s.fLayRadius),
  fZLength(s.fZLength),
  fSensorThick(s.fSensorThick),
  fLadderThick(s.fLadderThick),
  fLadderWidth(s.fLadderWidth),
  fLadderTilt(s.fLadderTilt),
  fNLadders(s.fNLadders),
  fNModules(s.fNModules),
  fDetTypeID(s.fDetTypeID),
  fIsTurbo(s.fIsTurbo),
  fBuildLevel(s.fBuildLevel),
  fStaveModel(s.fStaveModel)
{
  //
  // Copy constructor
  //
}

//________________________________________________________________________
AliITSUv1Layer& AliITSUv1Layer::operator=(const AliITSUv1Layer &s)
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
  fLadderThick = s.fLadderThick;
  fLadderWidth = s.fLadderWidth;
  fLadderTilt  = s.fLadderTilt;
  fNLadders    = s.fNLadders;
  fNModules    = s.fNModules;
  fIsTurbo     = s.fIsTurbo;
  fDetTypeID   = s.fDetTypeID;
  fBuildLevel  = s.fBuildLevel;
  fStaveModel  = s.fStaveModel;

  return *this;
}

//________________________________________________________________________
AliITSUv1Layer::~AliITSUv1Layer() {
  //
  // Destructor
  //
}

//________________________________________________________________________
void AliITSUv1Layer::CreateLayer(TGeoVolume *moth){
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
  if (fLayRadius <= 0) AliFatal(Form("Wrong layer radius (%f)",fLayRadius));
  if (fZLength   <= 0) AliFatal(Form("Wrong layer length (%f)",fZLength));
  if (fNLadders  <= 0) AliFatal(Form("Wrong number of ladders (%d)",fNLadders));
  if (fNModules  <= 0) AliFatal(Form("Wrong number of modules (%d)",fNModules));

  if (fLadderThick <= 0) {
    AliInfo(Form("Ladder thickness wrong or not set (%f), using default (%f)",
		 fLadderThick,fgkDefaultLadderThick));
    fLadderThick = fgkDefaultLadderThick;
  }

  if (fSensorThick <= 0) {
    AliInfo(Form("Sensor thickness wrong or not set (%f), using default (%f)",
		 fSensorThick,fgkDefaultSensorThick));
    fSensorThick = fgkDefaultSensorThick;
  }

  if (fSensorThick > fLadderThick) {
    AliWarning(Form("Sensor thickness (%f) is greater than ladder thickness (%f), fixing",
		 fSensorThick,fLadderThick));
    fSensorThick = fLadderThick;
  }


  // If a Turbo layer is requested, do it and exit
  if (fIsTurbo) {
    CreateLayerTurbo(moth);
    return;
  }


  // First create the ladder container
  alpha = (360./(2*fNLadders))*DegToRad();

  //  fLadderWidth = fLayRadius*Tan(alpha);

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSLayerPattern(),fLayerNumber);
  TGeoVolume *layVol = new TGeoVolumeAssembly(volname);
  layVol->SetUniqueID(fDetTypeID);

//  layVol->SetVisibility(kFALSE);
  layVol->SetVisibility(kTRUE);
  layVol->SetLineColor(1);

  TGeoVolume *laddVol = CreateLadder();


  // Now build up the layer
  alpha = 360./fNLadders;
  Double_t r = fLayRadius + ((TGeoBBox*)laddVol->GetShape())->GetDY();
  for (Int_t j=0; j<fNLadders; j++) {
    Double_t phi = j*alpha + fPhi0;
    xpos = r*CosD(phi);// r*SinD(-phi);
    ypos = r*SinD(phi);// r*CosD(-phi);
    zpos = 0.;
    phi += 90;
    layVol->AddNode(laddVol, j, new TGeoCombiTrans( xpos, ypos, zpos,
						    new TGeoRotation("",phi,0,0)));
  }


  // Finally put everything in the mother volume
  moth->AddNode(layVol, 1, 0);


  // Upgrade geometry is served
  return;
}

//________________________________________________________________________
void AliITSUv1Layer::CreateLayerTurbo(TGeoVolume *moth){
//
// Creates the actual Layer and places inside its mother volume
// A so-called "turbo" layer is a layer where ladders overlap in phi
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
  if (fLadderWidth <= 0)
    AliFatal(Form("Wrong ladder width (%f)",fLadderWidth));
  if (Abs(fLadderTilt) > 45)
    AliWarning(Form("Ladder tilt angle (%f) greater than 45deg",fLadderTilt));


  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSLayerPattern(), fLayerNumber);
  TGeoVolume *layVol = new TGeoVolumeAssembly(volname);
  layVol->SetUniqueID(fDetTypeID);
  layVol->SetVisibility(kTRUE);
  layVol->SetLineColor(1);
  TGeoVolume *laddVol = CreateLadder();


  // Now build up the layer
  alpha = 360./fNLadders;
  Double_t r = fLayRadius /* +module thick ?! */;
  for (Int_t j=0; j<fNLadders; j++) {
    Double_t phi = j*alpha + fPhi0;
    xpos = r*CosD(phi);// r*SinD(-phi);
    ypos = r*SinD(phi);// r*CosD(-phi);
    zpos = 0.;
    phi += 90;
    layVol->AddNode(laddVol, j, new TGeoCombiTrans( xpos, ypos, zpos,
						    new TGeoRotation("", phi-fLadderTilt,0,0)));
  }


  // Finally put everything in the mother volume
  moth->AddNode(layVol, 1, 0);

  return;
}

//________________________________________________________________________
TGeoVolume* AliITSUv1Layer::CreateLadder(const TGeoManager * /*mgr*/){
//
// Creates the actual Ladder
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
//

  char volname[30];
  Double_t xLenO = 5.79/2;
 
  Double_t xlenI, ylen, zlen;
  Double_t xpos, ypos, zpos, zmod;
  Double_t alpha;


  // First create all needed shapes
  alpha = (360./(2*fNLadders))*DegToRad();

  // The ladder
  xlenI = fLayRadius*Tan(alpha);
  if (fIsTurbo) xlenI = 0.5*fLadderWidth;
  ylen = 0.5*fLadderThick;
  zlen = 0.5*fZLength;

  Double_t yplus = 0.46;
  TGeoXtru *ladder = new TGeoXtru(2); //z sections
  Double_t xv[5] = {xlenI,xlenI,0,-xlenI,-xlenI};
  Double_t yv[5] = {ylen+0.09,-0.15,-yplus-fSensorThick,-0.15,ylen+0.09};    
  ladder->DefinePolygon(5,xv,yv);
  ladder->DefineSection(0,-zlen,0,0,1.);
  ladder->DefineSection(1,+zlen,0,0,1.);

  // We have all shapes: now create the real volumes

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSLadderPattern(), fLayerNumber);
//  TGeoVolume *laddVol = new TGeoVolume(volname, ladder, medAir);
  TGeoVolume *laddVol = new TGeoVolumeAssembly(volname);

  //  laddVol->SetVisibility(kFALSE);
  laddVol->SetVisibility(kTRUE);
  laddVol->SetLineColor(2);
  TGeoVolume *modVol = 0;
  TGeoVolume *mechLaddVol = 0;

  // Now build up the ladder
  if (fLayerNumber<3) {
    modVol = CreateModuleInnerB(xlenI,ylen,zlen);
    zmod = ((TGeoBBox*)modVol->GetShape())->GetDZ();
    for (Int_t j=0; j<fNModules; j++) {
      xpos = 0.;
      ypos = 0.021;  // Remove small overlap - M.S: 21may13
      zpos = -ladder->GetDZ() + j*2*zmod + zmod;
      laddVol->AddNode(modVol, j, new TGeoTranslation(xpos, ypos, zpos));
    }
 
  // put mechanical stave structure, only inner barrel up to now
    mechLaddVol = CreateStaveStructInnerB(xlenI,zlen); 
    if (mechLaddVol)
      laddVol->AddNode(mechLaddVol, fNModules, new TGeoCombiTrans(0, -0.15-ylen, 0, new TGeoRotation("",0, 0, 180)));
  }

  else{
    if (fStaveModel == AliITSUv1::kOBModel0) { // Create simplified stave struct as in v0
      modVol = CreateModuleInnerB(xlenI,ylen,zlen);
  printf("?????? %f %f %f\n",xlenI,ylen,zlen);
      zmod = ((TGeoBBox*)modVol->GetShape())->GetDZ();
      for (Int_t j=0; j<fNModules; j++) {
	xpos = 0.;
	ypos = 0.021;  // Remove small overlap - M.S: 21may13
	zpos = -ladder->GetDZ() + j*2*zmod + zmod;
	laddVol->AddNode(modVol, j, new TGeoTranslation(xpos, ypos, zpos));
      }
    } else { // (if fStaveModel) Create new stave struct as in TDR
      modVol = CreateStaveOuterB(xLenO);
      laddVol->AddNode(modVol, 1, new TGeoTranslation(0, 2.5, 0));

      mechLaddVol = CreateSpaceFrameOuterB(xLenO); 
      if (mechLaddVol)
	laddVol->AddNode(mechLaddVol, 1,
			 new TGeoCombiTrans(0, 0, 0,
					    new TGeoRotation("", 180, 0, 0)));
    } // if (fStaveModel)
  }
  

  // Done, return the ladder
  return laddVol;
}

//________________________________________________________________________
TGeoVolume* AliITSUv1Layer::CreateStaveStructInnerB(const Double_t xlad,
						    const Double_t zlad,
						    const TGeoManager *mgr){
//
// Create the mechanical stave structure
//
// Input:
//         xlad : X length
//         zlad : Z length
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      22 Mar 2013  Chinorat Kobdaj
// Updated:      26 Apr 2013  Mario Sitta
//

  TGeoVolume *mechLaddVol = 0;

  switch (fStaveModel) {
    case AliITSUv1::kIBModelDummy:
      mechLaddVol = CreateStaveModelInnerBDummy(xlad,zlad,mgr);
      break;
    case AliITSUv1::kIBModel0:
      mechLaddVol = CreateStaveModelInnerB0(xlad,zlad,mgr);
      break;
    case AliITSUv1::kIBModel1:
      mechLaddVol = CreateStaveModelInnerB1(xlad,zlad,mgr);
      break;
    case AliITSUv1::kIBModel21:
      mechLaddVol = CreateStaveModelInnerB21(xlad,zlad,mgr);
      break;
    case AliITSUv1::kIBModel22:
      mechLaddVol = CreateStaveModelInnerB22(xlad,zlad,mgr);
      break;
    case AliITSUv1::kIBModel3:
      mechLaddVol = CreateStaveModelInnerB3(xlad,zlad,mgr);
      break;
    default:
      AliFatal(Form("Unknown stave model %d",fStaveModel));
      break;
  }

  return mechLaddVol; 
}


//________________________________________________________________________
TGeoVolume* AliITSUv1Layer::CreateStaveModelInnerBDummy(const Double_t ,
							const Double_t ,
							const TGeoManager *) const {
//
// Create dummy stave
//
// Input:
//         xlad : X length
//         zlad : Z length
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
TGeoVolume* AliITSUv1Layer::CreateStaveModelInnerB0(const Double_t xlad,
						    const Double_t zlad,
						    const TGeoManager *mgr){
//
// Create the mechanical stave structure for Model 0 of TDR
//
// Input:
//         xlad : X length
//         zlad : Z length
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      22 Mar 2013  Chinorat Kobdaj
// Updated:      26 Apr 2013  Mario Sitta
//
  
  // Materials defined in AliITSUv1
  TGeoMedium *medAir    = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medWater  = mgr->GetMedium("ITS_WATER$");

  TGeoMedium *medM60J3K    = mgr->GetMedium("ITS_M60J3K$"); 
  TGeoMedium *medKapton    = mgr->GetMedium("ITS_KAPTON(POLYCH2)$");
  TGeoMedium *medGlue      = mgr->GetMedium("ITS_GLUE$");
  TGeoMedium *medFlexCable = mgr->GetMedium("ITS_FLEXCABLE$");

  // Local parameters
  Double_t kConeOutRadius = 0.15/2;
  Double_t kConeInRadius = 0.1430/2;
  Double_t kStaveLength = zlad*2;
  Double_t kStaveWidth = xlad*2-kConeOutRadius*2;
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
  snprintf(volname, 30, "%s%d_StaveStruct", AliITSUGeomTGeo::GetITSLadderPattern(), fLayerNumber);

  Double_t z=0, y=-0.011+0.0150, x=0;

   TGeoVolume *mechLaddVol = 0;

  if (fBuildLevel < 5) {

    // world (trapezoid)
    TGeoXtru *mechStruct = new TGeoXtru(2); //z sections
    Double_t xv[5] = {kStaveWidth/2+0.1,kStaveWidth/2+0.1,0,-kStaveWidth/2-0.1,-kStaveWidth/2-0.1};
    Double_t yv[5] = {-kConeOutRadius*2-0.07,0,kStaveHeight,0,-kConeOutRadius*2-0.07};    
    mechStruct->DefinePolygon(5,xv,yv);
    mechStruct->DefineSection(0,-kStaveLength-0.1,0,0,1.);
    mechStruct->DefineSection(1,kStaveLength+0.1,0,0,1.);

    mechLaddVol = new TGeoVolume(volname, mechStruct, medAir);
    mechLaddVol->SetLineColor(12);
    mechLaddVol->SetFillColor(12); 
    mechLaddVol->SetVisibility(kTRUE);
      
    // detailed structure ++++++++++++++
    //Pipe Kapton grey-35
    TGeoTube *coolTube = new TGeoTube(kConeInRadius,kConeOutRadius,kStaveLength/2);
    TGeoVolume *volCoolTube= new TGeoVolume("pipe", coolTube, medKapton);
    volCoolTube->SetFillColor(35);
    volCoolTube->SetLineColor(35);
    mechLaddVol->AddNode(volCoolTube,0,new TGeoTranslation(x+(kStaveWidth/2),y-(kHeight-kConeOutRadius),0));
    mechLaddVol->AddNode(volCoolTube,1,new TGeoTranslation(x-(kStaveWidth/2),y-(kHeight-kConeOutRadius),0));
  }

  if (fBuildLevel < 4) {
    TGeoTube *coolTubeW = new TGeoTube(0.,kConeInRadius,kStaveLength/2);
    TGeoVolume *volCoolTubeW= new TGeoVolume("pipeWater", coolTubeW, medWater);
    volCoolTubeW->SetFillColor(4);
    volCoolTubeW->SetLineColor(4);
    mechLaddVol->AddNode(volCoolTubeW,0,new TGeoTranslation(x+(kStaveWidth/2),y-(kHeight-kConeOutRadius),0));
    mechLaddVol->AddNode(volCoolTubeW,1,new TGeoTranslation(x-(kStaveWidth/2),y-(kHeight-kConeOutRadius),0));
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
      mechLaddVol->AddNode(volT2,4*i+0,
				  new TGeoCombiTrans(x+kWidth,y+(2*kConeOutRadius),z-kStaveLength/2+(i*(4/n)*kL1)+kS1/2,//z-14.25+(i*2*kL1),
						     new TGeoRotation("volT2",90,90-kAlpha,90-kBeta)));
      mechLaddVol->AddNode(volT2,4*i+1,
				  new TGeoCombiTrans(x-kWidth,y+(2*kConeOutRadius),z-kStaveLength/2+(i*(4/n)*kL1)+kS1/2,//z-14.25+(i*2*kL1),
						     new TGeoRotation("volT2",90,-90+kAlpha,-90+kBeta)));
      mechLaddVol->AddNode(volT2,4*i+2,
				  new TGeoCombiTrans(x+kWidth,y+(2*kConeOutRadius),z-kStaveLength/2+(i*(4/n)*kL1)+kS1/2,//z-14.25+(i*2*kL1),
						     new TGeoRotation("volT2",90,-90+kAlpha,90-kBeta)));
      mechLaddVol->AddNode(volT2,4*i+3,
				  new TGeoCombiTrans(x-kWidth,y+(2*kConeOutRadius),z-kStaveLength/2+(i*(4/n)*kL1)+kS1/2,//z-14.25+(i*2*kL1),  
						     new TGeoRotation("volT2",90,90-kAlpha,-90+kBeta)));
    }


    //Bottom CFRP Filament black-12 Carbon structure  TGeoBBox (thickness,width,length)
    TGeoBBox *t1=new TGeoBBox(0.007/2,0.15/2,kS1);//(0.002,0.02,kS1);
    TGeoVolume *volT1=new TGeoVolume("CFRPBottom", t1, medM60J3K);
    volT1->SetLineColor(12);
    volT1->SetFillColor(12); 

    for(int i=1;i<loop;i++){
      mechLaddVol->AddNode(volT1,4*i+0,
				  new TGeoCombiTrans(x+kWidth,y-kHeight,z-kStaveLength/2+((4/n)*kL1*i)+kS1/2, //z-14.25+(i*2*kL1),  
						     new TGeoRotation("volT1",-90,kAlpha,0)));
      mechLaddVol->AddNode(volT1,4*i+1,
				  new TGeoCombiTrans(x-kWidth,y-kHeight,z-kStaveLength/2+((4/n)*kL1*i)+kS1/2,  //z-14.25+(i*2*kL1), 
						     new TGeoRotation("volT1",90,kAlpha,0)));
      mechLaddVol->AddNode(volT1,4*i+2,
				  new TGeoCombiTrans(x+kWidth,y-kHeight,z-kStaveLength/2+(i*(4/n)*kL1)+kS1/2, //z-14.25+(i*2*kL1), 
						     new TGeoRotation("volT1",-90,-kAlpha,0)));
      mechLaddVol->AddNode(volT1,4*i+3,
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
      mechLaddVol->AddNode(volTG,4*i+0,
				  new TGeoCombiTrans(x+kWidth,y-0.16,z-kStaveLength/2+((4/n)*kL1*i)+kS1/2, //z-14.25+(2*kL1*i), 
						     new TGeoRotation("volTG",-90,kAlpha,0)));
      mechLaddVol->AddNode(volTG,4*i+1,
				  new TGeoCombiTrans(x-kWidth,y-0.16,z-kStaveLength/2+((4/n)*kL1*i)+kS1/2, //z-14.25+(2*kL1*i), 
						     new TGeoRotation("volTG",90,kAlpha,0)));
      mechLaddVol->AddNode(volTG,4*i+2,
				  new TGeoCombiTrans(x+kWidth,y-0.16,z-kStaveLength/2+((4/n)*i*kL1)+kS1/2, //z-14.25+(i*2*kL1), 
						     new TGeoRotation("volTG",-90,-kAlpha,0)));
      mechLaddVol->AddNode(volTG,4*i+3,
				  new TGeoCombiTrans(x-kWidth,y-0.16,z-kStaveLength/2+(i*(4/n)*kL1)+kS1/2, //z-14.25+(i*2*kL1), 
						     new TGeoRotation("volTG",-90,+kAlpha,0)));
    }

    TGeoBBox *glue = new TGeoBBox(xlad, 0.005/2, zlad);
    TGeoVolume *volGlue=new TGeoVolume("Glue2", glue, medGlue);
    volGlue->SetLineColor(5);
    volGlue->SetFillColor(5); 
    //mechLaddVol->AddNode(volGlue, 0, new TGeoCombiTrans(x, y-0.16, z, new TGeoRotation("",0, 0, 0)));
    mechLaddVol->AddNode(volGlue, 1, new TGeoCombiTrans(x, y-0.165-fSensorThick-0.005, z, new TGeoRotation("",0, 0, 0)));
  }

  if (fBuildLevel < 1) {
    //Flex cable brown-28 TGeoBBox(width,thickness,length); 
    TGeoBBox *kapCable = new TGeoBBox(xlad, 0.01/2, zlad);
    TGeoVolume *volCable=new TGeoVolume("FlexCable", kapCable, medFlexCable);
    volCable->SetLineColor(28);
    volCable->SetFillColor(28); 
    mechLaddVol->AddNode(volCable, 0, new TGeoCombiTrans(x, y-0.165-fSensorThick-0.005-0.01, z, new TGeoRotation("",0, 0, 0)));
 }

  // Done, return the stave structur
  return mechLaddVol;
}


//________________________________________________________________________
TGeoVolume* AliITSUv1Layer::CreateStaveModelInnerB1(const Double_t xlad,
						    const Double_t zlad,
						    const TGeoManager *mgr){
//
// Create the mechanical stave structure for Model 1 of TDR
//
// Input:
//         xlad : X length
//         zlad : Z length
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      22 Mar 2013  Chinorat Kobdaj
// Updated:      26 Apr 2013  Mario Sitta
//
  
  // Materials defined in AliITSUv1
  TGeoMedium *medAir    = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medWater  = mgr->GetMedium("ITS_WATER$");

  TGeoMedium *medM60J3K    = mgr->GetMedium("ITS_M60J3K$"); 
  TGeoMedium *medKapton    = mgr->GetMedium("ITS_KAPTON(POLYCH2)$");
  TGeoMedium *medGlue      = mgr->GetMedium("ITS_GLUE$");
  TGeoMedium *medFlexCable = mgr->GetMedium("ITS_FLEXCABLE$");

  // Local parameters
  Double_t kConeOutRadius = 0.15/2;
  //    Double_t kConeInRadius = 0.1430/2;
  Double_t kStaveLength = zlad*2;
  //    Double_t kStaveWidth = xlad*2-kConeOutRadius*2;
  Double_t kStaveWidth = xlad*2;
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
  

  TGeoVolume *mechLaddVol = 0;

  char volname[30];
  snprintf(volname, 30, "%s%d_StaveStruct", AliITSUGeomTGeo::GetITSLadderPattern(), fLayerNumber);
    

  // detailed structure ++++++++++++++
  Double_t z=0, y=-0.011+0.0150, x=0;

  // Polimide micro channels numbers
  Double_t yMC = y-kHeight+0.01;
  Int_t nb = (Int_t)(kStaveWidth/0.1)+1;
  Double_t xladMC = (nb*0.1-0.08)/2;


  if (fBuildLevel < 5) {
    // world (trapezoid)
    TGeoXtru *mechStruct = new TGeoXtru(2); //z sections
    Double_t xv[5] = {kStaveWidth/2+0.1,kStaveWidth/2+0.1,0,-kStaveWidth/2-0.1,-kStaveWidth/2-0.1};
    Double_t yv[5] = {-kConeOutRadius*2-0.07,0,kStaveHeight,0,-kConeOutRadius*2-0.07};    
    mechStruct->DefinePolygon(5,xv,yv);
    mechStruct->DefineSection(0,-kStaveLength-0.1,0,0,1.);
    mechStruct->DefineSection(1,kStaveLength+0.1,0,0,1.);

    mechLaddVol = new TGeoVolume(volname, mechStruct, medAir);
    mechLaddVol->SetLineColor(12);
    mechLaddVol->SetFillColor(12); 
    mechLaddVol->SetVisibility(kTRUE);
      
    // Polimide micro channels numbers
    TGeoBBox *tM0=new TGeoBBox(xladMC, 0.005/2, zlad);
    TGeoVolume *volTM0=new TGeoVolume("MicroChanCover", tM0, medKapton);
    volTM0->SetLineColor(35);
    volTM0->SetFillColor(35); 
    mechLaddVol->AddNode(volTM0, 0, new TGeoCombiTrans(x,-0.0125+yMC, z, new TGeoRotation("",0, 0, 0)));
    mechLaddVol->AddNode(volTM0, 1, new TGeoCombiTrans(x,+0.0125+yMC, z, new TGeoRotation("",0, 0, 0)));
      
    TGeoBBox *tM0b=new TGeoBBox(0.02/2, 0.02/2, zlad);
    TGeoVolume *volTM0b=new TGeoVolume("MicroChanWalls", tM0b, medKapton);
    volTM0b->SetLineColor(35);
    volTM0b->SetFillColor(35); 
    for (Int_t ib=0;ib<nb;ib++) {
      mechLaddVol->AddNode(volTM0b, ib, new TGeoCombiTrans(x+ib*0.1-xladMC+0.01,yMC, z, new TGeoRotation("",0, 0, 0)));
    }
      
  }
    
  if (fBuildLevel < 4) {
    // Water in Polimide micro channels
    TGeoBBox *water=new TGeoBBox(0.08/2, 0.02/2, zlad+0.1);
    TGeoVolume *volWater=new TGeoVolume("Water", water, medWater);
    volWater->SetLineColor(4);
    volWater->SetFillColor(4); 
    for (Int_t ib=0;ib<(nb-1);ib++) {
      mechLaddVol->AddNode(volWater, ib, new TGeoCombiTrans(x+ib*0.1-xladMC+0.06,yMC, z, new TGeoRotation("",0, 0, 0)));
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
      mechLaddVol->AddNode(volT1,4*i+0,
				  new TGeoCombiTrans(x+kWidth,y-kHeight+0.04+filHeight/2,z-kStaveLength/2+(4*kL1)+kS1/2, 
						     new TGeoRotation("volT1",-90,kAlpha,0)));
      mechLaddVol->AddNode(volT1,4*i+1,
				  new TGeoCombiTrans(x-kWidth,y-kHeight+0.04+filHeight/2,z-kStaveLength/2+(4*kL1*i)+kS1/2, 
						     new TGeoRotation("volT1",90,kAlpha,0)));
      mechLaddVol->AddNode(volT1,4*i+2,
				  new TGeoCombiTrans(x+kWidth,y-kHeight+0.04+filHeight/2,z-kStaveLength/2+2*kL1+(i*4*kL1)+kS1/2, 
						     new TGeoRotation("volT1",-90,-kAlpha,0)));
      mechLaddVol->AddNode(volT1,4*i+3,
				  new TGeoCombiTrans(x-kWidth,y-kHeight+0.04+filHeight/2,z-kStaveLength/2+2*kL1+(i*4*kL1)+kS1/2, 
						     new TGeoRotation("volT1",-90,+kAlpha,0)));
    }

      // Top filament CFRP black-12 Carbon structure TGeoBBox (length,thickness,width)
    TGeoBBox *t2=new TGeoBBox(kS2,filHeight/2,filWidth/2);
    TGeoVolume *volT2=new TGeoVolume("CFRPTop", t2, medM60J3K);
    volT2->SetLineColor(12);
    volT2->SetFillColor(12); 
    for(int i=0;i<loop;i++){ //i<30;i++){
      mechLaddVol->AddNode(volT2,4*i+0,
				  new TGeoCombiTrans(x+kWidth,y+0.04+filHeight/2,z-kStaveLength/2+(i*4*kL1)+kS1/2,
						     new TGeoRotation("volT2",90,90-kAlpha,90-kBeta)));
      mechLaddVol->AddNode(volT2,4*i+1,
				  new TGeoCombiTrans(x-kWidth,y+0.04+filHeight/2,z-kStaveLength/2+(i*4*kL1)+kS1/2,
						     new TGeoRotation("volT2",90,-90+kAlpha,-90+kBeta)));
      mechLaddVol->AddNode(volT2,4*i+2,
				  new TGeoCombiTrans(x+kWidth,y+0.04+filHeight/2,z-kStaveLength/2+2*kL1+(i*4*kL1)+kS1/2,
						     new TGeoRotation("volT2",90,-90+kAlpha,90-kBeta)));
      mechLaddVol->AddNode(volT2,4*i+3,
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
      mechLaddVol->AddNode(volT3,4*i+0,
				  new TGeoCombiTrans(x+kWidth,y-kHeight+0.0325,z-kStaveLength/2+(4*kL1*i)+kS1/2, 
						     new TGeoRotation("volT1",-90,kAlpha,0)));
      mechLaddVol->AddNode(volT3,4*i+1,
				  new TGeoCombiTrans(x-kWidth,y-kHeight+0.0325,z-kStaveLength/2+(4*kL1*i)+kS1/2, 
						     new TGeoRotation("volT1",90,kAlpha,0)));
      mechLaddVol->AddNode(volT3,4*i+2,
				  new TGeoCombiTrans(x+kWidth,y-kHeight+0.0325,z-kStaveLength/2+2*kL1+(i*4*kL1)+kS1/2, 
						     new TGeoRotation("volT1",-90,-kAlpha,0)));
      mechLaddVol->AddNode(volT3,4*i+3,
				  new TGeoCombiTrans(x-kWidth,y-kHeight+0.0325,z-kStaveLength/2+2*kL1+(i*4*kL1)+kS1/2, 
						     new TGeoRotation("volT1",-90,+kAlpha,0)));
    }
      
    // Glue microchannel and sensor
    TGeoBBox *glueM = new TGeoBBox(xlad, 0.01/2, zlad);
    TGeoVolume *volGlueM=new TGeoVolume("MicroChanGlue", glueM, medGlue);
    volGlueM->SetLineColor(5);
    volGlueM->SetFillColor(5); 
    mechLaddVol->AddNode(volGlueM, 0, new TGeoCombiTrans(x, y-0.16, z, new TGeoRotation("",0, 0, 0)));

    // Glue sensor and kapton
    TGeoBBox *glue = new TGeoBBox(xlad, 0.005/2, zlad);
    TGeoVolume *volGlue=new TGeoVolume("SensorGlue", glue, medGlue);
    volGlue->SetLineColor(5);
    volGlue->SetFillColor(5); 
    mechLaddVol->AddNode(volGlue, 1, new TGeoCombiTrans(x, y-0.165-fSensorThick-0.005, z, new TGeoRotation("",0, 0, 0)));
  }

  if (fBuildLevel < 1) {
    TGeoBBox *kapCable = new TGeoBBox(xlad, 0.01/2, zlad);
    TGeoVolume *volCable=new TGeoVolume("FlexCable", kapCable, medFlexCable);
    volCable->SetLineColor(28);
    volCable->SetFillColor(28); 
    mechLaddVol->AddNode(volCable, 0, new TGeoCombiTrans(x, y-0.165-fSensorThick-0.005-0.01, z, new TGeoRotation("",0, 0, 0)));
  }
    
  // Done, return the stave structur
  return mechLaddVol;

}

//________________________________________________________________________
TGeoVolume* AliITSUv1Layer::CreateStaveModelInnerB21(const Double_t xlad,
						     const Double_t zlad,
						     const TGeoManager *mgr){
//
// Create the mechanical stave structure for Model 2.1 of TDR
//
// Input:
//         xlad : X length
//         zlad : Z length
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      22 Mar 2013  Chinorat Kobdaj
// Updated:      26 Apr 2013  Mario Sitta
//
  
  // Materials defined in AliITSUv1
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
  Double_t kStaveLength = zlad;
  Double_t kStaveWidth = xlad*2;
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
  snprintf(volname, 30, "%s%d_StaveStruct", AliITSUGeomTGeo::GetITSLadderPattern(), fLayerNumber);

  Double_t z=0, y=-(kConeOutRadius+0.03)+0.0385, x=0;

  TGeoVolume *mechLaddVol = 0;

  if (fBuildLevel < 5) {
    // world (trapezoid)
    TGeoXtru *mechStruct = new TGeoXtru(2); //z sections
    Double_t xv[5] = {kStaveWidth/2+0.1,kStaveWidth/2+0.1,0,-kStaveWidth/2-0.1,-kStaveWidth/2-0.1};
    Double_t yv[5] = {-kConeOutRadius*2-0.07,0,kStaveHeigth,0,-kConeOutRadius*2-0.07};    
    mechStruct->DefinePolygon(5,xv,yv);
    mechStruct->DefineSection(0,-kStaveLength-0.1,0,0,1.);
    mechStruct->DefineSection(1,kStaveLength+0.1,0,0,1.);

    mechLaddVol = new TGeoVolume(volname, mechStruct, medAir);
    mechLaddVol->SetLineColor(12);
    mechLaddVol->SetFillColor(12); 
    mechLaddVol->SetVisibility(kTRUE);  
      
    //Pipe Kapton grey-35 
    TGeoCone *cone1 = new TGeoCone(kStaveLength,kConeInRadius,kConeOutRadius,kConeInRadius,kConeOutRadius);
    TGeoVolume *volCone1= new TGeoVolume("PolyimidePipe", cone1, medKapton);
    volCone1->SetFillColor(35);
    volCone1->SetLineColor(35);
    mechLaddVol->AddNode(volCone1,1,new TGeoTranslation(x+0.25,y,z));
    mechLaddVol->AddNode(volCone1,2,new TGeoTranslation(x-0.25,y,z));
  }

  if (fBuildLevel < 4) {
    
    TGeoTube *coolTubeW = new TGeoTube(0.,kConeInRadius,kStaveLength);
    TGeoVolume *volCoolTubeW= new TGeoVolume("Water", coolTubeW, medWater);
    volCoolTubeW->SetFillColor(4);
    volCoolTubeW->SetLineColor(4);
    mechLaddVol->AddNode(volCoolTubeW,0,new TGeoTranslation(x-0.25,y,z));
    mechLaddVol->AddNode(volCoolTubeW,1,new TGeoTranslation(x+0.25,y,z));
  }

  if (fBuildLevel < 3) {
    //top fillament
    // Top filament M60J black-12 Carbon structure TGeoBBox (length,thickness,width)
    TGeoBBox *t2=new TGeoBBox(kS2,0.02/2,0.04/2); //TGeoBBox *t2=new TGeoBBox(kS2,0.01,0.02);
    TGeoVolume *volT2=new TGeoVolume("TopFilament", t2, medM60J3K);
    volT2->SetLineColor(12);
    volT2->SetFillColor(12); 
    for(int i=0;i<loop;i++){// i<28;i++){
      mechLaddVol->AddNode(volT2,i*4+1,new TGeoCombiTrans(x+kWidth,y+kHeight+(0.12/2)-0.014+0.007,z-kStaveLength+(i*4*kL1)+kS1/2, new TGeoRotation("volT2",90,90-kAlpha,90-kBeta)));
      mechLaddVol->AddNode(volT2,i*4+2,new TGeoCombiTrans(x-kWidth,y+kHeight+(0.12/2)-0.014+0.007,z-kStaveLength+(i*4*kL1)+kS1/2, new TGeoRotation("volT2",90,-90+kAlpha,-90+kBeta)));
      mechLaddVol->AddNode(volT2,i*4+3,new TGeoCombiTrans(x+kWidth,y+kHeight+(0.12/2)-0.014+0.007,z-kStaveLength+2*kL1+(i*4*kL1)+kS1/2, new TGeoRotation("volT2",90,-90+kAlpha,90-kBeta)));
      mechLaddVol->AddNode(volT2,i*4+4,new TGeoCombiTrans(x-kWidth,y+kHeight+(0.12/2)-0.014+0.007,z-kStaveLength+2*kL1+(i*4*kL1)+kS1/2, new TGeoRotation("volT2",90,90-kAlpha,-90+kBeta)));
//    mechLaddVol->AddNode(volT2,i*4+1,new TGeoCombiTrans(x+kWidth+0.0036,y+kHeight-(0.12/2)+0.072,z+kStaveLength+(i*4*kL1)+kS1/2, new TGeoRotation("volT2",90,90-kAlpha,90-kBeta)));

 }
 
    //wall side structure out
    TGeoBBox *box4 = new TGeoBBox(0.03/2,0.12/2,kStaveLength-0.50);
    TGeoVolume *plate4 = new TGeoVolume("WallOut",box4,medM60J3K);
    plate4->SetFillColor(35);
    plate4->SetLineColor(35);
    mechLaddVol->AddNode(plate4,1,new TGeoCombiTrans(x+(2*kStaveWidth/4)-(0.03/2),y-0.0022-kConeOutRadius+0.12/2+0.007,z,new TGeoRotation("plate4",0,0,0)));
    mechLaddVol->AddNode(plate4,2,new TGeoCombiTrans(x-(2*kStaveWidth/4)+(0.03/2),y-0.0022-kConeOutRadius+0.12/2+0.007,z,new TGeoRotation("plate4",0,0,0)));
    //wall side in
    TGeoBBox *box5 = new TGeoBBox(0.015/2,0.12/2,kStaveLength-0.50);
    TGeoVolume *plate5 = new TGeoVolume("WallIn",box5,medM60J3K);
    plate5->SetFillColor(12);
    plate5->SetLineColor(12);
    mechLaddVol->AddNode(plate5,1,new TGeoCombiTrans(x+(2*kStaveWidth/4)-0.03-0.015/2,y-0.0022-kConeOutRadius+0.12/2+0.007,z,new TGeoRotation("plate5",0,0,0)));
    mechLaddVol->AddNode(plate5,2,new TGeoCombiTrans(x-(2*kStaveWidth/4)+0.03+0.015/2,y-0.0022-kConeOutRadius+0.12/2+0.007,z,new TGeoRotation("plate5",0,0,0)));

     //Amec Thermasol red-2 cover tube FGS300
    TGeoConeSeg *cons1 = new TGeoConeSeg(kStaveLength-0.50,kConeOutRadius,kConeOutRadius+kLay1,kConeOutRadius,kConeOutRadius+kLay1,0,180);
    TGeoVolume *cone11 = new TGeoVolume("ThermasolPipeCover",cons1,medFGS003);
    cone11->SetFillColor(2);
    cone11->SetLineColor(2);
    mechLaddVol->AddNode(cone11,1,new TGeoCombiTrans(x+0.25,y,z,new TGeoRotation("Cone11",0,0,0)));
    mechLaddVol->AddNode(cone11,2,new TGeoCombiTrans(x-0.25,y,z,new TGeoRotation("Cone11",0,0,0)));

    TGeoBBox *box2 = new TGeoBBox((0.50-(2*kConeOutRadius))/2,kLay1/2,kStaveLength-0.50);
    TGeoVolume *plate2 = new TGeoVolume("ThermasolMiddle",box2,medFGS003);
    plate2->SetFillColor(2);
    plate2->SetLineColor(2);
    mechLaddVol->AddNode(plate2,1,new TGeoCombiTrans(x,y-kConeOutRadius+(kLay1/2),z,new TGeoRotation("plate2",0,0,0)));

    TGeoBBox *box21 = new TGeoBBox((0.75-0.25-kConeOutRadius-kLay1)/2,kLay1/2,kStaveLength-0.50);
    TGeoVolume *plate21 = new TGeoVolume("ThermasolLeftRight",box21,medFGS003);
    plate21->SetFillColor(2);
    plate21->SetLineColor(2);
    mechLaddVol->AddNode(plate21,1,new TGeoCombiTrans(x+0.25+kConeOutRadius+(0.75-0.25-kConeOutRadius)/2-(kLay1/2),y-kConeOutRadius+(kLay1/2),z,new TGeoRotation("plate21",0,0,0)));
    mechLaddVol->AddNode(plate21,2,new TGeoCombiTrans(x-0.25-kConeOutRadius-(0.75-0.25-kConeOutRadius)/2+(kLay1/2),y-kConeOutRadius+(kLay1/2),z,new TGeoRotation("plate21",0,0,0)));

    TGeoBBox *box22 = new TGeoBBox((kLay1/2),kConeOutRadius/2,kStaveLength-0.50);
    TGeoVolume *plate22 = new TGeoVolume("ThermasolVertical",box22,medFGS003);
    plate22->SetFillColor(2);
    plate22->SetLineColor(2);
    mechLaddVol->AddNode(plate22,1,new TGeoCombiTrans(x+0.25+kConeOutRadius+(kLay1/2),y-kConeOutRadius/2,z,new TGeoRotation("plate22",0,0,0)));
    mechLaddVol->AddNode(plate22,2,new TGeoCombiTrans(x+0.25-kConeOutRadius-(kLay1/2),y-kConeOutRadius/2,z,new TGeoRotation("plate22",0,0,0)));
    mechLaddVol->AddNode(plate22,3,new TGeoCombiTrans(x-0.25+kConeOutRadius+(kLay1/2),y-kConeOutRadius/2,z,new TGeoRotation("plate22",0,0,0)));
    mechLaddVol->AddNode(plate22,4,new TGeoCombiTrans(x-0.25-kConeOutRadius-(kLay1/2),y-kConeOutRadius/2,z,new TGeoRotation("plate22",0,0,0)));

    //C Fleece
    TGeoConeSeg *cons2 = new TGeoConeSeg(kStaveLength-0.50,kConeOutRadius+kLay1,kConeOutRadius+kLay1+kLay2,kConeOutRadius+kLay1,kConeOutRadius+kLay1+kLay2,0,180); 
    TGeoVolume *cone12 = new TGeoVolume("CFleecePipeCover",cons2,medCarbonFleece);
    cone12->SetFillColor(28);
    cone12->SetLineColor(28);
    mechLaddVol->AddNode(cone12,1,new TGeoCombiTrans(x+0.25,y,z,new TGeoRotation("Cone12",0,0,0)));
    mechLaddVol->AddNode(cone12,2,new TGeoCombiTrans(x-0.25,y,z,new TGeoRotation("Cone12",0,0,0)));

    TGeoBBox *box3 = new TGeoBBox((0.50-(2*(kConeOutRadius+kLay1)))/2,kLay2/2,kStaveLength-0.50);
    TGeoVolume *plate3 = new TGeoVolume("CFleeceMiddle",box3,medCarbonFleece);
    plate3->SetFillColor(28);
    plate3->SetLineColor(28);
    mechLaddVol->AddNode(plate3,1,new TGeoCombiTrans(x,y-kConeOutRadius+kLay1+(kLay2/2),z,new TGeoRotation("plate3",0,0,0)));

    TGeoBBox *box31 = new TGeoBBox((0.75-0.25-kConeOutRadius-kLay1)/2,kLay2/2,kStaveLength-0.50);
    TGeoVolume *plate31 = new TGeoVolume("CFleeceLeftRight",box31,medCarbonFleece);
    plate31->SetFillColor(28);
    plate31->SetLineColor(28);
    mechLaddVol->AddNode(plate31,1,new TGeoCombiTrans(x+0.25+kConeOutRadius+kLay1+(0.75-0.25-kConeOutRadius-kLay1)/2,y-kConeOutRadius+kLay1+(kLay2/2),z,new TGeoRotation("plate31",0,0,0)));
    mechLaddVol->AddNode(plate31,2,new TGeoCombiTrans(x-0.25-kConeOutRadius-kLay1-(0.75-0.25-kConeOutRadius-kLay1)/2,y-kConeOutRadius+kLay1+(kLay2/2),z,new TGeoRotation("plate31",0,0,0)));

    TGeoBBox *box32 = new TGeoBBox((kLay2/2),(kConeOutRadius-kLay1)/2,kStaveLength-0.50);
    TGeoVolume *plate32 = new TGeoVolume("CFleeceVertical",box32,medCarbonFleece);
    plate32->SetFillColor(28);
    plate32->SetLineColor(28);
    mechLaddVol->AddNode(plate32,1,new TGeoCombiTrans(x+0.25+kConeOutRadius+kLay1+(kLay2/2),y+(kLay1-kConeOutRadius)/2,z,new TGeoRotation("plate32",0,0,0)));
    mechLaddVol->AddNode(plate32,2,new TGeoCombiTrans(x+0.25-kConeOutRadius-kLay1-(kLay2/2),y+(kLay1-kConeOutRadius)/2,z,new TGeoRotation("plate32",0,0,0)));
    mechLaddVol->AddNode(plate32,3,new TGeoCombiTrans(x-0.25+kConeOutRadius+kLay1+(kLay2/2),y+(kLay1-kConeOutRadius)/2,z,new TGeoRotation("plate32",0,0,0)));
    mechLaddVol->AddNode(plate32,4,new TGeoCombiTrans(x-0.25-kConeOutRadius-kLay1-(kLay2/2),y+(kLay1-kConeOutRadius)/2,z,new TGeoRotation("plate32",0,0,0)));


    //K13D2U carbon plate
    TGeoBBox *box1 = new TGeoBBox(2*kWidth,kLay3/2,kStaveLength-0.50);
    TGeoVolume *plate1 = new TGeoVolume("CarbonPlate",box1,medK13D2U2k);
    plate1->SetFillColor(5);
    plate1->SetLineColor(5);
    mechLaddVol->AddNode(plate1,1,new TGeoCombiTrans(x,y-(kConeOutRadius+(kLay3/2)),z,new TGeoRotation("plate1",0,0,0)));

    //C Fleece bottom plate 
    TGeoBBox *box6 = new TGeoBBox(2*kWidth,kLay2/2,kStaveLength-0.50);
    TGeoVolume *plate6 = new TGeoVolume("CFleeceBottom",box6,medCarbonFleece);
    plate6->SetFillColor(2);
    plate6->SetLineColor(2);
    mechLaddVol->AddNode(plate6,1,new TGeoCombiTrans(x,y-(kConeOutRadius+kLay3+(kLay2/2)),z,new TGeoRotation("plate1",0,0,0)));
      
      
  }

  if (fBuildLevel < 2) {
    //Glue layers and kapton
    TGeoBBox *glue = new TGeoBBox(kStaveWidth/2, 0.005/2, zlad);
    TGeoVolume *volGlue=new TGeoVolume("Glue", glue, medGlue);
    volGlue->SetLineColor(5);
    volGlue->SetFillColor(5); 
    mechLaddVol->AddNode(volGlue, 0, new TGeoCombiTrans(x,y-(kConeOutRadius+kLay3+(kLay2/2)+(0.01/2)), z, new TGeoRotation("",0, 0, 0)));
    mechLaddVol->AddNode(volGlue, 1, new TGeoCombiTrans(x,y-(kConeOutRadius+kLay3+(kLay2/2)+0.01+fSensorThick+(0.01/2)), z, new TGeoRotation("",0, 0, 0)));
  }

  if (fBuildLevel < 1) {
    TGeoBBox *kapCable = new TGeoBBox(kStaveWidth/2, 0.01/2, zlad);
    TGeoVolume *volCable=new TGeoVolume("FlexCable", kapCable, medFlexCable);
    volCable->SetLineColor(28);
    volCable->SetFillColor(28); 
    mechLaddVol->AddNode(volCable, 0, new TGeoCombiTrans(x, y-(kConeOutRadius+kLay3+(kLay2/2)+0.01+fSensorThick+0.01+(0.01/2)), z, new TGeoRotation("",0, 0, 0)));
  }
    

  // Done, return the stave structure
  return mechLaddVol;
  
}
// new model22
//________________________________________________________________________
TGeoVolume* AliITSUv1Layer::CreateStaveModelInnerB22(const Double_t xlad,
						     const Double_t zlad,
						     const TGeoManager *mgr){
//
// Create the mechanical stave structure for Model 2.2 of TDR
//
// Input:
//         xlad : X length
//         zlad : Z length
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
  
  // Materials defined in AliITSUv1
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
  Double_t kConeOutRadius =0.107/2;//0.107/2;
  Double_t kConeInRadius = 0.1015/2;//0.10105/2
  Double_t kStaveLength = zlad;
  Double_t kStaveWidth = xlad*2;
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
  snprintf(volname, 30, "%s%d_StaveStruct", AliITSUGeomTGeo::GetITSLadderPattern(), fLayerNumber);

  Double_t z=0, y=-(2*kConeOutRadius)+klay1+klay2+fSensorThick/2-0.0004, x=0;

  TGeoVolume *mechLaddVol = 0;

  if (fBuildLevel < 5) {
    // world (trapezoid)
    TGeoXtru *mechStruct = new TGeoXtru(2); //z sections
    Double_t xv[6] = {kStaveWidth/2,kStaveWidth/2,0.012,-0.012,-kStaveWidth/2,-kStaveWidth/2}; 
    /* Double_t yv[6] = {-2*(kConeOutRadius+klay1+1.5*klay2+klay3+klay4+fSensorThick+klay5),
                        0-0.02,kStaveHeight+0.01,kStaveHeight+0.01,0-0.02,
			-2*(kConeOutRadius+klay1+1.5*klay2+klay3+klay4+fSensorThick+klay5)};  // (kConeOutRadius*2)-0.0635 */
    Double_t yv[6] = {-(kConeOutRadius*2)-0.06395,0-0.02,kStaveHeight+0.01,kStaveHeight+0.01,0-0.02,-(kConeOutRadius*2)-0.06395};  // (kConeOutRadius*2)-0.064
    mechStruct->DefinePolygon(6,xv,yv);
    mechStruct->DefineSection(0,-kStaveLength,0,0,1.);
    mechStruct->DefineSection(1,kStaveLength,0,0,1.);

    mechLaddVol = new TGeoVolume(volname, mechStruct, medAir);
    mechLaddVol->SetLineColor(12);
    mechLaddVol->SetFillColor(12); 
    mechLaddVol->SetVisibility(kTRUE);  
      
    //Polyimide Pipe Kapton grey-35 
    TGeoCone *cone1 = new TGeoCone(kStaveLength,kConeInRadius,kConeOutRadius-0.0001,kConeInRadius,kConeOutRadius-0.0001);
    TGeoVolume *volCone1= new TGeoVolume("PolyimidePipe", cone1, medKapton);
    volCone1->SetFillColor(35);
    volCone1->SetLineColor(35);
    mechLaddVol->AddNode(volCone1,1,new TGeoTranslation(x+0.25,y,z));
    mechLaddVol->AddNode(volCone1,2,new TGeoTranslation(x-0.25,y,z));
    }

  if (fBuildLevel < 4) {
    TGeoTube *coolTubeW = new TGeoTube(0.,kConeInRadius-0.0001,kStaveLength);
    TGeoVolume *volCoolTubeW= new TGeoVolume("Water", coolTubeW, medWater);
    volCoolTubeW->SetFillColor(4);
    volCoolTubeW->SetLineColor(4);
    mechLaddVol->AddNode(volCoolTubeW,0,new TGeoTranslation(x-0.25,y,z));
    mechLaddVol->AddNode(volCoolTubeW,1,new TGeoTranslation(x+0.25,y,z));
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
       mechLaddVol->AddNode(volT2,i*4+1,new TGeoCombiTrans(x+kWidth+0.0036,y+kHeight+0.01,z-kStaveLength+0.1+(i*4*kL1)+kS1/2, new TGeoRotation("volT2",90,90-kAlpha,90-kBeta)));
      // 2) Front Right Top Filament
      mechLaddVol->AddNode(volT2,i*4+2,new TGeoCombiTrans(x-kWidth-0.0036,y+kHeight+0.01,z-kStaveLength+0.1+(i*4*kL1)+kS1/2, new TGeoRotation("volT2",90,-90+kAlpha,-90+kBeta)));
      // 3) Back Left  Top Filament
      mechLaddVol->AddNode(volT2,i*4+3,new TGeoCombiTrans(x+kWidth+0.0036,y+kHeight+0.01,z-kStaveLength+0.1+2*kL1+(i*4*kL1)+kS1/2, new TGeoRotation("volT2",90,-90+kAlpha,90-kBeta)));
      // 4) Back Right Top Filament
      mechLaddVol->AddNode(volT2,i*4+4,new TGeoCombiTrans(x-kWidth-0.0036,y+kHeight+0.01,z-kStaveLength+0.1+2*kL1+(i*4*kL1)+kS1/2, new TGeoRotation("volT2",90,90-kAlpha,-90+kBeta)));
   }
 
     //Vertex  structure 

      //top ver trd1
      TGeoTrd1 *trd1 = new TGeoTrd1(0,kTopVertexMaxWidth/2,kStaveLength,kTopVertexHeight/2);
      TGeoVolume *ibdv = new TGeoVolume("TopVertex",trd1,medM60J3K);
      ibdv->SetFillColor(12);
      ibdv->SetLineColor(12);
      mechLaddVol->AddNode(ibdv,1,new TGeoCombiTrans(x,y+kStaveHeight+0.03,z,new TGeoRotation("ibdv",0.,-90,0)));//y+kStaveHeight+0.056

      //left trd2
      TGeoTrd1 *trd2 = new TGeoTrd1(0,kSideVertexMWidth/2,kStaveLength, kSideVertexHeight/2);
      TGeoVolume *ibdv2 = new TGeoVolume("LeftVertex",trd2,medM60J3K);
      ibdv2->SetFillColor(12);
      ibdv2->SetLineColor(12);
      mechLaddVol->AddNode(ibdv2,1,new TGeoCombiTrans(x+kStaveWidth/2-0.06,y-0.0355,z,new TGeoRotation("ibdv2",-103.3,90,0))); //x-kStaveWidth/2-0.09 old Config.C y-0.0355,

      //right trd3
      TGeoTrd1 *trd3 = new TGeoTrd1(0,kSideVertexMWidth/2,kStaveLength, kSideVertexHeight/2);
      TGeoVolume *ibdv3 = new TGeoVolume("RightVertex",trd3,medM60J3K);
      ibdv3->SetFillColor(12);
      ibdv3->SetLineColor(12);
      mechLaddVol->AddNode(ibdv3,1,new TGeoCombiTrans(x-kStaveWidth/2+0.06,y-0.0355,z,new TGeoRotation("ibdv3",103.3,90,0))); //x-kStaveWidth/2+0.09 old Config.C
      
     //Carbon Fleece
      TGeoConeSeg *cons2 = new TGeoConeSeg(zlad,kConeOutRadius+klay1,kConeOutRadius+klay1+klay2,kConeOutRadius+klay1,kConeOutRadius+klay1+klay2,0,180); 
      TGeoVolume *cone12 = new TGeoVolume("CarbonFleecePipeCover",cons2,medCarbonFleece);
      cone12->SetFillColor(28);
      cone12->SetLineColor(28);
      mechLaddVol->AddNode(cone12,1,new TGeoCombiTrans(x+0.25,y,z,new TGeoRotation("cone12",0,0,0)));
      mechLaddVol->AddNode(cone12,2,new TGeoCombiTrans(x-0.25,y,z,new TGeoRotation("cone12",0,0,0)));

      TGeoBBox *box3 = new TGeoBBox((0.50-(2*(kConeOutRadius+klay1)))/2,klay2/2,zlad);//kStaveLength-0.50);
      TGeoVolume *plate3 = new TGeoVolume("CarbonFleeceMiddle",box3,medCarbonFleece);
      plate3->SetFillColor(28);
      plate3->SetLineColor(28);
      mechLaddVol->AddNode(plate3,1,new TGeoCombiTrans(x,y-kConeOutRadius+klay1+(klay2/2),z,new TGeoRotation("plate3",0,0,0)));

      TGeoBBox *box31 = new TGeoBBox((0.75-0.25-kConeOutRadius-klay1)/2+0.0025,klay2/2,zlad);
      TGeoVolume *plate31 = new TGeoVolume("CarbonFleeceLeftRight",box31,medCarbonFleece);
      plate31->SetFillColor(28);
      plate31->SetLineColor(28);
      mechLaddVol->AddNode(plate31,1,new TGeoCombiTrans(x+0.25+kConeOutRadius+klay1+(0.75-0.25-kConeOutRadius-klay1)/2,y-kConeOutRadius+klay1+(klay2/2),z,new TGeoRotation("plate31",0,0,0)));
      mechLaddVol->AddNode(plate31,2,new TGeoCombiTrans(x-0.25-kConeOutRadius-klay1-(0.75-0.25-kConeOutRadius-klay1)/2,y-kConeOutRadius+klay1+(klay2/2),z,new TGeoRotation("plate31",0,0,0)));

      TGeoBBox *box32 = new TGeoBBox((klay2/2),(kConeOutRadius-klay1)/2,zlad);
      TGeoVolume *plate32 = new TGeoVolume("CarbonFleeceVertical",box32,medCarbonFleece);
      plate32->SetFillColor(28);
      plate32->SetLineColor(28);
      mechLaddVol->AddNode(plate32,1,new TGeoCombiTrans(x+0.25+kConeOutRadius+klay1+(klay2/2),y+(klay1-kConeOutRadius)/2,z,new TGeoRotation("plate32",0,0,0)));
      mechLaddVol->AddNode(plate32,2,new TGeoCombiTrans(x+0.25-kConeOutRadius-klay1-(klay2/2),y+(klay1-kConeOutRadius)/2,z,new TGeoRotation("plate32",0,0,0)));
      mechLaddVol->AddNode(plate32,3,new TGeoCombiTrans(x-0.25+kConeOutRadius+klay1+(klay2/2),y+(klay1-kConeOutRadius)/2,z,new TGeoRotation("plate32",0,0,0)));
      mechLaddVol->AddNode(plate32,4,new TGeoCombiTrans(x-0.25-kConeOutRadius-klay1-(klay2/2),y+(klay1-kConeOutRadius)/2,z,new TGeoRotation("plate32",0,0,0)));

     //Amec Thermasol red-2 cover tube FGS300 or Carbon Paper
      TGeoConeSeg *cons1 = new TGeoConeSeg(zlad,kConeOutRadius,kConeOutRadius+klay1-0.0001,kConeOutRadius,kConeOutRadius+klay1-0.0001,0,180);//kConeOutRadius+klay1-0.0001
      TGeoVolume *cone11 = new TGeoVolume("ThermasolPipeCover",cons1,medFGS003);
      cone11->SetFillColor(2);
      cone11->SetLineColor(2);
      mechLaddVol->AddNode(cone11,1,new TGeoCombiTrans(x+0.25,y,z,new TGeoRotation("cone11",0,0,0)));
      mechLaddVol->AddNode(cone11,2,new TGeoCombiTrans(x-0.25,y,z,new TGeoRotation("cone11",0,0,0)));

      TGeoBBox *box2 = new TGeoBBox((0.50-(2*kConeOutRadius))/2,(klay1/2),zlad);//kStaveLength-0.50);
      TGeoVolume *plate2 = new TGeoVolume("ThermasolMiddle",box2,medFGS003);
      plate2->SetFillColor(2);
      plate2->SetLineColor(2);
      mechLaddVol->AddNode(plate2,1,new TGeoCombiTrans(x,y-kConeOutRadius+(klay1/2),z,new TGeoRotation("plate2",0,0,0)));

      TGeoBBox *box21 = new TGeoBBox((0.75-0.25-kConeOutRadius-klay1)/2+0.0025,(klay1/2),zlad);
      TGeoVolume *plate21 = new TGeoVolume("ThermasolLeftRight",box21,medFGS003);
      plate21->SetFillColor(2);
      plate21->SetLineColor(2);
      mechLaddVol->AddNode(plate21,1,new TGeoCombiTrans(x+0.25+kConeOutRadius+(0.75-0.25-kConeOutRadius)/2-(klay1/2)+0.0025,y-kConeOutRadius+(klay1/2),z,new TGeoRotation("plate21",0,0,0)));
      mechLaddVol->AddNode(plate21,2,new TGeoCombiTrans(x-0.25-kConeOutRadius-(0.75-0.25-kConeOutRadius)/2+(klay1/2)-0.0025,y-kConeOutRadius+(klay1/2),z,new TGeoRotation("plate21",0,0,0)));

      TGeoBBox *box22 = new TGeoBBox((klay1/2),kConeOutRadius/2,zlad);
      TGeoVolume *plate22 = new TGeoVolume("ThermasolVertical",box22,medFGS003);
      plate22->SetFillColor(2);
      plate22->SetLineColor(2);
      mechLaddVol->AddNode(plate22,1,new TGeoCombiTrans(x+0.25+kConeOutRadius+(klay1/2),y-kConeOutRadius/2,z,new TGeoRotation("plate22",0,0,0)));
      mechLaddVol->AddNode(plate22,2,new TGeoCombiTrans(x+0.25-kConeOutRadius-(klay1/2),y-kConeOutRadius/2,z,new TGeoRotation("plate22",0,0,0)));
      mechLaddVol->AddNode(plate22,3,new TGeoCombiTrans(x-0.25+kConeOutRadius+(klay1/2),y-kConeOutRadius/2,z,new TGeoRotation("plate22",0,0,0)));
      mechLaddVol->AddNode(plate22,4,new TGeoCombiTrans(x-0.25-kConeOutRadius-(klay1/2),y-kConeOutRadius/2,z,new TGeoRotation("plate22",0,0,0)));

     //K13D2U CF plate
      TGeoBBox *box1 = new TGeoBBox(2*kWidth,(klay3)/2,zlad);
      TGeoVolume *plate1 = new TGeoVolume("CFPlate",box1,medK13D2U2k);
      plate1->SetFillColor(5);
      plate1->SetLineColor(5);
      mechLaddVol->AddNode(plate1,1,new TGeoCombiTrans(x,y-(kConeOutRadius+(klay3/2)),z,new TGeoRotation("plate1",0,0,0)));

     //C Fleece bottom plate 
      TGeoBBox *box6 = new TGeoBBox(2*kWidth,(klay2)/2,zlad);
      TGeoVolume *plate6 = new TGeoVolume("CarbonFleeceBottom",box6,medCarbonFleece);
      plate6->SetFillColor(2);
      plate6->SetLineColor(2);
      mechLaddVol->AddNode(plate6,1,new TGeoCombiTrans(x,y-(kConeOutRadius+klay3+(klay2/2)),z,new TGeoRotation("plate6",0,0,0)));

    }
   if (fBuildLevel < 2) {
      //Glue klayers and kapton
     TGeoBBox *glue = new TGeoBBox(kStaveWidth/2, (klay4)/2, zlad);
      TGeoVolume *volGlue=new TGeoVolume("Glue", glue, medGlue);
      volGlue->SetLineColor(5);
      volGlue->SetFillColor(5); 
      // mechLaddVol->AddNode(volGlue, 0, new TGeoCombiTrans(x,y-(kConeOutRadius+klay3+klay2+(klay4/2)), z, new TGeoRotation("",0, 0, 0)));
      mechLaddVol->AddNode(volGlue, 0, new TGeoCombiTrans(x,y-(kConeOutRadius+klay3+klay2+(klay4)/2)+0.00005, z, new TGeoRotation("",0, 0, 0)));
    }

     if (fBuildLevel < 1) {
     //Flex Cable or Bus
      TGeoBBox *kapCable = new TGeoBBox(kStaveWidth/2, klay5/2, zlad);//klay5/2
      TGeoVolume *volCable=new TGeoVolume("FlexCable", kapCable, medFlexCable);
      volCable->SetLineColor(28);
      volCable->SetFillColor(28); 
      //      mechLaddVol->AddNode(volCable, 0, new TGeoCombiTrans(x, y-(kConeOutRadius+klay3+klay2+klay4+fSensorThick+(klay5)/2)+0.0002, z, new TGeoRotation("",0, 0, 0)));
      mechLaddVol->AddNode(volCable, 0, new TGeoCombiTrans(x, y-(kConeOutRadius+klay3+klay2+klay4+fSensorThick+(klay5)/2)+0.01185, z, new TGeoRotation("",0, 0, 0)));
      }
    // Done, return the stave structe
    return mechLaddVol;
}

// model3
//________________________________________________________________________
TGeoVolume* AliITSUv1Layer::CreateStaveModelInnerB3(const Double_t xlad,
						    const Double_t zlad,
						    const TGeoManager *mgr){
//
// Create the mechanical stave structure for Model 3 of TDR
//
// Input:
//         xlad : X length
//         zlad : Z length
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
  
  // Materials defined in AliITSUv1
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
    Double_t kStaveLength = zlad*2;
    Double_t kStaveWidth = xlad*2;
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
       Int_t modules = 4;
       Double_t headWidth=0.25;
       Double_t smcLength=kStaveLength/modules-2*headWidth;//6.25;
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
    snprintf(volname, 30, "%s%d_StaveStruct", AliITSUGeomTGeo::GetITSLadderPattern(), fLayerNumber);
    
    // detailed structure ++++++++++++++
    Double_t z=0, y=0-0.007, x=0;

    // Polimide micro channels numbers
    Double_t yMC = y-h+0.01;
    Int_t nb = (Int_t)(kStaveWidth/0.1)+1;
    Double_t xladMC = (nb*0.1-0.08)/2;


    TGeoVolume *mechLaddVol = 0;
    if (fBuildLevel < 5) {
      // world (trapezoid)
      TGeoXtru *mechStruct = new TGeoXtru(2); //z sections
      Double_t xv[5] = {kStaveWidth/2+0.1,kStaveWidth/2+0.1,0,-kStaveWidth/2-0.1,-kStaveWidth/2-0.1};
      Double_t yv[5] = {-kConeOutRadius*2-0.07,0,staveHeight,0,-kConeOutRadius*2-0.07};    
      mechStruct->DefinePolygon(5,xv,yv);
      mechStruct->DefineSection(0,-kStaveLength-0.1,0,0,1.);
      mechStruct->DefineSection(1,kStaveLength+0.1,0,0,1.);
      mechLaddVol = new TGeoVolume(volname, mechStruct, medAir);
      mechLaddVol->SetLineColor(12);
      mechLaddVol->SetFillColor(12); 
      mechLaddVol->SetVisibility(kTRUE);

       // Silicon micro channels numbers
      
      TGeoBBox *tM0a=new TGeoBBox(smcWidth/2, 0.003/2, headWidth/2);
      TGeoVolume *volTM0a=new TGeoVolume("microChanTop1", tM0a, medKapton);
      volTM0a->SetLineColor(35);
      volTM0a->SetFillColor(35); 

      for(Int_t  mo=1; mo<=modules; mo++) {
      mechLaddVol->AddNode(volTM0a, 0, new TGeoCombiTrans(x,yMC+0.03, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth+smcLength/2+(headWidth/2), new TGeoRotation("",ang1, ang2, ang3)));
      mechLaddVol->AddNode(volTM0a, 1, new TGeoCombiTrans(x,yMC+0.03, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth-smcLength/2-(headWidth/2), new TGeoRotation("",ang1, ang2, ang3)));
      }
      TGeoBBox *tM0c=new TGeoBBox(0.3/2, 0.003/2,smcLength/2);
      TGeoVolume *volTM0c=new TGeoVolume("microChanTop2", tM0c, medKapton);
      volTM0c->SetLineColor(35);
      volTM0c->SetFillColor(35); 
      for(Int_t  mo=1; mo<=modules; mo++) {
      mechLaddVol->AddNode(volTM0c, 0, new TGeoCombiTrans(x+(smcWidth/2)-(0.3/2),yMC+0.03, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));
      mechLaddVol->AddNode(volTM0c, 1, new TGeoCombiTrans(x-(smcWidth/2)+(0.3/2),yMC+0.03, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0c1=new TGeoBBox(0.2225/2, 0.003/2,smcLength/2);
      TGeoVolume *volTM0c1=new TGeoVolume("microChanBot1", tM0c1, medKapton);
      volTM0c1->SetLineColor(6);
      volTM0c1->SetFillColor(6); 
      for(Int_t  mo=1; mo<=modules; mo++) {
      mechLaddVol->AddNode(volTM0c1, 0, new TGeoCombiTrans(x+smcWidth/2-(smcSide1Thick)-(vaporThick)-(smcSide2Thick)-(smcSide3Thick)-(0.2225/2),yMC+0.03-hh-(0.003), z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      mechLaddVol->AddNode(volTM0c1, 1, new TGeoCombiTrans(x-smcWidth/2+(smcSide1Thick)+(liquidThick)+(smcSide2Thick)+(smcSide4Thick)+(0.2225/2),yMC+0.03-hh-(0.003), z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0c2=new TGeoBBox(0.072/2, 0.003/2,smcLength/2);
      TGeoVolume *volTM0c2=new TGeoVolume("microChanBot2", tM0c2, medKapton);
      volTM0c2->SetLineColor(35);
      volTM0c2->SetFillColor(35); 
      for(Int_t  mo=1; mo<=modules; mo++) {
      mechLaddVol->AddNode(volTM0c2, 0, new TGeoCombiTrans(x+smcWidth/2-(0.072/2),yMC+0.03-(0.035+0.0015)-(0.003)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0c2r=new TGeoBBox(0.068/2, 0.003/2,smcLength/2);
      TGeoVolume *volTM0c2r=new TGeoVolume("microChanBot3", tM0c2r, medKapton);
      volTM0c2r->SetLineColor(35);
      volTM0c2r->SetFillColor(35); 
      for(Int_t  mo=1; mo<=modules; mo++) {      
      mechLaddVol->AddNode(volTM0c2r, 0, new TGeoCombiTrans(x-smcWidth/2+(0.068/2),yMC+0.03-(0.035+0.0015)-(0.003)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0d=new TGeoBBox(smcSide1Thick/2, 0.035/2,smcLength/2);
      TGeoVolume *volTM0d=new TGeoVolume("microChanSide1", tM0d, medKapton);
      volTM0d->SetLineColor(12);
      volTM0d->SetFillColor(12); 
      for(Int_t  mo=1; mo<=modules; mo++) {
      mechLaddVol->AddNode(volTM0d, 0, new TGeoCombiTrans(x+smcWidth/2-(smcSide1Thick/2),yMC+0.03-0.0015-(0.035)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      mechLaddVol->AddNode(volTM0d, 1, new TGeoCombiTrans(x-smcWidth/2+(smcSide1Thick/2),yMC+0.03-0.0015-(0.035)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }

      TGeoBBox *tM0d1=new TGeoBBox(smcSide2Thick/2, 0.035/2,smcLength/2);
      TGeoVolume *volTM0d1=new TGeoVolume("microChanSide2", tM0d1, medKapton);
      volTM0d1->SetLineColor(12);
      volTM0d1->SetFillColor(12); 
      for(Int_t  mo=1; mo<=modules; mo++) {
      mechLaddVol->AddNode(volTM0d1, 0, new TGeoCombiTrans(x+smcWidth/2-(smcSide1Thick)-(vaporThick)-(smcSide2Thick/2),yMC+0.03-(0.003+0.035)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      mechLaddVol->AddNode(volTM0d1, 1, new TGeoCombiTrans(x-smcWidth/2+(smcSide1Thick)+(liquidThick)+(smcSide2Thick/2),yMC+0.03-(0.003+0.035)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0d2=new TGeoBBox(smcSide3Thick/2, (hh+0.003)/2, smcLength/2);
      TGeoVolume *volTM0d2=new TGeoVolume("microChanSide3", tM0d2, medKapton);
      volTM0d2->SetLineColor(12);
      volTM0d2->SetFillColor(12); 
      for(Int_t  mo=1; mo<=modules; mo++) {
      mechLaddVol->AddNode(volTM0d2, 0, new TGeoCombiTrans(x+smcWidth/2-(smcSide1Thick)-(vaporThick)-(smcSide2Thick)-(smcSide3Thick/2),yMC+0.03-(0.003+hh+0.003)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0d2r=new TGeoBBox(smcSide4Thick/2, (hh+0.003)/2, smcLength/2);
      TGeoVolume *volTM0d2r=new TGeoVolume("microChanSide4", tM0d2r, medKapton);
      volTM0d2r->SetLineColor(12);
      volTM0d2r->SetFillColor(12); 
      for(Int_t  mo=1; mo<=modules; mo++) {
      mechLaddVol->AddNode(volTM0d2r, 0, new TGeoCombiTrans(x-smcWidth/2+(smcSide1Thick)+(liquidThick)+(smcSide2Thick)+(smcSide4Thick/2),yMC+0.03-(0.003+hh+0.003)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0e=new TGeoBBox(smcSide5Thick/2, hh/2,smcLength/2);
      TGeoVolume *volTM0e=new TGeoVolume("microChanSide5", tM0e, medKapton);    
      volTM0e->SetLineColor(12);
      volTM0e->SetFillColor(12); 
      for(Int_t  mo=1; mo<=modules; mo++) {
      for (Int_t ie=0;ie<11;ie++) {
	mechLaddVol->AddNode(volTM0e, 0, new TGeoCombiTrans(x-(ie*(smcSpace+smcSide5Thick))+smcWidth/2-(smcSide1Thick)-(vaporThick)-(smcSide2Thick)-(smcSide3Thick)-smcSpace-(smcSide5Thick/2),yMC+0.03-(0.003+hh)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
	mechLaddVol->AddNode(volTM0e, 1, new TGeoCombiTrans(x+(ie*(smcSpace+smcSide5Thick))-smcWidth/2+(smcSide1Thick)+(liquidThick)+(smcSide2Thick)+(smcSide4Thick)+smcSpace+(smcSide5Thick/2),yMC+0.03-(0.003+hh)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
         }
      }
      
      TGeoBBox *tM0f=new TGeoBBox(0.02/2, hh/2, smcLength/2);
      TGeoVolume *volTM0f=new TGeoVolume("microChanTop3", tM0f, medKapton);
      //Double_t smcChannels=12;
      Double_t smcCloseWallvapor=smcWidth/2-smcSide1Thick-vaporThick-smcSide2Thick-smcSide3Thick-12*smcSpace-11*smcSide5Thick;
      Double_t smcCloseWallliquid=smcWidth/2-smcSide1Thick-liquidThick-smcSide2Thick-smcSide4Thick-12*smcSpace-11*smcSide5Thick;
      volTM0f->SetLineColor(12);
      volTM0f->SetFillColor(12); 
      for(Int_t  mo=1; mo<=modules; mo++) {
       mechLaddVol->AddNode(volTM0f, 0, new TGeoCombiTrans(x+smcCloseWallvapor-(0.02)/2,yMC+0.03-(0.003+hh)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
       mechLaddVol->AddNode(volTM0f, 1, new TGeoCombiTrans(x-smcCloseWallliquid+(0.02)/2,yMC+0.03-(0.003+hh)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      //Head(back) microchannel

      TGeoBBox *tM0hb=new TGeoBBox(smcWidth/2, 0.025/2, headWidth/2);
      TGeoVolume *volTM0hb=new TGeoVolume("microChanHeadBackBottom1", tM0hb, medKapton);
      volTM0hb->SetLineColor(4);
      volTM0hb->SetFillColor(4); 
      for(Int_t  mo=1; mo<=modules; mo++) {
      mechLaddVol->AddNode(volTM0hb, 0, new TGeoCombiTrans(x,yMC+0.03-0.0145-(0.025/2), z+(mo-3)*kStaveLength/4+smcLength/2+headWidth+smcLength/2+(headWidth/2), new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      mechLaddVol->AddNode(volTM0hb, 1, new TGeoCombiTrans(x,yMC+0.03-0.0145-(0.025)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth-smcLength/2-(headWidth/2), new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0h1=new TGeoBBox(smcWidth/2, 0.013/2, 0.05/2);
      TGeoVolume *volTM0h1=new TGeoVolume("microChanHeadBackBottom2", tM0h1, medKapton);
      volTM0h1->SetLineColor(5);
      volTM0h1->SetFillColor(5); 
      for(Int_t  mo=1; mo<=modules; mo++) {
      mechLaddVol->AddNode(volTM0h1, 0, new TGeoCombiTrans(x,yMC+0.03-0.0015-(0.013/2), z+(mo-3)*kStaveLength/4+smcLength/2+headWidth-smcLength/2-headWidth+(0.05/2), new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0h2=new TGeoBBox(smcWidth/2, 0.003/2, 0.18/2);
      TGeoVolume *volTM0h2=new TGeoVolume("microChanHeadBackBottom7", tM0h2, medKapton);
      volTM0h2->SetLineColor(6);
      volTM0h2->SetFillColor(6);
      for(Int_t  mo=1; mo<=modules; mo++) {
      mechLaddVol->AddNode(volTM0h2, 0, new TGeoCombiTrans(x,yMC+0.03-0.0015-0.01-(0.003/2), z+(mo-3)*kStaveLength/4+smcLength/2+headWidth-smcLength/2-0.02-(0.18/2), new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0h3=new TGeoBBox(smcWidth/2, 0.013/2, 0.02/2);
      TGeoVolume *volTM0h3=new TGeoVolume("microChanHeadBackBottom3", tM0h3, medKapton);
      volTM0h3->SetLineColor(5);
      volTM0h3->SetFillColor(5); 
      for(Int_t  mo=1; mo<=modules; mo++) {
      mechLaddVol->AddNode(volTM0h3, 0, new TGeoCombiTrans(x,yMC+0.03-0.0015-(0.013/2), z+(mo-3)*kStaveLength/4+smcLength/2+headWidth-smcLength/2-(0.02/2), new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0b1=new TGeoBBox(smcWidth/2, 0.013/2, 0.03/2);
      TGeoVolume *volTM0b1=new TGeoVolume("microChanHeadBackBottom4", tM0b1, medKapton);
      volTM0b1->SetLineColor(5);
      volTM0b1->SetFillColor(5); 
      for(Int_t  mo=1; mo<=modules; mo++) {
      mechLaddVol->AddNode(volTM0b1, 0, new TGeoCombiTrans(x,yMC+0.03-0.0015-(0.013/2), z+(mo-3)*kStaveLength/4+smcLength/2+headWidth+smcLength/2+headWidth-(0.03/2), new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0b2=new TGeoBBox(smcWidth/2, 0.003/2, 0.2/2);
      TGeoVolume *volTM0b2=new TGeoVolume("microChanHeadBackBottom5", tM0b2, medKapton);
      volTM0b2->SetLineColor(6);
      volTM0b2->SetFillColor(6); 
      for(Int_t  mo=1; mo<=modules; mo++) {
      mechLaddVol->AddNode(volTM0b2, 0, new TGeoCombiTrans(x,yMC+0.03-0.0015-0.01-(0.003/2), z+(mo-3)*kStaveLength/4+smcLength/2+headWidth+smcLength/2+0.02+(0.2/2), new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0b3=new TGeoBBox(smcWidth/2, 0.013/2, 0.02/2);
      TGeoVolume *volTM0b3=new TGeoVolume("microChanHeadBackBottom6", tM0b3, medKapton);
      volTM0b3->SetLineColor(5);
      volTM0b3->SetFillColor(5); 
      for(Int_t  mo=1; mo<=modules; mo++) {
      mechLaddVol->AddNode(volTM0b3, 0, new TGeoCombiTrans(x,yMC+0.03-0.0015-(0.013/2), z+(mo-3)*kStaveLength/4+smcLength/2+headWidth+smcLength/2+(0.02/2), new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
     
      TGeoBBox *tM0b=new TGeoBBox(0.02/2, 0.02/2, zlad);
      TGeoVolume *volTM0b=new TGeoVolume("microChanWalls", tM0b, medKapton);
      volTM0b->SetLineColor(35);
      volTM0b->SetFillColor(35); 
      for (Int_t ib=0;ib<nb;ib++) {
	//mechLaddVol->AddNode(volTM0b, ib, new TGeoCombiTrans(x+ib*0.1-xladMC+0.01,yMC, z, new TGeoRotation("",0, 0, 0)));
      }
      
      } 
    
    if (fBuildLevel < 4) {

      //**********cooling  inlet outlet

      TGeoBBox *tM0dv=new TGeoBBox(vaporThick/2, 0.035/2,smcLength/2);
      TGeoVolume *volTM0dv=new TGeoVolume("microChanVapor", tM0dv, medWater);
      volTM0dv->SetLineColor(2);
      volTM0dv->SetFillColor(2);
      for(Int_t  mo=1; mo<=modules; mo++) {
      mechLaddVol->AddNode(volTM0dv, 0, new TGeoCombiTrans(x+smcWidth/2-(smcSide1Thick)-(vaporThick/2),yMC+0.03-0.0015-(0.035)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
      }
      TGeoBBox *tM0dl=new TGeoBBox(liquidThick/2, 0.035/2,smcLength/2);
      TGeoVolume *volTM0dl=new TGeoVolume("microChanLiquid", tM0dl, medWater);
      volTM0dl->SetLineColor(3);
      volTM0dl->SetFillColor(3); 
      for(Int_t  mo=1; mo<=modules; mo++) {
      mechLaddVol->AddNode(volTM0dl, 0, new TGeoCombiTrans(x-smcWidth/2+(smcSide1Thick)+(liquidThick/2),yMC+0.03-0.0015-(0.035)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
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
      for(Int_t  mo=1; mo<=modules; mo++) {
      for (Int_t is=0;is<12;is++) {
	mechLaddVol->AddNode(volTM0dlq, 0, new TGeoCombiTrans(x+(is*(smcSpace+smcSide5Thick))-smcWidth/2+(smcSide1Thick)+(vaporThick)+(smcSide2Thick)+(smcSide3Thick)+smcSpace/2,yMC+0.03-(0.003+hh)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
	mechLaddVol->AddNode(volTM0dvp, 1, new TGeoCombiTrans(x-(is*(smcSpace+smcSide5Thick))+smcWidth/2-(smcSide1Thick)-(vaporThick)-(smcSide2Thick)-(smcSide3Thick)-smcSpace/2,yMC+0.03-(0.003+hh)/2, z+(mo-3)*kStaveLength/4+smcLength/2+headWidth, new TGeoRotation("",ang1, ang2, ang3)));//("",0, 0, 0)));
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
       	mechLaddVol->AddNode(volT1,4*i+0,
				    new TGeoCombiTrans(x+w,y-h+0.04+filHeight/2,z-kStaveLength/2+(4*l*i)+s1/2, 
						       new TGeoRotation("volT1",-90,alpha,0)));
	mechLaddVol->AddNode(volT1,4*i+1,
				    new TGeoCombiTrans(x-w,y-h+0.04+filHeight/2,z-kStaveLength/2+(4*l*i)+s1/2, 
						       new TGeoRotation("volT1",90,alpha,0)));
	mechLaddVol->AddNode(volT1,4*i+2,
				    new TGeoCombiTrans(x+w,y-h+0.04+filHeight/2,z-kStaveLength/2+2*l+(i*4*l)+s1/2, 
						       new TGeoRotation("volT1",-90,-alpha,0)));
	mechLaddVol->AddNode(volT1,4*i+3,
				    new TGeoCombiTrans(x-w,y-h+0.04+filHeight/2,z-kStaveLength/2+2*l+(i*4*l)+s1/2, 
						       new TGeoRotation("volT1",-90,+alpha,0)));
	}
 
     // Top filament CERP black-12 Carbon structure TGeoBBox (length,thickness,width)

      TGeoBBox *t2=new TGeoBBox(s2,filHeight/2,filWidth/2);
      TGeoVolume *volT2=new TGeoVolume("topFilament", t2, medM60J3K);
      volT2->SetLineColor(12);
      volT2->SetFillColor(12); 
      for(int i=0;i<loop;i++){ //i<30;i++){
       	mechLaddVol->AddNode(volT2,4*i+0,
				    new TGeoCombiTrans(x+w,y+0.04+filHeight/2,z-kStaveLength/2+(i*4*l)+s1/2,
						       new TGeoRotation("volT2",90,90-alpha,90-beta)));
	mechLaddVol->AddNode(volT2,4*i+1,
				    new TGeoCombiTrans(x-w,y+0.04+filHeight/2,z-kStaveLength/2+(i*4*l)+s1/2,
						       new TGeoRotation("volT2",90,-90+alpha,-90+beta)));
	mechLaddVol->AddNode(volT2,4*i+2,
				    new TGeoCombiTrans(x+w,y+0.04+filHeight/2,z-kStaveLength/2+2*l+(i*4*l)+s1/2,
						       new TGeoRotation("volT2",90,-90+alpha,90-beta)));
	mechLaddVol->AddNode(volT2,4*i+3,
				    new TGeoCombiTrans(x-w,y+0.04+filHeight/2,z-kStaveLength/2+2*l+(i*4*l)+s1/2, 
						       new TGeoRotation("volT2",90,90-alpha,-90+beta)));
	}
    }

    if (fBuildLevel < 2) {

      // Glue Filament and Silicon MicroChannel
      TGeoBBox *tM0=new TGeoBBox(xladMC/5, klay4/2, zlad);
      TGeoVolume *volTM0=new TGeoVolume("glueFM", tM0,medGlue );
      volTM0->SetLineColor(5);
      volTM0->SetFillColor(5); 
      mechLaddVol->AddNode(volTM0, 0, new TGeoCombiTrans(x-xlad/2-0.25,0.03+yMC, z, new TGeoRotation("",0, 0, 0)));
      mechLaddVol->AddNode(volTM0, 1, new TGeoCombiTrans(x+xlad/2+0.25,0.03+yMC, z, new TGeoRotation("",0, 0, 0)));

            
      // Glue microchannel and sensor
      TGeoBBox *glueM = new TGeoBBox(xladMC/5, klay4/2, zlad);
      TGeoVolume *volGlueM=new TGeoVolume("glueMS", glueM, medGlue);
      volGlueM->SetLineColor(5);
      volGlueM->SetFillColor(5); 
      mechLaddVol->AddNode(volGlueM, 0, new TGeoCombiTrans(x-xlad/2-0.25,yMC-0.01, z, new TGeoRotation("",0, 0, 0)));
      mechLaddVol->AddNode(volGlueM, 1, new TGeoCombiTrans(x+xlad/2+0.25,yMC-0.01, z, new TGeoRotation("",0, 0, 0)));
     
       // Glue sensor and kapton
      TGeoBBox *glue = new TGeoBBox(xlad, klay4/2, zlad);
      TGeoVolume *volGlue=new TGeoVolume("glueSensorBus", glue, medGlue);
      volGlue->SetLineColor(5);
      volGlue->SetFillColor(5); 
       mechLaddVol->AddNode(volGlue, 1, new TGeoCombiTrans(x, y-0.154-fSensorThick-klay4/2, z, new TGeoRotation("",0, 0, 0)));
    }      

    if (fBuildLevel < 1) {
      TGeoBBox *kapCable = new TGeoBBox(xlad, klay5/2, zlad);
      TGeoVolume *volCable=new TGeoVolume("Flexcable", kapCable, medFlexCable);
      volCable->SetLineColor(28);
      volCable->SetFillColor(28); 
      mechLaddVol->AddNode(volCable, 0, new TGeoCombiTrans(x, y-0.154-fSensorThick-klay4-klay5/2, z, new TGeoRotation("",0, 0, 0)));
    }

  // Done, return the stave structur
    return mechLaddVol;
 }

//________________________________________________________________________
TGeoVolume* AliITSUv1Layer::CreateStaveOuterB(const Double_t xlad,
					      const TGeoManager *mgr){
//
// Create the module stave for the Outer Barrel
//
// Input:
//         xlad : X length
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      20 Dec 2013  Mario Sitta
//

  TGeoVolume *mechLaddVol = 0;

  switch (fStaveModel) {
    case AliITSUv1::kOBModelDummy:
      mechLaddVol = CreateStaveModelOuterBDummy(xlad,mgr);
      break;
    case AliITSUv1::kOBModel0:
      mechLaddVol = CreateStaveModelOuterB0(xlad,mgr);
      break;
    case AliITSUv1::kOBModel1:
      mechLaddVol = CreateStaveModelOuterB1(xlad,mgr);
      break;
    default:
      AliFatal(Form("Unknown stave model %d",fStaveModel));
      break;
  }

  return mechLaddVol; 
}

//________________________________________________________________________
TGeoVolume* AliITSUv1Layer::CreateStaveModelOuterBDummy(const Double_t ,
							const TGeoManager *) const {
//
// Create dummy stave
//
// Input:
//         xlad : X length
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
TGeoVolume* AliITSUv1Layer::CreateStaveModelOuterB0(const Double_t ,
						    const TGeoManager *) const {
//
// Creation of the mechanical stave structure for the Outer Barrel as in v0
// is done directly in CreateLadder, so this method does nothing
// (doing it there is simpler, since all needed dimensions are known)
//
// Input:
//         xlad : X length
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
TGeoVolume* AliITSUv1Layer::CreateStaveModelOuterB1(const Double_t xlad,
						    const TGeoManager *mgr){
//
// Create the mechanical stave structure for the Outer Barrel as in TDR
//
// Input:
//         xlad : X length
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      20 Nov 2013  Anastasia Barbano
// Updated:      16 Jan 2014  Mario Sitta
//

  // Materials defined in AliITSUv0
  TGeoMedium *medAluminum     = mgr->GetMedium("ITS_ALUMINUM$");
  TGeoMedium *medCarbon       = mgr->GetMedium("ITS_CARBON$");
  TGeoMedium *medKapton       = mgr->GetMedium("ITS_KAPTON(POLYCH2)$");
  TGeoMedium *medGlue         = mgr->GetMedium("ITS_GLUE$");
  TGeoMedium *medWater        = mgr->GetMedium("ITS_WATER$");
  TGeoMedium *medCarbonFleece = mgr->GetMedium("ITS_CarbonFleece$"); 
  TGeoMedium *medFGS003       = mgr->GetMedium("ITS_FGS003$"); //amec thermasol
  TGeoMedium *medAir          = mgr->GetMedium("ITS_AIR$");


  // Local parameters
  Double_t modGap             = fgkOBModuleGap;
  Double_t yFlex1             = fgkOBFlexCable1Thick;
  Double_t yFlex2             = fgkOBFlexCable2Thick;
  Double_t yBus1              = fgkOBBusCable1Thick;
  Double_t yBus2              = fgkOBBusCable2Thick;
  Double_t xModPlate          = fgkOBHalfStaveWidth;
  Double_t yModPlate          = fgkOBCarbonPlateThick;
  Double_t xCPlate            = fgkOBHalfStaveWidth;
  Double_t yCPlate            = fgkOBColdPlateThick;
  Double_t yGlue              = fgkOBGlueThick;
  Double_t flexOverlap        = 5;
  Double_t deltaY             = 0.176;
  Double_t xOverlap           = 0.23;       //overlapping of the halfStaves to cover the dead zone of sensors
  Double_t zMod               = fgkOBModuleZLength;
  Double_t xHalfSt            = fgkOBHalfStaveWidth/2;
  Double_t xPos               = xOverlap/2 - xHalfSt;
  Double_t xlen               = xlad;
  Double_t rMin               = 0.267/2;
  Double_t rMax               = rMin + 0.0065;
  Double_t kLay1              = 0.004;      //carbon fleece
  Double_t kLay2              = 0.003;      //carbon paper
  Double_t yPos               = kLay1+kLay2;
  Double_t ylen,zact;
  Double_t zpos, zpos5cm;
  Double_t ymod;
  Double_t zbus;
  Double_t zlen;


  if (fIsTurbo) xlen = 0.5*fLadderWidth;
  //ylen    = 0.5*fLadderThick;
  ymod    = 0.005/2;//0.5*fSensorThick;
  ylen    = 0.5*(2*kLay1+2*kLay2+2*rMax+yCPlate+yGlue+ yModPlate + ymod + 2*yFlex1 + 2*yFlex2 + yBus1 + yBus2 + deltaY);
  zact    = fNModules*zMod; //active area 
  zbus    = zact + (fNModules-1)*modGap; 
  zlen    = zbus/2;


  // First create all needed shapes

  TGeoTube *coolTube   = new TGeoTube("CoolingTube",rMin,rMax,zbus/2);
  TGeoTube *coolTubeW  = new TGeoTube("CoolingTubeWater",0,rMin,zbus/2);
  TGeoBBox *coldPlate  = new TGeoBBox("ColdPlate",xCPlate/2,yCPlate/2,zbus/2);
  TGeoBBox *glue       = new TGeoBBox("Glue",xCPlate/2,yGlue/2,zbus/2);
  TGeoBBox *modPlate   = new TGeoBBox("CarbonPlate",xModPlate/2,yModPlate/2,zbus/2);
  TGeoBBox *flex1      = new TGeoBBox("Flex1MV",xHalfSt,yFlex1/2,zMod/2);
  TGeoBBox *flex2      = new TGeoBBox("Flex2MV",xHalfSt,yFlex2/2,zMod/2);
  TGeoBBox *flex1_5cm  = new TGeoBBox("Flex1MV_5cm",xHalfSt,yFlex1/2,flexOverlap/2);
  TGeoBBox *flex2_5cm  = new TGeoBBox("Flex2MV_5cm",xHalfSt,yFlex2/2,flexOverlap/2);
  TGeoBBox *bus1       = new TGeoBBox("Bus1HV",xHalfSt,yBus1/2,zbus/2);
  TGeoBBox *bus2       = new TGeoBBox("Bus2HV",xHalfSt,yBus2/2,zbus/2); 
  TGeoTubeSeg *cone1   = new TGeoTubeSeg(rMax +kLay2,rMax+kLay1+kLay2,zlen,180.,360.);  //Carbon Fleece 
  TGeoTubeSeg *cone2   = new TGeoTubeSeg(rMax ,rMax+kLay2,zlen,180.,360.);  //Graphite paper 
  TGeoBBox *box11      = new TGeoBBox((0.95-kLay2-rMax)/2,kLay1/2,zlen);
  TGeoBBox *box12      = new TGeoBBox((1.11-2*kLay2-2*rMax)/2,kLay1/2,zlen);
  TGeoBBox *box13      = new TGeoBBox(kLay1/2,(rMax-(kLay1+kLay2))/2,zlen);
  TGeoBBox *box21      = new TGeoBBox((0.95-rMax)/2,kLay2/2,zlen);
  TGeoBBox *box22      = new TGeoBBox((1.11-2*rMax)/2,kLay2/2,zlen);
  TGeoBBox *box23      = new TGeoBBox(kLay2/2,(rMax-kLay2)/2,zlen);
  TGeoBBox *mechStruct = new TGeoBBox("mechanicalStructure",xlen, ylen, 0.5*fZLength); 


  TGeoVolume *modVol       = CreateModuleOuterB(xHalfSt, ymod, zMod);
 
  TGeoVolume *coolTubeVol  = new TGeoVolume("CoolingTubeVol",coolTube,medKapton);
  TGeoVolume *coolTubeWVol = new TGeoVolume("CoolingTubeWaterVol",coolTubeW,medWater);
  TGeoVolume *coldPlateVol = new TGeoVolume("ColdPlateVol",coldPlate,medCarbon);
  TGeoVolume *glueVol      = new TGeoVolume("GlueVol",glue,medGlue);
  TGeoVolume *modPlateVol  = new TGeoVolume("CarbonPlateVol",modPlate,medCarbon);
  TGeoVolume *flex1Vol     = new TGeoVolume("Flex1Vol",flex1,medAluminum);
  TGeoVolume *flex2Vol     = new TGeoVolume("Flex2Vol",flex2,medKapton);
  TGeoVolume *flex1_5cmVol = new TGeoVolume("Flex1Vol5cm",flex1_5cm,medAluminum);
  TGeoVolume *flex2_5cmVol = new TGeoVolume("Flex2Vol5cm",flex2_5cm,medKapton);
  TGeoVolume *bus1Vol      = new TGeoVolume("Bus1Vol",bus1,medAluminum);
  TGeoVolume *bus2Vol      = new TGeoVolume("Bus2Vol",bus2,medKapton);
  TGeoVolume *cone1Vol     = new TGeoVolume("CarbonFleecePipeCover",cone1,medCarbonFleece);
  TGeoVolume *cone2Vol     = new TGeoVolume("GraphitePaperPipeCover",cone2,medFGS003);
  TGeoVolume *plate11Vol   = new TGeoVolume("CarbonFleeceLR1",box11,medCarbonFleece);
  TGeoVolume *plate12Vol   = new TGeoVolume("CarbonFleeceMiddle1",box12,medCarbonFleece); 
  TGeoVolume *plate13Vol   = new TGeoVolume("CarbonFleeceVertical1",box13,medCarbonFleece); 
  TGeoVolume *plate21Vol   = new TGeoVolume("CarbonFleeceLR2",box21,medFGS003);
  TGeoVolume *plate22Vol   = new TGeoVolume("CarbonFleeceMiddle2",box22,medFGS003);
  TGeoVolume *plate23Vol   = new TGeoVolume("CarbonFleeceVertical2",box23,medFGS003);
  TGeoVolume *mechLaddVol  = new TGeoVolume("mechLadderVolume",mechStruct,medAir);

  mechLaddVol->SetLineColor(12);
  mechLaddVol->SetFillColor(12); 
  mechLaddVol->SetVisibility(kTRUE);

  modVol->SetVisibility(kTRUE);
  flex1_5cmVol->SetLineColor(kRed);
  flex2_5cmVol->SetLineColor(kGreen);
  modPlateVol->SetLineColor(kMagenta-8);
  coolTubeVol->SetLineColor(kGray);
  coolTubeWVol->SetLineColor(kBlue);
  coldPlateVol->SetLineColor(kYellow-3);
  glueVol->SetLineColor(kBlack);
  flex1Vol->SetLineColor(kRed);
  flex2Vol->SetLineColor(kGreen);
  bus1Vol->SetLineColor(kCyan);
  bus2Vol->SetLineColor(kBlue);
  cone1Vol->SetFillColor(kViolet);
  plate11Vol->SetFillColor(kViolet);
  plate12Vol->SetLineColor(kViolet);
  plate13Vol->SetLineColor(kViolet);
  cone2Vol->SetLineColor(kGreen);
  plate22Vol->SetFillColor(kGreen);
  plate21Vol->SetLineColor(kGreen);
  plate23Vol->SetLineColor(kGreen);


  //Carbon Fleece

  mechLaddVol->AddNode(plate11Vol,1,new TGeoTranslation(xPos -(1.11/2+rMax+box11->GetDX()+kLay2),-ylen + yPos +2*rMax-kLay2-kLay1/2,0));
  mechLaddVol->AddNode(plate11Vol,2,new TGeoTranslation(xPos +(1.11/2+rMax+box11->GetDX()+kLay2),-ylen + yPos +2*rMax-kLay2-kLay1/2,0));
  mechLaddVol->AddNode(plate11Vol,3,new TGeoTranslation(-xPos -(1.11/2+rMax+box11->GetDX()+kLay2),-ylen + yPos +2*rMax-kLay2-kLay1/2 +deltaY,0));
  mechLaddVol->AddNode(plate11Vol,4,new TGeoTranslation(-xPos +(1.11/2+rMax+box11->GetDX()+kLay2),-ylen + yPos +2*rMax-kLay2-kLay1/2 +deltaY,0));
  mechLaddVol->AddNode(plate12Vol,1,new TGeoTranslation(xPos ,-ylen + yPos +2*rMax-kLay2-kLay1/2,0));
  mechLaddVol->AddNode(plate12Vol,2,new TGeoTranslation(-xPos ,-ylen + yPos +2*rMax-kLay2-kLay1/2 + deltaY,0));
  mechLaddVol->AddNode(plate13Vol,1,new TGeoTranslation(xPos -(1.11/2+rMax+kLay2+box13->GetDX()),-ylen + yPos +2*rMax-kLay1-kLay2-box13->GetDY(),0));
  mechLaddVol->AddNode(plate13Vol,2,new TGeoTranslation(xPos -1.11/2+rMax+kLay2+box13->GetDX(),-ylen + yPos +2*rMax-kLay1-kLay2-box13->GetDY(),0));
  mechLaddVol->AddNode(plate13Vol,3,new TGeoTranslation(xPos +(1.11/2+rMax+kLay2+box13->GetDX()),-ylen + yPos +2*rMax-kLay1-kLay2-box13->GetDY(),0));
  mechLaddVol->AddNode(plate13Vol,4,new TGeoTranslation(xPos +1.11/2-rMax-kLay2-box13->GetDX(),-ylen + yPos +2*rMax-kLay1-kLay2-box13->GetDY(),0));
  mechLaddVol->AddNode(plate13Vol,5,new TGeoTranslation(-xPos -(1.11/2+rMax+kLay2+box13->GetDX()),-ylen + yPos +2*rMax-kLay1-kLay2-box13->GetDY() +deltaY,0));
  mechLaddVol->AddNode(plate13Vol,6,new TGeoTranslation(-xPos -1.11/2+rMax+kLay2+box13->GetDX(),-ylen + yPos +2*rMax-kLay1-kLay2-box13->GetDY() +deltaY,0));
  mechLaddVol->AddNode(plate13Vol,7,new TGeoTranslation(-xPos +(1.11/2+rMax+kLay2+box13->GetDX()),-ylen + yPos +2*rMax-kLay1-kLay2-box13->GetDY() +deltaY,0));
  mechLaddVol->AddNode(plate13Vol,8,new TGeoTranslation(-xPos +1.11/2-rMax-kLay2-box13->GetDX(),-ylen + yPos +2*rMax-kLay1-kLay2-box13->GetDY() +deltaY,0));

  mechLaddVol->AddNode(cone1Vol,1,new TGeoTranslation(xPos - 0.555,-ylen + yPos + rMax,0));
  mechLaddVol->AddNode(cone1Vol,2,new TGeoTranslation(xPos + 0.555,-ylen + yPos + rMax,0));
  mechLaddVol->AddNode(cone1Vol,3,new TGeoTranslation(-xPos - 0.555,-ylen + yPos + rMax + deltaY,0));
  mechLaddVol->AddNode(cone1Vol,4,new TGeoTranslation(-xPos + 0.555,-ylen + yPos + rMax + deltaY,0));


  //Carbon Paper

  mechLaddVol->AddNode(plate21Vol,1,new TGeoTranslation(xPos -(1.11/2+rMax+box21->GetDX()),-ylen + yPos +2*rMax-kLay2/2,0));
  mechLaddVol->AddNode(plate21Vol,2,new TGeoTranslation(xPos +(1.11/2+rMax+box21->GetDX()) ,-ylen + yPos +2*rMax-kLay2/2,0));
  mechLaddVol->AddNode(plate21Vol,3,new TGeoTranslation(-xPos -(1.11/2+rMax+box21->GetDX()) ,-ylen + yPos +2*rMax-kLay2/2 +deltaY,0));
  mechLaddVol->AddNode(plate21Vol,4,new TGeoTranslation(-xPos +(1.11/2+rMax+box21->GetDX()) ,-ylen + yPos +2*rMax-kLay2/2 +deltaY,0));
  mechLaddVol->AddNode(plate22Vol,1,new TGeoTranslation(xPos ,-ylen + yPos +2*rMax-kLay2/2,0));
  mechLaddVol->AddNode(plate22Vol,2,new TGeoTranslation(-xPos ,-ylen + yPos +2*rMax-kLay2/2 + deltaY,0));
  mechLaddVol->AddNode(plate23Vol,1,new TGeoTranslation(xPos -(1.11/2+rMax+box23->GetDX()),-ylen + yPos +2*rMax-kLay2-(rMax-kLay2)/2,0));
  mechLaddVol->AddNode(plate23Vol,2,new TGeoTranslation(xPos -1.11/2+rMax+box23->GetDX(),-ylen + yPos +2*rMax-kLay2-(rMax-kLay2)/2,0));
  mechLaddVol->AddNode(plate23Vol,3,new TGeoTranslation(xPos +(1.11/2+rMax+box23->GetDX()),-ylen + yPos +2*rMax-kLay2-(rMax-kLay2)/2,0));
  mechLaddVol->AddNode(plate23Vol,4,new TGeoTranslation(xPos +1.11/2-rMax-box23->GetDX(),-ylen + yPos +2*rMax-kLay2-(rMax-kLay2)/2,0));
  mechLaddVol->AddNode(plate23Vol,5,new TGeoTranslation(-xPos -(1.11/2+rMax+box23->GetDX()),-ylen + yPos +2*rMax-kLay2-(rMax-kLay2)/2+deltaY,0));
  mechLaddVol->AddNode(plate23Vol,6,new TGeoTranslation(-xPos -1.11/2+rMax+box23->GetDX(),-ylen + yPos +2*rMax-kLay2-(rMax-kLay2)/2+deltaY,0));
  mechLaddVol->AddNode(plate23Vol,7,new TGeoTranslation(-xPos +(1.11/2+rMax+box23->GetDX()),-ylen + yPos +2*rMax-kLay2-(rMax-kLay2)/2+deltaY,0));
  mechLaddVol->AddNode(plate23Vol,8,new TGeoTranslation(-xPos +1.11/2-rMax-box23->GetDX(),-ylen + yPos +2*rMax-kLay2-(rMax-kLay2)/2+deltaY,0));

  mechLaddVol->AddNode(cone2Vol,1,new TGeoTranslation(xPos - 0.555,-ylen + yPos + rMax,0));
  mechLaddVol->AddNode(cone2Vol,2,new TGeoTranslation(xPos + 0.555,-ylen + yPos + rMax,0));
  mechLaddVol->AddNode(cone2Vol,3,new TGeoTranslation(-xPos - 0.555,-ylen + yPos + rMax + deltaY,0));
  mechLaddVol->AddNode(cone2Vol,4,new TGeoTranslation(-xPos + 0.555,-ylen + yPos + rMax + deltaY,0));

  //Cooling Tubes + water

  mechLaddVol->AddNode(coolTubeVol,1,new TGeoTranslation(xPos - 0.555,-ylen + yPos + rMax,0));
  mechLaddVol->AddNode(coolTubeWVol,1,new TGeoTranslation(xPos - 0.555,-ylen + yPos + rMax,0));
  mechLaddVol->AddNode(coolTubeVol,2,new TGeoTranslation(xPos + 0.555,-ylen + yPos + rMax,0));
  mechLaddVol->AddNode(coolTubeWVol,2,new TGeoTranslation(xPos + 0.555,-ylen + yPos + rMax,0));
  mechLaddVol->AddNode(coolTubeVol,3,new TGeoTranslation(-xPos - 0.555,-ylen + yPos + rMax + deltaY,0));
  mechLaddVol->AddNode(coolTubeWVol,3,new TGeoTranslation(-xPos - 0.555,-ylen + yPos + rMax + deltaY,0));
  mechLaddVol->AddNode(coolTubeVol,4,new TGeoTranslation(-xPos + 0.555,-ylen + yPos + rMax + deltaY,0));
  mechLaddVol->AddNode(coolTubeWVol,4,new TGeoTranslation(-xPos + 0.555,-ylen + yPos + rMax + deltaY,0));

  //Cold Plate

  mechLaddVol->AddNode(coldPlateVol,1,new TGeoTranslation(xPos,-ylen + yPos + 2*rMax + yCPlate/2,0));
  mechLaddVol->AddNode(coldPlateVol,2,new TGeoTranslation(-xPos,-ylen + yPos + 2*rMax + yCPlate/2 + deltaY,0));

  //Glue

  mechLaddVol->AddNode(glueVol,1,new TGeoTranslation(xPos,-ylen + yPos + 2*rMax + yCPlate + yGlue/2,0));
  mechLaddVol->AddNode(glueVol,2,new TGeoTranslation(-xPos,-ylen + yPos + 2*rMax + yCPlate + yGlue/2 + deltaY,0));

  //Module Carbon Plate

  mechLaddVol->AddNode(modPlateVol,1,new TGeoTranslation(xPos, -ylen + yPos + 2*rMax + yCPlate + yGlue + yModPlate/2,0));
  mechLaddVol->AddNode(modPlateVol,2,new TGeoTranslation(-xPos, -ylen + yPos + 2*rMax + yCPlate + yGlue + yModPlate/2 + deltaY,0));

  //Bus

  mechLaddVol->AddNode(bus1Vol,1,new TGeoTranslation(xPos, -ylen + yPos + 2*rMax + yCPlate + yGlue + yModPlate + 2*ymod + 2*yFlex1 + 2*yFlex2 + yBus1/2,0));
  mechLaddVol->AddNode(bus1Vol,2,new TGeoTranslation(-xPos, -ylen + yPos + 2*rMax + yCPlate + yGlue + yModPlate + 2*ymod + 2*yFlex1 + 2*yFlex2 + yBus1/2 + deltaY,0));
  mechLaddVol->AddNode(bus2Vol,1,new TGeoTranslation(xPos, -ylen + yPos + 2*rMax + yCPlate + yGlue + yModPlate + 2*ymod + 2*yFlex1 + 2*yFlex2 + yBus1 + yBus2/2,0));
  mechLaddVol->AddNode(bus2Vol,2,new TGeoTranslation(-xPos, -ylen + yPos + 2*rMax + yCPlate + yGlue + yModPlate + 2*ymod + 2*yFlex1 + 2*yFlex2 + yBus1 + yBus2/2 + deltaY,0));
  
  //FPC + modules

  for (Int_t j=0; j<fNModules; j++) {

    zpos = -(zact + (fNModules-1)*modGap)/2 + j*(zMod + modGap) + zMod/2;
    zpos5cm = -(zact + (fNModules-1)*modGap)/2 + (j+1)*(zMod + modGap) + flexOverlap/2 ;

    mechLaddVol->AddNode(modVol, j, new TGeoTranslation(xPos, -ylen + yPos + 2*rMax + yCPlate + yGlue + yModPlate + ymod, zpos));
    mechLaddVol->AddNode(modVol, fNModules+j, new TGeoTranslation(-xPos, -ylen + yPos + 2*rMax + yCPlate + yGlue + yModPlate + ymod +deltaY, zpos));
    mechLaddVol->AddNode(flex1Vol,j,new TGeoTranslation(xPos, -ylen + yPos + 2*rMax + yCPlate + yGlue + yModPlate + 2*ymod + yFlex1/2,zpos));
    mechLaddVol->AddNode(flex1Vol,fNModules+j,new TGeoTranslation(-xPos, -ylen + yPos + 2*rMax + yCPlate + yGlue + yModPlate + 2*ymod + yFlex1/2 + deltaY,zpos));
    mechLaddVol->AddNode(flex2Vol,j,new TGeoTranslation(xPos, -ylen + yPos + yModPlate + 2*rMax + yCPlate + yGlue + 2*ymod + yFlex1 + yFlex2/2,zpos));
    mechLaddVol->AddNode(flex2Vol,fNModules+j,new TGeoTranslation(-xPos, -ylen + yPos + 2*rMax + yCPlate + yGlue + yModPlate + 2*ymod + yFlex1 + yFlex2/2 + deltaY,zpos));

    if((j+1)!=fNModules){
      mechLaddVol->AddNode(flex1_5cmVol,j,new TGeoTranslation(xPos,-ylen + yPos + 2*rMax + yCPlate + yGlue + yModPlate + 2*ymod + yFlex1 + yFlex2 + yFlex1/2,zpos5cm));
      mechLaddVol->AddNode(flex1_5cmVol,fNModules+j,new TGeoTranslation(-xPos,-ylen + yPos + 2*rMax + yCPlate + yGlue + yModPlate + 2*ymod + yFlex1 + yFlex2 + yFlex1/2 +deltaY,zpos5cm));
      mechLaddVol->AddNode(flex2_5cmVol,j,new TGeoTranslation(xPos,-ylen + yPos + 2*rMax + yCPlate + yGlue + yModPlate + 2*ymod + 2*yFlex1 + 3*yFlex2/2,zpos5cm));
      mechLaddVol->AddNode(flex2_5cmVol,fNModules+j,new TGeoTranslation(-xPos,-ylen + yPos + 2*rMax + yCPlate + yGlue + yModPlate + 2*ymod + 2*yFlex1 + 3*yFlex2/2 +deltaY,zpos5cm));
    }
    else {
      mechLaddVol->AddNode(flex1_5cmVol,j,new TGeoTranslation(xPos,-ylen + yPos + 2*rMax + yCPlate + yGlue + yModPlate + 2*ymod + yFlex1/2,zpos5cm-modGap));
      mechLaddVol->AddNode(flex1_5cmVol,fNModules+j,new TGeoTranslation(-xPos,-ylen + yPos + 2*rMax + yCPlate + yGlue + yModPlate + 2*ymod + yFlex1/2 +deltaY,zpos5cm-modGap));
      mechLaddVol->AddNode(flex2_5cmVol,j,new TGeoTranslation(xPos,-ylen + yPos + 2*rMax + yCPlate + yGlue + yModPlate +2*ymod + yFlex1 + yFlex2/2,zpos5cm-modGap));
      mechLaddVol->AddNode(flex2_5cmVol,fNModules+j,new TGeoTranslation(-xPos,-ylen + yPos + 2*rMax + yCPlate + yGlue + yModPlate + 2*ymod + yFlex1 + yFlex2/2 +deltaY,zpos5cm-modGap));

      }
  }
  

  // Done, return the stave structur
  return mechLaddVol;
}

//________________________________________________________________________
TGeoVolume* AliITSUv1Layer::CreateSpaceFrameOuterB(const Double_t xlad,
						   const TGeoManager *mgr){
//
// Create the space frame for the Outer Barrel
//
// Input:
//         xlad : X length
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
//

  TGeoVolume *mechLaddVol = 0;

  switch (fStaveModel) {
    case AliITSUv1::kOBModelDummy:
      mechLaddVol = CreateSpaceFrameOuterBDummy(xlad,mgr);
      break;
    case AliITSUv1::kOBModel1:
      mechLaddVol = CreateSpaceFrameOuterB0(xlad,mgr);
      break;
    default:
      AliFatal(Form("Unknown stave model %d",fStaveModel));
      break;
  }

  return mechLaddVol; 
}

//________________________________________________________________________
TGeoVolume* AliITSUv1Layer::CreateSpaceFrameOuterBDummy(const Double_t ,
							const TGeoManager *) const {
//
// Create dummy stave
//
// Input:
//         xlad : X length
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
TGeoVolume* AliITSUv1Layer::CreateSpaceFrameOuterB0(const Double_t xlen,
						    const TGeoManager *mgr){
//
// Create the space frame for the Outer Barrel (Model 0)
//
// Input:
//         xlen : X length
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//         a TGeoVolume with the Space Frame of a stave
//
// Created:      20 Dec 2013  Anastasia Barbano
// Updated:      15 Jan 2014  Mario Sitta
//


  // Materials defined in AliITSUv0
  TGeoMedium *medCarbon       = mgr->GetMedium("ITS_CARBON$");
  TGeoMedium *medAir          = mgr->GetMedium("ITS_AIR$");


  // Local parameters
  Double_t ladderWidth        = 4.2;
  Double_t ladderHeight       = 4.2;
//  Double_t ladderSegBoxDW     = 7.5;
//  Double_t ladderSegBoxDH     = 7.1;
  Double_t ladderBeamRadius   = 0.06;
  Double_t ladderLa           = 0.3;
  Double_t ladderHa           = 0.0721979;
  Double_t ladderLb           = 0.37;
  Double_t ladderHb           = 0.0890428;
  Double_t ladderl            = 0.025;
  Double_t beamSidePhi        = 65;
  Double_t bottomBeamAngle    = 56.5;
//  Double_t dy                 = ladderSegBoxDH/2;
  Double_t triangleHeight     = ladderHeight - ladderBeamRadius;
  Double_t halfTheta          = TMath::ATan( 0.5*ladderWidth/triangleHeight );
  Double_t alpha              = TMath::Pi()*3./4. - halfTheta/2.;
  Double_t beta               = (TMath::Pi() - 2.*halfTheta)/4.;
//  Double_t dYTranslation      = (ladderHeight/2. -0.5*ladderWidth*TMath::Tan(beta)-ladderBeamRadius);
  Double_t distCenterSideDown =  0.5*ladderWidth/TMath::Cos(beta);
  Double_t zact;
  Double_t zbus;
  Double_t zlen;
  Double_t seglen;    


  zact    = fNModules*fgkOBModuleZLength; //active area 
  zbus    = zact + (fNModules-1)*fgkOBModuleGap; 
  zlen    = zbus/2;
  seglen  = zlen/10;

  // First create all needed shapes and volumes

  TGeoBBox *spaceFrame = new TGeoBBox("CarbonFrame",xlen, 2.2, zlen);
  TGeoBBox *segment    = new TGeoBBox(ladderWidth/2,ladderHeight/2,seglen/2);

  TGeoVolume *spaceFrameVol = new TGeoVolume("CarbonFrameVolume",
					     spaceFrame, medAir);
  spaceFrameVol->SetVisibility(kTRUE);

  TGeoVolume *segmentVol    = new TGeoVolume("segmentVol",segment,medAir);

  //SpaceFrame

  //--- the top V of the Carbon Fiber Ladder (segment)
  TGeoArb8 *cfLaddTop1 = CreateLadderSide("CFladdTopCornerVol1shape", seglen/2., halfTheta, -1, ladderLa, ladderHa, ladderl);
  TGeoVolume *cfLaddTopVol1 = new TGeoVolume("ITSsddCFladdTopCornerVol1", cfLaddTop1,medCarbon);
  TGeoArb8 *cfLaddTop2 = CreateLadderSide( "CFladdTopCornerVol2shape", seglen/2., halfTheta, 1, ladderLa, ladderHa, ladderl);
  TGeoVolume *cfLaddTopVol2 = new TGeoVolume("ITSsddCFladdTopCornerVol2",cfLaddTop2,medCarbon );

  //TGeoTranslation *trTop1 = new TGeoTranslation(0, fgkLadderHeight/2-dy, 0);
  TGeoTranslation *trTop1 = new TGeoTranslation(0, ladderHeight/2, 0);
  
  //--- the 2 side V
  TGeoArb8 *cfLaddSide1 = CreateLadderSide( "CFladdSideCornerVol1shape", seglen/2., beta, -1,ladderLb, ladderHb, ladderl);
  TGeoVolume *cfLaddSideVol1 = new TGeoVolume( "ITSsddCFladdSideCornerVol1", cfLaddSide1,medCarbon);
  TGeoArb8 *cfLaddSide2 = CreateLadderSide( "CFladdSideCornerVol2shape", seglen/2., beta, 1, ladderLb, ladderHb, ladderl);
  TGeoVolume *cfLaddSideVol2 = new TGeoVolume( "ITSsddCFladdSideCornerVol2", cfLaddSide2,medCarbon );

  
  TGeoCombiTrans *ctSideR = CreateCombiTrans("", distCenterSideDown, 0,alpha*TMath::RadToDeg());
  //AddTranslationToCombiTrans(ctSideR, 0, -dYTranslation-dy, 0);
  AddTranslationToCombiTrans(ctSideR, 0, ladderHeight/2-2.85/*2.765250/*triangleHeight*/, 0);
  TGeoCombiTrans *ctSideL = CreateCombiTrans("", distCenterSideDown,0,-alpha*TMath::RadToDeg());
  //AddTranslationToCombiTrans(ctSideL, 0, -dYTranslation-dy, 0);
  AddTranslationToCombiTrans(ctSideL, 0, ladderHeight/2-2.85/*triangleHeight*/, 0);

  segmentVol->AddNode(cfLaddTopVol1,1,trTop1);
  segmentVol->AddNode(cfLaddTopVol2,1,trTop1);
  segmentVol->AddNode(cfLaddSideVol1,1,ctSideR);
  segmentVol->AddNode(cfLaddSideVol1,2,ctSideL);
  segmentVol->AddNode(cfLaddSideVol2,1,ctSideR);
  segmentVol->AddNode(cfLaddSideVol2,2,ctSideL);


  //--- The beams
  // Beams on the sides
  Double_t beamPhiPrime = TMath::ASin(1./TMath::Sqrt( (1+TMath::Sin(2*beta)*TMath::Sin(2*beta)/(TanD(beamSidePhi)*TanD(beamSidePhi))) ));
  Double_t beamLength = TMath::Sqrt( ladderHeight*ladderHeight/( TMath::Sin(beamPhiPrime)*TMath::Sin(beamPhiPrime))+ ladderWidth*ladderWidth/4.)-ladderLa/2-ladderLb/2;
  TGeoTubeSeg *sideBeamS = new TGeoTubeSeg(0, ladderBeamRadius,beamLength/2.,0, 180);
  TGeoVolume *sideBeam = new TGeoVolume("ITSsddCFSideBeamVol", sideBeamS,medCarbon);

  //Euler rotation : about Z, then new X, then new Z
  TGeoRotation *beamRot1 = new TGeoRotation("", 90-2.*beta*TMath::RadToDeg(),-beamPhiPrime*TMath::RadToDeg(),-90);
  TGeoRotation *beamRot2 = new TGeoRotation("", 90-2.*beta*TMath::RadToDeg(), beamPhiPrime*TMath::RadToDeg(), -90);
  TGeoRotation *beamRot3 = new TGeoRotation("", 90+2.*beta*TMath::RadToDeg(), beamPhiPrime*TMath::RadToDeg(), -90);
  TGeoRotation *beamRot4 = new TGeoRotation("", 90+2.*beta*TMath::RadToDeg(),-beamPhiPrime*TMath::RadToDeg(),-90);

  TGeoCombiTrans *beamTransf[8];
  beamTransf[0] = new TGeoCombiTrans( 0.5*triangleHeight*TMath::Tan(halfTheta),ladderBeamRadius/2. ,-3*seglen/8, beamRot1);

  beamTransf[1] = new TGeoCombiTrans( 0.5*triangleHeight*TMath::Tan(halfTheta),ladderBeamRadius/2. ,-3*seglen/8, beamRot1);
  AddTranslationToCombiTrans(beamTransf[1], 0, 0, seglen/2);

  beamTransf[2] = new TGeoCombiTrans(0.5*triangleHeight*TMath::Tan(halfTheta),ladderBeamRadius/2. ,-seglen/8, beamRot2);

  beamTransf[3] = new TGeoCombiTrans(0.5*triangleHeight*TMath::Tan(halfTheta),ladderBeamRadius/2. ,-seglen/8, beamRot2);
  AddTranslationToCombiTrans(beamTransf[3], 0, 0, seglen/2);

  beamTransf[4] = new TGeoCombiTrans(-0.5*triangleHeight*TMath::Tan(halfTheta),ladderBeamRadius/2. ,-3*seglen/8, beamRot3);

  beamTransf[5] = new TGeoCombiTrans(-0.5*triangleHeight*TMath::Tan(halfTheta),ladderBeamRadius/2. ,-3*seglen/8, beamRot3);
  AddTranslationToCombiTrans(beamTransf[5], 0, 0, seglen/2);

  beamTransf[6] = new TGeoCombiTrans(-0.5*triangleHeight* TMath::Tan(halfTheta),ladderBeamRadius/2., -seglen/8,beamRot4);
  beamTransf[7] = new TGeoCombiTrans(-0.5*triangleHeight* TMath::Tan(halfTheta),ladderBeamRadius/2.,3*seglen/8,beamRot4);

  //--- Beams of the bottom
  TGeoTubeSeg *bottomBeam1 = new TGeoTubeSeg(0, ladderBeamRadius,ladderWidth/2.-ladderLb/3, 0, 180);
  TGeoVolume *bottomBeam1Vol = new TGeoVolume("ITSsddBottomBeam1Vol", bottomBeam1, medCarbon);
  TGeoTubeSeg *bottomBeam2 = new TGeoTubeSeg(0, ladderBeamRadius,ladderWidth/2.-ladderLb/3, 0, 90);
  TGeoVolume *bottomBeam2Vol = new TGeoVolume("ITSsddBottomBeam2Vol",bottomBeam2, medCarbon);
  TGeoTubeSeg *bottomBeam3 = new TGeoTubeSeg(0, ladderBeamRadius,0.5*ladderWidth/SinD(bottomBeamAngle) - ladderLb/3, 0, 180);
  TGeoVolume *bottomBeam3Vol = new TGeoVolume("ITSsddBottomBeam3Vol", bottomBeam3, medCarbon);
  TGeoRotation *bottomBeamRot1 = new TGeoRotation("", 90, 90,  90);
  TGeoRotation *bottomBeamRot2 = new TGeoRotation("",-90, 90, -90);

  TGeoCombiTrans *bottomBeamTransf1 = new TGeoCombiTrans("",0,-(ladderHeight/2-ladderBeamRadius),0, bottomBeamRot1);
  TGeoCombiTrans *bottomBeamTransf2 = new TGeoCombiTrans(0,-(ladderHeight/2-ladderBeamRadius),-seglen/2, bottomBeamRot1);
  TGeoCombiTrans *bottomBeamTransf3 = new TGeoCombiTrans(0,-(ladderHeight/2-ladderBeamRadius), seglen/2, bottomBeamRot2);
  // be careful for beams #3: when "reading" from -z to +z and 
  // from the bottom of the ladder, it should draw a Lambda, and not a V
  TGeoRotation *bottomBeamRot4 = new TGeoRotation("", -90, bottomBeamAngle, -90);
  TGeoRotation *bottomBeamRot5 = new TGeoRotation("" ,-90,-bottomBeamAngle, -90);
  TGeoCombiTrans *bottomBeamTransf4 = new TGeoCombiTrans(0,-(ladderHeight/2-ladderBeamRadius),-seglen/4,bottomBeamRot4);
  TGeoCombiTrans *bottomBeamTransf5 = new TGeoCombiTrans(0,-(ladderHeight/2-ladderBeamRadius),seglen/4, bottomBeamRot5);

  cfLaddTopVol1->SetLineColor(35);
  cfLaddTopVol2->SetLineColor(35);
  cfLaddSideVol1->SetLineColor(35);
  cfLaddSideVol2->SetLineColor(35);
  sideBeam->SetLineColor(35);
  bottomBeam1Vol->SetLineColor(35);
  bottomBeam2Vol->SetLineColor(35);
  bottomBeam3Vol->SetLineColor(35);


  segmentVol->AddNode(sideBeam,1, beamTransf[0]);
  segmentVol->AddNode(sideBeam,2, beamTransf[1]);
  segmentVol->AddNode(sideBeam,3, beamTransf[2]);
  segmentVol->AddNode(sideBeam,4, beamTransf[3]);
  segmentVol->AddNode(sideBeam,5, beamTransf[4]);
  segmentVol->AddNode(sideBeam,6, beamTransf[5]);
  segmentVol->AddNode(sideBeam,7, beamTransf[6]);
  segmentVol->AddNode(sideBeam,8, beamTransf[7]);
  segmentVol->AddNode(bottomBeam1Vol,1,bottomBeamTransf1);
  segmentVol->AddNode(bottomBeam2Vol,1,bottomBeamTransf2);
  segmentVol->AddNode(bottomBeam2Vol,1,bottomBeamTransf3);
  segmentVol->AddNode(bottomBeam3Vol,1,bottomBeamTransf4);
  segmentVol->AddNode(bottomBeam3Vol,1,bottomBeamTransf5);

  for(Int_t i=0;i<10;i++){
    spaceFrameVol->AddNode(segmentVol,i,new TGeoTranslation(0,0,seglen*(0.5+i)));
    spaceFrameVol->AddNode(segmentVol,11+i,new TGeoTranslation(0,0,-seglen*(0.5+i)));
  }


  // Done, return the stave structur
  return spaceFrameVol;
}

//________________________________________________________________________
TGeoVolume* AliITSUv1Layer::CreateModuleInnerB(const Double_t xlad,
					       const Double_t ylad,   
					       const Double_t zlad,
					       const TGeoManager *mgr){
//
// Creates the actual Module
//
// Input:
//         xlad,zlad : the ladder dimensions
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

  // The module
  TGeoBBox *module = new TGeoBBox(xlad,  ylad, zlad/fNModules);

  // The sensor
  xlen = module->GetDX();
  ylen = 0.5*fSensorThick;
  zlen = module->GetDZ();
  TGeoBBox *sensor = new TGeoBBox(xlen, ylen, zlen);


  // We have all shapes: now create the real volumes
  //TGeoMedium *medAir = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medSi  = mgr->GetMedium("ITS_SI$");

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSModulePattern(), fLayerNumber);
  // TGeoVolume *modVol = new TGeoVolume(volname, module, medAir);
  TGeoVolume *modVol = new TGeoVolume(volname, module, medSi);
  modVol->SetVisibility(kTRUE);
  modVol->SetLineColor(1);

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSSensorPattern(), fLayerNumber);
  TGeoVolume *sensVol = new TGeoVolume(volname, sensor, medSi);
  sensVol->SetVisibility(kTRUE);
  sensVol->SetLineColor(8);
  sensVol->SetLineWidth(1);
  sensVol->SetFillColor(sensVol->GetLineColor());
  sensVol->SetFillStyle(4000); // 0% transparent


  // Now build up the module
  xpos = 0.;
  ypos = -module->GetDY() + sensor->GetDY();
  zpos = 0.;

  modVol->AddNode(sensVol, 1, new TGeoTranslation(xpos, ypos, zpos));

  // Done, return the module
  return modVol;
}

//________________________________________________________________________
TGeoVolume* AliITSUv1Layer::CreateModuleOuterB(const Double_t xlad,
					       const Double_t ylad,   
					       const Double_t zmod,
					       const TGeoManager *mgr){
//
// Creates the actual Module
//
// Input:
//         xlad,ylad,zlad : the half stave dimensions
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      18 Dec 2013  M. Sitta, A. Barbano
//


  char volname[30];

  Double_t xGap  = 0.01; 
  Double_t zGap  = 0.01;

  Double_t xlen, ylen, zlen;
  Double_t xpos, ypos, zpos;
  
  // First create all needed shapes

  // The module
  TGeoBBox *module = new TGeoBBox(xlad,  ylad, zmod/2);

  // The sensor
  xlen = 0.5*(module->GetDX()-xGap/2);
  //xlen = 0.5*1.5;
  ylen = ylad;
  zlen = (2*module->GetDZ()-6*zGap)/14;
  TGeoBBox *sensor = new TGeoBBox(xlen, ylen, zlen);


  // We have all shapes: now create the real volumes
 
  TGeoMedium *medSi  = mgr->GetMedium("ITS_SI$");
  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSModulePattern(), fLayerNumber);
  TGeoVolume *modVol = new TGeoVolume(volname, module, medSi);
  modVol->SetVisibility(kTRUE);

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSSensorPattern(), fLayerNumber);
  TGeoVolume *sensVol = new TGeoVolume(volname, sensor, medSi);
  


  // Now build up the module
  xpos = -module->GetDX() + sensor->GetDX();
  //xpos = -xGap/2 -sensor->GetDX();
  ypos = -module->GetDY() + sensor->GetDY();
  for(Int_t k=0;k<7;k++)   //put 7x2 chip into one module
    {
      zpos = -module->GetDZ() + sensor->GetDZ() + k*(2*sensor->GetDZ() + zGap);
      modVol->AddNode(sensVol, k+1, new TGeoTranslation(xpos, ypos, zpos));
      modVol->AddNode(sensVol, k+2, new TGeoTranslation(-xpos, ypos, zpos));
    }
  
  //sensVol->SetVisibility(kTRUE);
  sensVol->SetLineColor(kYellow);
  //sensVol->SetLineWidth(1);
  //sensVol->SetTransparency(30);
  sensVol->SetFillColor(sensVol->GetLineColor());
  sensVol->SetFillStyle(4000); // 0% transparent
  // Done, return the module
  return modVol;
}

//________________________________________________________________________
Double_t AliITSUv1Layer::RadiusOfTurboContainer(){
//
// Computes the inner radius of the air container for the Turbo configuration
// as the radius of either the circle tangent to the ladder or the circle
// passing for the ladder's lower vertex
//
// Input:
//         none (all needed parameters are class members)
//
// Output:
//
// Return:
//        the radius of the container if >0, else flag to use the lower vertex
//
// Created:      08 Mar 2012  Mario Sitta
//

  Double_t rr, delta, z, lladd, rladd;

  if (fLadderThick > 89.) // Very big angle: avoid overflows since surely
    return -1;            // the radius from lower vertex is the right value

  rladd = fLayRadius + 0.5*fLadderThick;
  delta = (0.5*fLadderThick)/CosD(fLadderTilt);
  z     = (0.5*fLadderThick)*TanD(fLadderTilt);

  rr = rladd - delta;
  lladd = (0.5*fLadderWidth) - z;

  if ( (rr*SinD(fLadderTilt) < lladd) )
    return (rr*CosD(fLadderTilt));
  else
    return -1;
}

//________________________________________________________________________
void AliITSUv1Layer::SetLadderTilt(const Double_t t)
{
//
// Sets the Ladder tilt angle (for turbo layers only)
//
// Input:
//         t :  the ladder tilt angle
//
// Output:
//
// Return:
//
// Created:      08 Jul 2011  Mario Sitta
//

  if (fIsTurbo)
    fLadderTilt = t;
  else
    AliError("Not a Turbo layer");

}

//________________________________________________________________________
void AliITSUv1Layer::SetLadderWidth(const Double_t w){
//
// Sets the Ladder width (for turbo layers only)
//
// Input:
//         w :  the ladder width
//
// Output:
//
// Return:
//
// Created:      08 Jul 2011  Mario Sitta
//

  if (fIsTurbo)
    fLadderWidth = w;
  else
    AliError("Not a Turbo layer");

}

//________________________________________________________________________
TGeoArb8 *AliITSUv1Layer::CreateLadderSide(const char *name,
                         Double_t dz, Double_t angle, Double_t xSign,
                         Double_t L, Double_t H, Double_t l) {
//
// Creates the V-shaped sides of the OB space frame
// (from a similar method with same name and function
// in AliITSv11GeometrySDD class by L.Gaudichet)
//

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
TGeoCombiTrans *AliITSUv1Layer::CreateCombiTrans(const char *name,
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
void AliITSUv1Layer::AddTranslationToCombiTrans(TGeoCombiTrans* ct,
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
