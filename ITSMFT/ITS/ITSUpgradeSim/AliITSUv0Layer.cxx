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


/* $Id: AliITSUv0Layer.cxx  */
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
#include "AliITSUv0Layer.h"
#include "AliITSUGeomTGeo.h"
#include <TGeoBBox.h>
#include <TGeoShape.h>
#include <TGeoTrd1.h>

#include <AliLog.h>

using namespace TMath;

const Double_t AliITSUv0Layer::fgkDefaultSensorThick = 300*1.0E-4;
const Double_t AliITSUv0Layer::fgkDefaultStaveThick =   1.;

ClassImp(AliITSUv0Layer)

#define SQ(A) (A)*(A)

//________________________________________________________________________
AliITSUv0Layer::AliITSUv0Layer(): 
  TObject(),
  fLayerNumber(0),
  fPhi0(0),
  fLayRadius(0),
  fZLength(0),
  fSensorThick(0),
  fStaveThick(0),
  fStaveWidth(0),
  fStaveTilt(0),
  fNStaves(0),
  fNChips(0),
  fChipTypeID(0),
  fIsTurbo(0),
  fBuildLevel(0),
  fStaveModel(AliITSUv0::kModelDummy)
{
  //
  // Standard constructor
  //
}

//________________________________________________________________________
AliITSUv0Layer::AliITSUv0Layer(Int_t debug): 
  TObject(),
  fLayerNumber(0),
  fPhi0(0),
  fLayRadius(0),
  fZLength(0),
  fSensorThick(0),
  fStaveThick(0),
  fStaveWidth(0),
  fStaveTilt(0),
  fNStaves(0),
  fNChips(0),
  fChipTypeID(0),
  fIsTurbo(0),
  fBuildLevel(0),
  fStaveModel(AliITSUv0::kModelDummy)
{
  //
  // Constructor setting debugging level
  //
}

//________________________________________________________________________
AliITSUv0Layer::AliITSUv0Layer(Int_t lay, Int_t debug): 
  TObject(),
  fLayerNumber(lay),
  fPhi0(0),
  fLayRadius(0),
  fZLength(0),
  fSensorThick(0),
  fStaveThick(0),
  fStaveWidth(0),
  fStaveTilt(0),
  fNStaves(0),
  fNChips(0),
  fChipTypeID(0),
  fIsTurbo(0),
  fBuildLevel(0),
  fStaveModel(AliITSUv0::kModelDummy)
{
  //
  // Constructor setting layer number and debugging level
  //
}

//________________________________________________________________________
AliITSUv0Layer::AliITSUv0Layer(Int_t lay, Bool_t turbo, Int_t debug): 
  TObject(),
  fLayerNumber(lay),
  fPhi0(0),
  fLayRadius(0),
  fZLength(0),
  fSensorThick(0),
  fStaveThick(0),
  fStaveWidth(0),
  fStaveTilt(0),
  fNStaves(0),
  fNChips(0),
  fChipTypeID(0),
  fIsTurbo(turbo),
  fBuildLevel(0),
  fStaveModel(AliITSUv0::kModelDummy)
{
  //
  // Constructor setting layer number and debugging level
  // for a "turbo" layer (i.e. where staves overlap in phi)
  //
}

//________________________________________________________________________
AliITSUv0Layer::AliITSUv0Layer(const AliITSUv0Layer &s):
  TObject(),
  fLayerNumber(s.fLayerNumber),
  fPhi0(s.fPhi0),
  fLayRadius(s.fLayRadius),
  fZLength(s.fZLength),
  fSensorThick(s.fSensorThick),
  fStaveThick(s.fStaveThick),
  fStaveWidth(s.fStaveWidth),
  fStaveTilt(s.fStaveTilt),
  fNStaves(s.fNStaves),
  fNChips(s.fNChips),
  fChipTypeID(s.fChipTypeID),
  fIsTurbo(s.fIsTurbo),
  fBuildLevel(s.fBuildLevel),
  fStaveModel(s.fStaveModel)
{
  //
  // Copy constructor
  //
}

//________________________________________________________________________
AliITSUv0Layer& AliITSUv0Layer::operator=(const AliITSUv0Layer &s)
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
  fStaveThick = s.fStaveThick;
  fStaveWidth = s.fStaveWidth;
  fStaveTilt  = s.fStaveTilt;
  fNStaves    = s.fNStaves;
  fNChips    = s.fNChips;
  fIsTurbo     = s.fIsTurbo;
  fChipTypeID   = s.fChipTypeID;
  fBuildLevel  = s.fBuildLevel;
  fStaveModel  = s.fStaveModel;

  return *this;
}

//________________________________________________________________________
AliITSUv0Layer::~AliITSUv0Layer() {
  //
  // Destructor
  //
}

//________________________________________________________________________
void AliITSUv0Layer::CreateLayer(TGeoVolume *moth){
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
  if (fNStaves  <= 0) AliFatal(Form("Wrong number of staves (%d)",fNStaves));
  if (fNChips  <= 0) AliFatal(Form("Wrong number of chips (%d)",fNChips));

  if (fStaveThick <= 0) {
    AliInfo(Form("Stave thickness wrong or not set (%f), using default (%f)",
		 fStaveThick,fgkDefaultStaveThick));
    fStaveThick = fgkDefaultStaveThick;
  }

  if (fSensorThick <= 0) {
    AliInfo(Form("Sensor thickness wrong or not set (%f), using default (%f)",
		 fSensorThick,fgkDefaultSensorThick));
    fSensorThick = fgkDefaultSensorThick;
  }

  if (fSensorThick > fStaveThick) {
    AliWarning(Form("Sensor thickness (%f) is greater than stave thickness (%f), fixing",
		 fSensorThick,fStaveThick));
    fSensorThick = fStaveThick;
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

  layVol->GetShape()->ComputeBBox(); //RS Enforce recomputing BBox

  // Finally put everything in the mother volume
  moth->AddNode(layVol, 1, 0);


  // Upgrade geometry is served
  return;
}

//________________________________________________________________________
void AliITSUv0Layer::CreateLayerTurbo(TGeoVolume *moth){
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
  layVol->SetUniqueID(fChipTypeID);
  layVol->SetVisibility(kTRUE);
  layVol->SetLineColor(1);
  TGeoVolume *stavVol = CreateStave();


  // Now build up the layer
  alpha = 360./fNStaves;
  Double_t r = fLayRadius /* +chip thick ?! */;
  for (Int_t j=0; j<fNStaves; j++) {
    Double_t phi = j*alpha + fPhi0;
    xpos = r*CosD(phi);// r*SinD(-phi);
    ypos = r*SinD(phi);// r*CosD(-phi);
    zpos = 0.;
    phi += 90;
    layVol->AddNode(stavVol, j, new TGeoCombiTrans( xpos, ypos, zpos,
						    new TGeoRotation("", phi-fStaveTilt,0,0)));
  }

  layVol->GetShape()->ComputeBBox(); //RS Enforce recomputing BBox

  // Finally put everything in the mother volume
  moth->AddNode(layVol, 1, 0);

  return;
}

//________________________________________________________________________
TGeoVolume* AliITSUv0Layer::CreateStave(const TGeoManager * /*mgr*/){
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
//

  char volname[30];
  Double_t xlen, ylen, zlen;
  Double_t xpos, ypos, zpos, zmod;
  Double_t alpha;


  // First create all needed shapes
  alpha = (360./(2*fNStaves))*DegToRad();

  // The stave
  xlen = fLayRadius*Tan(alpha);
  if (fIsTurbo) xlen = 0.5*fStaveWidth;
  ylen = 0.5*fStaveThick;
  zlen = 0.5*fZLength;

  Double_t yplus = 0.46;
  //  Double_t origin[]={0,-yplus/2,0};
  //  TGeoBBox *stave = new TGeoBBox(xlen, ylen+yplus/2, zlen, origin);
  TGeoXtru *stave = new TGeoXtru(2); //z sections
  Double_t xv[5] = {xlen,xlen,0,-xlen,-xlen};
  Double_t yv[5] = {ylen+0.09,-0.15,-yplus-fSensorThick,-0.15,ylen+0.09};    
  stave->DefinePolygon(5,xv,yv);
  stave->DefineSection(0,-zlen,0,0,1.);
  stave->DefineSection(1,+zlen,0,0,1.);
  

  // We have all shapes: now create the real volumes
//  TGeoMedium *medAir = mgr->GetMedium("ITS_AIR$");

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSStavePattern(), fLayerNumber);
//  TGeoVolume *stavVol = new TGeoVolume(volname, stave, medAir);
  TGeoVolume *stavVol = new TGeoVolumeAssembly(volname);

  //  stavVol->SetVisibility(kFALSE);
  stavVol->SetVisibility(kTRUE);
  stavVol->SetLineColor(2);
  TGeoVolume *modVol = CreateChip(xlen,ylen, zlen);


  // Now build up the stave
  zmod = ((TGeoBBox*)modVol->GetShape())->GetDZ();
  for (Int_t j=0; j<fNChips; j++) {
    xpos = 0.;
    ypos = 0.021;  // Remove small overlap - M.S: 21may13
    zpos = -stave->GetDZ() + j*2*zmod + zmod;
    stavVol->AddNode(modVol, j, new TGeoTranslation(xpos, ypos, zpos));
  }
 

  // put mechanical stave structure, only inner barrel up to now
  if (fLayerNumber<3) {
    TGeoVolume *mechStavVol = CreateStaveStruct(xlen,zlen); 
    if (mechStavVol)
      stavVol->AddNode(mechStavVol, fNChips, new TGeoCombiTrans(0, -0.15-ylen, 0, new TGeoRotation("",0, 0, 180)));
  }

  stavVol->GetShape()->ComputeBBox(); //RS Enforce recomputing BBox  

  // Done, return the stave
  return stavVol;
}

//________________________________________________________________________
TGeoVolume* AliITSUv0Layer::CreateStaveStruct(const Double_t xsta,
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
    case AliITSUv0::kModelDummy:
      mechStavVol = CreateStaveModelDummy(xsta,zsta,mgr);
      break;
    case AliITSUv0::kModel0:
      mechStavVol = CreateStaveModel0(xsta,zsta,mgr);
      break;
    case AliITSUv0::kModel1:
      mechStavVol = CreateStaveModel1(xsta,zsta,mgr);
      break;
    case AliITSUv0::kModel21:
      mechStavVol = CreateStaveModel21(xsta,zsta,mgr);
      break;
    case AliITSUv0::kModel22:
      mechStavVol = CreateStaveModel22(xsta,zsta,mgr);
      break;
    case AliITSUv0::kModel3:
      mechStavVol = CreateStaveModel3(xsta,zsta,mgr);
      break;
    default:
      AliFatal(Form("Unknown stave model %d",fStaveModel));
      break;
  }

  return mechStavVol; 
}


//________________________________________________________________________
TGeoVolume* AliITSUv0Layer::CreateStaveModelDummy(const Double_t ,
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
TGeoVolume* AliITSUv0Layer::CreateStaveModel0(const Double_t xsta,
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
  
  // Materials defined in AliITSUv0
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
TGeoVolume* AliITSUv0Layer::CreateStaveModel1(const Double_t xsta,
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
  
  // Materials defined in AliITSUv0
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
TGeoVolume* AliITSUv0Layer::CreateStaveModel21(const Double_t xsta,
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
  
  // Materials defined in AliITSUv0
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
TGeoVolume* AliITSUv0Layer::CreateStaveModel22(const Double_t xsta,
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
  
  // Materials defined in AliITSUv0
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
    Double_t yv[6] = {-(kConeOutRadius*2)-0.06395,0-0.02,kStaveHeight+0.01,kStaveHeight+0.01,0-0.02,-(kConeOutRadius*2)-0.06395};  // (kConeOutRadius*2)-0.064
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
      mechStavVol->AddNode(ibdv2,1,new TGeoCombiTrans(x+kStaveWidth/2-0.06,y-0.0355,z,new TGeoRotation("ibdv2",-103.3,90,0))); //x-kStaveWidth/2-0.09 old Config.C y-0.0355,

      //right trd3
      TGeoTrd1 *trd3 = new TGeoTrd1(0,kSideVertexMWidth/2,kStaveLength, kSideVertexHeight/2);
      TGeoVolume *ibdv3 = new TGeoVolume("RightVertex",trd3,medM60J3K);
      ibdv3->SetFillColor(12);
      ibdv3->SetLineColor(12);
      mechStavVol->AddNode(ibdv3,1,new TGeoCombiTrans(x-kStaveWidth/2+0.06,y-0.0355,z,new TGeoRotation("ibdv3",103.3,90,0))); //x-kStaveWidth/2+0.09 old Config.C
      
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
      mechStavVol->AddNode(volCable, 0, new TGeoCombiTrans(x, y-(kConeOutRadius+klay3+klay2+klay4+fSensorThick+(klay5)/2)+0.01185, z, new TGeoRotation("",0, 0, 0)));
      }
    // Done, return the stave structe
    return mechStavVol;
}

// model3
//________________________________________________________________________
TGeoVolume* AliITSUv0Layer::CreateStaveModel3(const Double_t xsta,
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
  
  // Materials defined in AliITSUv0
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
    Int_t  loop = (kStaveLength/(2*l))/2;
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

//________________________________________________________________________
TGeoVolume* AliITSUv0Layer::CreateChip(const Double_t xsta,
                                         const Double_t ysta,   
					 const Double_t zsta,
					 const TGeoManager *mgr){
//
// Creates the actual Chip
//
// Input:
//         xsta,zsta : the stave dimensions
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
  TGeoBBox *chip = new TGeoBBox(xsta,  ysta, zsta/fNChips);

  // The sensor
  xlen = chip->GetDX();
  ylen = 0.5*fSensorThick;
  zlen = chip->GetDZ();
  TGeoBBox *sensor = new TGeoBBox(xlen, ylen, zlen);


  // We have all shapes: now create the real volumes
  //TGeoMedium *medAir = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medSi  = mgr->GetMedium("ITS_SI$");

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSChipPattern(), fLayerNumber);
  // TGeoVolume *modVol = new TGeoVolume(volname, chip, medAir);
  TGeoVolume *modVol = new TGeoVolume(volname, chip, medSi);
  modVol->SetVisibility(kTRUE);
  modVol->SetLineColor(1);

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSSensorPattern(), fLayerNumber);
  TGeoVolume *sensVol = new TGeoVolume(volname, sensor, medSi);
  sensVol->SetVisibility(kTRUE);
  sensVol->SetLineColor(8);
  sensVol->SetLineWidth(1);
  sensVol->SetFillColor(sensVol->GetLineColor());
  sensVol->SetFillStyle(4000); // 0% transparent


  // Now build up the chip
  xpos = 0.;
  ypos = -chip->GetDY() + sensor->GetDY();
  zpos = 0.;

  modVol->AddNode(sensVol, 1, new TGeoTranslation(xpos, ypos, zpos));

  // Done, return the chip
  return modVol;
}

//________________________________________________________________________
Double_t AliITSUv0Layer::RadiusOfTurboContainer(){
//
// Computes the inner radius of the air container for the Turbo configuration
// as the radius of either the circle tangent to the stave or the circle
// passing for the stave's lower vertex
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

  Double_t rr, delta, z, lstav, rstav;

  if (fStaveThick > 89.) // Very big angle: avoid overflows since surely
    return -1;            // the radius from lower vertex is the right value

  rstav = fLayRadius + 0.5*fStaveThick;
  delta = (0.5*fStaveThick)/CosD(fStaveTilt);
  z     = (0.5*fStaveThick)*TanD(fStaveTilt);

  rr = rstav - delta;
  lstav = (0.5*fStaveWidth) - z;

  if ( (rr*SinD(fStaveTilt) < lstav) )
    return (rr*CosD(fStaveTilt));
  else
    return -1;
}

//________________________________________________________________________
void AliITSUv0Layer::SetStaveTilt(const Double_t t)
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
void AliITSUv0Layer::SetStaveWidth(const Double_t w){
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
