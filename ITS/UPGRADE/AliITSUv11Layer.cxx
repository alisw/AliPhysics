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
//*************************************************************************


/* $Id: AliITSUv11Layer.cxx  */
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
#include "AliITSUv11Layer.h"
#include "AliITSUGeomTGeo.h"
using namespace TMath;

const Double_t AliITSUv11Layer::fgkDefaultSensorThick = 300*fgkmicron;
const Double_t AliITSUv11Layer::fgkDefaultLadderThick =   1*fgkcm;

ClassImp(AliITSUv11Layer)

#define SQ(A) (A)*(A)

//________________________________________________________________________
AliITSUv11Layer::AliITSUv11Layer(): 
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
  fIsTurbo(0)
{
  //
  // Standard constructor
  //
}

//________________________________________________________________________
AliITSUv11Layer::AliITSUv11Layer(Int_t debug): 
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
  fIsTurbo(0)
{
  //
  // Constructor setting debugging level
  //
}

//________________________________________________________________________
AliITSUv11Layer::AliITSUv11Layer(Int_t lay, Int_t debug): 
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
  fIsTurbo(0)
{
  //
  // Constructor setting layer number and debugging level
  //
}

//________________________________________________________________________
AliITSUv11Layer::AliITSUv11Layer(Int_t lay, Bool_t turbo, Int_t debug): 
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
  fIsTurbo(turbo)
{
  //
  // Constructor setting layer number and debugging level
  // for a "turbo" layer (i.e. where ladders overlap in phi)
  //
}

//________________________________________________________________________
AliITSUv11Layer::AliITSUv11Layer(const AliITSUv11Layer &s):
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
  fDetTypeID(0),
  fIsTurbo(s.fIsTurbo)
{
  //
  // Copy constructor
  //
}

//________________________________________________________________________
AliITSUv11Layer& AliITSUv11Layer::operator=(const AliITSUv11Layer &s)
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
  return *this;
}

//________________________________________________________________________
AliITSUv11Layer::~AliITSUv11Layer() {
  //
  // Destructor
  //
}

//________________________________________________________________________
void AliITSUv11Layer::CreateLayer(TGeoVolume *moth,const TGeoManager *mgr){
//
// Creates the actual Layer and places inside its mother volume
//
// Input:
//         moth : the TGeoVolume owing the volume structure
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      17 Jun 2011  Mario Sitta
// Updated:      08 Jul 2011  Mario Sitta
//
  // Local variables
  char volname[30];
  Double_t rmin, rmax;
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
    CreateLayerTurbo(moth, mgr);
    return;
  }


  // First create the ladder container
  alpha = (360./(2*fNLadders))*DegToRad();

  //  fLadderWidth = fLayRadius*Tan(alpha);

  rmin = 0.98*fLayRadius;
  rmax = 1.02*Sqrt( fLadderWidth*fLadderWidth +
			  (rmin+fLadderThick)*(rmin+fLadderThick) );

  TGeoTube *layer = new TGeoTube(rmin, rmax, 0.5*fZLength);


  // We have all shapes: now create the real volumes
  TGeoMedium *medAir = mgr->GetMedium("ITS_AIR$");

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSLayerPattern(),fLayerNumber);
  TGeoVolume *layVol = new TGeoVolume(volname, layer, medAir);
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
void AliITSUv11Layer::CreateLayerTurbo(TGeoVolume *moth,
					  const TGeoManager *mgr){
//
// Creates the actual Layer and places inside its mother volume
// A so-called "turbo" layer is a layer where ladders overlap in phi
// User can set width and tilt angle, no check is performed here
// to avoid volume overlaps
//
// Input:
//         moth : the TGeoVolume owing the volume structure
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      08 Jul 2011  Mario Sitta
// Updated:      08 Mar 2012  Mario Sitta  Correct way to compute container R
//


  // Local variables
  char volname[30];
  Double_t rmin, rmax, rladd, rcont, d;
  Double_t xpos, ypos, zpos;
  Double_t alpha, gamma;


  // Check if the user set the proper (remaining) parameters
  if (fLadderWidth <= 0)
    AliFatal(Form("Wrong ladder width (%f)",fLadderWidth));
  if (Abs(fLadderTilt) > 45)
    AliWarning(Form("Ladder tilt angle (%f) greater than 45deg",fLadderTilt));


  // First create the ladder container
  // d is half the diagonal of the ladder section
  // rladd is the radius at the ladder's center-of-gravity
  // alpha here is the angle between the diagonal and rladd
  d = 0.5*Sqrt(fLadderThick*fLadderThick + fLadderWidth*fLadderWidth);
  alpha = ACos(0.5*fLadderThick/d)*RadToDeg();
  gamma = alpha - fLadderTilt;
  rladd = fLayRadius + 0.5*fLadderThick;

  // rcont is the radius of the air container
  rcont = RadiusOfTurboContainer();
  
  if (rcont > 0)
    rmin = 0.98*rcont;
  else
    rmin = 0.98*Sqrt( rladd*rladd + d*d - 2*rladd*d*CosD(gamma) );
  
  rmax = 1.02*Sqrt( rladd*rladd + d*d + 2*rladd*d*CosD(gamma) );
  
  TGeoTube *layer = new TGeoTube(rmin, rmax, 0.5*fZLength);


  // We have all shapes: now create the real volumes
  TGeoMedium *medAir = mgr->GetMedium("ITS_AIR$");

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSLayerPattern(), fLayerNumber);
  TGeoVolume *layVol = new TGeoVolume(volname, layer, medAir);
  layVol->SetUniqueID(fDetTypeID);
  layVol->SetVisibility(kTRUE);
  layVol->SetLineColor(1);
  TGeoVolume *laddVol = CreateLadder();


  // Now build up the layer


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
						    new TGeoRotation("", phi-fLadderTilt,0,0)));
  }


  // Finally put everything in the mother volume
  moth->AddNode(layVol, 1, 0);

  return;
}

//________________________________________________________________________
TGeoVolume* AliITSUv11Layer::CreateLadder(const TGeoManager *mgr){
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
//

  char volname[30];
  Double_t xlen, ylen, zlen;
  Double_t xpos, ypos, zpos, zmod;
  Double_t alpha;


  // First create all needed shapes
  alpha = (360./(2*fNLadders))*DegToRad();

  // The ladder
  xlen = fLayRadius*Tan(alpha);
  if (fIsTurbo) xlen = 0.5*fLadderWidth;
  ylen = 0.5*fLadderThick;
  zlen = 0.5*fZLength;

  TGeoBBox *ladder = new TGeoBBox(xlen, ylen, zlen);


  // We have all shapes: now create the real volumes
  TGeoMedium *medAir = mgr->GetMedium("ITS_AIR$");

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSLadderPattern(), fLayerNumber);
  TGeoVolume *laddVol = new TGeoVolume(volname, ladder, medAir);

//  laddVol->SetVisibility(kFALSE);
  laddVol->SetVisibility(kTRUE);
  laddVol->SetLineColor(2);
  TGeoVolume *modVol = CreateModule(ladder->GetDX(), ladder->GetDY(),
				    ladder->GetDZ());


  // Now build up the ladder
  zmod = ((TGeoBBox*)modVol->GetShape())->GetDZ();
  for (Int_t j=0; j<fNModules; j++) {
    xpos = 0.;
    ypos = 0.;
    zpos = -ladder->GetDZ() + j*2*zmod + zmod;
    laddVol->AddNode(modVol, j, new TGeoTranslation(xpos, ypos, zpos));
  }


  // Done, return the ladder
  return laddVol;
}

//________________________________________________________________________
TGeoVolume* AliITSUv11Layer::CreateModule(const Double_t xlad,
						   const Double_t ylad,
						   const Double_t zlad,
						   const TGeoManager *mgr){
//
// Creates the actual Module
//
// Input:
//         xlad,ylad,zlad : the ladder dimensions
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
  TGeoBBox *module = new TGeoBBox(xlad, ylad, zlad/fNModules);

  // The sensor
  xlen = module->GetDX();
  ylen = 0.5*fSensorThick;
  zlen = module->GetDZ();
  TGeoBBox *sensor = new TGeoBBox(xlen, ylen, zlen);


  // We have all shapes: now create the real volumes
  //  TGeoMedium *medAir = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medSi  = mgr->GetMedium("ITS_SI$");

  snprintf(volname, 30, "%s%d", AliITSUGeomTGeo::GetITSModulePattern() ,fLayerNumber);
  //  TGeoVolume *modVol = new TGeoVolume(volname, module, medAir);
  TGeoVolume *modVol = new TGeoVolume(volname, module, medSi);
  modVol->SetVisibility(kFALSE);
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
Double_t AliITSUv11Layer::RadiusOfTurboContainer(){
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
void AliITSUv11Layer::SetLadderTilt(const Double_t t)
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
void AliITSUv11Layer::SetLadderWidth(const Double_t w){
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
