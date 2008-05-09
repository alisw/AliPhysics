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
/** @file    AliFMDGeometryBuilder.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:41:17 2006
    @brief   Class to build the FMD geometry 
*/
//____________________________________________________________________
//                                                                          
// Builder of FMD geometry. 
//
// This class takes care of actually building the geometry using the 
// TGeo classes.  Various parameters are fecthed from the
// AliFMDGeometry manager.  
// Forward Multiplicity Detector based on Silicon wafers. This class
// contains the base procedures for the Forward Multiplicity detector
// Detector consists of 3 sub-detectors FMD1, FMD2, and FMD3, each of
// which has 1 or 2 rings of silicon sensors. 
//                                                       
// 

#include <TArrayD.h>		// ROOT_TArrayD
#include <TGeoManager.h>	// ROOT_TGeoManager
#include <TGeoMatrix.h>	        // ROOT_TGeoMatrix
#include <TGeoTube.h>		// ROOT_TGeoTube
#include <TGeoTrd1.h>		// ROOT_TGeoTrd1
#include <TGeoCone.h>		// ROOT_TGeoTrd1
#include <TGeoVolume.h>		// ROOT_TGeoVolume
#include <TGeoXtru.h>		// ROOT_TGeoXtru
#include <TGeoCompositeShape.h>
#include <TMath.h>
#include <TVector2.h>		// ROOT_TVector2
//#include <TGeoMaterial.h>	// ROOT_TGeoMaterial
//#include <TGeoMedium.h>		// ROOT_TGeoMedium
//#include <TGeoPcon.h>		// ROOT_TGeoPcon
//#include <TGeoPolygon.h>	// ROOT_TGeoPolygon

#include "AliFMDGeometryBuilder.h"	// ALIFMDGEOSIMULATOR_H
#include "AliFMDGeometry.h"	// ALIFMDGEOMETRY_H
#include "AliFMDDetector.h"	// ALIFMDDETECTOR_H
#include "AliFMDRing.h"		// ALIFMDRING_H
#include "AliFMD1.h"		// ALIFMD1_H
#include "AliFMD2.h"		// ALIFMD2_H
#include "AliFMD3.h"		// ALIFMD3_H
// #include "AliFMD.h"		// ALIFMD_H
#include "AliFMDDebug.h"		// ALILOG_H

//====================================================================
ClassImp(AliFMDGeometryBuilder)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
const Char_t* AliFMDGeometryBuilder::fgkActiveName	= "F%cAC";
const Char_t* AliFMDGeometryBuilder::fgkSectorName	= "F%cSC";
const Char_t* AliFMDGeometryBuilder::fgkStripName	= "F%cST";
const Char_t* AliFMDGeometryBuilder::fgkSensorName	= "F%cSE";
const Char_t* AliFMDGeometryBuilder::fgkPCBName	        = "F%cPB";
const Char_t* AliFMDGeometryBuilder::fgkCuName	        = "F%cCU";
const Char_t* AliFMDGeometryBuilder::fgkChipName	= "F%cCH";
const Char_t* AliFMDGeometryBuilder::fgkLongLegName	= "F%cLL";
const Char_t* AliFMDGeometryBuilder::fgkShortLegName	= "F%cSL";
const Char_t* AliFMDGeometryBuilder::fgkFrontVName	= "F%cFH";
const Char_t* AliFMDGeometryBuilder::fgkBackVName	= "F%cBH";
const Char_t* AliFMDGeometryBuilder::fgkRingTopName	= "F%cTV";
const Char_t* AliFMDGeometryBuilder::fgkRingBotName	= "F%cBV";
const Char_t* AliFMDGeometryBuilder::fgkHCName		= "F%dH%c";
const Char_t* AliFMDGeometryBuilder::fgkIHCName		= "F%dI%c";
const Char_t* AliFMDGeometryBuilder::fgkNoseName        = "F3SN";
const Char_t* AliFMDGeometryBuilder::fgkBackName        = "F%dSB";
const Char_t* AliFMDGeometryBuilder::fgkTopName         = "F%dSU";
const Char_t* AliFMDGeometryBuilder::fgkBeamName        = "F%dSL";
const Char_t* AliFMDGeometryBuilder::fgkFlangeName      = "F%dSF";
const Char_t* AliFMDGeometryBuilder::fgkFMDDCuName      = "F%cDC";
const Char_t* AliFMDGeometryBuilder::fgkFMDDPCBName     = "F%cDP";
const Char_t* AliFMDGeometryBuilder::fgkFMDDChipName    = "F%cDI";
const Char_t* AliFMDGeometryBuilder::fgkFMDDName        = "F%cDD";
const Char_t* AliFMDGeometryBuilder::fgkFMDName         = "F%dM%c";

//____________________________________________________________________
AliFMDGeometryBuilder::AliFMDGeometryBuilder() 
  : TTask("FMD", "Geomtry builder"),
    fActiveId(0),
    fDetailed(kTRUE),
    fUseAssembly(kTRUE),
    fSectorOff(0),
    fModuleOff(0),
    fRingOff(0),
    fDetectorOff(0),
    fSi(0),
    fC(0),
    fAl(0),
    fPCB(0),
    fChip(0),
    fAir(0),
    fPlastic(0),
    fCopper(0),
    fSteel(0)
{
  // Default constructor
  fActiveId.Set(2);
}

//____________________________________________________________________
AliFMDGeometryBuilder::AliFMDGeometryBuilder(Bool_t detailed) 
  : TTask("FMD", "Geometry builder"),
    fActiveId(0),
    fDetailed(detailed),
    fUseAssembly(kTRUE),
    fSectorOff(0),
    fModuleOff(0),
    fRingOff(0),
    fDetectorOff(0),
    fSi(0),
    fC(0),
    fAl(0),
    fPCB(0),
    fChip(0),
    fAir(0),
    fPlastic(0),
    fCopper(0),
    fSteel(0)
{
  // Normal constructor
  // 
  // Parameters: 
  // 
  //      fmd		Pointer to AliFMD object 
  //      detailed      Whether to make a detailed simulation or not 
  // 
  fActiveId.Set(2);
}


//____________________________________________________________________
TGeoVolume*
AliFMDGeometryBuilder::RingGeometry(AliFMDRing* r) 
{
  // Setup the geometry of a ring.    The defined TGeoVolume is
  // returned, and should be used when setting up the rest of the
  // volumes. 
  // 
  // 
  // Parameters:
  //
  //     r		Pointer to ring geometry object 
  // 
  // Returns:
  //    pointer to ring volume 
  //
  if (!r) { 
    AliError("Didn't get a ring object");
    return 0;
  }
  Char_t      id       = r->GetId();
  Double_t    siThick  = r->GetSiThickness();
  const Int_t knv      = r->GetNVerticies();
  TVector2*   a        = r->GetVertex(5);
  TVector2*   b        = r->GetVertex(3);
  TVector2*   c        = r->GetVertex(4);
  Double_t    theta    = r->GetTheta();
  Double_t    off      = (TMath::Tan(TMath::Pi() * theta / 180) 
			  * r->GetBondingWidth());
  Double_t    rmax     = b->Mod();
  Double_t    rmin     = r->GetLowR();
  Double_t    pcbThick = r->GetPrintboardThickness();
  Double_t    cuThick  = r->GetCopperThickness();
  Double_t    chipThick= r->GetChipThickness();
  Double_t    modSpace = r->GetModuleSpacing();
  Double_t    legr     = r->GetLegRadius();
  Double_t    legl     = r->GetLegLength();
  Double_t    legoff   = r->GetLegOffset();
  Int_t       ns       = r->GetNStrips();
  Double_t    stripoff = a->Mod();
  Double_t    dstrip   = (rmax - stripoff) / ns;
  Double_t    space    = r->GetSpacing();
  TArrayD xs(knv);
  TArrayD ys(knv);
  for (Int_t i = 0; i < knv; i++) {
    // Reverse the order 
    TVector2* vv = r->GetVertex(knv - 1 - i);
    if (!vv) {
      AliError(Form("Failed to get vertex # %d", knv - 1 - i));
      continue;
    }
    xs[i] = vv->X();
    ys[i] = vv->Y();
  }
  
  // Shape of actual sensor 
  TGeoXtru* sensorShape = new TGeoXtru(2);
  sensorShape->DefinePolygon(knv, xs.fArray, ys.fArray);
  sensorShape->DefineSection(0, - siThick/2);
  sensorShape->DefineSection(1, siThick/2);
  TGeoVolume* sensorVolume = new TGeoVolume(Form(fgkSensorName, id), 
					    sensorShape, fSi);
  sensorVolume->VisibleDaughters(kFALSE);
  Int_t sid = sensorVolume->GetNumber();
  fSectorOff   = -1;
  fModuleOff   = 1;
  fRingOff     = 2;
  fDetectorOff = 3;
  if (fDetailed) {
    fSectorOff   = 1;
    fModuleOff   = 4;
    fRingOff     = 5;
    fDetectorOff = 6;
    // Virtual volume shape to divide - This volume is only defined if
    // the geometry is set to be detailed. 
    TGeoTubeSeg* activeShape = new TGeoTubeSeg(rmin, rmax, siThick/2, 
					       - theta, theta);
    TGeoVolume* activeVolume = new TGeoVolume(Form(fgkActiveName, id),
					      activeShape,fSi);
    TGeoVolume* sectorVolume = activeVolume->Divide(Form(fgkSectorName,id), 
						      2, 2, -theta,0,0,"N");
    TGeoVolume* stripVolume  = sectorVolume->Divide(Form(fgkStripName, id), 
						    1, ns, stripoff, dstrip, 
						    0, "SX");
    sid = stripVolume->GetNumber();
    sensorVolume->AddNodeOverlap(activeVolume, 0);
  }
  
  switch (id) {
  case 'i': case 'I': fActiveId[0] = sid; break;
  case 'o': case 'O': fActiveId[1] = sid; break;
  }

  // Shape of Printed circuit Board 
  for (Int_t i = 0;       i < knv / 2; i++) ys[i] -= off;
  for (Int_t i = knv / 2; i < knv;     i++) ys[i] += off;
  TGeoXtru* pcbShape         = new TGeoXtru(2);
  pcbShape->DefinePolygon(knv, xs.fArray, ys.fArray);
  pcbShape->DefineSection(0, - pcbThick/2);
  pcbShape->DefineSection(1, pcbThick/2);
  TGeoVolume* pcbVolume      = new TGeoVolume(Form(fgkPCBName, id), 
					      pcbShape, fPCB);

  // Copper layer
  TGeoXtru* cuShape       = new TGeoXtru(2);
  cuShape->DefinePolygon(6, xs.fArray, ys.fArray);
  cuShape->DefineSection(0, - cuThick/2);
  cuShape->DefineSection(1, cuThick/2);
  TGeoVolume* cuVolume    = new TGeoVolume(Form(fgkCuName,id),cuShape,fCopper);

  // Chip layer
  TGeoXtru*   chipShape   = new TGeoXtru(2);
  chipShape->DefinePolygon(6, xs.fArray, ys.fArray);
  chipShape->DefineSection(0, - chipThick/2);
  chipShape->DefineSection(1, chipThick/2);
  TGeoVolume* chipVolume = new TGeoVolume(Form(fgkChipName,id),
					  chipShape,fChip);

  // Short leg shape 
  TGeoTube*   shortLegShape  = new TGeoTube(0, legr, legl / 2);
  TGeoVolume* shortLegVolume = new TGeoVolume(Form(fgkShortLegName, id), 
					      shortLegShape, fCopper);

  // Long leg shape
  TGeoTube*   longLegShape   = new TGeoTube(0, legr, (legl + modSpace) / 2);
  TGeoVolume* longLegVolume  = new TGeoVolume(Form(fgkLongLegName, id), 
					      longLegShape, fCopper);
  
  
  // Back container volume 
  TGeoVolume* backVolume     = new TGeoVolumeAssembly(Form(fgkBackVName, id));
  Double_t x = 0;
  Double_t y = 0;
  Double_t z = siThick + space + pcbThick / 2;
  backVolume->AddNode(pcbVolume, 0, new TGeoTranslation(x,y,z));
  z          += (pcbThick + cuThick) / 2;
  backVolume->AddNode(cuVolume, 0, new TGeoTranslation(0, 0, z));
  z          += (cuThick + chipThick) / 2;
  backVolume->AddNode(chipVolume, 0, new TGeoTranslation(0, 0, z));
  x          =  a->X() + legoff + legr;
  y          =  0;
  z          += pcbThick / 2 + legl / 2;
  backVolume->AddNode(shortLegVolume, 0, new TGeoTranslation(x,y,z));
  x          =  c->X();
  y          =  c->Y() - legoff - legr - off;
  backVolume->AddNode(shortLegVolume, 1, new TGeoTranslation(x,y,z));
  y          =  -y;
  backVolume->AddNode(shortLegVolume, 2, new TGeoTranslation(x,y,z));

  // Front container volume 
  TGeoVolume* frontVolume    = new TGeoVolumeAssembly(Form(fgkFrontVName, id));
  x         =  0;
  y         =  0;
  z         =  siThick + space + pcbThick / 2;
  frontVolume->AddNode(pcbVolume, 1, new TGeoTranslation(x,y,z));
  z          += (pcbThick + cuThick) / 2;
  frontVolume->AddNode(cuVolume, 0, new TGeoTranslation(0, 0, z));
  z          += (cuThick + chipThick) / 2;
  frontVolume->AddNode(chipVolume, 0, new TGeoTranslation(0, 0, z));
  x         =  a->X() + legoff + legr;
  y         =  0;
  z         += pcbThick / 2 + (legl + modSpace)/ 2;
  frontVolume->AddNode(longLegVolume, 0, new TGeoTranslation(x,y,z));
  x         =  c->X();
  y         =  c->Y() - legoff - legr - off;
  frontVolume->AddNode(longLegVolume, 1, new TGeoTranslation(x,y,z));
  y         =  -y;
  frontVolume->AddNode(longLegVolume, 2, new TGeoTranslation(x,y,z));


  // FMDD 
  Double_t ddlr = r->GetFMDDLowR();
  Double_t ddhr = r->GetFMDDHighR();
  Double_t ddpt = r->GetFMDDPrintboardThickness();
  Double_t ddct = r->GetFMDDCopperThickness();
  Double_t ddit = r->GetFMDDChipThickness();
  Double_t ddt  = ddpt + ddct + ddit;
  
  TGeoShape* fmddPcbShape  = new TGeoTubeSeg(ddlr, ddhr, ddpt/2,0,180);
  TGeoShape* fmddCuShape   = new TGeoTubeSeg(ddlr, ddhr, ddct/2,0,180);
  TGeoShape* fmddChipShape = new TGeoTubeSeg(ddlr, ddhr, ddit/2,0,180);
  fmddPcbShape->SetName(Form(fgkFMDDPCBName, id));
  fmddCuShape->SetName(Form(fgkFMDDCuName, id));
  fmddChipShape->SetName(Form(fgkFMDDChipName, id));
  if (id == 'O' || id == 'o') { 
    TString pcbName(fmddPcbShape->GetName());
    TString cuName(fmddCuShape->GetName());
    TString chipName(fmddChipShape->GetName());
    
    fmddPcbShape->SetName(Form("%s_inner",  pcbName.Data()));
    fmddCuShape->SetName(Form("%s_inner",   cuName.Data()));
    fmddChipShape->SetName(Form("%s_inner", chipName.Data()));
    new TGeoBBox(Form("%s_clip",  pcbName.Data()), ddlr+3, ddhr/2, ddpt);
    new TGeoBBox(Form("%s_clip",  cuName.Data()),  ddlr+3, ddhr/2, ddpt);
    new TGeoBBox(Form("%s_clip",  chipName.Data()),ddlr+3, ddhr/2, ddpt);
    TGeoTranslation* trans = new TGeoTranslation(Form("%s_trans",
						      pcbName.Data()), 
						 0, ddhr/2, 0);
    trans->RegisterYourself();
    fmddPcbShape = new TGeoCompositeShape(pcbName.Data(), 
					  Form("%s_inner*%s_clip:%s_trans",
					       pcbName.Data(), 
					       pcbName.Data(), 
					       pcbName.Data())); 
    fmddCuShape = new TGeoCompositeShape(cuName.Data(), 
					 Form("%s_inner*%s_clip:%s_trans",
					      cuName.Data(), 
					      cuName.Data(), 
					      pcbName.Data()));
    fmddChipShape = new TGeoCompositeShape(chipName.Data(), 
					   Form("%s_inner*%s_clip:%s_trans",
						chipName.Data(), 
						chipName.Data(), 
						pcbName.Data()));
  }

  TGeoVolume*  fmddPcbVolume = new TGeoVolume(Form(fgkFMDDPCBName, id),
					      fmddPcbShape, fPCB);
  TGeoVolume*  fmddCuVolume  = new TGeoVolume(Form(fgkFMDDCuName, id),
					      fmddCuShape, fCopper);
  TGeoVolume*  fmddChipVolume= new TGeoVolume(Form(fgkFMDDChipName, id),
					      fmddChipShape, fChip);
  // Half ring mother volumes. 
  TGeoVolume* ringTopVolume = new TGeoVolumeAssembly(Form(fgkRingTopName,id));
  TGeoVolume* ringBotVolume = new TGeoVolumeAssembly(Form(fgkRingBotName,id));
  TGeoVolume* halfRing      = ringTopVolume;

  // Adding modules to half-rings
  Int_t    nmod =  r->GetNModules();
  AliFMDDebug(10, ("making %d modules in ring %c", nmod, id));
  for (Int_t i = 0; i < nmod; i++) {
    if (i == nmod / 2) halfRing = ringBotVolume;
    Bool_t      front =  (i % 2 == 0);
    TGeoVolume* vol   =  (front ? frontVolume : backVolume);
    vol->AddNode(sensorVolume, i, new TGeoTranslation(0,0,siThick/2));
    Double_t    z1    =  (i % 2) * modSpace;
    Double_t    th    =  (2 * i + 1) * theta;
    TGeoMatrix* mat1  =  new TGeoCombiTrans(0,0,z1,0); 
    mat1->RotateZ(th);
    halfRing->AddNode(vol, i, mat1);
#if 0
    Double_t    z2    =  z1 + siThick / 2 + space;
    Double_t    th    =  (2 * i + 1) * theta;
    AliFMDDebug(20, ("Placing copy %d of %s and %s in %s at z=%f and %f, "
		      "and theta=%f", i, sensorVolume->GetName(), 
		      vol->GetName(), halfRing->GetName(), z1, z2, th));
    TGeoMatrix* mat1  =  new TGeoCombiTrans(0,0,z1,0); 
    mat1->RotateZ(th);
    halfRing->AddNode(sensorVolume, i, mat1);
    TGeoMatrix* mat2  =  new TGeoCombiTrans(0,0,z2,0); 
    mat2->RotateZ(th);
    halfRing->AddNode(vol, i, mat2);
#endif
  }

  // Add the FMDD 
  Double_t zi = r->GetFullDepth() - ddt;
  Int_t    n  = 2;
  for (Int_t i = 0; i  < n; i++) {
    TGeoVolume*   halfRing = (i == 0 ? ringTopVolume : ringBotVolume);
    Double_t      phi    = 360. / n * i;
    TGeoRotation* rot    = new TGeoRotation(Form("FMDD%c rotation %d", id, i));
    rot->RotateZ(phi);
    z         =  zi + ddpt / 2;
    halfRing->AddNode(fmddPcbVolume, i, new TGeoCombiTrans(0,0,z,rot));
    z          += (ddpt + ddct) / 2;
    halfRing->AddNode(fmddCuVolume, i, new TGeoCombiTrans(0,0,z,rot));
    z          += (ddct + ddit) / 2;
    halfRing->AddNode(fmddChipVolume, i, new TGeoCombiTrans(0,0,z,rot));
  }
  

  return 0;
}

//____________________________________________________________________
TGeoVolume*
AliFMDGeometryBuilder::DetectorGeometry(AliFMDDetector* d, 
					TGeoVolume* topMother, 
					TGeoVolume* botMother, 
					Double_t    zMother, 
					TGeoVolume* innerTop, 
					TGeoVolume* innerBot, 
					TGeoVolume* outerTop, 
					TGeoVolume* outerBot) 
{
  // Common stuff for setting up the FMD1, FMD2, and FMD3 geometries.
  // This includes putting the Honeycomb support plates and the rings
  // into the mother volumes.   
  // 
  // Parameeters:
  //	d	  The detector geometry to use 
  //	mother	  The mother volume of the detector 
  //    zmother	  The midpoint in global coordinates of detector vol.
  //	inner	  Pointer to inner ring volume 
  //    outer	  Pointer to outer ring volume
  //
  // Returns:
  //    Pointer to mother (detector volume) 
  // 
  if (!d) return 0;
  // Loop over the defined rings 
  for (int i = 0; i < 2; i++) {
    AliFMDRing* r     = 0;
    Double_t    lowr  = 0;
    Double_t    highr = 0;
    Double_t    rz    = 0;
    TGeoVolume* tvol  = 0;
    TGeoVolume* bvol  = 0;
    switch (i) {
    case 0: 
      r      = d->GetInner();
      lowr   = d->GetInnerHoneyLowR();
      highr  = d->GetInnerHoneyHighR();
      rz     = d->GetInnerZ();
      tvol   = innerTop;
      bvol   = innerBot;
      break;
    case 1: 
      r      = d->GetOuter();
      lowr   = d->GetOuterHoneyLowR();
      highr  = d->GetOuterHoneyHighR();
      rz     = d->GetOuterZ();
      tvol   = outerTop;
      bvol   = outerBot;
      break;
    }
    if (!r) continue;
    Char_t   c       = r->GetId();
    Int_t    id      = d->GetId();
    Double_t hcThick = r->GetHoneycombThickness();
    Double_t alThick = r->GetAlThickness();
    Double_t z       = TMath::Abs(rz - zMother);

    // Place ring in mother volume
    // TGeoMatrix*matrix=new TGeoTranslation(Form("FMD%d%c trans",id,c),0,0,0);
    AliFMDDebug(1, ("Placing volumes %s and %s in %s and %s at z=%f", 
		     tvol->GetName(), bvol->GetName(), 
		     topMother->GetName(), botMother->GetName(), z));
    topMother->AddNode(tvol, Int_t(c), new TGeoTranslation(0,0,z));
    botMother->AddNode(bvol, Int_t(c), new TGeoTranslation(0,0,z));

    // Top of Honeycomb
    TGeoTubeSeg* hcSha = new TGeoTubeSeg(lowr, highr, hcThick/2, 0, 180);
    TGeoVolume*  hcVol = new TGeoVolume(Form(fgkHCName,id,c),hcSha,fAl);
    // Air in top of honeycomb
    TGeoTubeSeg* ihcSha = new TGeoTubeSeg(lowr+alThick, highr - alThick, 
					     (hcThick-alThick)/2, 0, 180);
    TGeoVolume*  ihcVol = new TGeoVolume(Form(fgkIHCName,id,c),ihcSha,fAir);
    hcVol->AddNode(ihcVol, 0);
    hcVol->VisibleDaughters(kFALSE);    
    hcVol->SetVisibility(kTRUE);
    
    z += (r->GetSiThickness() + 
	  r->GetSpacing() + 
	  r->GetPrintboardThickness() + 
	  r->GetCopperThickness() + 
	  r->GetChipThickness() + 
	  r->GetModuleSpacing() +
	  r->GetLegLength() + 
	  r->GetHoneycombThickness() + 
	  r->GetFMDDPrintboardThickness() - 
	  hcThick / 2); 

    AliFMDDebug(15, ("Placing a copy of %s in %s and %s at z=%f", 
		      hcVol->GetName(), topMother->GetName(), 
		      botMother->GetName(), z));
    // Add to top 
    topMother->AddNode(hcVol, 0, new TGeoTranslation(0, 0, z));

    // Add to bottom
    TGeoMatrix*   bhcMatrix = new TGeoCombiTrans(0,0,z,0);
    bhcMatrix->RotateZ(180);
    botMother->AddNode(hcVol, 1, bhcMatrix);
  }
  return 0;
}

//____________________________________________________________________
TGeoVolume*
AliFMDGeometryBuilder::FMD1Geometry(AliFMD1* fmd1, 
				    TGeoVolume* innerTop, 
				    TGeoVolume* innerBot) 
{
  // Setup the FMD1 geometry.  The FMD1 only has one ring, and no
  // special support as it is at the momement. 
  // 
  // See also AliFMDGeometryBuilder::DetectorGeometry 
  // 
  if (!fmd1 || !innerTop || !innerBot) return 0;
  AliFMDRing* r             = fmd1->GetInner();
  Double_t    z             = fmd1->GetInnerZ();  
  Double_t    disce         = 2;
  Double_t    backlr        = fmd1->GetInnerHoneyHighR();
  Double_t    backhr        = fmd1->GetInnerHoneyHighR()+5;
  Double_t    backth        = 0.2;
  Double_t    toplr         = r->GetLowR();
  Double_t    tophr         = fmd1->GetInnerHoneyHighR()+disce;
  Double_t    wallbh        = (r->GetFullDepth() + disce);
  Double_t    wallth        = wallbh+0.1;
  
  TGeoVolume* fmd1TopVolume = new TGeoVolumeAssembly(Form(fgkFMDName, 
							  fmd1->GetId(), 'T'));
  TGeoVolume* fmd1BotVolume = new TGeoVolumeAssembly(Form(fgkFMDName, 
							  fmd1->GetId(), 'B'));
  
  // Basic detector geometry 
  DetectorGeometry(fmd1, fmd1TopVolume, fmd1BotVolume, z, 
		   innerTop, innerBot, 0, 0);


  // Back
  TGeoTubeSeg* backShape  = new TGeoTubeSeg(backlr, backhr, backth / 2, 0, 180);
  TGeoTubeSeg* wallbShape = new TGeoTubeSeg(backlr, backlr + backth, 
					    wallbh/2, 0, 180);
  TGeoTubeSeg* topShape   = new TGeoTubeSeg(toplr, tophr, backth / 2, 0, 180);
  TGeoTubeSeg* walltShape = new TGeoTubeSeg(tophr, tophr + backth, 
					    wallth/2, 0, 180);
  TGeoVolume*  backVolume = new TGeoVolume(Form(fgkBackName, fmd1->GetId()), 
					   backShape, fC);
  TGeoVolume*  wallbVolume= new TGeoVolume(Form(fgkFlangeName, fmd1->GetId()), 
					   wallbShape, fC);
  TGeoVolume*  topVolume  = new TGeoVolume(Form(fgkTopName, fmd1->GetId()), 
					   topShape, fC);
  TGeoVolume*  walltVolume= new TGeoVolume(Form(fgkBeamName, fmd1->GetId()), 
					   walltShape, fC);
  backVolume->SetFillColor(kGray);
  topVolume->SetFillColor(kGray);
  wallbVolume->SetFillColor(kGray);
  walltVolume->SetFillColor(kGray);
  
  // Place volumes
  Double_t zb = TMath::Abs(fmd1->GetInnerZ() - z);
  Double_t zi = zb;
  Int_t    n  = 2;
  
  // Place top cover
  zi -= disce / 2 + backth / 2;
  zb =  zi;
  for (Int_t i = 0; i  < 2; i++) {
    TGeoVolume*   mother = (i == 0 ? fmd1TopVolume : fmd1BotVolume);
    Double_t      phi    = 360. / n * i;
    TGeoRotation* rot    = new TGeoRotation(Form("FMD1 top rotation %d",
						 i));
    rot->RotateZ(phi);
    TGeoMatrix* matrix   = new TGeoCombiTrans(Form("FMD1 top wall trans %d", 
						   i),
					    0, 0, zi, rot);
    mother->AddNode(topVolume, i, matrix);    
  }
  // Place outer wall
  zi += wallth / 2 + backth / 2;
  for (Int_t i = 0; i  < 2; i++) {
    TGeoVolume*   mother = (i == 0 ? fmd1TopVolume : fmd1BotVolume);
    Double_t      phi    = 360. / n * i;
    TGeoRotation* rot    = new TGeoRotation(Form("FMD1 outer wall rotation %d",
						 i));
    rot->RotateZ(phi);
    TGeoMatrix* matrix   = new TGeoCombiTrans(Form("FMD1 outer wall trans %d", 
						   i),
					    0, 0, zi, rot);
    mother->AddNode(walltVolume, i, matrix);    
  }
  // Place back
  zi += wallth / 2 + backth / 2; // + disce / 2;
  for (Int_t i = 0; i  < 2; i++) {
    TGeoVolume*   mother = (i == 0 ? fmd1TopVolume : fmd1BotVolume);
    Double_t      phi    = 360. / n * i;
    TGeoRotation* rot    = new TGeoRotation(Form("FMD1 back rotation %d", i));
    rot->RotateZ(phi);
    TGeoMatrix* matrix   = new TGeoCombiTrans(Form("FMD1 back trans %d", i),
					     0, 0, zi, rot);
    mother->AddNode(backVolume, i, matrix);    
  }
  // Place inner wall
  zi -= wallbh / 2 + backth / 2; // + disce / 2;
  for (Int_t i = 0; i  < 2; i++) {
    TGeoVolume*   mother = (i == 0 ? fmd1TopVolume : fmd1BotVolume);
    Double_t      phi    = 360. / n * i;
    TGeoRotation* rot    = new TGeoRotation(Form("FMD1 inner wall rotation %d",
						 i)); 
    rot->RotateZ(phi);
    TGeoMatrix*   matrix = new TGeoCombiTrans(Form("FMD1 inner wall trans %d", 
						   i),
					      0, 0, zi, rot);
    mother->AddNode(wallbVolume, i, matrix);    
  }


  // Must add this after filling the assembly.
  TGeoVolume* top    = gGeoManager->GetVolume("ALIC");
  // TGeoMatrix* matrix = new TGeoTranslation("FMD1 trans", 0, 0, z);
  TGeoRotation* rot = new TGeoRotation("FMD1 rotatation");
  rot->RotateZ(-90);
  TGeoMatrix* matrix = new TGeoCombiTrans("FMD1 trans", 0, 0, z, rot);
  AliFMDDebug(5, ("Placing volumes %s and %s in ALIC at z=%f", 
		   fmd1TopVolume->GetName(), fmd1BotVolume->GetName(), z));
  top->AddNode(fmd1TopVolume, fmd1->GetId(), matrix);
  top->AddNode(fmd1BotVolume, fmd1->GetId(), matrix);

  return 0;
}

//____________________________________________________________________
TGeoVolume*
AliFMDGeometryBuilder::FMD2Geometry(AliFMD2* fmd2, 
				    TGeoVolume* innerTop, 
				    TGeoVolume* innerBot, 
				    TGeoVolume* outerTop,
				    TGeoVolume* outerBot) 
{
  // Setup the FMD2 geometry.  The FMD2 has no
  // special support as it is at the momement. 
  // 
  // See also AliFMDGeometryBuilder::DetectorGeometry 
  // 
  if (!fmd2 || !innerTop || !innerBot || !outerTop || !outerBot) return 0;
  AliFMDRing* r             = fmd2->GetOuter();
  Double_t    z             = fmd2->GetOuterZ();  
  Double_t    framelr       = fmd2->GetOuterHoneyHighR()+0.5;
  Double_t    framehr       = fmd2->GetOuterHoneyHighR()+1.8;
  Double_t    framelz       = -1;
  Double_t    framehz       = (fmd2->GetInnerZ()-z) + r->GetFullDepth() + 1;
  Double_t    framel        = framehz - framelz;
  Double_t    coverlr       = fmd2->GetInner()->GetLowR()+1;
  Double_t    backth        = 0.05;

  TGeoVolume* fmd2TopVolume = new TGeoVolumeAssembly(Form(fgkFMDName, 
							  fmd2->GetId(), 'T'));
  TGeoVolume* fmd2BotVolume = new TGeoVolumeAssembly(Form(fgkFMDName, 
							  fmd2->GetId(), 'B'));
  
  DetectorGeometry(fmd2, fmd2TopVolume, fmd2BotVolume, z, 
		   innerTop, innerBot, outerTop, outerBot);

  TGeoShape*  cylinderShape   = new TGeoTubeSeg(framelr,framehr,framel/2,0,180);
  TGeoVolume* cylinderVolume  = new TGeoVolume(Form(fgkBackName, fmd2->GetId()),
					       cylinderShape, fC);
  TGeoShape*  coverShape      = new TGeoTubeSeg(coverlr,framehr,backth/2,0,180);
  TGeoVolume* coverVolume     = new TGeoVolume(Form(fgkTopName, fmd2->GetId()), 
					       coverShape, fC);
  cylinderVolume->SetTransparency(63);
  coverVolume->SetTransparency(63);
  
  for (Int_t i = 0; i  < 2; i++) {
    TGeoVolume*   mother = (i == 0 ? fmd2TopVolume : fmd2BotVolume);
    
    Double_t      phi    = 360. / 2 * i;
    TGeoRotation* rot    = new TGeoRotation(Form("FMD2 support rot %d",i)); 
    rot->RotateZ(phi);
    TGeoMatrix*   matrix = new TGeoCombiTrans(Form("FMD2 cyl trans %d", i),
					      0, 0, framelz+framel/2, rot);
    mother->AddNode(cylinderVolume, i, matrix);    
    matrix               = new TGeoCombiTrans(Form("FMD2 fcov trans %d", i),
					      0, 0, framelz-backth/2, rot);
    mother->AddNode(coverVolume, 2*i+0, matrix);    
    matrix               = new TGeoCombiTrans(Form("FMD2 bcov trans %d", i),
					      0, 0, framelz+framel+backth/2, 
					      rot);
    mother->AddNode(coverVolume, 2*i+1, matrix);    
  }


  Double_t    f1l           = 10;
  Double_t    f1w           = 6;
  Double_t    f1d           = 1.2;
  
  TGeoBBox*   flange1Shape  = new TGeoBBox(f1l/2, f1w/2, f1d/2);
  TGeoVolume* flange1Volume = new TGeoVolume(Form(fgkFlangeName, fmd2->GetId()),
					     flange1Shape, fAl);
  TGeoBBox*   flange2Shape  = new TGeoBBox(f1w/2, f1d/2, (framel+backth)/2);
  TGeoVolume* flange2Volume = new TGeoVolume(Form("F%dSG", fmd2->GetId()),
					     flange2Shape, fAl);
  flange1Volume->SetTransparency(42);
  for (Int_t i = 0; i  < 4; i++) {
    TGeoVolume*   mother = (i < 2 ? fmd2TopVolume : fmd2BotVolume);
    
    Double_t      phi    = 360. / 4 * i - 45;
    Double_t      rphi   = TMath::Pi()*phi/180;
    Double_t      x      = (framelr + f1l/2) * TMath::Sin(rphi);
    Double_t      y      = (framelr + f1l/2) * TMath::Cos(rphi);
    TGeoRotation* rot    = new TGeoRotation(Form("FMD2 support rot %d",i)); 
    rot->RotateZ(phi);
    TGeoMatrix*   matrix = new TGeoCombiTrans(Form("FMD2 flange 1 trans %d", i),
					      x,y, framelz-backth-f1d/2, rot);
    mother->AddNode(flange1Volume, 2*i+0, matrix);    
    matrix               = new TGeoCombiTrans(Form("FMD2 flange 2 trans %d", i),
					      x,y,framelz+framel+backth+f1d/2, 
					      rot);
    mother->AddNode(flange1Volume, 2*i+1, matrix);    
    Double_t x1 = x - (f1w-f1d) / 2 * TMath::Cos(rphi); 
    Double_t y1 = y + (f1w-f1d) / 2 * TMath::Sin(rphi);
    matrix               = new TGeoCombiTrans(Form("FMD2 flange 3 trans %d", i),
					      x1,y1,framelz+framel/2, rot);
    mother->AddNode(flange2Volume, 2*i+0, matrix);    
    Double_t x2 = x + (f1w-f1d) / 2 * TMath::Cos(rphi); 
    Double_t y2 = y - (f1w-f1d) / 2 * TMath::Sin(rphi);
    matrix               = new TGeoCombiTrans(Form("FMD2 flange 4 trans %d", i),
					      x2,y2,framelz+framel/2, rot);
    mother->AddNode(flange2Volume, 2*i+1, matrix);    
  }
  
  

  // Must be done after filling the assemblies 
  TGeoVolume* top = gGeoManager->GetVolume("ALIC");
  TGeoMatrix* matrix = new TGeoTranslation("FMD2 trans", 0, 0, z);
  AliFMDDebug(5, ("Placing volumes %s and %s in ALIC at z=%f", 
		   fmd2TopVolume->GetName(), fmd2BotVolume->GetName(), z));
  top->AddNode(fmd2TopVolume, fmd2->GetId(), matrix);
  top->AddNode(fmd2BotVolume, fmd2->GetId(), matrix);


  return 0;
}
  
//____________________________________________________________________
TGeoVolume*
AliFMDGeometryBuilder::FMD3Geometry(AliFMD3* fmd3, 
				    TGeoVolume* innerTop, 
				    TGeoVolume* innerBot, 
				    TGeoVolume* outerTop,
				    TGeoVolume* outerBot) 
{
  // Setup the FMD3 geometry.  The FMD2 has a rather elaborate support
  // structure, as the support will also support the vacuum
  // beam-pipe. 
  // 
  // See also AliFMDGeometryBuilder::DetectorGeometry 
  // 
  if (!fmd3 || !innerTop || !innerBot || !outerTop || !outerBot) return 0;
  Double_t nlen    = fmd3->GetNoseLength();
  Double_t nz      = fmd3->GetNoseZ();
  Double_t noser1  = fmd3->GetNoseLowR();
  Double_t noser2  = fmd3->GetNoseHighR();
  Double_t conet   = fmd3->GetBeamThickness();
  Double_t conel   = fmd3->GetConeLength();
  Double_t backl   = fmd3->GetBackLength();
  Double_t backr1  = fmd3->GetBackLowR();
  Double_t backr2  = fmd3->GetBackHighR();
  Double_t zdist   = conel -  backl - nlen;
  Double_t tdist   = backr2 - noser2;
  Double_t beaml   = TMath::Sqrt(zdist * zdist + tdist * tdist);
  Double_t theta   = -180. * TMath::ATan2(tdist, zdist) / TMath::Pi();
  Double_t flanger = fmd3->GetFlangeR();
  Double_t z       = fmd3->GetInnerZ(); // fmd3->GetZ();
  Double_t zi;

  TGeoVolume* fmd3TopVolume = new TGeoVolumeAssembly(Form(fgkFMDName, 
							  fmd3->GetId(), 'T'));
  TGeoVolume* fmd3BotVolume = new TGeoVolumeAssembly(Form(fgkFMDName, 
							  fmd3->GetId(), 'B'));

  
  DetectorGeometry(fmd3, fmd3TopVolume, fmd3BotVolume, z, 
		   innerTop, innerBot, outerTop, outerBot);

  
  TGeoVolumeAssembly* support = new TGeoVolumeAssembly("F3SU");
  
  // Nose volume 
  TGeoTubeSeg* noseShape  = new TGeoTubeSeg(noser1, noser2, nlen / 2, 0, 180);
  TGeoVolume*  noseVolume = new TGeoVolume(fgkNoseName, noseShape, fC);
  support->AddNode(noseVolume, 0, new TGeoTranslation(0, 0, nlen/2));
  
  // Steel bolts 
  TGeoTube*       boltShape  = new TGeoTube("F3SB", 0, 0.3, conet / 2);
  TGeoVolume*     boltVolume = new TGeoVolume("F3SB", boltShape, fSteel);
  Double_t        z1         = -10;
  Double_t        x1         = (fmd3->ConeR(nz+z1));
  TGeoRotation*   r1         = new TGeoRotation();
  r1->RotateY(theta);
  TGeoCombiTrans* t          = new TGeoCombiTrans("F3SB1",x1,0,-z1,r1);
  support->AddNode(boltVolume, 1, t);
  z1                         = -20;
  x1                         = (fmd3->ConeR(nz+z1));
  t                          = new TGeoCombiTrans("F3SB2",x1,0,-z1,r1);
  support->AddNode(boltVolume, 2, t);

  // Cooling plates
  TGeoTrd1*   plateShape  = new TGeoTrd1(2, 8, 0.1, (conel-2-2)/2-.1);
  TGeoVolume* plateVolume = new TGeoVolume("F3CO", plateShape, fAl);

  // Shape for carbon half-cone
  new TGeoConeSeg("F3SC_inner", conel/2,noser2-conet, noser2, 
		  backr2-conet, backr2, 0., 180.);
  new TGeoTrd1("F3SC_hole",2,8,conet*3,(conel-2-2)/2);
  Double_t        holeAng   = TMath::ATan2(backr2 - noser2, conel);
  Double_t        holeX     = ((conel-2) / 2 * TMath::Sin(holeAng) +
			       conet     * TMath::Cos(holeAng) +
			       noser2);
  TGeoRotation*   holeRot   = new TGeoRotation();
  holeRot->RotateZ(90);
  holeRot->RotateY(holeAng*180./TMath::Pi());
  TGeoCombiTrans* holeTrans = new TGeoCombiTrans(holeX, 0, -2, holeRot);

  // Build-up the composite shape for the cone, and add cooling plates
  // at the same time. 
  TString coneExp("F3SC_inner-(");
  for (int i = 0; i < 4; i++) { 
    Double_t        thisAng   = 360. / 8 * (i + .5);
    TGeoCombiTrans* thisTrans = new TGeoCombiTrans(*holeTrans);
    thisTrans->RotateZ(thisAng);
    thisTrans->SetName(Form("F3SC_rot%d", i));
    thisTrans->RegisterYourself();
    coneExp.Append(Form("F3SC_hole:F3SC_rot%d+", i));

    const Double_t* tt         = thisTrans->GetTranslation();
    Double_t        x          = tt[0]+1*TMath::Cos(thisAng*TMath::Pi()/180);
    Double_t        y          = tt[1]+1*TMath::Sin(thisAng*TMath::Pi()/180);
    TGeoCombiTrans* plateTrans = new TGeoCombiTrans(x,y,tt[2]-1+nlen+conel/2,
						    thisTrans->GetRotation());
    support->AddNode(plateVolume, i, plateTrans);
  }
  // Remove bolt holes 
  coneExp.Append("F3SB:F3SB1+F3SB:F3SB2)");

  // Finalize the half-cone shape and add volume
  TGeoCompositeShape* coneShape  = new TGeoCompositeShape(coneExp.Data());
  TGeoVolume*         coneVolume = new TGeoVolume("F3SC", coneShape, fC);
  support->AddNode(coneVolume,1,new TGeoTranslation(0,0,nlen+conel/2));
  
  // The flanges 
  TGeoBBox* flangeShape    = new TGeoBBox((flanger - backr2) / 2, 
					  fmd3->GetBeamWidth() / 2,
					  backl / 2);
  TGeoVolume* flangeVolume = new TGeoVolume(Form(fgkFlangeName, fmd3->GetId()),
					    flangeShape, fC);
  Int_t    n               = fmd3->GetNFlange();
  Double_t r               = backr2 + (flanger - backr2) / 2;
  for (Int_t i = 0; i  < n/2; i++) {
    Double_t phi       = 360. / n * i + 180. / n;
    Double_t x         = r * TMath::Cos(TMath::Pi() / 180 * phi);
    Double_t y         = r * TMath::Sin(TMath::Pi() / 180 * phi);
    TGeoRotation* rot  = new TGeoRotation;
    rot->RotateZ(phi);
    TGeoMatrix* matrix = new TGeoCombiTrans(x, y, nlen+conel-backl/2, rot);
    support->AddNode(flangeVolume, i, matrix);
  }

  // Place support volumes in half-detector volumes 
  z                          = fmd3->GetInnerZ();
  z1                         = z-nz;
  fmd3TopVolume->AddNode(support, 1, new TGeoTranslation(0,0,z1));
  r1                         = new TGeoRotation();
  r1->RotateZ(180);
  t                          = new TGeoCombiTrans(0,0,z1,r1);
  fmd3BotVolume->AddNode(support, 2, t);

  TGeoRotation*   rot        = new TGeoRotation("FMD3 rotatation");
  rot->RotateY(180);
  TGeoVolume*     top        = gGeoManager->GetVolume("ALIC");
  TGeoMatrix* mmatrix        = new TGeoCombiTrans("FMD3 trans", 0, 0, z, rot);
  AliFMDDebug(5, ("Placing volumes %s and %s in ALIC at z=%f", 
		   fmd3TopVolume->GetName(), fmd3BotVolume->GetName(), z));
  top->AddNode(fmd3TopVolume, fmd3->GetId(), mmatrix);
  top->AddNode(fmd3BotVolume, fmd3->GetId(), mmatrix);

  return 0;
}

//____________________________________________________________________
void
AliFMDGeometryBuilder::Exec(Option_t*) 
{
  // Setup up the FMD geometry. 
  AliFMDDebug(1, ("\tGeometry options: %s",
		    (fDetailed  ? "divided into strips" : "one volume")));
  if (!gGeoManager) {
    AliFatal("No TGeoManager defined");
    return;
  }

  fSi      = gGeoManager->GetMedium("FMD_Si$");
  fC       = gGeoManager->GetMedium("FMD_Carbon$");
  fAl      = gGeoManager->GetMedium("FMD_Aluminum$");
  fChip    = gGeoManager->GetMedium("FMD_Si Chip$");
  fAir     = gGeoManager->GetMedium("FMD_Air$");
  fPCB     = gGeoManager->GetMedium("FMD_PCB$");
  fPlastic = gGeoManager->GetMedium("FMD_Plastic$");
  fCopper  = gGeoManager->GetMedium("FMD_Copper$");
  fSteel   = gGeoManager->GetMedium("FMD_Steel$");

  if (!fSi||!fC||!fAl||!fChip||!fAir||!fPCB||!fPlastic||!fCopper||!fSteel) {
    AliError("Failed to get some or all tracking mediums");
    return;
  }    
  AliFMDGeometry* fmd = AliFMDGeometry::Instance();
  AliFMDRing* inner = fmd->GetInner();
  AliFMDRing* outer = fmd->GetOuter();
  RingGeometry(inner);
  RingGeometry(outer);
  TGeoVolume* innerTop = gGeoManager->GetVolume(Form(fgkRingTopName, 
						     inner->GetId()));
  TGeoVolume* innerBot = gGeoManager->GetVolume(Form(fgkRingBotName, 
						     inner->GetId()));
  TGeoVolume* outerTop = gGeoManager->GetVolume(Form(fgkRingTopName, 
						     outer->GetId()));
  TGeoVolume* outerBot = gGeoManager->GetVolume(Form(fgkRingBotName, 
						     outer->GetId()));
  
  FMD1Geometry(fmd->GetFMD1(), innerTop, innerBot);
  FMD2Geometry(fmd->GetFMD2(), innerTop, innerBot, outerTop, outerBot);
  FMD3Geometry(fmd->GetFMD3(), innerTop, innerBot, outerTop, outerBot);
#ifndef USE_PRE_MOVE
  fmd->SetSectorOff(fSectorOff);
  fmd->SetModuleOff(fModuleOff);
  fmd->SetRingOff(fRingOff);
  fmd->SetDetectorOff(fDetectorOff);
  fmd->SetActive(fActiveId.fArray, fActiveId.fN);
#endif
  // fmd->ExtractGeomInfo();
  
}


//____________________________________________________________________
//
// EOF
//
