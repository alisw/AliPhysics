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
#include <TGeoPcon.h>		// ROOT_TGeoPcon
#include <TGeoCompositeShape.h>
#include <TMath.h>
#include <TVector2.h>		// ROOT_TVector2
#include <TVector3.h>		// ROOT_TVector3
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
#include <iostream>

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
  Char_t        id       = r->GetId();
  const Char_t* lName    = (id == 'i' || id == 'I' ? "inner" : "outer");
  Double_t      siThick  = r->GetSiThickness();
  const Int_t   knv      = r->GetNVerticies();
  TVector2*     a        = r->GetVertex(5);
  TVector2*     b        = r->GetVertex(3);
  TVector2*     c        = r->GetVertex(4);
  Double_t      theta    = r->GetTheta();
  Double_t      off      = (TMath::Tan(TMath::Pi() * theta / 180) 
			    * r->GetBondingWidth());
  Double_t      rmax     = b->Mod();
  Double_t      rmin     = r->GetLowR();
  Double_t      pcbThick = r->GetPrintboardThickness();
  Double_t      cuThick  = r->GetCopperThickness();
  Double_t      chipThick= r->GetChipThickness();
  Double_t      modSpace = r->GetModuleSpacing();
  Double_t      legr     = r->GetLegRadius();
  Double_t      legl     = r->GetLegLength();
  Double_t      legoff   = r->GetLegOffset();
  Int_t         ns       = r->GetNStrips();
  Double_t      stripoff = a->Mod();
  Double_t      dstrip   = (rmax - stripoff) / ns;
  Double_t      space    = r->GetSpacing();
  TArrayD       xs(knv);
  TArrayD       ys(knv);
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
  sensorShape->SetName(Form(fgkSensorName, id));
  sensorShape->SetTitle(Form("FMD %s Sensor", lName));
  TGeoVolume* sensorVolume = new TGeoVolume(Form(fgkSensorName, id), 
					    sensorShape, fSi);
  sensorVolume->SetTitle(Form("FMD %s Sensor", lName));
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
    activeShape->SetName(Form(fgkActiveName, id));
    activeShape->SetTitle(Form("FMD %s active area", lName));
    TGeoVolume* activeVolume = new TGeoVolume(Form(fgkActiveName, id),
					      activeShape,fSi);
    activeVolume->SetTitle(Form("FMD %s active area", lName));
    TGeoVolume* sectorVolume = activeVolume->Divide(Form(fgkSectorName,id), 
						      2, 2, -theta,0,0,"N");
    sectorVolume->SetTitle(Form("FMD %s sector", lName));
    TGeoVolume* stripVolume  = sectorVolume->Divide(Form(fgkStripName, id), 
						    1, ns, stripoff, dstrip, 
						    0, "SX");
    stripVolume->SetTitle(Form("FMD %s strip", lName));
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
  pcbShape->SetName(Form(fgkPCBName, id));
  pcbShape->SetTitle(Form("FMD %s hybrid PCB", lName));
  TGeoVolume* pcbVolume      = new TGeoVolume(Form(fgkPCBName, id), 
					      pcbShape, fPCB);
  pcbVolume->SetTitle(Form("FMD %s hybrid PCB", lName));

  // Copper layer
  TGeoXtru* cuShape       = new TGeoXtru(2);
  cuShape->DefinePolygon(6, xs.fArray, ys.fArray);
  cuShape->DefineSection(0, - cuThick/2);
  cuShape->DefineSection(1, cuThick/2);
  cuShape->SetTitle(Form("FMD %s hybrid copper", lName));
  TGeoVolume* cuVolume    = new TGeoVolume(Form(fgkCuName,id),cuShape,fCopper);
  cuVolume->SetTitle(Form("FMD %s hybrid copper", lName));

  // Chip layer
  TGeoXtru*   chipShape   = new TGeoXtru(2);
  chipShape->DefinePolygon(6, xs.fArray, ys.fArray);
  chipShape->DefineSection(0, - chipThick/2);
  chipShape->DefineSection(1, chipThick/2);
  chipShape->SetTitle(Form("FMD %s hybrid chip", lName));
  TGeoVolume* chipVolume = new TGeoVolume(Form(fgkChipName,id),
					  chipShape,fChip);
  chipVolume->SetTitle(Form("FMD %s hybrid chip", lName));

  // Short leg shape 
  TGeoTube*   shortLegShape  = new TGeoTube(0, legr, legl / 2);
  shortLegShape->SetName(Form(fgkShortLegName, id));
  shortLegShape->SetTitle(Form("FMD %s short support foot", lName));
  TGeoVolume* shortLegVolume = new TGeoVolume(Form(fgkShortLegName, id), 
					      shortLegShape, fCopper);
  shortLegVolume->SetTitle(Form("FMD %s short support foot", lName));
  // Long leg shape
  TGeoTube*   longLegShape   = new TGeoTube(0, legr, (legl + modSpace) / 2);
  longLegShape->SetName(Form(fgkLongLegName, id));
  longLegShape->SetTitle(Form("FMD %s long support foot", lName));
  TGeoVolume* longLegVolume  = new TGeoVolume(Form(fgkLongLegName, id), 
					      longLegShape, fCopper);
  longLegVolume->SetTitle(Form("FMD %s long support foot", lName));
  
  
  // Back container volume 
  TGeoVolume* backVolume     = new TGeoVolumeAssembly(Form(fgkBackVName, id));
  backVolume->SetTitle(Form("FMD %s back module", lName));
  Double_t x = 0;
  Double_t y = 0;
  Double_t z = siThick / 2;
  backVolume->AddNode(sensorVolume, 0, new TGeoTranslation(x, y, z));
  z          += siThick / 2 + space + pcbThick / 2;
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
  frontVolume->SetTitle(Form("FMD %s front module", lName));
  x         =  0;
  y         =  0;
  z         = siThick / 2;
  frontVolume->AddNode(sensorVolume, 0, new TGeoTranslation(x, y, z));
  z          += siThick / 2 + space + pcbThick / 2;
  frontVolume->AddNode(pcbVolume, 0, new TGeoTranslation(x,y,z));
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
  
  TString    pcbName(Form(fgkFMDDPCBName, id));
  TString    cuName(Form(fgkFMDDCuName, id));
  TString    chipName(Form(fgkFMDDChipName, id));
  new TGeoTubeSeg(Form("%s_inner", pcbName.Data()),  ddlr, ddhr, ddpt/2,0,180);
  new TGeoTubeSeg(Form("%s_inner", cuName.Data()),   ddlr, ddhr, ddct/2,0,180);
  new TGeoTubeSeg(Form("%s_inner", chipName.Data()), ddlr, ddhr, ddit/2,0,180);
  
  Double_t clipWX = 0;
  Double_t clipWY = 0;
  Double_t clipY  = 1;
  
  if (id == 'I' || id == 'i') { 
    clipWX = ddhr;
    clipWY = ddhr/2;
  }
  else { 
    clipWX = ddlr+3;
    clipWY = ddhr/2;
  }
  
  new TGeoBBox(Form("%s_clip",  pcbName.Data()), clipWX, clipWY, ddpt);
  new TGeoBBox(Form("%s_clip",  cuName.Data()),  clipWX, clipWY, ddct);
  new TGeoBBox(Form("%s_clip",  chipName.Data()),clipWX, clipWY, ddit);
  TGeoTranslation* trans = new TGeoTranslation(Form("%s_trans",
						    pcbName.Data()), 
					       0, clipWY+clipY, 0);
  trans->RegisterYourself();
  TGeoShape* fmddPcbShape = 
    new TGeoCompositeShape(pcbName.Data(), 
			   Form("%s_inner*%s_clip:%s_trans",
				pcbName.Data(), 
				pcbName.Data(), 
				pcbName.Data())); 
  TGeoShape* fmddCuShape = 
    new TGeoCompositeShape(cuName.Data(), 
			   Form("%s_inner*%s_clip:%s_trans",
				cuName.Data(), 
				cuName.Data(), 
				pcbName.Data()));
  TGeoShape* fmddChipShape = 
    new TGeoCompositeShape(chipName.Data(), 
			   Form("%s_inner*%s_clip:%s_trans",
				chipName.Data(), 
				chipName.Data(), 
				pcbName.Data()));
  fmddPcbShape->SetTitle(Form("FMD %s digitiser PCB", lName));
  fmddCuShape->SetTitle(Form("FMD %s digitiser copper", lName));
  fmddChipShape->SetTitle(Form("FMD %s digitiser chip", lName));

  TGeoVolume*  fmddPcbVolume = new TGeoVolume(Form(fgkFMDDPCBName, id),
					      fmddPcbShape, fPCB);
  TGeoVolume*  fmddCuVolume  = new TGeoVolume(Form(fgkFMDDCuName, id),
					      fmddCuShape, fCopper);
  TGeoVolume*  fmddChipVolume= new TGeoVolume(Form(fgkFMDDChipName, id),
					      fmddChipShape, fChip);
  fmddPcbVolume->SetTitle(Form("FMD %s digitiser PCB", lName));
  fmddCuVolume->SetTitle(Form("FMD %s digitiser copper", lName));
  fmddChipVolume->SetTitle(Form("FMD %s digitiser chip", lName));

  // Half ring mother volumes. 
  TGeoVolume* ringTopVolume = new TGeoVolumeAssembly(Form(fgkRingTopName,id));
  TGeoVolume* ringBotVolume = new TGeoVolumeAssembly(Form(fgkRingBotName,id));
  TGeoVolume* halfRing      = ringTopVolume;
  ringTopVolume->SetTitle(Form("FMD %s top half-ring", lName));
  ringBotVolume->SetTitle(Form("FMD %s bottom half-ring", lName));
  
  // Adding modules to half-rings
  Int_t    nmod =  r->GetNModules();
  AliFMDDebug(10, ("making %d modules in ring %c", nmod, id));
  for (Int_t i = 0; i < nmod; i++) {
    if (i == nmod / 2) halfRing = ringBotVolume;
    Bool_t      front =  (i % 2 == 0);
    TGeoVolume* vol   =  (front ? frontVolume : backVolume);
    // vol->AddNode(sensorVolume, i, new TGeoTranslation(0,0,siThick/2));
    Double_t    z1    =  (i % 2) * modSpace;
    Double_t    th    =  (2 * i + 1) * theta;
    TGeoMatrix* mat1  =  new TGeoCombiTrans(0,0,z1,0); 
    mat1->RotateZ(th);
    mat1->SetName(Form("FMD%c_module_%02d", id, i));
    mat1->SetTitle(Form("FMD %s module %2d matrix", lName, i));
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
    halfRing             = (i == 0 ? ringTopVolume : ringBotVolume);
    Double_t      phi    = 360. / n * i;
    TGeoRotation* rot    = new TGeoRotation(Form("FMDD%c rotation %d", id, i));
    rot->RotateZ(phi);
    rot->SetTitle(Form("FMD %s digitiser rotation %2d", lName, i));
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
TGeoShape*
AliFMDGeometryBuilder::HoneycombShape(Int_t id, Char_t ring,
				      double r1, double r2, 
				      double w, double t, double c)
{
  // Make a honey comb shape from passed parameters.
  // Parameters: 
  //   id	Detector identifier (1,2, or 3)
  //   ring	Ring identifier ('I' or 'O')
  //   r1       Inner radius
  //   r2       Outer radius
  //   w        width 
  //   t        Thickness of material 
  //   c        Clearing from horizontal. 
  // Return 
  //   Pointer to newly allocated composite shape. 
  TString      form  = Form("FMD%d%c_%%c_%%c", id, ring);
  double       a1    = TMath::ATan2(c, r1) * 180  / TMath::Pi();

  TString      fn    = Form(form.Data(),'F','1');
  TString      bn    = Form(form.Data(),'B','1');
  TString      cn    = Form(form.Data(),'C','O');
  TString      in    = Form(form.Data(),'R','I');
  TString      on    = Form(form.Data(),'R','O');
  TString      en    = Form(form.Data(),'E','X');
  double       y     = c;
  double       x     = r1 * TMath::Cos(TMath::Pi()*a1/180);
  new TGeoTubeSeg(fn.Data(),r1,r2,t/2,0,180);
  new TGeoTubeSeg(bn.Data(),r1,r2,t/2,0,180);
  new TGeoBBox(cn.Data(),(r2-r1)/2,t/2,w/2);
  new TGeoTubeSeg(in.Data(),r1,r1+t,w/2,0,180);
  new TGeoTubeSeg(on.Data(),r2-t,r2,w/2,0,180);
  new TGeoBBox(en.Data(),r2+.005,c/2+.005,w/2+.005);
    
  TString          ftn = Form(form.Data(),'F','T');
  TString          btn = Form(form.Data(),'F','B');
  TString          ltn = Form(form.Data(),'C','L');
  TString          rtn = Form(form.Data(),'C','R');
  TString          etn = Form(form.Data(),'E','X');
  (new TGeoTranslation(ftn.Data(),0,0,+w/2-t/2))->RegisterYourself();
  (new TGeoTranslation(btn.Data(),0,0,-w/2+t/2))->RegisterYourself();
  (new TGeoTranslation(ltn.Data(),-(x+(r2-r1)/2), y+t/2,0))->RegisterYourself();
  (new TGeoTranslation(rtn.Data(),(x+(r2-r1)/2), y+t/2,0))->RegisterYourself();
  (new TGeoTranslation(etn.Data(),0, c/2,0))->RegisterYourself();
  
  TString comp(Form("(%s:%s+%s:%s+%s+%s+%s:%s+%s:%s)-%s:%s", 
		    fn.Data(),ftn.Data(),
		    bn.Data(),btn.Data(),
		    in.Data(),on.Data(),
		    cn.Data(),ltn.Data(),
		    cn.Data(),rtn.Data(),
		    en.Data(),etn.Data()));
  TGeoCompositeShape* shape = new TGeoCompositeShape(comp.Data());
  shape->SetName(Form(fgkHCName,id,ring));
  shape->SetTitle(Form("FMD%d%c Honeycomb shape", id, ring));
  return shape;
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

    // Honeycomp 
    TGeoShape*   hcSha = HoneycombShape(id, c, lowr, highr, hcThick, alThick);
    TGeoVolume*  hcVol = new TGeoVolume(Form(fgkHCName,id,c),hcSha,fAl);
    hcVol->SetTitle(Form("FMD%d%c honeycomb shell", id, c));
    
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
    bhcMatrix->SetName(Form("FMD%d%c_honeycomp", id, c));
    bhcMatrix->SetTitle(Form("FMD%d%c honeycomp", id, c));
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
  fmd1TopVolume->SetTitle("FMD1 top half");
  TGeoVolume* fmd1BotVolume = new TGeoVolumeAssembly(Form(fgkFMDName, 
							  fmd1->GetId(), 'B'));
  fmd1BotVolume->SetTitle("FMD1 bottom half");
  
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
  backShape->SetName(Form(fgkBackName, fmd1->GetId()));
  wallbShape->SetName(Form(fgkFlangeName, fmd1->GetId()));
  topShape->SetName(Form(fgkTopName, fmd1->GetId()));
  walltShape->SetName(Form(fgkBeamName, fmd1->GetId()));
  backShape->SetTitle("FMD1 back saucer rim");
  wallbShape->SetTitle("FMD1 back saucer wall");
  topShape->SetTitle("FMD1 top lid");
  walltShape->SetTitle("FMD1 top lid wall");
  backVolume->SetFillColor(kGray);
  topVolume->SetFillColor(kGray);
  wallbVolume->SetFillColor(kGray);
  walltVolume->SetFillColor(kGray);
  backVolume->SetTitle("FMD1 back saucer rim");
  wallbVolume->SetTitle("FMD1 back saucer wall");
  topVolume->SetTitle("FMD1 top lid");
  walltVolume->SetTitle("FMD1 top lid wall");
  
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
    TGeoRotation* rot    = new TGeoRotation(Form("FMD1 top rotation %d",i));
    rot->RotateZ(phi);
    TGeoMatrix* matrix   = new TGeoCombiTrans(Form("FMD1 top wall trans %d", i),
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
						   i), 0, 0, zi, rot);
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
						   i), 0, 0, zi, rot);
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
  Double_t    framelz       = -.5;
  Double_t    framehz       = (fmd2->GetInnerZ()-z) + r->GetFullDepth() + .5;
  Double_t    framel        = framehz - framelz;
  Double_t    coverlr       = fmd2->GetInner()->GetLowR()+1;
  Double_t    backth        = 0.05;

  TGeoVolume* fmd2TopVolume = new TGeoVolumeAssembly(Form(fgkFMDName, 
							  fmd2->GetId(), 'T'));
  TGeoVolume* fmd2BotVolume = new TGeoVolumeAssembly(Form(fgkFMDName, 
							  fmd2->GetId(), 'B'));
  fmd2TopVolume->SetTitle("FMD2 top half");
  fmd2BotVolume->SetTitle("FMD2 bottom half");
  
  DetectorGeometry(fmd2, fmd2TopVolume, fmd2BotVolume, z, 
		   innerTop, innerBot, outerTop, outerBot);

  TGeoShape*  cylinderShape   = new TGeoTubeSeg(framelr,framehr,framel/2,0,180);
  TGeoVolume* cylinderVolume  = new TGeoVolume(Form(fgkBackName, fmd2->GetId()),
					       cylinderShape, fC);
  TGeoShape*  coverShape      = new TGeoTubeSeg(coverlr,framehr,backth/2,0,180);
  TGeoVolume* coverVolume     = new TGeoVolume(Form(fgkTopName, fmd2->GetId()), 
					       coverShape, fC);
  cylinderShape->SetName(Form(fgkBackName, fmd2->GetId()));
  cylinderShape->SetTitle("FMD2 cylinder");
  cylinderVolume->SetTitle("FMD2 cylinder");
  cylinderVolume->SetTransparency(63);
  coverShape->SetName(Form(fgkTopName, fmd2->GetId()));
  coverShape->SetTitle("FMD2 cover");
  coverVolume->SetTitle("FMD2 cover");
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
  Double_t    f1d           = 1;
  
  TGeoBBox*   flange1Shape  = new TGeoBBox(f1l/2, f1w/2, f1d/2);
  TGeoVolume* flange1Volume = new TGeoVolume(Form(fgkFlangeName, fmd2->GetId()),
					     flange1Shape, fAl);
  TGeoBBox*   flange2Shape  = new TGeoBBox(f1w/2, f1d/2, (framel+backth)/2);
  TGeoVolume* flange2Volume = new TGeoVolume(Form("F%dSG", fmd2->GetId()),
					     flange2Shape, fAl);
  flange1Shape->SetName(Form(fgkFlangeName, fmd2->GetId()));
  flange1Shape->SetTitle("FMD2 vertical flange");
  flange1Volume->SetTitle("FMD2 vertical flange");
  flange2Shape->SetName(Form("F%dSG", fmd2->GetId()));
  flange2Shape->SetTitle("FMD2 horizontal flange");
  flange2Volume->SetTitle("FMD2 horizontal flange ");
  
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
  
#if 1
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

  //__________________________________________________________________
  // Basic detector set-up.
  TGeoVolume* fmd3TopVolume = new TGeoVolumeAssembly(Form(fgkFMDName, 
							  fmd3->GetId(), 'T'));
  TGeoVolume* fmd3BotVolume = new TGeoVolumeAssembly(Form(fgkFMDName, 
							  fmd3->GetId(), 'B'));
  fmd3TopVolume->SetTitle("FMD3 top half");
  fmd3BotVolume->SetTitle("FMD3 bottom half");
  DetectorGeometry(fmd3, fmd3TopVolume, fmd3BotVolume, fmd3->GetInnerZ(), 
		   innerTop, innerBot, outerTop, outerBot);

  //__________________________________________________________________
  // Mother for all support material
  TGeoVolumeAssembly* support = new TGeoVolumeAssembly("F3SU");
  support->SetTitle("FMD3 support");

  //__________________________________________________________________
  // Base of cone
  const TObjArray& radii    = fmd3->ConeRadii();
  Int_t            nRadii   = radii.GetEntriesFast();
  TGeoPcon*        coneBase = new TGeoPcon("FMD3_cone_base", 0., 180., nRadii);
  TVector3*        r5       = 0;
  TVector3*        r4       = 0;
  for (Int_t i = 0; i < nRadii; i++) { 
    TVector3* v = static_cast<TVector3*>(radii.At(i));
    coneBase->DefineSection(i, v->X(), v->Y(), v->Z());
    if      (i == 5) r5 = v;
    else if (i == 4) r4 = v;
  }
  TString          coneComb("(FMD3_cone_base");

  //__________________________________________________________________
  // Flanges 
  double    flangeDepth    = fmd3->GetFlangeDepth() / 2;
  double    flangeLength   = fmd3->GetFlangeLength() / 2;
  double    flangeWidth    = fmd3->GetFlangeWidth() / 2;
  new TGeoBBox("FMD3_flange_base", flangeLength, flangeWidth, flangeDepth);

  // Fiducial holes 
  const TObjArray& fiducialHoles  = fmd3->FiducialHoles();
  double           fiducialRadius = fmd3->GetFiducialRadius();
  TGeoTube*        fiducialShape  = new TGeoTube("FMD3_fiducial_hole", 0, 
						 fiducialRadius, 
						 flangeDepth+.1);
  Int_t            nFiducialHoles = fiducialHoles.GetEntriesFast();
  double           flangeAngle    = TMath::Pi() / 4;
  double           flangeX        = r5->Y()+flangeLength;
  TVector2         flangeC(flangeX * TMath::Cos(flangeAngle), 
			   flangeX * TMath::Sin(flangeAngle));
  TString          flangeComb("FMD3_flange_base-(");
#if 0// For debugging geometry 
  TGeoVolume* fiducialVolume = new TGeoVolume("FMD3_fiducial", fiducialShape);
  fiducialVolume->SetLineColor(kGreen);
#else
  (void*)fiducialShape;
#endif
  for (Int_t i = 0; i < nFiducialHoles; i++) { 
    TVector2&        v  =  *(static_cast<TVector2*>(fiducialHoles.At(i)));
    v                   -= flangeC;
    TVector2         r  =  v.Rotate(-flangeAngle);
    TGeoTranslation* t1 =  new TGeoTranslation(r.X(),  r.Y(), 0);
    TGeoTranslation* t2 =  new TGeoTranslation(r.X(), -r.Y(), 0);
    t1->SetName(Form("FMD3_fiducial_hole_rot%d", 2*i+0));
    t2->SetName(Form("FMD3_fiducial_hole_rot%d", 2*i+1));
    t1->RegisterYourself();
    t2->RegisterYourself();
    flangeComb.Append(Form("FMD3_fiducial_hole:FMD3_fiducial_hole_rot%d+"
			   "FMD3_fiducial_hole:FMD3_fiducial_hole_rot%d%c",
			   2*i+0, 2*i+1, (i == nFiducialHoles-1 ? ')' : '+')));
#if 1 // For debugging geometry 
    // support->AddNode(fiducialVolume, 2*i+0, t1);
    // support->AddNode(fiducialVolume, 2*i+1, t2);
#endif
  }
  
  // Final flange shape, and at to full shape 
  TGeoCompositeShape* flangeShape = new TGeoCompositeShape(flangeComb.Data());
  flangeShape->SetName("FMD3_flange");
  for (Int_t i = 0; i < 2; i++) { 
    TGeoRotation* rot = new TGeoRotation();
    rot->RotateZ((i+.5)*90);
    TVector2 v(flangeX, 0);
    TVector2 w = v.Rotate((i+.5) * 2 * flangeAngle);
    TGeoCombiTrans* trans = new TGeoCombiTrans(w.X(),w.Y(),
					       r4->X()+flangeDepth, rot);
    trans->SetName(Form("FMD3_flange_matrix%d", i));
    trans->RegisterYourself();
    coneComb.Append(Form("+FMD3_flange:FMD3_flange_matrix%d", i));
  }
  coneComb.Append(")-(");
  
  //__________________________________________________________________
  // Holes 
  Double_t holeL  = (fmd3->GetHoleLength()-1)/2;
  Double_t holeD  = fmd3->GetHoleDepth()/2;
  Double_t holeLW = fmd3->GetHoleLowWidth()/2;
  Double_t holeHW = fmd3->GetHoleHighWidth()/2;
  Double_t holeA  = fmd3->GetConeOuterAngle() - 1 * TMath::Pi() / 180;
  Double_t holeZ  = (fmd3->GetHoleOffset() 
		     + holeL * TMath::Cos(holeA)
		     - holeD * TMath::Sin(holeA));
  Double_t holeX  = (fmd3->ConeR(-holeZ + fmd3->GetInnerZ() + fmd3->GetNoseZ()) 
		     - holeD * TMath::Cos(holeA));
  Double_t plateZ  = (fmd3->GetHoleOffset() 
		     + holeL * TMath::Cos(holeA)
		     - 0.033 * TMath::Sin(holeA));
  Double_t plateX = (fmd3->ConeR(-plateZ + fmd3->GetInnerZ()+fmd3->GetNoseZ()) 
		     - 0.033 * TMath::Cos(holeA));
  TGeoTrd1* holeShape = new TGeoTrd1("FMD3_cone_hole", 
				     holeLW, holeHW, holeD, holeL);
  TGeoTrd1* plateShape = new TGeoTrd1("FMD3_cooling_plate", 
				      holeLW, holeHW, .033, holeL);
  TGeoRotation* holeRot = new TGeoRotation();
  holeRot->SetName("FMD3_cone_hole_rotation");
  holeRot->RotateZ(90);
  holeRot->RotateY(holeA*180/TMath::Pi());
  TGeoCombiTrans* holeBaseTrans = new TGeoCombiTrans(holeX, 0, holeZ, holeRot);
  holeBaseTrans->SetName("FMD3_cone_hole_base_matrix");
  TGeoCombiTrans* plateBaseTrans = new TGeoCombiTrans(plateX, 0,plateZ,holeRot);
  (void*)holeShape;
  TGeoVolume* plateVolume = new TGeoVolume("FMD3_cooling_plate", 
					   plateShape, fAl);
  plateShape->SetTitle("FMD3 cooling plate");
  plateVolume->SetTitle("FMD3 cooling plate");
  for (Int_t i = 0; i < 4; i++) { 
    Double_t        ang   = 360. / 8 * (i + .5);
    TGeoCombiTrans* trans = new TGeoCombiTrans(*holeBaseTrans);
    trans->RotateZ(ang);
    trans->SetName(Form("FMD3_cone_hole_matrix%d", i));
    trans->RegisterYourself();
    trans = new TGeoCombiTrans(*plateBaseTrans);
    trans->RotateZ(ang);
    trans->SetName(Form("FMD3_cooling_plate_matrix%d", i));
    coneComb.Append(Form("FMD3_cone_hole:FMD3_cone_hole_matrix%d+", i));
    support->AddNode(plateVolume, i, trans);
  }
  
  //__________________________________________________________________
  // Bolts
  Double_t boltRadius = fmd3->GetBoltRadius();
  Double_t boltLength = fmd3->GetBoltLength() / 2;
  Double_t boltZ1     = fmd3->GetInnerZ()+fmd3->GetNoseZ()-10;
  Double_t boltZ2     = fmd3->GetInnerZ()+fmd3->GetNoseZ()-20;
  Double_t boltXE     = 2*boltLength*TMath::Cos(fmd3->GetConeOuterAngle());
  Double_t boltX1     = (fmd3->ConeR(boltZ1) - boltXE);
  Double_t boltX2     = (fmd3->ConeR(boltZ2) - boltXE);
  
  new TGeoTube("FMD3_bolt_hole", 0, boltRadius, boltLength+.2);
  TGeoTube* boltShape = new TGeoTube("FMD3_bolt", 0, boltRadius, boltLength);
  TGeoRotation* boltRot = new TGeoRotation();
  boltRot->RotateY(-fmd3->GetConeOuterAngle()*180/TMath::Pi());
  TGeoCombiTrans* boltTrans1 = new TGeoCombiTrans(boltX1, 0, 10, boltRot);
  TGeoCombiTrans* boltTrans2 = new TGeoCombiTrans(boltX2, 0, 20, boltRot);
  TGeoCombiTrans* boltTrans3 = new TGeoCombiTrans(*boltTrans1);
  TGeoCombiTrans* boltTrans4 = new TGeoCombiTrans(*boltTrans2);
  boltTrans3->RotateZ(180);
  boltTrans4->RotateZ(180);
  boltTrans1->SetName("FMD3_bolt_matrix1");
  boltTrans2->SetName("FMD3_bolt_matrix2");
  boltTrans3->SetName("FMD3_bolt_matrix3");
  boltTrans4->SetName("FMD3_bolt_matrix4");
  boltTrans1->RegisterYourself();
  boltTrans2->RegisterYourself();
  boltTrans3->RegisterYourself();
  boltTrans4->RegisterYourself();
  coneComb.Append("FMD3_bolt_hole:FMD3_bolt_matrix1"
		  "+FMD3_bolt_hole:FMD3_bolt_matrix2"
		  "+FMD3_bolt_hole:FMD3_bolt_matrix3"
		  "+FMD3_bolt_hole:FMD3_bolt_matrix4)");
  TGeoVolume*     boltVolume = new TGeoVolume("FMD3_bolt", boltShape, fSteel);
  support->AddNode(boltVolume, 1, boltTrans1);
  support->AddNode(boltVolume, 2, boltTrans2);
  boltShape->SetTitle("FMD3 steering bolt");
  boltVolume->SetTitle("FMD3 steering bolt");
  
  //__________________________________________________________________
  // Final cone
  TGeoCompositeShape* coneShape = new TGeoCompositeShape(coneComb.Data());
  coneShape->SetName("FMD3_cone");
  coneShape->SetTitle("FMD3 cone");
  TGeoVolume*  coneVolume = new TGeoVolume("FMD3_Cone", coneShape, fC);
  coneVolume->SetLineColor(kRed);
  support->AddNode(coneVolume, 0, new TGeoTranslation(0, 0, 0));

  //__________________________________________________________________
  // Tension boxes. 
  TGeoBBox* tensionOuter = new TGeoBBox("FMD3_tension_outer", .5, 3, 5);
  new TGeoBBox("FMD3_tension_inner", .51, 2.5, 4.6);
  TString tensionExpr("FMD3_tension_outer-FMD3_tension_inner");
  TGeoCompositeShape* tensionShape = new TGeoCompositeShape(tensionExpr.Data());
  tensionShape->SetName("FMD3_tension_box");
  tensionShape->SetTitle("FMD3 tension box");
  TGeoVolume*      tensionFrame = new TGeoVolume("FMD3_tension_frame",
						  tensionShape, fAl);
  TGeoTube*        springShape  = new TGeoTube("FMD3_tension_spring", 
					       0, .3, 4.3/2);
  TGeoVolume*      springVolume = new TGeoVolume("FMD3_tension_spring", 
						 springShape, fSteel);
  TGeoVolume*      tensionBox   = new TGeoVolume("FMD3_tension_box", 
						 tensionOuter, fAir);
  tensionBox->AddNode(tensionFrame, 0);
  tensionBox->AddNode(springVolume, 0, new TGeoTranslation(0,0,4.2/2));
  
  Double_t         tensionD     = 5*TMath::Cos(fmd3->GetConeOuterAngle());
  Double_t         tensionZ     = (r4->Z() - 2 * tensionD - 5 -
				   2*.5*TMath::Cos(fmd3->GetConeOuterAngle()));
  Double_t         tensionX     = (fmd3->ConeR(fmd3->GetInnerZ()
					       +fmd3->GetNoseZ() 
					       -tensionZ) + 
				   2*.5*TMath::Cos(fmd3->GetConeOuterAngle())); 
  TGeoRotation*    tensionRot   = new TGeoRotation();
  tensionRot->RotateY(180/TMath::Pi()*fmd3->GetConeOuterAngle());
  TGeoCombiTrans*  tensionBase  = new TGeoCombiTrans(tensionX, 0, tensionZ, 
						     tensionRot);
  
  Double_t         wireT        = .1;
  Double_t         wireR1       = fmd3->ConeR(fmd3->GetInnerZ()
					      +fmd3->GetNoseZ()) + wireT;
  Double_t         wireR2       = fmd3->ConeR(fmd3->GetInnerZ()
					      +fmd3->GetNoseZ()-
					      tensionZ+tensionD) + wireT;
  Double_t         wireL        = TMath::Sqrt(TMath::Power(wireR1-wireR2,2)+ 
					      TMath::Power(tensionZ-
							   tensionD,2));
  Double_t         wireAngle    = TMath::ATan2(wireR2-wireR1,tensionZ-tensionD);
  TGeoTube*        wireShape    = new TGeoTube("FMD3_wire", 0, wireT, wireL/2);
  TGeoVolume*      wireVolume   = new TGeoVolume("FMD3_wire", wireShape,fSteel);
  TGeoRotation*    wireRot      = new TGeoRotation();
  wireRot->RotateY(180/TMath::Pi()*wireAngle);
  TGeoCombiTrans*  wireBase     = new TGeoCombiTrans((wireR2-wireR1)/2+wireR1
						     +.1*TMath::Cos(wireAngle),
						     0,(tensionZ-tensionD)/2,
						     wireRot);
  for (Int_t i = 0; i < 2; i++) { 
    Double_t        thisAngle = (i+.5) * 90;
    TGeoCombiTrans* thisTrans = new TGeoCombiTrans(*tensionBase);
    thisTrans->RotateZ(thisAngle);
    support->AddNode(tensionBox, i, thisTrans);
    thisTrans = new TGeoCombiTrans(*wireBase);
    thisTrans->RotateZ(thisAngle);
    support->AddNode(wireVolume, i, thisTrans);
  }

  //__________________________________________________________________
  // Place support volumes in half-detector volumes 
  Double_t         z  = fmd3->GetInnerZ();
  TGeoTranslation* t1 = new TGeoTranslation(0, 0, -fmd3->GetNoseZ());
  fmd3TopVolume->AddNode(support, 1, t1);
  TGeoCombiTrans*  t2 = new TGeoCombiTrans(*t1);
  t2->RotateZ(180);
  fmd3BotVolume->AddNode(support, 2, t2);

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

#else
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
  // Double_t backr1  = fmd3->GetBackLowR();
  Double_t backr2  = fmd3->GetBackHighR();
  Double_t zdist   = conel -  backl - nlen;
  Double_t tdist   = backr2 - noser2;
  // Double_t beaml   = TMath::Sqrt(zdist * zdist + tdist * tdist);
  Double_t theta   = -180. * TMath::ATan2(tdist, zdist) / TMath::Pi();
  Double_t flanger = fmd3->GetFlangeR();
  Double_t z       = fmd3->GetInnerZ(); // fmd3->GetZ();

  TGeoVolume* fmd3TopVolume = new TGeoVolumeAssembly(Form(fgkFMDName, 
							  fmd3->GetId(), 'T'));
  TGeoVolume* fmd3BotVolume = new TGeoVolumeAssembly(Form(fgkFMDName, 
							  fmd3->GetId(), 'B'));
  fmd3TopVolume->SetTitle("FMD3 top half");
  fmd3BotVolume->SetTitle("FMD3 bottom half");
  
  
  DetectorGeometry(fmd3, fmd3TopVolume, fmd3BotVolume, z, 
		   innerTop, innerBot, outerTop, outerBot);

  
  TGeoVolumeAssembly* support = new TGeoVolumeAssembly("F3SU");
  support->SetTitle("FMD3 support");
  
  // Cone shape 
  TGeoPcon*        coneBase = new TGeoPcon("FMD3 cone base", 0, 180, 6);
  const TObjArray& radii    = fmd3.ConeRadii();
  TVector3*        v1       = 0;
  TVector3*        v4       = 0;
  for (Int_t i = 0; i < radii.GetEntriesFast(); i++) { 
    TVector3* v = static_cast<TVector3*>(radii.At(i));
    coneBase->DefineSection(i,  v->X(), v->Y(), v->Z());
    if (i == 1) v1 = v;
    if (i == 4) v4 = v;
    
  }
  Double_t  holeL = TMath::Sqrt(TMath::Power(v4->Z()-v1->Z(),2) + 
				TMath::Power(v4->X()-v1->X(),2));
  
  TGeoTrd1*       coneHole  = new TGeoTrd1("F3SC_hole",2,8,conet*3,
					   (conel-2-2)/2);
  


  // Nose volume 
  TGeoTubeSeg* noseShape  = new TGeoTubeSeg(noser1, noser2, nlen / 2, 0, 180);
  TGeoVolume*  noseVolume = new TGeoVolume(fgkNoseName, noseShape, fC);
  support->AddNode(noseVolume, 0, new TGeoTranslation(0, 0, nlen/2));
  noseShape->SetName(fgkNoseName);
  noseShape->SetTitle("FMD3 nose");
  noseVolume->SetTitle("FMD3 nose");
  
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
  boltShape->SetTitle("FMD3 steering bolt");
  boltVolume->SetTitle("FMD3 steering bolt");

  // Cooling plates
  TGeoTrd1*   plateShape  = new TGeoTrd1(2, 8, 0.1, (conel-2-2)/2-.1);
  TGeoVolume* plateVolume = new TGeoVolume("F3CO", plateShape, fAl);
  plateShape->SetName("F3C0");
  plateShape->SetTitle("FMD3 cooling plate");
  plateVolume->SetTitle("FMD3 cooling plate");

  // Shape for carbon half-cone
  TGeoConeSeg*    innerCone = new TGeoConeSeg("F3SC_inner", conel/2,
					      noser2-conet, noser2, 
					      backr2-conet, backr2, 0., 180.);
  innerCone->SetTitle("FMD3 cone inner");
  TGeoTrd1*       coneHole  = new TGeoTrd1("F3SC_hole",2,8,conet*3,
					   (conel-2-2)/2);
  coneHole->SetTitle("FMD3 cone hole");
  Double_t        holeAng   = TMath::ATan2(backr2 - noser2, conel);
  Double_t        holeX     = ((conel-2) / 2 * TMath::Sin(holeAng) +
			       conet     * TMath::Cos(holeAng) +
			       noser2);
  TGeoRotation*   holeRot   = new TGeoRotation();
  holeRot->SetName("FMD3 cone hole rotation");
  holeRot->RotateZ(90);
  holeRot->RotateY(holeAng*180./TMath::Pi());
  TGeoCombiTrans* holeTrans = new TGeoCombiTrans(holeX, 0, -2, holeRot);
  holeRot->SetName("FMD3 cone hole");

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
  coneShape->SetName("F3SC");
  coneShape->SetTitle("FMD3 cone");
  coneVolume->SetTitle("FMD3 cone");
  support->AddNode(coneVolume,1,new TGeoTranslation(0,0,nlen+conel/2));
  
  // The flanges 
  TGeoBBox* flangeShape    = new TGeoBBox((flanger - backr2) / 2, 
					  fmd3->GetBeamWidth() / 2,
					  backl / 2);
  TGeoVolume* flangeVolume = new TGeoVolume(Form(fgkFlangeName, fmd3->GetId()),
					    flangeShape, fC);
  flangeShape->SetName(Form(fgkFlangeName, fmd3->GetId()));
  flangeShape->SetTitle("FMD3 flange");
  flangeVolume->SetTitle("FMD3 flange");
  
  Int_t    n               = fmd3->GetNFlange();
  Double_t r               = backr2 + (flanger - backr2) / 2;
  for (Int_t i = 0; i  < n/2; i++) {
    Double_t phi       = 360. / n * i + 180. / n;
    Double_t x         = r * TMath::Cos(TMath::Pi() / 180 * phi);
    Double_t y         = r * TMath::Sin(TMath::Pi() / 180 * phi);
    TGeoRotation* rot  = new TGeoRotation;
    rot->RotateZ(phi);
    TGeoMatrix* matrix = new TGeoCombiTrans(x, y, nlen+conel-backl/2, rot);
    matrix->SetName(Form("FMD3_flange_%02d", i));
    matrix->SetTitle(Form("FMD3_flange_%2d", i));
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
#endif

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
