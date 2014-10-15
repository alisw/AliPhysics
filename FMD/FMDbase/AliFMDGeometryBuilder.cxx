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
#include <TGeoTorus.h>		// ROOT_TGeoTorus
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
TGeoShape*
AliFMDGeometryBuilder::MakeXTRU(const TObjArray& verticies, 
				Double_t thick) const
{
  // 
  // Make a polygonic extrusion shape based on verticies passed in @a
  // verticies 
  // 
  // Parameters:
  //    verticies List of verticies
  //    thick     Thickness
  // 
  // Return:
  //    newly allocated polygonic extrusion shape
  //
  TArrayD xs(6);
  TArrayD ys(6);
  for (Int_t i = 0; i < 3; i++) { 
    TVector2* v = static_cast<TVector2*>(verticies.At(i+1));
    xs[i]     =  v->Y();
    ys[i]     = -v->X();
    xs[6-1-i] =  v->Y();
    ys[6-1-i] =  v->X();
  }
  TGeoXtru* shape = new TGeoXtru(2);
  shape->DefinePolygon(xs.fN, xs.fArray, ys.fArray);
  shape->DefineSection(0, -thick/2);
  shape->DefineSection(1, +thick/2);
  
  return shape;
}

//____________________________________________________________________
TGeoVolume*
AliFMDGeometryBuilder::RingGeometry(const AliFMDRing* r) 
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
  Char_t        rng      = toupper(id);
  const Char_t* lName    = (rng == 'I' ? "inner" : "outer");
  Double_t siThick  = r->GetSiThickness();
  Double_t pcbThick = r->GetPrintboardThickness();
  Double_t cuThick  = r->GetCopperThickness();
  Double_t chipThick= r->GetChipThickness();
  Double_t modSpace = r->GetModuleSpacing();
  Double_t theta    = r->GetTheta();
  
  //------------------------------------------------------------------
  // Sensor
  // Physical sensor
  TGeoShape* sensorShape = MakeXTRU(r->GetSensorVerticies(), siThick);
  sensorShape->SetName(Form("FMD%c_physical_sensor", id));
  sensorShape->SetTitle(Form("FMD %s physical sensor", lName));
  TString sensorName = TString::Format(fgkSensorName, id);
  TGeoVolume* sensorVolume = new TGeoVolume(sensorName, sensorShape, fSi);
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
    TGeoTubeSeg* activeShape = new TGeoTubeSeg(r->GetLowR(),
					       r->GetHighR(), 
					       siThick / 2, 
					       - theta, 
					       + theta);
    activeShape->SetName(Form(fgkActiveName, id));
    activeShape->SetTitle(Form("FMD %s active area", lName));
    TString activeName = TString::Format(fgkActiveName, id);
    TGeoVolume* activeVolume = new TGeoVolume(activeName,activeShape,fSi);
    activeVolume->SetTitle(Form("FMD %s active area", lName));
    TString sectorName = TString::Format(fgkSectorName,id);
    TGeoVolume* sectorVolume = activeVolume->Divide(sectorName, 2, 2, -theta,
						    0,0,"N");

    Int_t    ns       = r->GetNStrips();
    Double_t stripoff = r->GetLowR(); // 0; // a->Mod();
    Double_t dstrip   = (r->GetHighR() - stripoff) / ns;

    sectorVolume->SetTitle(Form("FMD %s sector", lName));
    TString stripName = TString::Format(fgkStripName, id);
    TGeoVolume* stripVolume  = sectorVolume->Divide(stripName, 
						    1, ns, stripoff, dstrip, 
						    0, "SX");
    stripVolume->SetTitle(Form("FMD %s strip", lName));
    sid = stripVolume->GetNumber();
    sensorVolume->AddNodeOverlap(activeVolume, 0);
  }
  
  switch (rng) {
  case 'I': fActiveId[0] = sid; break;
  case 'O': fActiveId[1] = sid; break;
  }

  //------------------------------------------------------------------
  // Hybrid
  // PCB layer of hybrid 
  TGeoShape* pcbShape = MakeXTRU(r->GetHybridVerticies(), pcbThick);
  pcbShape->SetName(Form("FMD%c_hybrid_pcb", id));
  pcbShape->SetTitle(Form("FMD %s hybrid PCB", lName));
  TString pcbName = TString::Format(fgkPCBName, id);
  TGeoVolume* pcbVolume = new TGeoVolume(pcbName, pcbShape, fPCB);
  pcbVolume->SetTitle(Form("FMD %s hybrid PCB", lName));

  // Copper layer
  TGeoShape* cuShape = MakeXTRU(r->GetHybridVerticies(), cuThick);
  cuShape->SetName(Form("FMD%c_hybrid_copper", id));
  cuShape->SetTitle(Form("FMD %s hybrid copper", lName));
  TString cuName = TString::Format(fgkCuName,id);
  TGeoVolume* cuVolume = new TGeoVolume(cuName,cuShape,fCopper);
  cuVolume->SetTitle(Form("FMD %s hybrid copper", lName));

  // Chip layer
  TGeoShape* chipShape = MakeXTRU(r->GetHybridVerticies(), chipThick);
  chipShape->SetName(Form("FMD%c_hybrid_chip", id));
  chipShape->SetTitle(Form("FMD %s hybrid chip", lName));
  TString chipName = TString::Format(fgkChipName,id);
  TGeoVolume* chipVolume = new TGeoVolume(chipName,chipShape,fChip);
  chipVolume->SetTitle(Form("FMD %s hybrid chip", lName));

  //------------------------------------------------------------------
  // Legs
  Double_t      legr     = r->GetLegRadius();
  Double_t      legl     = r->GetLegLength();
  Double_t      lege     = .05;

  // Short leg shape 
  TGeoTube*   shortLegShape  = new TGeoTube(0, legr, (legl-lege) / 2);
  shortLegShape->SetName(Form(fgkShortLegName, id));
  shortLegShape->SetTitle(Form("FMD %s short support foot", lName));
  TString shortLegName = TString::Format(fgkShortLegName, id);
  TGeoVolume* shortLegVolume = new TGeoVolume(shortLegName, 
					      shortLegShape, fCopper);
  shortLegVolume->SetTitle(Form("FMD %s short support foot", lName));
  // Long leg shape
  TGeoTube*   longLegShape   = new TGeoTube(0, legr, 
					    (legl - lege + modSpace) / 2);
  longLegShape->SetName(Form(fgkLongLegName, id));
  longLegShape->SetTitle(Form("FMD %s long support foot", lName));
  TString longLegName = TString::Format(fgkLongLegName, id);
  TGeoVolume* longLegVolume = new TGeoVolume(longLegName, 
					      longLegShape, fCopper);
  longLegVolume->SetTitle(Form("FMD %s long support foot", lName));
  


  //------------------------------------------------------------------
  // Placement of module volumes in assemblies 
  TArrayD xfs(3);
  TArrayD yfs(3);
  for (Int_t i = 0; i < 3; i++) { 
    TVector2* vv = r->GetFootPosition(i);
    // TVector2  uu = vv->Rotate(TMath::Pi()/2);
    xfs[i]       = vv->Y();
    yfs[i]       = vv->X();
  }

  // Back container volume 
  TGeoVolume* backVolume     = new TGeoVolumeAssembly(Form(fgkBackVName, id));
  backVolume->SetTitle(Form("FMD %s back module", lName));
  TGeoVolume* frontVolume    = new TGeoVolumeAssembly(Form(fgkFrontVName, id));
  frontVolume->SetTitle(Form("FMD %s front module", lName));

  Double_t space    = r->GetSpacing();
  Double_t x        = 0;
  Double_t y        = 0;
  Double_t zb       = siThick / 2;
  Double_t zf       = siThick / 2;
  backVolume->AddNode(sensorVolume, 0, new TGeoTranslation(x, y, zb));
  frontVolume->AddNode(sensorVolume, 0, new TGeoTranslation(x, y, zf));
  zb         += siThick / 2 + space + pcbThick / 2;
  zf         += siThick / 2 + space + pcbThick / 2;
  backVolume->AddNode(pcbVolume, 0, new TGeoTranslation(x, y, zb));
  frontVolume->AddNode(pcbVolume, 0, new TGeoTranslation(x, y, zf));
  zb         += (pcbThick + cuThick) / 2;
  zf         += (pcbThick + cuThick) / 2;
  backVolume->AddNode(cuVolume, 0, new TGeoTranslation(0, 0, zf));
  frontVolume->AddNode(cuVolume, 0, new TGeoTranslation(0, 0, zb));
  zb         += (cuThick + chipThick) / 2;
  zf         += (cuThick + chipThick) / 2;
  backVolume->AddNode(chipVolume, 0, new TGeoTranslation(0, 0, zb));
  frontVolume->AddNode(chipVolume, 0, new TGeoTranslation(0, 0, zf));
  zb         += pcbThick / 2 + (legl)/ 2  - lege;
  zf         += pcbThick / 2 + (legl + modSpace)/ 2 - lege;
  for (Int_t i = 0; i < 3; i++) { 
    x          =  xfs[i]; // a->X() + legoff + legr;
    y          =  yfs[i]; // 0;
    backVolume->AddNode(shortLegVolume, i, new TGeoTranslation(x,y,zb));
    frontVolume->AddNode(longLegVolume, i, new TGeoTranslation(x,y,zf));
  }

  //------------------------------------------------------------------
  // FMDD 
  Double_t ddlr = r->GetFMDDLowR();
  Double_t ddhr = r->GetFMDDHighR();
  Double_t ddpt = r->GetFMDDPrintboardThickness();
  Double_t ddct = r->GetFMDDCopperThickness();
  Double_t ddit = r->GetFMDDChipThickness();
  Double_t ddt  = ddpt + ddct + ddit;
  
  TString    pcbdName(Form(fgkFMDDPCBName, id));
  TString    cudName(Form(fgkFMDDCuName, id));
  TString    chipdName(Form(fgkFMDDChipName, id));
  new TGeoTubeSeg(Form("%s_inner", pcbdName.Data()),  ddlr, ddhr, ddpt/2,0,180);
  new TGeoTubeSeg(Form("%s_inner", cudName.Data()),   ddlr, ddhr, ddct/2,0,180);
  new TGeoTubeSeg(Form("%s_inner", chipdName.Data()), ddlr, ddhr, ddit/2,0,180);
  
  Double_t clipWX = 0;
  Double_t clipWY = 0;
  Double_t clipY  = 1;
  
  if (rng == 'I') { 
    clipWX = ddhr;
    clipWY = ddhr/2;
  }
  else { 
    clipWX = ddlr+3;
    clipWY = ddhr/2;
  }
  
  new TGeoBBox(Form("%s_clip",  pcbdName.Data()), clipWX, clipWY, ddpt);
  new TGeoBBox(Form("%s_clip",  cudName.Data()),  clipWX, clipWY, ddct);
  new TGeoBBox(Form("%s_clip",  chipdName.Data()),clipWX, clipWY, ddit);
  TGeoTranslation* trans = new TGeoTranslation(Form("%s_trans",
						    pcbdName.Data()), 
					       0, clipWY+clipY, 0);
  trans->RegisterYourself();
  TGeoShape* fmddPcbShape = 
    new TGeoCompositeShape(pcbdName.Data(), 
			   Form("%s_inner*%s_clip:%s_trans",
				pcbdName.Data(), 
				pcbdName.Data(), 
				pcbdName.Data())); 
  TGeoShape* fmddCuShape = 
    new TGeoCompositeShape(cudName.Data(), 
			   Form("%s_inner*%s_clip:%s_trans",
				cudName.Data(), 
				cudName.Data(), 
				pcbdName.Data()));
  TGeoShape* fmddChipShape = 
    new TGeoCompositeShape(chipdName.Data(), 
			   Form("%s_inner*%s_clip:%s_trans",
				chipdName.Data(), 
				chipdName.Data(), 
				pcbdName.Data()));
  fmddPcbShape->SetTitle(Form("FMD %s digitiser PCB", lName));
  fmddCuShape->SetTitle(Form("FMD %s digitiser copper", lName));
  fmddChipShape->SetTitle(Form("FMD %s digitiser chip", lName));

  TString fmddPcbName = TString::Format(fgkFMDDPCBName, id);
  TGeoVolume* fmddPcbVolume = new TGeoVolume(fmddPcbName,
					      fmddPcbShape, fPCB);
  TString fmddCuName = TString::Format(fgkFMDDCuName, id);
  TGeoVolume* fmddCuVolume = new TGeoVolume(fmddCuName,
					      fmddCuShape, fCopper);
  TString fmddChipName = TString::Format(fgkFMDDChipName, id);
  TGeoVolume* fmddChipVolume = new TGeoVolume(fmddChipName,
					      fmddChipShape, fChip);
  fmddPcbVolume->SetTitle(Form("FMD %s digitiser PCB", lName));
  fmddCuVolume->SetTitle(Form("FMD %s digitiser copper", lName));
  fmddChipVolume->SetTitle(Form("FMD %s digitiser chip", lName));

  //------------------------------------------------------------------
  // Half ring mother volumes. 
  TGeoVolume* ringTopVolume = new TGeoVolumeAssembly(Form(fgkRingTopName,id));
  TGeoVolume* ringBotVolume = new TGeoVolumeAssembly(Form(fgkRingBotName,id));
  TGeoVolume* halfRing      = ringTopVolume;
  ringTopVolume->SetTitle(Form("FMD %s top half-ring", lName));
  ringBotVolume->SetTitle(Form("FMD %s bottom half-ring", lName));
  
  //------------------------------------------------------------------
  // Adding modules to half-rings
  Int_t    nmod =  r->GetNModules();
  AliFMDDebug(10, ("making %d modules in ring %c", nmod, id));
  for (Int_t i = 0; i < nmod; i++) {
    if (i == nmod / 2) halfRing = ringBotVolume;
    Bool_t      front =  (i % 2 == (rng == 'I' ? 1 : 0));
    TGeoVolume* vol   =  (front ? frontVolume : backVolume);
    // vol->AddNode(sensorVolume, i, new TGeoTranslation(0,0,siThick/2));
    Double_t    z1    =  (front ? -1 : 1) * modSpace / 2;
    // Double_t z1    =  (front ? 0 : modSpace);
    Double_t    th    =  (2 * i + 1) * theta;
    TGeoMatrix* mat1  =  new TGeoCombiTrans(0,0,z1,0); 
    mat1->RotateZ(th);
    mat1->SetName(Form("FMD%c_module_%02d", id, i));
    mat1->SetTitle(Form("FMD %s module %2d matrix", lName, i));
    halfRing->AddNode(vol, i, mat1);
  }

  //------------------------------------------------------------------
  // Add the FMDD 
  Double_t zi = r->GetFullDepth() - ddt;
  Int_t    n  = 2;
  for (Int_t i = 0; i  < n; i++) {
    halfRing             = (i == 0 ? ringTopVolume : ringBotVolume);
    Double_t      phi    = 360. / n * i;
    TGeoRotation* rot    = new TGeoRotation(Form("FMDD%c rotation %d", id, i));
    rot->RotateZ(phi);
    rot->SetTitle(Form("FMD %s digitiser rotation %2d", lName, i));
    Double_t z =  zi + ddpt / 2;
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
AliFMDGeometryBuilder::TensionBox()
{
  // 
  // Get the tension box volume
  // 
  // 
  // Return:
  //    
  //
  static TGeoVolumeAssembly* tensionBox = 0;
  if (tensionBox) return tensionBox;
  
  TGeoBBox* tensionEndS = new TGeoBBox("FMD_tension_end", .6, 3,  .25);
  TGeoBBox* tensionTopS = new TGeoBBox("FMD_tension_top", .1, .5, 3.5);
  TGeoVolume* tensionEndV = new TGeoVolume("FMD_tension_end", tensionEndS,fAl);
  TGeoVolume* tensionTopV = new TGeoVolume("FMD_tension_top", tensionTopS,fAl);
  tensionBox = new TGeoVolumeAssembly("FMD_tension_box");
  tensionBox->AddNode(tensionEndV, 1, new TGeoTranslation(.6, 0,   -3.75));
  tensionBox->AddNode(tensionEndV, 2, new TGeoTranslation(.6, 0,   +3.75));
  tensionBox->AddNode(tensionTopV, 1, new TGeoTranslation(0.1, +2.5, 0));
  tensionBox->AddNode(tensionTopV, 2, new TGeoTranslation(0.1, -2.5, 0));
  tensionBox->AddNode(tensionTopV, 3, new TGeoTranslation(1.1, +2.5, 0));
  tensionBox->AddNode(tensionTopV, 4, new TGeoTranslation(1.1, -2.5, 0));
  return tensionBox;
}


//____________________________________________________________________
TGeoVolume*
AliFMDGeometryBuilder::DetectorGeometry(const AliFMDDetector* d, 
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
    AliFMDDebug(2, ("Placing volumes %s and %s in %s and %s at z=%f", 
		     tvol->GetName(), bvol->GetName(), 
		     topMother->GetName(), botMother->GetName(), z));
    topMother->AddNode(tvol, Int_t(c), new TGeoTranslation(0,0,z));
    botMother->AddNode(bvol, Int_t(c), new TGeoTranslation(0,0,z));

    // Honeycomp 
    TGeoShape*   hcSha = HoneycombShape(id, c, lowr, highr, hcThick, alThick);
    TGeoVolume*  hcVol = new TGeoVolume(Form(fgkHCName,id,c),hcSha,fAl);
    hcVol->SetTitle(Form("FMD%d%c honeycomb shell", id, c));
    
    z += (r->GetModuleDepth() 
	  + r->GetModuleSpacing() / 2
	  + r->GetHoneycombThickness() / 2);

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
AliFMDGeometryBuilder::FMD1Geometry(const AliFMD1* fmd1, 
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
  
  // `Top' or `Outside' master volume
  TString fmd1TopName = TString::Format(fgkFMDName, fmd1->GetId(), 'T');
  TGeoVolume* fmd1TopVolume = new TGeoVolumeAssembly(fmd1TopName);
  fmd1TopVolume->SetTitle("FMD1 top half");

  // `Bottom' or `Inside' master volume
  TString fmd1BotName = TString::Format(fgkFMDName, fmd1->GetId(), 'B');
  TGeoVolume* fmd1BotVolume = new TGeoVolumeAssembly(fmd1BotName);
  fmd1BotVolume->SetTitle("FMD1 bottom half");
  
  // Basic detector geometry 
  DetectorGeometry(fmd1, fmd1TopVolume, fmd1BotVolume, z, 
		   innerTop, innerBot, 0, 0);

  Double_t lidP[][3] = { {  0.00,  4.20, 20.95 }, 
			 {  0.15,  4.20, 20.95 }, 
			 {  0.15, 20.80, 20.95 }, 
			 {  3.00, 20.80, 20.95 }, 
			 {  3.00, 20.80, 22.30 }, 
			 {  3.15, 20.80, 22.30 }, 
			 {  3.15, 20.95, 24.65 },
			 {  3.30, 20.95, 24.65 }, 
			 {  3.30, 24.50, 24.65 }, 
			 {  6.80, 24.50, 24.65 },
			 {  6.80, 24.50, 26.00 },
			 {  6.95, 24.50, 26.00 } };
  Double_t  lidZStart = lidP[11][0];
  TGeoPcon* lidBaseS  = new TGeoPcon("FMD1_lid_base", 0, 180, 12);
  for (size_t i = 0; i < 12; i++) 
    lidBaseS->DefineSection(i, lidP[i][0] - lidZStart, lidP[i][1], lidP[i][2]);
  
  
  Double_t lidH[][2] = { {  7.84903, 24.15680  }, 
			 { 20.54900, 14.92970  },
			 { 21.99700, 12.70000  },
			 { 25.26090,  2.65502  } };
  Double_t lidHR = .53 / 2;
  Double_t lidHL = 0.16;
  
  new TGeoTube("FMD1_lid_hole", 0, lidHR, lidHL/2);
  TString lidComp("FMD1_lid_base-(");
  TGeoTranslation* trans = 0;
  for (size_t i = 0; i < 4; i++) { 
    trans = new TGeoTranslation(-lidH[i][0], lidH[i][1], /*6.95*/-lidHL/2);
    trans->SetName(Form("FMD1_lid_hole_mat%d", int(2*i+0)));
    trans->RegisterYourself();
    trans = new TGeoTranslation(+lidH[i][0], lidH[i][1], /*6.95*/-lidHL/2);
    trans->SetName(Form("FMD1_lid_hole_mat%d", int(2*i+1)));
    trans->RegisterYourself();
    lidComp.Append(Form("FMD1_lid_hole:FMD1_lid_hole_mat%d+" 
			"FMD1_lid_hole:FMD1_lid_hole_mat%d%c", 
			int(2 * i), int(2 * i + 1), int(i == 3 ? ')' : '+')));
  }
  TGeoCompositeShape* lidS = new TGeoCompositeShape(lidComp.Data());
  lidS->SetName("FMD1_lid");
  TGeoVolume* lidV = new TGeoVolume("FMD1_lid", lidS, fC);
  lidV->SetTransparency(63);
  
  // Place top cover
  Double_t lidZ = (lidZStart - 
		   (3.3 - r->GetModuleDepth() - r->GetModuleSpacing() / 2));
  AliFMDDebug(1, ("FMD1 lid offset in Z=%f", lidZ));

  for (Int_t i = 0; i  < 2; i++) {
    TGeoVolume*   mother = (i == 0 ? fmd1TopVolume : fmd1BotVolume);
    Double_t      phi    = 360. / 2 * i;
    TGeoRotation* rot    = new TGeoRotation(Form("FMD1_lid_rot%d",i));
    rot->RotateZ(phi);
    TGeoMatrix* matrix   = new TGeoCombiTrans(Form("FMD1_lid_mat%d", i),
					      0, 0, lidZ, rot);
    mother->AddNode(lidV, i, matrix);    
  }

  // Must add this after filling the assembly.
  TGeoVolume* top    = gGeoManager->GetVolume("ALIC");
  // TGeoMatrix* matrix = new TGeoTranslation("FMD1 trans", 0, 0, z);
  TGeoRotation* rot = new TGeoRotation("FMD1 rotatation");
  rot->RotateZ(90);
  TGeoMatrix* matrix = new TGeoCombiTrans("FMD1 trans", 0, 0, z, rot);

  AliFMDDebug(5, ("Placing volumes %s and %s in ALIC at z=%f", 
		   fmd1TopVolume->GetName(), fmd1BotVolume->GetName(), z));
  top->AddNode(fmd1TopVolume, fmd1->GetId(), matrix);
  top->AddNode(fmd1BotVolume, fmd1->GetId(), matrix);


  // Survey points on V0A (screw holes for the FMD) 
  const Double_t icb[] = { +12.700, -21.997, 324.670 };
  const Double_t ict[] = { +12.700, +21.997, 324.670 };
  const Double_t ocb[] = { -12.700, -21.997, 324.670 };
  const Double_t oct[] = { -12.700, +21.997, 324.670 };

  TGeoTube* surveyShape = new TGeoTube("FMD1_survey_marker", 
					0, .2, .001);

  TGeoMatrix* outMat = matrix;
#if 0
  if (gGeoManager->cd("/ALIC_1/F1MT_1")) 
    outMat = gGeoManager->GetCurrentMatrix();
  else 
    AliWarning("Couldn't cd to /ALIC_1/F1MT_1");
#endif

  Double_t loct[3], locb[3];
  outMat->MasterToLocal(oct, loct);
  outMat->MasterToLocal(ocb, locb);
  TGeoVolume* vOct = new TGeoVolume("V0L_OCT", surveyShape, fPlastic);
  TGeoVolume* vOcb = new TGeoVolume("V0L_OCB", surveyShape, fPlastic);
  
  fmd1TopVolume->AddNode(vOct, 1, new TGeoTranslation(loct[0],loct[1],loct[2]));
  fmd1TopVolume->AddNode(vOcb, 1, new TGeoTranslation(locb[0],locb[1],locb[2]));
    
  
  TGeoMatrix* inMat = matrix;
#if 0
  if (gGeoManager->cd("/ALIC_1/F1MT_1")) 
    inMat = gGeoManager->GetCurrentMatrix();
  else 
    AliWarning("Couldn't cd to /ALIC_1/F1MT_1");
#endif

  Double_t lict[3], licb[3];
  inMat->MasterToLocal(ict, lict);
  inMat->MasterToLocal(icb, licb);
  TGeoVolume* vIct = new TGeoVolume("V0L_ICT", surveyShape, fPlastic);
  TGeoVolume* vIcb = new TGeoVolume("V0L_ICB", surveyShape, fPlastic);
  
  fmd1BotVolume->AddNode(vIct, 1, new TGeoTranslation(lict[0],lict[1],lict[2]));
  fmd1BotVolume->AddNode(vIcb, 1, new TGeoTranslation(licb[0],licb[1],licb[2]));

  return 0;
}

//____________________________________________________________________
TGeoVolume*
AliFMDGeometryBuilder::FMD2Geometry(const AliFMD2* fmd2, 
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
  AliFMDRing* ring          = fmd2->GetOuter();
  Double_t    z             = fmd2->GetOuterZ();  
  Double_t    framelr       = 32.01;  // fmd2->GetOuterHoneyHighR()+0.5;
  Double_t    framehr       = 33.611; // fmd2->GetOuterHoneyHighR()+1.8;
  Double_t    framel        = 14.8; // framehz - framelz;
  // Double_t    backth        = 0.3;
  Double_t    backth        = 0.03;
  Double_t    framelz       = -(2.38 
				- ring->GetModuleDepth() 
				- ring->GetModuleSpacing() / 2);
  // Double_t    framelz       = -0.8;
  // Double_t    framehz       = framelz + backth + framel;
  Double_t    coverlr       = 4.3; // fmd2->GetInner()->GetLowR()+1;
  Double_t    coverhr       = framehr; //  - 1;
  
  TString fmd2TopName = TString::Format(fgkFMDName, fmd2->GetId(), 'T');
  TGeoVolume* fmd2TopVolume = new TGeoVolumeAssembly(fmd2TopName);
  TString fmd2BotName = TString::Format(fgkFMDName, fmd2->GetId(), 'B');
  TGeoVolume* fmd2BotVolume = new TGeoVolumeAssembly(fmd2BotName);
  fmd2TopVolume->SetTitle("FMD2 top half");
  fmd2BotVolume->SetTitle("FMD2 bottom half");
  
  DetectorGeometry(fmd2, fmd2TopVolume, fmd2BotVolume, z, 
		   innerTop, innerBot, outerTop, outerBot);

  TGeoVolumeAssembly* support = new TGeoVolumeAssembly("FMD2_support");
  TGeoShape*  cylinderShape   = new TGeoTubeSeg(framelr,framehr,framel/2,0,180);
  TGeoVolume* cylinderVolume  = new TGeoVolume(Form(fgkBackName, fmd2->GetId()),
					       cylinderShape, fC);
  TGeoShape*  coverShape      = new TGeoTubeSeg(coverlr,coverhr,backth/2,0,180);
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
  
  TGeoTranslation* trans = 0;
  support->AddNode(coverVolume,1, new TGeoTranslation(0,0,backth/2));
  support->AddNode(cylinderVolume, 1, new TGeoTranslation(0,0,backth+framel/2));
  

  Double_t    f1l            = 15.6085;
  Double_t    f1w            = 6;
  Double_t    f1d            = 1;
  Int_t       nFiducialHoles = 4;
  Double_t    precHoles[][2] = { { 32.4948, 29.6663 },
				 { 33.9104, 31.0819 },
				 { 34.8177, 33.4035 }, 
				 { 35.5028, 32.6744 } };
  Double_t    precRadius     = .25;
  Double_t    flangeA        = TMath::Pi()/4;
  
  new TGeoBBox("FMD2_flange_base", f1l/2, f1w/2, f1d/2);
  new TGeoTube("FMD2_fiducial_hole", 0, precRadius, f1d/2+.1);
  Double_t         flangeX        = framehr + f1l/2;
  TVector2         flangeC(flangeX * TMath::Cos(flangeA), 
			   flangeX * TMath::Sin(flangeA));
  TString          flangeComb("FMD2_flange_base-(");  
  new TGeoBBox("FMD2_flange_slit", 7./2, 1.5/2, f1d/2+.1);
  trans = new TGeoTranslation(-f1l/2+1+7./2, +.5+1.5/2, 0);
  trans->SetName("FMD2_flange_slit_mat1");
  trans->RegisterYourself();
  trans = new TGeoTranslation(-f1l/2+1+7./2, -.5-1.5/2, 0);
  trans->SetName("FMD2_flange_slit_mat2");
  trans->RegisterYourself();
  flangeComb.Append("FMD2_flange_slit:FMD2_flange_slit_mat1+"
		    "FMD2_flange_slit:FMD2_flange_slit_mat2+");
  for (Int_t i = 0; i < nFiducialHoles; i++) { 
    TVector2         v(precHoles[i][0], precHoles[i][1]);
    v                   -= flangeC;
    TVector2         r  =  v.Rotate(-flangeA);
    TGeoTranslation* t1 =  new TGeoTranslation(r.X(),  r.Y(), 0);
    TGeoTranslation* t2 =  new TGeoTranslation(r.X(), -r.Y(), 0);
    t1->SetName(Form("FMD2_fiducial_hole_rot%d", 2*i+0));
    t2->SetName(Form("FMD2_fiducial_hole_rot%d", 2*i+1));
    t1->RegisterYourself();
    t2->RegisterYourself();
    flangeComb.Append(Form("FMD2_fiducial_hole:FMD2_fiducial_hole_rot%d+"
			   "FMD2_fiducial_hole:FMD2_fiducial_hole_rot%d%c",
			   2*i+0, 2*i+1, (i == nFiducialHoles-1 ? ')' : '+')));
  }
  // Final flange shape, and at to full shape 
  TGeoCompositeShape* flangeS = new TGeoCompositeShape(flangeComb.Data());
  flangeS->SetName("FMD2_flange");
  TGeoVolume* flangeV = new TGeoVolume("FMD2_flange", flangeS, fAl);
  
  Double_t f2l = 7;
  Double_t f2d = 12.5;
  Double_t f2w = 1;

  new TGeoBBox("FMD2_flange_spacer_base", f2l/2, f2w/2, f2d/2);
  new TGeoTube("FMD2_flange_spacer_hole", 0, 2.5, f2w/2+.1);
  TGeoRotation* holeRot = new TGeoRotation();
  holeRot->RotateY(90);
  holeRot->RotateZ(90);
  TGeoCombiTrans* combo = 0;
  combo = new TGeoCombiTrans(0, 0, f2d/2-.5-2.5, holeRot);
  combo->SetName("FMD2_flange_spacer_hole_mat1");
  combo->RegisterYourself();
  combo = new TGeoCombiTrans(0, 0, -f2d/2+.5+2.5, holeRot);
  combo->SetName("FMD2_flange_spacer_hole_mat2");
  combo->RegisterYourself();
  TString spacerComp("FMD2_flange_spacer_base-("
		     "FMD2_flange_spacer_hole:FMD2_flange_spacer_hole_mat1+"
		     "FMD2_flange_spacer_hole:FMD2_flange_spacer_hole_mat2)");
  TGeoCompositeShape* spacerS = new TGeoCompositeShape(spacerComp.Data());
  TGeoVolume*         spacerV = new TGeoVolume("FMD2_flange_spacer",
					       spacerS, fAl);

  Double_t            extraL  = framehr-framelr;
  TGeoBBox*           extraS  = new TGeoBBox("FMD2_flange_extra", 
					     extraL/2, f1w/2, f1d/2);
  TGeoVolume*         extraV  = new TGeoVolume("FMD2_flange_extra", extraS,fAl);
  TGeoVolumeAssembly* wingV   = new TGeoVolumeAssembly("FMD2_wing");
  TGeoVolume*         tension = TensionBox();
  TGeoTube*           wireS   = new TGeoTube(0, .05, (framehr-coverlr)/2);
  TGeoVolume*         wireV   = new TGeoVolume("FMD2_tension_wire", 
					       wireS, fSteel);
  wingV->AddNode(flangeV, 1, new TGeoTranslation(f1l/2,    0, f1d/2));
  wingV->AddNode(flangeV, 2, new TGeoTranslation(f1l/2,    0, -f2d-f1d/2));
  wingV->AddNode(extraV, 1, new TGeoCombiTrans(-extraL/2, 0, f1d/2, 0));
  wingV->AddNode(spacerV, 1, new TGeoTranslation(1+f2l/2,-f2w/2+f1w/2,
						 -f2d/2));
  wingV->AddNode(spacerV, 2, new TGeoTranslation(1+f2l/2,+f2w/2-f1w/2,
						 -f2d/2));
  TGeoRotation* tensionR = new TGeoRotation;
  tensionR->RotateY(90);
  wingV->AddNode(tension, 1, new TGeoCombiTrans(4, 0, f1d+1.2, tensionR));
  TGeoRotation* wireR = new TGeoRotation;
  wireR->RotateY(90);
  wingV->AddNode(wireV, 1, new TGeoCombiTrans(-(framehr-coverlr)/2, 0, f1d+1,
					      wireR));
  
  TGeoCombiTrans* extraM1 = new TGeoCombiTrans(coverhr-extraL/2,0,0,0);
  extraM1->RotateZ(45);
  extraM1->RegisterYourself();
  extraM1->SetName("FMD2_back_cover_slit1");
  TGeoCombiTrans* extraM2 = new TGeoCombiTrans(coverhr-extraL/2,0,0,0);
  extraM2->RotateZ(135);
  extraM2->RegisterYourself();
  extraM2->SetName("FMD2_back_cover_slit2");
  TString coverComp(Form(fgkTopName, fmd2->GetId()));
  coverComp.Append("-(FMD2_flange_extra:FMD2_back_cover_slit1"
		   "+FMD2_flange_extra:FMD2_back_cover_slit2)");
  TGeoCompositeShape* cover2Shape = new TGeoCompositeShape(coverComp.Data());
  cover2Shape->SetName("FMD2_back_cover");
  TGeoVolume* cover2Volume = new TGeoVolume("FMD2_back_cover", cover2Shape,fC);
  support->AddNode(cover2Volume,2, 
		   new TGeoTranslation(0,0,backth+framel+backth/2));

  TGeoCombiTrans* trans1 = new TGeoCombiTrans(framehr, 0, backth+framel, 0);
  TGeoCombiTrans* trans2 = new TGeoCombiTrans(framehr, 0, backth+framel, 0);
  trans1->RotateZ(45);
  trans2->RotateZ(135);
  support->AddNode(wingV, 1, trans1);
  support->AddNode(wingV, 2, trans2);
  AliFMDDebug(1, ("FMD2 support offset is %f", framelz));
  
  for (Int_t i = 0; i  < 2; i++) {
    TGeoVolume*   mother = (i < 1 ? fmd2TopVolume : fmd2BotVolume);
    
    Double_t      phi    = 360. / 2 * i;
    TGeoRotation* rot    = new TGeoRotation(Form("FMD2 support rot %d",i)); 
    rot->RotateZ(phi);
    TGeoMatrix*   matrix = new TGeoCombiTrans(0, 0, framelz, rot);
    mother->AddNode(support, i, matrix);    
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
AliFMDGeometryBuilder::FMD3Geometry(const AliFMD3* fmd3, 
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
  TString fmd3TopName = TString::Format(fgkFMDName, fmd3->GetId(), 'T');
  TGeoVolume* fmd3TopVolume = new TGeoVolumeAssembly(fmd3TopName);
  TString fmd3BotName = TString::Format(fgkFMDName, fmd3->GetId(), 'B');
  TGeoVolume* fmd3BotVolume = new TGeoVolumeAssembly(fmd3BotName);
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
#if 0
  TGeoTube*        fiducialShape  = 
#endif
    new TGeoTube("FMD3_fiducial_hole", 0, fiducialRadius, flangeDepth+.1);
  Int_t            nFiducialHoles = fiducialHoles.GetEntriesFast();
  double           flangeAngle    = TMath::Pi() / 4;
  double           flangeX        = r5->Y()+flangeLength;
  TVector2         flangeC(flangeX * TMath::Cos(flangeAngle), 
			   flangeX * TMath::Sin(flangeAngle));
  TString          flangeComb("FMD3_flange_base-(");
#if 0// For debugging geometry 
  TGeoVolume* fiducialVolume = new TGeoVolume("FMD3_fiducial", fiducialShape);
  fiducialVolume->SetLineColor(kGreen);
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
#if 0 // For debugging geometry 
    support->AddNode(fiducialVolume, 2*i+0, t1);
    support->AddNode(fiducialVolume, 2*i+1, t2);
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
  Double_t holeL  = fmd3->GetHoleLength()/2;
  Double_t holeD  = fmd3->GetHoleDepth()/2;
  Double_t holeLW = fmd3->GetHoleLowWidth()/2;
  Double_t holeHW = fmd3->GetHoleHighWidth()/2;
  Double_t holeA  = fmd3->GetConeOuterAngle();
  Double_t holeA2 = TMath::Pi() - fmd3->GetConeOuterAngle();
  Double_t holeO  = fmd3->GetHoleOffset();
  Double_t holeZ  = (holeO
		     + holeL * TMath::Cos(holeA)
		     - holeD * TMath::Sin(holeA2));
  Double_t holeX  = (fmd3->ConeR(-holeZ + fmd3->GetInnerZ() + fmd3->GetNoseZ())
		     - holeD * TMath::Sin(holeA2));
  new TGeoTrd1("FMD3_cone_hole", holeLW, holeHW, holeD, holeL);
  TGeoTrd1* plateShape = new TGeoTrd1("FMD3_cooling_plate", 
				      holeLW, holeHW, .033, holeL);
  TGeoRotation* holeRot = new TGeoRotation();
  holeRot->SetName("FMD3_cone_hole_rotation");
  holeRot->RotateZ(90);
  holeRot->RotateY(holeA*180/TMath::Pi());
  TGeoCombiTrans* holeBaseTrans = new TGeoCombiTrans(holeX, 0, holeZ, holeRot);
  holeBaseTrans->SetName("FMD3_cone_hole_base_matrix");
  // TGeoRotation* plateRot = new TGeoRotation();
  // plateRot->SetName("FMD3_cone_plate_rotation");
  // plateRot->RotateZ(90);
  // plateRot->RotateY(plateA*180/TMath::Pi());
  // TGeoCombiTrans* plateBaseTrans = new 
  //                  TGeoCombiTrans(plateX,0,plateZ,plateRot);
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
    trans = new TGeoCombiTrans(*holeBaseTrans);
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
		  "+FMD3_bolt_hole:FMD3_bolt_matrix4");
  TGeoVolume*     boltVolume = new TGeoVolume("FMD3_bolt", boltShape, fSteel);
  support->AddNode(boltVolume, 1, boltTrans1);
  support->AddNode(boltVolume, 2, boltTrans2);
  boltShape->SetTitle("FMD3 steering bolt");
  boltVolume->SetTitle("FMD3 steering bolt");
  
  //__________________________________________________________________
  // Cut-outs for tension wheel sheeve 
  new TGeoBBox("FMD3_sheeve_hole", .55, .75, 1.16);
  Double_t        sheeveHoleZ = fmd3->GetInnerZ() + fmd3->GetNoseZ() - .75;
  Double_t        sheeveHoleR = fmd3->ConeR(sheeveHoleZ) - .55 + .2572222;
  TGeoCombiTrans* sheeveMat1  = new TGeoCombiTrans(sheeveHoleR,0,1.15,0);
  TGeoCombiTrans* sheeveMat2  = new TGeoCombiTrans(sheeveHoleR,0,1.15,0);
  sheeveMat1->RotateZ(45);
  sheeveMat2->RotateZ(135);
  sheeveMat1->SetName("FMD3_sheeve_hole_matrix1");
  sheeveMat2->SetName("FMD3_sheeve_hole_matrix2");
  sheeveMat1->RegisterYourself();
  sheeveMat2->RegisterYourself();
  coneComb.Append("+FMD3_sheeve_hole:FMD3_sheeve_hole_matrix1"
		  "+FMD3_sheeve_hole:FMD3_sheeve_hole_matrix2)");
  
  //__________________________________________________________________
  // Sheeve boxes 
  Double_t       sheeveL     = 1.15;
  TGeoBBox*      sheeveSideS = new TGeoBBox("FMD3_sheeve_side",
					   .55, .25, 1.15);
  TGeoBBox*      sheeveBackS = new TGeoBBox("FMD3_sheeve_back", 
					    .55, .25, .15);
  TGeoBBox*      sheeveWingS = new TGeoBBox("FMD3_sheeve_wing", 
					    .15, .15, 1.15);
  TGeoPcon*      sheeveWheelS = new TGeoPcon("FMD3_sheeve_wheel", 0, 360, 9);
  Double_t       sheeveInnerR = 0; // .2;
  Double_t       sheeveR      = .875;
  Double_t       sheeveWheelZ = .95;
  sheeveWheelS->DefineSection(0, -.25,   sheeveInnerR, 1);
  sheeveWheelS->DefineSection(1, -.125,  sheeveInnerR, 1);
  sheeveWheelS->DefineSection(2, -.125,  sheeveInnerR, sheeveWheelZ);
  sheeveWheelS->DefineSection(3, -.0625, sheeveInnerR, sheeveR+.02);
  sheeveWheelS->DefineSection(4, 0.000,  sheeveInnerR, sheeveR);
  sheeveWheelS->DefineSection(5, +.0625, sheeveInnerR, sheeveR+.02);
  sheeveWheelS->DefineSection(6, +.125,  sheeveInnerR, sheeveWheelZ);
  sheeveWheelS->DefineSection(7, +.125,  sheeveInnerR, 1);
  sheeveWheelS->DefineSection(8, +.25,   sheeveInnerR, 1);
  TGeoVolume*    sheeveSideV = new TGeoVolume("FMD3_sheeve_side", 
					      sheeveSideS, fPlastic);
  TGeoVolume*    sheeveBackV = new TGeoVolume("FMD3_sheeve_back", 
					      sheeveBackS, fPlastic);
  TGeoVolume*    sheeveWingV = new TGeoVolume("FMD3_sheeve_wing", 
					      sheeveWingS, fPlastic);
  TGeoVolume*    sheeveWheelV= new TGeoVolume("FMD3_sheeve_wheel", 
					      sheeveWheelS, fPlastic);
  TGeoVolumeAssembly* sheeveBox = new TGeoVolumeAssembly("FMD3_sheeve_box");
  sheeveBox->AddNode(sheeveSideV, 1, new TGeoTranslation(0, -.5, 0));
  sheeveBox->AddNode(sheeveSideV, 2, new TGeoTranslation(0, +.5, 0));
  sheeveBox->AddNode(sheeveBackV, 1, new TGeoTranslation(0, 0, 2.0+.15-1.15));
  sheeveBox->AddNode(sheeveWingV, 1, new TGeoTranslation(.55-.15, -.90, 0));
  sheeveBox->AddNode(sheeveWingV, 2, new TGeoTranslation(.55-.15, +.90, 0));
  TGeoRotation*   sheeveWheelR = new TGeoRotation;
  sheeveWheelR->RotateX(90);
  TGeoCombiTrans* sheeveWheelM = new TGeoCombiTrans(0, 0, sheeveWheelZ-sheeveL,
						    sheeveWheelR);
  sheeveBox->AddNode(sheeveWheelV, 1, sheeveWheelM);
  support->AddNode(sheeveBox, 1, sheeveMat1);
  support->AddNode(sheeveBox, 2, sheeveMat2);
  
  

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
  TGeoVolume*     tensionBox = TensionBox();
  Double_t        tensionH  = .6;
  Double_t        tensionL  = 4;
  Double_t        tensionZ  = 23.654;
  Double_t        tensionR  = fmd3->ConeR(fmd3->GetInnerZ() + fmd3->GetNoseZ() 
					  -  tensionZ);
  Double_t 	  tensionAr = fmd3->GetConeOuterAngle();
  Double_t 	  tensionA  = tensionAr * 180 / TMath::Pi();
  TGeoRotation*   tensionQ  = new TGeoRotation;
  tensionQ->RotateY(tensionA);
  TGeoCombiTrans* tensionM1 = new TGeoCombiTrans(tensionR,0,tensionZ, tensionQ);
  TGeoCombiTrans* tensionM2 = new TGeoCombiTrans(tensionR,0,tensionZ, tensionQ);
  tensionM1->RotateZ(45);
  tensionM2->RotateZ(135);
  support->AddNode(tensionBox, 1, tensionM1);
  support->AddNode(tensionBox, 2, tensionM2);
  
  // Double_t         tensionHR    = 0.15;
  Double_t         wireT        = .1/2;
  Double_t         wireZ1       = (tensionZ
				   - tensionL * TMath::Cos(tensionAr) 
				   - tensionH * TMath::Sin(tensionAr));
  Double_t         wireR1       = (tensionR 
				   - tensionL * TMath::Sin(tensionAr) 
				   + tensionH * TMath::Cos(tensionAr));
  AliFMDDebug(10, ("Wire Z1: %f=%f-%f*cos(%f)-%f*sin(%f)", 
		  wireZ1, tensionZ, tensionL, tensionAr, tensionH, tensionAr));
  AliFMDDebug(10, ("Wire R1: %f=%f-%f*sin(%f)-%f*cos(%f)", 
		  wireR1, tensionR, tensionL, tensionAr, tensionH, tensionAr));
  
  Double_t         wireStartA   = 42.3 * TMath::Pi() / 180;
  Double_t         wireZ2       = (sheeveWheelZ * (1 - TMath::Sin(wireStartA))
				   // - sheeveL - 
				   - wireT * TMath::Sin(wireStartA));
  /* (sheeveWheelZ * (1 - TMath::Sin(wireStartA))
				   - wireT * TMath::Sin(wireStartA) 
				   - sheeveL); */
  AliFMDDebug(10, ("wireZ2=%f=%f*(1-%f)", wireZ2, sheeveWheelZ, 
		  TMath::Sin(wireStartA)));
  Double_t         wireR2       = (sheeveHoleR + 
				   sheeveWheelZ * TMath::Cos(wireStartA) + 
				   wireT * TMath::Cos(wireStartA));
  Double_t         wireDR       = wireR1-wireR2;
  Double_t         wireDZ       = wireZ1-wireZ2;
  Double_t         wireL        = TMath::Sqrt(wireDR*wireDR+wireDZ*wireDZ)-.01;
  Double_t         wireAngle    = TMath::ATan2(wireDR,wireDZ);
  TGeoTube*        wireShape    = new TGeoTube("FMD3_wire", 0, wireT, wireL/2);
  TGeoVolume*      wireVolume   = new TGeoVolume("FMD3_wire", wireShape,fSteel);
  TGeoRotation*    wireRot      = new TGeoRotation();
  wireRot->RotateY(180/TMath::Pi()*wireAngle);
  Double_t         wireR        = wireR2 + wireDR / 2;
  Double_t         wireZ        = wireZ2 + wireDZ / 2;
  TGeoCombiTrans*  wireM1       = new TGeoCombiTrans(wireR, 0,wireZ, wireRot);
  TGeoCombiTrans*  wireM2       = new TGeoCombiTrans(wireR, 0,wireZ, wireRot);
  wireM1->RotateZ(45);
  wireM2->RotateZ(135);
  support->AddNode(wireVolume, 1, wireM1);
  support->AddNode(wireVolume, 2, wireM2);


  TGeoTorus*       wireTS  = new TGeoTorus(sheeveWheelZ+wireT, 0, wireT, 0, 
					   90-wireStartA*180/TMath::Pi());
  TGeoVolume*      wireTV  = new TGeoVolume("FMD3_bend_wire",wireTS,fSteel);
  TGeoRotation*    wireTR  = new TGeoRotation;
  wireTR->RotateY(90);
  wireTR->RotateZ(-90);
  Double_t         wireTZ  = sheeveWheelZ;
  TGeoCombiTrans*  wireTM1 = new TGeoCombiTrans(sheeveHoleR,0,wireTZ,wireTR);
  TGeoCombiTrans*  wireTM2 = new TGeoCombiTrans(sheeveHoleR,0,wireTZ,wireTR);
  wireTM1->RotateZ(45);
  wireTM2->RotateZ(135);
  support->AddNode(wireTV, 1, wireTM1);
  support->AddNode(wireTV, 2, wireTM2);

  Double_t         colarR = 4.05;
  Double_t         wireEL = sheeveHoleR - colarR;
  TGeoTube*        wireES = new TGeoTube("FMD3_end_wire", 0, wireT, wireEL/2);
  TGeoVolume*      wireEV = new TGeoVolume("FMD3_end_wire", wireES, fSteel);
  TGeoRotation*    wireER = new TGeoRotation;
  wireER->RotateY(90);
  TGeoCombiTrans*  wireEM1 = new TGeoCombiTrans(colarR+wireEL/2,0,
						-wireT,wireER);
  TGeoCombiTrans*  wireEM2 = new TGeoCombiTrans(colarR+wireEL/2,0,
						-wireT,wireER);
  wireEM1->RotateZ(45);
  wireEM2->RotateZ(135);
  support->AddNode(wireEV, 1, wireEM1);
  support->AddNode(wireEV, 2, wireEM2);
  
  

  
  //__________________________________________________________________
  // Place support volumes in half-detector volumes 
  Double_t         z  = fmd3->GetInnerZ();
  AliFMDDebug(1, ("FMD3 support at z=%f", -fmd3->GetNoseZ()));
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
