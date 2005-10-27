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

//____________________________________________________________________
//                                                                          
// Forward Multiplicity Detector based on Silicon wafers. This class
// contains the base procedures for the Forward Multiplicity detector
// Detector consists of 3 sub-detectors FMD1, FMD2, and FMD3, each of
// which has 1 or 2 rings of silicon sensors. 
//                                                       
// This is the base class for all FMD manager classes. 
//                    
// The actual code is done by various separate classes.   Below is
// diagram showing the relationship between the various FMD classes
// that handles the simulation
//
//      +--------+ 1     +-----------------+ 
//      | AliFMD |<>-----| AliFMDSimulator |
//      +--------+	 +-----------------+
//                               ^              
//                               |
//                 +-------------+-------------+
//                 |                           |	      
//        +--------------------+   +-------------------+
//        | AliFMDGeoSimulator |   | AliFMDG3Simulator | 
//        +--------------------+   +---------+---------+
//                                           ^
//                                           |
//                                +--------------------+
//				  | AliFMDOldSimulator |
//				  +--------------------+
//      
// *  AliFMD 
//    This defines the interface for the various parts of AliROOT that
//    uses the FMD, like AliFMDSimulator, AliFMDDigitizer, 
//    AliFMDReconstructor, and so on. 
//
// *  AliFMDSimulator
//    This is the base class for the FMD simulation tasks.   The
//    simulator tasks are responsible to implment the geoemtry, and
//    process hits. 
//                                                                          
// *  AliFMDGeoSimulator
//    This is a concrete implementation of the AliFMDSimulator that
//    uses the TGeo classes directly only.  This defines the active
//    volume as an ONLY XTRU shape with a divided MANY TUBS shape
//    inside to implement the particular shape of the silicon
//    sensors. 
//
// *  AliFMDG3Simulator
//    This is a concrete implementation of the AliFMDSimulator that
//    uses the TVirtualMC interface with GEANT 3.21-like messages.
//    This implements the active volume as a divided TUBS shape.  Hits
//    in the corners should be cut away at run time (but currently
//    isn't). 
//
// *  AliFMDOldSimulator
//    This is a concrete implementation of AliFMDSimulator.   It
//    approximates the of the rings as segmented disks. 
// 
#include "AliFMDGeoSimulator.h"	// ALIFMDGEOSIMULATOR_H
#include "AliFMDGeometry.h"	// ALIFMDGEOMETRY_H
#include "AliFMDDetector.h"	// ALIFMDDETECTOR_H
#include "AliFMDRing.h"		// ALIFMDRING_H
#include "AliFMD1.h"		// ALIFMD1_H
#include "AliFMD2.h"		// ALIFMD2_H
#include "AliFMD3.h"		// ALIFMD3_H
#include "AliFMD.h"		// ALIFMD_H
#include "AliLog.h"		// ALILOG_H
#include <TGeoVolume.h>		// ROOT_TGeoVolume
#include <TGeoTube.h>		// ROOT_TGeoTube
#include <TGeoPcon.h>		// ROOT_TGeoPcon
#include <TGeoMaterial.h>	// ROOT_TGeoMaterial
#include <TGeoMedium.h>		// ROOT_TGeoMedium
#include <TGeoXtru.h>		// ROOT_TGeoXtru
#include <TGeoPolygon.h>	// ROOT_TGeoPolygon
#include <TGeoTube.h>		// ROOT_TGeoTube
#include <TGeoManager.h>	// ROOT_TGeoManager
#include <TVector2.h>		// ROOT_TVector2
#include <TArrayD.h>		// ROOT_TArrayD

//====================================================================
ClassImp(AliFMDGeoSimulator)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDGeoSimulator::AliFMDGeoSimulator() 
  : fSi(0),
    fC(0),
    fAl(0),
    fPCB(0),
    fChip(0),
    fPlastic(0)
{
  // Default constructor
  fSectorOff   = 1;
  fModuleOff   = 4;
  fRingOff     = 5;
  fDetectorOff = 6;
}

//____________________________________________________________________
AliFMDGeoSimulator::AliFMDGeoSimulator(AliFMD* fmd, Bool_t detailed) 
  : AliFMDSimulator(fmd, detailed),
    fSi(0),
    fC(0),
    fAl(0),
    fPCB(0),
    fChip(0),
    fPlastic(0)
{
  // Normal constructor
  // 
  // Parameters: 
  // 
  //      fmd		Pointer to AliFMD object 
  //      detailed      Whether to make a detailed simulation or not 
  // 
  fSectorOff   = 1;
  fModuleOff   = 4;
  fRingOff     = 5;
  fDetectorOff = 6;
}

//____________________________________________________________________
void
AliFMDGeoSimulator::DefineMaterials() 
{
  // Define the materials and tracking mediums needed by the FMD
  // simulation.   These mediums are made by sending the messages
  // AliMaterial, AliMixture, and AliMedium to the passed AliModule
  // object module.   The defined mediums are 
  // 
  //	FMD Si$		Silicon (active medium in sensors)
  //	FMD C$		Carbon fibre (support cone for FMD3 and vacuum pipe)
  //	FMD Al$		Aluminium (honeycomb support plates)
  //	FMD PCB$	Printed Circuit Board (FEE board with VA1_ALICE)
  //	FMD Chip$	Electronics chips (currently not used)
  //	FMD Air$	Air (Air in the FMD)
  //	FMD Plastic$	Plastic (Support legs for the hybrid cards)
  //
  // Pointers to TGeoMedium objects are retrived from the TGeoManager
  // singleton.  These pointers are later used when setting up the
  // geometry 
  AliDebug(10, "\tCreating materials");

  if (!gGeoManager) {
    AliFatal("No TGeoManager defined");
    return;
  }
  AliFMDSimulator::DefineMaterials();
  fSi      = gGeoManager->GetMedium("FMD_Si$");
  fC       = gGeoManager->GetMedium("FMD_Carbon$");
  fAl      = gGeoManager->GetMedium("FMD_Aluminum$");
  fChip    = gGeoManager->GetMedium("FMD_Si Chip$");
  fAir     = gGeoManager->GetMedium("FMD_Air$");
  fPCB     = gGeoManager->GetMedium("FMD_PCB$");
  fPlastic = gGeoManager->GetMedium("FMD_Plastic$");
  fCopper  = gGeoManager->GetMedium("FMD_Copper$");
}

//____________________________________________________________________
TGeoVolume*
AliFMDGeoSimulator::RingGeometry(AliFMDRing* r) 
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
  const Int_t nv       = r->GetNVerticies();
  TVector2*   a        = r->GetVertex(5);
  TVector2*   b        = r->GetVertex(3);
  TVector2*   c        = r->GetVertex(4);
  Double_t    theta    = r->GetTheta();
  Double_t    off      = (TMath::Tan(TMath::Pi() * theta / 180) 
			  * r->GetBondingWidth());
  Double_t    rmax     = b->Mod();
  Double_t    rmin     = r->GetLowR();
  Double_t    pcbThick = r->GetPrintboardThickness();
  Double_t    modSpace = r->GetModuleSpacing();
  Double_t    legr     = r->GetLegRadius();
  Double_t    legl     = r->GetLegLength();
  Double_t    legoff   = r->GetLegOffset();
  Int_t       ns       = r->GetNStrips();
  Double_t    stripoff = a->Mod();
  Double_t    dstrip   = (rmax - stripoff) / ns;
  Double_t    space    = r->GetSpacing();
  TArrayD xs(nv);
  TArrayD ys(nv);
  for (Int_t i = 0; i < nv; i++) {
    // Reverse the order 
    TVector2* vv = r->GetVertex(nv - 1 - i);
    if (!vv) {
      AliError(Form("Failed to get vertex # %d", nv - 1 - i));
      continue;
    }
    xs[i] = vv->X();
    ys[i] = vv->Y();
  }
  
  // Shape of actual sensor 
  TGeoXtru* moduleShape = new TGeoXtru(2);
  moduleShape->DefinePolygon(nv, xs.fArray, ys.fArray);
  moduleShape->DefineSection(0, - siThick/2);
  moduleShape->DefineSection(1, siThick/2);
  TGeoVolume* moduleVolume = new TGeoVolume(Form(fgkModuleName, id), 
					    moduleShape, fSi);
  Int_t sid = moduleVolume->GetNumber();
  fSectorOff   = -1;
  fModuleOff   = 1;
  fRingOff     = 2;
  fDetectorOff = 3;
  if (fUseDivided) {
    fSectorOff   = 1;
    fModuleOff   = 4;
    fRingOff     = 5;
    fDetectorOff = 6;
    // Virtual volume shape to divide - This volume is only defined if
    // the geometry is set to be detailed. 
    TGeoVolume* activeVolume = 0;
    if (fDetailed) {
      TGeoTubeSeg* activeShape = 
	new TGeoTubeSeg(rmin, rmax, siThick/2, - theta, theta);
      activeVolume = new TGeoVolume(Form(fgkActiveName, id), activeShape, fSi);
      TGeoVolume* sectorVolume = activeVolume->Divide(Form(fgkSectorName, id), 
						      2, 2, -theta, 0, 0, "N");
      TGeoVolume* stripVolume = sectorVolume->Divide(Form(fgkStripName, id), 
						     1, ns, stripoff, dstrip, 
						     0, "SX");
      sid = stripVolume->GetNumber();
    }
    // Add divived MANY volume to the true shape of the module, but only
    // if a detailed simulation is reguested. 
    if (activeVolume) moduleVolume->AddNodeOverlap(activeVolume, 0);
  }
  
  switch (id) {
  case 'i':
  case 'I': fActiveId[0] = sid; break;
  case 'o':
  case 'O': fActiveId[2] = sid; break;
  }

  // Shape of Printed circuit Board 
  TGeoXtru* pcbShape = new TGeoXtru(2);
  for (Int_t i = 0;      i < nv / 2; i++) ys[i] -= off;
  for (Int_t i = nv / 2; i < nv;     i++) ys[i] += off;
  pcbShape->DefinePolygon(nv, xs.fArray, ys.fArray);
  pcbShape->DefineSection(0, - pcbThick/2);
  pcbShape->DefineSection(1, pcbThick/2);
  TGeoVolume* pcbVolume = new TGeoVolume(Form(fgkPCBName, id, 'B'), 
					 pcbShape, fPCB);

  // Short leg shape 
  TGeoTube*   shortLegShape  = new TGeoTube(0, legr, legl / 2);
  TGeoVolume* shortLegVolume = new TGeoVolume(Form(fgkShortLegName, id), 
					      shortLegShape, fPlastic);

  // Long leg shape
  TGeoTube*   longLegShape   = new TGeoTube(0, legr, (legl + modSpace) / 2);
  TGeoVolume* longLegVolume  = new TGeoVolume(Form(fgkLongLegName, id), 
					      longLegShape, fPlastic);
  
  TGeoMatrix* matrix = 0;
  // Back container volume 
  Double_t contThick     = siThick + pcbThick + legl;
  TGeoTubeSeg* backShape = new TGeoTubeSeg(rmin, rmax, contThick/2, 
					   - theta, theta);
  TGeoVolume* backVolume = new TGeoVolume(Form(fgkBackVName, id), 
					  backShape, fAir);
  Double_t x = 0;
  Double_t y = 0;
  Double_t z = -contThick / 2 + siThick / 2;
  matrix     = new TGeoTranslation(Form("FMD Ring  %c mod 1 transform", id), 
				   x, y, z);
  backVolume->AddNode(moduleVolume, 0, matrix);
  z          += siThick / 2 + space + pcbThick / 2;
  matrix     =  new TGeoTranslation(Form("FMD Ring %c pcb 1 transfrom", id), 
				    x, y, z);
  backVolume->AddNode(pcbVolume, 0, matrix);
  x          =  a->X() + legoff + legr;
  y          =  0;
  z          += pcbThick / 2 + legl / 2;
  matrix     = new TGeoTranslation(Form("FMD Ring %c leg 1 transfrom", id), 
				   x, y, z);
  backVolume->AddNode(shortLegVolume, 0, matrix);
  x          =  c->X();
  y          =  c->Y() - legoff - legr - off;
  matrix     =  new TGeoTranslation(Form("FMD Ring %c leg 2 transfrom", id), 
				    x, y, z);
  backVolume->AddNode(shortLegVolume, 1, matrix);
  y          =  -y;
  matrix     =  new TGeoTranslation(Form("FMD Ring %c leg 3 transfrom", id), 
				    x, y, z);
  backVolume->AddNode(shortLegVolume, 2, matrix);
  // backVolume->SetVisibility(kFALSE);
  // backVolume->VisibleDaughters(kTRUE);

  // Front container volume 
  contThick += modSpace;
  TGeoTubeSeg* frontShape = new TGeoTubeSeg(rmin, rmax, contThick/2, 
					    -theta, theta);
  TGeoVolume* frontVolume = new TGeoVolume(Form(fgkFrontVName, id),
					   frontShape, fAir);
  x         =  0;
  y         =  0;
  z         =  -contThick / 2 + siThick / 2 ;
  matrix    = new TGeoTranslation(Form("FMD Ring %c mod 2 transfrom", id), 
				  0, 0, z);
  frontVolume->AddNode(moduleVolume, 1, matrix);
  z         += siThick / 2 + space + pcbThick / 2;
  matrix    =  new TGeoTranslation(Form("FMD Ring %c pcb 2 transfrom", id), 
				   x, y, z);
  frontVolume->AddNode(pcbVolume, 1, matrix);
  x         =  a->X() + legoff + legr;
  y         =  0;
  z         += pcbThick / 2 + (legl + modSpace)/ 2;
  matrix    =  new TGeoTranslation(Form("FMD Ring %c leg 4 transfrom", id), 
				   x, y, z);
  frontVolume->AddNode(longLegVolume, 0, matrix);
  x         =  c->X();
  y         =  c->Y() - legoff - legr - off;
  matrix    =  new TGeoTranslation(Form("FMD Ring %c leg 4 transfrom", id), 
				   x, y, z);
  frontVolume->AddNode(longLegVolume, 1, matrix);
  y         =  -y;
  matrix    =  new TGeoTranslation(Form("FMD Ring %c leg 4 transfrom", id), 
				   x, y, z);
  frontVolume->AddNode(longLegVolume, 2, matrix);
  // frontVolume->SetVisibility(kFALSE);
  // frontVolume->VisibleDaughters(kTRUE);
  
  // Ring mother volume 
  TGeoTube* ringShape    = new TGeoTube(rmin, rmax, contThick / 2);
  TGeoVolume* ringVolume = new TGeoVolume(Form(fgkRingName, id), ringShape, 
					  fAir);

  Int_t nmod = r->GetNModules();
  AliDebug(10, Form("making %d modules in ring %c", nmod, id));
  for (Int_t i = 0; i < nmod; i++) {
    Bool_t isFront    = (i % 2 == 0);
    TGeoVolume* vol   = (isFront ? frontVolume : backVolume);
    TGeoRotation* rot = new TGeoRotation(Form("FMD Ring %c rotation %d",id,i));
    rot->RotateZ((i + .5) * 2 * theta);
    Double_t z = (isFront ? 0 : modSpace) / 2;
    matrix     = new TGeoCombiTrans(Form("FMD Ring %c transform %d", id, i), 
				    0, 0, z, rot);
    ringVolume->AddNode(vol, i, matrix);
  }

  ringVolume->SetVisibility(kFALSE);
  ringVolume->VisibleDaughters(kTRUE);
  return ringVolume;
}

//____________________________________________________________________
TGeoVolume*
AliFMDGeoSimulator::DetectorGeometry(AliFMDDetector* d, 
				     TGeoVolume* mother, 
				     Double_t    zmother, 
				     TGeoVolume* inner, 
				     TGeoVolume* outer) 
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
    TGeoVolume* rvol  = 0;
    switch (i) {
    case 0: 
      r      = d->GetInner();
      lowr   = d->GetInnerHoneyLowR();
      highr  = d->GetInnerHoneyHighR();
      rz     = d->GetInnerZ();
      rvol   = inner;
      break;
    case 1: 
      r      = d->GetOuter();
      lowr   = d->GetOuterHoneyLowR();
      highr  = d->GetOuterHoneyHighR();
      rz     = d->GetOuterZ();
      rvol   = outer;
      break;
    }
    if (!r) continue;
    Char_t   c       = r->GetId();
    Int_t    id      = d->GetId();
    Double_t hcThick = d->GetHoneycombThickness();
    Double_t alThick = d->GetAlThickness();
    Double_t z;
    if (zmother > 0) z = rz - zmother + r->GetRingDepth() / 2;
    else             z = zmother - rz + r->GetRingDepth() / 2;
    // Place ring in mother volume
    mother->AddNode(rvol, Int_t(c), 
		    new TGeoTranslation(Form("FMD%d%c transform", id, c), 
					0, 0, z));

    z += r->GetRingDepth() / 2 + hcThick / 2;
    // Top of Honeycomb
    TGeoTubeSeg* topHCShape  = new TGeoTubeSeg(lowr, highr, hcThick/2, 0, 180);
    TGeoVolume*  topHCVolume = new TGeoVolume(Form(fgkTopHCName, id, c), 
					      topHCShape, fAl);
    TGeoMatrix*  topHCMatrix = 
      new TGeoTranslation(Form("FMD%d%c top HC transform", id, c), 0, 0, z);
    mother->AddNode(topHCVolume, 0, topHCMatrix);

    // Air in top of honeycomb
    TGeoTubeSeg* topIHCShape = new TGeoTubeSeg(lowr+alThick, highr - alThick, 
					       (hcThick-alThick)/2, 0, 180);
    TGeoVolume*  topIHCVolume = new TGeoVolume(Form(fgkTopIHCName, id, c), 
					       topIHCShape, fAir);
    topHCVolume->AddNode(topIHCVolume, 0);
    topHCVolume->VisibleDaughters(kFALSE);    
    topHCVolume->SetVisibility(kTRUE);


    // Bottom of Honeycomb
    TGeoTubeSeg* botHCShape = new TGeoTubeSeg(lowr, highr, hcThick/2, 
					      180, 360);
    TGeoVolume*  botHCVolume = new TGeoVolume(Form(fgkBotHCName, id, c), 
					      botHCShape, fAl);
    TGeoMatrix*  botHCMatrix = 
      new TGeoTranslation(Form("FMD%d%c bottom HC transform", id, c), 0, 0, z);
    mother->AddNode(botHCVolume, 0, botHCMatrix);

    // Air in bot of honeycomb
    TGeoTubeSeg* botIHCShape = new TGeoTubeSeg(lowr+alThick, highr - alThick, 
					       (hcThick-alThick)/2, 180, 360);
    TGeoVolume*  botIHCVolume = new TGeoVolume(Form(fgkBotIHCName, id, c), 
					       botIHCShape, fAir);
    botHCVolume->AddNode(botIHCVolume, 0);
    botHCVolume->VisibleDaughters(kFALSE);    
    botHCVolume->SetVisibility(kTRUE);    
  }
  mother->SetVisibility(kFALSE);
  mother->VisibleDaughters(kTRUE);
  return mother;
}

//____________________________________________________________________
TGeoVolume*
AliFMDGeoSimulator::FMD1Geometry(AliFMD1* fmd1, TGeoVolume* inner) 
{
  // Setup the FMD1 geometry.  The FMD1 only has one ring, and no
  // special support as it is at the momement. 
  // 
  // See also AliFMDGeoSimulator::DetectorGeometry 
  // 
  if (!fmd1 || !inner) return 0;
  Double_t rmin    = fmd1->GetInner()->GetLowR();
  Double_t rmax    = fmd1->GetInnerHoneyHighR();
  Double_t hcThick = fmd1->GetHoneycombThickness();
  Double_t w       = fmd1->GetInner()->GetRingDepth() + hcThick;
  Double_t z       = fmd1->GetInnerZ() + w / 2;

  TGeoVolume* fmd1Volume = 0;
  if (!fUseAssembly) {
    TGeoTube* fmd1Shape = new TGeoTube(rmin, rmax, w / 2);
    fmd1Volume = new TGeoVolume(fmd1->GetName(), fmd1Shape, fAir);
  }
  else
    fmd1Volume = new TGeoVolumeAssembly(fmd1->GetName());
  
  TGeoVolume* top = gGeoManager->GetVolume("ALIC");
  TGeoMatrix* matrix = new TGeoTranslation("FMD1 transform", 0, 0, z);
  top->AddNode(fmd1Volume, fmd1->GetId(), matrix);

  return DetectorGeometry(fmd1, fmd1Volume, z, inner, 0);
}

//____________________________________________________________________
TGeoVolume*
AliFMDGeoSimulator::FMD2Geometry(AliFMD2* fmd2, 
				 TGeoVolume* inner, 
				 TGeoVolume* outer) 
{
  // Setup the FMD2 geometry.  The FMD2 has no
  // special support as it is at the momement. 
  // 
  // See also AliFMDGeoSimulator::DetectorGeometry 
  // 
  if (!fmd2 || !inner || !outer) return 0;
  Double_t rmin     = fmd2->GetInner()->GetLowR();
  Double_t rmax     = fmd2->GetOuterHoneyHighR();
  Double_t hcThick  = fmd2->GetHoneycombThickness();
  Double_t ow       = fmd2->GetInner()->GetRingDepth();
  Double_t iz       = fmd2->GetInnerZ();
  Double_t oz       = fmd2->GetOuterZ();
  Double_t w        = TMath::Abs(oz - iz) + ow + hcThick;
  Double_t z        = oz + w / 2;
  
  TGeoVolume* fmd2Volume = 0;
  if (!fUseAssembly) {
    TGeoTube* fmd2Shape = new TGeoTube(rmin, rmax, w / 2);
    fmd2Volume = new TGeoVolume(fmd2->GetName(), fmd2Shape, fAir);
  }
  else 
    fmd2Volume = new TGeoVolumeAssembly(fmd2->GetName());
  
  TGeoVolume* top = gGeoManager->GetVolume("ALIC");
  TGeoMatrix* matrix = new TGeoTranslation("FMD2 transform", 0, 0, z);
  top->AddNode(fmd2Volume, fmd2->GetId(), matrix);

  return DetectorGeometry(fmd2, fmd2Volume, z, inner, outer);
}
  
//____________________________________________________________________
TGeoVolume*
AliFMDGeoSimulator::FMD3Geometry(AliFMD3* fmd3, 
				 TGeoVolume* inner, 
				 TGeoVolume* outer) 
{
  // Setup the FMD3 geometry.  The FMD2 has a rather elaborate support
  // structure, as the support will also support the vacuum
  // beam-pipe. 
  // 
  // See also AliFMDGeoSimulator::DetectorGeometry 
  // 
  if (!fmd3 || !inner || !outer) return 0;
  Double_t nlen    = fmd3->GetNoseLength();
  Double_t nz      = fmd3->GetNoseZ();
  Double_t noser1  = fmd3->GetNoseLowR();
  Double_t noser2  = fmd3->GetNoseHighR();
  Double_t conel   = fmd3->GetConeLength();
  Double_t backl   = fmd3->GetBackLength();
  Double_t backr1  = fmd3->GetBackLowR();
  Double_t backr2  = fmd3->GetBackHighR();
  Double_t zdist   = conel -  backl - nlen;
  Double_t tdist   = backr2 - noser2;
  Double_t beaml   = TMath::Sqrt(zdist * zdist + tdist * tdist);
  Double_t theta   = -180. * TMath::ATan2(tdist, zdist) / TMath::Pi();
  Double_t innerZ  = fmd3->GetInnerZ();
  Double_t innerZh = (innerZ - fmd3->GetInner()->GetRingDepth() 
		      - fmd3->GetHoneycombThickness());
  Double_t outerZ  = fmd3->GetOuterZ();
  Double_t outerZh = (outerZ - fmd3->GetOuter()->GetRingDepth() 
		      - fmd3->GetHoneycombThickness());
  Double_t innerr1 = fmd3->GetInner()->GetLowR();
  // Double_t innerr2 = fmd3->GetInner()->GetHighR();
  Double_t outerr1 = fmd3->GetOuter()->GetLowR();
  // Double_t outerr2 = fmd3->GetOuter()->GetHighR();
  Double_t flanger = fmd3->GetFlangeR();
  Double_t minZ    = TMath::Min(nz - conel, outerZh);
  Double_t z       = fmd3->GetZ();
  Double_t zi;

  // FMD3 volume 
  TGeoVolume* fmd3Volume = 0;
  if (!fUseAssembly) {
    TGeoPcon* fmd3Shape = new TGeoPcon(0, 360, 8);
    zi = z - nz;
    fmd3Shape->DefineSection(0, zi,  noser1,   noser2);
    zi = z - (nz - nlen);
    fmd3Shape->DefineSection(1, zi,  noser1,   fmd3->ConeR(z - zi)+.15);
    zi = z - innerZ;
    fmd3Shape->DefineSection(2, zi,  innerr1,  fmd3->ConeR(z - zi)+.15);
    zi = z - innerZh;
    fmd3Shape->DefineSection(3, zi,  innerr1,  fmd3->ConeR(z - zi)+.15);
    fmd3Shape->DefineSection(4, zi,  outerr1,  fmd3->ConeR(z - zi)+.15);
    zi = z - nz + zdist + nlen;
    fmd3Shape->DefineSection(5, zi,  outerr1,  fmd3->ConeR(z - zi)+.15);
    zi = z - nz + nlen + zdist;
    fmd3Shape->DefineSection(6, zi,  outerr1,  flanger+1.5);
    zi = z - minZ;
    fmd3Shape->DefineSection(7, zi,  outerr1,  flanger+1.5);
    fmd3Volume = new TGeoVolume(fmd3->GetName(), fmd3Shape, fAir);
  }
  else 
    fmd3Volume = new TGeoVolumeAssembly(fmd3->GetName());
  
  TGeoRotation* rot = new TGeoRotation("FMD3 rotatation");
  rot->RotateY(180);
  TGeoVolume* top = gGeoManager->GetVolume("ALIC");
  TGeoMatrix* mmatrix = new TGeoCombiTrans("FMD3 transform", 0, 0, z, rot);
  top->AddNode(fmd3Volume, fmd3->GetId(), mmatrix);
  
  // Nose volume 
  TGeoTube* noseShape = new TGeoTube(noser1, noser2, nlen / 2);
  TGeoVolume* noseVolume = new TGeoVolume(fgkNoseName, noseShape, fC);
  zi = z - nz + nlen / 2;
  TGeoMatrix* nmatrix = new TGeoTranslation("FMD3 Nose translation", 0, 0, zi);
  // fmd3Volume->AddNodeOverlap(noseVolume, 0, nmatrix);
  fmd3Volume->AddNode(noseVolume, 0, nmatrix);
  
  // Back
  TGeoTube* backShape = new TGeoTube(backr1, backr2, backl / 2);
  TGeoVolume* backVolume = new TGeoVolume(fgkBackName, backShape, fC);
  zi = z - nz + conel - backl / 2;
  TGeoMatrix* bmatrix = new TGeoTranslation("FMD3 Back translation", 0, 0, zi);
  fmd3Volume->AddNode(backVolume, 0, bmatrix);
  
  Int_t n;
  Double_t r;
  // The flanges 
  TGeoBBox* flangeShape = new TGeoBBox((flanger - backr2) / 2, 
				       fmd3->GetBeamWidth() / 2,
				       backl / 2);
  TGeoVolume* flangeVolume = new TGeoVolume(fgkFlangeName, flangeShape, fC);
  n = fmd3->GetNFlange();
  r = backr2 + (flanger - backr2) / 2;
  for (Int_t i = 0; i  < n; i++) {
    Double_t phi = 360. / n * i + 180. / n;
    Double_t x   = r * TMath::Cos(TMath::Pi() / 180 * phi);
    Double_t y   = r * TMath::Sin(TMath::Pi() / 180 * phi);
    TGeoRotation* rot = new TGeoRotation(Form("FMD3 Flange rotation %d", i));
    rot->RotateZ(phi);
    TGeoMatrix* matrix = new TGeoCombiTrans(Form("FMD3 flange transform %d", 
						 i), x, y, zi, rot);
    // fmd3Volume->AddNodeOverlap(flangeVolume, i, matrix);
    fmd3Volume->AddNode(flangeVolume, i, matrix);
    
  }

  // The Beams 
  TGeoBBox* beamShape = new TGeoBBox(fmd3->GetBeamThickness() / 2, 
				     fmd3->GetBeamWidth() / 2 - .1,
				     beaml / 2);
  TGeoVolume* beamVolume = new TGeoVolume(fgkBeamName, beamShape, fC);
  n = fmd3->GetNBeam();
  r = noser2 + tdist / 2;
  zi = z - nz + nlen + zdist / 2;
  for (Int_t i = 0; i  < n; i++) {
    Double_t phi = 360. / n * i;
    Double_t x   = r * TMath::Cos(TMath::Pi() / 180 * phi);
    Double_t y   = r * TMath::Sin(TMath::Pi() / 180 * phi);
    TGeoRotation* rot = new TGeoRotation(Form("FMD3 beam rotation %d", i));
    // Order is important
    rot->RotateY(-theta);
    rot->RotateZ(phi);
    TGeoMatrix* matrix = new TGeoCombiTrans(Form("FMD3 beam transform %d", i),
					    x, y, zi, rot);
    fmd3Volume->AddNode(beamVolume, i, matrix);    
  }
  
  
  return DetectorGeometry(fmd3, fmd3Volume, z, inner, outer);
}

//____________________________________________________________________
void
AliFMDGeoSimulator::DefineGeometry() 
{
  // Setup up the FMD geometry. 
  AliDebug(10, "Setting up volume");

  AliFMDGeometry* fmd = AliFMDGeometry::Instance();
  TGeoVolume* inner = RingGeometry(fmd->GetInner());
  TGeoVolume* outer = RingGeometry(fmd->GetOuter());

  if (!inner || !outer) {
    AliError("Failed to create one of the ring volumes");
    return;
  }
  FMD1Geometry(fmd->GetFMD1(), inner);
  FMD2Geometry(fmd->GetFMD2(), inner, outer);
  FMD3Geometry(fmd->GetFMD3(), inner, outer);
}


//____________________________________________________________________
//
// EOF
//
