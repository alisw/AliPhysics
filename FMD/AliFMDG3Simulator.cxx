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
//        +--------------------+   +-------------------+
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
#include <math.h>
#include "AliFMDG3Simulator.h"	// ALIFMDG3SIMULATOR_H
#include "AliFMDGeometry.h"	// ALIFMDGEOMETRY_H
#include "AliFMDDetector.h"	// ALIFMDDETECTOR_H
#include "AliFMDRing.h"		// ALIFMDRING_H
#include "AliFMD1.h"		// ALIFMD1_H
#include "AliFMD2.h"		// ALIFMD2_H
#include "AliFMD3.h"		// ALIFMD3_H
#include "AliFMD.h"		// ALIFMD_H
#include <AliLog.h>		// ALILOG_H
#include <TVector2.h>		// ROOT_TVector2
#include <TVirtualMC.h>		// ROOT_TVirtualMC
#include <TArrayI.h>		// ROOT_TArrayI

//====================================================================
ClassImp(AliFMDG3Simulator)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDG3Simulator::AliFMDG3Simulator() 
{
  // Default constructor
  fSectorOff   = 1;
  fModuleOff   = 3;
  fRingOff     = 4;
  fDetectorOff = 5;
}

//____________________________________________________________________
AliFMDG3Simulator::AliFMDG3Simulator(AliFMD* fmd, Bool_t detailed) 
  : AliFMDSimulator(fmd, detailed)
{
  // Normal constructor
  // 
  // Parameters: 
  // 
  //      fmd		Pointer to AliFMD object 
  //      detailed      Whether to make a detailed simulation or not 
  // 
  fSectorOff   = 1;
  fModuleOff   = 3;
  fRingOff     = 4;
  fDetectorOff = 5;
}

//____________________________________________________________________
Bool_t
AliFMDG3Simulator::RingGeometry(AliFMDRing* r) 
{
  // Setup the geometry of a ring.    The defined TGeoVolume is
  // returned, and should be used when setting up the rest of the
  // volumes. 
  // 
  // Parameters:
  //
  //     r		Pointer to ring geometry object 
  // 
  // Returns:
  //    true on success 
  //
  if (!r) { 
    AliError("Didn't get a ring object");
    return kFALSE;
  }
  Char_t      id       = r->GetId();
  Double_t    siThick  = r->GetSiThickness();
  // const Int_t nv       = r->GetNVerticies();
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
  Double_t    space    = r->GetSpacing();
  Double_t    stripoff = a->Mod();
  Double_t    dstrip   = (rmax - stripoff) / ns;
  Double_t    par[10];
  TString     name;
  TString     name2;
  TVirtualMC* mc       = TVirtualMC::GetMC();
  
  Int_t siId  = fFMD->GetIdtmed()->At(kSiId);
  Int_t airId = fFMD->GetIdtmed()->At(kAirId);
  Int_t pcbId = fFMD->GetIdtmed()->At(kPcbId);
  Int_t plaId = fFMD->GetIdtmed()->At(kPlasticId);
  Int_t copId = fFMD->GetIdtmed()->At(kCopperId);
  Int_t chiId = fFMD->GetIdtmed()->At(kSiChipId);

  Double_t ringWidth  = r->GetRingDepth();
  Double_t x          = 0;
  Double_t y          = 0;
  Double_t z          = 0;
  Double_t backWidth  = siThick + pcbThick + legl + space;
  Double_t frontWidth = backWidth + modSpace;

  // Ring mother volume 
  par[0]     =  rmin;
  par[1]     =  rmax;
  par[2]     =  ringWidth / 2;
  name       =  Form(fgkRingName, id);
  mc->Gsvolu(name.Data(), "TUBE", airId, par, 3);

  // Back container volume 
  par[0] = rmin;
  par[1] = rmax;
  par[2] = backWidth / 2;
  par[3] = -theta;
  par[4] = +theta;
  TString backName(Form(fgkBackVName, id));
  mc->Gsvolu(backName.Data(), "TUBS", airId, par, 5);
  
  // Front container volume 
  par[2]     =  frontWidth / 2;
  TString frontName(Form(fgkFrontVName, id));
  mc->Gsvolu(frontName.Data(), "TUBS", airId, par, 5);
  
  Double_t topL = (b->X() - c->X());
  Double_t botL = (c->X() - a->X());
  Int_t    rot;
  mc->Matrix(rot, 90, 90, 0, 90, 90, 0);  

  Double_t zFront = - frontWidth / 2 + siThick / 2;
  Double_t zBack  = - backWidth / 2 + siThick / 2;
  if (fUseDivided) {
    fSectorOff   = 1;
    fModuleOff   = 3;
    fRingOff     = 4;
    fDetectorOff = 5;

    // Virtual volume shape to divide - This volume is only defined if
    // the geometry is set to be detailed. 
    par[0]     = rmin;
    par[1]     = rmax;
    par[2]     = siThick / 2;
    par[3]     = -theta;
    par[4]     = theta;
    name       = Form(fgkActiveName, id);
    mc->Gsvolu(name.Data(), "TUBS", (fDetailed ? airId : siId), par, 5);

    mc->Gspos(name.Data(), 0, backName.Data(), x, y, zBack, 0, "ONLY");
    mc->Gspos(name.Data(), 0, frontName.Data(), x, y, zFront, 0, "ONLY");
    
    Int_t sid = -1;
    if (fDetailed) {
      // Divide the volume into sectors
      name2 = name;
      name  = Form(fgkSectorName, id);
      mc->Gsdvn2(name.Data(), name2.Data(), 2, 2, -theta, siId);
      
      // Divide the volume into strips
      name2 = name;
      name  = Form(fgkStripName, id);
      mc->Gsdvt2(name.Data(), name2.Data(), dstrip, 1, stripoff, siId, ns);
      sid = mc->VolId(name.Data());
      AliDebug(10, Form("Got volume id %d for volume %s", sid, name.Data()));
    }
  
    switch (id) {
    case 'i': case 'I': fActiveId[0] = sid; break;
    case 'o': case 'O': fActiveId[2] = sid; break;
    }
  }
  else {
    fSectorOff   = -1;
    fModuleOff   = 1;
    fRingOff     = 2;
    fDetectorOff = 3;

    // Create top of module shape 
    par[0]    = c->Y();
    par[1]    = b->Y();
    par[2]    = siThick / 2;
    par[3]    = topL / 2;
    name      = Form(fgkModuleName, id);
    name[3]   = 'T';
    mc->Gsvolu(name.Data(), "TRD1", siId, par, 4);
    Int_t tid = mc->VolId(name.Data());
    x         =  rmin + botL + topL / 2;
    mc->Gspos(name.Data(), 0, backName.Data(), x, y, zBack, rot, "ONLY");
    mc->Gspos(name.Data(), 0, frontName.Data(), x, y, zFront, rot, "ONLY");


    // Create bottom of module shape 
    par[0]    = a->Y();
    par[1]    = c->Y();
    par[3]    = botL / 2;
    name      = Form(fgkModuleName, id);
    name[3]   = 'B';
    mc->Gsvolu(name.Data(), "TRD1", siId, par, 4);
    Int_t bid = mc->VolId(name.Data());
    x         =  rmin + botL / 2;
    z         =  - backWidth / 2 + siThick / 2;
    mc->Gspos(name.Data(), 0, backName.Data(), x, y, zBack, rot, "ONLY");
    mc->Gspos(name.Data(), 0, frontName.Data(), x, y, zFront, rot, "ONLY");
  
    switch (id) {
    case 'i': case 'I': fActiveId[0] = tid; fActiveId[1] = bid; break;
    case 'o': case 'O': fActiveId[2] = tid; fActiveId[3] = bid; break;
    }
  }
  
    
  // Shape of Printed circuit Board 
  // Top
  par[0] =  c->Y() - off;
  par[1] =  b->Y() - off;
  par[2] =  pcbThick / 2;
  par[3] =  topL / 2;
  x      =  rmin + botL + topL / 2;
  zBack  += siThick / 2 + space + pcbThick / 2;
  zFront += siThick / 2 + space + pcbThick / 2;
  name   = Form(fgkPCBName, id, 'T');
  mc->Gsvolu(name.Data(), "TRD1", pcbId, par, 4);
  mc->Gspos(name.Data(), 0, backName.Data(), x, y, zBack, rot, "ONLY");
  mc->Gspos(name.Data(), 0, frontName.Data(), x, y, zFront, rot, "ONLY");
  
  // Bottom
  par[0] = a->Y() - off;
  par[1] = c->Y() - off;
  par[3] = botL / 2;
  name   = Form(fgkPCBName, id, 'B');
  x      =  rmin + botL / 2;
  mc->Gsvolu(name.Data(), "TRD1", pcbId, par, 4);
  mc->Gspos(name.Data(), 0, backName.Data(), x, y, zBack, rot, "ONLY");
  mc->Gspos(name.Data(), 0, frontName.Data(), x, y, zFront, rot, "ONLY");

  Double_t x1, y1;
  // Short leg volume 
  par[0] =  legr - .1;
  par[1] =  legr;
  par[2] =  legl / 2;
  x      =  a->X() + legoff + legr;
  x1     =  c->X();
  y1     =  c->Y() - legoff - legr - off;
  zBack  += pcbThick / 2 + legl / 2;
  zFront += pcbThick / 2 + legl / 2 + modSpace / 2;
  name   = Form(fgkShortLegName, id);
  mc->Gsvolu(name.Data(),  "TUBE",  plaId, par, 3);
  mc->Gspos(name.Data(), 0, backName.Data(), x,    y, zBack, 0, "ONLY");
  mc->Gspos(name.Data(), 1, backName.Data(), x1,  y1, zBack, 0, "ONLY");
  mc->Gspos(name.Data(), 2, backName.Data(), x1, -y1, zBack, 0, "ONLY");
  
  // Long leg volume 
  par[2] += modSpace / 2;
  name   = Form(fgkLongLegName, id);
  mc->Gsvolu(name.Data(),  "TUBE",  plaId, par, 3);
  mc->Gspos(name.Data(), 0, frontName.Data(), x,    y, zFront, 0, "ONLY");
  mc->Gspos(name.Data(), 1, frontName.Data(), x1,  y1, zFront, 0, "ONLY");
  mc->Gspos(name.Data(), 2, frontName.Data(), x1, -y1, zFront, 0, "ONLY");
  
  // Place modules+pcb+legs in ring volume 
  Int_t nmod = r->GetNModules();
  name2      =  Form(fgkRingName, id);
  AliDebug(10, Form("making %d modules in ring %c", nmod, id));
  for (Int_t i = 0; i < nmod; i++) {
    Double_t th      = (i + .5) * 2 * theta;
    Bool_t   isFront = (i % 2 == 0);
    name             = (isFront ? frontName : backName);
    z                = (isFront ? 0 : modSpace) / 2;
    mc->Matrix(rot, 90, th, 90, fmod(90 + th, 360), 0, 0);
    mc->Gspos(name.Data(), i, name2.Data(), 0, 0, z, rot, "ONLY");
  }
  
  return kTRUE;
}

//____________________________________________________________________
Bool_t
AliFMDG3Simulator::DetectorGeometry(AliFMDDetector* d, Double_t zmother) 
{
  // Common stuff for setting up the FMD1, FMD2, and FMD3 geometries.
  // This includes putting the Honeycomb support plates and the rings
  // into the mother volumes.   
  // 
  // Parameeters:
  //	d	  The detector geometry to use 
  //    zmother	  The midpoint in global coordinates of detector vol.
  //
  // Returns:
  //    true on success
  // 
  if (!d) return kFALSE;

  TString     name;
  TString     name2;
  TVirtualMC* mc       = TVirtualMC::GetMC();

  // Loop over the defined rings 
  for (int i = 0; i < 2; i++) {
    AliFMDRing* r     = 0;
    Double_t    lowr  = 0;
    Double_t    highr = 0;
    Double_t    rz    = 0;
    switch (i) {
    case 0: 
      r      = d->GetInner();
      lowr   = d->GetInnerHoneyLowR();
      highr  = d->GetInnerHoneyHighR();
      rz     = d->GetInnerZ();
      break;
    case 1: 
      r      = d->GetOuter();
      lowr   = d->GetOuterHoneyLowR();
      highr  = d->GetOuterHoneyHighR();
      rz     = d->GetOuterZ();
      break;
    }
    if (!r) continue;
    Char_t   c       = r->GetId();
    Int_t    id      = d->GetId();
    Int_t    airId   = (fFMD->GetIdtmed()->At(kAirId));
    Int_t    alId    = (fFMD->GetIdtmed()->At(kAlId));
    Double_t hcThick = d->GetHoneycombThickness();
    Double_t alThick = d->GetAlThickness();
    Double_t par[10];
    Double_t z;
    // Place ring in mother volume
    if (zmother > 0) z = rz - zmother + r->GetRingDepth() / 2;
    else             z = zmother - rz + r->GetRingDepth() / 2;
    name  = Form(fgkRingName, c);
    name2 = d->GetName();
    mc->Gspos(name.Data(), Int_t(c), name2.Data(), 0, 0, z, 0, "ONLY");
    
    // Place Top Honeycomb in mother volume 
    z += + r->GetRingDepth() / 2 + hcThick / 2;
    // Top of Honeycomb
    par[0] = lowr;
    par[1] = highr;
    par[2] = hcThick / 2;
    par[3] = 0;
    par[4] = 180;
    name   = Form(fgkTopHCName, id, c);
    mc->Gsvolu(name.Data(), "TUBS", alId, par, 5);
    mc->Gspos(name.Data(), 0, name2.Data(), 0, 0, z, 0, "ONLY");

    par[0] += alThick;
    par[1] -= alThick;
    par[2] -= alThick / 2;
    name2  =  name;
    name   =  Form(fgkTopIHCName, id, c);
    mc->Gsvolu(name.Data(), "TUBS", airId, par, 5);
    mc->Gspos(name.Data(), 0, name2.Data(), 0, 0, 0, 0, "ONLY");
    
    // Bot of Honeycomb
    par[0] = lowr;
    par[1] = highr;
    par[2] = hcThick / 2;
    par[3] = 180;
    par[4] = 360;
    name2  = d->GetName();
    name   = Form(fgkBotHCName, id, c);
    mc->Gsvolu(name.Data(), "TUBS", alId, par, 5);
    mc->Gspos(name.Data(), 0, name2.Data(), 0, 0, z, 0, "ONLY");

    par[0] += alThick;
    par[1] -= alThick;
    par[2] -= alThick / 2;
    name2  =  name;
    name   =  Form(fgkBotIHCName, id, c);
    mc->Gsvolu(name.Data(), "TUBS", airId, par, 5);
    mc->Gspos(name.Data(), 0, name2.Data(), 0, 0, 0, 0, "ONLY");
  }
  return kTRUE;
}

//____________________________________________________________________
Bool_t
AliFMDG3Simulator::FMD1Geometry(AliFMD1* fmd1) 
{
  // Setup the FMD1 geometry.  The FMD1 only has one ring, and no
  // special support as it is at the momement. 
  // 
  // See also AliFMDG3Simulator::DetectorGeometry 
  // 
  if (!fmd1) return kFALSE;
  Double_t rmin    = fmd1->GetInner()->GetLowR();
  Double_t rmax    = fmd1->GetInnerHoneyHighR();
  Double_t hcThick = fmd1->GetHoneycombThickness();
  Double_t w       = fmd1->GetInner()->GetRingDepth() + hcThick;
  Double_t z       = fmd1->GetInnerZ() + w / 2;
  TVirtualMC* mc   = TVirtualMC::GetMC();
  Int_t    airId   = (fFMD->GetIdtmed()->At(kAirId));

  Double_t par[3];
  par[0] = rmin;
  par[1] = rmax;
  par[2] = w / 2;
  mc->Gsvolu(fmd1->GetName(), "TUBE", airId, par, 3);
  mc->Gspos(fmd1->GetName(), fmd1->GetId(), "ALIC", 0, 0, z, 0, "ONLY");

  return DetectorGeometry(fmd1, z);
}

//____________________________________________________________________
Bool_t
AliFMDG3Simulator::FMD2Geometry(AliFMD2* fmd2) 
{
  // Setup the FMD2 geometry.  The FMD2 has no
  // special support as it is at the momement. 
  // 
  // See also AliFMDG3Simulator::DetectorGeometry 
  // 
  if (!fmd2) return kFALSE;
  Double_t rmin     = fmd2->GetInner()->GetLowR();
  Double_t rmax     = fmd2->GetOuterHoneyHighR();
  Double_t hcThick  = fmd2->GetHoneycombThickness();
  Double_t ow       = fmd2->GetInner()->GetRingDepth();
  Double_t iz       = fmd2->GetInnerZ();
  Double_t oz       = fmd2->GetOuterZ();
  Double_t w        = TMath::Abs(oz - iz) + ow + hcThick;
  Double_t z        = oz + w / 2;
  
  TVirtualMC* mc   = TVirtualMC::GetMC();
  Int_t    airId   = (fFMD->GetIdtmed()->At(kAirId));

  Double_t par[3];
  par[0] = rmin;
  par[1] = rmax;
  par[2] = w / 2;
  mc->Gsvolu(fmd2->GetName(), "TUBE", airId, par, 3);
  mc->Gspos(fmd2->GetName(), fmd2->GetId(), "ALIC", 0, 0, z, 0, "ONLY");

  return DetectorGeometry(fmd2, z);
}
  
//____________________________________________________________________
Bool_t
AliFMDG3Simulator::FMD3Geometry(AliFMD3* fmd3) 
{
  // Setup the FMD3 geometry.  The FMD2 has a rather elaborate support
  // structure, as the support will also support the vacuum
  // beam-pipe. 
  // 
  // See also AliFMDG3Simulator::DetectorGeometry 
  // 
  if (!fmd3) return kFALSE;
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
  TVirtualMC* mc   = TVirtualMC::GetMC();
  Int_t    airId   = (fFMD->GetIdtmed()->At(kAirId));
  Int_t    cId     = (fFMD->GetIdtmed()->At(kCarbonId));
  Double_t par[27];
  
  // FMD3 volume 
  par[0]  = 0;
  par[1]  = 360;
  par[2]  = 8;
  // First
  par[3]  = z - nz;
  par[4]  = noser1;
  par[5]  = noser2;
  // Second
  par[6]  = z - (nz - nlen);
  par[7]  = noser1;
  par[8]  = fmd3->ConeR(z - par[6])+.15;
  // Third
  par[9]  = z - innerZ;
  par[10] = innerr1;
  par[11] = fmd3->ConeR(z - par[9])+.15;
  // Fourth
  par[12] = z - innerZh;
  par[13] = innerr1;
  par[14] = fmd3->ConeR(z - par[12])+.15;
  // Fifth
  par[15] = par[12];
  par[16] = outerr1;
  par[17] = fmd3->ConeR(z - par[15])+.15;
  // Sixth
  par[18] = z - nz + zdist + nlen;
  par[19] = outerr1;
  par[20] = fmd3->ConeR(z - par[18])+.15;
  // Seventh
  par[21] = z - nz + nlen + zdist;
  par[22] = outerr1;
  par[23] = flanger+1.5;
  // Eight
  par[24] = z - minZ;
  par[25] = outerr1;
  par[26] = flanger+1.5;
  mc->Gsvolu(fmd3->GetName(), "PCON", airId, par, 27);
  
  Int_t id;
  mc->Matrix(id, 270, 180, 90, 90, 180, 0);
  mc->Gspos(fmd3->GetName(), fmd3->GetId(), "ALIC", 0, 0, z, id, "ONLY");
  
  // Nose volume 
  par[0] = noser1;
  par[1] = noser2;
  par[2] = nlen / 2;
  zi = z - nz + nlen / 2;
  mc->Gsvolu(fgkNoseName, "TUBE", cId, par, 3);
  mc->Gspos(fgkNoseName, 0, fmd3->GetName(), 0, 0, zi, 0, "MANY");
  
  // Back
  par[0] = backr1;
  par[1] = backr2;
  par[2] = backl / 2;
  zi = z - nz + conel - backl / 2;
  mc->Gsvolu(fgkBackName, "TUBE", cId, par, 3);
  mc->Gspos(fgkBackName, 0, fmd3->GetName(), 0, 0, zi, 0, "ONLY");
  
  Int_t n;
  Double_t r;
  // The flanges 
  par[0] = (flanger - backr2) / 2;
  par[1] = fmd3->GetBeamWidth() / 2;
  par[2] = backl / 2;
  mc->Gsvolu(fgkFlangeName, "BOX", cId, par, 3);
  n = fmd3->GetNFlange();
  r = backr2 + (flanger - backr2) / 2;
  for (Int_t i = 0; i  < n; i++) {
    Double_t phi = 360. / n * i + 180. / n;
    Double_t x   = r * TMath::Cos(TMath::Pi() / 180 * phi);
    Double_t y   = r * TMath::Sin(TMath::Pi() / 180 * phi);
    Int_t    id;
    mc->Matrix(id, 90, phi, 90, 90 + phi, 0, 0);
    mc->Gspos(fgkFlangeName, i, fmd3->GetName(), x, y, zi, id, "ONLY");
  }

  // The Beams 
  par[0] = fmd3->GetBeamThickness() / 2;
  par[1] = fmd3->GetBeamWidth() / 2;	
  par[2] = beaml / 2;		   
  mc->Gsvolu(fgkBeamName, "BOX", cId, par, 3);
  n      = fmd3->GetNBeam();
  r      = noser2 + tdist / 2;
  zi     = z - nz + nlen + zdist / 2;
  for (Int_t i = 0; i  < n; i++) {
    Double_t phi = 360. / n * i;
    Double_t x   = r * TMath::Cos(TMath::Pi() / 180 * phi);
    Double_t y   = r * TMath::Sin(TMath::Pi() / 180 * phi);
    Int_t    id;
    (void)theta;
    mc->Matrix(id, 90-theta, phi, 90, 90 + phi, 360 - theta, phi);
    mc->Gspos(fgkBeamName, i, fmd3->GetName(), x, y, zi, id, "MANY");
  }
  
  return DetectorGeometry(fmd3, z);
}

//____________________________________________________________________
void
AliFMDG3Simulator::DefineGeometry() 
{
  // Setup up the FMD geometry. 
  AliDebug(10, "Setting up volume");

  AliFMDGeometry* fmd   = AliFMDGeometry::Instance();
  if (!RingGeometry(fmd->GetInner())) {
    AliError("Failed to create inner ring volume");
    return;
  }
  if (!RingGeometry(fmd->GetOuter())) {
    AliError("Failed to create outer ring volume");
    return;
  }
  FMD1Geometry(fmd->GetFMD1());
  FMD2Geometry(fmd->GetFMD2());
  FMD3Geometry(fmd->GetFMD3());
}

//____________________________________________________________________
//
// EOF
//
