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
//                                 +----------------------+
//				   | AliFMDG3OldSimulator |
//				   +----------------------+
//      
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
// *  AliFMDG3OldSimulator
//    This is a concrete implementation of the AliFMDSimulator that
//    uses the TVirtualMC interface with GEANT 3.21-like messages.
//    This implements the active volume as a divided TUBS shape.  Hits
//    in the corners should be cut away at run time (but currently
//    isn't). 
//
#include <math.h>
#include "AliFMDGeoOldSimulator.h" // ALIFMDG3OLDSIMULATOR_H
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
#include <TGeoVolume.h>		// ROOT_TGeoVolume
#include <TGeoTube.h>		// ROOT_TGeoTube
#include <TGeoMatrix.h>         // ROOT_TGeoMatrix
 
//====================================================================
ClassImp(AliFMDGeoOldSimulator)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDGeoOldSimulator::AliFMDGeoOldSimulator() 
{
  // Default constructor
  fSectorOff   = 1;
  fModuleOff   = -1;
  fRingOff     = 3;
  fDetectorOff = 4;
  fUseDivided  = kTRUE;
}

//____________________________________________________________________
AliFMDGeoOldSimulator::AliFMDGeoOldSimulator(AliFMD* fmd, Bool_t detailed) 
  : AliFMDGeoSimulator(fmd, detailed)
{
  // Normal constructor
  // 
  // Parameters: 
  // 
  //      fmd		Pointer to AliFMD object 
  //      detailed      Whether to make a detailed simulation or not 
  // 
  fSectorOff   = 1;
  fModuleOff   = -1;
  fRingOff     = 3;
  fDetectorOff = 4;
  fUseDivided  = kTRUE;
}

//____________________________________________________________________
TGeoVolume*
AliFMDGeoOldSimulator::RingGeometry(AliFMDRing* r) 
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
  Char_t      id          = r->GetId();
  Double_t    siThick     = r->GetSiThickness();
  // const Int_t nv       = r->GetNVerticies();
  TVector2*   a           = r->GetVertex(5);
  TVector2*   b           = r->GetVertex(3);
  TVector2*   c           = r->GetVertex(4);
  Double_t    theta       = r->GetTheta();
  Double_t    off         = (TMath::Tan(TMath::Pi() * theta / 180) 
			     * r->GetBondingWidth());
  Double_t    rmax        = b->Mod();
  Double_t    rmin        = r->GetLowR();
  Double_t    pcbThick    = r->GetPrintboardThickness();
  Double_t    copperThick = r->GetCopperThickness(); // .01;
  Double_t    chipThick   = r->GetChipThickness(); // .01;
  Double_t    modSpace    = r->GetModuleSpacing();
  Double_t    legr        = r->GetLegRadius();
  Double_t    legl        = r->GetLegLength();
  Double_t    legoff      = r->GetLegOffset();
  Int_t       ns          = r->GetNStrips();
  Int_t       nsec        = Int_t(360 / theta);
  Double_t    stripoff    = a->Mod();
  Double_t    dstrip      = (rmax - stripoff) / ns;
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

  Double_t ringWidth = (siThick + 2 * (pcbThick + copperThick + chipThick));
  // Virtual volume shape to divide - This volume is only defined if
  // the geometry is set to be detailed. 
  // Ring mother volume 
  TGeoShape* ringShape   = new TGeoTube(rmin, rmax, ringWidth / 2);
  name                   =  Form(fgkRingName, id);
  TGeoVolume* ringVolume = new TGeoVolume(name.Data(), ringShape, fAir);

  TGeoShape* activeShape   = new TGeoTube(rmin, rmax, siThick / 2);
  name                     = Form(fgkActiveName, id);
  Double_t z               = - ringWidth / 2 + siThick / 2;
  TGeoVolume* activeVolume = new TGeoVolume(name.Data(), activeShape, fSi);
  ringVolume->AddNode(activeVolume, 1, new TGeoTranslation(0, 0, z));
  
  Int_t sid = activeVolume->GetNumber();
  if (fUseDivided) {
    name                     = Form(fgkSectorName, id);
    TGeoVolume* sectorVolume = activeVolume->Divide(name.Data(), 2, nsec,
						    0, 0, 0, "N");
    
    name                     = Form(fgkStripName, id);
    TGeoVolume* stripVolume  = sectorVolume->Divide(name.Data(), 1, ns,
						    0, 0, 0, "N");
    sid                      = stripVolume->GetNumber();
    AliDebug(10, Form("Got volume id %d for volume %s", sid, name.Data()));
  }
  switch (id) {
  case 'i':
  case 'I': fActiveId[0] = sid; break;
  case 'o':
  case 'O': fActiveId[2] = sid; break;
  }

  // Shape of Printed circuit Board 
  Double_t boardThick = (pcbThick + copperThick + chipThick);
  TGeoShape*  boardShape  =  new TGeoTube(rmin + .1, rmax - .1, boardThick/ 2);
  name                    =  Form(fgkPCBName, id, 'B');
  TGeoVolume* boardVolume =  new TGeoVolume(name.Data(), boardShape, fAir);
  z                       += siThick / 2 + boardThick / 2;
  ringVolume->AddNode(boardVolume,  0, new TGeoTranslation(0, 0, z));
  ringVolume->AddNode(boardVolume,  1, new TGeoTranslation(0,0,z+boardThick));

  // PCB
  TGeoShape*  pcbShape    = new TGeoTube(rmin+.1, rmax-.1, pcbThick / 2);
  name                    =  Form("F%cPC", id);
  z                       =  -boardThick / 2 + pcbThick / 2;
  TGeoVolume* pcbVolume   = new TGeoVolume(name.Data(), pcbShape, fPCB);
  boardVolume->AddNode(pcbVolume, 0, new TGeoTranslation(0, 0, z));

  // Copper
  TGeoShape*  cuShape     =  new TGeoTube(rmin+.1, rmax-.1, copperThick / 2);
  name                    =  Form("F%cCO", id);
  z                       += -pcbThick / 2 + copperThick / 2;
  TGeoVolume* cuVolume    =  new TGeoVolume(name.Data(), cuShape, fCopper);
  boardVolume->AddNode(cuVolume, 0, new TGeoTranslation(0, 0, z));

  // Chip
  TGeoShape*  chipShape   = new TGeoTube(rmin+.1, rmax-.1, chipThick / 2);
  name                    =  Form("F%cCH", id);
  z                       =  -copperThick / 2 + chipThick / 2;
  TGeoVolume* chipVolume  = new TGeoVolume(name.Data(), chipShape, fChip);
  boardVolume->AddNode(chipVolume, 0, new TGeoTranslation(0, 0, z));

  return ringVolume;
}

//____________________________________________________________________
//
// EOF
//
