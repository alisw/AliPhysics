/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved.      *
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

//////////////////////////////////////////////////////////////////////////////
//
// Utility class to help implement the FMD geometry.  This provides
// the interface for the concrete geometry implementations of the FMD
// sub-detectors. 
//
// The AliFMD object owns the AliFMDSubDetector objects
//
// Latest changes by Christian Holm Christensen
//
//////////////////////////////////////////////////////////////////////////////

#include "AliFMDSubDetector.h"	// ALIFMDSUBDETECTOR_H
#include "AliFMDRing.h"		// ALIFMDRING_H
#include <AliLog.h>		// ALILOG_H
#include <TVirtualMC.h>		// ROOT_TVirtualMC
#include <TList.h>		// ROOT_TList
#include <TString.h>		// ROOT_TString

ClassImp(AliFMDSubDetector);

//____________________________________________________________________
AliFMDSubDetector::AliFMDSubDetector(Int_t n)
  : fId(n), 
    fInnerZ(0), 
    fOuterZ(0), 
    fInner(0), 
    fOuter(0)
{}

//____________________________________________________________________
void 
AliFMDSubDetector::Draw(Option_t* /* opt */) const 
{
  // DebugGuard guard("AliFMDSubDetector::Draw");
  AliDebug(10, "AliFMDSubDetector::Draw");
}

//____________________________________________________________________
Bool_t
AliFMDSubDetector::CheckHit(Char_t ring, Int_t module, Double_t x, Double_t y) 
{
  // Check if a hit (x,y) in module module of ring ring is within the
  // actual shape. 
  Bool_t ret = kFALSE;
  switch (ring) {
  case 'i':
  case 'I': 
    if (!fInner) break;
    ret = fInner->IsWithin(module, x, y);
    break;
  case 'o':
  case 'O': 
    if (!fOuter) break;
    ret = fOuter->IsWithin(module, x, y);
    break;
  }
  return ret;
}

//____________________________________________________________________
void 
AliFMDSubDetector::SimpleGeometry(TList* nodes, 
				  TNode* mother, 
				  Int_t colour, 
				  Double_t zMother) 
{
  // Make a simplified geometry for event display 
  // 
  // Parameters
  // 
  //    nodes     List of nodes to register all create nodes in 
  //    mother    Mother node to put the ring in. 
  //    colour    Colour of the nodes 
  //    zMother   Z position of the node in the mother volume 
  // 
  for (int i = 0; i < 2; i++) {
    AliFMDRing* r = 0;
    Double_t z = 0;
    switch (i) {
    case 0: 
      r     = fInner;
      z     = fInnerZ;
      break;
    case 1: 
      r     =  fOuter;
      z     =  fOuterZ;
      break;
    }
    if (!r) continue;

    // Make the coordinates relative to the mother volume.   If we're
    // on the positive side, then we need to flip the z-coordinate, as
    // we'll rotate the whole sub-detector afterwards. 
    z -= zMother;
    if (zMother > 0) z *= -1;
    
    r->SimpleGeometry(nodes, mother, colour, z, fId);
  }
}

  
//____________________________________________________________________
void 
AliFMDSubDetector::SetupGeometry(Int_t airId, Int_t kaptionId) 
{
  // Set up the geometry of this particular detector. 
  // 
  // In this class, it defines the support honey comp calls and
  // nothing else. 
  // 
  // Parameters
  //   airId           Medium of inactive virtual volumes 
  //   kaptionId       Medium of honeycomb
  // 
  // DebugGuard guard("AliFMDSubDetector::SetupGeometry");
  AliDebug(10, "AliFMDSubDetector::SetupGeometry");
  TString name;
  Double_t par[5];

  for (int i = 0; i < 2; i++) {
    AliFMDRing* r       = 0;
    char  c = '\0';
    switch (i) {
    case 0: 
      r      = fInner;
      c      = 'I';
      par[0] = fInnerHoneyLowR;
      par[1] = fInnerHoneyHighR;
      break;
    case 1: 
      r     = fOuter;
      c     = 'O';
      par[0] = fOuterHoneyLowR;
      par[1] = fOuterHoneyHighR;
      break;
    }
    if (!r) continue;
    // Top of honeycomb, inner ring 
    par[2] = fHoneycombThickness / 2;
    par[3] = 0;
    par[4] = 180;
    name   = Form("HT%c%d", c, fId);
    gMC->Gsvolu(name.Data(), "TUBS", kaptionId, par, 5);
    
    // Bottom of honeycomb, inner ring 
    par[3] = 180;
    par[4] = 360;
    name   = Form("HB%c%d", c, fId);
    gMC->Gsvolu(name.Data(), "TUBS", kaptionId, par, 5);
    
    // Air in top of honeycomb, inner ring 
    par[0] += fKaptionThickness;
    par[1] -= fKaptionThickness;
    par[2] -= fKaptionThickness;
    par[3] = 0;
    par[4] = 180;
    name   = Form("GT%c%d", c, fId);
    gMC->Gsvolu(name.Data(), "TUBS", airId, par, 5);
    
    // Air in bottom of honeycomb, inner ring 
    par[3] = 180;
    par[4] = 360;
    name   = Form("GB%c%d", c, fId);
    gMC->Gsvolu(name.Data(), "TUBS", airId, par, 5);
  }
}

//____________________________________________________________________
void 
AliFMDSubDetector::Geometry(const char* mother, Int_t pbRotId, Int_t idRotId, 
		      Double_t zMother) 
{
  // Place the volume inside mother volume. 
  // 
  // Parameters
  // 
  //     mother     Volume to position this detector in 
  //     pbRotId    Print board rotation matrix, 
  //     idRotId    Identity rotation matrix 
  //     zMother    The Z passed in, is the position of the middle
  //                point of the mother volume. 
  // 
  // In this base class, it asks the contained rings to position
  // themselves in the mother volume, and positions the honey comb
  // support in the mother volume 
  // 
  // DebugGuard guard("AliFMDSubDetector::Geometry");
  AliDebug(10, "AliFMDSubDetector::Geometry");

  Double_t  ringW;
  Double_t  z = 0;
  // Double_t* b = 0;
  TString name;
  TString name2;
  
  for (int i = 0; i < 2; i++) {
    AliFMDRing* r = 0;
    char  c = '\0';
    switch (i) {
    case 0: 
      r     = fInner;
      c     = 'I';
      z     = fInnerZ;
      break;
    case 1: 
      r     =  fOuter;
      c     =  'O';
      z     =  fOuterZ;
      break;
    }
    if (!r) continue;

    // Make the coordinates relative to the mother volume.   If we're
    // on the positive side, then we need to flip the z-coordinate, as
    // we'll rotate the whole sub-detector afterwards. 
    z -= zMother;
    if (zMother > 0) z *= -1;
    
    r->Geometry(mother, fId, z, pbRotId, idRotId);
    ringW =  r->GetRingDepth();
    z     -= ringW + fHoneycombThickness / 2;

    // Top of honeycomb
    name = Form("HT%c%d", c, fId);
    gMC->Gspos(name.Data(), 1, mother, 0, 0, z, idRotId);

    // Air in top of honeycomb
    name2 = name;
    name  = Form("GT%c%d", c, fId);
    gMC->Gspos(name.Data(), 1, name2.Data(),0,fKaptionThickness,0,idRotId);
    
    // Bottom of honeycomb
    name = Form("HB%c%d", c, fId);
    gMC->Gspos(name.Data(), 1, mother, 0, 0, z, idRotId);

    // Air in bottom of honeycomb
    name2 = name;
    name  = Form("GB%c%d", c, fId);
    gMC->Gspos(name.Data(),1,name2.Data(),0,-fKaptionThickness,0,idRotId);
  }
}

//____________________________________________________________________
void
AliFMDSubDetector::Gsatt() 
{
  // Set drawing attributes for the detector 
  // 
  // DebugGuard guard("AliFMDSubDetector::Gsatt");
  AliDebug(10, "AliFMDSubDetector::Gsatt");
  TString name(Form("FMD%d", fId));
  gMC->Gsatt(name.Data(), "SEEN", 0);
}


//____________________________________________________________________
// 
// EOF
//
