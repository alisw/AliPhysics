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
/** @file    AliFMDDetector.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:36:27 2006
    @brief   Sub-detector base class implementation
    @ingroup FMD_base
*/

//____________________________________________________________________
//
// AliFMDDetector.   
//
// Base class for concrete FMD detectors, like AliFMD1, AliFMD2,
// AliFMD3. 
// Utility class to help implement the FMD geometry.  This provides
// the interface for the concrete geometry implementations of the FMD
// sub-detectors. 
//
// The AliFMDGeometry object owns the AliFMDDetector objects
//
// Latest changes by Christian Holm Christensen
//

#include <TGeoManager.h>	// ROOT_TGeoManager 
#include <TGeoMatrix.h>		// ROOT_TGeoMatrix 
#include <TMath.h>              // ROOT_TMath

#include "AliFMDDetector.h"	// ALIFMDSUBDETECTOR_H
#include "AliFMDRing.h"		// ALIFMDRING_H
#include "AliLog.h"             // ALILOG_H

//====================================================================
ClassImp(AliFMDDetector)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDDetector::AliFMDDetector(Int_t id, AliFMDRing* inner, AliFMDRing* outer) 
  : TNamed(Form("FMD%d", id), "Forward multiplicity ring"), 
    fId(id), 
    fInnerZ(0.),
    fOuterZ(0.),
    fHoneycombThickness(0.),
    fAlThickness(0.),
    fInnerHoneyLowR(0.),
    fInnerHoneyHighR(0.),
    fOuterHoneyLowR(0.),
    fOuterHoneyHighR(0.),
    fInner(inner),
    fOuter(outer), 
    fInnerTransforms(0),
    fOuterTransforms(0)
{
  // Constructor
  // 
  //   ID         Id of detector (1,2, or 3)
  //   INNER      Inner ring geometry 
  //   OUTER      Outer ring geometry (if any)
  // 
  SetHoneycombThickness();
  SetAlThickness();
  SetInnerHoneyLowR(0);
  SetInnerHoneyHighR(0);
  SetInnerZ(0);
  SetOuterZ(0);
  SetOuterHoneyLowR(0);
  SetOuterHoneyHighR(0);
}

//____________________________________________________________________
AliFMDDetector::AliFMDDetector(const AliFMDDetector& other)
  : TNamed(other), 
    fId(other.fId),
    fInnerZ(0.),
    fOuterZ(0.),
    fHoneycombThickness(0.),
    fAlThickness(0.),
    fInnerHoneyLowR(0.),
    fInnerHoneyHighR(0.),
    fOuterHoneyLowR(0.),
    fOuterHoneyHighR(0.),
    fInner(other.fInner),
    fOuter(other.fOuter),
    fInnerTransforms(other.fInnerTransforms),
    fOuterTransforms(other.fOuterTransforms)
{
  // Copy constructor 
  SetHoneycombThickness(other.GetHoneycombThickness());
  SetAlThickness(other.GetAlThickness());
  SetInnerHoneyLowR(other.GetInnerHoneyLowR());
  SetInnerHoneyHighR(other.GetInnerHoneyHighR());
  SetInnerZ(other.GetInnerZ());
  SetOuterZ(other.GetOuterZ());
  SetOuterHoneyLowR(other.GetOuterHoneyLowR());
  SetOuterHoneyHighR(other.GetOuterHoneyHighR());
}

//____________________________________________________________________
AliFMDDetector&
AliFMDDetector::operator=(const AliFMDDetector& other)
{
  // Assignment operator
  SetName(other.GetName());
  SetTitle(other.GetTitle());
  fId              = other.fId;
  fInner           = other.fInner;
  fOuter           = other.fOuter;
  fInnerTransforms = other.fInnerTransforms;
  fOuterTransforms = other.fOuterTransforms;
  SetHoneycombThickness(other.GetHoneycombThickness());
  SetAlThickness(other.GetAlThickness());
  SetInnerHoneyLowR(other.GetInnerHoneyLowR());
  SetInnerHoneyHighR(other.GetInnerHoneyHighR());
  SetInnerZ(other.GetInnerZ());
  SetOuterZ(other.GetOuterZ());
  SetOuterHoneyLowR(other.GetOuterHoneyLowR());
  SetOuterHoneyHighR(other.GetOuterHoneyHighR());
  return *this;
}

//____________________________________________________________________
void
AliFMDDetector::Init()
{
  // Initialize. 
  if (fInner) {
    SetInnerHoneyLowR(fInner->GetLowR() + 1.);
    SetInnerHoneyHighR(fInner->GetHighR() + 1.);
  }
  if (fOuter) {
    SetOuterHoneyLowR(fOuter->GetLowR() + 1.);
    SetOuterHoneyHighR(fOuter->GetHighR() + 1.);
  }  
}

//____________________________________________________________________
Bool_t
AliFMDDetector::HasAllTransforms(Char_t ring) const
{
  // Check if we got all transformations for a given ring.  Return
  // true in that case. 
  AliFMDRing* r = GetRing(ring);
  if (!r) return kTRUE;
  TObjArray* matricies = (r == fInner ? fInnerTransforms : fOuterTransforms);
  if (!matricies) return kTRUE;
  if (matricies->GetEntries() == r->GetNModules()) return kTRUE;
  return kFALSE;
}

#define IS_NODE_THIS(name) \
  (name[0] == 'F' && name[2] == 'M' && name[1] == Char_t(48+fId) && \
   (name[3] == 'T' || name[3] == 'B'))
#define IS_NODE_SENSOR(name) \
  (name[0] == 'F' && name[2] == 'S' && name[3] == 'E')
#define IS_NODE_HALF(name) \
  (name[0] == 'F' && name[2] == 'M' && (name[3] == 'B' || name[3] == 'T'))

//____________________________________________________________________
void
AliFMDDetector::InitTransformations()
{
  // Find all local<->global transformations for this detector. 
  if ((!fInner || (fInner && fInnerTransforms)) && 
      (!fOuter || (fOuter && fOuterTransforms))) {
    AliDebug(5, Form("Transforms for FMD%d already registered", fId));
    return;
  }
  AliDebug(5, Form("Initializing transforms for FMD%d", fId));
  if (!gGeoManager) {
    AliFatal("No TGeoManager defined");
    return;
  }
  TGeoVolume* topVolume = gGeoManager->GetTopVolume();
  if (!topVolume) {
    AliFatal("No top-level volume defined");
    return;
  }
  // Make container of transforms 
  if (fInner && !fInnerTransforms) 
    fInnerTransforms = new TObjArray(fInner->GetNModules());
  if (fOuter && !fOuterTransforms) 
    fOuterTransforms = new TObjArray(fOuter->GetNModules());
  
  // Make an iterator
  TGeoIterator next(topVolume);
  TGeoNode* node = 0;
  
  // Find the node corresponding to this detector, and then find the
  // sensor volumes 
  Bool_t thisNodeFound = kFALSE;
  Bool_t allInners     = HasAllTransforms('I');
  Bool_t allOuters     = HasAllTransforms('O');
  
  while ((node = static_cast<TGeoNode*>(next())) 
	 && !(allInners && allOuters)) {
    // Get nodes names 
    const Char_t* name = node->GetName();
    if (!name) continue;
    AliDebug(50, Form("Got volume %s", name));
    // Check if this node is this detector 
    // The base offset for numbers in the ASCII table is 48
    if (IS_NODE_THIS(name)) {
      AliDebug(20, Form("Found detector node '%s' for FMD%d", name, fId));
      thisNodeFound = kTRUE;
    }
    // if the detector was found, then we're on that branch, and we
    // check if this node represents a module in that branch.
    if (thisNodeFound && IS_NODE_SENSOR(name)) {
      AliDebug(20, Form("Found sensor node '%s' for FMD%d", name, fId));
      // Get the ring Id.
      Char_t ringid = name[1];

      // Get the approprate ring
      AliFMDRing* ring = GetRing(ringid);
      if (!ring) continue;

      // Check whether we have all the modules we need for this ring,
      // and if so, go on to the next node. 
      Bool_t& done = (ring == fInner ? allInners : allOuters);
      if ((done = HasAllTransforms(ringid))) {
	AliDebug(20, Form("Already has all module transforms for ring %c", 
			 ringid));
	continue;
      }

      // Get the approprate container
      TObjArray* matricies = (ringid == 'i' || ringid == 'I' 
			      ? fInnerTransforms : fOuterTransforms);

      // Get the copy (module) number, and check that it hasn't
      // already been added to the container. 
      Int_t copy  = node->GetNumber();
      if (matricies->At(copy)) {
	AliWarning(Form("Have a transformation for module %d in ring %c", 
			copy, ringid));
	continue;
      }

      // Get the global transformation matrix, and store it. 
      TGeoMatrix*  trans = new TGeoHMatrix(*(next.GetCurrentMatrix()));
      matricies->AddAt(trans, copy);

    }
  }
}

//____________________________________________________________________
void
AliFMDDetector::SetAlignableVolumes() const
{
  AliDebug(10, Form("Making alignable volumes for FMD%d", fId));
  if (!gGeoManager) {
    AliFatal("No TGeoManager defined");
    return;
  }
  TGeoVolume* topVolume = gGeoManager->GetTopVolume();
  if (!topVolume) {
    AliFatal("No top-level volume defined");
    return;
  }

  // Make an iterator
  TGeoIterator next(topVolume);
  next.Reset(topVolume);
  next.SetTopName(Form("/%s_1", topVolume->GetName()));
  TGeoNode* node = 0;
  
  Int_t nInnerSensor = (fInner ? fInner->GetNModules() : 0);
  Int_t nOuterSensor = (fOuter ? fOuter->GetNModules() : 0);
  // Find the node corresponding to this detector, and then find the
  // sensor volumes 
  Bool_t thisNodeFound = kFALSE;
  Char_t thisHalf      = '\0';
  Int_t  iInnerSensor  = 0;
  Int_t  iOuterSensor  = 0;
  Bool_t hasTop        = false;
  Bool_t hasBottom     = false;
  
  TString path, align;
  while ((node = static_cast<TGeoNode*>(next())) 
	 && (iInnerSensor < nInnerSensor || iOuterSensor < nOuterSensor
	     || !hasBottom || !hasTop)) {
    // Get nodes names 
    const Char_t* name = node->GetName();
    if (!name) continue;
    AliDebug((name[0] == 'F' ? 40 : 50), Form("Got volume %s", name));
    // Check if this node is this detector 
    // The base offset for numbers in the ASCII table is 48
    if (IS_NODE_THIS(name)) {
      AliDebug(20, Form("Found detector node '%s' for FMD%d", name, fId));
      thisNodeFound = kTRUE;
    }

    // if a half ring is found, then we're on that branch, and we
    // check if this node represents a half ring on that branch 
    if (thisNodeFound && IS_NODE_HALF(name)) {
      AliDebug(30, Form("Found half node '%s' for FMD%d", name, fId));
      // Get the half Id.
      thisHalf = name[3];

      // Check if we're done 
      Bool_t done = (thisHalf == 'T' ? hasTop : hasBottom);
      if (done) {
	AliDebug(20,Form("Already has all halves for detector %c",name[1]));
	continue;
      }

      switch (thisHalf) {
      case 'T': hasTop = true; break;
      case 'B': hasBottom = true; break;
      default:  
	AliWarning(Form("Unknown part '%c' of FMD%d", fId));
	continue; // because the node is unknown. 
      }
      
      // Get the node path 
      next.GetPath(path);
      align = Form("FMD/FMD%d_%c", fId, thisHalf);
    }
    
    // if the detector was found, then we're on that branch, and we
    // check if this node represents a module in that branch.
    if (thisNodeFound && thisHalf && IS_NODE_SENSOR(name)) {
      AliDebug(30, Form("Found sensor node '%s' for FMD%d", name, fId));
      // Get the ring Id.
      Char_t ringid = name[1];

      // check that the ring is valid 
      if (!GetRing(ringid)) {
	AliWarning(Form("Invalid ring %c for FMD%d", ringid, fId));
	continue;
      }

      // Check if we're done
      Bool_t done = false;
      switch (ringid) {
      case 'I': done = iInnerSensor >= nInnerSensor; break;
      case 'O': done = iOuterSensor >= nOuterSensor; break;
      default: continue;
      }
      if (done) {
	AliDebug(20,Form("Already has all sensor volumes for ring %c",ringid));
	continue;
      }
      // Get the copy (module) number, and check that it hasn't
      // already been added to the container. 
      Int_t copy  = node->GetNumber();
      next.GetPath(path);
      // path.Replace("ALIC", "/ALIC_1");
      align = Form("FMD/FMD%d_%c/FMD%c_%02d", fId, thisHalf, ringid, copy);
      
      switch (ringid) {
      case 'I': iInnerSensor++; break;
      case 'O': iOuterSensor++; break;
      }
    }
    if (!align.IsNull() && !path.IsNull()) {
      AliDebug(20, Form("Got %s -> %s", path.Data(), align.Data()));
      TGeoPNEntry* entry = 
	gGeoManager->SetAlignableEntry(align.Data(),path.Data());
      if(!entry)
	AliFatal(Form("Alignable entry %s not created. "
		      "Volume path %s not valid", 
			align.Data(),path.Data()));
#ifdef MAKE_ALIGNABLE_PHYSICAL
      TGeoPhysicalNode* phys = gGeoManager->MakeAlignablePN(entry);
      if (!phys) 
	AliWarning(Form("Physical node entry %s not created. "
			"Volume path %s not valid", 
			align.Data(),path.Data()));
#endif
      align = "";
    }
    AliDebug(20, Form("FMD%d: top: %d bottom: %d Inner: %d/%d Outer %d/%d", 
		      fId, hasTop, hasBottom, iInnerSensor,  nInnerSensor, 
		      iOuterSensor, nOuterSensor));
  }
}

  

//____________________________________________________________________
AliFMDRing*
AliFMDDetector::GetRing(Char_t id) const
{
  // Get the specified ring 
  // 
  //   ID      Id of ring ('I' or 'O')
  // 
  switch (id) {
  case 'i':
  case 'I': return GetInner();
  case 'o':
  case 'O': return GetOuter();
  }
  return 0;
}

//____________________________________________________________________
Double_t
AliFMDDetector::GetRingZ(Char_t id) const
{
  // Get the z-coordinate specified ring 
  // 
  //   ID      Id of ring ('I' or 'O')
  // 
  switch (id) {
  case 'i':
  case 'I': return GetInnerZ();
  case 'o':
  case 'O': return GetOuterZ();
  }
  return 0;
}

//____________________________________________________________________
TGeoMatrix*
AliFMDDetector::FindTransform(Char_t ring, UShort_t sector) const 
{
  // Find the transformation that corresponds to sector sector in ring
  // ring. 
  TObjArray* matricies = 0;
  switch (ring) {
  case 'i': case 'I': matricies = fInnerTransforms; break;
  case 'o': case 'O': matricies = fOuterTransforms; break;
  }
  if (!matricies) { 
    AliWarning(Form("Unknown ring %c of FMD%d", ring, fId));
    return 0;
  }
  UInt_t module = sector / 2;
  TGeoMatrix* m = static_cast<TGeoMatrix*>(matricies->At(module));
  if (!m) {
    AliWarning(Form("No matrix found for sector %d in FMD%d%c", 
		    sector, fId, ring));
    return 0;
  }
  return m;
}

  
//____________________________________________________________________
void
AliFMDDetector::Detector2XYZ(Char_t   ring, 
			     UShort_t sector,
			     UShort_t strip, 
			     Double_t& x, 
			     Double_t& y, 
			     Double_t& z) const
{
  // Translate detector coordinates (this,ring,sector,strip) into
  // (x,y,z) coordinates (in global reference frame)
  AliFMDRing* r = GetRing(ring);
  if (!r) return;
  TGeoMatrix* m = FindTransform(ring, sector);
  if (!m) return;
  Double_t rho      = r->GetStripRadius(strip);
  Double_t phi      = ((sector % 2) - .5) * r->GetTheta();
  Double_t siThick  = r->GetSiThickness();
  Double_t modThick = (siThick
		       + r->GetPrintboardThickness()
		       + r->GetCopperThickness()
		       + r->GetChipThickness()
		       + r->GetSpacing());
  AliDebug(30,Form("Rho %7.3f, angle %7.3f", rho, phi));
# define DEGRAD TMath::Pi() / 180. 
  Double_t local[]  = { rho * TMath::Cos(phi * DEGRAD), 
		        rho * TMath::Sin(phi * DEGRAD), 
		        -modThick + siThick / 2 };
  Double_t master[3];
  AliDebug(30, Form("Local (%7.3f,%7.3f,%7.3f)",local[0], local[1], local[2]));
  m->LocalToMaster(local, master);
  AliDebug(30, Form("Master (%7.3f,%7.3f,%7.3f)",
		    master[0],master[1],master[2]));
  x = master[0];
  y = master[1];
  z = master[2];
}

//____________________________________________________________________
Bool_t
AliFMDDetector::XYZ2Detector(Double_t  x,
			     Double_t  y,
			     Double_t  z,
			     Char_t&   ring, 
			     UShort_t& sector,
			     UShort_t& strip) const
{
  // Translate (x,y,z) coordinates (in global reference frame) into 
  // detector coordinates (this,ring,sector,strip).
  AliFMDRing* rng = 0;
  ring = -1;
  for (int j = 0; j < 2; j++) {
    rng = GetRing(j == 0 ? 'I'  : 'O');
    if (!rng) continue;
    Double_t ringZ    = GetRingZ(j == 0 ? 'I'  : 'O');
    Double_t modSpace = TMath::Sign(rng->GetModuleSpacing(), ringZ);
    if (TMath::Abs(z - ringZ) < 0.01 || 
	TMath::Abs(z - ringZ + modSpace) < 0.01) break;
    rng = 0;
  }
  if (rng && rng->XYZ2Detector(x, y, z - GetRingZ(rng->GetId()),
			       sector, strip)) {
    ring = rng->GetId();
    return kTRUE;
  }
  return kFALSE;
}

  

//____________________________________________________________________
// 
// EOF
//
