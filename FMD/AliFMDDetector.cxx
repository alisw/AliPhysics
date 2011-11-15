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
#include <TGeoPhysicalNode.h>   // ROOT_TGeoPhysicalNode
#include <TGeoMatrix.h>		// ROOT_TGeoMatrix 
#include <TMath.h>              // ROOT_TMath

#include "AliFMDDetector.h"	// ALIFMDSUBDETECTOR_H
#include "AliFMDRing.h"		// ALIFMDRING_H
#include "AliFMDDebug.h"		// ALIFMDDEBUG_H ALILOG_H

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
  if (&other == this) return *this; 
  SetName(other.GetName());
  SetTitle(other.GetTitle());
  fId              = other.fId;
  fInner           = other.fInner;
  fOuter           = other.fOuter;
  fInnerTransforms = other.fInnerTransforms;
  fOuterTransforms = other.fOuterTransforms;
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
#define IS_NODE_SENSOR(name)				\
  (name[0] == 'F' && (name[2] == 'B' || name[2] == 'F') && name[3] == 'H')
//#define IS_NODE_SENSOR(name)				
//  (name[0] == 'F' && name[2] == 'S' && name[3] == 'E')
#define IS_NODE_HALF(name) \
  (name[0] == 'F' && name[2] == 'M' && (name[3] == 'B' || name[3] == 'T'))
#define HALF_FORMAT   "FMD/FMD%d_%c"
#define SENSOR_FORMAT "FMD/FMD%d_%c/FMD%c_%02d"

//____________________________________________________________________
void
AliFMDDetector::InitTransformations()
{
  // Find all local<->global transformations for this detector. 
  if ((!fInner || (fInner && fInnerTransforms)) && 
      (!fOuter || (fOuter && fOuterTransforms))) {
    AliFMDDebug(5, ("Transforms for FMD%d already registered", fId));
    return;
  }
  AliFMDDebug(5, ("Initializing transforms for FMD%d", fId));
  if (!gGeoManager) {
    AliFatal("No TGeoManager defined");
    return;
  }

  // Implementation using alignable volume names. 
  // Make container of transforms 
  if (fInner && !fInnerTransforms) 
    fInnerTransforms = new TObjArray(fInner->GetNModules());
  if (fOuter && !fOuterTransforms) 
    fOuterTransforms = new TObjArray(fOuter->GetNModules());
  
  // Loop over bottom/top 
  for (size_t ihalf = 0; ihalf < 2; ihalf++) {
    char  half = (ihalf == 0 ? 'T' : 'B');
    TString path(Form(HALF_FORMAT, fId, half));
    TGeoPNEntry* entry = gGeoManager->GetAlignableEntry(path.Data());
    if (!entry) {
      AliError(Form("Alignable entry for half-detector \"%s\" not found!", 
		    path.Data()));
      continue;
    }
    TGeoPhysicalNode* pn = entry->GetPhysicalNode();
    if (!pn) {
      AliWarning(Form("Making physical volume for \"%s\"", path.Data()));
      pn = gGeoManager->MakeAlignablePN(entry);
      if (!pn) {
	AliError(Form("No physical node for \"%s\"", path.Data()));
	continue;
      }
    }
  }
  
  // Loop over rings 
  for (size_t iring = 0; iring < 2; iring++) {
    char ring = (iring == 0 ? 'I' : 'O');
    TObjArray*  trans = 0;
    AliFMDRing* r     = 0; 
    switch (ring) {
    case 'I': r = fInner; trans = fInnerTransforms; break;
    case 'O': r = fOuter; trans = fOuterTransforms; break; 
    }
    if (!r || !trans) continue;

    Int_t nModules = r->GetNModules();
    if (nModules <= 0) continue;

    // Loop over bottom/top 
    for (size_t ihalf = 0; ihalf < 2; ihalf++) {
      char  half = (ihalf == 0 ? 'T' : 'B');
      Int_t base = (half == 'T' ? 0 : nModules / 2);
      
      // Loop over modules in this half ring 
      for (Int_t imod = 0; imod < nModules / 2; imod++) {
	// Find physical node entry
	TString path(Form(SENSOR_FORMAT, fId, half, ring, base+imod));
	TGeoPNEntry* entry = gGeoManager->GetAlignableEntry(path.Data());
	if (!entry) {
	  AliError(Form("Alignable entry for sensor \"%s\" not found!", 
			path.Data()));
	  continue;
	}
	TGeoPhysicalNode* pn = entry->GetPhysicalNode();
	if (!pn) {
	  AliWarning(Form("Making physical volume for \"%s\"", path.Data()));
	  pn = gGeoManager->MakeAlignablePN(entry);
	  if (!pn) {
	    AliError(Form("No physical node for \"%s\"", path.Data()));
	    continue;
	  }
	}
	
	const TGeoMatrix* pm = pn->GetMatrix();
	if (!pm) {
	  AliError(Form("No matrix for path \"%s\"", path.Data()));
	  continue;
	}
	// Get transformation matrix for this node, and store it. 
	TGeoMatrix*  t = new TGeoHMatrix(*pm);
	trans->AddAt(t, base+imod);
	AliFMDDebug(5, ("Found matrix for path \"%s\": %p",path.Data(),pm));
      }
    }
  }
  if (HasAllTransforms('I') && HasAllTransforms('O')) return;

  // Alternative implementation using TGeoIter. 
  TGeoVolume* topVolume = gGeoManager->GetTopVolume();
  if (!topVolume) {
    AliFatal("No top-level volume defined");
    return;
  }
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
    AliFMDDebug(50, ("Got volume %s", name));
    // Check if this node is this detector 
    // The base offset for numbers in the ASCII table is 48
    if (IS_NODE_THIS(name)) {
      AliFMDDebug(20, ("Found detector node '%s' for FMD%d", name, fId));
      thisNodeFound = kTRUE;
    }
    // if the detector was found, then we're on that branch, and we
    // check if this node represents a module in that branch.
    if (thisNodeFound && IS_NODE_SENSOR(name)) {
      AliFMDDebug(20, ("Found sensor node '%s' for FMD%d", name, fId));
      // Get the ring Id.
      Char_t ringid = name[1];

      // Get the approprate ring
      AliFMDRing* ring = GetRing(ringid);
      if (!ring) continue;

      // Check whether we have all the modules we need for this ring,
      // and if so, go on to the next node. 
      Bool_t& done = (ring == fInner ? allInners : allOuters);
      if ((done = HasAllTransforms(ringid))) {
	AliFMDDebug(20, ("Already has all module transforms for ring %c", 
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
  // Set alignable volumes. 
  // This will define the alignable volumes. 
  // That is currently, the modules and the half-rings. 
  
  AliFMDDebug(10, ("Making alignable volumes for FMD%d", fId));
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
    AliFMDDebug((name[0] == 'F' ? 40 : 50), ("Got volume %s", name));
    // Check if this node is this detector 
    // The base offset for numbers in the ASCII table is 48
    if (IS_NODE_THIS(name)) {
      AliFMDDebug(20, ("Found detector node '%s' for FMD%d", name, fId));
      thisNodeFound = kTRUE;
    }

    // if a half ring is found, then we're on that branch, and we
    // check if this node represents a half ring on that branch 
    if (thisNodeFound && IS_NODE_HALF(name)) {
      AliFMDDebug(30, ("Found half node '%s' for FMD%d", name, fId));
      // Get the half Id.
      thisHalf = name[3];

      // Check if we're done 
      Bool_t done = (thisHalf == 'T' ? hasTop : hasBottom);
      if (done) {
	AliFMDDebug(20, ("Already has all halves for detector %c",name[1]));
	continue;
      }

      switch (thisHalf) {
      case 'T': hasTop = true; break;
      case 'B': hasBottom = true; break;
      default:  
	AliWarning(Form("Unknown part '%c' of FMD%d", thisHalf, fId));
	continue; // because the node is unknown. 
      }
      
      // Get the node path 
      next.GetPath(path);
      align = Form(HALF_FORMAT, fId, thisHalf);
    }
    
    // if the detector was found, then we're on that branch, and we
    // check if this node represents a module in that branch.
    if (thisNodeFound && thisHalf && IS_NODE_SENSOR(name)) {
      AliFMDDebug(30, ("Found sensor node '%s' for FMD%d", name, fId));
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
	AliFMDDebug(20, ("Already has all sensor volumes for ring %c",ringid));
	continue;
      }
      // Get the copy (module) number, and check that it hasn't
      // already been added to the container. 
      Int_t copy  = node->GetNumber();
      next.GetPath(path);
      // path.Replace("ALIC", "/ALIC_1");
      align = Form(SENSOR_FORMAT, fId, thisHalf, ringid, copy);
      
      switch (ringid) {
      case 'I': iInnerSensor++; break;
      case 'O': iOuterSensor++; break;
      }
    }
    if (!align.IsNull() && !path.IsNull()) {
      AliFMDDebug(20, ("Got %s -> %s", path.Data(), align.Data()));
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
    AliFMDDebug(20, ("FMD%d: top: %d bottom: %d Inner: %d/%d Outer %d/%d", 
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
  if (!r) { 
    AliWarning(Form("No such ring FMD%d%c ", fId, ring));
    return;
  }
  TGeoMatrix* m = FindTransform(ring, sector);
  if (!m) { 
    AliWarning(Form("No transfrmation found for FMD%d%c[%02d]", 
		    fId, ring, sector));
    return;
  }
  Double_t rho      = r->GetStripRadius(strip);
  Double_t phi      = ((sector % 2) - .5) * r->GetTheta();
  Double_t siThick  = r->GetSiThickness();
#if 0 
  Double_t modThick = (siThick
		       + r->GetPrintboardThickness()
		       + r->GetCopperThickness()
		       + r->GetChipThickness()
		       + r->GetSpacing());
#endif
  AliFMDDebug(30, ("Rho %7.3f, angle %7.3f", rho, phi));
  Double_t local[]  = { rho * TMath::Cos(phi * TMath::DegToRad()), 
		        rho * TMath::Sin(phi * TMath::DegToRad()), 
		        /* -modThick + */ siThick / 2 };
  Double_t master[3];
  AliFMDDebug(30, ("Local (%7.3f,%7.3f,%7.3f)",local[0], local[1], local[2]));
  m->LocalToMaster(local, master);
  AliFMDDebug(30, ("Master (%7.3f,%7.3f,%7.3f)",
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
