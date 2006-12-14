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
/** @file    AliFMDGeometry.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:40:37 2006
    @brief   Geometry mananger for the FMD
*/
//____________________________________________________________________
//                                                                          
// Forward Multiplicity Detector based on Silicon wafers. 
//
// This class is a singleton that handles the geometry parameters of
// the FMD detectors.  
//                                                       
// The actual code is done by various separate classes.   Below is
// diagram showing the relationship between the various FMD classes
// that handles the geometry 
//
//                               +------------+ 
//                            +- | AliFMDRing |
// 			   2  |  +------------+
//      +----------------+<>--+        |				
//      | AliFMDGeometry |             ^                       	
//      +----------------+<>--+        V 1..2                     	
//           		   3  | +----------------+ 		
//            		      +-| AliFMDDetector | 		
//             		        +----------------+		
//                                     ^
//                                     |
//                       +-------------+-------------+
//                       |             |             |	      
//                  +---------+   +---------+   +---------+
//                  | AliFMD1 |   | AliFMD2 |   | AliFMD3 |
//                  +---------+   +---------+   +---------+
//      
//
// *  AliFMDRing 
//    This class contains all stuff needed to do with a ring.  It's
//    used by the AliFMDDetector objects to instantise inner and
//    outer rings.  The AliFMDRing objects are shared by the
//    AliFMDDetector objects, and owned by the AliFMDv1 object. 
//
// *  AliFMD1, AliFMD2, and AliFMD3 
//    These are specialisation of AliFMDDetector, that contains the
//    particularities of each of the sub-detector system.  It is
//    envisioned that the classes should also define the support
//    volumes and material for each of the detectors.                          
//                                                                          
//
#include "AliFMDGeometry.h"	// ALIFMDGEOMETRY_H
#include "AliFMDRing.h"		// ALIFMDRING_H
#include "AliFMD1.h"		// ALIFMD1_H
#include "AliFMD2.h"		// ALIFMD2_H
#include "AliFMD3.h"		// ALIFMD2_H
#include "AliRecPoint.h"	// ALIRECPOINT_H
#include "AliLog.h"		// ALIRECPOINT_H
#include <TVector3.h>           // ROOT_TVector3
// #include <TMatrix.h>            // ROOT_TMatrix
// #include <TParticle.h>          // ROOT_TParticle
#include <Riostream.h>
#include "AliFMDGeometryBuilder.h"
// #include <TArrayI.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoNode.h>
static Int_t FindNodeDepth(const char* name, const char* volname);


//====================================================================
ClassImp(AliFMDGeometry)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDGeometry* AliFMDGeometry::fgInstance = 0;

//____________________________________________________________________
AliFMDGeometry* 
AliFMDGeometry::Instance() 
{
  // Return (newly created) singleton instance 
  if (!fgInstance) fgInstance = new AliFMDGeometry;
  return fgInstance;
}

//____________________________________________________________________
AliFMDGeometry::AliFMDGeometry() 
  : AliGeometry("FMD", "Forward multiplicity"), 
    fIsInitialized(kFALSE), 
    fInner(0),
    fOuter(0),
    fFMD1(0),
    fFMD2(0),
    fFMD3(0),
    fUseFMD1(kFALSE),
    fUseFMD2(kFALSE),
    fUseFMD3(kFALSE),
    fBuilder(0),
    fDetectorOff(0),
    fModuleOff(0),  
    fRingOff(0),
    fSectorOff(0),
    fActive(2),
    fDetailed(kFALSE),       
    fUseAssembly(kFALSE)
{
  // PROTECTED
  // Default constructor 
  fUseFMD1     = kTRUE;
  fUseFMD2     = kTRUE;
  fUseFMD3     = kTRUE;  
  fDetailed    = kTRUE;
  fUseAssembly = kTRUE;
  fInner = new AliFMDRing('I');
  fOuter = new AliFMDRing('O');
  fFMD1  = new AliFMD1(fInner);
  fFMD2  = new AliFMD2(fInner, fOuter);
  fFMD3  = new AliFMD3(fInner, fOuter);
  fIsInitialized = kFALSE;
  fActive.Reset(-1);
}

//____________________________________________________________________
AliFMDGeometry::AliFMDGeometry(const AliFMDGeometry& other) 
  : AliGeometry(other),
    fIsInitialized(other.fIsInitialized),
    fInner(other.fInner), 
    fOuter(other.fOuter), 
    fFMD1(other.fFMD1), 
    fFMD2(other.fFMD2), 
    fFMD3(other.fFMD3), 
    fUseFMD1(other.fUseFMD1), 
    fUseFMD2(other.fUseFMD2), 
    fUseFMD3(other.fUseFMD3), 
    fBuilder(other.fBuilder),
    fDetectorOff(other.fDetectorOff),
    fModuleOff(other.fModuleOff),  
    fRingOff(other.fRingOff),
    fSectorOff(other.fSectorOff),
    fActive(other.fActive),
    fDetailed(other.fDetailed),
    fUseAssembly(other.fUseAssembly)
{
  // PROTECTED
  // Copy constructor
}



//____________________________________________________________________
AliFMDGeometry&
AliFMDGeometry::operator=(const AliFMDGeometry& other) 
{
  // PROTECTED
  // Assignment operator 
  fUseFMD1		= other.fUseFMD1; 
  fUseFMD2		= other.fUseFMD2; 
  fUseFMD3		= other.fUseFMD3; 
  fFMD1			= other.fFMD1; 
  fFMD2			= other.fFMD2; 
  fFMD3			= other.fFMD3; 
  fInner		= other.fInner; 
  fOuter		= other.fOuter; 
  fIsInitialized	= other.fIsInitialized;
  return *this;
}

//____________________________________________________________________
void
AliFMDGeometry::Init()
{
  // Initialize the the singleton if not done so already 
  if (fIsInitialized) return;
  fInner->Init();
  fOuter->Init();
  fFMD1->Init();
  fFMD2->Init();
  fFMD3->Init();
}

//____________________________________________________________________
void
AliFMDGeometry::InitTransformations()
{
  // Find all local <-> global transforms 
  if (!gGeoManager) {
    AliError("No TGeoManager defined");
    return;
  }
  AliDebug(0, "Initialising transforms for FMD geometry");
  if (fFMD1) fFMD1->InitTransformations();
  if (fFMD2) fFMD2->InitTransformations();
  if (fFMD3) fFMD3->InitTransformations();
}

//____________________________________________________________________
void
AliFMDGeometry::Build()
{
  // Build the geometry 
  if (!fBuilder) fBuilder = new AliFMDGeometryBuilder(fDetailed);
  fBuilder->SetDetailed(fDetailed);
  fBuilder->UseAssembly(fUseAssembly);
  fBuilder->Exec();
}

//____________________________________________________________________
void
AliFMDGeometry::SetActive(Int_t* active, Int_t n) 
{
  // Set active volumes 
  fActive.Set(n);
  for (Int_t i = 0; i < n; i++) { 
    AliDebug(1, Form("Active vol id # %d: %d", i, active[i]));
    fActive[i] = active[i];
  }
}

//____________________________________________________________________
void
AliFMDGeometry::AddActive(Int_t active)
{
  // Add an active volume 
  Int_t n = fActive.fN;
  fActive.Set(n+1);
  fActive[n] = active;
}

//____________________________________________________________________
Bool_t
AliFMDGeometry::IsActive(Int_t vol) const
{
  // Check if a volume is active 
  for (Int_t i = 0; i < fActive.fN; i++) 
    if (fActive[i] == vol) return kTRUE;
  return kFALSE;
}
  
//____________________________________________________________________
AliFMDDetector*
AliFMDGeometry::GetDetector(Int_t i) const
{
  // Get the ith detector.   i should be one of 1, 2, or 3.  If an
  // invalid value is passed, 0 (NULL) is returned. 
  switch (i) {
  case 1: return fUseFMD1 ? static_cast<AliFMDDetector*>(fFMD1) : 0;
  case 2: return fUseFMD2 ? static_cast<AliFMDDetector*>(fFMD2) : 0;
  case 3: return fUseFMD3 ? static_cast<AliFMDDetector*>(fFMD3) : 0;
  }
  return 0;
}

//____________________________________________________________________
AliFMDRing*
AliFMDGeometry::GetRing(Char_t i) const
{
  // Get the ith ring.  i should be one of 'I' or 'O' (case
  // insensitive).  If an invalid parameter is passed, 0 (NULL) is
  // returned. 
  switch (i) {
  case 'I':
  case 'i': return fInner;
  case 'O':
  case 'o': return fOuter;
  }
  return 0;
}

//____________________________________________________________________
void
AliFMDGeometry::Enable(Int_t i)
{
  // Enable the ith detector.  i should be one of 1, 2, or 3
  switch (i) {
  case 1: fUseFMD1 = kTRUE; break;
  case 2: fUseFMD2 = kTRUE; break;
  case 3: fUseFMD3 = kTRUE; break;
  }
}

//____________________________________________________________________
void
AliFMDGeometry::Disable(Int_t i)
{
  // Disable the ith detector.  i should be one of 1, 2, or 3
  switch (i) {
  case 1: fUseFMD1 = kFALSE; break;
  case 2: fUseFMD2 = kFALSE; break;
  case 3: fUseFMD3 = kFALSE; break;
  }
}

//____________________________________________________________________
void
AliFMDGeometry::Detector2XYZ(UShort_t  detector, 
			     Char_t    ring, 
			     UShort_t  sector, 
			     UShort_t  strip, 
			     Double_t& x, 
			     Double_t& y, 
			     Double_t& z) const
{
  // Translate detector coordinates (detector, ring, sector, strip) to
  // spatial coordinates (x, y, z) in the master reference frame of
  // ALICE. 
  AliFMDDetector* det = GetDetector(detector);
  if (!det) { 
    AliWarning(Form("Unknown detector %d", detector));
    return;
  }
  det->Detector2XYZ(ring, sector, strip, x, y, z);
}

//____________________________________________________________________
Bool_t
AliFMDGeometry::XYZ2Detector(Double_t  x, 
			     Double_t  y, 
			     Double_t  z,
			     UShort_t& detector, 
			     Char_t&   ring, 
			     UShort_t& sector, 
			     UShort_t& strip) const
{
  // Translate spatial coordinates (x,y,z) in the master reference frame of
  // ALICE to the detector coordinates (detector, ring, sector,
  // strip).  Note, that if this method is to be used in
  // reconstruction or the like, then the input z-coordinate should be
  // corrected for the events interactions points z-coordinate, like 
  // geom->XYZ2Detector(x,y,z-ipz,d,r,s,t);
  AliFMDDetector* det = 0;
  detector = 0;
  for (int i = 1; i <= 3; i++) {
    det = GetDetector(i);
    if (!det) continue;
    if (det->XYZ2Detector(x, y, z, ring, sector, strip)) {
      detector = det->GetId();
      return kTRUE;
    }
  }
  return kFALSE;
}


//____________________________________________________________________
void
AliFMDGeometry::GetGlobal(const AliRecPoint* p, 
			  TVector3& pos, 
			  TMatrixF& /* mat */) const 
{
  // Get the global coordinates cooresponding to the reconstructed
  // point p.  The coordiates is returned in the 3-vector pos passed
  // to this member function.  The matrix mat is used for rotations. 
  GetGlobal(p, pos);
}

//____________________________________________________________________
void
AliFMDGeometry::GetGlobal(const AliRecPoint* p, TVector3& pos) const 
{
  // Get the global coordinates cooresponding to the reconstructed
  // point p.  The coordiates is returned in the 3-vector pos passed
  // to this member function. Note, as AliRecPoint only has places for
  // 3 indicies, it is assumed that the ring hit is an inner ring -
  // which obviously needn't be the case. This makes the member
  // function pretty darn useless. 
  // FIXME: Implement this function to work with outer rings too. 
  Double_t x, y, z;
  TVector3 local;
  p->GetLocalPosition(local);
  UShort_t detector = UShort_t(local.X());
  UShort_t sector   = UShort_t(local.Y());
  UShort_t strip    = UShort_t(local.Z());
  Detector2XYZ(detector, 'I', sector, strip, x, y, z);
  pos.SetXYZ(x, y, z);
}

//____________________________________________________________________
Bool_t
AliFMDGeometry::Impact(const TParticle* /* particle */) const 
{ 
  // Return true, if the particle will hit the active detector
  // elements, and false if not.  Should be used for fast
  // simulations.  Note, that the function currently return false
  // always.  
  // FIXME: Implement this function. 
  return kFALSE; 
}

//____________________________________________________________________	
void  
AliFMDGeometry::SetAlignableVolumes() const
{
  // Declare alignable volumes
  for (Int_t d = 1; d <= 3; d++) 
    if (GetDetector(d)) GetDetector(d)->SetAlignableVolumes();
}


//____________________________________________________________________	
void  
AliFMDGeometry::ExtractGeomInfo()
{
  // Check the volume depth of some nodes, get the active volume
  // numbers, and so forth. 
  // 
  // TODO: Here, we should actually also get the parameters of the
  // shapes, like the verticies of the polygon shape that makes up the
  // silicon sensor, the strip pitch, the ring radii, the z-positions,
  // and so on - that is, all the geometric information we need for
  // futher processing, such as simulation, digitization,
  // reconstruction, etc. 
  Int_t detectorDepth = FindNodeDepth("FMD1_1", "ALIC");
  Int_t ringDepth     = FindNodeDepth(Form("FMDI_%d", Int_t('I')), "ALIC");
  Int_t moduleDepth   = FindNodeDepth("FIFV_0", "ALIC");
  Int_t sectorDepth   = FindNodeDepth("FISE_1", "ALIC");
  fActive.Reset(-1);
  AliDebug(1, Form("Geometry depths:\n"
		   "   Sector:     %d\n"
		   "   Module:     %d\n"
		   "   Ring:       %d\n"
		   "   Detector:   %d", 
		   sectorDepth, moduleDepth, ringDepth, detectorDepth));
  if (sectorDepth < 0 && moduleDepth < 0) {
    fDetailed    = kFALSE;
    fSectorOff   = -1;
    fModuleOff   = -1;
    fRingOff     = 0;
    fDetectorOff = (ringDepth - detectorDepth);
    TGeoVolume* actiVol = gGeoManager->GetVolume("FIAC");
    TGeoVolume* actoVol = gGeoManager->GetVolume("FOAC");
    if (actiVol) AddActive(actiVol->GetNumber());
    if (actiVol) AddActive(actoVol->GetNumber());
  }
  else if (sectorDepth < 0) {
    fDetailed    = kFALSE;
    fSectorOff   = -1;
    fModuleOff   = 1;
    fRingOff     = (moduleDepth - ringDepth) + 1;
    fDetectorOff = (moduleDepth - detectorDepth) + 1;
    TGeoVolume* modiVol = gGeoManager->GetVolume("FIMO");
    TGeoVolume* modoVol = gGeoManager->GetVolume("FOMO");
    if (modiVol) AddActive(modiVol->GetNumber());
    if (modoVol) AddActive(modoVol->GetNumber());
  }
  else {
    Int_t stripDepth    = FindNodeDepth("FIST_1", "ALIC");
    fDetailed    = kTRUE;
    fSectorOff   = (stripDepth - sectorDepth);
    fModuleOff   = (moduleDepth >= 0 ? (stripDepth - moduleDepth) : -1);
    fRingOff     = (stripDepth - ringDepth);
    fDetectorOff = (stripDepth - detectorDepth );
    TGeoVolume* striVol = gGeoManager->GetVolume("FIST");
    TGeoVolume* stroVol = gGeoManager->GetVolume("FOST");
    if (striVol) AddActive(striVol->GetNumber());
    if (stroVol) AddActive(stroVol->GetNumber());
  }    
  AliDebug(1, Form("Geometry offsets:\n"
		   "   Sector:     %d\n"
		   "   Module:     %d\n"
		   "   Ring:       %d\n"
		   "   Detector:   %d", 
		   fSectorOff, fModuleOff, fRingOff, fDetectorOff));
}

  
//____________________________________________________________________	
static Int_t 
CheckNodes(TGeoNode* node, const char* name, Int_t& lvl)
{
  // If there's no node here. 
  if (!node) return -1;
  // Check if it this one 
  TString sname(name);
  if (sname == node->GetName()) return lvl;

  // Check if the node is an immediate daugther 
  TObjArray* nodes = node->GetNodes();
  if (!nodes) return -1;
  // Increase the level, and search immediate sub nodes. 
  lvl++;
  TGeoNode*  found = static_cast<TGeoNode*>(nodes->FindObject(name));
  if (found) return lvl;

  // Check the sub node, if any of their sub-nodes match.
  for (Int_t i = 0; i < nodes->GetEntries(); i++) {
    TGeoNode* sub = static_cast<TGeoNode*>(nodes->At(i));
    if (!sub) continue;
    // Recurive check 
    if (CheckNodes(sub, name, lvl) >= 0) return lvl;
  }
  // If not found, decrease the level 
  lvl--;
  return -1;
}
//____________________________________________________________________	
Int_t 
FindNodeDepth(const char* name, const char* volname) 
{
  // Find the depth of a node 
  TGeoVolume* vol  = gGeoManager->GetVolume(volname);
  if (!vol) {
    std::cerr << "No top volume defined" << std::endl;
    return -1;
  }
  TObjArray* nodes = vol->GetNodes();
  if (!nodes) { 
    std::cerr << "No nodes in top volume" << std::endl;
    return -1;
  }
  TIter next(nodes);
  TGeoNode* node = 0;
  Int_t lvl = 0;
  while ((node = static_cast<TGeoNode*>(next()))) 
    if (CheckNodes(node, name, lvl) >= 0) return lvl;
  return -1;
}

//____________________________________________________________________
//
// EOF
//
