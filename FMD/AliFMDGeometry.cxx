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
#include <TMatrix.h>            // ROOT_TMatrix
#include <TParticle.h>          // ROOT_TParticle
#include <Riostream.h>

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
  if (!fgInstance) fgInstance = new AliFMDGeometry;
  return fgInstance;
}

//____________________________________________________________________
AliFMDGeometry::AliFMDGeometry() 
  : AliGeometry("FMD", "Forward multiplicity")
{
  fUseFMD1 = kTRUE;
  fUseFMD2 = kTRUE;
  fUseFMD3 = kTRUE;  
  fInner = new AliFMDRing('I');
  fOuter = new AliFMDRing('O');
  fFMD1  = new AliFMD1(fInner);
  fFMD2  = new AliFMD2(fInner, fOuter);
  fFMD3  = new AliFMD3(fInner, fOuter);
  fIsInitialized = kFALSE;
}

//____________________________________________________________________
void
AliFMDGeometry::Init()
{
  if (fIsInitialized) return;
  fInner->Init();
  fOuter->Init();
  fFMD1->Init();
  fFMD2->Init();
  fFMD3->Init();
}

//____________________________________________________________________
AliFMDDetector*
AliFMDGeometry::GetDetector(Int_t i) const
{
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
  AliFMDDetector* det = GetDetector(detector);
  if (!det) return;
  det->Detector2XYZ(ring, sector, strip, x, y, z);
}


//____________________________________________________________________
void
AliFMDGeometry::GetGlobal(const AliRecPoint* p, 
			  TVector3& pos, 
			  TMatrix& /* mat */) const 
{
  GetGlobal(p, pos);
}

//____________________________________________________________________
void
AliFMDGeometry::GetGlobal(const AliRecPoint* p, TVector3& pos) const 
{
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
  return kFALSE; 
}

//____________________________________________________________________
//
// EOF
//
