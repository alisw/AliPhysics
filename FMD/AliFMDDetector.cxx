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

//____________________________________________________________________
//
// Utility class to help implement the FMD geometry.  This provides
// the interface for the concrete geometry implementations of the FMD
// sub-detectors. 
//
// The AliFMDGeometry object owns the AliFMDDetector objects
//
// Latest changes by Christian Holm Christensen
//
#include "AliFMDDetector.h"	// ALIFMDSUBDETECTOR_H
#include "AliFMDRing.h"		// ALIFMDRING_H
#include <AliLog.h>		// ALILOG_H

//====================================================================
ClassImp(AliFMDDetector)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDDetector::AliFMDDetector(Int_t id, AliFMDRing* inner, AliFMDRing* outer) 
  : TNamed(Form("FMD%d", id), "Forward multiplicity ring"), 
    fId(id), 
    fInner(inner),
    fOuter(outer)
{
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
void
AliFMDDetector::Init()
{
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
AliFMDRing*
AliFMDDetector::GetRing(Char_t id) const
{
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
  switch (id) {
  case 'i':
  case 'I': return GetInnerZ();
  case 'o':
  case 'O': return GetOuterZ();
  }
  return 0;
}
//____________________________________________________________________
void
AliFMDDetector::Detector2XYZ(Char_t ring, 
			     UShort_t sector,
			     UShort_t strip, 
			     Double_t& x, 
			     Double_t& y, 
			     Double_t& z) const
{
  AliFMDRing* r = GetRing(ring);
  if (!r) return;
  z = GetRingZ(ring);
  r->Detector2XYZ(sector, strip, x, y, z);
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
