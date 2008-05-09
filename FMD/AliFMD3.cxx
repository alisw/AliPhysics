/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved. *
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
/** @file    AliFMD3.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:26:12 2006
    @brief   Concrete implementation of AliFMDDetector for FMD3
*/
//____________________________________________________________________
//                                                                          
// Concrete implementation of AliFMDDetector 
//
// This implements the geometry for FMD3.
// This has 2 rings.
// The support of the FMD3 is a carbon-fibre cone, attached to the ITS
// support via flanges.  The cone also supports the beam-pipe.
// The support is a special cone of carbon-fibre made by a Danish
// Yacht company.
//

#include <TMath.h>		// ROOT_TMath

#include "AliFMD3.h"		// ALIFMD3_H 
#include "AliFMDDebug.h"		// ALIFMDDEBUG_H ALILOG_H
#include "AliFMDRing.h"		// ALIFMDRING_H 

//====================================================================
ClassImp(AliFMD3)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMD3::AliFMD3(AliFMDRing* inner, AliFMDRing* outer) 
  : AliFMDDetector(3, inner, outer),
    fNoseZ(0),
    fNoseLowR(0),
    fNoseHighR(0),
    fNoseLength(0),
    fBackLowR(0),
    fBackHighR(0),
    fBackLength(0),
    fBeamThickness(0),
    fBeamWidth(0),
    fConeLength(0),
    fFlangeR(0),
    fZ(0),
    fAlpha(0), 
    fNBeam(0), 
    fNFlange(0)
{
  // Constructor. 
  SetInnerZ(-62.8);
  SetOuterZ(-75.2);
  SetNoseZ();
  SetNoseLowR();
  SetNoseHighR();
  SetNoseLength();
  SetBackLowR();
  SetBackHighR();
  SetBackLength();
  SetBeamThickness();
  SetBeamWidth();
  SetConeLength();
  SetFlangeR();
  SetNBeam();
  SetNFlange();
}

//____________________________________________________________________
void
AliFMD3::Init() 
{
  // Initialize 
  AliFMDDetector::Init();
  SetInnerHoneyHighR(GetOuterHoneyHighR());
  Double_t zdist   = fConeLength;
  Double_t tdist   = fBackHighR - fNoseHighR;
  Double_t innerZh = (fInnerZ - fInner->GetRingDepth());
  Double_t outerZh = (fOuterZ - fOuter->GetRingDepth() - 
		      fOuter->GetHoneycombThickness());
  Double_t minZ    = TMath::Min(fNoseZ - fConeLength, outerZh);
  fAlpha           = tdist / zdist;
  fZ               = fNoseZ + (minZ - fNoseZ) / 2;
  fInnerHoneyHighR = ConeR(innerZh,"I");
  fOuterHoneyHighR = GetBackLowR();
}

//____________________________________________________________________
Double_t
AliFMD3::ConeR(Double_t z, Option_t* opt) const
{
  // Calculate the cone radius at Z
  if (fAlpha < 0) {
    AliWarning(Form("alpha not set: %lf", fAlpha));
    return -1;
  }
  if (z > fNoseZ) {
    AliWarning(Form("z=%lf is before start of cone %lf", z, fNoseZ));
    return -1;
  }
  if (z < fOuterZ - fOuter->GetFullDepth()) {
    AliWarning(Form("z=%lf is after end of cone %lf", z, 
		    fOuterZ - fOuter->GetFullDepth()));
    return -1;
  }
  Bool_t   inner = opt[0] == 'I' || opt[1] == 'i';
  Double_t off1  = (inner ? fNoseLowR : fNoseHighR);
  Double_t off2  = (inner ? fBackLowR : fBackHighR);
  Double_t off3  = (inner ? 0         : fBeamThickness/fAlpha);
  if (z > fNoseZ - fNoseLength)               return off1;
  if (z < fNoseZ - fConeLength + fBackLength) return off2;
  return (off1 + off3 + fAlpha * TMath::Abs(z - fNoseZ + fNoseLength));
}


//____________________________________________________________________
//
// EOF
//
