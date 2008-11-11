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
#include <TVector3.h>

//====================================================================
ClassImp(AliFMD3)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMD3::AliFMD3(AliFMDRing* inner, AliFMDRing* outer) 
  : AliFMDDetector(3, inner, outer),
    // fNoseZ(16.54667)
    fNoseZ(18.13 - inner->GetModuleDepth()-inner->GetModuleSpacing()/2),    // From drawing
    fFlangeDepth(0),
    fFlangeHighR(49.25), // From drawing 
    fFlangeLength(0),
    fFlangeWidth(6),     // From drawing 
    fFiducialRadius(.25),
    fConeInnerAngle(0),
    fConeOuterAngle(0),
    fHoleOffset(6.9),    // From drawing
    fHoleDepth(2),       // What's needed
    fHoleLength(31.2),   // From drawing
    fHoleLowWidth(3), // 4),    // What's needed
    fHoleHighWidth(18.5), // 17.5), // 18),  // What's needed
    fBoltLength(1),      // Guessed
    fBoltRadius(0.15),   // Estimate
    fConeRadii(6),    
    fFiducialHoles(4)
{
  // Constructor. 
  Double_t off = -0.39615-0.10185; // -0.25;
  if (off != 0) 
    AliWarning(Form("Z position of FMD3 rings may be off by %fcm!", off));

  SetInnerZ(-62.8+off);             // By design
  SetOuterZ(-75.2+off);             // By design

  SetInnerHoneyLowR(4.18207);   // From drawing
  SetInnerHoneyHighR(19.74922); // From drawing
  SetOuterHoneyLowR(13.4776);   // From drawing
  SetOuterHoneyHighR(31.01964); // From drawing
  
  // These are from the drawings
  fConeRadii.Add(new TVector3( 0,       5.55,  6.25));
  fConeRadii.Add(new TVector3( 2.35,    5.55,  6.25));
  fConeRadii.Add(new TVector3( 2.9935,  5.55,  6.88479));
  fConeRadii.Add(new TVector3(28.9435, 31.50, 32.75850));
  fConeRadii.Add(new TVector3(29.5,    31.50, 33.4));
  fConeRadii.Add(new TVector3(30.9,    31.50, 33.4));

  // These are from the drawings
  fFiducialHoles.Add(new TVector2(29.666, 32.495));
  fFiducialHoles.Add(new TVector2(31.082, 33.910));
  fFiducialHoles.Add(new TVector2(32.674, 35.503));
  fFiducialHoles.Add(new TVector2(33.403, 34.818));
}

//____________________________________________________________________
void
AliFMD3::Init() 
{
  // Initialize 
  AliFMDDetector::Init();
  // TVector3& v0 = *(static_cast<TVector3*>(fConeRadii.At(0)));
  TVector3& v1 = *(static_cast<TVector3*>(fConeRadii.At(1)));
  TVector3& v2 = *(static_cast<TVector3*>(fConeRadii.At(2)));
  TVector3& v3 = *(static_cast<TVector3*>(fConeRadii.At(3)));
  TVector3& v4 = *(static_cast<TVector3*>(fConeRadii.At(4)));
  TVector3& v5 = *(static_cast<TVector3*>(fConeRadii.At(5)));
  
  fFlangeDepth     = v5.X() - v4.X();
  fFlangeLength    = fFlangeHighR - v5.Y();
  
  fConeInnerAngle  = TMath::ATan2(v4.Z()-v1.Z(), v4.X()-v1.X());
  fConeOuterAngle  = TMath::ATan2(v3.Y()-v2.Y(), v3.X()-v2.X());
  
#if 0
  Double_t    hz1  = -fHoleOffset+fInnerZ+fNoseZ;
  fHoleLength      = TMath::Sqrt(TMath::Power(v4.Z()-ConeR(hz1),2) + 
				 TMath::Power(v4.X()-fHoleOffset,2));  
#endif
}

//____________________________________________________________________
Double_t
AliFMD3::ConeR(Double_t z, Option_t* opt) const
{
  // Calculate the cone radius at Z
  // TVector3& v0 = *(static_cast<TVector3*>(fConeRadii.At(0)));
  TVector3& v1 = *(static_cast<TVector3*>(fConeRadii.At(1)));
  TVector3& v2 = *(static_cast<TVector3*>(fConeRadii.At(2)));
  TVector3& v3 = *(static_cast<TVector3*>(fConeRadii.At(3)));
  TVector3& v4 = *(static_cast<TVector3*>(fConeRadii.At(4)));
  TVector3& v5 = *(static_cast<TVector3*>(fConeRadii.At(5)));

  if (z > fInnerZ + fNoseZ) {
    AliWarning(Form("z=%lf is before start of cone %lf", z, fInnerZ + fNoseZ));
    return -1;
  }
  if (z < fInnerZ + fNoseZ - v5.Z()) {
    AliWarning(Form("z=%lf is after end of cone %lf", z, 
		    fInnerZ + fNoseZ - v5.Z()));
    return -1;
  }
  Double_t rz    = -(z-fInnerZ-fNoseZ);
  Bool_t   inner = opt[0] == 'I' || opt[1] == 'i';
  if (inner  && rz <= v2.X()) return v2.Y();
  if (!inner && rz <= v1.X()) return v1.Z();
  if (inner  && rz >  v3.X()) return v3.Y();
  if (!inner && rz >  v4.X()) return v4.Z();
  
  rz             -= (inner ? v2.X() : v1.X());
  Double_t sr    =  (inner ? v2.Y() : v1.Z());
  Double_t ang   =  (inner ? fConeInnerAngle : fConeOuterAngle);
  return sr + rz * TMath::Tan(ang);
}


//____________________________________________________________________
//
// EOF
//
