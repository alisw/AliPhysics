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
/** @file    AliFMDRing.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:47:43 2006
    @brief   FMD ring geometry parameters 
*/
//__________________________________________________________________
//
// Utility class to help implement collection of FMD modules into
// rings.  This is used by AliFMDDetector and AliFMDGeometry.
// The AliFMDGeometry object owns the AliFMDRing objects, and the
// AliFMDDetector objects reference these.  That is, the AliFMDRing
// objects are share amoung the AliFMDDetector objects.
//
// Latest changes by Christian Holm Christensen
//

#include <TMath.h>		// ROOT_TMath
#include <TVector2.h>		// ROOT_TVector2

// #include <AliLog.h>		// ALILOG_H
#include "AliFMDRing.h"		// ALIFMDRING_H

//====================================================================
ClassImp(AliFMDRing)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDRing::AliFMDRing(Char_t id) 
  : TNamed(Form("FMD%c", id), "Forward multiplicity ring"), 
    fId(id), 
    fBondingWidth(0),
    fWaferRadius(0),
    fSiThickness(0),
    fLowR(0),
    fHighR(0),
    fMinR(0),
    fMaxR(0),
    fTheta(0),
    fNStrips(0),
    fRingDepth(0),
    fLegRadius(0),
    fLegLength(0),
    fLegOffset(0),
    fModuleSpacing(0),
    fPrintboardThickness(0),
    fCopperThickness(0),
    fChipThickness(0),
    fSpacing(0),
    fVerticies(0)
{
  // CTOR
  SetBondingWidth();
  SetWaferRadius();
  SetSiThickness();
  SetLegRadius();
  SetLegLength();
  SetLegOffset();
  SetModuleSpacing();
  SetPrintboardThickness();
  SetCopperThickness();
  SetChipThickness();
  SetSpacing();
  
  if (fId == 'I' || fId == 'i') {
    SetLowR(4.3);
    SetHighR(17.2);
    SetTheta(36/2);
    SetNStrips(512);
  }
  else if (fId == 'O' || fId == 'o') {
    SetLowR(15.6);
    SetHighR(28.0);
    SetTheta(18/2);
    SetNStrips(256);
  }
}

//____________________________________________________________________
void
AliFMDRing::Init()
{
  // Initialize 
  Double_t tanTheta  = TMath::Tan(fTheta * TMath::Pi() / 180.);
  Double_t tanTheta2 = TMath::Power(tanTheta,2);
  Double_t r2        = TMath::Power(fWaferRadius,2);
  Double_t yA        = tanTheta * fLowR;
  Double_t lr2       = TMath::Power(fLowR, 2);
  Double_t hr2       = TMath::Power(fHighR,2);
  Double_t xD        = fLowR + TMath::Sqrt(r2 - tanTheta2 * lr2);
  Double_t xD2       = TMath::Power(xD,2);
  Double_t yB        = TMath::Sqrt(r2 - hr2 + 2 * fHighR * xD - xD2);
  Double_t xC        = ((xD + TMath::Sqrt(-tanTheta2 * xD2 + r2
					  + r2 * tanTheta2)) 
			/ (1 + tanTheta2));
  Double_t yC        = tanTheta * xC;
  
  fVerticies.Expand(6);
  fVerticies.AddAt(new TVector2(fLowR,  -yA), 0);
  fVerticies.AddAt(new TVector2(xC,     -yC), 1);
  fVerticies.AddAt(new TVector2(fHighR, -yB), 2);
  fVerticies.AddAt(new TVector2(fHighR,  yB), 3);
  fVerticies.AddAt(new TVector2(xC,      yC), 4);
  fVerticies.AddAt(new TVector2(fLowR,   yA), 5);  

  // A's length. Corresponds to distance from nominal beam line to the
  // cornor of the active silicon element. 
  fMinR = GetVertex(5)->Mod();
  // A's length. Corresponds to distance from nominal beam line to the
  // cornor of the active silicon element. 
  fMaxR = fHighR;

  fRingDepth = (fSiThickness + fPrintboardThickness 
		+ fCopperThickness + fChipThickness 
		+ fLegLength + fModuleSpacing + fSpacing);
}

//____________________________________________________________________
TVector2*
AliFMDRing::GetVertex(Int_t i) const
{
  // Get the i'th vertex of polygon shape
  return static_cast<TVector2*>(fVerticies.At(i));
}

//____________________________________________________________________
Double_t
AliFMDRing::GetStripRadius(UShort_t strip) const
{
  // Return the nominal strip radius 
  Double_t rmax     = GetMaxR();
  Double_t stripoff = GetMinR();
  Double_t dstrip   = (rmax - stripoff) / GetNStrips();
  return (strip + .5) * dstrip + stripoff; // fLowR
}
 
//____________________________________________________________________
void
AliFMDRing::Detector2XYZ(UShort_t sector,
			 UShort_t strip, 
			 Double_t& x, 
			 Double_t& y, 
			 Double_t& z) const
{
  // Translate detector coordinates (this,sector,strip) to global
  // coordinates (x,y,z)
  if (sector >= GetNSectors()) {
    Error("Detector2XYZ", "Invalid sector number %d (>=%d) in ring %c", 
	  sector, GetNSectors(), fId);
    return;
  }
  if (strip >= GetNStrips()) {
    Error("Detector2XYZ", "Invalid strip number %d (>=%d)", 
	  strip, GetNStrips(), fId);
    return;
  }
  Double_t phi = Float_t(sector + .5) / GetNSectors() * 2 * TMath::Pi();
  Double_t r   = Float_t(strip + .5) / GetNStrips() * (fHighR - fLowR) + fLowR;
  x = r * TMath::Cos(phi);
  y = r * TMath::Sin(phi);
  if (((sector / 2) % 2) == 1) 
    z += TMath::Sign(fModuleSpacing, z);
}

//____________________________________________________________________
Bool_t
AliFMDRing::XYZ2Detector(Double_t  x, 
			 Double_t  y, 
			 Double_t  z,
			 UShort_t& sector,
			 UShort_t& strip) const
{
  // Translate global coordinates (x,y,z) to detector coordinates
  // (this,sector,strip)
  sector = strip = 0;
  Double_t r = TMath::Sqrt(x * x + y * y);
  Int_t str = Int_t((r - fMinR) / GetPitch());
  if (str < 0 || str >= GetNStrips()) return kFALSE;

  Double_t phi = TMath::ATan2(y, x) * 180. / TMath::Pi();
  if (phi < 0) phi = 360. + phi;
  Int_t sec = Int_t(phi / fTheta);
  if (sec < 0 || sec >= GetNSectors()) return kFALSE;
  if ((sec / 2) % 2 == 1) {
    if (TMath::Abs(z - TMath::Sign(fModuleSpacing, z)) >= 0.01)
      return kFALSE;
  }
  else if (TMath::Abs(z) >= 0.01) return kFALSE;

  strip  = str;
  sector = sec;
  return kTRUE;
}


//
// EOF
//
