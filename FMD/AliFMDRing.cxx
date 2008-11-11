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
    fHoneycombThickness(0.),
    fAlThickness(0.),
    fVerticies(0), 
    fSensorVerticies(0),
    fHybridVerticies(0),
    fFeetPositions(0)
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
  SetHoneycombThickness();
  SetAlThickness();
  
  if (fId == 'I' || fId == 'i') {
    SetLowR(4.3);
    SetHighR(17.2);
    SetTheta(36/2);
    SetNStrips(512);
    Double_t base = 0; // 4.1915;
    fFeetPositions.Add(new TVector2( 0.0551687, 8.0534-base));
    fFeetPositions.Add(new TVector2( 2.9993,   12.9457-base));
    fFeetPositions.Add(new TVector2(-2.9062,   12.9508-base));
    
    fHybridVerticies.Add(new TVector2(0.0000,  4.1700));
    fHybridVerticies.Add(new TVector2(1.0574,  4.1700));
    fHybridVerticies.Add(new TVector2(4.6614, 15.2622));
    fHybridVerticies.Add(new TVector2(0.9643, 17.4000));
    fHybridVerticies.Add(new TVector2(0.0000, 17.4000));

    fSensorVerticies.Add(new TVector2(0.0000,  4.1915));
    fSensorVerticies.Add(new TVector2(1.5793,  4.1915));
    fSensorVerticies.Add(new TVector2(5.2293, 15.4251));
    fSensorVerticies.Add(new TVector2(1.9807, 17.3035));
    fSensorVerticies.Add(new TVector2(0.0000, 17.3035));

    fVerticies.Add(new TVector2(0.0000,  4.3000));
    fVerticies.Add(new TVector2(1.3972,  4.3000));
    fVerticies.Add(new TVector2(4.9895, 15.3560));
    fVerticies.Add(new TVector2(1.8007, 17.2000));
    fVerticies.Add(new TVector2(0.0000, 17.2000));
  }
  else if (fId == 'O' || fId == 'o') {
    SetLowR(15.6);
    SetHighR(28.0);
    SetTheta(18/2);
    SetNStrips(256);
    Double_t base = 0; // 14.9104;
    fFeetPositions.Add(new TVector2(-1.72540000, 20.6267-base));
    fFeetPositions.Add(new TVector2( 1.72900000, 20.6267-base));
    fFeetPositions.Add(new TVector2( 0.00177616, 26.6007-base));

    fHybridVerticies.Add(new TVector2(0.0000, 14.9104));
    fHybridVerticies.Add(new TVector2(2.0783, 14.9104));
    fHybridVerticies.Add(new TVector2(3.9202, 26.5395));
    fHybridVerticies.Add(new TVector2(0.6784, 28.2500));
    fHybridVerticies.Add(new TVector2(0.0000, 28.2500));

    fSensorVerticies.Add(new TVector2(0.0000, 15.0104));
    fSensorVerticies.Add(new TVector2(2.5799, 15.0104));
    fSensorVerticies.Add(new TVector2(4.4439, 26.7766));
    fSensorVerticies.Add(new TVector2(1.8350, 28.1500));
    fSensorVerticies.Add(new TVector2(0.0000, 28.1500));

    fVerticies.Add(new TVector2(0.0000, 15.2104));
    fVerticies.Add(new TVector2(2.4091, 15.2104));
    fVerticies.Add(new TVector2(4.2231, 26.6638));
    fVerticies.Add(new TVector2(1.8357, 27.9500));
    fVerticies.Add(new TVector2(0.0000, 27.9500));
  }
}

//____________________________________________________________________
void
AliFMDRing::Init()
{
  // Initialize 
#if 0
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
#endif

  // A's length. Corresponds to distance from nominal beam line to the
  // cornor of the active silicon element. 
  fMinR = GetVertex(1)->Mod(); // GetVertex(5)->Mod();
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
TVector2*
AliFMDRing::GetSensorVertex(Int_t i) const
{
  // Get the i'th vertex of polygon shape
  return static_cast<TVector2*>(fSensorVerticies.At(i));
}

//____________________________________________________________________
TVector2*
AliFMDRing::GetHybridVertex(Int_t i) const
{
  // Get the i'th vertex of polygon shape
  return static_cast<TVector2*>(fHybridVerticies.At(i));
}

//____________________________________________________________________
TVector2*
AliFMDRing::GetFootPosition(Int_t i) const
{
  // Get the i'th vertex of polygon shape
  return static_cast<TVector2*>(fFeetPositions.At(i));
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
Double_t
AliFMDRing::GetModuleDepth() const
{
  return (GetSiThickness() 
	  + GetSpacing() 
	  + GetPrintboardThickness()
	  + GetCopperThickness()
	  + GetChipThickness()
	  + GetLegLength());
  
}

//____________________________________________________________________
Double_t
AliFMDRing::GetFullDepth() const
{
  return (GetModuleDepth() 
	  + GetModuleSpacing()
	  + GetHoneycombThickness()
	  + GetFMDDPrintboardThickness()
	  + GetFMDDCopperThickness()
	  + GetFMDDChipThickness() 
	  + 0.5);
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
