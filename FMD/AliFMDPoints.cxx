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
/** @file    AliFMDPoints.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Tue Apr 11 12:42:50 2006
    @brief   Specialised class for drawing hits in the FMD. 
    @ingroup FMD_sim    
*/
// Specialised class for drawing hits in the FMD.   The normal
// AliPoints class isn't really suited for the FMD, as it connects the
// dots between hits from the same particle.  However, in the FMD, all
// hits are not identified by a track, so this has little meaning.
// What is more interresting, is the actual size of the hit. 
#include <TMath.h>		// ROOT_TMath
#include <TVector3.h>           // ROOT_TVector3
#include <TMarker3DBox.h>       // ROOT_TMarker3DBox
#include "AliFMDHit.h"          // ALIFMDHIT_H 
#include "AliFMDPoints.h"	// ALIFMDPOINTS_H

//____________________________________________________________________
AliFMDPoints::AliFMDPoints(AliFMDHit* hit, UInt_t color) 
  : AliPoints(1), fMarker(0)
{
  // Constructor 
  // Params: 
  //   hit     Hit to represent 
  //   color   Color of marker
  //
  if (!hit) return;
  Float_t  size  = TMath::Min(TMath::Max(hit->Edep() * .1, .1), 1.);
  TVector3 p(hit->Px(), hit->Py(), hit->Pz());
  fMarker = new TMarker3DBox(hit->X(), hit->Y(), hit->Z(), size, size, size,
			     p.Theta(), p.Phi());
  fMarker->SetLineColor(color);
  fMarker->SetRefObject(this);
  fP[0] = hit->X();
  fP[1] = hit->Y();
  fP[2] = hit->Z();
}

//____________________________________________________________________
AliFMDPoints::AliFMDPoints(const AliFMDPoints& other) 
  : AliPoints(other), fMarker(other.fMarker)
{
  // Copy constructor 
}

//____________________________________________________________________
AliFMDPoints&
AliFMDPoints::operator=(const AliFMDPoints& other) 
{
  // Assignment operator 
  fMarker = other.fMarker;
  return *this;
}


//____________________________________________________________________
AliFMDPoints::~AliFMDPoints() 
{
  // Destructor 
  // if (fMarker) delete  fMarker;
}
//____________________________________________________________________
void 
AliFMDPoints::SetXYZ(Double_t x, Double_t y, Double_t z)
{
  // Set (x,y,z) position of marker 
  // Params
  //   X     X position 
  //   Y     Y position 
  //   Z     Z position 
  if (fMarker) fMarker->SetPosition(x, y, z);
}

//____________________________________________________________________
Int_t
AliFMDPoints::DistancetoPrimitive(Int_t px, Int_t py) 
{
  // Calculate distance from (px,py) to this 
  // Params: 
  //   Px      X-coordinate of marker 
  //   Py      Y-coordinate of marker 
  // Return 
  //   Distance to this 
  return fMarker->DistancetoPrimitive(px, py);
}
//____________________________________________________________________
void 
AliFMDPoints::Draw(Option_t* option) 
{
  // Draw on pad 
  // Params 
  //   option   See TMarker3DBox::Draw
  if (fMarker) fMarker->Draw(option);
}
//____________________________________________________________________
void 
AliFMDPoints::Paint(Option_t* option)
{
  // Draw on pad 
  // Params 
  //   option   See TMarker3DBox::Paint
  if (fMarker) fMarker->Paint(option);
}

//____________________________________________________________________
void 
AliFMDPoints::SetMarkerColor(Color_t colour) 
{
  // Set the marker color
  // Params 
  //   colour   Colour of marker 
  if (fMarker) fMarker->SetLineColor(colour);
}

//____________________________________________________________________
//
// EOF
//
