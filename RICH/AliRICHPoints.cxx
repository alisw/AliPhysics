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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  This class contains the points for the ALICE event display               //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliRICHPointsClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "AliRICHPoints.h"
#include "AliRICHDisplay.h"
#include "AliRICHChamber.h"
#include "AliRICH.h"
#include <TPad.h>
#include <TTree.h>
#include <TView.h>
#include <TMath.h>
#include <TPolyMarker3D.h>
#include <TMarker3DBox.h>
#include <TParticle.h>
#include <AliRun.h>
#include <AliMC.h>
#include <TRotMatrix.h>

const Int_t kMaxNipx=400, kMaxNipy=800;
 
ClassImp(AliRICHPoints)

//_____________________________________________________________________________
AliRICHPoints::AliRICHPoints()
{
  //
  // Default constructor
  //
  fHitIndex = 0;
  fTrackIndex = 0;
  fDigitIndex = 0;
  fMarker[0] = fMarker[1] = fMarker[2]=0;
}

//_____________________________________________________________________________
AliRICHPoints::AliRICHPoints(Int_t npoints)
  :AliPoints(npoints)
{
  //
  // Standard constructor
  //
  fHitIndex = 0;
  fTrackIndex = 0;
  fDigitIndex = 0;
  fMarker[0] = fMarker[1] = fMarker[2]=0;
}
	 
//_____________________________________________________________________________
AliRICHPoints::~AliRICHPoints()
{
  //
  // Default destructor
  //
  fHitIndex = 0;
  fTrackIndex = 0;
  fDigitIndex = 0;
}

//_____________________________________________________________________________
void AliRICHPoints::DumpHit()
{
  //
  //   Dump hit corresponding to this point
  //
  AliRICHhit *hit = GetHit();
  if (hit) hit->Dump();
}

//_____________________________________________________________________________
void AliRICHPoints::DumpDigit()
{
  //
  //   Dump digit corresponding to this point
  //
  AliRICHdigit *digit = GetDigit();
  if (digit) digit->Dump();
}

//_____________________________________________________________________________
void AliRICHPoints::InspectHit()
{
  //
  //   Inspect hit corresponding to this point
  //
  AliRICHhit *hit = GetHit();
  if (hit) hit->Inspect();
}

//_____________________________________________________________________________
void AliRICHPoints::InspectDigit()
{
  //
  //   Inspect digit corresponding to this point
  //
  AliRICHdigit *digit = GetDigit();
  if (digit) digit->Inspect();
}

//_____________________________________________________________________________
Int_t AliRICHPoints::GetTrackIndex()
{
  //
  //   Dump digit corresponding to this point
  //
  printf("GetTrackIndex - fTrackIndex %d \n",fTrackIndex);
  this->Inspect();
  return fTrackIndex;
}
//_____________________________________________________________________________
TParticle *AliRICHPoints::GetParticle() const
{
  //
  //   Returns pointer to particle index in AliRun::fParticles
  //
  if (fIndex < 0 || fIndex >= gAlice->GetMCApp()->GetNtrack()) return 0;
  return gAlice->GetMCApp()->Particle(fIndex);
}

//_____________________________________________________________________________
AliRICHhit *AliRICHPoints::GetHit() const
{
  //
  //   Returns pointer to hit index in AliRun::fParticles
  //
  AliRICH *pRICH  = (AliRICH*)gAlice->GetDetector("RICH");
  pRICH->TreeH()->GetEvent(fTrackIndex);
  TClonesArray *pRICHhits  = pRICH->Hits();
  Int_t nhits = pRICHhits->GetEntriesFast();
  if (fHitIndex < 0 || fHitIndex >= nhits) return 0;
  return (AliRICHhit*)pRICHhits->UncheckedAt(fHitIndex);
}

//_____________________________________________________________________________
AliRICHdigit *AliRICHPoints::GetDigit() const
{
  //
  //   Returns pointer to digit index in AliRun::fParticles
  //

  AliRICHDisplay *display=(AliRICHDisplay*)gAlice->Display();
  Int_t chamber=display->GetChamber();
  Int_t cathode=display->GetCathode();
   
  AliRICH *pRICH  = (AliRICH*)gAlice->GetDetector("RICH");
  TClonesArray *pRICHdigits  = pRICH->Digits(chamber);
  gAlice->TreeD()->GetEvent(cathode);
  Int_t ndigits = pRICHdigits->GetEntriesFast();
  if (fDigitIndex < 0 || fDigitIndex >= ndigits) return 0;
  return (AliRICHdigit*)pRICHdigits->UncheckedAt(fDigitIndex);
}
//----------------------------------------------------------------------------
void AliRICHPoints::ShowRing(Int_t highlight)
{

//
// Highlights all pads generated by the same mother particle

  highlight++;   
  

}

//_____________________________________________________________________________
const Text_t *AliRICHPoints::GetName() const
{
  //
  // Return name of the Geant3 particle corresponding to this point
  //
  TParticle *particle = GetParticle();
  if (!particle) return "Particle";
  return particle->GetName();
}

//_____________________________________________________________________________
Text_t *AliRICHPoints::GetObjectInfo(Int_t, Int_t)
{
  //
  //   Redefines TObject::GetObjectInfo.
  //   Displays the info (particle,etc
  //   corresponding to cursor position px,py
  //
  static char info[64];
  sprintf(info,"%s %d",GetName(),fIndex);
  return info;
}



