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

/*
$Log$
Revision 1.4  1999/09/29 09:24:29  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  This class contains the points for the ALICE event display               //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliPointsClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliPoints.h"
#include "AliRun.h"
#include "AliDetector.h"
#include "TPad.h"
#include "TView.h"
 
ClassImp(AliPoints)

//_____________________________________________________________________________
AliPoints::AliPoints()
{
  //
  // Default constructor
  //
  fDetector = 0;	
  fIndex    = 0;
}

//_____________________________________________________________________________
AliPoints::AliPoints(const AliPoints &pts)
{
  //
  // Copy constructor
  //
  pts.Copy(*this);
}

//_____________________________________________________________________________
AliPoints::AliPoints(Int_t nhits)
  :TPolyMarker3D(nhits)
{
  //
  // Standard constructor
  //
  fDetector = 0;	
  fIndex    = 0;
  ResetBit(kCanDelete);
}
	 
//_____________________________________________________________________________
AliPoints::~AliPoints()
{
  //
  // Default constructor
  //
  fDetector = 0;	
  fIndex    = 0;
}

//_____________________________________________________________________________
void AliPoints::Copy(AliPoints &pts) const
{
  //
  // Copy *this onto pts
  //
  if(this != &pts) {
    ((TPolyMarker3D*)this)->Copy((TPolyMarker3D&)pts);
    pts.fGLList = fGLList;
    pts.fLastPoint = fLastPoint;
    pts.fDetector = fDetector;
    pts.fIndex = fIndex;
  }
}


//_____________________________________________________________________________
Int_t AliPoints::DistancetoPrimitive(Int_t px, Int_t py)
{
  //
  //*-*-*-*-*-*-*Compute distance from point px,py to a 3-D polymarker*-*-*-*-*
  //*-*          =====================================================
  //*-*
  //*-*  Compute the closest distance of approach from point
  //*-*  px,py to each segment
  //*-*  of the polyline.
  //*-*  Returns when the distance found is below DistanceMaximum.
  //*-*  The distance is computed in pixels units.
  //*-*
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  //const Int_t inaxis = 7;
  //Int_t dist = 9999;
  return TPolyMarker3D::DistancetoPrimitive(px,py);
}

//_____________________________________________________________________________
void AliPoints::DumpParticle()
{
  //
  //   Dump particle corresponding to this point
  //
  TParticle *particle = GetParticle();
  if (particle) particle->Dump();
}

//_____________________________________________________________________________
void AliPoints::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
  //
  //*-*-*-*-*-*-*-*-*-*Execute action corresponding to one event*-*-*-*-*-*-*-*
  //*-*                =========================================
  //*-*
  //*-*  This member function must be implemented to realize the action
  //*-*  corresponding to the mouse click on the object in the window
  //*-*
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  gPad->SetCursor(kCross);
  
  if (gPad->GetView())
    gPad->GetView()->ExecuteRotateView(event, px, py);

}

//_____________________________________________________________________________
const Text_t *AliPoints::GetName() const
{
  //
  // Return name of the Geant3 particle corresponding to this point
  //
  TParticle *particle = GetParticle();
  if (!particle) return "Particle";
  return particle->GetName();
}

//_____________________________________________________________________________
Text_t *AliPoints::GetObjectInfo(Int_t, Int_t)
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

//_____________________________________________________________________________
TParticle *AliPoints::GetParticle() const
{
  //
  //   Returns pointer to particle index in AliRun::fParticles
  //
  TClonesArray *particles = gAlice->Particles();
  Int_t nparticles = particles->GetEntriesFast();
  if (fIndex < 0 || fIndex >= nparticles) return 0;
  return (TParticle*)particles->UncheckedAt(fIndex);
}

//_____________________________________________________________________________
void AliPoints::InspectParticle()
{
  //
  //   Inspect particle corresponding to this point
  //
  TParticle *particle = GetParticle();
  if (particle) particle->Inspect();
}

//_____________________________________________________________________________
AliPoints & AliPoints::operator=(const AliPoints &pts)
{
  //
  // Assignment operator
  //
  pts.Copy(*this);
  return (*this);
}

//_____________________________________________________________________________
void AliPoints::Propagate()
{
  //
  //   Set attributes of all detectors to be the attributes of this point
  //
  Int_t ntracks,track;
  TObjArray *points;
  AliPoints *pm;
  //  
  TIter next(gAlice->Detectors());
  AliDetector *detector;
  while((detector = (AliDetector*)next())) {
    if (!detector->IsActive()) continue;
    points = detector->Points();
    if (!points) continue;
    ntracks = points->GetEntriesFast();
    for (track=0;track<ntracks;track++) {
      pm = (AliPoints*)points->UncheckedAt(track);
      if (!pm) continue;
      if (fIndex == pm->GetIndex()) {
	pm->SetMarkerColor(GetMarkerColor());
	pm->SetMarkerSize(GetMarkerSize());
	pm->SetMarkerStyle(GetMarkerStyle());
      }
    }
  }
}
