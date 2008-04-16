// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

/***********************************************************************
*  This code defines the reconstructed v0 visualized with EVE
*
* Ludovic Gaudichet (gaudichet@to.infn.it)
************************************************************************/

#include "AliEveV0.h"

#include <TEveTrack.h>
#include <TEveTrackPropagator.h>
#include <TEveManager.h>

#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TColor.h>

#include <vector>


/***********************************************************************
*
*  AliEveV0 class
*
************************************************************************/

ClassImp(AliEveV0)

AliEveV0::AliEveV0() :
  TEvePointSet(),

  fRecBirthV(),
  fRecDecayV(),
  fRecDecayP(),
  fNegTrack(0),
  fPosTrack(0),
  fRnrStyle(0),
  fPolyLineV0(),
  fESDIndex(-1),
  fDaughterDCA(999),
  fChi2V0(-1)
{}


AliEveV0::AliEveV0(TEveRecTrack* tNeg, TEveRecTrack* tPos,
		   TEveRecV0* v0, TEveTrackPropagator* rs) :
  TEvePointSet(),

  fRecBirthV(v0->fV0Birth),
  fRecDecayV(v0->fVCa),
  fRecDecayP(v0->fPNeg + v0->fPPos),

  fNegTrack(new TEveTrack(tNeg, rs)),
  fPosTrack(new TEveTrack(tPos, rs)),

  fRnrStyle(rs),
  fPolyLineV0(),
  fESDIndex(-1),
  fDaughterDCA(999),
  fChi2V0(-1)
{
  fPolyLineV0.SetLineColor(fMarkerColor);
 
  fPosTrack->SetLineColor(2);  // red
  fNegTrack->SetLineColor(7);  // light blue

  fMainColorPtr = &fMarkerColor;
  fMarkerStyle = 20;
  fMarkerColor = 5;
  fMarkerSize  = 0.3;

  AddElement(fNegTrack);
  AddElement(fPosTrack);
}

AliEveV0::~AliEveV0()
{}


void AliEveV0::Reset(TPolyLine3D* polyLine)
{
  //polyLine->SetPolyLine(n_points);
  polyLine->SetPolyLine(0);
}

//______________________________________________________________________________
void AliEveV0::MakeV0path()
{
  fPolyLineV0.SetPoint(0, fRecBirthV.fX, fRecBirthV.fY, fRecBirthV.fZ);
  fPolyLineV0.SetPoint(1, fRecDecayV.fX, fRecDecayV.fY, fRecDecayV.fZ);
}


//______________________________________________________________________________
void AliEveV0::MakeV0()
{
  SetPoint(0, fRecDecayV.fX, fRecDecayV.fY, fRecDecayV.fZ);

  fNegTrack->MakeTrack();
  fPosTrack->MakeTrack();
  MakeV0path();
}


/***********************************************************************
*
*  AliEveV0List class
*
************************************************************************/

ClassImp(AliEveV0List)

//______________________________________________________________________________
AliEveV0List::AliEveV0List() :
  TEveElementList(),
  fTitle(),
  fRnrStyle(0),
  fRnrDaughters(kTRUE),
  fRnrV0vtx(kTRUE),
  fRnrV0path(kTRUE),
  fNegColor(0),
  fPosColor(0)
{
  fChildClass = AliEveV0::Class(); // override member from base TEveElementList
}

//______________________________________________________________________________
AliEveV0List::AliEveV0List(TEveTrackPropagator* rs) :
  TEveElementList(),
  fTitle(),
  fRnrStyle(rs),
  fRnrDaughters(kTRUE),
  fRnrV0vtx(kTRUE),
  fRnrV0path(kTRUE),
  fNegColor(0),
  fPosColor(0)
{
  fChildClass = AliEveV0::Class(); // override member from base TEveElementList

  Init();
}

//______________________________________________________________________________
AliEveV0List::AliEveV0List(const Text_t* name, TEveTrackPropagator* rs) :
  TEveElementList(),
  fTitle(),
  fRnrStyle(rs),
  fRnrDaughters(kTRUE),
  fRnrV0vtx(kTRUE),
  fRnrV0path(kTRUE),
  fNegColor(0),
  fPosColor(0)
{
  fChildClass = AliEveV0::Class(); // override member from base TEveElementList

  Init();
  SetName(name);
}

//______________________________________________________________________________
void AliEveV0List::Init()
{
  if (fRnrStyle== 0) fRnrStyle = new TEveTrackPropagator;

}

//______________________________________________________________________________
AliEveV0List::~AliEveV0List()
{

}

//______________________________________________________________________________
void AliEveV0List::Paint(Option_t* option)
{
  if(fRnrSelf) {

    if(fRnrV0vtx) {
      for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
	if((*i)->GetRnrSelf()) {
	  ((AliEveV0*)(*i))->Paint(option);
	}
      }
    }

    if(fRnrDaughters) {
      for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
	if((*i)->GetRnrSelf()) {
	  ((AliEveV0*)(*i))->PaintDaughters(option);
	}
      }
    }

    if(fRnrV0path) {
      for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
	if((*i)->GetRnrSelf()) {
	  ((AliEveV0*)(*i))->PaintPath(option);
	}
      }
    }
  }
}


//______________________________________________________________________________

void AliEveV0List::SetRnrV0vtx(Bool_t rnr)
{
  fRnrV0vtx = rnr;
  gEve->Redraw3D();
}

void AliEveV0List::SetRnrV0path(Bool_t rnr)
{
  fRnrV0path = rnr;
  gEve->Redraw3D();
}

void AliEveV0List::SetRnrDaughters(Bool_t rnr)
{
  fRnrDaughters = rnr;
  gEve->Redraw3D();
}


//______________________________________________________________________________

void AliEveV0List::MakeV0s()
{
  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    ((AliEveV0*)(*i))->MakeV0();
  }
  gEve->Redraw3D();
}


void AliEveV0List::MakeMarkers()
{
  gEve->Redraw3D();
}
