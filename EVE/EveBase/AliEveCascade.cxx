// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveCascade.h"

#include <TEveTrack.h>
#include <TEveTrackPropagator.h>
#include <TEveManager.h>

#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TColor.h>

#include <vector>


/***********************************************************************
*
*  AliEveCascade class
*
************************************************************************/

ClassImp(AliEveCascade)

//______________________________________________________________________________
AliEveCascade::AliEveCascade() :
  TEvePointSet(),

  fRecBirthV(),
  fRecDecayV(),
  fRecDecayP(),
  fRecDecayV0(),
  fBacTrack(0),
  fRnrStyle(0),
  fPointingCurve(0),
  fV0Path(0),
  fESDIndex(-1),
  fDaughterDCA(999),
  fChi2Cascade(-1)
{
  // Default constructor.

  // Override from TEveElement.
  fPickable = kTRUE;
  fMainColorPtr = &fMarkerColor;
}

//______________________________________________________________________________
AliEveCascade::AliEveCascade(TEveRecTrack* tBac, TEveRecV0* v0, TEveRecCascade* cascade, TEveTrackPropagator* rs) :
  TEvePointSet(),

  fRecBirthV(cascade->fCascadeBirth),
  fRecDecayV(cascade->fCascadeVCa),
  fRecDecayP(cascade->fPBac + v0->fPNeg + v0->fPPos),
  fRecDecayV0(v0->fVCa),
  fBacTrack(new TEveTrack(tBac, rs)),

  fRnrStyle(rs),
  fPointingCurve(new TEveLine("Pointing Curve")),
  fV0Path(new TEveLine("V0 Path")),
  fESDIndex(-1),
  fDaughterDCA(999),
  fChi2Cascade(-1)
{
  // Constructor with full Cascade specification.

  // Override from TEveElement.
  fPickable = kTRUE;
  fMainColorPtr = &fMarkerColor;

  fMarkerStyle = 2;
  fMarkerColor = kViolet;
  fMarkerSize  = 1;

  fPointingCurve->SetLineColor(fMarkerColor);
  fPointingCurve->SetLineWidth(2);
  fPointingCurve->IncDenyDestroy();
  AddElement(fPointingCurve);

  fV0Path->SetLineColor(fMarkerColor);
  fV0Path->SetLineStyle(3);
  fV0Path->SetLineWidth(2);
  fV0Path->IncDenyDestroy();
  AddElement(fV0Path);

  fBacTrack->SetLineColor(6);
  fBacTrack->SetStdTitle();

  fBacTrack->IncDenyDestroy();
  AddElement(fBacTrack);
}

//______________________________________________________________________________
AliEveCascade::~AliEveCascade()
{
  // Destructor. Dereferences bachelor track and pointing-line objects.

  fBacTrack->DecDenyDestroy();
  fPointingCurve->DecDenyDestroy();
  fV0Path->DecDenyDestroy();
}

//______________________________________________________________________________
void AliEveCascade::MakeCascade()
{
  // Set all dependant components for drawing.

  SetPoint(0, fRecDecayV.fX, fRecDecayV.fY, fRecDecayV.fZ);

  fBacTrack->MakeTrack();

  fPointingCurve->SetPoint(0, fRecBirthV.fX, fRecBirthV.fY, fRecBirthV.fZ);
  fPointingCurve->SetPoint(1, fRecDecayV.fX, fRecDecayV.fY, fRecDecayV.fZ);

  fV0Path->SetPoint(0, fRecDecayV.fX, fRecDecayV.fY, fRecDecayV.fZ);
  fV0Path->SetPoint(1, fRecDecayV0.fX, fRecDecayV0.fY, fRecDecayV0.fZ);
}


/***********************************************************************
*
*  AliEveCascadeList class
*
************************************************************************/

ClassImp(AliEveCascadeList)

//______________________________________________________________________________
AliEveCascadeList::AliEveCascadeList() :
  TEveElementList(),
  fTitle(),
  fRnrStyle(0),
  fRnrDaughters(kTRUE),
  fRnrCascadevtx(kTRUE),
  fRnrCascadepath(kTRUE),
  fBacColor(0),
  fMinRCut(0),
  fMaxRCut(250),
  fMinDaughterDCA(0),
  fMaxDaughterDCA(1),
  fMinPt(0),
  fMaxPt(20)
{
  // Default constructor.

  fChildClass = AliEveCascade::Class(); // override member from base TEveElementList
}

//______________________________________________________________________________
AliEveCascadeList::AliEveCascadeList(TEveTrackPropagator* rs) :
  TEveElementList(),
  fTitle(),
  fRnrStyle(rs),
  fRnrDaughters(kTRUE),
  fRnrCascadevtx(kTRUE),
  fRnrCascadepath(kTRUE),
  fBacColor(0),
  fMinRCut(0),
  fMaxRCut(250),
  fMinDaughterDCA(0),
  fMaxDaughterDCA(1),
  fMinPt(0),
  fMaxPt(20)
{
  // Constructor with given track-propagator..

  fChildClass = AliEveCascade::Class(); // override member from base TEveElementList

  Init();
}

//______________________________________________________________________________
AliEveCascadeList::AliEveCascadeList(const Text_t* name, TEveTrackPropagator* rs) :
  TEveElementList(),
  fTitle(),
  fRnrStyle(rs),
  fRnrDaughters(kTRUE),
  fRnrCascadevtx(kTRUE),
  fRnrCascadepath(kTRUE),
  fBacColor(0),
  fMinRCut(0),
  fMaxRCut(100),
  fMinDaughterDCA(0),
  fMaxDaughterDCA(1),
  fMinPt(0),
  fMaxPt(20)
{
  // Standard constructor.

  fChildClass = AliEveCascade::Class(); // override member from base TEveElementList

  Init();
  SetName(name);
}

//______________________________________________________________________________
void AliEveCascadeList::Init()
{
  // Initialize members needed for drawing operations.

  if (fRnrStyle== 0) fRnrStyle = new TEveTrackPropagator;
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveCascadeList::MakeCascades()
{
  // Call MakeCascade() for all elements.

  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    ((AliEveCascade*)(*i))->MakeCascade();
  }
  gEve->Redraw3D();
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveCascadeList::FilterByRadius(Float_t minR, Float_t maxR)
{
  // Select visibility of elements based on their axial radius.

  fMinRCut = minR;
  fMaxRCut = maxR;

  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
  {
    AliEveCascade* cascade = (AliEveCascade*) *i;
    Float_t  rad = cascade->GetRadius();
    Bool_t  show = rad >= fMinRCut && rad <= fMaxRCut;
    cascade->SetRnrState(show);
  }
  ElementChanged();
  gEve->Redraw3D();
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveCascadeList::FilterByDaughterDCA(Float_t minDaughterDCA, Float_t maxDaughterDCA)
{
  // Select visibility of elements based on the DCA between daughters.

  fMinDaughterDCA = minDaughterDCA;
  fMaxDaughterDCA = maxDaughterDCA;

  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
  {
    AliEveCascade* cascade = (AliEveCascade*) *i;
    Float_t  dca = cascade->GetDaughterDCA();
    Bool_t  show = dca >= fMinDaughterDCA && dca <= fMaxDaughterDCA;
    cascade->SetRnrState(show);
  }
  ElementChanged();
  gEve->Redraw3D();
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveCascadeList::FilterByPt(Float_t minPt, Float_t maxPt)
{
  // Select visibility of elements based on the Cascade pt.

  fMinPt = minPt;
  fMaxPt = maxPt;

  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
  {
    AliEveCascade* cascade = (AliEveCascade*) *i;
    Float_t  pt = cascade->GetPt();
    Bool_t  show = pt >= fMinPt && pt <= fMaxPt;
    cascade->SetRnrState(show);
  }
  ElementChanged();
  gEve->Redraw3D();
}
