// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveV0.h"

#include <TEveTrack.h>
#include <TEveTrackPropagator.h>
#include <TEveManager.h>

#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TColor.h>

#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include <vector>


/***********************************************************************
*
*  AliEveV0 class
*
************************************************************************/

ClassImp(AliEveV0)

//______________________________________________________________________________
AliEveV0::AliEveV0() :
  TEvePointSet(),

  fRecBirthV(),
  fRecDecayV(),
  fRecDecayP(),
  fNegTrack(0),
  fPosTrack(0),
  fRnrStyle(0),
  fPointingLine(0),
  fESDIndex(-1),
  fOnFlyStatus(kFALSE),
  fDaughterDCA(999),
  fChi2V0(-1)
{
  // Default constructor.

  // Override from TEveElement.
  fPickable = kTRUE;
  fMainColorPtr = &fMarkerColor;
}

//______________________________________________________________________________
AliEveV0::AliEveV0(TEveRecTrack* tNeg, TEveRecTrack* tPos,
		   TEveRecV0* v0, TEveTrackPropagator* rs) :
  TEvePointSet(),

  fRecBirthV(v0->fV0Birth),
  fRecDecayV(v0->fVCa),
  fRecDecayP(v0->fPNeg + v0->fPPos),

  fNegTrack(new TEveTrack(tNeg, rs)),
  fPosTrack(new TEveTrack(tPos, rs)),

  fRnrStyle(rs),
  fPointingLine(new TEveLine("Pointing line")),
  fESDIndex(-1),
  fOnFlyStatus(kFALSE),
  fDaughterDCA(999),
  fChi2V0(-1)
{
  // Constructor with full V0 specification.

  // Override from TEveElement.
  fPickable = kTRUE;
  fMainColorPtr = &fMarkerColor;

  fMarkerStyle = 2;
  fMarkerColor = kSpring + 6;
  fMarkerSize  = 1;

  fPointingLine->SetLineColor(fMarkerColor);
  fPointingLine->SetLineWidth(2);
  fPointingLine->IncDenyDestroy();
  AddElement(fPointingLine);

  fPosTrack->SetLineColor(2);  // red
  fPosTrack->SetStdTitle();
  fNegTrack->SetLineColor(7);  // light blue
  fNegTrack->SetStdTitle();

  fNegTrack->IncDenyDestroy();
  AddElement(fNegTrack);
  fPosTrack->IncDenyDestroy();
  AddElement(fPosTrack);
}

//______________________________________________________________________________
AliEveV0::~AliEveV0()
{
  // Destructor. Dereferences pos/neg tracks and pointing-line objects.

  fNegTrack->DecDenyDestroy();
  fPosTrack->DecDenyDestroy();
  fPointingLine->DecDenyDestroy();
}

//______________________________________________________________________________
Float_t AliEveV0::GetInvMass(Float_t nPdgCode, Float_t pPdgCode) const
{
  // Returns Invariant Mass assuming the masses of the daughter particles
  TEveVector lNegMomentum = fNegTrack->GetMomentum();
  // Does not work properly because momenta at the primary vertex !!!!!!!
  TEveVector lPosMomentum = fPosTrack->GetMomentum();
  Double_t nMass=TDatabasePDG::Instance()->GetParticle(nPdgCode)->Mass();
  Double_t pMass=TDatabasePDG::Instance()->GetParticle(pPdgCode)->Mass();

  printf("\n check the mass of the particle negative %.5f positive %.5f \n",nMass,pMass);

  Double_t eNeg = TMath::Sqrt(nMass*nMass + lNegMomentum.Mag2());
  Double_t ePos = TMath::Sqrt(pMass*pMass + lPosMomentum.Mag2());

  return TMath::Sqrt( (eNeg+ePos)*(eNeg+ePos) - fRecDecayP.Mag2() );
}

//______________________________________________________________________________
void AliEveV0::MakeV0()
{
  // Set all dependant components for drawing.

  SetPoint(0, fRecDecayV.fX, fRecDecayV.fY, fRecDecayV.fZ);

  fNegTrack->MakeTrack();
  fPosTrack->MakeTrack();

  fPointingLine->SetPoint(0, fRecBirthV.fX, fRecBirthV.fY, fRecBirthV.fZ);
  fPointingLine->SetPoint(1, fRecDecayV.fX, fRecDecayV.fY, fRecDecayV.fZ);
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
  fPosColor(0),
  fMinRCut(0),
  fMaxRCut(250),
  fMinDaughterDCA(0),
  fMaxDaughterDCA(1),
  fMinPt(0),
  fMaxPt(20)
{
  // Default constructor.

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
  fPosColor(0),
  fMinRCut(0),
  fMaxRCut(250),
  fMinDaughterDCA(0),
  fMaxDaughterDCA(1),
  fMinPt(0),
  fMaxPt(20)
{
  // Constructor with given track-propagator..

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
  fPosColor(0),
  fMinRCut(0),
  fMaxRCut(100),
  fMinDaughterDCA(0),
  fMaxDaughterDCA(1),
  fMinPt(0),
  fMaxPt(20)
{
  // Standard constructor.

  fChildClass = AliEveV0::Class(); // override member from base TEveElementList

  Init();
  SetName(name);
}

//______________________________________________________________________________
void AliEveV0List::Init()
{
  // Initialize members needed for drawing operations.

  if (fRnrStyle== 0) fRnrStyle = new TEveTrackPropagator;
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveV0List::MakeV0s()
{
  // Call MakeV0() for all elements.

  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    ((AliEveV0*)(*i))->MakeV0();
  }
  gEve->Redraw3D();
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveV0List::FilterByRadius(Float_t minR, Float_t maxR)
{
  // Select visibility of elements based on their axial radius.

  fMinRCut = minR;
  fMaxRCut = maxR;

  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
  {
    AliEveV0* v0 = (AliEveV0*) *i;
    Float_t  rad = v0->GetRadius();
    Bool_t  show = rad >= fMinRCut && rad <= fMaxRCut;
    v0->SetRnrState(show);
  }
  ElementChanged();
  gEve->Redraw3D();
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveV0List::FilterByDaughterDCA(Float_t minDaughterDCA, Float_t maxDaughterDCA)
{
  // Select visibility of elements based on the DCA between daughters.

  fMinDaughterDCA = minDaughterDCA;
  fMaxDaughterDCA = maxDaughterDCA;

  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
  {
    AliEveV0* v0 = (AliEveV0*) *i;
    Float_t  dca = v0->GetDaughterDCA();
    Bool_t  show = dca >= fMinDaughterDCA && dca <= fMaxDaughterDCA;
    v0->SetRnrState(show);
  }
  ElementChanged();
  gEve->Redraw3D();
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveV0List::FilterByPt(Float_t minPt, Float_t maxPt)
{
  // Select visibility of elements based on the V0 pt.

  fMinPt = minPt;
  fMaxPt = maxPt;

  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
  {
    AliEveV0* v0 = (AliEveV0*) *i;
    Float_t  pt = v0->GetPt();
    Bool_t  show = pt >= fMinPt && pt <= fMaxPt;
    v0->SetRnrState(show);
  }
  ElementChanged();
  gEve->Redraw3D();
}
