// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveV0.h"

#include <TEveTrackPropagator.h>
#include <TEveManager.h>

#include <TPolyLine3D.h>
#include <TColor.h>

#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include <vector>

//______________________________________________________________________________
// Representation of a reconstructed V0.

ClassImp(AliEveV0)

//______________________________________________________________________________
AliEveV0::AliEveV0() :
  TEvePointSet(),

  fRecBirthV(),
  fRecDecayV(),
  fRecDecayP(),
  fNegTrack(0),
  fPosTrack(0),
  fRnrStyleNeg(0),
  fRnrStylePos(0),
  fPointingLine(0),
  fESDIndex(-1),
  fOnFlyStatus(kFALSE),
  fDaughterDCA(999),
  fChi2V0(-1),
  fNegMaxProbPdg(0),
  fPosMaxProbPdg(0),
  fNegMaxProbPid(0),
  fPosMaxProbPid(0)
{
  // Default constructor.

  // Override from TEveElement.
  fPickable = kTRUE;
  fMainColorPtr = &fMarkerColor;
}

//______________________________________________________________________________
AliEveV0::AliEveV0(TEveRecTrack* tNeg, TEveRecTrack* tPos,
		   TEveRecV0* v0, TEveTrackPropagator* rsNeg, TEveTrackPropagator* rsPos) :
  TEvePointSet(),

  fRecBirthV(v0->fV0Birth),
  fRecDecayV(v0->fVCa),
  fRecDecayP(v0->fPNeg + v0->fPPos),

  fNegTrack(new AliEveTrack(tNeg, rsNeg)),
  fPosTrack(new AliEveTrack(tPos, rsPos)),

  fRnrStyleNeg(rsNeg),
  fRnrStylePos(rsPos),
  fPointingLine(new TEveLine("Pointing line")),
  fESDIndex(-1),
  fOnFlyStatus(kFALSE),
  fDaughterDCA(999),
  fChi2V0(-1),
  fNegMaxProbPdg(0),
  fPosMaxProbPdg(0),
  fNegMaxProbPid(0),
  fPosMaxProbPid(0)
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
void AliEveV0::SetMaxProbPdgPid(Int_t iDaughter, Int_t rPdg, Float_t rPid)
{
  // Sets the maximum probability Pdg value and Pid for one daughter
  // Should be moved to AliEveTrack property eventually (or AliEveTrack creation)

  if(iDaughter==0){
    fNegMaxProbPdg = rPdg;
    fNegMaxProbPid = rPid;
  }
  else if (iDaughter==1){
    fPosMaxProbPdg = rPdg;
    fPosMaxProbPid = rPid;
  }

}
//______________________________________________________________________________
Float_t AliEveV0::GetInvMass(Int_t nPdgCode, Int_t pPdgCode) const
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


//==============================================================================
//==============================================================================

//______________________________________________________________________________
//
// Container for AliEveV0s.
//
// Allows runtime selection by pT and DCA of daughters, radius of V0
// creation and PID priobabilities

ClassImp(AliEveV0List)

//______________________________________________________________________________
AliEveV0List::AliEveV0List() :
  TEveElementList(),
  fTitle(),
  fRnrStyleNeg(0),
  fRnrStylePos(0),
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
  fMaxPt(20),
  fNegCheckedPid(211),
  fPosCheckedPid(211),
  fNegCheckedProb(0.5),
  fPosCheckedProb(0.5),
  fMinInvariantMass(0),
  fMaxInvariantMass(1.2)
{
  // Default constructor.

  fChildClass = AliEveV0::Class(); // override member from base TEveElementList
}

//______________________________________________________________________________
AliEveV0List::AliEveV0List(TEveTrackPropagator* rsNeg,TEveTrackPropagator* rsPos) :
  TEveElementList(),
  fTitle(),
  fRnrStyleNeg(rsNeg),
  fRnrStylePos(rsPos),
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
  fMaxPt(20),
  fNegCheckedPid(211),
  fPosCheckedPid(211),
  fNegCheckedProb(0.5),
  fPosCheckedProb(0.5),
  fMinInvariantMass(0),
  fMaxInvariantMass(1.2)
{
  // Constructor with given track-propagator..

  fChildClass = AliEveV0::Class(); // override member from base TEveElementList

  Init();
}

//______________________________________________________________________________
AliEveV0List::AliEveV0List(const Text_t* name, TEveTrackPropagator* rsNeg,TEveTrackPropagator* rsPos) :
  TEveElementList(),
  fTitle(),
  fRnrStyleNeg(rsNeg),
  fRnrStylePos(rsPos),
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
  fMaxPt(20),
  fNegCheckedPid(211),
  fPosCheckedPid(211),
  fNegCheckedProb(0.5),
  fPosCheckedProb(0.5),
  fMinInvariantMass(0),
  fMaxInvariantMass(1.2)
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

  if (fRnrStyleNeg== 0) fRnrStyleNeg = new TEveTrackPropagator;
  if (fRnrStylePos== 0) fRnrStylePos = new TEveTrackPropagator;
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

//______________________________________________________________________________
void AliEveV0List::FilterByCheckedPidMinProb(Int_t rFlag,Int_t rDaughter, Int_t rPid, Float_t rProb)
{
  // Select visibility of elements based on one of the V0 daughters PID

  if (!rDaughter){
    fNegCheckedPid  = rPid;
    fNegCheckedProb = rProb;
  }
  else {
    fPosCheckedPid  = rPid;
    fPosCheckedProb = rProb;
  }

  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
  {
    AliEveV0* v0 = (AliEveV0*) *i;
    Int_t   pid = 0;
    Float_t prob = 0.0;
    Bool_t  show = 0;
    if     (!rDaughter) {// Negative daughter checked
      pid  = v0->GetNegMaxProbPdg();
      prob = v0->GetNegMaxProbPid();
      show = (pid == fNegCheckedPid && prob > fNegCheckedProb) || !rFlag ; 
    }
    else if (rDaughter) {// Positive daughter checked
      pid = v0->GetPosMaxProbPdg();
      prob = v0->GetPosMaxProbPid();
      show = (pid == fPosCheckedPid && prob > fPosCheckedProb) || !rFlag ; 
    }
    v0->SetRnrState(show);
  }
  ElementChanged();
  gEve->Redraw3D();
}

//______________________________________________________________________________
void AliEveV0List::FilterByInvariantMass(Float_t minInvariantMass, Float_t maxInvariantMass, Int_t nPdgCode, Int_t pPdgCode)
{
  // Select visibility of elements based on the V0 invariant mass.

  fMinInvariantMass = minInvariantMass;
  fMaxInvariantMass = maxInvariantMass;

  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
  {
    AliEveV0* v0 = (AliEveV0*) *i;
    Float_t  invMass = v0->GetInvMass(nPdgCode, pPdgCode);
    Bool_t  show = invMass >= fMinInvariantMass && invMass <= fMaxInvariantMass;
    v0->SetRnrState(show);
  }
  ElementChanged();
  gEve->Redraw3D();
}
