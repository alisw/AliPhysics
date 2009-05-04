// $Id$
// Main authors: Paraskevi Ganoti: 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveKink.h"

#include <TEveTrack.h>
#include <TEveTrackPropagator.h>

#include <TEveManager.h>

/***********************************************************************
*
*  AliEveKink class
*
************************************************************************/

ClassImp(AliEveKink)

//______________________________________________________________________________
AliEveKink::AliEveKink() :
  TEvePointSet(),

  fRecKinkPosition(),
  fMotherMomentum(),
  fDaughterMomentum(),
  fMotherTrack(0),
  fDaughterTrack(0),
  fRnrStyleMoth(0),
  fRnrStyleDaugh(0),
  fESDKinkIndex(-1),
  fDaugMaxProbPdg(0),
  fDaugMaxProbPid(0)
{
  // Default constructor.

  // Override from TEveElement.
  fKinkAngle[0]=fKinkAngle[1]=fKinkAngle[2]=0.0;
  fPickable = kTRUE;
}

//______________________________________________________________________________
AliEveKink::AliEveKink(TEveRecTrack* tMoth, TEveRecTrack* tDaug,
		   TEveRecKink* kink, TEveTrackPropagator* rsMoth, TEveTrackPropagator* rsDaugh) :
  TEvePointSet(),
  fRecKinkPosition(kink->fVKink),
  fMotherMomentum(kink->fPMother),
  fDaughterMomentum(kink->fPDaughter),
  fMotherTrack(new TEveTrack(tMoth, rsMoth)),
  fDaughterTrack(new TEveTrack(tDaug, rsDaugh)),
  fRnrStyleMoth(rsMoth),
  fRnrStyleDaugh(rsDaugh),
  fESDKinkIndex(-1),
  fDaugMaxProbPdg(0),
  fDaugMaxProbPid(0)
{
  // Constructor with full kink specification.

  // Override from TEveElement.
  fPickable = kTRUE;
  fKinkAngle[0]=fKinkAngle[1]=fKinkAngle[2]=0;

  fMotherTrack->SetLineColor(2);  // red
  fMotherTrack->SetStdTitle();
  fMotherTrack->SetLineWidth(3);
  TEvePathMark* pmM = 0;
  pmM = new TEvePathMark(TEvePathMark::kDecay);
  pmM->fV.Set(kink->fVKink);
  fMotherTrack->AddPathMark(*pmM);

  fDaughterTrack->SetLineColor(7);  // light blue
  fDaughterTrack->SetLineWidth(3);  //
  fDaughterTrack->SetStdTitle();

  fMotherTrack->IncDenyDestroy();
  AddElement(fMotherTrack);
  fDaughterTrack->IncDenyDestroy();
  AddElement(fDaughterTrack);
}

//______________________________________________________________________________
AliEveKink::~AliEveKink()
{
  // Destructor. Dereferences mother-daughter tracks

  fMotherTrack->DecDenyDestroy();
  fDaughterTrack->DecDenyDestroy();
}
//______________________________________________________________________________
void AliEveKink::SetMaxProbPdgPid(Int_t rPdg, Float_t rPid)
{
  // Sets the maximum probability Pdg value and Pid for the daughter
  // Should be moved to TEveTrack property eventually (or AliEveTrack creation)

    fDaugMaxProbPdg = rPdg;
    fDaugMaxProbPid = rPid;

}

//______________________________________________________________________________
Float_t AliEveKink::GetInvMass(Int_t dPdgCode) const
{
  // Returns Invariant Mass assuming the mass of the daughter particle
  TEveVector lDaugMomentum = fDaughterTrack->GetMomentum();
  TEveVector lMothMomentum = fMotherTrack->GetMomentum();
  Double_t dMass=TDatabasePDG::Instance()->GetParticle(dPdgCode)->Mass();

  Double_t eDaug = TMath::Sqrt(dMass*dMass + lDaugMomentum.Mag2());
  Double_t otherDaug=TMath::Sqrt(
  (lMothMomentum.fX-lDaugMomentum.fX)*(lMothMomentum.fX-lDaugMomentum.fX)+
  (lMothMomentum.fY-lDaugMomentum.fY)*(lMothMomentum.fY-lDaugMomentum.fY)+
  (lMothMomentum.fZ-lDaugMomentum.fZ)*(lMothMomentum.fZ-lDaugMomentum.fZ));

  return TMath::Sqrt( (eDaug+otherDaug)*(eDaug+otherDaug) - lMothMomentum.Mag2());
}

//______________________________________________________________________________
Float_t AliEveKink::GetQt() const
{
  // Returns the kink daughter Qt

  TEveVector lDaugMomentum = fDaughterTrack->GetMomentum();
  Float_t daugMom=TMath::Sqrt(
  (lDaugMomentum.fX)*(lDaugMomentum.fX)+
  (lDaugMomentum.fY)*(lDaugMomentum.fY)+
  (lDaugMomentum.fZ)*(lDaugMomentum.fZ));

  return TMath::Sin(fKinkAngle[2])*daugMom ;
}

//______________________________________________________________________________
void AliEveKink::MakeKink()
{
  // Set all dependant components for drawing.

  SetPoint(0, fRecKinkPosition.fX, fRecKinkPosition.fY, fRecKinkPosition.fZ);

  fMotherTrack->MakeTrack();
  fDaughterTrack->MakeTrack();

}


/***********************************************************************
*
*  AliEveKinkList class
*
************************************************************************/

ClassImp(AliEveKinkList)

//______________________________________________________________________________
AliEveKinkList::AliEveKinkList() :
  TEveElementList(),
  fRnrStyleMoth(0),
  fRnrStyleDaugh(0),
  fRnrKinkDaughter(kTRUE),
  fRnrKinkvtx(kTRUE),
  fMothColor(0),
  fDaugColor(0),
  fMinRCut(0),
  fMaxRCut(300),
  fMinKinkAngle(0),
  fMaxKinkAngle(300),
  fMinPt(0),
  fMaxPt(20),
  fMinInvariantMass(0),
  fMaxInvariantMass(1.0),
  fDaugCheckedPid(0),
  fDaugCheckedProb(0)
{
  // Default constructor.

  fChildClass = AliEveKink::Class(); // override member from base TEveElementList
}

//______________________________________________________________________________
AliEveKinkList::AliEveKinkList(TEveTrackPropagator* rsMoth,TEveTrackPropagator* rsDaugh) :
  TEveElementList(),
  fRnrStyleMoth(rsMoth),
  fRnrStyleDaugh(rsDaugh),
  fRnrKinkDaughter(kTRUE),
  fRnrKinkvtx(kTRUE),
  fMothColor(0),
  fDaugColor(0),
  fMinRCut(0),
  fMaxRCut(300),
  fMinKinkAngle(0),
  fMaxKinkAngle(300),
  fMinPt(0),
  fMaxPt(20),
  fMinInvariantMass(0),
  fMaxInvariantMass(1.0),
  fDaugCheckedPid(0),
  fDaugCheckedProb(0)
{
  // Constructor with given track-propagator..

  fChildClass = AliEveKink::Class(); // override member from base TEveElementList

  Init();
}

//______________________________________________________________________________
AliEveKinkList::AliEveKinkList(const Text_t* name, TEveTrackPropagator* rsMoth, TEveTrackPropagator* rsDaugh) :
  TEveElementList(),
  fRnrStyleMoth(rsMoth),
  fRnrStyleDaugh(rsDaugh),
  fRnrKinkDaughter(kTRUE),
  fRnrKinkvtx(kTRUE),
  fMothColor(0),
  fDaugColor(0),
  fMinRCut(0),
  fMaxRCut(300),
  fMinKinkAngle(0),
  fMaxKinkAngle(300),
  fMinPt(0),
  fMaxPt(20),
  fMinInvariantMass(0),
  fMaxInvariantMass(1.0),
  fDaugCheckedPid(0),
  fDaugCheckedProb(0)
{
  // Standard constructor.

  fChildClass = AliEveKink::Class(); // override member from base TEveElementList

  Init();
  SetName(name);
}

//______________________________________________________________________________
void AliEveKinkList::Init()
{
  // Initialize members needed for drawing operations.

  if (fRnrStyleMoth== 0) fRnrStyleMoth = new TEveTrackPropagator;
  if (fRnrStyleDaugh== 0) fRnrStyleDaugh = new TEveTrackPropagator;
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveKinkList::MakeKinks()
{
  // Call MakeKink() for all elements.

  for(List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
    ((AliEveKink*)(*i))->MakeKink();
  }
  gEve->Redraw3D();
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveKinkList::FilterByRadius(Float_t minR, Float_t maxR)
{
  // Select visibility of elements based on their axial radius.

  fMinRCut = minR;
  fMaxRCut = maxR;

  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
  {
    AliEveKink* kink = (AliEveKink*) *i;
    Float_t  rad = kink->GetKinkRadius();
    Bool_t  show = rad >= fMinRCut && rad <= fMaxRCut;
    kink->SetRnrState(show);
  }
  ElementChanged();
  gEve->Redraw3D();
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveKinkList::FilterByKinkAngle(Float_t minKinkAngle, Float_t maxKinkAngle)
{
  // Select visibility of elements based on the kink angle.

  fMinKinkAngle = minKinkAngle;
  fMaxKinkAngle = maxKinkAngle;

  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
  {
    AliEveKink* kink = (AliEveKink*) *i;
    Float_t  angle = TMath::RadToDeg()*kink->GetKinkAngle(2);
    Bool_t  show = angle >= fMinKinkAngle && angle <= fMaxKinkAngle;
    kink->SetRnrState(show);
  }
  ElementChanged();
  gEve->Redraw3D();
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveKinkList::FilterByPt(Float_t minPt, Float_t maxPt)
{
  // Select visibility of elements based on the kink pt.

  fMinPt = minPt;
  fMaxPt = maxPt;

  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
  {
    AliEveKink* kink = (AliEveKink*) *i;
    Float_t  pt = kink->GetKinkPMotherPerp();
    Bool_t  show = pt >= fMinPt && pt <= fMaxPt;
    kink->SetRnrState(show);
  }
  ElementChanged();
  gEve->Redraw3D();
}

//______________________________________________________________________________
void AliEveKinkList::FilterByInvariantMass(Float_t minInvariantMass, Float_t maxInvariantMass, Int_t dPdgCode)
{
  // Select visibility of elements based on the kink invariant mass.

  fMinInvariantMass = minInvariantMass;
  fMaxInvariantMass = maxInvariantMass;

  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
  {
    AliEveKink* kink = (AliEveKink*) *i;
    Float_t  invMass = kink->GetInvMass(dPdgCode);
    Bool_t  show = invMass >= fMinInvariantMass && invMass <= fMaxInvariantMass;
    kink->SetRnrState(show);
  }
  ElementChanged();
  gEve->Redraw3D();
}

//______________________________________________________________________________
void AliEveKinkList::FilterByCheckedPidMinProb(Int_t rFlag, Int_t rPid, Float_t rProb)
{
  // Select visibility of elements based on the kink daughter PID.

  fDaugCheckedPid  = rPid;
  fDaugCheckedProb = rProb;

  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
  {
    AliEveKink* kink = (AliEveKink*) *i;
    Int_t   pid  = 0;
    Float_t prob = 0;
    Bool_t  show = kFALSE;

    pid  = kink->GetDaugMaxProbPdg();
    prob = kink->GetDaugMaxProbPid();
    show = (pid == fDaugCheckedPid && prob > fDaugCheckedProb) || !rFlag ;

    kink->SetRnrState(show);
  }
  ElementChanged();
  gEve->Redraw3D();
}
