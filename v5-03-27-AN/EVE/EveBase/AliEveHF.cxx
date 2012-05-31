// $Id$
// Main author: Davide Caffarri 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * ee http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveHF.h"
#include "AliAODRecoDecay.h"

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
*  AliEveHF class
*
************************************************************************/


ClassImp(AliEveHF)

//______________________________________________________________________________
AliEveHF::AliEveHF():
  TEvePointSet(),

  fAODobj(0x0),
  fRecBirthHF(0,0,0),
  fRecDecayHF(0,0,0),
  fRecDecayMomHF(0,0,0),
  fPointingAngleHF(0),

  fNegTrack(0),
  fPosTrack(0),
  fRnrStyle(0),
  fPointingLine(0),

  fnProng(0),
  fAODIndex(-1),
  fChi2SecondVtx(-1),

  fProngDCA(0x0),
  fProngd0(0x0),
  fProngMaxProbPdg(0x0),
  fProngMaxProbPid(0x0),
  fInvariantMassPart(0),
  fInvariantMassAntiPart(0),
  fDecay(0)

{
  // Default constructor.

  // Override from TEveElement.
  fPickable = kTRUE;
  fMainColorPtr = &fMarkerColor;
}

//______________________________________________________________________________
AliEveHF::AliEveHF(TEveRecTrack* tNeg, TEveRecTrack* tPos, Double_t primVtx[3], AliAODRecoDecay* aodObj, TEveTrackPropagator* rs) :

  TEvePointSet(),

  fAODobj(aodObj),
  fRecBirthHF(primVtx[0], primVtx[1], primVtx[2]),
  fRecDecayHF(aodObj->GetSecVtxX(), aodObj->GetSecVtxY(), aodObj->GetSecVtxZ()),
  fRecDecayMomHF(aodObj->Px(),aodObj->Py(),aodObj->Pz()),
  fPointingAngleHF(aodObj->CosPointingAngle(primVtx)),

  fNegTrack(new TEveTrack(tNeg, rs)),
  fPosTrack(new TEveTrack(tPos, rs)),
  fRnrStyle(rs),
  fPointingLine(new TEveLine("Pointing line")),

  fnProng(aodObj->GetNProngs()),
  fAODIndex(-1),
  fChi2SecondVtx(aodObj->GetReducedChi2()),

  fProngDCA(0),
  fProngd0(0),
  fProngMaxProbPdg(0),
  fProngMaxProbPid(0),
  fInvariantMassPart(0),
  fInvariantMassAntiPart(0),
  fDecay(0)
{

  // Constructor with full HF specification.
  // Override from TEveElement.
  fPickable = kTRUE;
  fMainColorPtr = &fMarkerColor;

  fMarkerStyle = 2;
  fMarkerColor = kSpring + 6;
  fMarkerSize  = 0.3;

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
AliEveHF::~AliEveHF()
{
  // Destructor. Dereferences pos/neg tracks and pointing-line objects.

  fNegTrack->DecDenyDestroy();
  fPosTrack->DecDenyDestroy();
  fPointingLine->DecDenyDestroy();
}

//_____________________________________________________________________________
void AliEveHF::SetProngDCA() const
{
  for(Int_t ip=0; ip<fnProng; ip++)
    fProngDCA[ip] = fAODobj->GetDCA(ip);
}

//______________________________________________________________________________
void AliEveHF::Setd0Prong() const
{
  for(Int_t ip=0; ip<fnProng; ip++)
    fProngd0[ip] = fAODobj->Getd0Prong(ip);
}


//______________________________________________________________________________
void AliEveHF::SetMaxProbPdgPid()
{
  // Sets the maximum probability Pdg value and Pid for one daughter
  // Should be moved to TEveTrack property eventually (or AliEveTrack creation)
  Double_t pid[5];
  Int_t    pos = -1;

  for (Int_t ip=0; ip<fnProng; ip++){
    fAODobj->GetPIDProng(ip, pid);

    fProngMaxProbPid[ip]=pid[0]; 
    for (Int_t pp=1; pp<5; pp++)
      if (pid[pp]>pid[pp-1]) {
	fProngMaxProbPid[ip]=pid[pp];
	pos = pp;}
    switch (pos)
      {
      case 0:
	fProngMaxProbPdg[ip] = -11;
	break;
      case 1:
	fProngMaxProbPdg[ip] = -13;
    break;
      case 2:
	fProngMaxProbPdg[ip] = 211;
    break;
      case 3:
	fProngMaxProbPdg[ip] = 321;
    break;
      case 4:
	fProngMaxProbPdg[ip] = 2212;
    break;
      }
  }

}

//________________________________________________________________________________________
void AliEveHF::CalculateInvMass(Int_t decay)
{

  //Method to calculate the invariant mass of the particle and the antiparticle. 
  UInt_t   pdg2[2];
  Double_t mPDG,minv;

  switch (decay)
    {
    case 0:                  // D0->Kpi
      pdg2[0]=211; pdg2[1]=321;
      mPDG=TDatabasePDG::Instance()->GetParticle(421)->Mass();
      minv = fAODobj->InvMass(fnProng,pdg2);
      fInvariantMassPart=minv;

      pdg2[0]=321; pdg2[1]=211;
      minv = fAODobj->InvMass(fnProng,pdg2);
      fInvariantMassAntiPart=minv;
      break;
    }

}

//______________________________________________________________________________
Bool_t AliEveHF::SelectInvMass(Int_t decay, Float_t decayCuts)
{

  //Member fuction to select particles using the invariant mass cuts. 
  UInt_t   pdg2[2];
  Double_t mPDG,minv;

  Bool_t retval=kFALSE;
  switch (decay)
    {
    case 0:                  // D0->Kpi
      pdg2[0]=211; pdg2[1]=321;
      mPDG=TDatabasePDG::Instance()->GetParticle(421)->Mass();
      minv = fAODobj->InvMass(fnProng,pdg2);
      fInvariantMassPart=minv;
      if(TMath::Abs(minv-mPDG)<decayCuts) retval=kTRUE;
      pdg2[0]=321; pdg2[1]=211;
      minv = fAODobj->InvMass(fnProng,pdg2);
      fInvariantMassAntiPart=minv;
      if(TMath::Abs(minv-mPDG)<decayCuts) retval=kTRUE;
      break;
    }
  return retval;

}


//______________________________________________________________________________
void AliEveHF::MakeHF()
{
  // Set all dependant components for drawing.

  SetPoint(0, fRecDecayHF.fX, fRecDecayHF.fY, fRecDecayHF.fZ);

  fNegTrack->MakeTrack();
  fPosTrack->MakeTrack();

  fPointingLine->SetPoint(0, fRecBirthHF.fX, fRecBirthHF.fY, fRecBirthHF.fZ);
  fPointingLine->SetPoint(1, fRecDecayHF.fX, fRecDecayHF.fY, fRecDecayHF.fZ);
}

/***********************************************************************
*
*  AliEveHFList class
*
************************************************************************/

ClassImp(AliEveHFList)

//______________________________________________________________________________
AliEveHFList::AliEveHFList() :
  TEveElementList(),
  fTitle(),
  fRnrStyle(0),
  fRnrDaughters(kTRUE),
  fRnrHFvtx(kTRUE),
  fRnrHFpath(kTRUE),
  fProngColor(0),
  fMinRCut(0),
  fMaxRCut(250),
  fMinDaughterDCA(0),
  fMaxDaughterDCA(1),
  fMinPt(0),
  fMaxPt(20),
  fMinCosPointingAngle(0),
  fMaxCosPointingAngle(1),
  fMind0(-3),
  fMaxd0(3),
  fProngCheckedPid(0x0),
  fProngCheckedProb(0x0),
  fDeltaInvariantMass(0),
  fDecay(0)
{
  // Default constructor.

  fChildClass = AliEveHF::Class(); // override member from base TEveElementList
}

//______________________________________________________________________________
AliEveHFList::AliEveHFList(TEveTrackPropagator* rs) :
  TEveElementList(),
  fTitle(),
  fRnrStyle(rs),
  fRnrDaughters(kTRUE),
  fRnrHFvtx(kTRUE),
  fRnrHFpath(kTRUE),
  fProngColor(0),
  fMinRCut(0),
  fMaxRCut(250),
  fMinDaughterDCA(0),
  fMaxDaughterDCA(1),
  fMinPt(0),
  fMaxPt(20),
  fMinCosPointingAngle(0),
  fMaxCosPointingAngle(1),
  fMind0(-3),
  fMaxd0(3),
  fProngCheckedPid(0x0),
  fProngCheckedProb(0x0),
  fDeltaInvariantMass(1),
  fDecay(0)
{
  // Constructor with given track-propagator..

  fChildClass = AliEveHF::Class(); // override member from base TEveElementList

  Init();
}

//______________________________________________________________________________
AliEveHFList::AliEveHFList(const Text_t* name, TEveTrackPropagator* rs) :
  TEveElementList(),
  fTitle(),
  fRnrStyle(rs),
  fRnrDaughters(kTRUE),
  fRnrHFvtx(kTRUE),
  fRnrHFpath(kTRUE),
  fProngColor(0),
  fMinRCut(0),
  fMaxRCut(100),
  fMinDaughterDCA(0),
  fMaxDaughterDCA(1),
  fMinPt(0),
  fMaxPt(20),
  fMinCosPointingAngle(0),
  fMaxCosPointingAngle(1),
  fMind0(-3),
  fMaxd0(3),
  fProngCheckedPid(0x0),
  fProngCheckedProb(0x0),
  fDeltaInvariantMass(1),
  fDecay(0)
{
  // Standard constructor.

  fChildClass = AliEveHF::Class(); // override member from base TEveElementList

  Init();
  SetName(name);
}

//______________________________________________________________________________
void AliEveHFList::Init()
{
  // Initialize members needed for drawing operations.

  if (fRnrStyle== 0) fRnrStyle = new TEveTrackPropagator;
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveHFList::MakeHFs()
{
  // Call MakeHF() for all elements.

  for(List_i i=fChildren.begin(); i!=fChildren.end(); i++) {
    ((AliEveHF*)(*i))->MakeHF();
  }
  gEve->Redraw3D();
}

/******************************************************************************/


//______________________________________________________________________________
void AliEveHFList::FilterByPt(Float_t minPt, Float_t maxPt)
{

  //Select visibility of elements based on their Pt

  fMinPt = minPt;
  fMaxPt = maxPt;

  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
    {
      AliEveHF *hf = (AliEveHF*) *i;
      Float_t  pt  = hf->GetPt();
      Bool_t  show = ((pt >= fMinPt) && (pt <= fMaxPt));
      hf->SetRnrState(show);
    }
  ElementChanged();
  gEve->Redraw3D();
}
//______________________________________________________________________________
void AliEveHFList::FilterByRadius(Float_t minR, Float_t maxR)
{
  // Select visibility of elements based on their axial radius.

  fMinRCut = minR;
  fMaxRCut = maxR;

  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
    {
      AliEveHF *hf = (AliEveHF*) *i;
      Float_t  rad = hf->GetRadius();
      Bool_t  show = rad >= fMinRCut && rad <= fMaxRCut;
      hf->SetRnrState(show);
    }
  ElementChanged();
  gEve->Redraw3D();
}

/******************************************************************************/


//______________________________________________________________________________
void AliEveHFList::FilterByCosPointingAngle(Float_t minCosPointingAngle, Float_t maxCosPointingAngle)
{
  // Select visibility of elements based on the HF cosine of the pointing angle

  fMinCosPointingAngle = minCosPointingAngle;
  fMaxCosPointingAngle = maxCosPointingAngle;


  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
    {
      AliEveHF *hf = (AliEveHF*) *i;
      Float_t  cosPointingAngle = hf->GetCosPointingAngle();
      Bool_t  show = cosPointingAngle >= fMinCosPointingAngle && cosPointingAngle <= fMaxCosPointingAngle;
      hf->SetRnrState(show);

      ElementChanged();
      gEve->Redraw3D();
    }
}
//______________________________________________________________________________
void AliEveHFList::FilterByd0(Float_t mind0, Float_t maxd0)
{
  // Select visibility of elements based on the HF impact parameter.

  fMind0 = mind0;
  fMaxd0 = maxd0;

  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
    {
      AliEveHF *hf = (AliEveHF*) *i;
      Int_t nProng=(Int_t)hf->GetAODobj()->GetNProngs();
      for(Int_t ip = 1; ip < nProng; ++ip) {
	Double_t  d0 = hf->Getd0Prong(ip);
	Bool_t  show = d0 >= fMind0 && d0 <= fMaxd0;
	hf->SetRnrState(show);
      }
      ElementChanged();
      gEve->Redraw3D();
    }
}


 //______________________________________________________________________________
void AliEveHFList::FilterByDCA(Float_t minDaughterDCA, Float_t maxDaughterDCA)
{
  // Select visibility of elements based on the DCA between daughters.

  fMinDaughterDCA = minDaughterDCA;
  fMaxDaughterDCA = maxDaughterDCA;

  Double_t dca;
  AliEveHF *hf();


  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
  {
    AliEveHF* hf = (AliEveHF*) *i;
    Int_t nProng=(Int_t)hf->GetAODobj()->GetNProngs();
    for(Int_t ip = 1; ip < nProng; ++ip) {
      dca = hf->GetProngDCA(ip);
      Bool_t  show = dca >= fMinDaughterDCA && dca <= fMaxDaughterDCA;
      hf->SetRnrState(show);
    }
  }
  ElementChanged();
  gEve->Redraw3D();
}

/******************************************************************************/
//______________________________________________________________________________
/*void AliEveHFList::FilterByCheckedPidMinProb(Int_t rFlag,Int_t rDaughter, Int_t rPid, Float_t rProb)
{

   if (!rDaughter){
    fNegCheckedPid  = rPid;
    fNegCheckedProb = rProb;
  }
  else {
    fPosCheckedPid  = rPid;
    fPosCheckedProb = rProb;
  }

  // Select visibility of elements based on one of the V0 daughters PID
  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
  {
    AliEveHF* hf = (AliEveHF*) *i;
    Int_t   pid = 0;
    Float_t prob = 0.0;
    Bool_t  show = 0;
    if     (!rDaughter) {// Negative daughter checked
      pid  = hf->GetNegMaxProbPdg();
      prob = hf->GetNegMaxProbPid();
      show = (pid == fNegCheckedPid && prob > fNegCheckedProb) || !rFlag ;
    }
    else if (rDaughter) {// Positive daughter checked
      pid = hf->GetPosMaxProbPdg();
      prob = hf->GetPosMaxProbPid();
      show = (pid == fPosCheckedPid && prob > fPosCheckedProb) || !rFlag ;
    }
    hf->SetRnrState(show);
exit  }
  ElementChanged();
  gEve->Redraw3D();
  }*/

//______________________________________________________________________________

void AliEveHFList::FilterByInvariantMass(Int_t decay, Float_t deltaInvariantMass)
{
  // Select visibility of elements based on the HF  invariant mass.

  fDeltaInvariantMass = deltaInvariantMass;
  fDecay = decay;

  for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
  {
    AliEveHF* hf = (AliEveHF*) *i;
    Bool_t show  = hf->SelectInvMass(decay, fDeltaInvariantMass);
    // if (show) show = (invMass >= (fInvariantMass-fDeltaInvariantMass)) && (invMass <= (fInvariantMass+fDeltaInvariantMass));
    hf->SetRnrState(show);
  }
  ElementChanged();
  gEve->Redraw3D();
}
