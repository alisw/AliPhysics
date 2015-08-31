// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

//------------------------------------------------------------------------
//
// AliEveCascade class
//
//------------------------------------------------------------------------

//#include <TPolyLine3D.h>
//#include <TColor.h>

//#include <vector>

#include <TDatabasePDG.h>
#include <TEveTrackPropagator.h>
#include <TEveManager.h>

#include "AliEveCascade.h"




ClassImp(AliEveCascade)

//______________________________________________________________________________
AliEveCascade::AliEveCascade() :
  TEvePointSet(),

  fRecBirthV(),
  fRecDecayV(),
  fRecDecayP(),
  fRecDecayV0(),
  fBacTrack(0),
  fNegTrack(0),
  fPosTrack(0),
  fRnrStyleBac(0),
  fRnrStyleNeg(0),
  fRnrStylePos(0),
  fPointingCurve(0),
  fV0Path(0),
  fESDIndex(-1),
  fDaughterDCA(999),
  fChi2Cascade(-1),
  fLambdaP(0.,0.,0.),
  fBachP(0.,0.,0.)
{
  // Default constructor.

  // Override from TEveElement.
  fPickable = kTRUE;
  fMainColorPtr = &fMarkerColor;
}

//______________________________________________________________________________
AliEveCascade::AliEveCascade(TEveRecTrack* tBac, TEveRecTrack* tNeg, TEveRecTrack* tPos, TEveRecV0* v0, TEveRecCascade* cascade, TEveTrackPropagator* rsBac, TEveTrackPropagator* rsNeg, TEveTrackPropagator* rsPos) :
  TEvePointSet(),

  fRecBirthV(cascade->fCascadeBirth),
  fRecDecayV(cascade->fCascadeVCa),
  fRecDecayP(cascade->fPBac + v0->fPNeg + v0->fPPos),
  fRecDecayV0(v0->fVCa),
  fBacTrack(new AliEveTrack(tBac, rsBac)),
  fNegTrack(new AliEveTrack(tNeg, rsNeg)),
  fPosTrack(new AliEveTrack(tPos, rsPos)),
  

  fRnrStyleBac(rsBac),
  fRnrStyleNeg(rsNeg),
  fRnrStylePos(rsPos),
  fPointingCurve(new TEveLine("Pointing Curve")),
  fV0Path(new TEveLine("V0 Path")),
  fESDIndex(-1),
  fDaughterDCA(999), // DCA between Xi daughters is properly set in esd_make_cascade via the setter
  fChi2Cascade(-1),
  fLambdaP(0.,0.,0.), // Momemtum of Lambda is properly set in esd_make_cascade via the setter
  fBachP(0.,0.,0.)    // Momemtum of Lambda is properly set in esd_make_cascade via the setter
{
  // Constructor with full Cascade specification.

  // Override from TEveElement.
  fPickable = kTRUE;
  fMainColorPtr = &fMarkerColor;

  fMarkerStyle = 2;
  fMarkerColor = kMagenta+2;
  fMarkerSize  = 1;

  fPointingCurve->SetTitle("Xi pointing curve");
  fPointingCurve->SetLineColor(fMarkerColor);
  fPointingCurve->SetLineStyle(9);
  fPointingCurve->SetLineWidth(2);  
  fPointingCurve->IncDenyDestroy();
  AddElement(fPointingCurve);

  fV0Path->SetTitle("V0 path");
  fV0Path->SetLineColor(fMarkerColor);
  fV0Path->SetLineStyle(3);
  fV0Path->SetLineWidth(2);
  fV0Path->IncDenyDestroy();
  AddElement(fV0Path);

  fBacTrack->SetLineColor(kMagenta);
  fBacTrack->SetStdTitle();
  fBacTrack->IncDenyDestroy();
  AddElement(fBacTrack);
  
  fNegTrack->SetLineColor(kCyan+2);  // in V0 = kCyan
  fNegTrack->SetStdTitle();
  fNegTrack->IncDenyDestroy();
  AddElement(fNegTrack);
  
  fPosTrack->SetLineColor(kRed+2);  // in Eve V0 = kRed
  fPosTrack->SetStdTitle();
  fPosTrack->IncDenyDestroy();
  AddElement(fPosTrack);  
}

//______________________________________________________________________________
AliEveCascade::~AliEveCascade()
{
  // Destructor. Dereferences bachelor track and pointing-line objects.

  fBacTrack->DecDenyDestroy();
  fNegTrack->DecDenyDestroy();
  fPosTrack->DecDenyDestroy();
  fPointingCurve->DecDenyDestroy();
  fV0Path->DecDenyDestroy();
}

//______________________________________________________________________________
void AliEveCascade::MakeCascade()
{
  // Set all dependant components for drawing.

  SetPoint(0, fRecDecayV.fX, fRecDecayV.fY, fRecDecayV.fZ);

  fBacTrack->MakeTrack();  
  fNegTrack->MakeTrack();
  fPosTrack->MakeTrack();

  fPointingCurve->SetPoint(0, fRecBirthV.fX, fRecBirthV.fY, fRecBirthV.fZ);
  fPointingCurve->SetPoint(1, fRecDecayV.fX, fRecDecayV.fY, fRecDecayV.fZ);

  fV0Path->SetPoint(0, fRecDecayV.fX, fRecDecayV.fY, fRecDecayV.fZ);
  fV0Path->SetPoint(1, fRecDecayV0.fX, fRecDecayV0.fY, fRecDecayV0.fZ);
}



//______________________________________________________________________________
Double_t AliEveCascade::GetInvMass(Int_t cascadePdgCodeHyp) const
{
  // Returns Invariant Mass assuming the masses of the daughter particles
	TEveVector lBachMomentum = fBacTrack->GetMomentum();
  // Does not work properly because momenta at the primary vertex for V0 daughters (same for AliEVEv0) + Bach
	
	Double_t lBachMass = 0.0;
	if(cascadePdgCodeHyp == kXiMinus || cascadePdgCodeHyp == -kXiMinus )
		lBachMass = TDatabasePDG::Instance()->GetParticle("pi-")->Mass();
	else if(cascadePdgCodeHyp == kOmegaMinus || cascadePdgCodeHyp == -kOmegaMinus )
		lBachMass = TDatabasePDG::Instance()->GetParticle("K-")->Mass();
	
	Double_t lLambdaMass = TDatabasePDG::Instance()->GetParticle("Lambda0")->Mass();
	
	printf("\n Cascade : check the mass of the bachelor (%.5f) - the Lambda (%.5f)\n",lBachMass, lLambdaMass);
	
	Double_t eBach   = TMath::Sqrt(lBachMass*lBachMass     + lBachMomentum.Mag2());
	//Double_t eBach   = TMath::Sqrt(lBachMass  *lBachMass   + fBachP.Mag2());
	Double_t eLambda = TMath::Sqrt(lLambdaMass*lLambdaMass + fLambdaP.Mag2());
	
	TVector3 lCascadeMom(fBachP + fLambdaP);
	/*
	printf("Ptot(Lambda) : %f GeV/c \n", fLambdaP.Mag());
	printf("Ptot(Bach)   : %f GeV/c \n", lBachMomentum.Mag());
	printf("Ptot(BachESD): %f GeV/c \n", fBachP.Mag());
	printf("Ptot(Casc)   : %f GeV/c \n", fRecDecayP.Mag());
	printf("Ptot(CascESD): %f GeV/c \n", lCascadeMom.Mag());
	*/
	//return TMath::Sqrt( (eBach+eLambda)*(eBach+eLambda) - fRecDecayP.Mag2() );	// numerically instable ...
	//return TMath::Sqrt( ((eBach+eLambda) - lCascadeMom.Mag()) * ((eBach+eLambda) + lCascadeMom.Mag()) );
	return TMath::Sqrt( ((eBach+eLambda) - fRecDecayP.Mag()) * ((eBach+eLambda) + fRecDecayP.Mag()) );
}



//------------------------------------------------------------------------
//
// AliEveCascadeList class
//
//------------------------------------------------------------------------

ClassImp(AliEveCascadeList)

//______________________________________________________________________________
AliEveCascadeList::AliEveCascadeList() :
  TEveElementList(),
  fTitle(),
  fRnrStyleBac(0),
  fRnrStyleNeg(0),
  fRnrStylePos(0),
  fRnrDaughters(kTRUE),
  fRnrCascadevtx(kTRUE),
  fRnrCascadepath(kTRUE),
  fBacColor(0),
  fMinRCut(0),
  fMaxRCut(250),
  fMinDaughterDCA(0),
  fMaxDaughterDCA(1),
  fMinPt(0),
  fMaxPt(20),
  fInvMassHyp(kXiMinus),
  fMinInvariantMass(1.0),
  fMaxInvariantMass(6.0)
{
  // Default constructor.

  fChildClass = AliEveCascade::Class(); // override member from base TEveElementList
}

//______________________________________________________________________________
AliEveCascadeList::AliEveCascadeList(TEveTrackPropagator* rsBac,TEveTrackPropagator* rsNeg,TEveTrackPropagator* rsPos) :
  TEveElementList(),
  fTitle(),
  fRnrStyleBac(rsBac),
  fRnrStyleNeg(rsNeg),
  fRnrStylePos(rsPos),
  fRnrDaughters(kTRUE),
  fRnrCascadevtx(kTRUE),
  fRnrCascadepath(kTRUE),
  fBacColor(0),
  fMinRCut(0),
  fMaxRCut(250),
  fMinDaughterDCA(0),
  fMaxDaughterDCA(1),
  fMinPt(0),
  fMaxPt(20),
  fInvMassHyp(kXiMinus),
  fMinInvariantMass(1.0),
  fMaxInvariantMass(6.0)
{
  // Constructor with given track-propagator..

  fChildClass = AliEveCascade::Class(); // override member from base TEveElementList

  Init();
}

//______________________________________________________________________________
AliEveCascadeList::AliEveCascadeList(const Text_t* name, TEveTrackPropagator* rsBac, TEveTrackPropagator* rsNeg, TEveTrackPropagator* rsPos) :
  TEveElementList(),
  fTitle(),
  fRnrStyleBac(rsBac),
  fRnrStyleNeg(rsNeg),
  fRnrStylePos(rsPos),
  fRnrDaughters(kTRUE),
  fRnrCascadevtx(kTRUE),
  fRnrCascadepath(kTRUE),
  fBacColor(0),
  fMinRCut(0),
  fMaxRCut(100),
  fMinDaughterDCA(0),
  fMaxDaughterDCA(1),
  fMinPt(0),
  fMaxPt(20),
  fInvMassHyp(kXiMinus),
  fMinInvariantMass(1.0),
  fMaxInvariantMass(6.0)
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

  if (fRnrStyleBac== 0) fRnrStyleBac = new TEveTrackPropagator;
  if (fRnrStyleNeg== 0) fRnrStyleNeg = new TEveTrackPropagator;
  if (fRnrStylePos== 0) fRnrStylePos = new TEveTrackPropagator;
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

/******************************************************************************/

//______________________________________________________________________________
void AliEveCascadeList::FilterByInvariantMass(Float_t minInvariantMass, Float_t maxInvariantMass, Int_t cascadePdgCodeHyp)
{
  // Select visibility of elements based on the Cascade invariant mass.

	fMinInvariantMass = minInvariantMass;
	fMaxInvariantMass = maxInvariantMass;

	for(List_i i = fChildren.begin(); i != fChildren.end(); ++i)
	{
		AliEveCascade* cascade = (AliEveCascade*) *i;
		Float_t  invMass = (Float_t)cascade->GetInvMass(cascadePdgCodeHyp);
		Bool_t   show    = invMass >= fMinInvariantMass && invMass <= fMaxInvariantMass;
		cascade->SetRnrState(show);
	}
	ElementChanged();
	gEve->Redraw3D();
}

