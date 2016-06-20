// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTrackFitter.h"
#include "AliEveTrack.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TQObject.h"

#include "AliRieman.h"
#include "AliExternalTrackParam.h"

#include <TEveTrackPropagator.h>
#include <TEveVSDStructs.h>
#include <TEveManager.h>


//==============================================================================
//==============================================================================
// AliEveTrackFitter
//==============================================================================

//______________________________________________________________________________
//
// AliEveTrackFitter is an interface to TEvePointSet allowing AliRieman fit.
// It builds a list of points by listening to selection signal of any
// object of type TEvePointSet. After selection the list is feeded to
// AliRieman fitter, which returns helix parameters visualized with
// AliEveTrack.
//

ClassImp(AliEveTrackFitter)

AliEveTrackFitter::AliEveTrackFitter(const Text_t* name, Int_t nPoints) :
    TEvePointSet   (name, nPoints),

    fAlpha         (0),
    fRieman        (0),

    fConnected     (kFALSE),
    fSPMap         (),
    fTrackList     (0),

    fGraphPicked   (0),
    fGraphHelix    (0)
{
  // Constructor.

  SetMarkerColor(3);
  SetOwnIds(kFALSE);

  fTrackList = new TEveTrackList("Tracks");
  fTrackList->IncDenyDestroy();
  fTrackList->SetLineWidth(2);
  fTrackList->SetLineColor(8);
  fTrackList->GetPropagator()->SetEditPathMarks(kTRUE);
  AddElement(fTrackList);

  fGraphPicked = new TGraph();
  fGraphPicked->SetName("Selected points");
  fGraphPicked->SetMarkerColor(4);
  fGraphPicked->SetMarkerStyle(4);
  fGraphPicked->SetMarkerSize(2);

  fGraphHelix = new TGraphErrors();
  fGraphHelix->SetName("Fitted points");
  fGraphHelix->SetMarkerColor(2);
}

AliEveTrackFitter::~AliEveTrackFitter()
{
  // Destructor.

  if (fRieman) delete fRieman;

  fTrackList->DecDenyDestroy();
}

/******************************************************************************/

void AliEveTrackFitter::DestroyElements()
{
  // Virtual method of base class TEveElement.
  // Preserves TEveTrackList object for fitted helices.

  TEveElement::DestroyElements();

  // fTrackList is destroyed because DenyDestroy is set.
  gEve->AddElement(fTrackList, this);
  fTrackList->DestroyElements();
}

/******************************************************************************/

void AliEveTrackFitter::Start()
{
  // Clear existing point selection and maintain connection to the
  // TEvePointSet signal.

  Reset();
  if (fConnected == kFALSE)
  {
    TQObject::Connect("TEvePointSet", "PointSelected(Int_t)",
		      "AliEveTrackFitter", this, "AddFitPoint(Int_t)");
    fConnected = kTRUE;
  }
}

void AliEveTrackFitter::Stop()
{
  // Stop adding points for the fit.

  if (fConnected)
  {
    TQObject::Disconnect("TEvePointSet", "AddFitPoint(Int_t)");
    fConnected = kFALSE;
  }
}

void AliEveTrackFitter::Reset(Int_t nPoints, Int_t nIntIds)
{
  // Reset selection.

  if (fRieman) fRieman->Reset();
  TEvePointSet::Reset(nPoints, nIntIds);
  fSPMap.clear();
}

/******************************************************************************/

void AliEveTrackFitter::AddFitPoint(Int_t pointId)
{
  // Add or remove given point depending if exists in the map.

  Float_t x, y, z;

  TEvePointSet* ps = static_cast<TEvePointSet*>((TQObject*) gTQSender);

  PointMap_t::iterator g = fSPMap.find(Point_t(ps, pointId));
  if (g != fSPMap.end())
  {
    Int_t idx = g->second;
    if (idx != fLastPoint)
    {
      GetPoint(fLastPoint, x, y, z);
      SetPoint(idx, x, y, z);
    }
    fSPMap.erase(g);
    --fLastPoint;
  }
  else
  {
    fSPMap[Point_t(ps, pointId)] = Size();
    ps->GetPoint(pointId, x, y, z);
    SetNextPoint(x, y, z);
  }

  ResetBBox();
  ElementChanged(kTRUE, kTRUE);
}

/******************************************************************************/

void AliEveTrackFitter::FitTrack()
{
  // Fit selected points with AliRieman fitter.

  using namespace TMath;

  if (fRieman) delete fRieman;
  fRieman = new AliRieman(Size());

  Float_t x, y, z;
  Int_t alphaIdx = 0;
  GetPoint(alphaIdx, x, y, z);
  Float_t minR2=x*x + y*y;
  for (Int_t i=0; i<=fLastPoint; i++)
  {
    GetPoint(i, x, y, z);
    Float_t cR2 = x*x + y*y;
    if (minR2 > cR2)
    {
      minR2 = cR2;
      alphaIdx = i;
    }
  }
  GetPoint(alphaIdx, x, y, z);
  fAlpha = ATan2(y, x);
  Float_t sin = Sin(-fAlpha);
  Float_t cos = Cos(-fAlpha);
  for (Int_t i = 0; i <= fLastPoint; ++i)
  {
    GetPoint(i, x, y, z);
    fRieman->AddPoint(cos*x - sin*y, cos*y + sin*x, z, 1, 1);
  }
  fRieman->Update();

  Double_t r = Sqrt(minR2);
  Double_t param[5];
  Double_t cov[15];
  fRieman->GetExternalParameters(r, param, cov);
  // curvature to pt
  param[4] /= TEveTrackPropagator::fgDefMagField*TEveTrackPropagator::fgkB2C;
  // sign in tang
  if (param[4] < 0) param[3] = -param[3];
  AliExternalTrackParam trackParam(r, fAlpha, param, cov);
  trackParam.Print();

  // make track
  Double_t v0[3];
  trackParam.GetXYZAt(r, TEveTrackPropagator::fgDefMagField, v0);
  Double_t p0[3];
  trackParam.GetPxPyPzAt(r, TEveTrackPropagator::fgDefMagField, p0);
  TEveRecTrack rc;
  rc.fV.Set(v0);
  rc.fP.Set(p0);
  rc.fSign = trackParam.Charge();

  AliEveTrack* track = new AliEveTrack(&rc, fTrackList->GetPropagator());
  track->SetName(Form("track %f", fAlpha));

  track->MakeTrack();
  track->SetAttLineAttMarker(fTrackList);
  fTrackList->AddElement(track);
}


/******************************************************************************/

void AliEveTrackFitter::DrawDebugGraph()
{
  // Draw graph of picked points and helix points.

  static const TEveException kEH("AliEveTrackFitter::DrawRiemanGraph ");

  if (fRieman == 0)
    throw(kEH + "fitter not set.");

  Int_t nR = fRieman->GetN();
  fGraphPicked->Set(nR);
  fGraphHelix->Set(nR);

  Double_t* x  =  fRieman->GetX();
  Double_t* y  =  fRieman->GetY();
  Double_t* sy =  fRieman->GetSy();
  for (Int_t i=0; i<nR; i++)
  {
    fGraphPicked->SetPoint(i, x[i], y[i]);
    fGraphHelix->SetPoint (i, x[i], fRieman->GetYat(x[i]));
    fGraphHelix->SetPointError(i, 0.1, sy[i]); // now faked
  }

  if (gPad)
    gPad->Clear();

  fGraphPicked->Draw("AP");
  fGraphHelix->Draw("SAME P");
  gPad->GetCanvas()->SetTitle(Form("AliRieman alpha: %f", fAlpha));
  gPad->Modified();
  gPad->Update();
}
