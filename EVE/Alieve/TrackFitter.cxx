// $Header$

#include "TrackFitter.h"

//#include "TClass.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "AliCluster.h"
#include "AliRieman.h"
#include "AliExternalTrackParam.h"

#include <Reve/Track.h>
#include <Reve/PODs.h>
#include <Reve/ReveManager.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// TrackFitter
//
//  TrackFitter is an interface to helix fit. It creates a set of points, listening to signal 
//  PointCtrlClicked() of any Reve::PointSet. Via editor it fits selected points and creates a 
//  reconstructed track.
// 

ClassImp(TrackFitter)

TrackFitter::TrackFitter(const Text_t* name, Int_t n_points, TreeVarType_e tv_type) :
    Reve::PointSet(name, n_points, tv_type),

    fGraphSelected(0),
    fGraphFitted(0),

    fRieman(0),
    fConnected(kFALSE),
    fTrackList(0)
{
  // Constructor.

  SetMarkerColor(3);
  SetOwnIds(kFALSE);

  fGraphSelected = new TGraph();
  fGraphSelected->SetName("Selected points");
  fGraphSelected->SetMarkerColor(4);
  fGraphSelected->SetMarkerStyle(4);  
  fGraphSelected->SetMarkerSize(2);

  fGraphFitted = new TGraphErrors();
  fGraphFitted->SetName("Fitted points");
  fGraphFitted->SetMarkerColor(2);

  fTrackList = new TrackList("Tracks");
  fTrackList->SetLineWidth(2);
  fTrackList->SetLineColor(8);
  fTrackList->IncDenyDestroy();
  fTrackList->GetRnrStyle()->SetEditPathMarks(kTRUE);
  gReve->AddRenderElement(fTrackList, this);
  UpdateItems();
}

TrackFitter::~TrackFitter()
{
  // Destructor.

  if(fRieman) delete fRieman;
  fTrackList->DecDenyDestroy();
  delete fTrackList;
}

/**************************************************************************/
void TrackFitter::DestroyElements()
{
  // Virtual method of base class Reve::RenderElement.
  // It preserves track list to have coomon track propagator attributes.

  RenderElement::DestroyElements();
  gReve->AddRenderElement(fTrackList, this);
  fTrackList->DestroyElements();
  UpdateItems();
}

/**************************************************************************/
void TrackFitter::Start()
{
  // Start selection of points.

  Reset();
  if(fConnected == kFALSE)
  {
    TQObject::Connect("Reve::PointSet", "PointCtrlClicked(Reve::PointSet*,Int_t)",
		      "Alieve::TrackFitter", this, "AddFitPoint(Reve::PointSet*,Int_t)");

    fConnected = kTRUE;
  }
}

void TrackFitter::Stop()
{
  // Stop selection of points.

  if(fConnected)
  {
    TQObject::Disconnect("Reve::PointSet", "AddFitPoint(Reve::PointSet*,Int_t)");
    fConnected = kFALSE;
  }
}

/**************************************************************************/

void TrackFitter::AddFitPoint(Reve::PointSet* ps, Int_t n)
{ 
  // Add/remove given point depending if exists in the fMapPS.
 
  Float_t x, y, z;

  std::map<Point_t, Int_t>::iterator g = fMapPS.find(Point_t(ps, n));
  if(g != fMapPS.end())
  {
    Int_t idx = g->second;
    if(idx != fLastPoint)
    {
      GetPoint(fLastPoint, x, y, z);
      SetPoint(idx, x, y, z);
    }
    fMapPS.erase(g);
    fLastPoint--;
  }
  else 
  {
    fMapPS[Point_t(ps, n)] = Size();
    ps->GetPoint(n, x, y, z);
    SetNextPoint(x, y, z); 
    SetPointId(ps->GetPointId(n));
  }
  ResetBBox();
  ElementChanged(kTRUE, kTRUE);
}

/**************************************************************************/

void TrackFitter::FitTrack()
{
  // Fit selected points with AliRieman fitter.

  using namespace TMath;

  if(fRieman) delete fRieman;
  fRieman = new AliRieman(Size());

  Float_t x, y, z;
  Int_t alphaIdx = 0;
  GetPoint(alphaIdx, x, y, z);
  Float_t minR2=x*x + y*y;
  for (Int_t i=1; i<=fLastPoint; i++)
  {
    GetPoint(i, x, y, z);
    Float_t cR2 = x*x + y*y;
    if ( minR2 > cR2 )
    {
      minR2 = cR2;
      alphaIdx = i;
    }
  }
  GetPoint(alphaIdx, x, y, z);
  fAlpha = ATan2(y, x);
  Float_t sin = Sin(-fAlpha);
  Float_t cos = Cos(-fAlpha);  
  for (Int_t i=0; i<=fLastPoint; i++) { 
    GetPoint(i, x, y, z);
    fRieman->AddPoint(cos*x - sin*y, cos*y + sin*x, z, 1, 1);
  }
  fRieman->Update();

  Double_t r = Sqrt(minR2);
  Double_t param[5];
  Double_t cov[15];
  fRieman->GetExternalParameters(r, param, cov);
  // curvature to pt
  param[4] /= TrackRnrStyle::fgDefMagField*TrackRnrStyle::fgkB2C;
  // sign in tang
  if(param[4] < 0) param[3] *= -1;
  AliExternalTrackParam trackParam(r, fAlpha, param, cov);
  trackParam.Print();

  // make track
  Double_t V0[3];
  trackParam.GetXYZAt(r, TrackRnrStyle::fgDefMagField, V0);
  Double_t P0[3];
  trackParam.GetPxPyPzAt(r, TrackRnrStyle::fgDefMagField, P0);
  RecTrack rc;
  rc.V.Set(V0); 
  rc.P.Set(P0);
  rc.sign = trackParam.Charge();

  Track* track = new Track(&rc, fTrackList->GetRnrStyle());
  track->SetName(Form("track %f", fAlpha));
  PathMark* pm = new PathMark(PathMark::Daughter);
  for(Int_t i=0; i==fLastPoint; i++)
  {
    GetPoint(i, x, y, z);
    pm->V.Set(x, y, z);
    pm->P.Set(P0); 
    track->AddPathMark(pm);
  }
  track->MakeTrack();
  track->SetAttLineAttMarker(fTrackList); 
  gReve->AddRenderElement(track, fTrackList);
}


void TrackFitter::Reset(Int_t n, Int_t ids)
{
  // Reset selection.

  if(fRieman) fRieman->Reset();
  PointSet::Reset(n, ids);
  fMapPS.clear();
}

/**************************************************************************/
void TrackFitter::DrawRiemanGraph()
{
  // Draw graph of rieman fit.

   static const Exc_t eH("TrackFitter::DrawRiemanGraph ");

  if(fRieman == 0)
    throw(eH + "fitter not set.");

  Int_t nR = fRieman->GetN();
  fGraphSelected->Set(nR);
  fGraphFitted->Set(nR);

  Double_t* x =  fRieman->GetX();  
  Double_t* y =  fRieman->GetY();
  Double_t* sy =  fRieman->GetSy();
  for (Int_t i=0; i<nR; i++)
  {
    fGraphSelected->SetPoint(i, x[i], y[i]);
    fGraphFitted->SetPoint(i, x[i], fRieman->GetYat(x[i]));
    fGraphFitted->SetPointError(i, 0.1, sy[i]);
  }
  
  if (gPad) gPad->Clear();
  fGraphSelected->Draw("AP");
  fGraphFitted->Draw("SAME P");
  gPad->GetCanvas()->SetTitle(Form("AliRieman alpha: %f", fAlpha));
  gPad->Modified();
  gPad->Update();
}
