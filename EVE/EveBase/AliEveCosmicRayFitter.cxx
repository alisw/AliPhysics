// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "Riostream.h"

#include "AliEveCosmicRayFitter.h"

#include "TCanvas.h"
#include "TPad.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"

#include "TEveElement.h"

#include "TVirtualFitter.h"
#include "TLinearFitter.h"

#include "TEveLine.h"

#include "AliLog.h"

#include <TEveManager.h>
#include <Math/Vector3D.h>

//______________________________________________________________________
// AliEveCosmicRayFitter
//
//  Same schema used in AliEveTrackFitter class.
//
//  AliEveCosmicRayFitter is a TEvePointSet interface to allow
//  straight line fit. It creates a set of points by listening to
//  selection signal of any TEvePointSet object. After selection, the
//  list is feeded to TVirtualFitter class, which returns straight
//  line parameters. The fit result is visualized with a TEveLine
//  object.
// 
//  Thanks to L. Moneta for the algorithm prepared to make a 3D
//  straight line fit.
//
//  Author: A. De Caro
//

Double_t Distance3D(Double_t x, Double_t y, Double_t z, Double_t *p);
void SumDistance3D(Int_t &, Double_t *, Double_t & sum, Double_t * par, Int_t );

ClassImp(AliEveCosmicRayFitter)

AliEveCosmicRayFitter::AliEveCosmicRayFitter(const Text_t* name, Int_t n_points) :
    TEvePointSet    (name, n_points),

    fLine3DFitter   (0),
    fLineFitter1    (0),
    fLineFitter2    (0),

    fConnected      (kFALSE),
    fSPMap          (),
    fTrackList      (0),

    fStartingPointX (390.),

    fGraphPicked1   (0),
    fGraphLinear1   (0),
    fGraphPicked2   (0),
    fGraphLinear2   (0),
    fGraphPicked3   (),
    fGraphLinear3   ()

{
  // Constructor.

  fLine3DFitter = TVirtualFitter::Fitter(0, 4);
  fLineFitter1 = new TLinearFitter(1, "hyp1");
  fLineFitter2 = new TLinearFitter(1, "hyp1");

  SetMarkerColor(3);
  SetOwnIds(kFALSE);

  fTrackList = new TEveElementList("Cosmic ray");
  fTrackList->SetTitle("muons");
  fTrackList->SetMainColor((Color_t)8);
  gEve->AddElement(fTrackList);
  UpdateItems();

  fGraphPicked1 = new TGraph();
  fGraphPicked1->SetName("Selected points");
  fGraphPicked1->SetMarkerColor(4);
  fGraphPicked1->SetMarkerStyle(4);  
  fGraphPicked1->SetMarkerSize(2);

  fGraphLinear1 = new TGraphErrors();
  fGraphLinear1->SetName("Fitted points");
  fGraphLinear1->SetMarkerColor(2);

  fGraphPicked2 = new TGraph();
  fGraphPicked2->SetName("Selected points");
  fGraphPicked2->SetMarkerColor(4);
  fGraphPicked2->SetMarkerStyle(4);  
  fGraphPicked2->SetMarkerSize(2);

  fGraphLinear2 = new TGraphErrors();
  fGraphLinear2->SetName("Fitted points");
  fGraphLinear2->SetMarkerColor(2);

  fGraphPicked3 = new TGraph2D();
  fGraphPicked3->SetName("Selected points");
  fGraphPicked3->SetMarkerColor(4);
  fGraphPicked3->SetMarkerStyle(4);  
  fGraphPicked3->SetMarkerSize(2);

  fGraphLinear3 = new TGraph2DErrors();
  fGraphLinear3->SetName("Fitted points");
  fGraphLinear3->SetMarkerColor(2);

}

AliEveCosmicRayFitter::~AliEveCosmicRayFitter()
{
  // Destructor.

  if(fLine3DFitter) delete fLine3DFitter;
  if(fLineFitter1) delete fLineFitter1;
  if(fLineFitter2) delete fLineFitter2;

  fTrackList->DecDenyDestroy();
  delete fTrackList;

}

/**************************************************************************/
void AliEveCosmicRayFitter::Start()
{
  // Clear existing point selection and maintain connection to the
  // TEvePointSet signal.

  Reset();
  if(fConnected == kFALSE)
  {
    TQObject::Connect("TEvePointSet", "PointSelected(Int_t)",
		      "AliEveCosmicRayFitter", this, "AddFitPoint(Int_t)");
    fConnected = kTRUE;
  }
}

/**************************************************************************/
void AliEveCosmicRayFitter::Stop()
{
  // Stop selection of points.

  if(fConnected)
  {
    TQObject::Disconnect("TEvePointSet", "AddFitPoint(Int_t)");
    fConnected = kFALSE;
  }
}

/**************************************************************************/

void AliEveCosmicRayFitter::Reset(Int_t n, Int_t ids)
{
  // Reset selection.

  if(fGraphPicked3) fGraphPicked3->Clear();
  if(fLine3DFitter) fLine3DFitter->Clear();
  if(fLineFitter1) fLineFitter1->Clear();
  if(fLineFitter2) fLineFitter2->Clear();
  TEvePointSet::Reset(n, ids);
  fSPMap.clear();
}

/**************************************************************************/

//void AliEveCosmicRayFitter::AddFitPoint(TEvePointSet* ps, Int_t n)
void AliEveCosmicRayFitter::AddFitPoint(Int_t n)
{ 
  // Add/remove given point depending if exists in the fSPMap.
 
  Float_t x, y, z;
  TEvePointSet* ps = dynamic_cast<TEvePointSet*>((TQObject*) gTQSender);

  std::map<Point_t, Int_t>::iterator g = fSPMap.find(Point_t(ps, n));
  if(g != fSPMap.end())
  {
    Int_t idx = g->second;
    if(idx != fLastPoint)
    {
      GetPoint(fLastPoint, x, y, z);
      SetPoint(idx, x, y, z);
    }
    fSPMap.erase(g);
    fLastPoint--;
  }
  else 
  {
    fSPMap[Point_t(ps, n)] = Size();
    ps->GetPoint(n, x, y, z);
    SetNextPoint(x, y, z); 
    //SetPointId(ps->GetPointId(n));
  }
  ResetBBox();
  ElementChanged(kTRUE, kTRUE);
}

/**************************************************************************/

void AliEveCosmicRayFitter::FitTrack()
{
  // Fit selected points with TLinearFitter fitter.

  using namespace TMath;

  if(fLine3DFitter) delete fLine3DFitter;
  if(fLineFitter1) delete fLineFitter1;
  if(fLineFitter2) delete fLineFitter2;
  fLine3DFitter = TVirtualFitter::Fitter(0, 4);
  fLineFitter1 = new TLinearFitter(1, "hyp1");
  fLineFitter2 = new TLinearFitter(1, "hyp1");

  AliDebug(2, Form(" Selected point number %3i   %3i", fLastPoint+1, Size()));

  Double_t *posXpoint = new Double_t[fLastPoint+1];
  Double_t *posYpoint = new Double_t[fLastPoint+1];
  Double_t *posZpoint = new Double_t[fLastPoint+1];
  Double_t *errYpoint = new Double_t[fLastPoint+1];
  Double_t *errZpoint = new Double_t[fLastPoint+1];

  Double_t x, y, z;

  for (Int_t i=0; i<=fLastPoint; i++) { 
    x=0., y=0., z=0.;
    GetPoint(i, x, y, z);
    posXpoint[i] = x;
    posYpoint[i] = y;
    posZpoint[i] = z;
    errYpoint[i] = 10.;//0.05*posYpoint[i];//arbitrary
    errZpoint[i] = 10.;//0.05*posZpoint[i];//arbitrary

    fGraphPicked3->SetPoint(i,x,y,z);
    //fGraphPicked3->SetPointError(i,errXpoint,errYpoint[i],errZpoint[i]);

    AliDebug(2, Form("  %f    %f    %f",
		     posXpoint[i], posYpoint[i], posZpoint[i]));
  }
  /*
  fLine3DFitter->SetObjectFit(fGraphPicked3);
  fLine3DFitter->SetFCN(SumDistance3D);
  Double_t arglist[10];
  arglist[0] = 3;
  fLine3DFitter->ExecuteCommand("SET PRINT",arglist,1);
  
  Double_t pStart[4] = {1,1,1,1};
  fLine3DFitter->SetParameter(0,"y0",pStart[0],0.01,0,0);
  fLine3DFitter->SetParameter(1,"Ay",pStart[1],0.01,0,0);
  fLine3DFitter->SetParameter(2,"z0",pStart[2],0.01,0,0);
  fLine3DFitter->SetParameter(3,"Az",pStart[3],0.01,0,0);
    
  arglist[0] = 1000; // number of function calls 
  arglist[1] = 0.001; // tolerance 
  fLine3DFitter->ExecuteCommand("MIGRAD",arglist,2);

  //if (minos) fLine3DFitter->ExecuteCommand("MINOS",arglist,0);
  Int_t nvpar,nparx; 
  Double_t amin,edm, errdef;
  fLine3DFitter->GetStats(amin,edm,errdef,nvpar,nparx);
  fLine3DFitter->PrintResults(1,amin);
  */
  TEveLine* line3D = new TEveLine();
  line3D->SetLineColor(4);
  line3D->SetLineWidth(3);
  line3D->SetLineStyle(2);

  AliDebug(2,Form(" fStartingPointX = %f\n",fStartingPointX));

  Double_t parFit3D[4];
  for (int i = 0; i <4; ++i)
    parFit3D[i] = fLine3DFitter->GetParameter(i);

  Float_t xCoor = 0.;
  xCoor = fStartingPointX;
  /*
  Float_t radicando = parFit3D[1]*parFit3D[1]*parFit3D[0]*parFit3D[0] -
    (parFit3D[1]*parFit3D[1] + 1.)*(parFit3D[0]*parFit3D[0] -
				    fStartingPointX*fStartingPointX);
  if (radicando>=0.) xCoor = (- parFit3D[1]*parFit3D[0] + TMath::Sqrt(radicando))/
		       (parFit3D[1]*parFit3D[1] + 1.);
		       else xCoor = fStartingPointX;
  */

  Float_t yCoorP = parFit3D[0] + xCoor*parFit3D[1];
  Float_t zCoorP = parFit3D[2] + xCoor*parFit3D[3];
  Float_t yCoorM = parFit3D[0] - xCoor*parFit3D[1];
  Float_t zCoorM = parFit3D[2] - xCoor*parFit3D[3];

  line3D->SetNextPoint( xCoor, yCoorP, zCoorP);
  line3D->SetNextPoint(-xCoor, yCoorM, zCoorM);

  //gEve->AddElement(line3D, fTrackList);

  fLineFitter1->AssignData(fLastPoint+1, 1, posXpoint, posYpoint, errYpoint);
  fLineFitter2->AssignData(fLastPoint+1, 1, posXpoint, posZpoint, errZpoint);

  Int_t fitResult1 = fLineFitter1->Eval();
  Int_t fitResult2 = fLineFitter2->Eval();
  if (fitResult1 || fitResult2) {
    AliInfo(Form(" linear fit result is not ok (%1i, %1i)", fitResult1, fitResult2));
    return;
  }

  // make track
  TVectorD params1;
  TVectorD params2;
  fLineFitter1->GetParameters(params1);
  fLineFitter2->GetParameters(params2);

  TEveLine* line3Dbis = new TEveLine();
  line3Dbis->SetLineColor(4);
  line3Dbis->SetLineWidth(5);
  line3Dbis->SetLineStyle(2);

  xCoor = fStartingPointX;
  yCoorP = params1(0) + xCoor*params1(1);
  zCoorP = params2(0) + xCoor*params2(1);
  yCoorM = params1(0) - xCoor*params1(1);
  zCoorM = params2(0) - xCoor*params2(1);

  line3Dbis->SetNextPoint( xCoor, yCoorP, zCoorP);
  line3Dbis->SetNextPoint(-xCoor, yCoorM, zCoorM);

  gEve->AddElement(line3Dbis, fTrackList);

  delete [] posXpoint;
  delete [] posYpoint;
  delete [] posZpoint;
  delete [] errYpoint;
  delete [] errZpoint;

}

/**************************************************************************/
void AliEveCosmicRayFitter::DestroyElements()
{
  // Virtual method of base class Reve::RenderElement.
  // It preserves track list to have coomon track propagator attributes.

  TEveElement::DestroyElements();
  gEve->AddElement(fTrackList, this);
  fTrackList->DestroyElements();
  UpdateItems();
}

/**************************************************************************/
void AliEveCosmicRayFitter::DrawDebugGraph()
{
  //
  // Draw graph of Line fit.
  //

  static const TEveException eH("CosmicRayFitter::DrawDebugGraph ");

  if(fLineFitter1 == 0 || fLineFitter2 == 0)
    throw(eH + "fitter not set.");

  AliDebug(2,Form("           %3i   %3i\n", fLastPoint+1, Size()));

  Int_t nR = fLastPoint+1;
  fGraphPicked1->Set(nR);
  fGraphLinear1->Set(nR);
  fGraphPicked2->Set(nR);
  fGraphLinear2->Set(nR);

  TVectorD params1;
  TVectorD params2;
  fLineFitter1->GetParameters(params1);
  fLineFitter2->GetParameters(params2);
  /*
  TVectorD sy;
  fLineFitter1->GetErrors(sy);
  TVectorD sz;
  fLineFitter2->GetErrors(sz);
  */
  Double_t *x = new Double_t[nR];
  Double_t *y = new Double_t[nR];
  Double_t *z = new Double_t[nR];
  Float_t xx=0., yy=0., zz=0.;
  for (Int_t i=0; i<nR; i++) {
    GetPoint(i, xx, yy, zz);
    x[i] = (Double_t)xx;
    y[i] = (Double_t)yy;
    z[i] = (Double_t)zz;
    AliDebug(2,Form("  %f      %f      %f\n", x[i], y[i], z[i]));
  }

  fGraphLinear3->Set(nR);
  Double_t parFit3D[4];
   for (int i = 0; i <4; ++i) 
      parFit3D[i] = fLine3DFitter->GetParameter(i);

  Double_t yCoor = 0.;
  Double_t zCoor = 0.;
  for (Int_t i=0; i<nR; i++)
  {
    yCoor = params1(0) + x[i]*params1(1);
    zCoor = params2(0) + x[i]*params2(1);
    fGraphPicked1->SetPoint(i, x[i], y[i]);
    fGraphLinear1->SetPoint(i, x[i], yCoor);
    fGraphLinear1->SetPointError(i, TMath::Abs(x[i])*0.001, TMath::Abs(yCoor)*0.05);
    fGraphPicked2->SetPoint(i, x[i], z[i]);
    fGraphLinear2->SetPoint(i, x[i], zCoor);
    fGraphLinear2->SetPointError(i, TMath::Abs(x[i])*0.001, TMath::Abs(zCoor)*0.05);

    fGraphLinear3->SetPoint(i, x[i],
			    parFit3D[0]+x[i]*parFit3D[1],
			    parFit3D[2]+x[i]*parFit3D[3]);
    fGraphLinear3->SetPointError(i, TMath::Abs(x[i])*0.001,
				 TMath::Abs(parFit3D[0]+x[i]*parFit3D[1])*0.05,
				 TMath::Abs(parFit3D[2]+x[i]*parFit3D[3])*0.05);
  }
  
  TCanvas * canvas = 0;
  if (gPad) gPad->Clear();
  else if (gPad==0 || gPad->GetCanvas()->IsEditable() == kFALSE) {
    canvas = new TCanvas("canvas", "CosmicRayFitter", 800, 800);
    canvas->Clear();
  }
  canvas->SetHighLightColor(2);
  canvas->Range(0,0,1,1);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);

  canvas->Divide(1,2);

  TPad *canvas_1 = new TPad("canvas_1", "canvas_1",0.01,0.01,0.49,0.49);
  canvas_1->Draw();
  canvas_1->cd();
  canvas_1->SetFillColor(0);
  canvas_1->SetBorderSize(2);
  canvas_1->SetFrameBorderMode(0);
  fGraphPicked1->Draw("AP");
  fGraphLinear1->Draw("SAME P");

  canvas_1->Modified();
  canvas->cd();

  TPad *canvas_2 = new TPad("canvas_2", "canvas_2",0.51,0.01,0.99,0.49);
  canvas_2->Draw();
  canvas_2->cd();
  canvas_2->SetFillColor(0);
  canvas_2->SetBorderSize(2);
  canvas_2->SetFrameBorderMode(0);
  fGraphPicked2->Draw("AP");
  fGraphLinear2->Draw("SAME P");

  canvas_2->Modified();
  canvas->cd();

  TPad *canvas_3 = new TPad("canvas_2", "canvas_2",0.01,0.51,0.99,0.99);
  canvas_3->Draw();
  canvas_3->cd();
  canvas_3->SetFillColor(0);
  canvas_3->SetBorderSize(2);
  canvas_3->SetFrameBorderMode(0);
  fGraphPicked3->Draw("AP");
  fGraphLinear3->Draw("SAME P");

  canvas_3->Modified();
  canvas->cd();

  canvas->Modified();
  canvas->cd();

  delete x;
  delete y;
  delete z;

}

/**************************************************************************/
void SumDistance3D(Int_t &, Double_t *, Double_t & sum, Double_t * par, Int_t )
{
  //
  // Square sum of distances
  //

  TGraph2D * gr =
    dynamic_cast<TGraph2D*>((TVirtualFitter::GetFitter())->GetObjectFit() );
  assert(gr != 0);

  Double_t * x = gr->GetX();
  Double_t * y = gr->GetY();
  Double_t * z = gr->GetZ();
  Int_t npoints = gr->GetN();
  sum = 0;
  Double_t d = 0;
  for (Int_t i  = 0; i < npoints; ++i) { 
    d = Distance3D(x[i],y[i],z[i],par);
    sum += d;
#ifdef DEBUG
    std::cout << " point " << i << "\t"
	      << x[i] << "\t" 
	      << y[i] << "\t" 
	      << z[i] << "\t" 
	      << std::sqrt(d) << std::endl; 
    std::cout << " Total sum3D = " << sum << std::endl;
#endif

  }

}

/**************************************************************************/
Double_t Distance3D(Double_t x, Double_t y, Double_t z, Double_t *p)
{
  //
  // distance line point is D = | (xp-x0) cross  ux |
  // where ux is direction of line and x0 is a point in the line (like t = 0) 
  //

  using namespace ROOT::Math;

  XYZVector xp(x,y,z); 
  XYZVector x0(0., p[0], p[2]); 
  XYZVector x1(1., p[0] + p[1], p[2] + p[3]); 
  XYZVector u = (x1-x0).Unit(); 
  Double_t d2 = ((xp-x0).Cross(u)) .Mag2(); 
  return d2; 
}

/**************************************************************************/
