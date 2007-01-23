/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id$ */
/** @file    AliFMDPattern.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:39:09 2006
    @brief   FMD 2D Event display 
*/
//___________________________________________________________________
//
// The classes defined here, are utility classes for reading in data
// for the FMD.  They are  put in a seperate library to not polute the
// normal libraries.  The classes are intended to be used as base
// classes for customized class that do some sort of analysis on the
// various types of data produced by the FMD. 
//
// Latest changes by Christian Holm Christensen
//

#include <iostream>

#include <TApplication.h>
#include <TButton.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TMath.h>
#include <TPad.h>
#include <TRandom.h>
#include <TSlider.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVector2.h>
#include <TView.h>

#include "AliFMDPattern.h"	// ALIFMDDISPLAY_H
#include "AliFMDGeometry.h"	// ALIFMDGEOMETRY_H
#include "AliFMDParameters.h"	// ALIFMDPARAMETERS_H
#include "AliFMDRing.h"
#include "AliFMDDetector.h"
#include "AliFMDHit.h"
#include <AliLog.h>

//____________________________________________________________________
ClassImp(AliFMDPattern)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDPattern::Detector::Detector(UShort_t id) 
  : fId(id),
    fCounts(0), 
    fGraphs(0), 
    fFrame(0)
{}

//____________________________________________________________________
AliFMDPattern::Detector::~Detector()
{
  if (fFrame) delete fFrame;
}

//____________________________________________________________________
void
AliFMDPattern::Detector::DrawShape(TObjArray& a) 
{
  TIter next(&a);
  TGraph* g = 0;
  while ((g = static_cast<TGraph*>(next()))) {
    g->DrawClone("f same");
    g->DrawClone("l same");
  }
}

//____________________________________________________________________
void
AliFMDPattern::Detector::Begin(Int_t nlevel, Double_t r, 
			       TObjArray& inners, TObjArray& outers)
{
  if (nlevel < 1) nlevel = gStyle->GetNumberOfColors();
  fCounts.Set(nlevel);
  if (!fFrame) {
    fFrame = new TH2F(Form("fmd%dFrame", fId), Form("FMD%d", fId), 
		      10, -r, r, 10, -r, r);
    fFrame->SetStats(kFALSE);
    fFrame->Draw();
  }
  DrawShape(inners);
  if (fId != 1) DrawShape(outers);
  for (Int_t i = 0; i < nlevel; i++) { 
    TGraph* g = new TGraph;
    Int_t idx = Int_t(Float_t(i) / nlevel * gStyle->GetNumberOfColors());
    Int_t col = gStyle->GetColorPalette(idx);
    g->SetName(Form("FMD%d_L%02d", fId, i));
    g->SetMarkerColor(col);
    g->SetLineColor(col);
    g->SetFillColor(col);
    g->SetMarkerSize(i * .2 + .2);
    g->SetMarkerStyle(2);
    g->Draw("same p");
    fGraphs.AddAtAndExpand(g, i);
  }
  TIter   next(&fGraphs);
}

//____________________________________________________________________
void
AliFMDPattern::Detector::Clear() 
{
  fCounts.Reset(0);
}

//____________________________________________________________________
void
AliFMDPattern::Detector::End()
{
  TIter   next(&fGraphs);
  TGraph* g = 0;
  Int_t   i = 0;
  while ((g = static_cast<TGraph*>(next()))) g->Set(fCounts[i++]);
}
//____________________________________________________________________
void
AliFMDPattern::Detector::AddMarker(Double_t x, Double_t y, Float_t s, 
				   Float_t max)
{
  Int_t i = TMath::Min(Int_t(fCounts.fN * s / max), 
		       Int_t(fGraphs.GetEntries()-1));
  TGraph* g = static_cast<TGraph*>(fGraphs.At(i));
  if (!g) return;
  g->SetPoint(fCounts[i]++, x, y);
}


//____________________________________________________________________
AliFMDPattern::AliFMDPattern(const char* gAliceFile)
  : AliFMDDisplay(kTRUE, gAliceFile),
    fInnerMax(0), 
    fOuterMax(0),
    fFMD1Pad(0),
    fFMD1(1),
    fFMD2Pad(0),
    fFMD2(2),
    fFMD3Pad(0),
    fFMD3(3),
    fEvent(.1, .8, "Event #"),
    fFMD1Sum(.2, .7, "# in FMD1: "),
    fFMD2Sum(.2, .6, "# in FMD2: "),
    fFMD3Sum(.2, .5, "# in FMD3: "),
    fLine(.15, .47, .85, .47),
    fTotal(.2, .35, "Total:   ")
{
  // RemoveLoad(kGeometry);
  fEvent.SetBit(TLatex::kTextNDC);
  fFMD1Sum.SetBit(TLatex::kTextNDC);
  fFMD2Sum.SetBit(TLatex::kTextNDC);
  fFMD3Sum.SetBit(TLatex::kTextNDC);
  fLine.SetBit(TLine::kLineNDC);
  fTotal.SetBit(TLatex::kTextNDC);
}

//____________________________________________________________________
AliFMDPattern::~AliFMDPattern()
{
  fInners.Delete();
  fOuters.Delete();
}

    
//____________________________________________________________________
Bool_t 
AliFMDPattern::Init()
{
  // Initialize.  GEt transforms and such, 
  if (!AliFMDInput::Init()) return kFALSE;
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  geom->Init();
  geom->InitTransformations();
  
  Char_t rs[] = { 'I' , 'O', '\0' };
  Char_t *r   = rs;
  do {
    AliFMDRing* ring = geom->GetRing(*r);
    if (!ring) continue;
    const TObjArray& vs = ring->GetVerticies();
    TObjArray&       gs = (*r == 'I' ? fInners   : fOuters);
    Float_t&         mr = (*r == 'I' ? fInnerMax : fOuterMax);
    Int_t            nm = ring->GetNModules();
    AliInfo(Form("Making %d modules for %c", nm, *r));
    for (Int_t m = 0; m < nm; m++) {
      Int_t          nv = vs.GetEntries();
      Double_t       a  = TMath::Pi() / 180 * (m * 2 + 1) * ring->GetTheta();
      TGraph*        g  = new TGraph(nv+1);
      Double_t       x0 = 0, y0 = 0;
      gs.AddAtAndExpand(g, m);
      for (Int_t c = 0; c < nv; c++) {
	TVector2* v = static_cast<TVector2*>(vs.At(c));
	mr          = TMath::Max(mr, Float_t(v->Mod()));
	TVector2  w(v->Rotate(a));
	if (c == 0) { x0 = w.X(); y0 = w.Y(); }
	g->SetPoint(c, w.X(), w.Y());
      }
      g->SetName(Form("FMDX%c_%02d", *r, m));
      g->SetPoint(nv, x0, y0);
      g->SetFillColor((*rs == 'I' ? 
		       (m % 2 == 0 ? 18 : 17) :
		       (m % 2 == 0 ? 20 : 23)));
      g->SetFillStyle(3001);
      g->SetLineColor(1);
      g->SetLineWidth(1);
      g->SetLineStyle(2);
    }
  } while (*(++r));
    
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDPattern::Begin(Int_t event) 
{
  MakeAux();
  if (!fCanvas) {
    const char* which[] = { "Continue", "Redisplay", 0 };
    MakeCanvas(which);
    
    AliFMDGeometry* geom = AliFMDGeometry::Instance();
    AliFMDDetector* det;
    if ((det = geom->GetDetector(1))) {
      fPad->cd();
      fFMD1Pad = new TPad("FMD1", "FMD1", 0.0, 0.50, 0.5, 1.0, 0, 0);
      fFMD1Pad->Draw();
      fFMD1Pad->cd();
      fFMD1.Begin(-1, fInnerMax, fInners, fOuters);
    }
    if ((det = geom->GetDetector(2))) {
      fPad->cd();
      fFMD2Pad = new TPad("FMD2", "FMD2", 0.5, 0.50, 1.0, 1.0, 0, 0);
      fFMD2Pad->Draw();
      fFMD2Pad->cd();
      fFMD2.Begin(-1, fOuterMax, fInners, fOuters);
    }
    if ((det = geom->GetDetector(3))) {
      fPad->cd();
      fFMD3Pad = new TPad("FMD3", "FMD3", 0.0, 0.0, .5, .5, 0, 0);
      fFMD3Pad->Draw();
      fFMD3Pad->cd();
      fFMD3.Begin(-1, fOuterMax, fInners, fOuters);
    }
    fPad->cd();
    fSummary = new TPad("display", "Display", 0.5, 0.0, 1.0, 0.5, 0, 0);
    fSummary->Draw();
    fSummary->cd();
    fEvent.Draw();
    fFMD1Sum.Draw();
    fFMD2Sum.Draw();
    fFMD3Sum.Draw();
    fLine.Draw();
    fTotal.Draw();
  }
  fEvent.SetTitle(Form("Event # %6d", event));

  fCanvas->Modified();
  fCanvas->Update();
  fCanvas->cd();
  fFMD1.Clear();
  fFMD2.Clear();
  fFMD3.Clear();
  return AliFMDInput::Begin(event);
}

//____________________________________________________________________
void
AliFMDPattern::Redisplay()
{
  fFMD1.Clear();
  fFMD2.Clear();
  fFMD3.Clear();
  AliFMDDisplay::Redisplay();
}

//____________________________________________________________________
void
AliFMDPattern::AtEnd()
{
  DrawAux();
  
  Int_t total = 0;
  
  fFMD1.End();
  fFMD1Pad->Modified();
  fFMD1Sum.SetTitle(Form("# hits in FMD1: %5d", fFMD1.Total()));
  total += fFMD1.Total();

  fFMD2.End();
  fFMD2Pad->Modified();
  fFMD2Sum.SetTitle(Form("# hits in FMD2: %5d", fFMD2.Total()));
  total += fFMD2.Total();

  fFMD3.End();
  fFMD3Pad->Modified();
  fFMD3Sum.SetTitle(Form("# hits in FMD3: %5d", fFMD3.Total()));
  total += fFMD3.Total();

  fTotal.SetTitle(Form("Total:    %5d/51200 (%3d%%)", 
		      total, Int_t(100. / 51200 * total)));
  fSummary->Modified();
  fCanvas->Modified();
  fCanvas->Update();
  fCanvas->cd();
}

//____________________________________________________________________
Bool_t 
AliFMDPattern::ProcessHit(AliFMDHit* hit, TParticle*) 
{
  switch (hit->Detector()) {
  case 1: fFMD1.AddMarker(hit->X(), hit->Y(), hit->Edep(), 1); break;
  case 2: fFMD2.AddMarker(hit->X(), hit->Y(), hit->Edep(), 1); break;
  case 3: fFMD3.AddMarker(hit->X(), hit->Y(), hit->Edep(), 1); break;
  }
  return kTRUE;
}


//____________________________________________________________________
void
AliFMDPattern::AddMarker(UShort_t det, Char_t rng, UShort_t sec, UShort_t str,
			 TObject*, Float_t s, Float_t max)
{
  // Add a marker to the display
  //
  //    det 	Detector
  //    rng 	Ring
  //    sec 	Sector 
  //    str 	Strip
  //    o   	Object to refer to
  //    s   	Signal 
  //    max 	Maximum of signal 
  //
  Detector* d = 0;
  switch (det) {
  case 1: d = &fFMD1; break;
  case 2: d = &fFMD2; break;
  case 3: d = &fFMD3; break;
  }
  if (!d) return;
  AliFMDGeometry*   geom = AliFMDGeometry::Instance();
  Double_t x, y, z;
  geom->Detector2XYZ(det, rng, sec, str, x, y, z);
  if (true) {
    AliFMDRing* r  = geom->GetRing(rng);
    Double_t    t  = .9 * r->GetTheta() / 2;
    Double_t    a  = gRandom->Uniform(-t,t) * TMath::Pi() / 180;
    Double_t    x1 = x * TMath::Cos(a) - y * TMath::Sin(a);
    Double_t    y1 = x * TMath::Sin(a) + y * TMath::Cos(a);
    x = x1;
    y = y1;
  }
  d->AddMarker(x, y, s, max);
}

//____________________________________________________________________
//
// EOF
//
