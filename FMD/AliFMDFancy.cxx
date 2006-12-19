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
/** @file    AliFMDFancy.cxx
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
#include "AliFMDFancy.h"	// ALIFMDDISPLAY_H
#include "AliFMDGeometry.h"	// ALIFMDGEOMETRY_H
#include "AliFMDParameters.h"	// ALIFMDPARAMETERS_H
#include <AliLog.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TButton.h>
#include <TCanvas.h>
#include <TView.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TSystem.h>
#include <TVector2.h>
#include "AliFMDRing.h"
#include "AliFMDDetector.h"
#include "AliFMDHit.h"
#include <TRandom.h>
#include <TLatex.h>
#include <TLine.h>
#include <iostream>
#include <iomanip>

//____________________________________________________________________
ClassImp(AliFMDFancy)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDFancy::AliFMDFancy(const char* gAliceFile)
  : AliFMDDisplay(kTRUE, gAliceFile),
    fFMD1Pad(0),
    fFMD1(1),
    fFMD2Pad(0),
    fFMD2(2),
    fFMD3Pad(0),
    fFMD3(3),
    fEvent(.1, .8, "Event #"),
    fFMD1IHits(.2, .7, "# in FMD1I: "),
    fFMD2IHits(.2, .6, "# in FMD2I: "),
    fFMD2OHits(.2, .5, "# in FMD2O: "),
    fFMD3IHits(.2, .4, "# in FMD3I: "),
    fFMD3OHits(.2, .3, "# in FMD3O: "),
    fLine(.15, .27, .85, .27),
    fTotal(.2, .15, "Total:   ")
{
  fEvent.SetBit(TLatex::kTextNDC);
  fFMD1IHits.SetBit(TLatex::kTextNDC);
  fFMD2IHits.SetBit(TLatex::kTextNDC);
  fFMD2OHits.SetBit(TLatex::kTextNDC);
  fFMD3IHits.SetBit(TLatex::kTextNDC);
  fFMD3OHits.SetBit(TLatex::kTextNDC);
  fLine.SetBit(TLine::kLineNDC);
  fTotal.SetBit(TLatex::kTextNDC);
}

//____________________________________________________________________
AliFMDFancy::Detector::Detector(UShort_t id)
  : fId(id)
{
  fInnerHits.SetName(Form("FMD%dI", id));
  fInnerHits.SetMarkerStyle(1); // 20);
  fInnerHits.SetMarkerSize(.2);
  fInnerHits.SetMarkerColor(50); // 12);
  fOuterHits.SetName(Form("FMD%dO", id));
  fOuterHits.SetMarkerStyle(1); // 20);
  fOuterHits.SetMarkerSize(.2);
  fOuterHits.SetMarkerColor(50); // 12);
}

//____________________________________________________________________
AliFMDFancy::Detector::~Detector()
{
  fShapes.Delete();
  if (fFrame) delete fFrame;
}

//____________________________________________________________________
AliFMDFancy::~AliFMDFancy()
{}

//____________________________________________________________________
void
AliFMDFancy::Detector::AddHistogram(TGraph2D& g, const char* opt)
{
  TH2* h = g.GetHistogram(opt);
  if (!h) return;
  h->SetBins(1, -fMaxR, fMaxR, 1, -fMaxR, fMaxR);
  h->GetZaxis()->SetRangeUser(fMinZ, fMaxZ);
}

  

//____________________________________________________________________
void
AliFMDFancy::Detector::Init()
{
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  AliFMDDetector* det  = geom->GetDetector(fId);
  if (!det) return;
  Char_t   rs[]   = { 'I' , 'O', '\0' };
  Char_t*  rp     = rs;
  Char_t   r;
  Double_t maxR = 0;
  Double_t minZ = 10000;
  Double_t maxZ = -10000;
  Int_t    ns   = 0;
  while ((r = *(rp++))) {
    AliFMDRing* ring = det->GetRing(r);
    if (!ring) continue;
    // if (r == 'O') continue;
    const TObjArray& vs = ring->GetVerticies();
    Int_t            nm = ring->GetNModules();
    Double_t         zd = (r == 'I' ? det->GetInnerZ() : det->GetOuterZ());
    for (Int_t m = 0; m < nm; m++) {
      Int_t          nv = vs.GetEntries();
      Double_t       a  = TMath::Pi() / 180 * (m * 2 + 1) * ring->GetTheta();
      TGraph2D*      g  = new TGraph2D(nv);
      Double_t       x0 = 0, y0 = 0, z0 = 0;
      Double_t       z  = zd + (m % 2==0 ? 0 : 
				TMath::Sign(ring->GetModuleSpacing(), zd));
      minZ              = TMath::Min(minZ, z);
      maxZ              = TMath::Max(maxZ, z);
      g->SetName(Form("FMD%d%cM%02d", fId, r, m));
      fShapes.AddAtAndExpand(g, ns++);
      for (Int_t c = 0; c < nv; c++) {
	TVector2* v = static_cast<TVector2*>(vs.At(nv - 1 - c));
	TVector2  w(v->Rotate(a));
	if (c == 0) { x0 = w.X(); y0 = w.Y(); z0 = z; }
	g->SetPoint(c, w.X(), w.Y(), z);
	maxR        = TMath::Max(maxR, v->Mod());
      }
      //g->SetPoint(nv, x0, y0, z0);
      g->SetFillColor(2);
      g->SetFillStyle(3002);
      g->SetLineColor(2);
      g->SetLineWidth(1);
    }
  }
  fMaxR = 1.05 * maxR;
  fMinZ = (minZ  > 0 ? 0.95 * minZ : 1.05 * minZ);
  fMaxZ = (maxZ  > 0 ? 1.05 * maxZ : 0.95 * maxZ);

  TIter next(&fShapes);
  TGraph2D* g = 0;
  while ((g = static_cast<TGraph2D*>(next()))) AddHistogram(*g);
  if (det->GetInner()) AddHistogram(fInnerHits);
  if (det->GetOuter()) AddHistogram(fOuterHits);

  fFrame = new TH2F(Form("FMD%d", fId), Form("FMD%d", fId), 
		    1, -fMaxR, fMaxR, 1, -fMaxR, fMaxR);
  fFrame->SetStats(kFALSE);
  fFrame->GetXaxis()->SetTitle("x [cm]");
  fFrame->GetYaxis()->SetTitle("y [cm]");
  fFrame->GetZaxis()->SetTitle("z [cm]");
  fFrame->SetDirectory(0);
}

  
//____________________________________________________________________
Bool_t 
AliFMDFancy::Init()
{
  // Initialize.  GEt transforms and such, 
  if (!AliFMDInput::Init()) return kFALSE;
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  geom->Init();
  geom->InitTransformations();
  
  fFMD1.Init();
  fFMD2.Init();
  fFMD3.Init();
  return kTRUE;
}

//____________________________________________________________________
void
AliFMDFancy::Detector::Begin(Int_t /* event */)
{
  TIter next(&fShapes);
  TGraph2D* g = 0;
  fFrame->Draw("surf fb");
  fFrame->GetZaxis()->SetRangeUser(fMinZ, fMaxZ);
  while ((g = static_cast<TGraph2D*>(next()))) g->Draw("tri2 FB same");
}

//____________________________________________________________________
void
AliFMDFancy::Detector::Clear(Int_t /* event */)
{
  fNInnerHits = 0;
  fNOuterHits = 0;
}

  
//____________________________________________________________________
Bool_t 
AliFMDFancy::Begin(Int_t event) 
{
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
      fFMD1.Begin(event);
    }
    if ((det = geom->GetDetector(2))) {
      fPad->cd();
      fFMD2Pad = new TPad("FMD2", "FMD2", 0.5, 0.50, 1.0, 1.0, 0, 0);
      fFMD2Pad->Draw();
      fFMD2Pad->cd();
      fFMD2.Begin(event);
    }
    if ((det = geom->GetDetector(3))) {
      fPad->cd();
      fFMD3Pad = new TPad("FMD3", "FMD3", 0.0, 0.05, .5, .5, 0, 0);
      fFMD3Pad->Draw();
      fFMD3Pad->cd();
      fFMD3.Begin(event);
    }
    fPad->cd();
    fSummary = new TPad("display", "Display", 0.5, 0.05, 1.0, 0.5, 0, 0);
    fSummary->Draw();
    fSummary->cd();
    fEvent.Draw();
    fFMD1IHits.Draw();
    fFMD2IHits.Draw();
    fFMD2OHits.Draw();
    fFMD3IHits.Draw();
    fFMD3OHits.Draw();
    fLine.Draw();
    fTotal.Draw();
  }
  fEvent.SetTitle(Form("Event # %6d", event));
  // fEvent.Modify();
  fCanvas->Modified();
  fCanvas->Update();
  fCanvas->cd();
  fFMD1.Clear(event);
  fFMD2.Clear(event);
  fFMD3.Clear(event);
  return AliFMDInput::Begin(event);
}

//____________________________________________________________________
void
AliFMDFancy::Detector::End()
{
  Char_t  rs[] = { 'I', 'O', '\0' };
  Char_t* rp   = rs;
  Char_t  r;
  while ((r = *(rp++))) {
    TGraph2D& g = (r == 'I' ? fInnerHits : fOuterHits);
    Int_t&    n = (r == 'I' ? fNInnerHits : fNOuterHits);
    Int_t     m = (r == 'I' ? 512 * 10 * 2 : 256 * 20 * 2);
    if (n == 0) continue;
    for (Int_t i = n; i < g.GetN(); i++)  g.RemovePoint(i);
    std::cout << g.GetName() << " has " << std::setw(4) << n << "/" 
	      << std::setw(5) << m << " points" << std::endl;
    g.Draw("same fb p");
    AddHistogram(g, "empty");
  }
}

//____________________________________________________________________
Bool_t 
AliFMDFancy::End() 
{
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  AliFMDDetector* det;
  Int_t total = 0;
  if ((det = geom->GetDetector(1))) {
    fFMD1Pad->cd();
    fFMD1.End();
    fFMD1Pad->Modified();
    fFMD1IHits.SetTitle(Form("# hits in FMD1I:  %5d", fFMD1.fNInnerHits));
    total += fFMD1.fNInnerHits;
  }
  if ((det = geom->GetDetector(2))) {
    fFMD2Pad->cd();
    fFMD2.End();
    fFMD2Pad->Modified();
    fFMD2IHits.SetTitle(Form("# hits in FMD2I:  %5d", fFMD2.fNInnerHits));
    fFMD2OHits.SetTitle(Form("# hits in FMD2O: %5d", fFMD2.fNOuterHits));
    total += fFMD2.fNInnerHits;
    total += fFMD2.fNOuterHits;    
  }
  if ((det = geom->GetDetector(3))) {
    fFMD3Pad->cd();
    fFMD3.End();
    fFMD3Pad->Modified();
    fFMD3IHits.SetTitle(Form("# hits in FMD3I:  %5d", fFMD3.fNInnerHits));
    fFMD3OHits.SetTitle(Form("# hits in FMD3O: %5d", fFMD3.fNOuterHits));
    total += fFMD3.fNInnerHits;
    total += fFMD3.fNOuterHits;    
  }
  fTotal.SetTitle(Form("Total:    %5d/51200 (%3d%%)", 
		      total, Int_t(100. / 51200 * total)));
  fSummary->Modified();
  fCanvas->Modified();
  fCanvas->Update();
  fCanvas->cd();
  fWait = kTRUE;
  while (fWait) {
    gApplication->StartIdleing();
    gSystem->InnerLoop();
    gApplication->StopIdleing();
  }
  return AliFMDInput::End();
}

//____________________________________________________________________
Bool_t 
AliFMDFancy::ProcessHit(AliFMDHit* hit, TParticle*) 
{
  AddMarker(hit->Detector(), hit->Ring(), hit->Sector(), hit->Strip(), 
	    hit, hit->Edep(), 0);
  return kTRUE;
}
//____________________________________________________________________
void
AliFMDFancy::Detector::AddMarker(Char_t rng, UShort_t sec, UShort_t str,
				 Float_t, Float_t)
{
  AliFMDGeometry*   geom = AliFMDGeometry::Instance();
  Double_t x, y, z;
  geom->Detector2XYZ(fId, rng, sec, str, x, y, z);
  if (true) {
    AliFMDRing* r  = geom->GetRing(rng);
    Double_t    t  = .9 * r->GetTheta() / 2;
    Double_t    a  = gRandom->Uniform(-t,t) * TMath::Pi() / 180;
    Double_t    x1 = x * TMath::Cos(a) - y * TMath::Sin(a);
    Double_t    y1 = x * TMath::Sin(a) + y * TMath::Cos(a);
    x = x1;
    y = y1;
  }
  switch (rng) {
  case 'I':
  case 'i': fInnerHits.SetPoint(fNInnerHits++, x, y, z); break;
  case 'O': 
  case 'o': fOuterHits.SetPoint(fNOuterHits++, x, y, z); break;
  default:  return;
  }
}

    

//____________________________________________________________________
void
AliFMDFancy::AddMarker(UShort_t det, Char_t rng, UShort_t sec, UShort_t str,
		       TObject*, Float_t, Float_t)
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
  switch (det) {
  case 1: fFMD1.AddMarker(rng,sec,str,0,0); break;
  case 2: fFMD2.AddMarker(rng,sec,str,0,0); break;
  case 3: fFMD3.AddMarker(rng,sec,str,0,0); break;
  }
}

//____________________________________________________________________
//
// EOF
//
