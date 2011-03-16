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

// #include <TApplication.h>
// #include <TButton.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TMath.h>
#include <TPad.h>
#include <TRandom.h>
// #include <TSlider.h>
#include <TStyle.h>
// #include <TSystem.h>
#include <TVector2.h>
// #include <TView.h>
#include <TGraph.h>
#include "AliFMDPattern.h"	// ALIFMDDISPLAY_H
#include "AliFMDGeometry.h"	// ALIFMDGEOMETRY_H
//#include "AliFMDParameters.h"	// ALIFMDPARAMETERS_H
#include "AliFMDRing.h"
// #include "AliFMDDetector.h"
#include "AliFMDHit.h"
#include "AliMultiplicity.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
// #include <AliLog.h>
#include "AliFMDDebug.h" // Better debug macros
// #include "AliPhysicsSelection.h"
class AliFMDDetector;

//____________________________________________________________________
ClassImp(AliFMDPattern)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDPattern::AliFMDPatternDetector::AliFMDPatternDetector(UShort_t id) 
  : fId(id),
    fCounts(0), 
    fGraphs(0), 
    fFrame(0), 
    fInners(10), 
    fOuters(id == 1 ? 0 : 20)
{
  // CTOR 
  // 
  // Parameters: 
  // 
  //   ID 	Identifier 
}

//____________________________________________________________________
AliFMDPattern::AliFMDPatternDetector::~AliFMDPatternDetector()
{
  // DTOR 
  // Destructor - 
  // deletes mother frame 
  if (fFrame) delete fFrame;
}

//____________________________________________________________________
void
AliFMDPattern::AliFMDPatternDetector::DrawShape(const TObjArray& a) 
{
  // Draw all shapes. 
  // 
  // Paramters 
  // 
  //	a	Array of shapes 
  //
  TIter next(&a);
  TGraph* g = 0;
  while ((g = static_cast<TGraph*>(next()))) {
    g->DrawClone("f same");
    g->DrawClone("l same");
  }
}

//____________________________________________________________________
void
AliFMDPattern::AliFMDPatternDetector::CopyShapes(const TObjArray& src, 
						 TObjArray& dest, 
						 Double_t ang, 
						 Double_t fx, 
						 Double_t fy)
{
  // 
  // Copy shapes
  // 
  // Parameters:
  //    input  Source
  //    own    Ours
  //    ang    Angle 
  //    fx     Factor x
  //    fy     Factor y
  //
  TIter     next(&src);
  TGraph*   g = 0;
  while ((g = static_cast<TGraph*>(next()))) { 
    TGraph* gg = new TGraph(*g);
    Double_t* x  = gg->GetX();
    Double_t* y  = gg->GetY();
    for (Int_t i = 0; i < gg->GetN(); i++) { 
      Float_t xx = x[i] * TMath::Cos(ang) - y[i] * TMath::Sin(ang);
      Float_t yy = x[i] * TMath::Sin(ang) + y[i] * TMath::Cos(ang);
      gg->SetPoint(i, fx * xx, fy * yy);
    }
    gg->SetFillStyle(g->GetFillStyle());
    gg->SetFillColor(g->GetFillColor());
    gg->SetLineStyle(g->GetLineStyle());
    gg->SetLineColor(g->GetLineColor());
    gg->SetLineWidth(g->GetLineWidth());
    gg->SetMarkerStyle(g->GetMarkerStyle());
    gg->SetMarkerColor(g->GetMarkerColor());
    gg->SetMarkerSize(g->GetMarkerSize());
    TString name(g->GetName());
    name.ReplaceAll("X", Form("%d",fId));
    gg->SetName(name.Data());
    TString title(g->GetTitle());
    title.ReplaceAll("X", Form("%d",fId));
    gg->SetTitle(title.Data());
    dest.Add(gg);
  }
  dest.SetOwner();
}

//____________________________________________________________________
void
AliFMDPattern::AliFMDPatternDetector::Begin(Int_t      nlevel, 
					    Double_t   r, 
					    TObjArray& inners, 
					    TObjArray& outers)
{
  // Start of a run. 
  // 
  // Parameters 
  // 
  //	nlevel		Number of levels 
  //    r		Radius 
  //	inners		Array of inner shapes 
  //	outers		Array of outer shapes 
  //

  // To make code-checker shut up
  TStyle* style = gStyle;  
  if (nlevel < 1) nlevel = style->GetNumberOfColors();
  fCounts.Set(nlevel);
  Double_t rr = 1.05 * r;
  if (!fFrame) {
    // The code-checker thinks this is not using the declaration of
    // TH2F - what a morron!   
    fFrame = new TH2F(Form("fmd%dFrame", fId), Form("FMD%d", fId), 
		      100, -rr, rr, 100, -rr, rr);
    fFrame->SetStats(kFALSE);
    fFrame->Draw();
  }
  Double_t ang = (fId == 1 ? -TMath::Pi() / 2 : 0);
  Double_t fx  = (fId == 3 ? -1               : 1); // Flip around Y
  Double_t fy  = (fId == 1 ?  1               : 1); // Flip around X
  
  CopyShapes(inners, fInners, ang, fx, fy);
  DrawShape(fInners);
  if (fId != 1) { 
    CopyShapes(outers, fOuters, ang, fx, fy);
    DrawShape(fOuters);
  }

  for (Int_t i = 0; i < nlevel; i++) { 
    TGraph* g = new TGraph;
    Int_t idx = Int_t(Float_t(i) / nlevel * style->GetNumberOfColors());
    Int_t col = style->GetColorPalette(idx);
    g->SetName(Form("FMD%d_L%02d", fId, i));
    g->SetMarkerColor(col);
    g->SetLineColor(col);
    g->SetFillColor(col);
    g->SetMarkerSize(i * .2 + .2);
    g->SetMarkerStyle(2);
    g->SetEditable(kFALSE);
    g->Draw("same p");
    fGraphs.AddAtAndExpand(g, i);
  }
  // TIter   next(&fGraphs);
}

//____________________________________________________________________
void
AliFMDPattern::AliFMDPatternDetector::Clear() 
{
  // Clear this display.  
  // Simply reset counters to zero. 
  // Avoid deleting memory. 
  fCounts.Reset(0);
}

//____________________________________________________________________
void
AliFMDPattern::AliFMDPatternDetector::End()
{
  // Called when displaying the data.  
  // Simply resets number of points at each level to 
  // the seen number of hits at that level. 
  // Avoid deleting memory. 
  
 
  TIter   next(&fGraphs);
  TGraph* g = 0;
  Int_t   i = 0;
  while ((g = static_cast<TGraph*>(next()))) { 
    Int_t cnt = fCounts[i++];
    if (cnt > 0) { 
      g->Set(cnt);
      g->SetMarkerSize(i * .2 + .2);
    }
    else {
      g->SetPoint(0,0,0);
      g->SetMarkerSize(0);
    }
  }
  
}
//____________________________________________________________________
void
AliFMDPattern::AliFMDPatternDetector::AddMarker(Double_t x, 
						Double_t y, 
						Float_t  s, 
						Float_t  max)
{
  // Add a marker at (X,Y,Z).  The marker color and size is chosen
  // relative to the MAX argument. 
  // 
  // Parameters 
  // 
  //	X,Y,Z		Coordiantes 
  //	MAX		Maximum value. 
  // 
  // Sigh, for some odd reason, the code-checker does not recognise
  // this a usage of the TMath namespace declaration! Idiot 
  // 
  Int_t i = TMath::Min(Int_t(fCounts.fN * s / max),  
		       Int_t(fGraphs.GetEntries()-1));
  if (i < 0 || i >= fCounts.fN) { 
    std::cerr << "Graph index " << i << " out of bounds [0," 
	      << fCounts.fN << ") - " 
	      << fCounts.fN << " * " << s << " / " << max << std::endl;
    return;
  }
  TGraph* g = static_cast<TGraph*>(fGraphs.At(i));
  if (!g) return;
  g->SetPoint(fCounts[i]++, x, y);
}


//____________________________________________________________________
AliFMDPattern::AliFMDPattern(const char* gAliceFile)
  : AliFMDDisplay(kTRUE, gAliceFile),
    fInners(0), 
    fOuters(0),
    fInnerMax(0), 
    fOuterMax(0),
    fFMD1Pad(0),
    fFMD1(1),
    fFMD2Pad(0),
    fFMD2(2),
    fFMD3Pad(0),
    fFMD3(3),
    fSummary(0),
    fEvent(.1, .8, "Event #"),
    fFMD1Sum(.2, .7, "# in FMD1: "),
    fFMD2Sum(.2, .6, "# in FMD2: "),
    fFMD3Sum(.2, .5, "# in FMD3: "),
    fLine(.15, .47, .85, .47),
    fTotal(.2, .35, "Total:   "), 
    fFMD1Area(0),
    fFMD2Area(0),
    fFMD3Area(0)// ,fPhysicsSelection(0)
{
  // Constructor. 
  // 
  // Parameters 
  // 
  //   gAliceFile	The galice.root file to use - if any. 
  // 

  SetName("AliFMDPattern");
  SetName("2D display of FMD data");
  // fPhysicsSelection = new AliPhysicsSelection();
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
  // DTOR 
  // Free all allocated shapes. 
  // note, that most members are real objects, so we do not need to
  // deal with them here. 
  fInners.Delete();
  fOuters.Delete();
}

    
//____________________________________________________________________
Bool_t 
AliFMDPattern::Init()
{
  // Initialize.  Get transforms and such, 
  if (!AliFMDInput::Init()) return kFALSE;
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  if (!geom) return kFALSE;
  geom->Init();
  geom->InitTransformations();
  
  fFMD1Area = 0;
  fFMD2Area = 0;
  fFMD3Area = 0;

  Double_t innerArea = 0;
  Double_t outerArea = 0;

  Char_t rs[] = { 'I' , 'O', '\0' };
  Char_t *r   = rs;
  do {
    AliFMDRing* ring = geom->GetRing(*r);
    if (!ring) continue;

    Double_t rl   = ring->GetMinR();
    Double_t rh   = ring->GetMaxR();
    Double_t area = rh * rh * TMath::Pi() - rl * rl * TMath::Pi();
    if (*r == 'I') innerArea = area;
    else           outerArea = area;
      

    const TObjArray& vs = ring->GetVerticies();
    TObjArray&       gs = (*r == 'I' ? fInners   : fOuters);
    Float_t&         mr = (*r == 'I' ? fInnerMax : fOuterMax);
    Int_t            nm = ring->GetNModules();
    AliFMDDebug(1, ("Making %d modules for %c", nm, *r));
    for (Int_t m = 0; m < nm; m++) {
      Int_t          nv = 6; // vs.GetEntries();
      Double_t       a  = TMath::Pi() / 180 * (m * 2 + 1) * ring->GetTheta();
      TGraph*        g  = new TGraph(nv+1);
      Double_t       x0 = 0, y0 = 0;
      gs.AddAtAndExpand(g, m);
      for (Int_t c = 1; c < 4; c++) {
	TVector2* v = static_cast<TVector2*>(vs.At(c));
	mr          = TMath::Max(mr, Float_t(v->Mod()));
	TVector2  w(v->Rotate(a));
	if (c == 1) { x0 = w.X(); y0 = w.Y(); }
	g->SetPoint(c-1, w.X(), w.Y());
      }
      for (Int_t c = 3; c > 0; c--) {
	TVector2* v = static_cast<TVector2*>(vs.At(c));
	TVector2  u(-v->X(), v->Y());
	mr          = TMath::Max(mr, Float_t(u.Mod()));
	TVector2  w(u.Rotate(a));
	g->SetPoint(3+(3-c), w.X(), w.Y());
      }
      g->SetName(Form("FMDX%c_%02d%02d", *r, 2*m,2*m+1));
      g->SetTitle(Form("FMDX%c, sectors %d and %d", *r, 2*m,2*m+1));
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

  fFMD1Area = innerArea;
  fFMD2Area = innerArea + outerArea;
  fFMD3Area = innerArea + outerArea;
  
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDPattern::Begin(Int_t event) 
{
  // Called at the begining of an event. 
  // 
  // Parameters 
  //
  //	EVENT		The event number 
  // 
  MakeAux();
  if (!fCanvas) {
    const char* which[] = { "Continue", "Start", "Pause", "Redisplay", 0 };
    MakeCanvas(which);
    
    AliFMDGeometry* geom = AliFMDGeometry::Instance();
    // AliFMDDetector* det;
    if ((/* det = */ geom->GetDetector(1))) {
      fPad->cd();
      fFMD1Pad = new TPad("FMD1", "FMD1", 0.0, 0.50, 0.5, 1.0, 0, 0);
      fFMD1Pad->Draw();
      fFMD1Pad->cd();
      fFMD1.Begin(-1, fInnerMax, fInners, fOuters);
    }
    if ((/* det = */ geom->GetDetector(2))) {
      fPad->cd();
      fFMD2Pad = new TPad("FMD2", "FMD2", 0.5, 0.50, 1.0, 1.0, 0, 0);
      fFMD2Pad->Draw();
      fFMD2Pad->cd();
      fFMD2.Begin(-1, fOuterMax, fInners, fOuters);
    }
    if ((/* det = */ geom->GetDetector(3))) {
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
  
#if 0
  TString triggers = fESDEvent->GetFiredTriggerClasses();
  const AliESDVertex* vertex = fESDEvent->GetPrimaryVertexSPD();
  Double_t vertexXYZ[3];
  vertex->GetXYZ(vertexXYZ);
  const AliMultiplicity* mult = fESDEvent->GetMultiplicity();
  Int_t nTrackLets = mult->GetNumberOfTracklets();
  std::cout<<triggers.Data()<<"  "<<fPhysicsSelection->IsCollisionCandidate(fESDEvent)<<"    "<<nTrackLets<<"   "<<vertexXYZ[0]<<"   "<<vertexXYZ[1]<<"   "<<vertexXYZ[2]<<std::endl;
#endif   
  
  return AliFMDInput::Begin(event);
}

//____________________________________________________________________
void
AliFMDPattern::Redisplay()
{
  // Redraw the displayu 
  fFMD1.Clear();
  fFMD2.Clear();
  fFMD3.Clear();
  AliFMDDisplay::Redisplay();
}

//____________________________________________________________________
void
AliFMDPattern::AtEnd()
{
  // Called at the end of an event. 
  DrawAux();
  
  Int_t total = 0;
  
  fFMD1.End();
  fFMD1Pad->Modified();
  fFMD1Sum.SetTitle(Form("# hits in FMD1: %5d   (%4.2f /cm^{2})", 
			 fFMD1.Total(), fFMD1.Total()/fFMD1Area));
  total += fFMD1.Total();

  fFMD2.End();
  fFMD2Pad->Modified();
  fFMD2Sum.SetTitle(Form("# hits in FMD2: %5d   (%4.2f /cm^{2})", 
			 fFMD2.Total(), fFMD2.Total()/fFMD2Area));
  total += fFMD2.Total();

  fFMD3.End();
  fFMD3Pad->Modified();
  fFMD3Sum.SetTitle(Form("# hits in FMD3: %5d   (%4.2f /cm^{2})", 
			 fFMD3.Total(), fFMD3.Total()/fFMD3Area));
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
  // Process a hit. 
  // 
  // Parameters 
  // 
  //	HIT		The hit to process. 
  // 
  // The TParticle argument is never used. 
  static const Float_t rMin  = fgkEdepRange.fLow;
  static const Float_t rMax  = fgkEdepRange.fHigh;

  if (!hit) { AliError("No hit");   return kFALSE; }
  // if (!p)   { AliError("No track"); return kFALSE; }
  Float_t  edep  = hit->Edep();

  if (fHits)                        fHits->Add(hit);
  if (fSpec)                        fSpec->Fill(edep);
  if (InsideCut(edep, rMin, rMax) && fSpecCut) fSpecCut->Fill(edep);

  switch (hit->Detector()) {
  case 1: fFMD1.AddMarker(hit->X(), hit->Y(), hit->Edep(), rMax); break;
  case 2: fFMD2.AddMarker(hit->X(), hit->Y(), hit->Edep(), rMax); break;
  case 3: fFMD3.AddMarker(hit->X(), hit->Y(), hit->Edep(), rMax); break;
  }
  return kTRUE;
}


//____________________________________________________________________
void
AliFMDPattern::AddMarker(UShort_t det, Char_t rng, 
			 UShort_t sec, UShort_t str,
			 TObject*, Float_t s, Float_t /*min*/, Float_t max)
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
  AliFMDPatternDetector* d = 0;
  switch (det) {
  case 1: d = &fFMD1; break;
  case 2: d = &fFMD2; break;
  case 3: d = &fFMD3; break;
  }
  if (!d) return;
  AliFMDGeometry*   geom = AliFMDGeometry::Instance();
  Double_t x, y, z;
  geom->Detector2XYZ(det, rng, sec, str, x, y, z);
  // Make code-checker shut the f**k up 
  TRandom* rand = gRandom;
  if (false) {
    AliFMDRing* r  = geom->GetRing(rng);
    Double_t    t  = .9 * r->GetTheta() / 2;
    Double_t    a  = rand->Uniform(-t,t) * TMath::Pi() / 180;
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
