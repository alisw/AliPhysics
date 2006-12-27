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
/** @file    AliFMDDisplay.cxx
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:39:09 2006
    @brief   FMD Event display 
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
#include "AliFMDDisplay.h"	// ALIFMDDISPLAY_H
#include "AliFMDHit.h"          // ALIFMDHIT_H
#include "AliFMDDigit.h"        // ALIFMDDIGIT_H
#include "AliFMDRecPoint.h"     // ALIFMDRECPOINT_H
#include "AliFMDGeometry.h"	// ALIFMDGEOMETRY_H
#include "AliFMDParameters.h"	// ALIFMDPARAMETERS_H
#include <AliESDFMD.h>          // ALIESDFMD_H
#include <AliLog.h>
#include <TStyle.h>
// #include <TArrayF.h>
#include <TMarker3DBox.h>
#include <TGeoManager.h>
// #include <TMath.h>
#include <TApplication.h>
#include <TButton.h>
// #include <TParticle.h>
#include <TCanvas.h>
#include <TView.h>
#include <TVirtualX.h>
#include <TSlider.h>
#include <TH1D.h>

//____________________________________________________________________
ClassImp(AliFMDDisplay)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDDisplay* AliFMDDisplay::fgInstance = 0;

//____________________________________________________________________
AliFMDDisplay* 
AliFMDDisplay::Instance()
{
  // Return static instance 
  return fgInstance;
}

//____________________________________________________________________
AliFMDDisplay::~AliFMDDisplay()
{
  if (fMarkers) {
    fMarkers->Delete();
    delete fMarkers;
  }
  if (fHits) {
    fHits->Clear();
    delete fHits;
  }
  if (fPad)     delete fPad;
  fButtons.Delete();
  if (fSlider)  delete fSlider;
  if (fCanvas)  delete fCanvas;
}
  
//____________________________________________________________________
AliFMDDisplay::AliFMDDisplay(Bool_t onlyFMD, const char* gAliceFile)
  : AliFMDInput(gAliceFile),
    fWait(kFALSE),
    fMarkers(0),
    fHits(0),
    fCanvas(0), 
    fPad(0), 
    fSlider(0),
    fZoomMode(kFALSE),
    fX0(0),
    fY0(0),
    fX1(0),
    fY1(0),
    fMultCut(0),
    fPedestalFactor(0),
    fXPixel(0),
    fYPixel(0),
    fOldXPixel(0),
    fOldYPixel(0),
    fLineDrawn(0),
    fOnlyFMD(onlyFMD)
{
  // Constructor of an FMD display object. 
  AddLoad(kGeometry);
  if (fgInstance) delete fgInstance;
  fgInstance = this;
  SetMultiplicityCut();
  SetPedestalFactor();
}

//____________________________________________________________________
void           
AliFMDDisplay::MakeCanvas(const char** which)
{
  gStyle->SetPalette(1);
  Double_t y1 = .10;
  Int_t    w  = 700;
  fCanvas = new TCanvas("display", "Display", w, Int_t(w / (1-y1)));
  fCanvas->SetFillColor(1);
  fCanvas->ToggleEventStatus();
  fCanvas->cd();
  fPad = new TPad("view", "3DView", 0.0, y1, 1.0, 1.0, 1, 0, 0);
  fPad->Draw();

  Double_t yb = 0;
  fCanvas->cd();
  if (TESTBIT(fTreeMask, kESD) || 
      TESTBIT(fTreeMask, kDigits) || 
      TESTBIT(fTreeMask, kRaw)) {
    yb = .05;
    fSlider = new TSlider("multCut", "Multiplicity cut", 0, 0, 1, yb);
    fSlider->SetMethod("AliFMDDisplay::Instance()->ChangeCut()");
    fSlider->Draw();
    fSlider->SetMinimum(TESTBIT(fTreeMask, kESD) ? fMultCut * 10 :
			fPedestalFactor * 10);
  }
  const char** p = which;
  const char*  m;
  Int_t        n  = 0;
  Int_t        j  = 0;
  while (*(p++)) n++;
  AliInfo(Form("Got %d buttons", n));
  Float_t      x0 = 0;
  Float_t      dx = 1. / n;
  p               = which;
  while ((m = *(p++))) {
    fCanvas->cd();
    AliInfo(Form("Adding button %s", m));
    TButton* b = new TButton(m, Form("AliFMDDisplay::Instance()->%s()", m),
			     x0, yb, x0 + dx, y1);
    b->Draw();
    fButtons.Add(b);
    x0 += dx;
    j++;
  }
}

//____________________________________________________________________
void           
AliFMDDisplay::ShowOnlyFMD()
{
  if (!fGeoManager) return;
  static bool once = false;
  if (once) return;
  once = true;
  AliInfo("Will only show the FMD");
  TGeoVolume* top = gGeoManager->GetTopVolume();
  top->InvisibleAll(kTRUE);
  TGeoIterator next(top);
  TGeoNode* node;
  TGeoVolume* v = 0;
  Bool_t hasFMD1 = kFALSE;
  Bool_t hasFMD2 = kFALSE;
  Bool_t hasFMD3 = kFALSE;
  TObjArray toshow;
  while ((node = static_cast<TGeoNode*>(next()))) {
    const char* name = node->GetName();
    if (!name) continue;
    if (!(v = node->GetVolume())) continue;

    if (name[0] == 'F') {
      if (name[2] == 'M' && (name[3] == 'T' || name[3] == 'B')) {
	// Virtual Master half-ring volume - top-level
	Int_t det = node->GetNumber();
	switch (det) {
	case 1: hasFMD1 = true; break;
	case 2: hasFMD2 = true; break;
	case 3: hasFMD3 = true; break;
	default: continue;
	}
	toshow.Add(v);
      }
      else if (name[3] == 'V' && (name[2] == 'T' || name[2] == 'B')) 
	toshow.Add(v); // Virtual Half-ring, bare detectors
      // else if (name[3] == 'H' && (name[2] == 'F' || name[2] == 'B')) 
      //  toshow.Add(v); // Virtual Hybrid container 
    }
    v->SetVisibility(kFALSE);
    v->SetVisDaughters(kFALSE);
    v->InvisibleAll(kTRUE);
  }
  TIter i(&toshow);
  while ((v = static_cast<TGeoVolume*>(i()))) {
    v->SetVisibility(kTRUE);
    v->SetVisDaughters(kTRUE);
    v->InvisibleAll(kFALSE);
  }  
}

    
//____________________________________________________________________
void           
AliFMDDisplay::ExecuteEvent(Int_t event, Int_t px, Int_t py) 
{
  // AliInfo(Form("Event %d, at (%d,%d)", px, py));
  if (px == 0 && py == 0) return;
  if (!fZoomMode && fPad->GetView()) {
    fPad->GetView()->ExecuteRotateView(event, px, py);
    return;
  }
  fPad->SetCursor(kCross);
  switch (event) {
  case kButton1Down: 
    fPad->TAttLine::Modify();
    fX0        = fPad->AbsPixeltoX(px);
    fY0        = fPad->AbsPixeltoY(py);
    fXPixel    = fOldXPixel = px;
    fYPixel    = fOldYPixel = py;
    fLineDrawn = kFALSE;
    return;
  case kButton1Motion:
    if (fLineDrawn) 
      gVirtualX->DrawBox(fXPixel, fYPixel, fOldXPixel, fOldYPixel, 
			 TVirtualX::kHollow);
    fOldXPixel = px;
    fOldYPixel = py;
    fLineDrawn = kTRUE;
    gVirtualX->DrawBox(fXPixel, fYPixel, fOldXPixel, fOldYPixel, 
		       TVirtualX::kHollow);
    return;
  case kButton1Up:
    fPad->GetCanvas()->FeedbackMode(kFALSE);
    if (px == fXPixel || py == fYPixel) return;
    fX1 = fPad->AbsPixeltoX(px);
    fY1 = fPad->AbsPixeltoY(py);
    if (fX1 < fX0) std::swap(fX0, fX1); 
    if (fY1 < fY0) std::swap(fY0, fY1); 
    fPad->Range(fX0, fY0, fX1, fY1);
    fPad->Modified();
    return;
  }
}

//____________________________________________________________________
Int_t          
AliFMDDisplay::DistancetoPrimitive(Int_t px, Int_t) 
{
  // AliInfo(Form("@ (%d,%d)", px, py));
  fPad->SetCursor(kCross);
  Float_t xmin = fPad->GetX1();
  Float_t xmax = fPad->GetX2();
  Float_t dx   = .02 * (xmax - xmin);
  Float_t x    = fPad->AbsPixeltoX(px);
  if (x < xmin + dx || x > xmax - dx) return 9999;
  return (fZoomMode ? 0 : 7);
}
//____________________________________________________________________
Bool_t 
AliFMDDisplay::Init()
{
  // Initialize.  GEt transforms and such, 
  if (!AliFMDInput::Init()) return kFALSE;
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  geom->Init();
  geom->InitTransformations();
  if (TESTBIT(fTreeMask, kDigits) || TESTBIT(fTreeMask, kRaw)) 
    AliFMDParameters::Instance()->Init();

  fMarkers = new TObjArray;
  fHits    = new TObjArray;
  fMarkers->SetOwner(kTRUE);
  fHits->SetOwner(kFALSE);
  return kTRUE;
}

//____________________________________________________________________
void
AliFMDDisplay::MakeAux()
{
  if ((TESTBIT(fTreeMask, kESD) || 
       TESTBIT(fTreeMask, kDigits) || 
       TESTBIT(fTreeMask, kRaw))) {
    if (!fAux) {
      fAux = new TCanvas("aux", "Aux");
      fAux->SetLogy();
      if (TESTBIT(fTreeMask, kESD)) 
	fSpec = new TH1D("spec", "Mult spectra", 500, 0, 10);
      else 
	fSpec = new TH1D("spec", "Adc spectra", 1024, -.5, 1023.5);
      fSpecCut = static_cast<TH1*>(fSpec->Clone("specCut"));
      fSpec->SetFillColor(2);
      fSpec->SetFillStyle(3001);
      fSpecCut->SetFillColor(4);
      fSpecCut->SetFillStyle(3001);
    }
    else {
      fSpec->Reset();
      fSpecCut->Reset();
    }
  }
}

//____________________________________________________________________
void
AliFMDDisplay::DrawAux()
{
  if (!fAux) return;
  fAux->cd();
  fAux->Clear();
  fSpec->Draw();
  fSpecCut->Draw("same");
  fAux->Modified();
  fAux->Update();
  fAux->cd();
}

//____________________________________________________________________
Bool_t 
AliFMDDisplay::Begin(Int_t event) 
{
  // Begin of event.  Make canvas is not already done 
  if (!fCanvas) {
    const char* m[] = { "Continue", "Zoom", "Pick", "Redisplay", 0 }; 
    MakeCanvas(m);
  }
  MakeAux();

  // AliInfo("Clearing canvas");
  // fCanvas->Clear();
  if (!fGeoManager) {
    Warning("End", "No geometry manager");
    return kFALSE;
  }
  AliInfo("Drawing geometry");
  fPad->cd();
  fGeoManager->GetTopVolume()->Draw();
  if (fOnlyFMD) ShowOnlyFMD();
  AliInfo("Adjusting view");
  Int_t irep;
  if (fPad->GetView()) {
    fPad->GetView()->SetView(-200, -40, 80, irep);
    fPad->GetView()->Zoom();
    fPad->Modified();
    fPad->cd();
  }
  return AliFMDInput::Begin(event);
}

//____________________________________________________________________
void
AliFMDDisplay::AtEnd() 
{
  fPad->cd();
  fMarkers->Draw();
  fPad->cd();
  AppendPad();
  fPad->cd();
  DrawAux();
}

//____________________________________________________________________
void
AliFMDDisplay::Idle() 
{
  fWait = kTRUE;
  while (fWait) {
    gApplication->StartIdleing();
    gSystem->InnerLoop();
    gApplication->StopIdleing();
  }
  AliInfo("After idle loop");
  if (fMarkers) fMarkers->Delete();
  if (fHits)    fHits->Clear();
  AliInfo("After clearing caches");
}

//____________________________________________________________________
Bool_t 
AliFMDDisplay::End()
{
  // End of event.  Draw everything 
  AtEnd();
  Idle();
  return AliFMDInput::End();
}

//____________________________________________________________________
Int_t
AliFMDDisplay::LookupColor(Float_t x, Float_t max) const
{
  // Look-up color 
  Int_t idx = Int_t(x / max * gStyle->GetNumberOfColors());
  return gStyle->GetColorPalette(idx);
}

//____________________________________________________________________
void
AliFMDDisplay::ChangeCut() 
{
  fMultCut        = fSlider->GetMinimum() * 10;
  fPedestalFactor = fSlider->GetMinimum() * 10;
  AliInfo(Form("Multiplicity cut: %7.5f, Pedestal factor: %7.4f (%6.5f)", 
	       fMultCut, fPedestalFactor, fSlider->GetMinimum()));
  Redisplay();
}

//____________________________________________________________________
void
AliFMDDisplay::Redisplay()
{
  if (fMarkers) fMarkers->Delete();
  if (fHits)    fHits->Clear();
  if (fSpec)    fSpec->Reset();
  if (fSpecCut) fSpecCut->Reset();
  Event();
  AtEnd();
}

//____________________________________________________________________
void
AliFMDDisplay::AddMarker(UShort_t det, Char_t rng, UShort_t sec, UShort_t str,
			 TObject* o, Float_t s, Float_t max)
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
  AliFMDGeometry*   geom = AliFMDGeometry::Instance();
  Double_t x, y, z;
  geom->Detector2XYZ(det, rng, sec, str, x, y, z);
  Float_t  size  = .1;
  Float_t  zsize = s / max * 10;
  Float_t  r     = TMath::Sqrt(x * x + y * y);
  Float_t  theta = TMath::ATan2(r, z);
  Float_t  phi   = TMath::ATan2(y, x);
  Float_t  rz    = z + (z < 0 ? 1 : -1) * zsize;
  TMarker3DBox* marker = new  TMarker3DBox(x,y,rz,size,size,zsize,theta,phi);
  if (o) marker->SetRefObject(o);
  marker->SetLineColor(LookupColor(s, max));
  fMarkers->Add(marker);
}
  

//____________________________________________________________________
Bool_t 
AliFMDDisplay::ProcessHit(AliFMDHit* hit, TParticle* /* p */) 
{
  // Process a hit 
  if (!hit) { AliError("No hit");   return kFALSE; }
  // if (!p)   { AliError("No track"); return kFALSE; }

  if (fHits) fHits->Add(hit);
  Float_t  size  = .1;
  Float_t  zsize = TMath::Sqrt(hit->Edep() * 20);
  Float_t  z     = hit->Z() + (hit->Z() < 0 ? 1 : -1) * zsize; 
  Float_t  pt    = TMath::Sqrt(hit->Py()*hit->Py()+hit->Px()*hit->Px());
  Float_t  theta = TMath::ATan2(pt, hit->Pz());
  Float_t  phi   = TMath::ATan2(hit->Py(), hit->Px());
  TMarker3DBox* marker = new  TMarker3DBox(hit->X(), hit->Y(), z,
					   size, size, zsize, theta, phi);
  marker->SetLineColor(LookupColor(hit->Edep(), 1));
  marker->SetRefObject(hit);
  fMarkers->Add(marker);
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDDisplay::ProcessDigit(AliFMDDigit* digit)
{
  // Process a digit 
  if (!digit) { AliError("No digit");   return kFALSE; }

  AliFMDParameters* parm = AliFMDParameters::Instance();
  UShort_t det           =  digit->Detector();
  Char_t   ring          =  digit->Ring();
  UShort_t sec           =  digit->Sector();
  UShort_t str           =  digit->Strip();
  Double_t ped           =  parm->GetPedestal(det,ring, sec, str);
  Double_t pedW          =  parm->GetPedestalWidth(det,ring, sec, str);
  Double_t threshold     =  ped + fPedestalFactor * pedW;
  Float_t  counts        =  digit->Counts();
  AliDebug(10, Form("FMD%d%c[%2d,%3d] ADC: %d > %d (=%4.2f+%4.2f*%4.2f)", 
		    digit->Detector(), digit->Ring(), digit->Sector(), 
		    digit->Strip(), Int_t(counts), Int_t(threshold), 
		    ped, fPedestalFactor, pedW));
  if (fSpec) fSpec->Fill(counts);
  if (counts < threshold) return kTRUE;
  if (fHits) fHits->Add(digit);
  if (fSpecCut) fSpecCut->Fill(counts);

  AddMarker(det, ring, sec, str, digit, counts, 1024);
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDDisplay::ProcessRaw(AliFMDDigit* digit)
{
  // PRocess raw data 
  return ProcessDigit(digit);
}

//____________________________________________________________________
Bool_t 
AliFMDDisplay::ProcessRecPoint(AliFMDRecPoint* recpoint)
{
  // Process reconstructed point 
  if (!recpoint) { AliError("No recpoint");   return kFALSE; }
  if (recpoint->Particles() < fMultCut) return kTRUE;
  if (fHits) fHits->Add(recpoint);
  AddMarker(recpoint->Detector(), recpoint->Ring(), recpoint->Sector(),  
	    recpoint->Strip(), recpoint, recpoint->Particles(), 20);
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDDisplay::ProcessESD(UShort_t det, Char_t rng, UShort_t sec, UShort_t str,
			  Float_t, Float_t mult)
{
  Double_t cmult = mult;
  if (fSpec) fSpec->Fill(cmult);
  if (cmult < fMultCut || cmult == AliESDFMD::kInvalidMult) return kTRUE;
  AddMarker(det,rng,sec,str, 0, cmult, 20);
  if (fSpecCut) fSpecCut->Fill(cmult);
  return kTRUE;
}

//____________________________________________________________________
//
// EOF
//
