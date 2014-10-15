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
/** 
 * @file    AliFMDDisplay.cxx
 * @author  Christian Holm Christensen <cholm@nbi.dk>
 * @date    Mon Mar 27 12:39:09 2006
 * @brief   FMD Event display 
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

#include <TSystem.h>
#include <TApplication.h>
#include <TButton.h>
#include <TCanvas.h>
#include <TGeoManager.h>
#include <TH1D.h>
#include <TMarker3DBox.h>
#include <TMath.h>
#include <TSlider.h>
#include <TSliderBox.h>
#include <TStyle.h>
#include <TView.h>
#include <TVirtualX.h>
#include <TVirtualViewer3D.h>
#include <TList.h>
// #include <TArrayF.h>
// #include <TParticle.h>

#include "AliFMDDisplay.h"	// ALIFMDDISPLAY_H
#include "AliFMDHit.h"          // ALIFMDHIT_H
#include "AliFMDDigit.h"        // ALIFMDDIGIT_H
#include "AliFMDSDigit.h"       // ALIFMDSDIGIT_H
#include "AliFMDRecPoint.h"     // ALIFMDRECPOINT_H
#include "AliFMDGeometry.h"	// ALIFMDGEOMETRY_H
#include "AliFMDParameters.h"	// ALIFMDPARAMETERS_H
#include "AliFMDRawReader.h"    // ALIFMDRAWREADER_H
#include <AliESDFMD.h>          // ALIESDFMD_H
// #include <AliLog.h>
#include "AliFMDDebug.h" // Better debug macros

//____________________________________________________________________
ClassImp(AliFMDDisplay)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDDisplay* AliFMDDisplay::fgInstance = 0;

//____________________________________________________________________
const AliFMDDisplay::Range_t AliFMDDisplay::fgkEdepRange = {  100, 0.,    2. };
const AliFMDDisplay::Range_t AliFMDDisplay::fgkAdcRange  = { 1024, 0., 1023. };
const AliFMDDisplay::Range_t AliFMDDisplay::fgkMultRange = {  500, 0.,   20. };

  
//____________________________________________________________________
AliFMDDisplay* 
AliFMDDisplay::Instance()
{
  // Return static instance 
  // If the instance does not exist
  // it is not created!
  return fgInstance;
}

//____________________________________________________________________
AliFMDDisplay::~AliFMDDisplay()
{
  // Destructor. 
  // Cleans 
  // up 
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
    fButtons(0),
    fSlider(0),
    fFactor(0),
    fZoomMode(kFALSE),
    fX0(0),
    fY0(0),
    fX1(0),
    fY1(0),
    fXPixel(0),
    fYPixel(0),
    fOldXPixel(0),
    fOldYPixel(0),
    fLineDrawn(0),
    fOnlyFMD(onlyFMD),
    fSpec(0), 
    fSpecCut(0),
    fAux(0), 
    fReturn(kFALSE), 
    fContinous(kFALSE),
    fTimeout("gApplication->StopIdleing()", 10), 
    fInitialMin(0), 
    fInitialMax(1), 
    fInitialFactor(3/10.)
{
  // Constructor of an FMD display object. 
  // Must be called 
  // before Instance 
  SetName("AliFMDDisplay");
  SetTitle("3D Display of various kinds of FMD data");
  AddLoad(kGeometry);
  if (fgInstance) delete fgInstance;
  fgInstance = this;
}

//____________________________________________________________________
void           
AliFMDDisplay::MakeCanvas(const char** which)
{
  // Make a canvas 
  // Parameters: 
  //   which   Which button to put up. 
  gStyle->SetPalette(1);
  // gStyle->SetCanvasPreferGL(kTRUE);
  Double_t y1 = .10;
  Int_t    w  = 700;
  fCanvas = new TCanvas(Form("gl%s", GetName()), 
			Form("%s - Display", GetTitle()), 
			w, Int_t(w / (1-y1)));
  fCanvas->SetFillColor(1);
  fCanvas->ToggleEventStatus();
  fCanvas->cd();
  fPad = new TPad("glview", "3DView", 0.0, y1, 1.0, 1.0, 1, 0, 0);
  fPad->Draw();
  
  const char** p = which;
  const char*  m;
  Int_t        n  = 0;
  Int_t        j  = 0;
  while (*(p++)) n++;
  AliFMDDebug(1, ("Got %d buttons", n));
  if (n <= 0) return;

  Double_t yb = 0;
  Double_t xb = 1;
  fCanvas->cd();
  if (TESTBIT(fTreeMask, kDigits) || 
      TESTBIT(fTreeMask, kRawCalib) || 
      TESTBIT(fTreeMask, kRaw)) { 
    yb = .05;
    xb = .66;
    fFactor = new TSlider("pedFactor", "Pedestal Factor", xb+.01, 0, 1, yb);
    fFactor->SetMethod("AliFMDDisplay::Instance()->ChangeFactor()");
    fFactor->SetRange(fInitialFactor, 1);
    fFactor->Draw();
    fFactor->SetMinimum(fInitialFactor);
    TSliderBox *sbox = 
      static_cast<TSliderBox*>(fFactor->GetListOfPrimitives()->
			       FindObject("TSliderBox"));
    if (sbox) { 
      sbox->SetToolTipText("Adjust the noise suppression factor by moving "
			   "lower limit");
    }
  }
  if (TESTBIT(fTreeMask, kHits)    || 
      TESTBIT(fTreeMask, kESD)     || 
      TESTBIT(fTreeMask, kDigits)  || 
      TESTBIT(fTreeMask, kSDigits) || 
      TESTBIT(fTreeMask, kRaw)     ||
      TESTBIT(fTreeMask, kRawCalib)) {
    yb = .05;
    fSlider = new TSlider("genCut", "Multiplicity cut", 0, 0, xb, yb);
    fSlider->SetMethod("AliFMDDisplay::Instance()->ChangeCut()");
    fSlider->SetRange(fInitialMin,fInitialMax);
    fSlider->Draw();
    TSliderBox *sbox = 
      static_cast<TSliderBox*>(fSlider->GetListOfPrimitives()->
			       FindObject("TSliderBox"));
    if (sbox) { 
      sbox->SetToolTipText("Adjust lower and upper limit on data signal");
    }
  }
  // fCanvas->Modified();
  // fCanvas->Update();
  // fCanvas->cd();
  Float_t      x0 = 0;
  Float_t      dx = 1. / n;
  p               = which;
  while ((m = *(p++))) {
    fCanvas->cd();
    AliFMDDebug(1, ("Adding button %s", m));
    TButton* b = new TButton(m, Form("AliFMDDisplay::Instance()->%s()", m),
			     x0, yb, TMath::Min(x0 + dx,.999F), y1);
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
  // Show only the FMD 
  // Do not show 
  // other volumes 
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
  // Bool_t hasFMD1 = kFALSE;
  // Bool_t hasFMD2 = kFALSE;
  // Bool_t hasFMD3 = kFALSE;
  AliFMDDebug(1, ("Getting material FMD_Si$"));
  TGeoMaterial* si   = gGeoManager->GetMaterial("FMD_Si$");      // kRed 
  AliFMDDebug(1, ("Getting material FMD_Carbon$"));
  TGeoMaterial* c    = gGeoManager->GetMaterial("FMD_Carbon$");  // kGray
  AliFMDDebug(1, ("Getting material FMD_Aluminum$"));
  TGeoMaterial* al   = gGeoManager->GetMaterial("FMD_Aluminum$");// kGray-2
  AliFMDDebug(1, ("Getting material FMD_Copper$"));
  TGeoMaterial* cu   = gGeoManager->GetMaterial("FMD_Copper$");	// kGreen-2
  AliFMDDebug(1, ("Getting material FMD_PCB$"));
  TGeoMaterial* pcb  = gGeoManager->GetMaterial("FMD_PCB$");	// kGreen+2
  AliFMDDebug(1, ("Getting material FMD_PCB$"));
  TGeoMaterial* chip = gGeoManager->GetMaterial("FMD_Si Chip$");// kGreen+2
  TObjArray     toshow;
  while ((node = static_cast<TGeoNode*>(next()))) {
    const char* name = node->GetName();
    if (!name) continue;
    if (!(v = node->GetVolume())) continue;

    if (name[0] == 'F') {
      TGeoMaterial* m   = (v->IsAssembly() ? 0 : v->GetMaterial());
      Int_t         col = -1;
      if      (m == si)   col = kRed;
      else if (m == c)	  col = kGray;
      else if (m == al)   col = kYellow+4;  
      else if (m == cu)   col = kRed+6;   
      else if (m == pcb)  col = kGreen+2; 
      else if (m == chip) col = kGreen+4; 
      if (col >= 0) { 
	v->SetLineColor(col);
	v->SetFillColor(col);
      }
      if (name[2] == 'M' && (name[3] == 'T' || name[3] == 'B')) {
	// Virtual Master half-ring volume - top-level
	// Int_t det = node->GetNumber();
	/* switch (det) {
	   case 1: hasFMD1 = true; break;
	   case 2: hasFMD2 = true; break;
	   case 3: hasFMD3 = true; break;
	   default: continue;
	   } */
	toshow.Add(v);
      }
      else if (name[3] == 'V' && (name[2] == 'T' || name[2] == 'B')) 
	toshow.Add(v); // Virtual Half-ring, bare detectors
      else if (name[3] == 'H' && (name[2] == 'F' || name[2] == 'B')) 
	toshow.Add(v); // Virtual Hybrid container 
      else if (name[2] == 'S' && name[3] == 'U') 
	toshow.Add(v); // Virtual support structre 
      // else if (name[3] == 'H' && (name[2] == 'F' || name[2] == 'B')) 
      //  toshow.Add(v); // Virtual Hybrid container 
    }
    v->SetVisibility(kFALSE);
    v->SetVisDaughters(kFALSE);
    v->InvisibleAll(kTRUE);
  }
  TIter i(&toshow);
  while ((v = static_cast<TGeoVolume*>(i()))) {
    if (!v->IsAssembly())
      v->SetVisibility(kTRUE);
    v->InvisibleAll(kFALSE);
    v->SetVisDaughters(kTRUE);
    
  }  
}

    
//____________________________________________________________________
void           
AliFMDDisplay::ExecuteEvent(Int_t event, Int_t px, Int_t py) 
{
  // Execute an event on canvas 
  // Parameters: 
  //   event   What happened
  //   px, py  Pixel coordinates 
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
Bool_t 
AliFMDDisplay::Init()
{
  // Initialize.  GEt transforms and such,
  // so that we can draw thins properly 
  // Returns true on success 
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
  // MAke the aux canvas 
  // This is used to display spectra
  // etc, 
  const Range_t* range = 0;
  if      (TESTBIT(fTreeMask, kESD))      range = &fgkMultRange;
  else if (TESTBIT(fTreeMask, kRawCalib)) range = &fgkMultRange;
  else if (TESTBIT(fTreeMask, kDigits))   range = &fgkAdcRange;
  else if (TESTBIT(fTreeMask, kSDigits))  range = &fgkAdcRange;
  else if (TESTBIT(fTreeMask, kRaw))      range = &fgkAdcRange;
  else if (TESTBIT(fTreeMask, kHits))     range = &fgkEdepRange;
  if (!range) return;
  
  if (!fAux) {
    fAux = new TCanvas(Form("aux_%s", GetName()), 
		       Form("Aux - %s", GetTitle()));
    fAux->SetLogy();
    fAux->SetFillColor(kWhite);
    fAux->SetBorderMode(0);
    fAux->SetBorderSize(0);
    Float_t dBin = (range->fHigh - range->fLow) / range->fNbins;
    fSpec = new TH1D("spec", "Spectra", range->fNbins, 
		     range->fLow-dBin/2, range->fHigh+dBin/2);
    fSpecCut = static_cast<TH1*>(fSpec->Clone("specCut"));
    TString xTitle((TESTBIT(fTreeMask, kRawCalib) || 
		    TESTBIT(fTreeMask, kESD)) ? "#Delta E/#Delta E_{mip}" :
		   (TESTBIT(fTreeMask, kDigits) ||
		    TESTBIT(fTreeMask, kSDigits) || 
		    TESTBIT(fTreeMask, kRaw)) ? "ADC [counts]" :
		   TESTBIT(fTreeMask, kHits) ? "Hits" : "signal");
    fSpec->SetXTitle(xTitle.Data());
    fSpec->SetYTitle("events");
    fSpec->SetFillColor(2);
    fSpec->SetFillStyle(3001);
    fSpecCut->SetXTitle(xTitle.Data());
    fSpecCut->SetYTitle("events");
    fSpecCut->SetFillColor(4);
    fSpecCut->SetFillStyle(3001);
  }
  else {
    fSpec->Reset();
    fSpecCut->Reset();
  }  
}

//____________________________________________________________________
void
AliFMDDisplay::DrawAux()
{
  // Draw in the Aux the canvas 
  // For example draw the spectra 
  // or such stuff 
  if (!fAux) return;
  fAux->cd();
  fAux->Clear();
  fAux->SetLogy(fSpec->GetMaximum() > 10);
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
  // Parameters: 
  //   event   The event number 
  if (!fCanvas) {
    const char* m[] = { "Continue", 
			"Break", 
			"Zoom", 
			"Pick", 
			"Redisplay", 
			"Render", 
			0 }; 
    MakeCanvas(m);
  }
  MakeAux();
  fReturn = kFALSE;
  
  // AliInfo("Clearing canvas");
  // fCanvas->Clear();
  if (!fGeoManager) {
    Warning("End", "No geometry manager");
    return kFALSE;
  }
  AliFMDDebug(1, ("Drawing geometry"));
  fPad->cd();
  fGeoManager->GetTopVolume()->Draw();
  if (fOnlyFMD) ShowOnlyFMD();
  AliFMDDebug(1, ("Adjusting view"));
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
  // Called at of the event. 
  // Draw stuff. 
  // Draw spectrum. 
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
  // Idle loop. 
  // Sends the ROOT loop into the idle loop,
  // so that we can go on. 
  fWait = kTRUE;
  if (fContinous) fTimeout.Start(10, kTRUE);
  while (fWait) {
    gApplication->StartIdleing();
    gSystem->InnerLoop();
    gApplication->StopIdleing();
    if (fContinous) break;
  }
  AliFMDDebug(3, ("After idle loop"));
  if (fMarkers) fMarkers->Delete();
  if (fHits)    fHits->Clear();
  AliFMDDebug(3, ("After clearing caches"));
}

//____________________________________________________________________
Bool_t 
AliFMDDisplay::End()
{
  // End of event.  Draw everything 
  AtEnd();
  Idle();
  if (fReturn) return kFALSE;
  return AliFMDInput::End();
}

//____________________________________________________________________
Int_t
AliFMDDisplay::LookupColor(Float_t x, Float_t min, Float_t max) const
{
  // Look-up color.  
  // Get a colour from the  current palette depending 
  // on the ratio x/max 
  Float_t range  = (max-min);
  Float_t l      = fSlider->GetMinimum();
  Float_t h      = fSlider->GetMaximum();
  if (l == h) { l = 0; h = 1; }
  Float_t cmin   = range * l;
  Float_t cmax   = range * h;
  Float_t crange = (cmax-cmin);
  Int_t   idx    = Int_t((x-cmin) / crange * gStyle->GetNumberOfColors());
  return gStyle->GetColorPalette(idx);
} 

//____________________________________________________________________
void
AliFMDDisplay::SetCut(Float_t l, Float_t h) 
{
  // Change the cut on the slider. 
  fInitialMin = l;
  fInitialMax = h;
  if (!fSlider) return;
  fSlider->SetMinimum(fInitialMin);
  fSlider->SetMaximum(fInitialMax);
  ChangeCut();
}

//____________________________________________________________________
void
AliFMDDisplay::SetFactor(Float_t f) 
{
  // Change the cut on the slider. 
  fInitialFactor = f / 10;
  if (!fFactor) return;
  fFactor->SetMinimum(fInitialFactor);
  ChangeFactor();
}

//____________________________________________________________________
void
AliFMDDisplay::ChangeCut() 
{
  // Change the cut on the slider. 
  // The factor depends on what is 
  // drawn in the AUX canvas
  AliInfo(Form("Range is now %7.5f - %7.5f", fSlider->GetMinimum(), 
	       fSlider->GetMaximum()));
  if ((TESTBIT(fTreeMask, kESD)     || 
       TESTBIT(fTreeMask, kDigits)  || 
       TESTBIT(fTreeMask, kSDigits) || 
       TESTBIT(fTreeMask, kRaw)     ||
       TESTBIT(fTreeMask, kRawCalib))) {
    Float_t l = fSlider->GetMinimum();
    Float_t h = fSlider->GetMaximum();
    l         = 1024 * l + 0;
    h         = 1024 * h + 0;
    AliInfo(Form("ADC range is now %4d - %4d", int(l), int(h)));
  }
  Redisplay();
}
//____________________________________________________________________
void
AliFMDDisplay::ChangeFactor() 
{
  // Change the cut on the slider. 
  // The factor depends on what is 
  // drawn in the AUX canvas
  AliInfo(Form("Noise factor is now %4.1f", 10 * fFactor->GetMinimum()));
  Redisplay();
}

//____________________________________________________________________
void
AliFMDDisplay::Redisplay()
{
  // Redisplay stuff. 
  // Redraw markers, hits, 
  // spectra 
  if (fMarkers) fMarkers->Delete();
  if (fHits)    fHits->Clear();
  if (fSpec)    fSpec->Reset();
  if (fSpecCut) fSpecCut->Reset();
  Event();
  AtEnd();
}
//____________________________________________________________________
void
AliFMDDisplay::Break()
{
  // Redisplay stuff. 
  // Redraw markers, hits, 
  // spectra 
  if (fMarkers) fMarkers->Delete();
  if (fHits)    fHits->Clear();
  if (fSpec)    fSpec->Reset();
  if (fSpecCut) fSpecCut->Reset();
  fReturn = kTRUE;
  fWait   = kFALSE;
}
//____________________________________________________________________
void
AliFMDDisplay::Render()
{
  fPad->cd();
  TVirtualViewer3D* viewer = fPad->GetViewer3D("ogl");
  if (!viewer) return;
}

//____________________________________________________________________
void
AliFMDDisplay::AddMarker(Float_t x, Float_t y, Float_t z,
			 TObject* o, Float_t s, Float_t min, Float_t max)
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
  Float_t  size  = .1;
  Float_t  zsize = (s - min) / (max-min) * 10;
  Float_t  r     = TMath::Sqrt(x * x + y * y);
  Float_t  theta = TMath::ATan2(r, z);
  Float_t  phi   = TMath::ATan2(y, x);
  Float_t  rz    = z + (z < 0 ? 1 : -1) * zsize;
  TMarker3DBox* marker = new  TMarker3DBox(x,y,rz,size,size,zsize,theta,phi);
  if (o) marker->SetRefObject(o);
  marker->SetLineColor(LookupColor(s, min, max));
  fMarkers->Add(marker);
}
//____________________________________________________________________
void
AliFMDDisplay::AddMarker(UShort_t det, Char_t rng, UShort_t sec, UShort_t str,
			 TObject* o, Float_t s, Float_t min, Float_t max)
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
  AddMarker(x,y,z,o,s,min,max);
}
  
//____________________________________________________________________
Bool_t 
AliFMDDisplay::InsideCut(Float_t val, const Float_t& min, 
			 const Float_t& max) const
{
  // 
  // Whether a point is inside 
  // 
  // Parameters:
  //    v   Point
  //    min Minimum
  //    max Maximum
  // 
  // Return:
  //    true if @a v is inside cut 
  //
  Float_t r = max - min;
  Float_t l = fSlider->GetMinimum();
  Float_t h = fSlider->GetMaximum();
  if (l == h) { l = 0; h = 1; }
  if (val < r * l + min || val > r * h + min) { 
    AliFMDDebug(2, ("Value %f is outside cut %f - %f (range %f - %f)", 
		    val, min+r*l, min+r*h, min, max));
    return kFALSE;
  }
  return kTRUE;
}


//____________________________________________________________________
Bool_t 
AliFMDDisplay::ProcessHit(AliFMDHit* hit, TParticle* /* p */) 
{
  // Process a hit. 
  // Parameters: 
  //   hit   Hit data
  static const Float_t rMin  = fgkEdepRange.fLow;
  static const Float_t rMax  = fgkEdepRange.fHigh;

  if (!hit) { AliError("No hit");   return kFALSE; }
  // if (!p)   { AliError("No track"); return kFALSE; }
  Float_t  edep  = hit->Edep();

  if (fHits)                        fHits->Add(hit);
  if (fSpec)                        fSpec->Fill(edep);
  if (!InsideCut(edep, rMin, rMax)) return kTRUE;
  if (fSpecCut)                     fSpecCut->Fill(edep);
  
  AddMarker(hit->X(), hit->Y(), hit->Z(), hit, edep, rMin, rMax);
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDDisplay::ProcessDigit(AliFMDDigit* digit)
{
  // Process a digit 
  // Parameters: 
  //   digit Digit information 
  static const Float_t rMin  = fgkAdcRange.fLow;
  static const Float_t rMax  = fgkAdcRange.fHigh;

  if (!digit) { AliError("No digit");   return kFALSE; }

  UShort_t det           =  digit->Detector();
  Char_t   ring          =  digit->Ring();
  UShort_t sec           =  digit->Sector();
  UShort_t str           =  digit->Strip();
  Double_t threshold     =  GetADCThreshold(det, ring, sec, str);
  Float_t  counts        =  digit->Counts();
  if (threshold < 0) { counts += -threshold; threshold = 0; }

  AliFMDDebug(10, ("FMD%d%c[%02d,%03d] counts %4d threshold %4d", 
		   det, ring, sec, str, Int_t(counts), Int_t(threshold)));
  if (fHits)                                    fHits->Add(digit);
  if (fSpec)                                    fSpec->Fill(counts);
  if (!InsideCut(counts-threshold, rMin, rMax)) return kTRUE;
  if (fSpecCut)                                 fSpecCut->Fill(counts);
  

  AddMarker(det, ring, sec, str, digit, counts, rMin, rMax);
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDDisplay::ProcessSDigit(AliFMDSDigit* sdigit)
{
  // Process a sdigit 
  // Parameters: 
  //   sdigit Digit information 
  static const Float_t rMin  = fgkAdcRange.fLow;
  static const Float_t rMax  = fgkAdcRange.fHigh;

  if (!sdigit) { AliError("No sdigit");   return kFALSE; }
  
  UShort_t det           =  sdigit->Detector();
  Char_t   ring          =  sdigit->Ring();
  UShort_t sec           =  sdigit->Sector();
  UShort_t str           =  sdigit->Strip();
  Float_t  counts        =  sdigit->Counts();

  if (fHits)                          fHits->Add(sdigit);
  if (fSpec)                          fSpec->Fill(counts);
  if (!InsideCut(counts, rMin, rMax)) return kTRUE;
  if (fSpecCut)                       fSpecCut->Fill(counts);
  
  AddMarker(det, ring, sec, str, sdigit, counts, rMin, rMax);
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDDisplay::ProcessRawDigit(AliFMDDigit* digit)
{
  // PRocess raw data 
  // Parameters: 
  //   digit Digit information 
  AliFMDDebug(50, ("Forwarding call of ProcessRaw to ProcessDigit "
		  "for FMD%d%c[%02d,%03d] %d", 
		  digit->Detector(), digit->Ring(), digit->Sector(), 
		  digit->Strip(), digit->Counts()));
  return ProcessDigit(digit);
}

//____________________________________________________________________
Bool_t 
AliFMDDisplay::ProcessRawCalibDigit(AliFMDDigit* digit)
{
  // Process a digit 
  // Parameters: 
  //   digit Digit information 
  static const Float_t rMin  = fgkMultRange.fLow;
  static const Float_t rMax  = fgkMultRange.fHigh;
  static const Float_t aMin  = fgkAdcRange.fLow;
  static const Float_t aMax  = fgkAdcRange.fHigh;

  if (!digit) { AliError("No digit");   return kFALSE; }

  AliFMDParameters* parm = AliFMDParameters::Instance();
  UShort_t det           =  digit->Detector();
  Char_t   ring          =  digit->Ring();
  UShort_t sec           =  digit->Sector();
  UShort_t str           =  digit->Strip();
  Double_t gain          =  parm->GetPulseGain(det, ring, sec, str);
  Double_t ped           =  parm->GetPedestal(det, ring, sec, str);
  Double_t threshold     =  GetADCThreshold(det, ring, sec, str);
  Float_t  counts        =  digit->Counts();
  if (threshold < 0) { counts += -threshold; threshold = 0; ped = 0; }

  //   Double_t edep = ((counts * parm->GetEdepMip() / 
  // 		    (gain * parm->GetDACPerMIP()));
  Double_t mult = (counts-ped) / (gain * parm->GetDACPerMIP());
  if (gain < 0.1 || gain > 10) mult = 0;
  

  AliFMDDebug(10, ("FMD%d%c[%02d,%03d] adc %4d "
		   "(threshold %4d, gain %6.3f) -> mult %7.4f", 
		   det, ring, sec, str, int(counts), int(threshold), 
		   gain, mult));
  if (fHits)                                    fHits->Add(digit);
  if (fSpec)                                    fSpec->Fill(mult);
  if (!InsideCut(counts-threshold, aMin, aMax)) return kTRUE;
  if (fSpecCut)                                 fSpecCut->Fill(mult);
  
  if (mult >= 0) 
    AddMarker(det, ring, sec, str, digit, mult, rMin, rMax);
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDDisplay::ProcessRecPoint(AliFMDRecPoint* recpoint)
{
  // Process reconstructed point 
  // Parameters: 
  //  recpoint  Reconstructed multiplicity/energy
  static const Float_t rMin  = fgkMultRange.fLow;
  static const Float_t rMax  = fgkMultRange.fHigh;


  if (!recpoint) { AliError("No recpoint");   return kFALSE; }

  if (!InsideCut(recpoint->Particles(), rMin, rMax)) return kTRUE;

  if (fHits) fHits->Add(recpoint);
  AddMarker(recpoint->Detector(), recpoint->Ring(), recpoint->Sector(),  
	    recpoint->Strip(), recpoint, recpoint->Particles(), rMin, rMax);
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDDisplay::ProcessESD(UShort_t det, Char_t rng, UShort_t sec, UShort_t str,
			  Float_t, Float_t mult)
{
  // Process data from ESD 
  // Parameters 
  //   det,rng,sec,str   Detector coordinates. 
  //   mult              Multiplicity. 
  static const Float_t rMin  = fgkMultRange.fLow;
  static const Float_t rMax  = fgkMultRange.fHigh;
  
  Double_t cmult = mult;
  if (fSpec) fSpec->Fill(cmult);
  if (!InsideCut(cmult, rMin, rMax) || cmult == AliESDFMD::kInvalidMult) 
    return kTRUE;

  AddMarker(det,rng,sec,str, 0, cmult, rMin, rMax);

  if (fSpecCut) fSpecCut->Fill(cmult);

  return kTRUE;
}

//____________________________________________________________________
Double_t
AliFMDDisplay::GetADCThreshold(UShort_t d, Char_t r, 
			       UShort_t s, UShort_t t) const
{
  // 
  // Get the ADC threshold
  // 
  // Parameters:
  //    d Detector
  //    r Ring 
  //    s Sector 
  //    t Strip 
  // 
  // Return:
  //    The threshold 
  //
  AliFMDParameters* parm = AliFMDParameters::Instance();
  Double_t ped           =  parm->GetPedestal(d,r, s, t);
  Double_t pedW          =  parm->GetPedestalWidth(d,r, s, t);
  Double_t threshold     = 0;
  if (fFMDReader && fFMDReader->IsZeroSuppressed(d-1))
    threshold = - fFMDReader->NoiseFactor(d-1) * pedW;
  else 
    threshold = ped + pedW * 10 * fFactor->GetMinimum();
  AliFMDDebug(10, ("FMD%d%c[%2d,%3d] ped: %f +/- %f [factor: %f-%f]", 
		   d, r, s, t, ped, pedW, 
		   fFactor->GetMinimum(), fFactor->GetMaximum()));
  if (threshold > fgkAdcRange.fHigh) threshold = fgkAdcRange.fHigh;
  return threshold;
}

//____________________________________________________________________
//
// EOF
//

