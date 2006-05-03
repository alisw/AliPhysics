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
AliFMDDisplay::AliFMDDisplay(const char* gAliceFile)
  : AliFMDInput(gAliceFile),
    fWait(kFALSE),
    fCanvas(0), 
    fPad(0), 
    fButton(0), 
    fZoom(0),
    fPick(0),
    fZoomMode(kFALSE)
{
  // Constructor of an FMD display object. 
  AddLoad(kGeometry);
  fMarkers = new TObjArray;
  fHits    = new TObjArray;
  fMarkers->SetOwner(kTRUE);
  fHits->SetOwner(kFALSE);
  fgInstance = this;
  SetMultiplicityCut();
  SetPedestalFactor();
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
  // AliFMDParameters* parm = AliFMDParameters::Instance();
  // parm->Init();
  return kTRUE;
}
//____________________________________________________________________
Bool_t 
AliFMDDisplay::Begin(Int_t event) 
{
  // Begin of event.  Make canvas is not already done 
  if (!fCanvas) {
    gStyle->SetPalette(1);
    fCanvas = new TCanvas("display", "Display", 700, 700);
    fCanvas->SetFillColor(1);
    fCanvas->ToggleEventStatus();
    fPad = new TPad("view3D", "3DView", 0.0, 0.05, 1.0, 1.0, 1, 0, 0);
    fCanvas->cd();
    fPad->Draw();
  }
  if (!fButton) {
    fCanvas->cd();
    fButton = new TButton("Continue", "AliFMDDisplay::Instance()->Continue()",
			  0, 0, .5, .05);
    fButton->Draw();
    fZoom = new TButton("Zoom", "AliFMDDisplay::Instance()->Zoom()",
			.5, 0, .75, .05);
    fZoom->Draw();
    fPick = new TButton("Pick", "AliFMDDisplay::Instance()->Pick()",
			.75, 0, 1, .05);
    fPick->Draw();
  }
  AliInfo("Clearing canvas");
  // fCanvas->Clear();
  if (!fGeoManager) {
    Warning("End", "No geometry manager");
    return kFALSE;
  }
  AliInfo("Drawing geometry");
  fPad->cd();
  fGeoManager->GetTopVolume()->Draw();
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
Bool_t 
AliFMDDisplay::End()
{
  // End of event.  Draw everything 
  fPad->cd();
  fMarkers->Draw();
  fPad->cd();
  AppendPad();
  // fPad->Update();
  fPad->cd();
  // fCanvas->Modified(kTRUE);
  //fCanvas->Update();
  // fCanvas->cd();
  // fPad->cd();
  fWait = kTRUE;
  while (fWait) {
    gApplication->StartIdleing();
    gSystem->InnerLoop();
    gApplication->StopIdleing();
  }
  AliInfo("After idle loop");
  fMarkers->Delete();
  fHits->Clear();
  AliInfo("After clearing caches");
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
AliFMDDisplay::ProcessHit(AliFMDHit* hit, TParticle* p) 
{
  // Process a hit 
  if (!hit) { AliError("No hit");   return kFALSE; }
  if (!p)   { AliError("No track"); return kFALSE; }

  fHits->Add(hit);
  Float_t  size  = .1;
  Float_t  zsize = hit->Edep() * 10;
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
  Double_t threshold     =  ped * fPedestalFactor * pedW;
  Float_t  counts        = digit->Counts();
  if (counts < threshold) return kTRUE;
  fHits->Add(digit);

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
  fHits->Add(recpoint);
  AddMarker(recpoint->Detector(), recpoint->Ring(), recpoint->Sector(),  
	    recpoint->Strip(), recpoint, recpoint->Particles(), 20);
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDDisplay::ProcessESD(AliESDFMD* esd)
{
  // Process event summary data
  for (UShort_t det = 1; det <= 3; det++) {
    Char_t rings[] = { 'I', (det == 1 ? '\0' : 'O'), '\0' };
    for (Char_t* rng = rings; *rng != '\0'; rng++) {
      UShort_t nsec = (*rng == 'I' ?  20 :  40);
      UShort_t nstr = (*rng == 'O' ? 512 : 256);
      for (UShort_t sec = 0; sec < nsec; sec++) {
	for (UShort_t str = 0; str < nstr; str++) {
	  Float_t mult = esd->Multiplicity(det,*rng,sec,str);
	  if (mult < fMultCut) continue;
	  AddMarker(det,*rng,sec,str, 0, mult, 20);
	}
      }
    }
  }
  return kTRUE;
}

//____________________________________________________________________
//
// EOF
//
