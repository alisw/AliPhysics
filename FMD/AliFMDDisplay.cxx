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
#include <TArrayF.h>
#include <TMarker3DBox.h>
#include <TGeoManager.h>
#include <TMath.h>
#include <TApplication.h>
#include <TButton.h>
#include <TParticle.h>
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
  Int_t idx = Int_t(x / max * gStyle->GetNumberOfColors());
  return gStyle->GetColorPalette(idx);
}


//____________________________________________________________________
Bool_t 
AliFMDDisplay::ProcessHit(AliFMDHit* hit, TParticle* p) 
{
  if (!hit) { AliError("No hit");   return kFALSE; }
  if (!p)   { AliError("No track"); return kFALSE; }

  fHits->Add(hit);
  Float_t  size  = .1;
  Float_t  pt    = TMath::Sqrt(hit->Py()*hit->Py()+hit->Px()*hit->Px());
  Float_t  theta = TMath::ATan2(pt, hit->Pz());
  Float_t  phi   = TMath::ATan2(hit->Py(), hit->Px());
  TMarker3DBox* marker = new  TMarker3DBox(hit->X(), hit->Y(), hit->Z(),
					   size, size, size, theta, phi);
  marker->SetLineColor(LookupColor(hit->Edep(), 1));
  marker->SetRefObject(hit);
  fMarkers->Add(marker);
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDDisplay::ProcessDigit(AliFMDDigit* digit)
{
  if (!digit) { AliError("No digit");   return kFALSE; }

  Double_t x, y, z;
  AliFMDGeometry*   geom = AliFMDGeometry::Instance();
  AliFMDParameters* parm = AliFMDParameters::Instance();
  Double_t threshold = (parm->GetPedestal(digit->Detector(), 
					  digit->Ring(), 
					  digit->Sector(), 
					  digit->Strip())
			+ 4 * parm->GetPedestalWidth(digit->Detector(), 
						     digit->Ring(), 
						     digit->Sector(), 
						     digit->Strip()));
  if (digit->Counts() < threshold) return kTRUE;
  fHits->Add(digit);
  geom->Detector2XYZ(digit->Detector(), digit->Ring(), digit->Sector(), 
		    digit->Strip(), x, y, z);
  Float_t  size  = .1;
  Float_t  r     = TMath::Sqrt(x * x + y * y);
  Float_t  theta = TMath::ATan2(r, z);
  Float_t  phi   = TMath::ATan2(y, x);
  TMarker3DBox* marker = new  TMarker3DBox(x,y,z,size,size,size,theta,phi);
  marker->SetRefObject(digit);
  marker->SetLineColor(LookupColor(digit->Counts(), 1024));
  fMarkers->Add(marker);
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDDisplay::ProcessRaw(AliFMDDigit* digit)
{
  return ProcessDigit(digit);
}

//____________________________________________________________________
Bool_t 
AliFMDDisplay::ProcessRecPoint(AliFMDRecPoint* recpoint)
{
  if (!recpoint) { AliError("No recpoint");   return kFALSE; }
  if (recpoint->Particles() < .1) return kTRUE;
  fHits->Add(recpoint);
  Double_t x, y, z;
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  geom->Detector2XYZ(recpoint->Detector(), recpoint->Ring(), 
		    recpoint->Sector(),  recpoint->Strip(), x, y, z);

  Float_t  size  = .1;
  Float_t  r     = TMath::Sqrt(x * x + y * y);
  Float_t  theta = TMath::ATan2(r, z);
  Float_t  phi   = TMath::ATan2(y, x);
  TMarker3DBox* marker = new  TMarker3DBox(x,y,z,size,size,size,theta,phi);
  marker->SetRefObject(recpoint);
  marker->SetLineColor(LookupColor(recpoint->Particles(), 20));
  fMarkers->Add(marker);
  return kTRUE;
}

//____________________________________________________________________
//
// EOF
//
