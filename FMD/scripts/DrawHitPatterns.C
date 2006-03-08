//
// $Id$
//
// Script that contains a class to draw hits, using the
// AliFMDInputHits class in the util library. 
//
// It draws the energy loss versus the p/(mq^2).  It can be overlayed
// with the Bethe-Bloc curve to show how the simulation behaves
// relative to the expected. 
//
// Use the script `Compile.C' to compile this class using ACLic. 
//
#include <TH2D.h>
#include <AliFMDHit.h>
#include <AliFMDInput.h>
#include <iostream>
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

class DrawHitPatterns : public AliFMDInputHits
{
private:
  Bool_t     fWait;
  Bool_t     fPrimary;
  Int_t      fPdg;
  TObjArray* fMarkers;
  TObjArray* fHits;
  TCanvas*   fCanvas;
  TPad*      fPad;
  TButton*   fButton;
  TButton*   fZoom;
  TButton*   fPick;
  Bool_t     fZoomMode;
public:
  static DrawHitPatterns* fgInstance;
  DrawHitPatterns(Bool_t primary=kTRUE, Int_t pdg=211) 
    : fPrimary(primary), fPdg(pdg), 
      fCanvas(0), fPad(0), fButton(0)
  { 
    AddLoad(kKinematics);
    AddLoad(kGeometry);
    fMarkers = new TObjArray;
    fHits    = new TObjArray;
    fMarkers->SetOwner(kTRUE);
    fHits->SetOwner(kFALSE);
    fgInstance = this;
  }
  void Continue() { fWait = kFALSE; }
  void Zoom() { fZoomMode = kTRUE; }
  void Pick() { fZoomMode = kFALSE; }
  void ExecuteEvent(Int_t event, Int_t px, Int_t py) 
  {
    Info("ExecuteEvent", "Event %d, at (%d,%d)", px, py);
    static Float_t x0, y0, x1, y1;
    static Int_t   px0, py0, pxOld,pyOld;
    static Bool_t  lineDraw;
    if (px == 0 && py == 0) return;
    if (!fZoomMode && fPad->GetView()) {
      fPad->GetView()->ExecuteRotateView(event, px, py);
      return;
    }
    gPad->SetCursor(kCross);
    switch (event) {
    case kButton1Down: 
      gPad->TAttLine::Modify();
      x0 = fPad->AbsPixeltoX(px);
      y0 = fPad->AbsPixeltoY(py);
      px0 = pxOld = px;
      py0 = pyOld = py;
      lineDraw = kFALSE;
      return;
    case kButton1Motion:
      if (lineDraw) 
	gVirtualX->DrawBox(px0, py0, pxOld, pyOld, TVirtualX::kHollow);
      pxOld = px;
      pyOld = py;
      lineDraw = kTRUE;
      gVirtualX->DrawBox(px0, py0, pxOld, pyOld, TVirtualX::kHollow);
      return;
    case kButton1Up:
      fPad->GetCanvas()->FeedbackMode(kFALSE);
      if (px == px0 || py == py0) return;
      x1 = gPad->AbsPixeltoX(px);
      y1 = gPad->AbsPixeltoY(py);
      if (x1 < x0) std::swap(x0, x1); 
      if (y1 < y0) std::swap(y0, y1); 
      gPad->Range(x0, y0, x1, y1);
      gPad->Modified();
      return;
    }
  }
  Int_t  DistanceToPrimitive(Int_t px, Int_t py) 
  {
    Info("DistanceToPrimitive", "@ (%d,%d)", px, py);
    gPad->SetCursor(kCross);
    Float_t xmin = gPad->GetX1();
    Float_t xmax = gPad->GetX2();
    Float_t dx   = .02 * (xmax - xmin);
    Float_t x    = gPad->AbsPixeltoX(px);
    if (x < xmin + dx || x > xmax - dx) return 9999;
    return (fZoomMode ? 0 : 7);
  }
  Bool_t Begin(Int_t event) 
  {
    if (!fCanvas) {
      fCanvas = new TCanvas("display", "Display", 700, 700);
      fCanvas->SetFillColor(1);
      fCanvas->ToggleEventStatus();
      fPad = new TPad("view3D", "3DView", 0.0, 0.1, 1.0, 1.0, 1, 0, 0);
      fCanvas->cd();
      fPad->Draw();
    }
    if (!fButton) {
      fCanvas->cd();
      fButton = new TButton("Continue", 
			    "DrawHitPatterns::fgInstance->Continue()",
			    0, 0, .5, .05);
      fButton->Draw();
      fZoom = new TButton("Zoom", 
			  "DrawHitPatterns::fgInstance->Zoom()",
			  .5, 0, .75, .05);
      fZoom->Draw();
      fPick = new TButton("Pick", 
			  "DrawHitPatterns::fgInstance->Pick()",
			  .75, 0, 1, .05);
      fPick->Draw();
    }
    Info("Begin", "Clearing canvas");
    // fCanvas->Clear();
    if (!fGeoManager) {
      Warning("End", "No geometry manager");
      return kFALSE;
    }
    Info("Begin", "Drawing geometry");
    fPad->cd();
    fGeoManager->GetTopVolume()->Draw();
    Info("Begin", "Adjusting view");
    Int_t irep;
    if (fPad->GetView()) {
      fPad->GetView()->SetView(-200, -40, 80, irep);
      fPad->GetView()->Zoom();
      fPad->Modified();
      fPad->cd();
    }
    Info("Begin", "Drawing button");

    return AliFMDInputHits::Begin(event);
  }
    
  Bool_t ProcessHit(AliFMDHit* hit, TParticle* p) 
  {
    if (!hit) {
      std::cout << "No hit" << std::endl;
      return kFALSE;
    }

    if (!p) {
      std::cout << "No track" << std::endl;
      return kFALSE;
    }
    if (fPrimary && !p->IsPrimary()) return kTRUE;
    if (hit->IsStop()) return kTRUE;
    if (fPdg != 0) {
      if (TMath::Abs(hit->Pdg()) != fPdg) return kTRUE;
    }
    fHits->Add(hit);
    Float_t  size  = TMath::Min(TMath::Max(hit->Edep() * .1, .1), 1.);
    Float_t  pt    = TMath::Sqrt(hit->Py()*hit->Py()+hit->Px()*hit->Px());
    Float_t  theta = TMath::ATan2(pt, hit->Pz());
    Float_t  phi   = TMath::ATan2(hit->Py(), hit->Px());
    TMarker3DBox* marker = new  TMarker3DBox(hit->X(), hit->Y(), hit->Z(),
					     size, size, size,
					     theta, phi);
    marker->SetLineColor((p->IsPrimary() ? 3 : 4));
    marker->SetRefObject(hit);
    fMarkers->Add(marker);
    return kTRUE;
  }
  Bool_t End() 
  {
    // fGeoManager->GetTopVolume()->Draw();
    fPad->cd();
    fMarkers->Draw();
    AppendPad();
    fPad->Modified(kTRUE);
    fPad->Update();
    fPad->cd();
    fCanvas->Modified(kTRUE);
    fCanvas->Update();
    fCanvas->cd();
    
    fWait = kTRUE;
    while (fWait) {
      gApplication->StartIdleing();
      gSystem->InnerLoop();
      gApplication->StopIdleing();
    }
    Info("End", "After idle loop");
    fMarkers->Delete();
    fHits->Clear();
    Info("End", "After clearing caches");
    return AliFMDInputHits::End();
  }
  ClassDef(DrawHitPatterns,0);
};

DrawHitPatterns* DrawHitPatterns::fgInstance = 0;

//____________________________________________________________________
//
// EOF
//
