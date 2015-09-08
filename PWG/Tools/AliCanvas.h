/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Canvas for ALICE figures
//
// Author: Jochen Klein <jochen.klein@cern.ch>

#ifndef ALICANVAS_H
#define ALICANVAS_H

#include "TCanvas.h"
#include "TColor.h"

class TObject;
class TString;
class TCanvas;
class TText;
class TASImage;
class TStyle;
class TLatex;

class AliCanvas : public TCanvas
{
 public:
  AliCanvas(const char* name = "", const char* title = "", Int_t ww = 800, Int_t wh = 600);
  ~AliCanvas();

  enum Status_t {
    kWorkInProgress = 0,
    kThisWork,
    kPerformance,
    kPreliminary,
    kFinal,
    kStatusLast
  };

  enum Mode_t {
    kPresentation,
    kPrint
  };

  enum Pos_t {
    kN = 0, kNE, kE, kSE, kS, kSW, kW, kNW, kCenter
  };

  void Draw(Option_t *option = "");
  void Paint(Option_t *option = "");
  void Clear(Option_t *option = "");
  void Update();

  void SetStatus(Status_t status);

  void SetLogoFilename(Status_t status, TString filename);
  void SetLogoPos(Float_t x, Float_t y);
  void SetLogoPos(Pos_t pos);
  void SetLogoSize(Float_t size);

  void SetCollSystem(TString txt);
  void SetCollSystemPos(Float_t x, Float_t y);

  void SetDataSample(TString txt);

  void SetTextSize(Float_t size);

  void SetDisabled(Bool_t disable = kTRUE) { fDisabled = disable; }
  Bool_t GetDisabled() const { return fDisabled; }

  void SetDisabledMargin(Bool_t disable = kTRUE) { fDisabledMargin = disable; }
  Bool_t GetDisabledMargin() const { return fDisabledMargin; }

  static TStyle* Style();

 protected:
  void UpdateLogo();
  void UpdateLogoPos();
  void UpdatePad(TPad *pad);

  Bool_t   fDisabled; // updating of canvas disabled
  Bool_t   fDisabledMargin; // updating of canvas margins disabled

  Bool_t   isMC; // set if Monte-Carlo figure
  Status_t fStatus; // status of plot (published, preliminary, ...)

  TString  fDateFormat;
  Float_t  fTextSize;
  Color_t  fTextColor;

  TString  fLogoFilename[kStatusLast];
  Float_t  fLogoPosX;
  Float_t  fLogoPosY;
  Float_t  fLogoHeight;

  TPad     *fLogoPad;
  TText    *fDate;
  TText    *fStatusPad;
  TLatex   *fCollSystem;
  TLatex   *fDataSample;
  TText    *fTag;
  TASImage *fLogo;

  const char* fStatusString[kStatusLast];

  static TStyle *fgStyle;

 private:
  AliCanvas(const AliCanvas& rhs); // not implemented
  AliCanvas& operator=(const AliCanvas& rhs); // not implemented

  ClassDef(AliCanvas, 1);
};

#endif
