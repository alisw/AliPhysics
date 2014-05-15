/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Canvas for ALICE figures
//
// Author: Jochen Klein <jochen.klein@cern.ch>

#ifndef ALIFIGURE_H
#define ALIFIGURE_H

#include "TCanvas.h"
#include "TColor.h"

class TObject;
class TString;
class TCanvas;
class TText;
class TASImage;
class TStyle;
class TLatex;

class AliFigure : public TCanvas
{
 public:
  AliFigure(const char* name = "", const char* title = "", Int_t ww = 800, Int_t wh = 600);
  ~AliFigure();

  enum Status_t {
    kWorkInProgress = 0,
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

  static TStyle* Style();

 protected:
  void UpdateLogo();
  void UpdateLogoPos();
  void UpdatePad(TPad *pad);

  Bool_t   isMC;
  Status_t fStatus;

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
  AliFigure(const AliFigure& rhs); // not implemented
  AliFigure& operator=(const AliFigure& rhs); // not implemented

  ClassDef(AliFigure, 1);
};

#endif
