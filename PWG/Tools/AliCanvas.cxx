/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

////////////////////////////////////////////////////////////////////////////
//
// Generic canvas for ALICE figures
//
// This class inherits from TCanvas and can be used instead of it. It
// places the ALICE logo, and other things according to the figure status.
//
// Author: Jochen Klein <jochen.klein@cern.ch>
//
////////////////////////////////////////////////////////////////////////////

#include "TCanvas.h"
#include "TList.h"
#include "TDatime.h"

#include "TStyle.h"
#include "TColor.h"
#include "TText.h"
#include "TASImage.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "AliCanvas.h"

TStyle *AliCanvas::fgStyle;

AliCanvas::AliCanvas(const char* name, const char* title, Int_t ww, Int_t wh) :
  TCanvas(name, title, ww, wh),
  fDisabled(kFALSE),
  fDisabledMargin(kFALSE),
  isMC(kFALSE),
  fStatus(kWorkInProgress),
  fDateFormat("%d/%m/%Y"),
  fTextSize(20),
  fTextColor(kBlack),
  fLogoFilename(),
  fLogoPosX(0.2),
  fLogoPosY(0.5),
  fLogoHeight(100),
  fLogoPad(0x0),
  fDate(0x0),
  fStatusPad(0x0),
  fCollSystem(0x0),
  fDataSample(0x0),
  fTag(0x0),
  fLogo(0x0)
{
  // default ctor

  fStatusString[kWorkInProgress] = "work in progress";
  fStatusString[kThisWork]       = "- this work -";
  fStatusString[kPerformance]    = "Performance";
  fStatusString[kPreliminary]    = "Preliminary";
  fStatusString[kFinal]          = "";

  fLogoFilename[kWorkInProgress] = "";
  fLogoFilename[kThisWork]       = "";
  fLogoFilename[kPerformance]    = "$ALICE_PHYSICS/PWG/Tools/fig/2011-Nov-24-ALICE_logo_WithoutStrapline.eps";
  // fLogoFilename[kPerformance]    = "$ALICE_PHYSICS/PWG/Tools/fig/2011-Nov-24-ALICE_PERFORMANCE_logo_BLACK_small_usage_design.eps";
  fLogoFilename[kPreliminary]    = "$ALICE_PHYSICS/PWG/Tools/fig/2011-Nov-24-ALICE_PRELIMINARY_logo_BLACK_small_usage_design.eps";
  fLogoFilename[kFinal]          = "";

  this->SetLeftMargin(0.15);
  this->SetRightMargin(0.04);
  this->SetRightMargin(0.04);
  this->SetBottomMargin(0.15);

  TDatime now;
  char date[30];
  time_t t = (time_t) now.Convert();
  struct tm* loctis = localtime(&t);
  strftime(date, 30, fDateFormat.Data(), loctis);

  fCollSystem = new TLatex(0.96, 0.95, "pp #sqrt{s} = 7 TeV");
  fCollSystem->SetNDC();
  fCollSystem->SetTextSize(fTextSize);
  fCollSystem->SetTextFont(43);
  fCollSystem->SetTextAlign(33);
  fCollSystem->SetTextColor(fTextColor);

  fDataSample = new TLatex(0.15, 0.05, "Run ######");
  fDataSample->SetNDC();
  fDataSample->SetTextSize(fTextSize);
  fDataSample->SetTextFont(43);
  fDataSample->SetTextAlign(11);
  fDataSample->SetTextColor(fTextColor);

  fDate    = new TText(0.5, 0., date);
  fDate->SetNDC();
  fDate->SetTextSize(fTextSize);
  fDate->SetTextFont(43);
  fDate->SetTextAlign(22);
  fDate->SetTextColor(fTextColor);

  fStatusPad  = new TText(0.5, 0.15, fStatusString[fStatus]);
  fStatusPad->SetNDC();
  fStatusPad->SetTextSize(fTextSize);
  fStatusPad->SetTextFont(43);
  fStatusPad->SetTextAlign(22);
  fStatusPad->SetTextColor(fTextColor);

  fTag = new TText(0.5, 0., "");
  fTag->SetNDC();
  fTag->SetTextSize(fTextSize);
  fTag->SetTextFont(43);
  fTag->SetTextAlign(22);
  fTag->SetTextColor(fTextColor);

//   fLogoPad->ResetBit(kCanDelete);
//   fCollSystem->ResetBit(kCanDelete);
//   fStatusPad->ResetBit(kCanDelete);
//   fDate->ResetBit(kCanDelete);

//   GetListOfPrimitives()->SetOwner(kFALSE);

  // now compose the logo with the additional text
  fLogoPad = new TPad("aliceLogo", "Pad for ALICE Logo", .2, .2, .4, .4);
  this->UpdateLogo();

  if (fLogo)
    fLogoPad->GetListOfPrimitives()->Add(fLogo);

  UpdateLogoPos();
}

AliCanvas::~AliCanvas()
{
  // destructor

};

void AliCanvas::Draw(Option_t *option)
{
  // just passing on to the TCanvas method

  TCanvas::Draw(option);
}

void AliCanvas::Paint(Option_t *option)
{
  // just passing on to the TCanvas method

  TCanvas::Paint(option);
}

void AliCanvas::Clear(Option_t *option)
{
  // clean up

  // now walk through all the children
  TList *listOfPrimitives = this->GetListOfPrimitives();
  TIter iter(listOfPrimitives);

  TObject *obj = 0x0;
  while ((obj = iter())) {
    if ((obj == fCollSystem) || (obj == fStatusPad) || (obj == fLogoPad) || (obj == fDate) || (obj == fTag) || (obj == fDataSample)) {
      listOfPrimitives->Remove(obj);
    }
  }

  TCanvas::Clear(option);
}

void AliCanvas::Update()
{
  // while updating make sure that the additional logo and text
  // get drawn last

  if (fDisabled)
    return;

  this->UpdateLogoPos();
  this->UpdatePad(this);

  if (fLogo && fLogo->IsValid()) {
    this->GetListOfPrimitives()->Add(fLogoPad);
  }
  this->GetListOfPrimitives()->Add(fCollSystem);
  if ((fStatus == kWorkInProgress) || (fStatus == kThisWork))
    this->GetListOfPrimitives()->Add(fStatusPad);
  this->GetListOfPrimitives()->Add(fDataSample);
  // this->GetListOfPrimitives()->Add(fDate);
  this->GetListOfPrimitives()->Add(fTag);

  this->Modified(kTRUE);
  TCanvas::Update();
}

void AliCanvas::SetStatus(Status_t status)
{
  // set the figure status, i.e. "work in progress", "performance",
  // "preliminary", or "final"

  fStatus = status;
  fStatusPad->SetText(fStatusPad->GetX(), fStatusPad->GetY(), fStatusString[fStatus]);

  UpdateLogo();

  this->Modified(kTRUE);
}

void AliCanvas::SetLogoFilename(Status_t status, TString filename)
{
  // set the filename for the logo for the given figure status

  fLogoFilename[status] = filename;
}

void AliCanvas::UpdateLogo()
{
  // update the logo

  if (fLogo)
    fLogoPad->GetListOfPrimitives()->Remove(fLogo);
  delete fLogo;

  if (fLogoFilename[fStatus].Length() > 0) {
    fLogo    = new TASImage(fLogoFilename[fStatus]);
    if (fLogo->IsValid()) {
      fLogo->SetImageQuality(TASImage::kImgBest);
      // fLogo->Crop(70, 93, 634, 600);
    }
  }
  else
    fLogo = 0x0;

  this->UpdateLogoPos();
}

void AliCanvas::SetLogoPos(Float_t x, Float_t y)
{
  // set the position of the logo

  fLogoPosX = x;
  fLogoPosY = y;

  this->UpdateLogoPos();
}

void AliCanvas::SetLogoPos(Pos_t pos)
{
  // set the logo position
  // as north (kN), north-east (kNE), ...

  switch (pos) {
  case kN:
    SetLogoPos(0.5, 0.8);
    break;
  case kNE:
    SetLogoPos(0.8, 0.8);
    break;
  case kE:
    SetLogoPos(0.8, 0.5);
    break;
  case kSE:
    SetLogoPos(0.8, 0.2);
    break;
  case kS:
    SetLogoPos(0.5, 0.2);
    break;
  case kSW:
    SetLogoPos(0.25, 0.2);
    break;
  case kW:
    SetLogoPos(0.25, 0.5);
    break;
  case kNW:
    SetLogoPos(0.25, 0.8);
    break;
  case kCenter:
    SetLogoPos(0.5, 0.5);
    break;
  default:
    this->Error(__FUNCTION__, "Unknown position specifier");
  }
}

void AliCanvas::SetLogoSize(Float_t size)
{
  // set the height of the logo as a fraction of the canvas height

  fLogoHeight = this->GetWh() * size;

  this->UpdateLogoPos();
}

void AliCanvas::SetCollSystem(TString txt)
{
  // set the text for the collision system

  fCollSystem->SetText(fCollSystem->GetX(), fCollSystem->GetY(), txt);
}

void AliCanvas::SetCollSystemPos(Float_t x, Float_t y)
{
  // set the position for the collision system

  fCollSystem->SetX(x);
  fCollSystem->SetY(y);
}

void AliCanvas::SetDataSample(TString txt)
{
  // set the description of the data sample

  fDataSample->SetText(fDataSample->GetX(), fDataSample->GetY(), txt);
}

TStyle* AliCanvas::Style()
{
  // return the figure style

  if (!fgStyle) {
    // we start from the current style
    // how can we start from a well-defined style, e.g. plain?
    fgStyle = new TStyle(*gStyle);
    fgStyle->SetName("alice");
    fgStyle->SetTitle("ALICE figure style");

    const int font = 43;
    fgStyle->SetFrameBorderMode(0);
    fgStyle->SetFrameFillColor(0);
    fgStyle->SetCanvasBorderMode(0);
    fgStyle->SetPadBorderMode(0);
    fgStyle->SetPadColor(10);
    fgStyle->SetCanvasColor(10);
    fgStyle->SetOptTitle(0);
    fgStyle->SetTitleFillColor(10);
    fgStyle->SetTitleBorderSize(0);
    if ((font % 10) == 3)
      fgStyle->SetTitleFontSize(30);
    else
      fgStyle->SetTitleFontSize(0.08);
    fgStyle->SetStatColor(10);
    fgStyle->SetStatBorderSize(1);
    // legend settings
    // fgStyle->SetLegendBorderSize(1);
    // fgStyle->SetLegendFont(43);
    // fgStyle->SetLegendFillColor(0);

    fgStyle->SetDrawBorder(0);
    fgStyle->SetTextFont(font);
    fgStyle->SetStatFont(font);
    fgStyle->SetStatFontSize(0.05);
    fgStyle->SetStatX(0.97);
    fgStyle->SetStatY(0.98);
    fgStyle->SetStatH(0.03);
    fgStyle->SetStatW(0.3);
    fgStyle->SetTickLength(0.02,"y");
    fgStyle->SetEndErrorSize(3);
    if ((font % 10) == 3)
      fgStyle->SetLabelSize(30, "xyz");
    else
      fgStyle->SetLabelSize(0.05,"xyz");
    fgStyle->SetLabelFont(font,"xyz");
    fgStyle->SetLabelOffset(0.01,"xyz");
    fgStyle->SetTitleFont(font,"xyz");
    fgStyle->SetTitleOffset(1.,"xyz");
    if ((font % 10) == 3)
      fgStyle->SetTitleSize(34, "xyz");
    else
      fgStyle->SetTitleSize(0.06,"xyz");
    fgStyle->SetMarkerSize(1);
    fgStyle->SetPalette(1,0);
    if (kFALSE) {
      fgStyle->SetOptStat(1111);
      fgStyle->SetOptFit(1111);
    }
    else {
      fgStyle->SetOptStat(0);
      fgStyle->SetOptFit(0);
    }

    fgStyle->SetLegendFont(font);
  }

  return fgStyle;
}

void AliCanvas::SetTextSize(Float_t size)
{
  // specify the text size

  fTextSize = size;
  fStatusPad->SetTextSize(fTextSize);

  this->UpdateLogoPos();
  this->UpdatePad(this);
}

void AliCanvas::UpdateLogoPos()
{
  // move the logo to the correct position

  Float_t ratio  = 1.;
  Float_t height = 0.;
  Float_t width  = 0.;

  if (fLogo && fLogo->IsValid()) {
    ratio = fLogo->GetWidth();
    ratio /= (Float_t) fLogo->GetHeight();
    ratio *= (Float_t) this->GetWh();
    ratio /= (Float_t) this->GetWw();
    height = fLogoHeight / this->GetWh();
    width  = height*ratio;
    fLogoPad->SetPad(fLogoPosX - width/2, fLogoPosY - height/2,
		     fLogoPosX + width/2, fLogoPosY + height/2);
  }

  Float_t offset = - (0.*fTextSize)/this->GetWh();

  if (fStatus == kPerformance) {
    fStatusPad->SetX(fLogoPosX);
    fStatusPad->SetY(fLogoPosY - height/2 + offset);
    offset -= (1.2*fTextSize)/this->GetWh();
  }
  else if ((fStatus == kThisWork) ||
      (fStatus == kWorkInProgress)) {
    fStatusPad->SetX(0.25);
    fStatusPad->SetY(0.95);
  }

  fDate->SetX(fLogoPosX);
  fDate->SetY(fLogoPosY - height/2 + offset);

  offset -= (1.2*fTextSize)/this->GetWh();

  fTag->SetX(fLogoPosX);
  fTag->SetY(fLogoPosY - height/2 + offset);

  this->Modified();
}

void AliCanvas::UpdatePad(TPad *pad)
{
  // enforce settings of the given pad
  // and walk through all children

  // set the geometry for the pad
  if (!fDisabledMargin) {
    pad->SetLeftMargin(0.15);
    pad->SetRightMargin(0.04);
    pad->SetBottomMargin(0.15);
  }

  fCollSystem->SetX(.96);

  // now walk through all the children
  TList *listOfPrimitives = pad->GetListOfPrimitives();
  TIter iter(listOfPrimitives);

  TObject *obj = 0x0;
  while ((obj = iter())) {
    if ((obj == fCollSystem) || (obj == fStatusPad) || (obj == fLogoPad) || (obj == fDate) || (obj == fTag) || (obj == fDataSample)) {
      listOfPrimitives->Remove(obj);
    }
    else {
      // make sure that the current style is used
      // obj->UseCurrentStyle();
      // if (obj->InheritsFrom("TGraphErrors")) {
      //   TGraphErrors *graph = (TGraphErrors*) obj;
      //   graph->SetMarkerSize(5);
      //   // graph->SetMarkerStyle(currentMarkerStyle);
      //   // graph->SetMarkerColor(currentMarkerColor);
      //   // graph->SetLineStyle(currentLineStyle);
      //   // graph->SetLineColor(currentLineColor);
      // }
      if (obj->InheritsFrom("TPad")) {
	// this->UpdatePad((TPad*) obj);
      }
      else if (obj->InheritsFrom("TH2")) {
	const TString drawOption = ((TH2*) obj)->GetDrawOption();
	if (drawOption.Contains("surf")) {
	  pad->SetLeftMargin(.2);
	  ((TH2*) obj)->SetTitleOffset(1.3,"z");
	  ((TH2*) obj)->SetTitleOffset(1.3,"y");
	  ((TH2*) obj)->GetXaxis()->CenterTitle();
	  ((TH2*) obj)->GetYaxis()->CenterTitle();
	}
	else {
	  pad->SetRightMargin(.18);
	  fCollSystem->SetX(.82);
	}
      }
      else if (obj->InheritsFrom("TLegend")) {
        TLegend *leg = (TLegend*) obj;
        leg->SetTextFont(43);
        leg->SetBorderSize(1);
        leg->SetFillStyle(0);
        leg->SetFillColor(1);
	leg->SetShadowColor(0);
        leg->SetMargin(0.25);
        leg->SetTextSize(fTextSize);
        leg->SetEntrySeparation(0.25);
      }
    }
  }
}
