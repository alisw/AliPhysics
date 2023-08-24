#ifndef __PLOTTING_HH__
#define __PLOTTING_HH__

//-------------------------------------------------------
// Plotting beutificaton
//-------------------------------------------------------

#include "CandRoot.h"


void BeutifyCanvas(TCanvas* cv, Double_t top,Double_t right,Double_t bottom,Double_t left,Bool_t setLogY = kFALSE, Bool_t setLogZ =kFALSE);
void BeutifyPad(TVirtualPad* pad,Double_t top,Double_t right,Double_t bottom,Double_t left, Bool_t setTicks, Bool_t divide = kFALSE);
void BeutifyLegend(TLegend* legend, Double_t fontSize, Int_t font=42,Int_t fillColor=0,Double_t borderSize=0);
void BeutifyLatex(TLatex* lat, Double_t size =0.04,Int_t font =42, Int_t color =1);
void BeutifyTH1(TH1* h, const Char_t* title, Double_t lineWidth, Int_t lineColor,
                        Int_t markerStyle=1, Int_t markerColor=1, Double_t markerSize=1);
void BeutifyTH2(TH2* h, Bool_t stats,const Char_t* title);
void BeutifyTAxis(TAxis* ax, Double_t rangeMin, Double_t rangeMax, Double_t titleOffset=1., Double_t titleSize=0.04, Int_t titleFont=42,
                                               Double_t labelSize=0.04, Int_t labelFont=42, Int_t nDivisions=510);

//______________________________________________________________________________________
void BeutifyCanvas(TCanvas* cv, Double_t top,Double_t right,Double_t bottom,Double_t left,Bool_t setLogX /*= kFALSE*/, Bool_t setLogY /*= kFALSE*/){
  cv->SetTopMargin(top);
  cv->SetRightMargin(right);
  cv->SetBottomMargin(bottom);
  cv->SetLeftMargin(left);
  cv->SetTickx();
  cv->SetTicky();
  if (setLogX) cv->SetLogx();
  if (setLogY) cv->SetLogy();
}

//_______________________________________________________________________________________
void BeutifyLegend(TLegend* legend, Double_t fontSize, Int_t font /*=42*/,Int_t fillColor /*=0*/, Double_t borderSize /*=0*/){
  legend->SetTextSize(fontSize);
  legend->SetTextFont(font);
  legend->SetFillColor(fillColor);
  legend->SetBorderSize(borderSize);
}

void BeutifyLatex(TLatex* lat, Double_t size /*=0.04*/,Int_t font /*=42*/, Int_t color /*=1*/){
  lat->SetNDC();
  lat->SetTextSize(size);
  lat->SetTextFont(font);
  lat->SetTextColor(color);
}


//________________________________________________________________________________________
 void BeutifyTH1(TH1* h, const Char_t* title, Double_t lineWidth, Int_t lineColor,
                                              Int_t markerStyle /*=1*/, Int_t markerColor /*=1*/, Double_t markerSize /*=1*/) {
 //
 // set drawing options for a TH1
 //
    if(!h) return;
    h->SetTitle(title);
    h->SetLineWidth(lineWidth);
    h->SetLineColor(lineColor);
    h->SetMarkerStyle(markerStyle);
    h->SetMarkerColor(markerColor);
    h->SetMarkerSize(markerSize);
}
//________________________________________________________________________________________
 void BeutifyTH2(TH2* h, Bool_t stats,const Char_t* title) {
 //
 // set drawing options for a TH2
 //
    if(!h) return;
    h->SetStats(stats);
    h->SetTitle(title);
}

//________________________________________________________________________________________
void BeutifyTH1(TH1* h, const Char_t* title, Double_t lineWidth, Int_t lineColor,
                                              Int_t markerStyle /*=1*/, Int_t markerColor /*=1*/, Double_t markerSize /*=1*/,Bool_t stats /*=kFALSE*/) {
 //
 // set drawing options for a TH1
 //

    if(!h) return;

    h->SetTitle(title);
    h->SetLineWidth(lineWidth);
    h->SetLineColor(lineColor);
    h->SetMarkerStyle(markerStyle);
    h->SetMarkerColor(markerColor);
    h->SetMarkerSize(markerSize);
    h->SetStats(stats);
}

//________________________________________________________________________________________
void BeutifyTAxis(TAxis* ax, Double_t rangeMin, Double_t rangeMax, Double_t titleOffset /*=1.*/, Double_t titleSize /*=0.04*/, Int_t titleFont /*=42*/,
                                               Double_t labelSize /*=0.04*/, Int_t labelFont /*=42*/, Int_t nDivisions /*=510*/) {
  //
  // set drawing options for TAxis
  //
  if(!ax) return;
  ax->SetRangeUser(rangeMin, rangeMax);
  ax->SetTitleSize(titleSize);
  ax->SetTitleOffset(titleOffset);
  ax->SetTitleFont(titleFont);
  ax->SetLabelSize(labelSize);
  ax->SetLabelFont(labelFont);
  ax->SetNdivisions(nDivisions);
}

#endif

