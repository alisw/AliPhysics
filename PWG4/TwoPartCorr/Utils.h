// $Id$

#ifndef Utils_h
#define Utils_h

#include <string>
#include <vector>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TString.h>
#include <TObjArray.h>

class PlotUtils : public TObject
{
public:
  static void set_hist_props(TObjArray* arr,
			     Int_t linecolor,            
			     Int_t fillcolor,
			     Int_t markercolor,
			     Int_t markerstyle,
			     Double_t markersize = 0.75);
  
  static void set_hist_props(TH1F* h,
			     Int_t linecolor,
			     Int_t fillcolor,            
			     Int_t markercolor,
			     Int_t markerstyle,
			     Double_t markersize = 0.75);
  
  static void set_tgraph_props(TGraphErrors* g,
			       Int_t linecolor,            
			       Int_t markercolor,
			       Int_t markerstyle,
			       Double_t markersize = 0.75);
  
  static void set_tgraph_props(TGraphAsymmErrors* g,
			       Int_t linecolor,            
			       Int_t markercolor,
			       Int_t markerstyle,
			       Double_t markersize = 0.75);
  
  static void shade(TF1* f_hi, TF1* f_lo);
  static void draw_errorbox(TH1F* mid, 
			    TH1F* hi, 
			    TH1F* lo, 
			    Int_t color = 46, Double_t dx = 0.06);
  static void draw_errorbox(TGraph* mid,
			    TGraph* hi,
			    TGraph* lo, 
			    Int_t color = 46, Double_t dx = 0.06);
  static void draw_errorbox(TGraphErrors* mid,
			    TGraphErrors* hi,
			    TGraphErrors* lo, 
			    Int_t color = 46, Double_t dx = 0.06);
  static void draw_errorbox(TGraphAsymmErrors* mid,
			    TGraphAsymmErrors* hi,
			    TGraphAsymmErrors* lo,
			    Int_t color = 46, Double_t dx = 0.06);
  
  static void draw_box_at_point(Double_t x, Double_t dx, Double_t y,
				Double_t dy_plus, Double_t dy_minus,
				Int_t color = 46);
  
  static Double_t maximum(Double_t, Double_t, Double_t);
  static Double_t minimum(Double_t, Double_t, Double_t);

  // Return position of mid, max, or mid argument (0, 1, or 2)
  // Return -1 if something is wrong.
  static Int_t midpos(Double_t x, Double_t y, Double_t z);
  static Int_t maxpos(Double_t x, Double_t y, Double_t z);
  static Int_t minpos(Double_t x, Double_t y, Double_t z);

  static void offset_x(TGraphErrors* g, Double_t xoff);
  static void offset_x(TGraphAsymmErrors* g, Double_t xoff);

  static void padsetup(TCanvas *c, Int_t nx=4, Int_t ny=3, std::string opt = "",
		       Double_t lmarg = 0.12, Double_t rmarg = 0.01,
		       Double_t topmarg = 0.01, Double_t botmarg = 0.12);
  
  static void multiplot(TCanvas* c, TObjArray* hArray,
			Int_t nx, Int_t ny, TString opt = "", Int_t iStart=0);
  
  static void set_ylimits(TObjArray* a1, TObjArray* a2,
			  Double_t hispace=0, Double_t lospace=0);

  static void set_ylimits(TH1F* h1, TH1F* h2, 
			  Double_t hispace=0, Double_t lospace=0);
  
  static void set_ylimits(TF1* f1, TF1* f2, 
			  Double_t topspace = 0.1, Double_t lowspace = 0.1);

  static void set_012(TH1F* h0, TH1F* h1, TH1F* h2);

  static void set_012(TGraphErrors* gr0, 
		      TGraphErrors* gr1, 
		      TGraphErrors* gr2);

  static void set_012(TGraphAsymmErrors* gr0, 
		      TGraphAsymmErrors* gr1, 
		      TGraphAsymmErrors* gr2);
  
  static void scale_axis_labels(TH1F* h, Double_t c);
  static void scale_axis_labels(TObjArray* a, Double_t c);

  static void make_nice_axes(TCanvas* can, TH1F* h, Double_t c);
  static void make_nice_axes(TCanvas* can, TH2F* h, Double_t c);

 private:
  // None...allows use as a static class

  ClassDef(PlotUtils,0) // Plot utilities class
};
#endif
