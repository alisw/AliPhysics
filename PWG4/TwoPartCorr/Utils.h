#ifndef Utils_h
#define Utils_h

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TObjArray.h"

#include <string>
#include <vector>

// Or maybe this should live in TreeClasses.h?
class Noti: public TObject
{
public:
  Noti() : fc(0) {;}
  virtual ~Noti() {;}
  Bool_t Notify()         { fc=1; return 1; }
  Bool_t Notified() const { return fc;      }
  void   Reset()          { fc=0;           }
protected:
  Bool_t fc; //=1 when file changed
  ClassDef (Noti,0)
};

class PlotUtils : public TObject
 {
 public:
  static void set_hist_props(TObjArray* arr,
			     int linecolor,            
			     int fillcolor,
			     int markercolor,
			     int markerstyle,
			     float markersize = 0.75);
  
  static void set_hist_props(TH1F* h,
			     int linecolor,
			     int fillcolor,            
			     int markercolor,
			     int markerstyle,
			     float markersize = 0.75);
  
  static void set_tgraph_props(TGraphErrors* g,
			       int linecolor,            
			       int markercolor,
			       int markerstyle,
			       float markersize = 0.75);
  
  static void set_tgraph_props(TGraphAsymmErrors* g,
			       int linecolor,            
			       int markercolor,
			       int markerstyle,
			       float markersize = 0.75);
  
  static void shade(TF1* f_hi, TF1* f_lo);
  static void draw_errorbox(TH1F* mid, 
			    TH1F* hi, 
			    TH1F* lo, 
			    int color = 46, double dx = 0.06);
  static void draw_errorbox(TGraph* mid,
			    TGraph* hi,
			    TGraph* lo, 
			    int color = 46, double dx = 0.06);
  static void draw_errorbox(TGraphErrors* mid,
			    TGraphErrors* hi,
			    TGraphErrors* lo, 
			    int color = 46, double dx = 0.06);
  static void draw_errorbox(TGraphAsymmErrors* mid,
			    TGraphAsymmErrors* hi,
			    TGraphAsymmErrors* lo,
			    int color = 46, double dx = 0.06);
  
  static void draw_box_at_point(double x, double dx, double y,
				double dy_plus, double dy_minus,
				int color = 46);
  
  static double maximum(double, double, double);
  static double minimum(double, double, double);

  // Return position of mid, max, or mid argument (0, 1, or 2)
  // Return -1 if something is wrong.
  static int midpos(double x, double y, double z);
  static int maxpos(double x, double y, double z);
  static int minpos(double x, double y, double z);

  static void offset_x(TGraphErrors* g, double xoff);
  static void offset_x(TGraphAsymmErrors* g, double xoff);

  static void padsetup(TCanvas *c, int nx=4, int ny=3, std::string opt = "",
		       double lmarg = 0.12, double rmarg = 0.01,
		       double topmarg = 0.01, double botmarg = 0.12);
  
  static void multiplot(TCanvas* c, TObjArray* hArray,
			int nx, int ny,	TString opt = "", int iStart=0);
  
  static void set_ylimits(TObjArray* a1, TObjArray* a2,
			  double hispace=0, double lospace=0);

  static void set_ylimits(TH1F* h1, TH1F* h2, 
			  double hispace=0, double lospace=0);
  
  static void set_ylimits(TF1* f1, TF1* f2, 
			  double topspace = 0.1, double lowspace = 0.1);

  static void set_012(TH1F* h0, TH1F* h1, TH1F* h2);

  static void set_012(TGraphErrors* gr0, 
		      TGraphErrors* gr1, 
		      TGraphErrors* gr2);

  static void set_012(TGraphAsymmErrors* gr0, 
		      TGraphAsymmErrors* gr1, 
		      TGraphAsymmErrors* gr2);
  
  static void scale_axis_labels(TH1F* h, double c);
  static void scale_axis_labels(TObjArray* a, double c);

  static void make_nice_axes(TCanvas* can, TH1F* h, double c);
  static void make_nice_axes(TCanvas* can, TH2F* h, double c);

 private:
  // None...allows use as a static class

  ClassDef(PlotUtils,0)  
};

#endif
