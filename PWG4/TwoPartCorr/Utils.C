// $Id$

#include "Utils.h"
#include "TMath.h"
#include "TLine.h"
#include "TList.h"

using std::cout;
using std::endl;

ClassImp(PlotUtils)

void PlotUtils::set_hist_props(TObjArray* arr,
			       Int_t linecolor,
			       Int_t fillcolor,
			       Int_t markercolor,
			       Int_t markerstyle,
			       Double_t markersize) 
{
  TH1F* h = 0;
  
  for (Int_t i=0; i<arr->GetEntries(); i++) {
    
    if (arr->At(i)->InheritsFrom("TH1F"))  {
      h = (TH1F*)arr->At(i);
    }
    
    h->SetLineColor(linecolor);
    h->SetFillColor(fillcolor);
    h->SetMarkerColor(markercolor);
    h->SetMarkerStyle(markerstyle);
    h->SetMarkerSize(markersize);
  }
}

void PlotUtils::set_hist_props(TH1F* h,
			       Int_t linecolor,
			       Int_t fillcolor,
			       Int_t markercolor,
			       Int_t markerstyle,
			       Double_t markersize) 
{
  h->SetLineColor(linecolor);
  h->SetFillColor(fillcolor);
  h->SetMarkerColor(markercolor);
  h->SetMarkerStyle(markerstyle);
  h->SetMarkerSize(markersize);
}

void PlotUtils::set_tgraph_props(TGraphErrors* h,
				 Int_t linecolor,
				 Int_t markercolor,
				 Int_t markerstyle,
				 Double_t markersize) 
{
  h->SetLineColor(linecolor);
  h->SetMarkerColor(markercolor);
  h->SetMarkerStyle(markerstyle);
  h->SetMarkerSize(markersize);
  h->SetLineWidth(2);
}

void PlotUtils::set_tgraph_props(TGraphAsymmErrors* h,
				 Int_t linecolor,
				 Int_t markercolor,
				 Int_t markerstyle,
				 Double_t markersize) 
{
  h->SetLineColor(linecolor);
  h->SetMarkerColor(markercolor);
  h->SetMarkerStyle(markerstyle);
  h->SetMarkerSize(markersize);
  h->SetLineWidth(2);
}

void PlotUtils::shade(TF1 * fup, TF1 * flow)
{
  TGraph * gr = new TGraph();
  gr->SetFillColor(fup->GetFillColor());
  gr->SetFillStyle(fup->GetFillStyle());
  Double_t xmin = fup->GetXmin();
  Double_t xmax = fup->GetXmax();

  //process top side function
  Int_t npx = flow->GetNpx();
  Int_t npoints=0;
  Double_t dx = (xmax-xmin)/npx;
  Double_t x = xmin+0.5*dx;
  Double_t y =0;
  while (x <= xmax) {
    Double_t yup = fup->Eval(x);
    Double_t ylow = flow->Eval(x);

    if ( yup > ylow )
      y=yup;
    else
      y=ylow;
    gr->SetPoint(npoints,x,y);
    npoints++;
    x += dx;
  }
  //process bottom side                                                                                         
  x = xmax-0.5*dx;
  while (x >= xmin) {
    Double_t ylow = flow->Eval(x);
    Double_t yup = fup->Eval(x);
    if ( yup < ylow )
      y=yup;
    else
      y=ylow;
    gr->SetPoint(npoints,x,y);
    npoints++;
    x -= dx;
  }
  gr->Draw("f");
}

void PlotUtils::draw_errorbox(TH1F* mid, TH1F* hi, TH1F* lo, Int_t color, Double_t dx) 
{
  Double_t x, y, dy_hi, dy_lo;

  for (Int_t j=1; j <= mid->GetNbinsX(); j++) {

    if ( ! mid->GetBinContent(j) ) continue;

    x  = mid->GetBinCenter(j);
    dx = 0.08;
    y  = mid->GetBinContent(j);
    dy_hi = fabs(hi->GetBinContent(j) - mid->GetBinContent(j) );
    dy_lo =   fabs(lo->GetBinContent(j) - mid->GetBinContent(j));

    draw_box_at_point(x, dx, y, dy_hi, dy_lo, color);
  }
}

void PlotUtils::draw_errorbox(TGraph* mid,
			      TGraph* hi,
			      TGraph* lo, Int_t color, Double_t dx) 
{
  Double_t* x = mid->GetX();
  Double_t* y =  mid->GetY();
  Double_t* y_hi= hi->GetY();
  Double_t* y_lo= lo->GetY();

  for (Int_t j=0; j < mid->GetN(); j++) {

    Double_t dy_hi = fabs(y_hi[j] - y[j]);
    Double_t dy_lo = fabs(y_lo[j] - y[j]);

    draw_box_at_point(x[j], dx, y[j], dy_hi, dy_lo, color);
  }
  return;
}

void PlotUtils::draw_errorbox(TGraphErrors* mid,
			      TGraphErrors* hi,
			      TGraphErrors* lo, Int_t color, Double_t dx) 
{
  Double_t* x = mid->GetX();
  Double_t* y =  mid->GetY();
  Double_t* y_hi= hi->GetY();
  Double_t* y_lo= lo->GetY();

  for (Int_t j=0; j < mid->GetN(); j++) {

    Double_t dy_hi = fabs(y_hi[j] - y[j]);
    Double_t dy_lo = fabs(y_lo[j] - y[j]);

    draw_box_at_point(x[j], dx, y[j], dy_hi, dy_lo, color);
  }
  return;
}

void PlotUtils::draw_errorbox(TGraphAsymmErrors* mid,
			      TGraphAsymmErrors* hi,
			      TGraphAsymmErrors* lo, Int_t color, Double_t dx) 
{
  Double_t* x = mid->GetX();
  Double_t* y =  mid->GetY();
  Double_t* y_hi= hi->GetY();
  Double_t* y_lo= lo->GetY();

  for (Int_t j=0; j < mid->GetN(); j++) {

    Double_t dy_hi = fabs(y_hi[j] - y[j]);
    Double_t dy_lo = fabs(y_lo[j] - y[j]);

    draw_box_at_point(x[j], dx, y[j], dy_hi, dy_lo, color);
  }
}


void PlotUtils::draw_box_at_point(Double_t x, Double_t dx, Double_t y,
                                 Double_t dy_plus, Double_t dy_minus,
                                 Int_t color) 
{
  TLine tl;
  tl.SetLineWidth(2);
  tl.SetLineColor(color);

  tl.DrawLine(x-dx, y-dy_minus, x+dx, y-dy_minus);
  tl.DrawLine(x-dx, y+dy_plus, x+dx, y+dy_plus);
  tl.DrawLine(x-dx, y+dy_plus, x-dx, y-dy_minus);
  tl.DrawLine(x+dx, y+dy_plus, x+dx, y-dy_minus);
}

// This function modifies the graphs so that 0, 1, 2 truly have the
// mid, high, and low points at each x value.
void PlotUtils::set_012(TH1F* h0, 
			TH1F* h1, 
			TH1F* h2) 
{
  for (Int_t j=0; j<h0->GetNbinsX(); j++) {
    Double_t  y[3], dy[3];
    Int_t kmid=0, kmax=1, kmin=2;

    TH1F* h[3]; 
    h[0] = h0; // k=0: set to mid 
    h[1] = h1; // k=1: set to high
    h[2] = h2; // k=2: set to low
    
    for (Int_t k=0; k<3; k++) {
      y[k] =  h[k]->GetBinContent(j);
      dy[k] = h[k]->GetBinError(j);
    }

    kmid = midpos(y[0], y[1], y[2]);      
    kmax = maxpos(y[0], y[1], y[2]);
    kmin = minpos(y[0], y[1], y[2]);
    
    h0->SetBinContent(j, y[kmid]);
    h1->SetBinContent(j, y[kmax]);
    h2->SetBinContent(j, y[kmin]);

    h0->SetBinError(j, dy[kmid]);
    h1->SetBinError(j, dy[kmax]);
    h2->SetBinError(j, dy[kmin]);
  } // j loop
}

// This function modifies the graphs so that 0, 1, 2 truly have the
// mid, high, and low points at each x value.
void PlotUtils::set_012(TGraphErrors* gr0, 
			TGraphErrors* gr1, 
			TGraphErrors* gr2) 
{
  // Assuming that gr0,1,2 have matching
  // x-values at each poInt_t j.
  for (Int_t j=0; j<gr0->GetN(); j++) {
    
    Double_t x, y[3], dy[3];
    Int_t kmid=0, kmax=1, kmin=2;

    TGraphErrors* g[3]; 
    g[0] = gr0; // k=0: set to mid 
    g[1] = gr1; // k=1: set to high
    g[2] = gr2; // k=2: set to low
    
    for (Int_t k=0; k<3; k++) {
      Double_t ytmp=0;
      g[k]->GetPoint(j, x, ytmp);
      y[k] = ytmp;
      dy[k] = g[k]->GetErrorY(j);
    }

    kmid = midpos(y[0], y[1], y[2]);      
    kmax = maxpos(y[0], y[1], y[2]);
    kmin = minpos(y[0], y[1], y[2]);
        
    gr0->SetPoint(j, x, y[kmid]);
    gr1->SetPoint(j, x, y[kmax]);
    gr2->SetPoint(j, x, y[kmin]);

    gr0->SetPointError(j, 0, dy[kmid]);
    gr1->SetPointError(j, 0, dy[kmax]);
    gr2->SetPointError(j, 0, dy[kmin]);
  } // j loop
}

// This function modifies the graphs so that 0, 1, 2 truly have the
// mid, high, and low points at each x value.
void PlotUtils::set_012(TGraphAsymmErrors* gr0, 
			TGraphAsymmErrors* gr1, 
			TGraphAsymmErrors* gr2) 
{
  // Assuming that gr0,1,2 have matching
  // x-values at each poInt_t j.
  for (Int_t j=0; j<gr0->GetN(); j++) {
    
    Double_t x, y[3], yhi[3], ylo[3];
    Int_t kmid=0, kmax=1, kmin=2;

    TGraphAsymmErrors* g[3]; 
    g[0] = gr0; // k=0: set to mid 
    g[1] = gr1; // k=1: set to high
    g[2] = gr2; // k=2: set to low
    
    for (Int_t k=0; k<3; k++) {
      Double_t ytmp=0;
      g[k]->GetPoint(j, x, ytmp);
      y[k] = ytmp;
      yhi[k] = g[k]->GetErrorYhigh(j);
      ylo[k] = g[k]->GetErrorYlow(j);
    }
    
    kmid = midpos(y[0], y[1], y[2]);      
    kmax = maxpos(y[0], y[1], y[2]);
    kmin = minpos(y[0], y[1], y[2]);
    
    gr0->SetPoint(j, x, y[kmid]);
    gr1->SetPoint(j, x, y[kmax]);
    gr2->SetPoint(j, x, y[kmin]);
    
    gr0->SetPointError(j, 0, 0, ylo[kmid], yhi[kmid]);
    gr1->SetPointError(j, 0, 0, ylo[kmax], yhi[kmax]);
    gr2->SetPointError(j, 0, 0, ylo[kmin], yhi[kmin]);
  } // j loop
}

Int_t PlotUtils::maxpos(Double_t x, Double_t y, Double_t z) 
{
  Double_t max = maximum(x,y,z);
  if(x==max) return 0;
  if(y==max) return 1;
  if(z==max) return 2;

  return -1;
}

Int_t PlotUtils::minpos(Double_t x, Double_t y, Double_t z) 
{
  Double_t min = minimum(x,y,z);
  if(x==min) return 0;
  if(y==min) return 1;
  if(z==min) return 2;

  return -1;
}

Int_t PlotUtils::midpos(Double_t x, Double_t y, Double_t z) 
{
  Int_t maxp = maxpos(x,y,z);
  Int_t minp = minpos(x,y,z);
  
  if (maxp==-1 || minp==-1){
    cout << "PlotUtils::midpos(): Bad parameters!" 
	 << endl;
    return -1;
  }

  if (maxp==minp) return 0;
  if ((maxp==1 && minp==2) || (maxp==2 && minp==1)) return 0;
  if ((maxp==0 && minp==2) || (maxp==2 && minp==0)) return 1;
  if ((maxp==0 && minp==1) || (maxp==1 && minp==0)) return 2;
  
  cout << "PlotUtils::midpos(): Did not return a middle position!" 
       << endl;
  return -1;
}

Double_t PlotUtils::maximum(Double_t x, Double_t y, Double_t z) 
{
  Double_t max = (x > y ? x : y);
  max = (max > z ? max : z);
  return max;
}

Double_t PlotUtils::minimum(Double_t x, Double_t y, Double_t z) 
{
  Double_t min = (x < y ? x : y);
  min = (min < z ? min : z);
  return min;
}

void PlotUtils::offset_x(TGraphErrors* g, Double_t xoff)
{
  Int_t npoints = g->GetN();
  Double_t* x = g->GetX();
  Double_t* y = g->GetY();
  for (Int_t j=0; j<npoints; j++){
    g->SetPoint(j, x[j] + xoff, y[j]);
  }
}

void PlotUtils::offset_x(TGraphAsymmErrors* g, Double_t xoff)
{
  Int_t npoints = g->GetN();
  Double_t* x = g->GetX();
  Double_t* y = g->GetY();
  for (Int_t j=0; j<npoints; j++){
    g->SetPoint(j, x[j] + xoff, y[j]);
  }
}

void PlotUtils::padsetup(TCanvas* c, Int_t nx, Int_t ny, std::string opt,
			 Double_t lmarg, Double_t rmarg, 
			 Double_t himarg, Double_t lomarg)
{
  c->Divide(nx,ny,0,0);

  Int_t pad[100][100];

  for (Int_t x=1; x<=nx; x++) {
    for (Int_t y=1; y<=ny; y++) {
      pad[x][y] = nx * (y-1) + x;
    }
  }

  if (opt=="mergex" || opt=="MERGEX") {
    for (Int_t x=1; x<=nx; x++) {
      for (Int_t y=1; y<=ny; y++) {
	// left edge...
	if (x == 1)
	  c->GetPad(pad[x][y])->SetLeftMargin(lmarg);
	else
	  c->GetPad(pad[x][y])->SetLeftMargin(0);

	// right edge...
	if (x==nx)
	  c->GetPad(pad[x][y])->SetRightMargin(rmarg);
	else
	  c->GetPad(pad[x][y])->SetRightMargin(0);
	// top edge...
	if (y==1)  c->GetPad(pad[x][y])->SetTopMargin(himarg);
	// bottom...
	if (y==ny) c->GetPad(pad[x][y])->SetBottomMargin(lomarg);
      }
    }
  }
  else {
    for (Int_t x=1; x<=nx; x++) {
      for (Int_t y=1; y<=ny; y++) {
	// left edge...
	c->GetPad(pad[x][y])->SetLeftMargin(lmarg);
	// right edge...
	c->GetPad(pad[x][y])->SetRightMargin(rmarg);
	// top edge...
	if (y==1)  c->GetPad(pad[x][y])->SetTopMargin(himarg);
	// bottom...
	if (y==ny) c->GetPad(pad[x][y])->SetBottomMargin(lomarg);
      }
    }
  }
  for (Int_t i=1; i<=nx*ny; i++){
    c->GetPad(i)->SetTickx();
    c->GetPad(i)->SetTicky();
    //    c->GetPad(i)->SetGridy();
    c->GetPad(i)->SetFrameLineWidth(2);
  }
}

void PlotUtils::set_ylimits(TObjArray* arr1, TObjArray* arr2,
			    Double_t topspace, Double_t lowspace) 
{
  Int_t n1 = arr1->GetEntries();
  Int_t n2 = arr2->GetEntries();
  if (n1 != n2) {
    cout << __FILE__ << " " << __LINE__ <<": " << endl;
    cout << "Error: Arrays have unequal length: " 
	 << n1 << " vs " << n2 << "." << endl;
    return;
  }

  for (Int_t i=0; i<n1; i++) {
    TH1F* h1 = (TH1F*)arr1->At(i);
    TH1F* h2 = (TH1F*)arr2->At(i);
    set_ylimits(h1, h2, topspace, lowspace);
  }
}


void PlotUtils::set_ylimits(TH1F* h1, TH1F* h2, 
			    Double_t topspace, Double_t lowspace) 
{
  Int_t maxbin1 = h1->GetMaximumBin(); 
  Int_t maxbin2 = h2->GetMaximumBin(); 
  Int_t minbin1 = h1->GetMinimumBin(); 
  Int_t minbin2 = h2->GetMinimumBin(); 
  
  Double_t max1 = h1->GetBinContent(maxbin1) + h1->GetBinError(maxbin1);
  Double_t max2 = h2->GetBinContent(maxbin2) + h2->GetBinError(maxbin2);
  
  Double_t min1 = h1->GetBinContent(minbin1) - h1->GetBinError(minbin1);
  Double_t min2 = h2->GetBinContent(minbin2) - h2->GetBinError(minbin2);
   
  Double_t max = max1 > max2 ? max1 : max2;
  Double_t min = min1 < min2 ? min1 : min2;

  if (topspace == 0) topspace = 0.1*max;
  if (lowspace == 0) lowspace = 0.1*max;

  Double_t upperSpace = (max - min)*topspace;
  Double_t lowerSpace = (max - min)*lowspace;

  h1->GetYaxis()->SetRangeUser(min - lowerSpace, max + upperSpace);
  h2->GetYaxis()->SetRangeUser(min - lowerSpace, max + upperSpace);
}

void PlotUtils::set_ylimits(TF1* f1, TF1* f2, 
			    Double_t topspace, Double_t lowspace)
{
  Double_t xlo1 = f1->GetXaxis()->GetXmin();
  Double_t xhi1 = f1->GetXaxis()->GetXmax();
  Double_t xlo2 = f2->GetXaxis()->GetXmin();
  Double_t xhi2 = f2->GetXaxis()->GetXmax();

  Double_t min1 = f1->GetMinimum(xlo1, xhi1);
  Double_t max1 = f1->GetMaximum(xlo1, xhi1);

  Double_t min2 = f2->GetMinimum(xlo2, xhi2);
  Double_t max2 = f2->GetMaximum(xlo2, xhi2);

  Double_t max = max1 > max2 ? max1 : max2;
  Double_t min = min1 < min2 ? min1 : min2;

  Double_t upperSpace = (max - min)*topspace;
  Double_t lowerSpace = (max - min)*lowspace;

  f1->GetYaxis()->SetRangeUser(min - lowerSpace, max + upperSpace);
  f2->GetYaxis()->SetRangeUser(min - lowerSpace, max + upperSpace);
}

void PlotUtils::multiplot(TCanvas* c, TObjArray *hArray,
                          Int_t nx, Int_t ny, TString opt, Int_t iStart)
{
  Int_t ipad = 1;
  
  c->cd();
  TList* prims = c->GetListOfPrimitives();
  if (prims->GetEntries() > 0) {
    opt.Append("same");
  }
  else padsetup(c, nx, ny);
  
  for (Int_t i = iStart; i < iStart+nx*ny; i++) {
    c->cd(ipad);
    
    TObject* obj = hArray->At(i);

    if (!obj) 
      cout << __FILE__ << " " << __LINE__ 
	   << ": Object in hArray is null." << endl;
    else if (obj->InheritsFrom("TH1F")) 
      ((TH1F*)hArray->At(i))->Draw(opt.Data());
    else if (obj->InheritsFrom("TH1D")) 
      ((TH1D*)hArray->At(i))->Draw(opt.Data());
    else if (obj->InheritsFrom("TF1")) 
      ((TF1*)hArray->At(i))->Draw(opt.Data());
    else
      (hArray->At(i))->Draw(opt.Data());
    
    ipad++;
  }
}

void PlotUtils::scale_axis_labels(TObjArray* a, Double_t c)
{
  for (Int_t i=0; i<a->GetEntries(); i++) {
    TObject* obj = a->At(i);
    if (obj->InheritsFrom("TH1F")) {
      TH1F* h = (TH1F*)obj;
      scale_axis_labels(h, c);
    }
  }
}

void PlotUtils::scale_axis_labels(TH1F* h, Double_t c)
{
  TAxis* x = h->GetXaxis();
  TAxis* y = h->GetYaxis();

  x->SetLabelSize(c*x->GetLabelSize());
  y->SetLabelSize(c*y->GetLabelSize());

  x->SetLabelOffset(c*x->GetLabelOffset());
  y->SetLabelOffset(c*y->GetLabelOffset());

  x->SetTitleSize(c*x->GetTitleSize());
  y->SetTitleSize(c*y->GetTitleSize());
  return;
}

void PlotUtils::make_nice_axes(TCanvas* can, TH1F* h, Double_t c)
{
  TAxis* x = h->GetXaxis();
  TAxis* y = h->GetYaxis();

  x->SetLabelSize(c*x->GetLabelSize());
  y->SetLabelSize(c*y->GetLabelSize());

  x->SetLabelOffset(c*x->GetLabelOffset());
  y->SetLabelOffset(c*y->GetLabelOffset());

  x->SetTitleSize(c*x->GetTitleSize());
  y->SetTitleSize(c*y->GetTitleSize());

  // new below
  x->SetTitleOffset(0.9*c*x->GetTitleOffset());
  y->SetTitleOffset(c*y->GetTitleOffset());

  // default left and bottom margins: 0.10
  // assume defaults, and add extra space
  can->SetLeftMargin(c*0.12);
  can->SetBottomMargin(c*0.12);

  x->CenterTitle();
  y->CenterTitle();
}

// I should really learn how to write templates...
void PlotUtils::make_nice_axes(TCanvas* can, TH2F* h, Double_t c)
{
  TAxis* x = h->GetXaxis();
  TAxis* y = h->GetYaxis();

  x->SetLabelSize(c*x->GetLabelSize());
  y->SetLabelSize(c*y->GetLabelSize());

  x->SetLabelOffset(c*x->GetLabelOffset());
  y->SetLabelOffset(c*y->GetLabelOffset());

  x->SetTitleSize(c*x->GetTitleSize());
  y->SetTitleSize(c*y->GetTitleSize());

  // new below
  x->SetTitleOffset(c*x->GetTitleOffset());
  y->SetTitleOffset(c*y->GetTitleOffset());

  // default left and bottom margins: 0.10
  // assume defaults, and add extra space
  can->SetLeftMargin(c*0.12);
  can->SetBottomMargin(c*0.12);

  x->CenterTitle();
  y->CenterTitle();
}
