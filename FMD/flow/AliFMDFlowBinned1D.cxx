/* Copyright (C) 2007 Christian Holm Christensen <cholm@nbi.dk>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA
 */
/** @file 
    @brief Implementation of a 1-dimensional Flow "histogram" */
//____________________________________________________________________ 
//
// A histogram of flow bins.  The axis can by anything
// (pseudo-rapidity, transvers momentum) - there's no assumption on
// what is the basis of the histogram.  The method Event can be used
// to calculate everything in one go.   Alternatively, one can use the
// methods AddToEventPlane and AddToHarmonic.  See also the example
// TestFlow.C 
#include "flow/AliFMDFlowBinned1D.h"
#include "flow/AliFMDFlowBin.h"
#include "flow/AliFMDFlowSplitter.h"
// #include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <TBrowser.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>

//====================================================================
AliFMDFlowBinned1D::AliFMDFlowBinned1D()
  : TNamed("", ""), 
    fXAxis(),
    fN(0),
    fBins(0),
    fSplitter(0)
{
  // Default CTOR - do not use
}

//____________________________________________________________________
AliFMDFlowBinned1D::AliFMDFlowBinned1D(const char* name, 
				       const char* title, 
				       UShort_t order, 
				       UShort_t k,
				       UShort_t nxbins, 
				       Double_t* xbins, 
				       AliFMDFlowSplitter* splitter) 
  : TNamed(name,title), 
    fXAxis(nxbins, xbins),
    fN(0),
    fBins(0), 
    fSplitter(splitter)
{
  // Constructor 
  // Parameters: 
  //   Order	Order 
  //   nxbins   Number of bins 
  //   xbins    Bin borders 
  UShort_t n = fXAxis.N();
  fN         = n;
  fBins      = new AliFMDFlowBin*[n];
  for (UInt_t i = 0; i < n; i++) fBins[i]= new AliFMDFlowBin(order,k);
  if (!fSplitter) fSplitter = new AliFMDFlowShuffle;
}
//____________________________________________________________________
AliFMDFlowBinned1D::AliFMDFlowBinned1D(const char*           name, 
				       const char*           title, 
				       UShort_t              order, 
				       UShort_t              k,
				       const AliFMDFlowAxis& xaxis,
				       AliFMDFlowSplitter*   splitter)
  : TNamed(name,title), 
    fXAxis(xaxis), 
    fN(0),
    fBins(0), 
    fSplitter(splitter)
{
  // Constructor 
  // Parameters: 
  //   Order	Order 
  //   xaxis    X axis object
  UShort_t n = fXAxis.N();
  fN         = n;
  fBins      = new AliFMDFlowBin*[n];
  for (UInt_t i = 0; i < n; i++) fBins[i]= new AliFMDFlowBin(order,k);
  if (!fSplitter) fSplitter = new AliFMDFlowShuffle;
    
}
//____________________________________________________________________
AliFMDFlowBinned1D::AliFMDFlowBinned1D(const AliFMDFlowBinned1D& o)
  : TNamed(o), 
    TAttLine(o),
    TAttFill(o),
    TAttMarker(o),
    fXAxis(o.fXAxis), 
    fN(0),
    fBins(0),
    fSplitter(0)
{
  // Copy constructor 
  // Parameters: 
  //   o   Object to copy from 
  UShort_t n = fXAxis.N();
  fSplitter  = new AliFMDFlowShuffle;
  fN         = n;
  fBins      = new AliFMDFlowBin*[n];
  for (UInt_t i = 0; i < n; i++) fBins[i]= new AliFMDFlowBin(*(o.fBins[i]));
}
//____________________________________________________________________
AliFMDFlowBinned1D&
AliFMDFlowBinned1D::operator=(const AliFMDFlowBinned1D& o)
{
  // Assignment operator
  // Parameters: 
  //   o Object to assign from 
  // 
  // Returns reference to this object 
  SetLineColor(o.GetLineColor());
  SetLineStyle(o.GetLineStyle());
  SetLineWidth(o.GetLineWidth());
  SetFillColor(o.GetFillColor());
  SetFillStyle(o.GetFillStyle());
  SetMarkerColor(o.GetMarkerColor());
  SetMarkerStyle(o.GetMarkerStyle());
  SetMarkerSize(o.GetMarkerSize());
  this->SetName(o.GetName());
  this->SetTitle(o.GetTitle());
  if (fBins) { 
    for (UInt_t i = 0; i < fXAxis.N(); i++) delete fBins[i];
    delete [] fBins;
  }
  if (fSplitter) delete fSplitter;
  fXAxis     = o.fXAxis;
  UShort_t n = fXAxis.N();
  fN         = n;
  fSplitter  = new AliFMDFlowShuffle;
  fBins      = new AliFMDFlowBin*[n];
  for (UInt_t i = 0; i < n; i++) fBins[i]= new AliFMDFlowBin(*(o.fBins[i]));
  return *this;
}
  
//____________________________________________________________________
AliFMDFlowBinned1D::~AliFMDFlowBinned1D()
{
  // Destructor 
  // Parameters: 
  //    none 
  if (fBins) { 
    for (UInt_t i = 0; i < fXAxis.N(); i++) delete fBins[i];
    delete [] fBins;
  }
  if (fSplitter) delete fSplitter;
}

//____________________________________________________________________
AliFMDFlowBin* 
AliFMDFlowBinned1D::GetBin(UShort_t i) const
{
  // Get the ith bin 
  // Parameters: 
  //    i  Bin number 
  // 
  // Return pointer to bin, or null. 
  if (i >= fXAxis.N()) return 0;
  return fBins[i];
}
//____________________________________________________________________
AliFMDFlowBin* 
AliFMDFlowBinned1D::GetBin(Double_t x) const
{
  // Get the bin that contains x
  // Parameters: 
  //    x  X axis value
  // 
  // Return pointer to bin, or null. 
  Int_t i = fXAxis.FindBin(x);
  if (i < 0) return 0;
  UShort_t j = i;
  return GetBin(j);
}
  
//____________________________________________________________________
void 
AliFMDFlowBinned1D::Begin()
{
  // Called at the beginning of an event
  // Parameters: 
  //   none
  for (UInt_t i = 0; i < fXAxis.N(); i++) fBins[i]->Begin();
  fSplitter->Begin();
}
//____________________________________________________________________
void 
AliFMDFlowBinned1D::End()
{
  // Called at the end of an event
  // Parameters: 
  //   none
  for (UInt_t i = 0; i < fXAxis.N(); i++) fBins[i]->End();
  fSplitter->End();
}
//____________________________________________________________________
void 
AliFMDFlowBinned1D::Finish()
{
  // Called at the end of an job
  // Parameters: 
  //   none
  for (UInt_t i = 0; i < fXAxis.N(); i++) fBins[i]->Finish();
}
//____________________________________________________________________
Bool_t 
AliFMDFlowBinned1D::AddToEventPlane(Double_t x, Double_t phi, 
				    Double_t w, Bool_t a)
{
  // Called to add a contribution to the event plane 
  // Parameters:
  //    x   Bin to fill into 
  //    w   Weight
  //    phi The angle phi in radians 
  //    a   If true, add to sub-event A, otherwise sub-event B
  // 
  // Return false if x falls outside the defined range, true otherwise
  AliFMDFlowBin* bin = GetBin(x);
  if (!bin) return kFALSE;
  bin->AddToEventPlane(phi, w, a);
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDFlowBinned1D::AddToHarmonic(Double_t x,  Double_t phi, 
				  Double_t wp, Double_t wh)
{
  // Called to add a contribution to the harmonic
  // Parameters: 
  //    x    Bin to fill into 
  //    phi  The angle phi in radians
  //    wp   weight of event plane
  //    wh   weight of harmonic

  // Return false if x falls outside the defined range, true otherwise
  AliFMDFlowBin* bin = GetBin(x);
  if (!bin) return kFALSE;
  bin->AddToHarmonic(phi, wp, wh);
  return kTRUE;
}

//____________________________________________________________________
UShort_t
AliFMDFlowBinned1D::Order() const 
{ 
  return GetBin(UShort_t(0))->Order(); 
}
//____________________________________________________________________
UShort_t
AliFMDFlowBinned1D::PsiOrder() const 
{ 
  return GetBin(UShort_t(0))->Order(); 
}

//____________________________________________________________________
void 
AliFMDFlowBinned1D::Event(ULong_t   n,  Double_t* phis, Double_t* xs, 
			  Double_t* wp, Double_t* wh)
{
  // Process a full event. 
  // Parameters: 
  //   phis   List of n phi=[0,2pi] angles 
  //   xs     List of n x values. 
  //   ws     Weights
  //   n      Size of phis and xs
  Begin();
  fSplitter->Event(phis, xs, n);
  for (ULong_t i = 0; i < n; i++) 
    AddToEventPlane(xs[i], phis[i], (wp ? wp[i] : 1), fSplitter->Select(i));
  for (ULong_t i = 0; i < n; i++) 
    AddToHarmonic(xs[i], phis[i], (wp ? wp[i] : 1), (wh ? wh[i] : 1));
  End();
}

//____________________________________________________________________
void 
AliFMDFlowBinned1D::Browse(TBrowser* b)
{
  // Browse this object. 
  b->Add(&fXAxis, "xaxis");
  if (fSplitter) b->Add(fSplitter, "Splitter");
  for (UInt_t i = 0; i < fXAxis.N(); i++) 
    b->Add(fBins[i], Form("bin_%03d", i));
}

//____________________________________________________________________
TH1*
AliFMDFlowBinned1D::MakeHistogram(UInt_t which, UInt_t what)
{
  // Make a histogram of some of the stuff in this object 

  // Some strings 
  const char* names[]   = { "Bare",  "Naive", "STAR", "TDR" };
  const char* hNames[]  = { "flow",  "resolution", "counts" }; 
  const char* hTitles[] = { "Flow",  "Resolution", "Counts" }; 
  const char* hUnits[]  = { "v_{%d}","<cos(%d(#Psi_{%d}-#Psi_{R}))>", "N" };
  AliFMDFlowBin::CorType types[] = { AliFMDFlowBin::kNone, 
				     AliFMDFlowBin::kNaive,
				     AliFMDFlowBin::kStar,
				     AliFMDFlowBin::kTdr };
  // Figure out how many things to draw
  UShort_t nm = 0;
  Short_t  sm = -1;
  for (UShort_t i = 0; i < 4; i++) { 
    if ((which & (1 << i)) != 0) {
      nm++;
      sm = i;
    }
  }
  if (what == AliFMDFlowBin::kCounts) { nm = 1; sm = 0; }

  TH1* h = 0;
  if (nm > 1) { 
    // Make 2D histogram 
    h = new TH2D(Form("%s_%s", GetName(), hNames[what]), 
		 Form("%s %s", GetTitle(), hTitles[what]), 
		 fXAxis.N(), fXAxis.Bins(), nm, 0, nm);

    // Set titles and such
    h->SetXTitle("x");
    h->SetYTitle("method");
    switch (what) { 
    case AliFMDFlowBin::kHarmonic:   
      h->SetZTitle(Form(hUnits[what], Order())); break;
    case AliFMDFlowBin::kResolution: 
      h->SetZTitle(Form(hUnits[what], Order(), PsiOrder())); break;
    default:
      h->SetZTitle(hUnits[what]); break;
    }
    h->GetYaxis()->SetNdivisions(nm+1, kFALSE);

    // Set Bin labels for the methods 
    UInt_t j = 0;
    for (UShort_t i = 0; i < 4; i++) 
      if (which & (1 << i)) h->GetYaxis()->SetBinLabel(++j, names[i]);
  }
  else {
    TString name(what == AliFMDFlowBin::kCounts ? 
		 Form("%s_%s", GetName(), hNames[what]) : 
		 Form("%s_%s_%s", GetName(), hNames[what], names[sm]));
    TString title(what == AliFMDFlowBin::kCounts ?
		  Form("%s %s", GetTitle(), hTitles[what]) : 
		  Form("%s %s %s", GetTitle(), hTitles[what], names[sm]));
    h = new TH1D(name.Data(), title.Data(), fXAxis.N(), fXAxis.Bins());
    h->SetXTitle("x");
    switch (what) { 
    case AliFMDFlowBin::kHarmonic:   
      h->SetYTitle(Form(hUnits[what], Order())); break;
    case AliFMDFlowBin::kResolution: 
      h->SetYTitle(Form(hUnits[what], Order(), PsiOrder())); break;
    default:
      h->SetYTitle(hUnits[what]); break;
    }
  }

  for (UShort_t i = 0; i < fXAxis.N(); i++) { 
    Double_t       x   = fXAxis.BinCenter(i);
    AliFMDFlowBin* bin = GetBin(x);
    Double_t       y=0, e2=0, dummy;
    if (bin->Counts() <= 0) continue;
    if (nm == 1) { 
      switch (what) { 
      case AliFMDFlowBin::kHarmonic:   
	if (bin->Correction(dummy, types[sm]) < .01) continue;
	y = bin->Value(e2, types[sm]); break;
      case AliFMDFlowBin::kResolution: 
	y = bin->Correction(e2, types[sm]); break;
      case AliFMDFlowBin::kCounts:     
	y = bin->Counts(); break;
      }
      h->SetBinContent(i+1, y);
      h->SetBinError(i+1, sqrt(e2));
      continue;
    }
    UInt_t j = 0;
    for (UShort_t k = 0; k < 4; k++)  { 
      if (!(which & (1 << k))) continue;
      switch (what) { 
      case AliFMDFlowBin::kHarmonic:   
	if (bin->Correction(dummy,types[k]) < .01) continue;
	y = bin->Value(e2, types[k]); break;
      case AliFMDFlowBin::kResolution: 
	y = bin->Correction(e2, types[k]); break;
      case AliFMDFlowBin::kCounts:     
	y = bin->Counts(); break;
      }
      h->SetBinContent(i+1, j+1, y);
      h->SetBinError(i+1, j+1, sqrt(e2));
      j++;
    }
  }
  h->SetLineColor(GetLineColor());
  h->SetLineWidth(GetLineWidth());
  h->SetLineStyle(GetLineStyle());
  h->SetFillColor(GetFillColor());
  h->SetFillStyle(GetFillStyle());
  h->SetMarkerColor(GetMarkerColor());
  h->SetMarkerStyle(GetMarkerStyle());
  h->SetMarkerSize(GetMarkerSize());
  h->SetDirectory(0);
  return h;
}
  

//____________________________________________________________________
void 
AliFMDFlowBinned1D::Draw(Option_t* option)
{
  // Draw the distribution of the harmonics or the event plane
  // resolution. 
  // Parameters: 
  //    option     String of options 
  // 
  // Options:
  //    b          Draw bare averages of cos(n(phi-Psi))
  //    n          Draw harmonic scaled by naive correction
  //    s          Draw harmonic scaled by STAR correction 
  //    t          Draw harmonic scaled by TDR correction (*)
  //    r          Draw resolution rather than harmonic 
  //    c          Draw counts rather than harmonic 
  //    h          Draw harmonics (*)
  //   
  // (*) This is the default value 
  // 
  // If more than one of b, n, s, or t is given, a 2D histogram is
  // drawn. 
  // 
  // Only one of r, c, h can be specified.  The first specified wins.
  //
  TString opt(option);
  opt.ToLower();
  TString dopt = "";
  Int_t idx = opt.Index(":");
  if (idx != kNPOS) {
    dopt = opt(idx+1, opt.Length()-idx-1);
    opt.Remove(idx, opt.Length()-idx);
  }
  UInt_t which = 0;
  if (opt.Contains("b")) which |= AliFMDFlowBin::kNone;
  if (opt.Contains("n")) which |= AliFMDFlowBin::kNaive;
  if (opt.Contains("s")) which |= AliFMDFlowBin::kStar;  
  if (opt.Contains("t")) which |= AliFMDFlowBin::kTdr;
  
  UInt_t what = AliFMDFlowBin::kHarmonic;
  if (opt.Contains("c")) what = AliFMDFlowBin::kCounts;
  if (opt.Contains("r")) what = AliFMDFlowBin::kResolution;

  TH1* h = MakeHistogram(which, what);
  h->Draw(dopt.Data());
}

  
  
//____________________________________________________________________
void 
AliFMDFlowBinned1D::Print(Option_t* option) const
{
  // Print information to standard output. 
  // Parameters: 
  //     option    String of options 
  // 
  // Options: 
  //     d         Print details 
  //     s         Print summary 
  // 
  TString opt(option);
  opt.ToLower();
  Bool_t det = opt.Contains("d");
  Bool_t sum = opt.Contains("s");
  if (det) { 
    for (UShort_t i = 0; i < fXAxis.N(); i++) { 
      Double_t x = fXAxis.BinCenter(i);
      std::streamsize         oldP = std::cout.precision(3);
      std::ios_base::fmtflags oldF = std::cout.setf(std::ios_base::fixed, 
						     std::ios_base::floatfield);
      std::cout << "x=" << std::setw(5) << x << std::endl;
      fBins[i]->Print();
      std::cout.precision(oldP);
      std::cout.setf(oldF, std::ios_base::floatfield);
    }
  }
  
  if (sum) { 
    UInt_t       nType = 4;
    const char*  names[] = { "Bare",    "Naive",    "STAR",    "TDR" };
    AliFMDFlowBin::CorType types[] = { AliFMDFlowBin::kNone, 
				       AliFMDFlowBin::kNaive, 
				       AliFMDFlowBin::kStar, 
				       AliFMDFlowBin::kTdr };
    std::cout << GetName() << " - " << GetTitle() << "\n"
	      << "    x";
    for (UInt_t i = 0; i < nType; i++) 
      std::cout << " | " << std::setw(6+6+5) << names[i];
    std::cout << "\n-----" << std::setfill('-');
    for (UInt_t i = 0; i < nType; i++) 
      std::cout << "-+-" <<  std::setw(6+6+5) << "-";
    std::cout << std::setfill(' ') << std::endl;
    
    std::streamsize         oldP = std::cout.precision(2);
    std::ios_base::fmtflags oldF = std::cout.setf(std::ios_base::fixed, 
						  std::ios_base::floatfield);
    for (UShort_t i = 0; i < fXAxis.N(); i++) { 
      Double_t x = fXAxis.BinCenter(i);
      std::cout << std::setprecision(2) << std::setw(5) << x << std::flush;
      for (UShort_t j = 0; j < nType; j++) { 
	Double_t e2v;
	Double_t v  = fBins[i]->Value(e2v, types[j]);
	std::cout << std::setprecision(3)    << " | "
		  << std::setw(6) << 100 * v << " +/- " 
		  << std::setw(6) << 100 * sqrt(e2v);
      }
      std::cout << std::endl;
    }
    std::cout.precision(oldP);
    std::cout.setf(oldF, std::ios_base::floatfield);
  }
}

    
//____________________________________________________________________
//
// EOF
//
