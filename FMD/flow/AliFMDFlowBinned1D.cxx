/** @file 
    @brief Implementation of a 1-dimensional Flow "histogram" */
#include "flow/AliFMDFlowBinned1D.h"
#include "flow/AliFMDFlowBin.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <TBrowser.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>

//====================================================================
AliFMDFlowBinned1D::AliFMDFlowBinned1D(UShort_t order, 
				       UShort_t nxbins, 
				       Double_t* xbins) 
  : fXAxis(nxbins, xbins),
    fBins(0)
{
  UShort_t n = fXAxis.N();
  fBins   = new AliFMDFlowBin*[n];
  for (UInt_t i = 0; i < n; i++) fBins[i]= new AliFMDFlowBin(order);
}
//____________________________________________________________________
AliFMDFlowBinned1D::AliFMDFlowBinned1D(UShort_t order, 
				       const AliFMDFlowAxis& xaxis)
  : fXAxis(xaxis)
{
  UShort_t n = fXAxis.N();
  fBins   = new AliFMDFlowBin*[n];
  for (UInt_t i = 0; i < n; i++) fBins[i]= new AliFMDFlowBin(order);
}
//____________________________________________________________________
AliFMDFlowBinned1D::AliFMDFlowBinned1D(const AliFMDFlowBinned1D& o)
  : TObject(o), 
    fXAxis(o.fXAxis)
{
  UShort_t n = fXAxis.N();
  fBins   = new AliFMDFlowBin*[n];
  for (UInt_t i = 0; i < n; i++) fBins[i]= new AliFMDFlowBin(*(o.fBins[i]));
}
//____________________________________________________________________
AliFMDFlowBinned1D&
AliFMDFlowBinned1D::operator=(const AliFMDFlowBinned1D& o)
{
  if (fBins) { 
    for (UInt_t i = 0; i < fXAxis.N(); i++) delete fBins[i];
    delete [] fBins;
  }
  fXAxis     = o.fXAxis;
  UShort_t n = fXAxis.N();
  fBins   = new AliFMDFlowBin*[n];
  for (UInt_t i = 0; i < n; i++) fBins[i]= new AliFMDFlowBin(*(o.fBins[i]));
  return *this;
}
  
//____________________________________________________________________
AliFMDFlowBinned1D::~AliFMDFlowBinned1D()
{
  if (fBins) { 
    for (UInt_t i = 0; i < fXAxis.N(); i++) delete fBins[i];
    delete [] fBins;
  }
}

//____________________________________________________________________
AliFMDFlowBin* 
AliFMDFlowBinned1D::GetBin(UShort_t i) const
{
  if (i >= fXAxis.N()) return 0;
  return fBins[i];
}
//____________________________________________________________________
AliFMDFlowBin* 
AliFMDFlowBinned1D::GetBin(Double_t x) const
{
  Int_t i = fXAxis.FindBin(x);
  if (i < 0) return 0;
  UShort_t j = i;
  return GetBin(j);
}
  
//____________________________________________________________________
void 
AliFMDFlowBinned1D::Begin()
{
  for (UInt_t i = 0; i < fXAxis.N(); i++) fBins[i]->Begin();
}
//____________________________________________________________________
void 
AliFMDFlowBinned1D::End()
{
  for (UInt_t i = 0; i < fXAxis.N(); i++) fBins[i]->End();
}
//____________________________________________________________________
void 
AliFMDFlowBinned1D::Finish()
{
  for (UInt_t i = 0; i < fXAxis.N(); i++) fBins[i]->Finish();
}
//____________________________________________________________________
Bool_t 
AliFMDFlowBinned1D::AddToEventPlane(Double_t x, Double_t phi, Double_t w, Bool_t a)
{
  AliFMDFlowBin* bin = GetBin(x);
  if (!bin) return kFALSE;
  bin->AddToEventPlane(phi, w, a);
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDFlowBinned1D::AddToHarmonic(Double_t x, Double_t phi)
{
  AliFMDFlowBin* bin = GetBin(x);
  if (!bin) return kFALSE;
  bin->AddToHarmonic(phi);
  return kTRUE;
}

//____________________________________________________________________
void 
AliFMDFlowBinned1D::Event(Double_t* phis, Double_t* xs, Double_t* ws, ULong_t n)
{
  Begin();
  for (UInt_t i = 0; i < n; i++) 
    AddToEventPlane(xs[i], phis[i], (ws ? ws[i] : 1), 
		    Float_t(rand()) / RAND_MAX > 0.5);
  for (UInt_t i = 0; i < n; i++) 
    AddToHarmonic(xs[i], phis[i]);
  End();
}

//____________________________________________________________________
void 
AliFMDFlowBinned1D::Browse(TBrowser* b)
{
  b->Add(&fXAxis, "xaxis");
  for (UInt_t i = 0; i < fXAxis.N(); i++) 
    b->Add(fBins[i], Form("bin_%03d", i));
}

//____________________________________________________________________
void 
AliFMDFlowBinned1D::Draw(Option_t* option)
{
  TString opt(option);
  opt.ToLower();
  const char* names[] = { "Bare", "Naive", "STAR", "TDR" };
  const AliFMDFlowBin::CorType types[] = { AliFMDFlowBin::none, 
					   AliFMDFlowBin::naive, 
					   AliFMDFlowBin::star, 
					   AliFMDFlowBin::tdr };
  Bool_t meths[] = { opt.Contains("b"), 
		     opt.Contains("n"),
		     opt.Contains("s"), 
		     opt.Contains("t") };
  Bool_t res     = opt.Contains("r");
  UShort_t nm = 0;
  Short_t  sm = -1;
  for (UShort_t i = 0; i < 4; i++) { if (meths[i]) { nm++; sm = i; } }
  TH1* h = 0;
  if (nm > 1) { 
    h = new TH2D((res ? "res" : "flow"), (res ? "Resolution" : "Flow"), 
		 fXAxis.N(), fXAxis.Bins(), nm, 0, nm);
    h->SetXTitle("x");
    h->SetYTitle("method");
    h->GetYaxis()->SetNdivisions(nm+1, kFALSE);
    h->SetZTitle((res ? "<cos(n(#Psi_{m}-#Psi_{R}))>" : "v"));
    UInt_t j = 0;
    for (UShort_t i = 0; i < 4; i++) 
      if (meths[i]) h->GetYaxis()->SetBinLabel(++j, names[i]);
  }
  else {
    h = new TH1D(Form("%s_%s", (res ? "res" : "flow"), names[sm]),
		 Form("%s_%s", (res ? "Resolution" : "Flow"), names[sm]),
		 fXAxis.N(), fXAxis.Bins());
    h->SetXTitle("x");
    h->SetYTitle((res ? "<cos(n(#Psi_{m}-#Psi_{R}))>" : "v"));
  }

  for (UShort_t i = 0; i < fXAxis.N(); i++) { 
    Double_t       x   = fXAxis.BinCenter(i);
    AliFMDFlowBin* bin = GetBin(x);
    Double_t       v, e2;
    if (nm == 1) { 
      if (res) v = bin->Correction(e2, types[sm]);
      else     v = bin->Value(e2, types[sm]);
      h->SetBinContent(i+1, v);
      h->SetBinError(i+1, sqrt(e2));
      continue;
    }
    UInt_t j = 0;
    for (UShort_t k = 0; k < 4; k++)  { 
      if (!meths[k]) continue;
      if (res) v = bin->Correction(e2, types[k]);
      else     v = bin->Value(e2, types[k]);
      h->SetBinContent(i+1, j+1, v);
      h->SetBinError(i+1, j+1, sqrt(e2));
      j++;
    }
  }
  h->Draw(option);
}

  
  
//____________________________________________________________________
void 
AliFMDFlowBinned1D::Print(Option_t* option) const
{
  TString opt(option);
  opt.ToLower();
  Bool_t det = opt.Contains("d");
  Bool_t sum = opt.Contains("s");
  if (det) { 
    for (UShort_t i = 0; i < fXAxis.N(); i++) { 
      Double_t x = fXAxis.BinCenter(i);
      std::streamsize         old_p = std::cout.precision(3);
      std::ios_base::fmtflags old_f = std::cout.setf(std::ios_base::fixed, 
						     std::ios_base::floatfield);
      std::cout << "x=" << std::setw(5) << x << std::endl;
      fBins[i]->Print();
      std::cout.precision(old_p);
      std::cout.setf(old_f, std::ios_base::floatfield);
    }
  }
  
  if (sum) { 
    UInt_t       nType = 4;
    const char*  names[] = { "Bare",    "Naive",    "STAR",    "TDR" };
    AliFMDFlowBin::CorType types[] = { AliFMDFlowBin::none, 
				       AliFMDFlowBin::naive, 
				       AliFMDFlowBin::star, 
				       AliFMDFlowBin::tdr };
    std::cout << "    x";
    for (UInt_t i = 0; i < nType; i++) 
      std::cout << " | " << std::setw(6+6+5) << names[i];
    std::cout << "\n-----" << std::setfill('-');
    for (UInt_t i = 0; i < nType; i++) 
      std::cout << "-+-" <<  std::setw(6+6+5) << "-";
    std::cout << std::setfill(' ') << std::endl;
    
    std::streamsize         old_p = std::cout.precision(2);
    std::ios_base::fmtflags old_f = std::cout.setf(std::ios_base::fixed, 
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
    std::cout.precision(old_p);
    std::cout.setf(old_f, std::ios_base::floatfield);
  }
}

    
//____________________________________________________________________
//
// EOF
//
