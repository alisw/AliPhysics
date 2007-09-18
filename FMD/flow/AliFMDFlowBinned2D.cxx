/** @file 
    @brief Implementation of a 2-dimensional Flow "histogram" */
#include "flow/AliFMDFlowBinned2D.h"
#include "flow/AliFMDFlowBin.h"
#include <cmath>
#include <cstdlib>
#include <TString.h>
#include <TBrowser.h>

//====================================================================
AliFMDFlowBinned2D::AliFMDFlowBinned2D(UShort_t order, 
			 UShort_t nxbins, Double_t* xbins,
			 UShort_t nybins, Double_t* ybins) 
  : fXAxis(nxbins, xbins),
    fYAxis(nybins, ybins),
    fBins(0)
{
  UInt_t n = fXAxis.N() * fYAxis.N();
  fBins   = new AliFMDFlowBin*[n];
  for (UInt_t i = 0; i < n; i++) fBins[i]= new AliFMDFlowBin(order);
}
//____________________________________________________________________
AliFMDFlowBinned2D::AliFMDFlowBinned2D(UShort_t order, 
			 const AliFMDFlowAxis&    xaxis, 
			 const AliFMDFlowAxis&    yaxis)
  : fXAxis(xaxis), 
    fYAxis(yaxis)
{
  UShort_t n = fXAxis.N() * fYAxis.N();
  fBins   = new AliFMDFlowBin*[n];
  for (UInt_t i = 0; i < n; i++) fBins[i]= new AliFMDFlowBin(order);
}
//____________________________________________________________________
AliFMDFlowBinned2D::AliFMDFlowBinned2D(const AliFMDFlowBinned2D& o)
  : TObject(o), 
    fXAxis(o.fXAxis), 
    fYAxis(o.fYAxis)
{
  UShort_t n = fXAxis.N() * fYAxis.N();
  fBins   = new AliFMDFlowBin*[n];
  for (UInt_t i = 0; i < n; i++) fBins[i]= new AliFMDFlowBin(*(o.fBins[i]));
}
//____________________________________________________________________
AliFMDFlowBinned2D&
AliFMDFlowBinned2D::operator=(const AliFMDFlowBinned2D& o)
{
  if (fBins) { 
    UInt_t n = fXAxis.N() * fYAxis.N();
    for (UInt_t i = 0; i < n; i++) delete fBins[i];
    delete [] fBins;
  }
  fXAxis     = o.fXAxis;
  UShort_t n = fXAxis.N() * fYAxis.N();
  fBins   = new AliFMDFlowBin*[n];
  for (UInt_t i = 0; i < n; i++) fBins[i]= new AliFMDFlowBin(*(o.fBins[i]));
  return *this;
}

//____________________________________________________________________
AliFMDFlowBinned2D::~AliFMDFlowBinned2D()
{
  if (fBins) { 
    UInt_t n = fXAxis.N() * fYAxis.N();
    for (UInt_t i = 0; i < n; i++) delete fBins[i];
    delete [] fBins;
  }
}
//____________________________________________________________________
AliFMDFlowBin* 
AliFMDFlowBinned2D::GetBin(UShort_t i, UShort_t j) const
{
  if (i >= fXAxis.N() || j >= fYAxis.N()) return 0;
  return fBins[i * fYAxis.N() + j];
}
//____________________________________________________________________
AliFMDFlowBin* 
AliFMDFlowBinned2D::GetBin(Double_t x, Double_t y) const
{
  Int_t i = fXAxis.FindBin(x);
  if (i < 0) return 0;
  Int_t j = fYAxis.FindBin(y);
  if (j < 0) return 0;
  UShort_t k = i;
  UShort_t l = j;
  return GetBin(k, l);
}
//____________________________________________________________________
void 
AliFMDFlowBinned2D::Begin()
{
  UInt_t n = fXAxis.N() * fYAxis.N();
  for (UInt_t i = 0; i < n; i++) fBins[i]->Begin();
}
//____________________________________________________________________
void 
AliFMDFlowBinned2D::End()
{
  UInt_t n = fXAxis.N() * fYAxis.N();
  for (UInt_t i = 0; i < n; i++) fBins[i]->End();
}
//____________________________________________________________________
void 
AliFMDFlowBinned2D::Finish()
{
  UInt_t n = fXAxis.N() * fYAxis.N();
  for (UInt_t i = 0; i < n; i++) fBins[i]->Finish();
}
//____________________________________________________________________
Bool_t 
AliFMDFlowBinned2D::AddToEventPlane(Double_t x, Double_t y, Double_t phi, 
				Double_t w, Bool_t a)
{
  AliFMDFlowBin* bin = GetBin(x, y);
  if (!bin) return kFALSE;
  bin->AddToEventPlane(phi, w, a);
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDFlowBinned2D::AddToHarmonic(Double_t x, Double_t y, Double_t phi)
{
  AliFMDFlowBin* bin = GetBin(x, y);
  if (!bin) return kFALSE;
  bin->AddToHarmonic(phi);
  return kTRUE;
}

//____________________________________________________________________
void 
AliFMDFlowBinned2D::Event(Double_t* phis, Double_t* xs, Double_t* ys, 
			  Double_t* ws, ULong_t n)
{
  Begin();
  for (UInt_t i = 0; i < n; i++) 
    AddToEventPlane(xs[i], ys[i], phis[i], (ws ? ws[i] : 1), 
		    Float_t(rand()) / RAND_MAX > 0.5);
  for (UInt_t i = 0; i < n; i++) 
    AddToHarmonic(xs[i], ys[i], phis[i]);
  End();
}

//____________________________________________________________________
void 
AliFMDFlowBinned2D::Browse(TBrowser* b)
{
  b->Add(&fXAxis, "xaxis");
  b->Add(&fYAxis, "yaxis");
  for (UInt_t i = 0; i < fXAxis.N(); i++) { 
    for (UInt_t j = 0; i < fYAxis.N(); j++) { 
      b->Add(fBins[i*fXAxis.N()+j], Form("bin_%03d_%03d", i, j));
    }
  }
}

//____________________________________________________________________
//
// EOF
//
