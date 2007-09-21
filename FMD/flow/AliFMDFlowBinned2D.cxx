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
    @brief Implementation of a 2-dimensional Flow "histogram" */
//____________________________________________________________________ 
//
// A histogram of flow bins.  The axis can by anything
// (pseudo-rapidity, transvers momentum) - there's no assumption on
// what is the basis of the histogram.  The method Event can be used
// to calculate everything in one go.   Alternatively, one can use the
// methods AddToEventPlane and AddToHarmonic.  See also the example
// TestFlow.C 
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
  // Constructor 
  // Parameters: 
  //   Order	Order 
  //   nxbins   Number of X bins 
  //   xbins    X Bin borders 
  //   nybins   Number of Y bins 
  //   ybins    Y Bin borders 
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
  // Constructor 
  // Parameters: 
  //   Order	Order 
  //   xaxis    X axis object
  //   yaxis    Y axis object
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
  // Copy constructor 
  // Parameters: 
  //   o   Object to copy from 
  UShort_t n = fXAxis.N() * fYAxis.N();
  fBins   = new AliFMDFlowBin*[n];
  for (UInt_t i = 0; i < n; i++) fBins[i]= new AliFMDFlowBin(*(o.fBins[i]));
}
//____________________________________________________________________
AliFMDFlowBinned2D&
AliFMDFlowBinned2D::operator=(const AliFMDFlowBinned2D& o)
{
  // Assignment operator
  // Parameters: 
  //   o Object to assign from 
  // 
  // Returns reference to this object 
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
  // Destructor 
  // Parameters: 
  //    none 
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
  // Get the ith,jth bin 
  // Parameters: 
  //    i  X Bin number 
  //    j  Y Bin number 
  // 
  // Return pointer to bin, or null. 
  if (i >= fXAxis.N() || j >= fYAxis.N()) return 0;
  return fBins[i * fYAxis.N() + j];
}
//____________________________________________________________________
AliFMDFlowBin* 
AliFMDFlowBinned2D::GetBin(Double_t x, Double_t y) const
{
  // Get the bin that contains x,y
  // Parameters: 
  //    x  X axis value
  //    y  X axis value
  // 
  // Return pointer to bin, or null. 
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
  // Called at the beginning of an event
  // Parameters: 
  //   none
  UInt_t n = fXAxis.N() * fYAxis.N();
  for (UInt_t i = 0; i < n; i++) fBins[i]->Begin();
}
//____________________________________________________________________
void 
AliFMDFlowBinned2D::End()
{
  // Called at the end of an event
  // Parameters: 
  //   none
  UInt_t n = fXAxis.N() * fYAxis.N();
  for (UInt_t i = 0; i < n; i++) fBins[i]->End();
}
//____________________________________________________________________
void 
AliFMDFlowBinned2D::Finish()
{
  // Called at the end of an job
  // Parameters: 
  //   none
  UInt_t n = fXAxis.N() * fYAxis.N();
  for (UInt_t i = 0; i < n; i++) fBins[i]->Finish();
}
//____________________________________________________________________
Bool_t 
AliFMDFlowBinned2D::AddToEventPlane(Double_t x, Double_t y, Double_t phi, 
				Double_t w, Bool_t a)
{
  // Called to add a contribution to the event plane 
  // Parameters:
  //    x   X Bin value to fill into 
  //    y   Y Bin value to fill into 
  //    w   Weight
  //    phi The angle phi in radians 
  //    a   If true, add to sub-event A, otherwise sub-event B
  // 
  // Return false if (x,y) falls outside the defined range, true otherwise
  AliFMDFlowBin* bin = GetBin(x, y);
  if (!bin) return kFALSE;
  bin->AddToEventPlane(phi, w, a);
  return kTRUE;
}

//____________________________________________________________________
Bool_t 
AliFMDFlowBinned2D::AddToHarmonic(Double_t x, Double_t y, Double_t phi)
{
  // Called to add a contribution to the harmonic
  // Parameters: 
  //    x   X Bin value to fill into 
  //    y   Y Bin value to fill into 
  //    phi The angle phi in radians
  // 
  // Return false if (x,y) falls outside the defined range, true otherwise
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
  // Process a full event. 
  // Parameters: 
  //   phis   List of n phi=[0,2pi] angles 
  //   xs     List of n x values. 
  //   ys     List of n y values. 
  //   ws     Weights
  //   n      Size of phis and xs
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
  // Browse this object
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
