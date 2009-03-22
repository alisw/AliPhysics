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
    @brief Implementation of an EventPlane class */
//____________________________________________________________________
//
// Class to determine the event plane 
// 
// The event plane is calculated as 
// 
//    Psi_n = 1/n * atan((sum_i(w_i sin(n phi_i)))
//                        sum_i(w_i cos(n phi_i))))
//
// where i runs over all observations of phi in an event, and 
// w_i is the weight of the ith observation of phi
#include "flow/AliFMDFlowEventPlane.h"
#include "flow/AliFMDFlowUtil.h"
#include <TMath.h>
#include <TBrowser.h>
// #include <iostream>
// #include <cmath>
#ifndef _GNU_SOURCE
//#warning Using private implementation of sincos
extern "C" 
{
  /** Function to caculate @f$ \sin(a), \cos(a)@f$ in one go.  Note,
      that with GCC, this is a built-in, and we only have this to
      serve as a replacement function 
      @ingroup utils */
  inline void sincos(Double_t a, Double_t* sina, Double_t* cosa) {
    *sina = sin(a);
    *cosa = cos(a);
  }
}
#endif

//====================================================================
AliFMDFlowEventPlane::AliFMDFlowEventPlane(UShort_t m) 
  : fSumSinMPhi(0), 
    fSumCosMPhi(0),
    fOrder(m),
    fCache(0), 
    fSum("sum", Form("#sumw#cos(%d#varphi) vs #sum#sin(%d#varphi)", m, m), 
	 100, -1.1, 1.1, 100, -1.1, 1.1), 
    fPsi("psi", Form("#Psi_{%d}", m), 80, 0, 2*TMath::Pi()),
    fScale(0)
{ 
  Clear(); 
  fSum.SetDirectory(0);
  fSum.SetXTitle(Form("#sum_{i}w_{i}cos(%d#varphi_{i})", fOrder));
  fSum.SetYTitle(Form("#sum_{i}w_{i}sin(%d#varphi_{i})", fOrder));
  fSum.SetMarkerStyle(20);
  fPsi.SetDirectory(0);
  fPsi.SetXTitle(Form("#Psi_{%d}", fOrder));
}

//____________________________________________________________________
AliFMDFlowEventPlane::AliFMDFlowEventPlane(const AliFMDFlowEventPlane& o)
  : TObject(o), 
    fSumSinMPhi(o.fSumSinMPhi),
    fSumCosMPhi(o.fSumCosMPhi),
    fOrder(o.fOrder),
    fCache(-1), 
    fSum(o.fSum),
    fPsi(o.fPsi),
    fScale(o.fScale)
{
  // copy cosntructor 
  // Parameters 
  //   o  Object to copy from. 
  fSum.SetDirectory(0);
  fSum.SetXTitle(Form("#sum_{i}w_{i}cos(%d#varphi_{i})", fOrder));
  fSum.SetYTitle(Form("#sum_{i}w_{i}sin(%d#varphi_{i})", fOrder));
  fSum.SetMarkerStyle(20);
  fPsi.SetDirectory(0);
  fPsi.SetXTitle(Form("#Psi_{%d}", fOrder));
}
//____________________________________________________________________
AliFMDFlowEventPlane&
AliFMDFlowEventPlane::operator=(const AliFMDFlowEventPlane& o)
{
  // Assignment operator. 
  // Parameters: 
  //  o Object to assign from 
  fSumSinMPhi = o.fSumSinMPhi;
  fSumCosMPhi = o.fSumCosMPhi;
  fOrder      = o.fOrder;
  fCache      = -1;
  fScale      = o.fScale;

  fSum.Reset();
  fSum.Add(&o.fSum);
  fPsi.Reset();
  fPsi.Add(&o.fPsi);

  return *this;
}

//____________________________________________________________________
void 
AliFMDFlowEventPlane::Clear(Option_t*) 
{ 
  // clear internal variables. 
  fSumSinMPhi = 0;
  fSumCosMPhi = 0;
  fCache      = -1;
  fScale      = 0;
}
//____________________________________________________________________
void 
AliFMDFlowEventPlane::Add(Double_t phi, Double_t weight) 
{ 
  // Add a data point 
  // Parameters: 
  //   phi     The angle phi in[0,2pi]
  //   weight  The weight 
  Double_t a = NormalizeAngle(fOrder * phi);
  Double_t s, c;
  sincos(a, &s, &c);
  if (TMath::IsNaN(s) || !TMath::Finite(s) || 
      TMath::IsNaN(c) || !TMath::Finite(s)) return;
  if (weight == 0) return;
  fScale += 1./weight;
  fSumSinMPhi += weight * s;
  fSumCosMPhi += weight * c;
}

//____________________________________________________________________
void
AliFMDFlowEventPlane::End() 
{
  Double_t r2 = fSumCosMPhi*fSumCosMPhi + fSumSinMPhi*fSumSinMPhi;
  Double_t r  = (r2 < 0 ? 1 : TMath::Sqrt(r2));
  fSum.Fill(fSumCosMPhi/(r!=0?r:1),fSumSinMPhi/(r!=0?r:1));
  fPsi.Fill(Psi());
}

//____________________________________________________________________
void
AliFMDFlowEventPlane::Browse(TBrowser* b)
{
  b->Add(&fSum);
  b->Add(&fPsi);
}

//____________________________________________________________________
Double_t 
AliFMDFlowEventPlane::Psi() const  
{ 
  // Get the event plane 
  // Parameters: 
  //   none
  if (fCache < 0) { 
    fCache = DoPsi(fSumSinMPhi, fSumCosMPhi);
  }
  return fCache;
}
//____________________________________________________________________
Double_t 
AliFMDFlowEventPlane::Psi(Double_t phi, Double_t w) const  
{ 
  // Get the event plane angle Psi_k disregarding the contribution
  // from the observation phi with weight w.  This is to avoid
  // auto-correlations  
  // 
  // Parameters: 
  //  phi   The observation  phi
  //  w     The weight w of the obervation. 
  // 
  // Returns The event plane angle Psi with out the contribution from
  // phi_i
  Double_t a = NormalizeAngle(fOrder * phi);
  Double_t s, c;
  sincos(a, &s, &c);
  if (TMath::IsNaN(s) || !TMath::Finite(s) || 
      TMath::IsNaN(c) || !TMath::Finite(s)) return Psi();
  Double_t psi = DoPsi(fSumSinMPhi - w * s, fSumCosMPhi - w * c);
  return psi;
}

//____________________________________________________________________
Double_t 
AliFMDFlowEventPlane::DoPsi(Double_t sumsin, Double_t sumcos) const
{
  // Calculate the event plane 
  // Parameters: 
  //   sumsin    Sum of sines 
  //   sumcos    Sum of cosines 
  Double_t psi = 0;
  // Make sure we get an angle everywhere 
  if      (sumcos != 0) psi =  atan2(sumsin, sumcos);
  else if (sumsin == 0) psi =  0;
  else if (sumsin >  0) psi =  M_PI / 2;
  else                  psi = -M_PI / 2;
  psi =  NormalizeAngle(psi);
  psi /= fOrder;
  return psi;
}

//____________________________________________________________________
//
// EOF
//

