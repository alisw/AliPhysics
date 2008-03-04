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
    @brief Implementation of a Harmonic class */
//____________________________________________________________________
//
// Calculate the nth order harmonic. 
// Input is the phis of the observations,
// and the resolution of the event plane. 
// The class derives from AliFMDFlowStat to easy calculating the mean
// and the square variance of the harmonic. 
#include "flow/AliFMDFlowHarmonic.h"
#include "flow/AliFMDFlowUtil.h"
#include <TBrowser.h>
#include <iostream>
// #include <cmath>

//====================================================================
AliFMDFlowHarmonic::AliFMDFlowHarmonic(UShort_t n) 
  : fOrder(n),
    fPhi("phi", "#varphi-#Psi", 40, 0, 2*TMath::Pi()),
    fNPhi("nphi",Form("%d(#varphi-#Psi)",n),40,0,2*TMath::Pi()), 
    fWeight("weight", Form("cos(%d(#varphi-#Psi))", n), 100, -1, 1),
    fContrib("contrib", Form("(#varphi-#psi) vs cos(%d(#varphi-#Psi))", n), 
	     40, 0, 2 * TMath::Pi(), 100, -1, 1)
{
  fContrib.SetDirectory(0);
  fContrib.Sumw2();
  fContrib.SetXTitle("#varphi-#Psi");
  fContrib.SetYTitle(Form("cos(%d(#varphi-#Psi))", n));
  fPhi.SetDirectory(0);
  fPhi.Sumw2();
  fPhi.SetXTitle("#varphi-#Psi");
  fNPhi.SetDirectory(0);
  fNPhi.Sumw2();
  fNPhi.SetXTitle(Form("%d(#varphi-#Psi)", n));
  fWeight.SetDirectory(0);
  fWeight.Sumw2();
  fWeight.SetXTitle(Form("cos(%d(#varphi-#Psi))", n));
  fWeight.SetYTitle("#sum_iw_i");
} 

//____________________________________________________________________
AliFMDFlowHarmonic::AliFMDFlowHarmonic(const AliFMDFlowHarmonic& o)
  : AliFMDFlowStat(o), 
    fOrder(o.fOrder), 
    fPhi(o.fPhi),
    fNPhi(o.fNPhi),
    fWeight(o.fWeight),
    fContrib(o.fContrib)
{
  // Copy constructor 
  // Parameters: 
  //   o   Object to copy from 
  fContrib.SetDirectory(0);
  fContrib.Sumw2();
  fContrib.SetXTitle(Form("w_{i}cos(%d(#varphi-#Psi))", fOrder));
  fPhi.SetDirectory(0);
  fPhi.SetXTitle("#varphi-#Psi");
  fPhi.Sumw2();
  fNPhi.SetDirectory(0);
  fNPhi.SetXTitle(Form("%d(#varphi-#Psi)", fOrder));
  fNPhi.Sumw2();
  fWeight.SetDirectory(0);
  fWeight.Sumw2();
  fWeight.SetXTitle(Form("cos(%d(#varphi-#Psi))", fOrder));
  fWeight.SetYTitle("#sum_iw_i");
}

//____________________________________________________________________
AliFMDFlowHarmonic&
AliFMDFlowHarmonic::operator=(const AliFMDFlowHarmonic& o)
{
  // Assignment operator 
  // Parameters: 
  //   o   Object to assign from 
  // Return reference to this object. 
  AliFMDFlowStat::operator=(o);
  fOrder = o.fOrder;

  fContrib.Reset();
  fContrib.Add(&o.fContrib);
  fPhi.Reset();
  fPhi.Add(&o.fPhi);
  fNPhi.Reset();
  fNPhi.Add(&o.fNPhi);
  fWeight.Reset();
  fWeight.Add(&o.fWeight);
  return *this;
}

//____________________________________________________________________
void 
AliFMDFlowHarmonic::Browse(TBrowser* b)
{
  // Browse this object 
  // Parameters 
  //   b	Browser to use 
  // Return 
  //   nothing
  b->Add(&fContrib);
  b->Add(&fPhi);
  b->Add(&fNPhi);
  b->Add(&fWeight);
}

//____________________________________________________________________
void 
AliFMDFlowHarmonic::Add(Double_t phi, Double_t psi, Double_t weight) 
{ 
  // Add a data point. 
  // Parameters: 
  //    phi    Angle of this observation. 
  //    psi    Event plane of this observation
  //    weight The weight of this observation
  Double_t a       = NormalizeAngle(fOrder * (phi - psi));
  Double_t cosa    = TMath::Cos(a);
  Double_t contrib = cosa; // weight * cosa;
  AliFMDFlowStat::Add(contrib);
  fPhi.Fill(NormalizeAngle(phi-psi));
  fNPhi.Fill(a); 
  fWeight.Fill(cosa, weight);
  fContrib.Fill(a/*NormalizeAngle(phi-psi)*/,contrib);
}
//____________________________________________________________________
Double_t 
AliFMDFlowHarmonic::Value(Double_t r, Double_t er2, Double_t& e2) const 
{
  // The corrected value is given by 
  //
  //          v_n^obs
  //    v_n = -------
  //             R
  // 
  // where 
  // 
  //              1
  //    v_n^obs = - \sum_i(cos(n(\phi_i - \Psi)))
  //              N 
  // 
  // and R is the resolution 
  // 
  // The error on the corrected value is given by 
  //
  //                dv_n                    dv_n
  //    d^2v_n = (--------)^2 d^2v_n^obs + (----)^2 d^2R
  //              dv_n^obs                   dR
  //
  //             d^2v_n^obs R^2 + d^2R v_n^obs^2
  //           = -------------------------------
  //                         R^4 
  // 
  Double_t a = fAverage;
  Double_t v = a / r;
  if (fN != 0) { 
    Double_t s = fSqVar / fN;
    e2         = (s * r * r + er2 * a * a) / pow(r, 4);
  }
  else e2      = 0;
  return v;
}
//____________________________________________________________________
void
AliFMDFlowHarmonic::Print(Option_t* /*option*/) const 
{
  Double_t e2, er2 = 1, r = 1;
  Double_t v = Value(r, er2, e2);
  std::cout << v << " +/- " << e2 << std::endl;
  
} 

//____________________________________________________________________
//
// EOF
//
