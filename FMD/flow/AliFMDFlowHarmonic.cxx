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
// #include <cmath>

//====================================================================
AliFMDFlowHarmonic::AliFMDFlowHarmonic(const AliFMDFlowHarmonic& o)
  : AliFMDFlowStat(o), 
    fOrder(o.fOrder)
{
  // Copy constructor 
  // Parameters: 
  //   o   Object to copy from 
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
  return *this;
}

//____________________________________________________________________
void 
AliFMDFlowHarmonic::Add(Double_t phi, Double_t psi) 
{ 
  // Add a data point. 
  // Parameters: 
  //    phi  Angle. 
  //    psi  Event plane 
  Double_t a       = NormalizeAngle(fOrder * (phi - psi));
  Double_t contrib = cos(a);
  AliFMDFlowStat::Add(contrib);
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
  Double_t s = fSqVar / fN;
  e2       = (s * r * r + er2 * a * a) / pow(r, 4);
  return v;
}
//____________________________________________________________________
//
// EOF
//
