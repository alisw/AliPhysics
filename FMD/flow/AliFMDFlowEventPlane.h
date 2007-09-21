// -*- C++ -*- 
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
    @brief Declaration of an EventPlane class */
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
#ifndef ALIFMDFLOWEVENTPLANE_H
#define ALIFMDFLOWEVENTPLANE_H
#include <TObject.h>

//______________________________________________________
/** @class AliFMDFlowEventPlane flow/AliFMDFlowEventPlane.h <flow/AliFMDFlowEventPlane.h>
    @brief Class to determine the event plane 

    The event plane is calculated as 
    @f[ 
    \Psi_n = \frac1n\tan^{-1}\left[\frac{\sum_i(w_i\sin(n\varphi_i))}
    {\sum_i(w_i\cos(n\varphi_i))}\right]
    @f]
    where @f$ i @f$ runs over all observations of @f$\varphi@f$ in an
    event, and @f$ w_i@f$ is the weight of the @f$ i@f$ observation of
    @f$ \varphi@f$

    @ingroup a_basic
*/
class AliFMDFlowEventPlane : public TObject
{
public:
  /** Constructor 
      @param m Harmonic number */
  AliFMDFlowEventPlane(UShort_t m=0) 
    : fSumSinMPhi(0), 
      fSumCosMPhi(0),
      fOrder(m),
      fCache(0)
  { Clear(); }
  /** Copy constructor. 
      @param o Object to copy from */ 
  AliFMDFlowEventPlane(const AliFMDFlowEventPlane& o);
  /** Assignement operator. 
      @param o Object to copy from 
      @return Reference to this */
  AliFMDFlowEventPlane& operator=(const AliFMDFlowEventPlane& o);
  /** Destructor */
  ~AliFMDFlowEventPlane() {} 
  /** Clear it */
  void Clear(Option_t* option="");
  /** Add a data point 
      @param phi The angle @f$\varphi\in[0,2\pi]@f$
      @param weight The weight */
  void Add(Double_t phi, Double_t weight=1);
  /** Get the event plane 
      @return @f$ \Psi_k@f$ */
  Double_t Psi() const;
  /** Get the event plane angle @f$ \Psi_k@f$ @e disregarding the
      contribution from the observation @f$ \varphi_i@f$ with weight
      @f$ w_i@f$.  This is to avoid auto-correlations 
      @param phi The observation  @f$ \varphi_i@f$
      @param w   The weight @f$ w_i@f$ of the obervation. 
      @return The event plane angle @f$ \Psi_k@f$ with out the
      contribution from @f$ \varphi_i@f$ */
  Double_t Psi(Double_t phi, Double_t w=1) const;
  /** Get the harmnic order 
      @return @f$ k@f$  */
  UShort_t Order() const { return fOrder; }
protected:
  /** Utility function to calculate @f$ \Psi@f$ from the sum of
      sines and cosines. 
      @param sumsin Sum of sines 
      @param sumcos Sum of cosine. 
      @return @f$ \Psi@f$ */
  Double_t DoPsi(Double_t sumsin, Double_t sumcos) const;
  /** @f$ \sum_i w_i \sin(k \varphi_i)@f$ */
  Double_t fSumSinMPhi; // Sum of contributions 
  /** @f$ \sum_i w_i \cos(k \varphi_i)@f$ */
  Double_t fSumCosMPhi; // Sum of contributions 
  /** Order */
  UShort_t fOrder; // Order 
  /** Cache of Psi */
  mutable Double_t fCache; // Cache of calculated value 
  /** Define for ROOT I/O */
  ClassDef(AliFMDFlowEventPlane,1); 
};


#endif
//
// EOF
//
