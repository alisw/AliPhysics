// -*- mode: C++ -*-
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
    @brief Declaration of Bessel functions */
#ifndef ALIFMDFLOWBESSEL_H
#define ALIFMDFLOWBESSEL_H
#include <cctype>
#include <Rtypes.h>

/** @defgroup z_bessel Bessel utility functions.  
    @brief This group contains a number of functions to calculate
    the modified Bessel function of the first kind @f$
    I_{\nu}(x)@f$ for integer and half-integer values of @f$
    \nu@f$, of any real number @f$ x@f$.
    
    The main entry of these functions is the function Inu, which
    in turn call I01, Iwhole, or Ihalf.  Inu always returns a list
    of the function values for the specified values of @f$ \nu@f$.
    
    The convinience functions I and DiffI returns a single value
    of @f$ I_{\nu}(x)@f$ or @f$ dI_{\nu}(x)/dx@f$ for a specified
    value of @f$ \nu@f$. 
    
    @example test_bessel.cxx */

/** Namespace for Bessel functions 
    @ingroup z_bessel */
namespace AliFMDFlowBessel
{
  /** @{ 
      @ingroup z_bessel */
  /** Compute the modified Bessel functions 
      @f[
      I_0(x) = \sum_{k=1}^\infty\frac{\left(\frac{x^2}{4}\right)^k}{(k!)^2}
      @f]
      and @f$ I_1(x)@f$, and their first derivatives  
      @param x   Argument @f$ x@f$ 
      @param bi0 On return, @f$ I_0(x)@f$ 
      @param di0 On return, @f$ I_0'(x)@f$ 
      @param bi1 On return, @f$ I_1(x)@f$ 
      @param di1 On return, @f$ I_1'(x)@f$ 
  */
  void I01(Double_t x, 
	   Double_t& bi0, Double_t& di0, Double_t& bi1, Double_t& di1);
  
  /** Compute the modified Bessel functions @f$ I_{n/2}(x)@f$ and their
      derivatives.  
      @param x  Argument. 
      @param n  Order
      @param bi On output, @f$ I_{1/2}(x), \ldots, I_{n/2}(x)@f$ for
      @f$ n > 0@f$ and @f$ I_{-n/2}(x), \ldots, I_{-1/2}(x)@f$ for
      @f$ n < 0@f$ 
      @param di On output, @f$ I_{1/2}'(x), \ldots, I_{n/2}'(x)@f$ for
      @f$ n > 0@f$ and @f$ I_{-n/2}'(x), \ldots, I_{-1/2}'(x)@f$ for
      @f$ n < 0@f$ 
      @return number of valid entries in @a bi and @a di */
  UInt_t Ihalf(Int_t n, Double_t x, Double_t* bi, Double_t* di);
  /** Compute the modified Bessel functions @f$ I_n(x)@f$ and their
      derivatives.  Note, that @f$ I_{-n}(x) = I_n(x)@f$ and 
      @f$ dI_{-n}(x)/dx = dI_{n}(x)/dx@f$ so this function can be used
      to evaluate @f$ I_n \forall n \in \mathcal{Z}@f$ 
      @param x  Argument. 
      @param n  Order
      @param bi On output, @f$ I_0(x), \ldots, I_{mn}(x)@f$ 
      @param di On output, @f$ I_0'(x), \ldots, I_{mn}'(x)@f$ 
      @return number of valid entries in @a bi and @a di */
  UInt_t Iwhole(UInt_t n, Double_t x, Double_t* bi, Double_t* di);

  /** Compute the modified Bessel functions @f$ I_{\nu}(x)@f$ and
      their derivatives for @f$ \nu = n_1, n_1+1, \ldots, n_2@f$ for
      and integer @f$ n_1, n_2@f$ or any half-integer @f$ n_1,
      n_2@f$. 
      @param x  Argument. 
      @param n1 Lower order
      @param n2 Upper order
      @param bi On output, @f$ I_{n_1}(x), \ldots, I_{n_2}(x)@f$ 
      @param di On output, @f$ I_{n_1}'(x), \ldots, I_{n_2}'(x)@f$ 
      @return number of valid entries in @a bi and @a di */
  UInt_t Inu(Double_t n1, Double_t n2, Double_t x, Double_t* bi, Double_t* di);

  /** Compute the modified Bessel function of the first kind 
      @f[
      I_n(x) = \left(\frac{x}{2}\right)^n
      \sum_{k=0}^\infty\frac{\left(\frac{x^2}{4}\right)^k}
      {k!\Gamma(n+k+1)}
      @f]
      for arbirary integer order @f$ n@f$ 
      @param n  Order
      @param x  Argument 
      @return @f$ I_n(x)@f$ */ 
  Double_t I(Double_t n, Double_t x);
      
  /** Compute the derivative of the modified Bessel function of the
      first kind 
      @f[ 
      \frac{dI_n(x)}{dx} = I_{n-1}(x) - \frac{n}{x} I_{n}(x)
      @f]
      for arbirary integer order @f$ n@f$
      @param n Order 
      @param x Argument 
      @return @f$ \frac{dI_n(x)}{dx}@f$ */ 
  Double_t DiffI(Double_t n, Double_t x);
  /** @} */
}

#endif
//
// EOF
//


