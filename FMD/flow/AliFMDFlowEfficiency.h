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
#ifndef ALIFMDFLOWEFFICIENCY_H
#define ALIFMDFLOWEFFICIENCY_H
#include <Rtypes.h>

/** @defgroup z_eff Efficiency calculations 
    @brief Functions to do efficiency calculations based on a
    Baysian analysis. 
*/
/** Namespace for efficency calculations 
    @ingroup z_eff */
namespace AliFMDFlowEfficiency 
{
  /** @{ 
      @ingroup z_eff */
  /** Calculate @f$ \log(\Gamma(z))\ \forall z>0@f$ 
      @param z Argument. 
      @return  @f$ \log(\Gamma(z))@f$ */ 
  Double_t LnGamma(Double_t z);
  
  /** Continued fraction evaluation by modified Lentz's method 
      used in calculation of incomplete  Beta function. 
      @param x argument. 
      @param a lower limit
      @param b upper limit
      @return incomplete Beta function evaluated at x */ 
  Double_t BetaCf(Double_t x, Double_t a, Double_t b);
  
  /** Calculates the incomplete Beta function @f$ I_x(a,b)@f$;
      this is the incomplete Beta function divided by the
      complete Beta function. 
      @param a Lower bound 
      @param b Upper bound 
      @param x Order 
      @return  @f$ I_x(a,b)@f$ */
  Double_t IBetaI(Double_t a, Double_t b, Double_t x);

  /** Calculates the fraction of the area under the curve 
      @f$ x^k (1-x)^{n-k}@f$ between @f$ x=a@f$ and @f$ x=b@f$ 
      @param a lower limit 
      @param b upper limit 
      @param k Parameter @f$ k@f$ 
      @param n Parameter @f$ n@f$ 
      @return The fraction under the curve */ 
  Double_t BetaAB(Double_t a, Double_t b, Int_t k, Int_t n);
  
  /** Integrates the Binomial distribution with parameters @a k and @a
      n, and determines the upper edge of the integration region, 
      starting at @a low, which contains probability content @a c.  If
      an upper limit is found, the value is returned. If no solution
      is found, -1 is returned. Check to see if there is any solution
      by verifying that the integral up to the maximum upper limit (1)
      is greater than c 
      @param low Where to start the integration from. 
      @param k   @a k parameter of the Binomial distribution. 
      @param n   @a N parameter of the Binomial distribution. 
      @param c   Wanted confidence limit (defaults to 68% - similar to
      @f$ 1\sigma@f$ of a Gaussian distribution)
      @return the upper limit of the confidence interval, or -1 in
      case of failure */
  Double_t SearchUpper(Double_t low, Int_t k, Int_t n, Double_t c=0.683);

  /** Integrates the Binomial distribution with parameters @a k and @a
      n, and determines the lower edge of the integration region, 
      ending at @a high, which contains probability content @a c.  If
      a lower limit is found, the value is returned. If no solution
      is found, -1 is returned. Check to see if there is any solution
      by verifying that the integral up to the maximum upper limit (1)
      is greater than c 
      @param high Where to end the integration at. 
      @param k    @a k parameter of the Binomial distribution. 
      @param n    @a N parameter of the Binomial distribution. 
      @param c    Wanted confidence limit (defaults to 68% - similar to
      @f$ 1\sigma@f$ of a Gaussian distribution) 
      @return the upper limit of the confidence interval, or -1 in
      case of failure */
  Double_t SearchLower(Double_t high, Int_t k, Int_t n, Double_t c=0.683);

  /** Numerical equation solver.  This includes root finding and
      minimum finding algorithms. Adapted from Numerical Recipes in
      C, 2nd edition. Translated to C++ by Marc Paterno
      @param ax   Left side of interval 
      @param bx   Middle of interval 
      @param cx   Right side of interval 
      @param tol  Tolerance 
      @param xmin On return, the value of @f$ x@f$ such that @f$
      f(x)@f$ i minimal. 
      @param k    Parameter of the Binomial distribution 
      @param n    Parameter of the Binomial distribution 
      @param c    Confidence level. 
      @return Minimum of @f$ f = f(x)@f$. */
  Double_t Brent(Double_t ax, Double_t bx, Double_t cx, Double_t& xmin, 
	       Int_t k, Int_t n, Double_t c=.683, Double_t tol=1e-9);
  
  /** Return the length of the interval starting at @a l that contains
      @a c of the @f$ x^k (1-x)^{n-k}@f$ distribution. If there is no
      sufficient interval starting at @a l, we return 2.0 
      @param l Lower bound 
      @param k Binomial parameter k 
      @param n Binomial parameter n 
      @param c Condifience level 
      @return Legnth of interval */
  Double_t Length(Double_t l, Int_t k, Int_t n, Double_t c=0.683);
    

  /** Calculate the shortest central confidence interval containing
      the required probability content @a c. Interval(low) returns the
      length of the interval starting at low that contains @a c
      probability. We use Brent's method, except in two special cases:
      when @a k=0, or when @a k=n 
      @author Marc Paterno
      @param k    Binomial parameter @a k 
      @param n    Binomial parameter @a n 
      @param low  On return, the lower limit
      @param high On return, the upper limit 
      @param c    Required confidence level (defaults to 68% - similar
      to @f$ 1\sigma@f$ of a Gaussian distribution)
      @return The mode */ 
  Double_t Interval(Int_t k, Int_t n, Double_t& low, Double_t& high, Double_t c=0.683);
  /** @} */
}


#endif
//
// EOF
//

