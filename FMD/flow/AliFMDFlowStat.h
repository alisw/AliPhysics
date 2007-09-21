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
#ifndef ALIFMDFLOWSTAT_H
#define ALIFMDFLOWSTAT_H
// #include <cmath>
#include <TMath.h>
#include <TObject.h>

//______________________________________________________
/** @class AliFMDFlowStat flow/AliFMDFlowStat.h <flow/AliFMDFlowStat.h>
    @brief Single variable statistics. 
    @ingroup u_utils 

    This class calculates the single variable statistics of a random
    variable @f$ x@f$. 
    
    The average @f$\bar{x}@f$ after @f$n@f$ observations is calculated as 
    @f[ 
    \bar{x}^{(n)} = \left\{\begin{array}{rl} 
    x^{(n)} & n = 1\\
    \left(1-\frac1n\right)\bar{x}^{(n-1)} + \frac{x^{(n)}}{n} & n > 1
    \end{array}\right. 
    @f]
    and the square variance @f$s^2@f$ after @f$ n@f$ observations is
    given by  
    @f[ 
    {s^{2}}^{(n)} = \left\{\begin{array}{rl}
    0 & n = 1\\ 
    \left(1-\frac1n\right){s^{2}}^{(n-1)} + 
    \frac{\left(x^{(n)}-\bar{x}^{(n)}\right)^2}{n-1} & n > 1
    \end{array}\right.
    @f] 
*/
class AliFMDFlowStat : public TObject
{
public:
  /** Constructor */
  AliFMDFlowStat() : fAverage(0), fSqVar(0), fN(0) {}
  /** Copy constructor 
      @param o Object to copy from. */
  AliFMDFlowStat(const AliFMDFlowStat& o);
  /** Assuignement operator
      @param o Object to assign from. 
      @return Reference to this */
  AliFMDFlowStat& operator=(const AliFMDFlowStat& o);
  /** Destructor */
  virtual ~AliFMDFlowStat() {}
  /** Reset */
  void Clear(Option_t* option="");
  /** Add another data point */
  void Add(Double_t x);
  /** Get the average */
  Double_t Average() const { return fAverage; }
  /** Get the square variance */
  Double_t SqVar() const { return fSqVar; } 
  /** Get the number of data points */
  ULong_t N() const { return fN; }
  /**  Given the two partial sample of the same stocastic variable, of 
       size @f$ n_1@f$ and @f$ n_2@f$, with averages @f$ m_1@f$ and
       @f$ m2@f$, and square variances, @f$ s_1^2@f$ and @f$ s_2^2@f$,
       calculate the full sample average and square variance: 
       @f[ 
       m = \frac{n_1 m_1 + n_2 m_2}{n_1 + n_2} 
       @f] 
       and 
       @f[ 
       s^2 = \frac{n_1 * s_1^2 + n_2 s_2^2}{n_1 + n_2} 
             + \frac{n_1 n_2 (m_1 - m_2)^2}{(n_1 + n_2)^2}
       @f]
       @param o Other statistics object to add 
       @return Reference to this */
  AliFMDFlowStat& operator+=(const AliFMDFlowStat& o);
protected:
  /** Average */
  Double_t        fAverage; // Running 
  /** Square variance */
  Double_t        fSqVar; // Square variance 
  /** Number of data points */
  ULong_t fN; // Number of data points
  /** Define for ROOT I/O */
  ClassDef(AliFMDFlowStat,1);
};

//__________________________________________________________________
inline 
AliFMDFlowStat::AliFMDFlowStat(const AliFMDFlowStat& o)
  : TObject(o), 
    fAverage(o.fAverage), 
    fSqVar(o.fSqVar), 
    fN(o.fN)
{}
//__________________________________________________________________
inline AliFMDFlowStat&
AliFMDFlowStat::operator=(const AliFMDFlowStat& o)
{
  fAverage = o.fAverage;
  fSqVar   = o.fSqVar;
  fN       = o.fN;
  return *this;
}
//__________________________________________________________________
inline AliFMDFlowStat&
AliFMDFlowStat::operator+=(const AliFMDFlowStat& o)
{
  if (fN + o.fN <= 0) return *this;
  
  Double_t m = (fN * o.fAverage + o.fN * o.fAverage)/(fN + o.fN);
  Double_t s = ((fN * o.fSqVar + o.fN * o.fSqVar)/(fN + o.fN)
		+ (TMath::Power(fAverage - o.fAverage, 2) * fN * o.fN) 
		/ TMath::Power(fN + o.fN, 2));
  fAverage = m;
  fSqVar   = s;
  fN       = fN + o.fN;
  return *this;
}
//__________________________________________________________________
inline void 
AliFMDFlowStat::Clear(Option_t*) 
{ 
  fAverage = fSqVar = 0; 
  fN = 0; 
}

//__________________________________________________________________
inline void 
AliFMDFlowStat::Add(Double_t x) 
{ 
  if (TMath::IsNaN(x) || !TMath::Finite(x)) return;
  fN++;
  if (fN == 1) { 
    fAverage = x;
    fSqVar   = 0;
    return;
  }
  Double_t c = (1. - 1. / fN);
  fAverage  *= c;
  fAverage  += x / fN;
  fSqVar    *= c;
  fSqVar    += TMath::Power(x - fAverage, 2) / (fN - 1);
}


#endif
//
// EOF
//
