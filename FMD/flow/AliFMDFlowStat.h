// -*- mode: C++ -*- 
#ifndef FLOW_STAT_H
#define FLOW_STAT_H
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
protected:
  /** Average */
  Double_t        fAverage;
  /** Square variance */
  Double_t        fSqVar;
  /** Number of data points */
  ULong_t fN;
  /** Define for ROOT I/O */
  ClassDef(AliFMDFlowStat,1);
};

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
