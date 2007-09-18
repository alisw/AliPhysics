// -*- C++ -*- 
/** @file 
    @brief Declaration of an EventPlane class */
#ifndef FLOW_EVENTPLANE_H
#define FLOW_EVENTPLANE_H
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
  AliFMDFlowEventPlane(UShort_t m) 
    : fSumSinMPhi(0), 
      fSumCosMPhi(0),
      fOrder(m),
      fCache(0)
  { Clear(); }
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
  Double_t fSumSinMPhi;
  /** @f$ \sum_i w_i \cos(k \varphi_i)@f$ */
  Double_t fSumCosMPhi;
  /** Order */
  UShort_t fOrder;
  /** Cache of Psi */
  mutable Double_t fCache;
  /** Define for ROOT I/O */
  ClassDef(AliFMDFlowEventPlane,1); 
};


#endif
//
// EOF
//
