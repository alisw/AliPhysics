// -*- mode: C++ -*-
/** @file 
    @brief Declaration and implementation of classes to deal with a
    well-known event plane @f$ \Psi_R@f$. */
#ifndef FLOW_TRUE_H
#define FLOW_TRUE_H
#include <flow/AliFMDFlowStat.h>
#include <flow/AliFMDFlowBin.h>
#include <flow/AliFMDFlowBinned1D.h>


/** @defgroup x_true Classes for handling Flow when the event plane is
    known. */
/** @class TrueBin flow/True.h <flow/True.h>
    @brief A specialised Bin of flow in case the event plane is 
    well-known. 
    @ingroup x_true */
struct AliFMDFlowTrueBin : public AliFMDFlowBin
{
  /** Constructor */ 
  AliFMDFlowTrueBin(UShort_t order) : 
    AliFMDFlowBin(order, 1), 
    fPsiR(0), 
    fResReal()
  {}
  /** Set the well-known event plane angle @f$ \Psi_R@f$ for this
      event. 
      @param psi @f$ \Psi_R@f$ */ 
  void SetPsi(Double_t psi) { fPsiR = psi; } 
  /** Should be called at the end of an event */ 
  virtual void End();
  /** Add a contribution @f$ \cos(n(\varphi-\Psi_R))@f$ where   
      @f$ \Psi_R@f$ is the previously set, well-known event plane
      angle. 
      @param phi @f$ \varphi@f$ */ 
  void AddToHarmonic(Double_t phi, Double_t) { fHarmonic.Add(phi, fPsiR); }
  /** Get the value in this bin 
      @param t  Which type of correction
      @return the value of the harmonic */
  virtual Double_t Value(CorType t=naive) const;
  /** Get the value in this bin 
      @param e2 On return, the square error. 
      @param t  Which type of correction
      @return the value of the harmonic */
  Double_t Value(Double_t& e2, CorType t=naive) const;
  /** Get the value in this bin 
      @param e2 On return, the square error. 
      @param t  Which type  of correction
      @return the value of the harmonic */
  Double_t Correction(Double_t& e2, CorType t=naive) const;
  /** Print to standard out. */ 
  void Print(Option_t* option="") const;
  /** The well-known event plane */ 
  Double_t fPsiR; 
  /** True resolution */ 
  AliFMDFlowStat fResReal;
  /** define for ROOT I/O */
  ClassDef(AliFMDFlowTrueBin,1);
  
}; 
/** @brief A "histogram" of objects of class TrueBin in case the
    event plane angle @f$ \Psi_R@f$ is well known. 
    @ingroup x_true */
struct AliFMDFlowTrue1D : public AliFMDFlowBinned1D
{
  /** Constructor */ 
  AliFMDFlowTrue1D(UShort_t order, const AliFMDFlowAxis& xaxis);
  /** Set the well-known event plane angle @f$ \Psi_R@f$ for this
      event. 
      @param psi @f$ \Psi_R@f$ */ 
  void SetPsi(Double_t psi);
  /** Print to standard out */ 
  virtual void Print(Option_t* option="") const;
  /** define for ROOT I/O */
  ClassDef(AliFMDFlowTrue1D,1);
};


#endif
//
// EOF
// 
