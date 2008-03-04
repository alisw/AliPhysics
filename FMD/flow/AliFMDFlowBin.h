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
    @brief Declaration of a Bin in a Flow "histogram" */
//____________________________________________________________________
//
// This contains an of class AliFMDFlowHarmonic and an object of
// class AliFMDFlowEventPlane to calculate v_n and \Psi_k.  It contain
// two objects of class AliFMDFlowEventPlane to calculate the
// sub-event event planes Psi_A, \Psi_B.  It also contain 3 objects of
// class AliFMDFlowResolution to calculate the event plane angle
// resolution. 
//
#ifndef ALIFMDFLOWBIN_H
#define ALIFMDFLOWBIN_H
#include <flow/AliFMDFlowEventPlane.h>
#include <flow/AliFMDFlowHarmonic.h>
#include <flow/AliFMDFlowResolution.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TObject.h>

//Forward declaration 
class TBrowser;

/** @defgroup c_binned Binned flow 
    @brief This group contains code for binned flow analysis.  Two
    kinds of "histograms" are defined - a 1 dimensional and a 2
    dimensional set of binned objects of class AliFMDFlowBin.   

    Objects of class AliFMDFlowBin contains all the code needed to compute
    flow in a given bin.   

    The class AliFMDFlowAxis encodes look-up of a object of class
    AliFMDFlowBin in a flow "Histogram"
*/
//______________________________________________________
/** @class AliFMDFlowBin flow/AliFMDFlowBin.h <flow/AliFMDFlowBin.h>
    @brief A bin of flow.   

    This contains an of class AliFMDFlowHarmonic and an object of
    class AliFMDFlowEventPlane to calculate @f$ v_n@f$ and
    @f$\Psi_k@f$.  It contain two objects of class
    AliFMDFlowEventPlane to calculate the sub-event event planes
    @f$\Psi_A, \Psi_B@f$.  It also contain 3 objects of class
    AliFMDFlowResolution to calculate the event plane angle
    resolution.

    @ingroup c_binned 
*/
class AliFMDFlowBin : public TObject
{
public:
  /** Correction type */
  enum CorType {
    /** No correction */
    kNone = 1, 
    /** Naive, using the formulas in Voloshins paper */
    kNaive = 2,
    /** STARs way */
    kStar = 4, 
    /** The way used in the TDR */
    kTdr = 8
  };
  /** Types of things to draw */
  enum DrawVar { 
    /** The harmonic */
    kHarmonic, 
    /** The resolution */
    kResolution, 
    /** The statistics */ 
    kCounts
  };
  /** Default Ctor - do not use */ 
  AliFMDFlowBin() 
    : fPsi(),
      fPsiA(),
      fPsiB(),
      fRes(),
      fResStar(),
      fResTdr(),
      fHarmonic(),
      fSplit(),
      fPhi(),
      fNA(0),
      fNB(0),
      fN(0),
      fAB()
  {}
  /** Constructor 
      @param order Order of harmonic. 
      @param k     Factor of event plane order=k * m */
  AliFMDFlowBin(UShort_t order, UShort_t k=1); 
  /** Copy constructor 
      @param o Object top copy from */ 
  AliFMDFlowBin(const AliFMDFlowBin& o);
  /** Assignment operator
      @param o Object to assign from 
      @return Reference to this object */
  AliFMDFlowBin& operator=(const AliFMDFlowBin& o);
  
  /** Destructor */
  virtual ~AliFMDFlowBin() {} 
  
  /** @{ 
      @name Processing */
  /** Should be called at the start of an event */ 
  virtual void Begin();
  /** Called to add a contribution to the event plane 
      @param phi The angle @f$ \varphi \in[0,2\pi]@f$ 
      @param w   Weight
      @param a   If true, add to sub-event A, otherwise to sub-event
      B. */
  virtual void AddToEventPlane(Double_t phi, Double_t w=1, Bool_t a=kTRUE);
  /** Called to add a contribution to the harmonic. 
      @param phi The angle @f$ \varphi \in[0,2\pi]@f$
      @param wp  Weight of @a phi (only used in the calculation of
                 the event plane). 
      @param wh  Weight if the @f$ \varphi@f$ observation */
  virtual void AddToHarmonic(Double_t phi, Double_t wp=1, Double_t wh=1);
  /** Should be called at the end of an event */ 
  virtual void End();
  /** Analyse events 
      @param n    Size of @a phis and possibly @a ws
      @param phis @f$ (\varphi_i, \ldots, \varphi_n)@f$ 
      @param wp   Weights for event plane (optional)
      @param wh   Weights for harmonic (optional - can be the same as @a ws) */
  virtual void Event(UInt_t n, Double_t* phis, Double_t* wp=0, Double_t* wh=0);
  /** Finish run */
  virtual void Finish();
  /** @} */

  /** @{ 
      @name Orders */
  /** Get the harmonic order 
      @return The harmonic order */ 
  virtual UShort_t Order() const { return fHarmonic.Order(); }
  /** Get the harmonic order 
      @return The harmonic order */ 
  virtual UShort_t PsiOrder() const { return fPsi.Order(); }
  /** @} */

  /** @{
      @name Values */
  /** Get the value in this bin 
      @param t  Which type of correction
      @return the value of the harmonic */
  virtual Double_t Value(CorType t=kTdr) const;
  /** Get the value in this bin 
      @param t  Which type of correction 
      @return the error on the value of the harmonic */
  virtual Double_t EValue(CorType t=kTdr) const;
  /** Get the value in this bin 
      @param e2 On return, the square error. 
      @param t  Which type of correction
      @return the value of the harmonic */
  virtual Double_t Value(Double_t& e2, CorType t=kTdr) const;
  /** Get the value in this bin 
      @param e2 On return, the square error. 
      @param t  Which type  of correction
      @return the value of the event plane correction */
  virtual Double_t Correction(Double_t& e2, CorType t=kTdr) const;
  /** Get the number of counts used in this. 
      @return The number of counts used in this bin */
  virtual ULong_t  Counts() const;
  /** @} */

  /** @{ 
      @name Utility */
  /** Print summary to standard output */ 
  virtual void Print(Option_t* option="") const; //*MENU*
  /** Return true */ 
  virtual Bool_t IsFolder() const { return kTRUE; } 
  /** Browse this item */ 
  virtual void Browse(TBrowser* b); 
  /** Get the event plane angle */
  virtual Double_t Psi() const { return fPsi.Psi(); } 
  /** Get the sub-event A plane angle */
  virtual Double_t PsiA() const { return fPsiA.Psi(); } 
  /** Get the sub-event B plane angle */
  virtual Double_t PsiB() const { return fPsiB.Psi(); } 
  /** @} */

  /** @{ 
      @name histograms */ 
  /** @return Split histogram */
  const TH2& SplitHistogram() const { return fSplit; }
  /** @return Phi histogram */
  const TH1& PhiHistogram() const { return fPhi; }
  /** @} */
protected:
  /** Major event plane */
  AliFMDFlowEventPlane fPsi; // Major event plane
  /** Sub-event A event plane */
  AliFMDFlowEventPlane fPsiA; // Sub-event A event plane
  /** Sub-event B event plane */
  AliFMDFlowEventPlane fPsiB; // Sub-event B event plane
  /** Resolution */
  AliFMDFlowResolution fRes; // Resolution
  /** Resolution */
  AliFMDFlowResolutionStar fResStar; // Resolution
  /** Resolution */
  AliFMDFlowResolutionTDR fResTdr; // Resolution
  /** The harmonic */
  AliFMDFlowHarmonic fHarmonic; // Harmonic
  /** Histogram of the relative split between A and B */ 
  TH2D fSplit;
  /** Histogram of relative phi's */ 
  TH1D fPhi;
  /** Counter of number of observations in sub-A */
  ULong_t fNA;
  /** Counter of number of observations in sub-B */
  ULong_t fNB;
  /** Counter of number of observations */
  ULong_t fN;
  /** Psi from A vs from B */
  TH2D    fAB;

  /** Define for ROOT I/O */
  ClassDef(AliFMDFlowBin,1); // A flow analysis 
};


#endif
//
// EOF
//
