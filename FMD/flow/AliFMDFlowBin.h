// -*- mode: C++ -*-
/** @file 
    @brief Declaration of a Bin in a Flow "histogram" */
#ifndef FLOW_BIN_H
#define FLOW_BIN_H
#include <flow/AliFMDFlowEventPlane.h>
#include <flow/AliFMDFlowHarmonic.h>
#include <flow/AliFMDFlowResolution.h>
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
    none, 
    /** Naive, using the formulas in Voloshins paper */
    naive,
    /** STARs way */
    star, 
    /** The way used in the TDR */
    tdr
  };
  /** Constructor */
  AliFMDFlowBin(UShort_t order, UShort_t k=1) 
    : fPsi(order / k), 
      fPsiA(order / k), 
      fPsiB(order / k), 
      fRes(order / k), 
      fResStar(order / k), 
      fResTdr(order / k),
      fHarmonic(order) 
  {}
  /** Destructor */
  virtual ~AliFMDFlowBin() {} 
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
      @param w   Weight of @a phi (only used in the calculation of
      the event plane). */
  virtual void AddToHarmonic(Double_t phi, Double_t w=1);
  /** Should be called at the end of an event */ 
  virtual void End();
  /** Analyse events 
      @param phis @f$ (\varphi_i, \ldots, \varphi_n)@f$ 
      @param ws   Weights (optional)
      @param n    Size of @a phis and possibly @a ws */
  virtual void Event(Double_t* phis, Double_t* ws, UInt_t n);
  /** Finish run */
  virtual void Finish();
  /** Get the value in this bin 
      @param t  Which type of correction
      @return the value of the harmonic */
  virtual Double_t Value(CorType t=naive) const;
  /** Get the value in this bin 
      @param t  Which type of correction 
      @return the error on the value of the harmonic */
  virtual Double_t EValue(CorType t=naive) const;
  /** Get the value in this bin 
      @param e2 On return, the square error. 
      @param t  Which type of correction
      @return the value of the harmonic */
  virtual Double_t Value(Double_t& e2, CorType t=naive) const;
  /** Get the value in this bin 
      @param e2 On return, the square error. 
      @param t  Which type  of correction
      @return the value of the harmonic */
  virtual Double_t Correction(Double_t& e2, CorType t=naive) const;
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

protected:
  /** Major event plane */
  AliFMDFlowEventPlane fPsi;
  /** Sub-event A event plane */
  AliFMDFlowEventPlane fPsiA;
  /** Sub-event B event plane */
  AliFMDFlowEventPlane fPsiB;
  /** Resolution */
  AliFMDFlowResolution fRes;
  /** Resolution */
  AliFMDFlowResolutionStar fResStar;
  /** Resolution */
  AliFMDFlowResolutionTDR fResTdr;
  /** The harmonic */
  AliFMDFlowHarmonic fHarmonic;
  /** Define for ROOT I/O */
  ClassDef(AliFMDFlowBin,1);
};


#endif
//
// EOF
//
