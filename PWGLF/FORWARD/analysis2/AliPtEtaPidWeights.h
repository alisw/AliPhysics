#ifndef ALIROOT_ALIPTETAPIDWEIGHTS
#define ALIROOT_ALIPTETAPIDWEIGHTS
/**
 * @file   AliPtEtaPidWeights.h
 * @author Marek Chojnacki <mchojnac#cern.ch>
 * @date   Wed Feb  4 00:30:56 2015
 * 
 * @brief  
 * 
 * 
 */
#include "AliBaseMCWeights.h"
#include <TArrayI.h>
#include <TList.h>
class TH2;

/**
 *  Class for weights in calculation of the secondary maps using @f$
 *  p_T, \eta@f$ PID weighting.
 *
 * @ingroup pwglf_forward_mc
 */

class AliPtEtaPidWeights : public AliBaseMCWeights
{
public:
  /** 
   * Default CTOR 
   */
  AliPtEtaPidWeights();
  /** 
   * Copy Constructor 
   * 
   * @param o Object to copy from 
   */
  AliPtEtaPidWeights(const AliPtEtaPidWeights& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this. 
   */
  AliPtEtaPidWeights& operator=(const AliPtEtaPidWeights& o);
  /** 
   * Destructor 
   */
  virtual ~AliPtEtaPidWeights() {}
  /** 
   * Calculate the weight of a single particle 
   * 
   * @param eta    Pseudo rapidity @f$\eta@f$ 
   * @param pt     Transvers momentum @f$p_{T}@f$ 
   * @param phi    Azimuth angle @f$\varphi@f$ 
   * @param id     Particle ID (PDG code)
   * @param phiR   Event plane angle @f$\Psi_R@f$ 
   * @param b      Event impact parameter @f$b@f$ in fermi-meter 
   * 
   * @return Particle weight 
   */
  virtual Double_t CalcWeight(Double_t eta,
			      Double_t pt,
			      Double_t phi,
			      Int_t    id,
			      Double_t phiR,
			      Double_t b) const;
  /** 
   * Initialize this object 
   * 
   * @param l List of output objects 
   */
  virtual void Init(TList* l);
  /** 
   * Add a PDG code to be weighted.  Note that this object will claim
   * ownership of the passed histogram.  The X axis of the histogram
   * is assumed to be the @f$\eta@f$ axis, while the Y axis is assumed
   * to be @f$p_T@f$ axis.
   * 
   * @param pdg Particle code 
   * @param weight The histogram of (eta,phi) dependent weight factors 
   */
  virtual void AddPDGCode(Int_t  pdg, TH2* weight);
  /** 
   * Print information 
   * 
   * @param option Not used 
   */
  virtual void Print(Option_t* option="") const;
private:
  TArrayI  fPdgs; // list of the particles to be weighted 	
  TList   fWeights; // the weights 
  ClassDef(AliPtEtaPidWeights,1);
};

#endif
// Local Variables:
//   mode: C++
// End:
