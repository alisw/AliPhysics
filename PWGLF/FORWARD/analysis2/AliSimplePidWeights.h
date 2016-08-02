#ifndef ALIROOT_ALISIMPLEPIDWEIGHTS
#define ALIROOT_ALISIMPLEPIDWEIGHTS
/**
 * @file   AliSimplePidWeights.h
 * @author Marek Chojnacki <mchojnac#cern.ch>
 * @date   Wed Feb  4 00:30:56 2015
 * 
 * @brief  
 * 
 * 
 */
#include "AliBaseMCWeights.h"
#include <TArrayI.h>
#include <TArrayD.h>

/**
 *  Class for weights in calculation of the secondary maps using
 *  simple PID weighting.
 *
 * @ingroup pwglf_forward_mc
 */

class AliSimplePidWeights : public AliBaseMCWeights
{
public:
  /** 
   * Default CTOR 
   */
  AliSimplePidWeights();
  /** 
   * Copy Constructor 
   * 
   * @param o Object to copy from 
   */
  AliSimplePidWeights(const AliSimplePidWeights& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this. 
   */
  AliSimplePidWeights& operator=(const AliSimplePidWeights& o);
  /** 
   * Destructor 
   */
  virtual ~AliSimplePidWeights() {}
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
   * Add a PDG code to be weighted
   * 
   * @param pdg    Particle code 
   * @param weight The factor 
   * @param anti   If true, also add for anti-particle 
   */
  virtual void AddPDGCode(Int_t  pdg, Double_t weight, Bool_t anti=false);
  /** 
   * Print information 
   * 
   * @param option Not used 
   */
  virtual void Print(Option_t* option="") const;
private:
  TArrayI  fPdgs; // list of the particles to be weighted 	
  TArrayD  fWeights; // the weights 
  ClassDef(AliSimplePidWeights,1);
};

#endif
// Local Variables:
//   mode: C++
// End:
