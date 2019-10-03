#ifndef ALIROOT_ALIBASEDMCWEIGHTS
#define ALIROOT_ALIBASEDMCWEIGHTS
/**
 * @file   AliBaseMCWeights.h
 * @author Marek Chojnacki <mchojnac#cern.ch>
 * @date   Wed Feb  4 00:30:56 2015
 * 
 * @brief  
 * 
 * 
 */
#include <TObject.h>
class TList;
class AliMCParticle;

/**
 * Base class for weights in calculation of the secondary maps
 *
 * @ingroup pwglf_forward_mc
 */
class AliBaseMCWeights : public TObject
{
public:
  /** 
   * Default CTOR 
   */
  AliBaseMCWeights() : TObject() {}
  /** 
   * Copy Constructor 
   * 
   * @param o Object to copy from 
   */
  AliBaseMCWeights(const AliBaseMCWeights& o) : TObject(o) {}
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this. 
   */
  AliBaseMCWeights& operator=(const AliBaseMCWeights& o);
  /** 
   * Destructor
   */
  virtual ~AliBaseMCWeights() {}
  /** 
   * Calculate the weight of a given particle.  In the default
   * implementation, this forwards to the same function but with
   * explicit eta, pt, phi, and id arguments.  
   * 
   * However, this is the function that is called by the referencing
   * code, so one can safely overload this and make the other a dummy
   * 
   * @param p         Particle 
   * @param isPrimary True if this refers to a primary particle 
   * @param phiR      Event plane angle @f$\Psi_R@f$ 
   * @param b         Event impact parameter @f$b@f$ in fermi-meter 
   * 
   * @return The particle weight
   */
  virtual Double_t CalcWeight(const AliMCParticle* p,
			      Bool_t               isPrimary,
			      Double_t             phiR,
			      Double_t             b) const;
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
			      Double_t b) const = 0;
  /** 
   * Initialize this object 
   * 
   * @param l Output list 
   */
  virtual void Init(TList* l); 
  /** 
   * Print information to standard out
   * 
   * @param option Not used
   */
  virtual void Print(Option_t* option="") const;
  
  ClassDef(AliBaseMCWeights,1);
};

#endif
// Local Variables:
//   mode: C++
// End:
