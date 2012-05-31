#ifndef ALIFMDESDREVERTEXER_H
#define ALIFMDESDREVERTEXER_H
#include <TObject.h>
class AliESDFMD;

//
// Class to recaculate quantities in an AliESDFMD object based on new
// a value for the z-coordinate of the primary interaction point. 
//
// This allows us, in case we get a better determination of the
// z-coordinate of the primary interaction, to recalibrate the signals
// in the FMD ESD object, without having to redo the reconstruction. 
//
class AliFMDESDRevertexer : public TObject
{
public:
  /** 
   * Constructor
   * 
   */
  AliFMDESDRevertexer();
  /** 
   * Destructor 
   * 
   */
  virtual ~AliFMDESDRevertexer() {}
  /** 
   * Revertex the passed ESD.   The passed ESD object will be modified
   * directly. 
   * 
   * @param fmdEsd ESD object to revertex. 
   * @param vz     New Z coordinate of primary vertex. 
   * 
   * @return @c true on success, @c false failure.
   */
  Bool_t Revertex(AliESDFMD* fmdEsd, Double_t vz) const;
  
  /** 
   * Calculate the physical coordinates (@a eta, @a phi) corresponding
   * to the detector coordinates (@a det, @a rng, @a sec, @a str).
   * 
   * @param det   The detector identifier 
   * @param rng   The ring identifier 
   * @param sec   The sector identifier 
   * @param str   The strip identifier 
   * @param vz    The z coordinate of the current primary interation vertex
   * @param eta   On return, the psuedo-rapidity
   * @param phi   On return, the azimuthal angle
   * @param r     On return, the radius
   * @param theta On return, the polar angle
   */
  Bool_t PhysicalCoordinates(UShort_t  det, 
			     Char_t    rng, 
			     UShort_t  sec, 
			     UShort_t  str,
			     Double_t  vz,
			     Double_t& eta, 
			     Double_t& phi,
			     Double_t& r,
			     Double_t& theta) const;

  /** 
   * Calculate the polar angle @f$ \theta@f$ corresponding to the
   * psuedo-rapidity @f$ \eta@f$ 
   * 
   * @param eta Psuedo rapidity @f$ \eta=-\log[\tan(\theta/2)]@f$ 
   * 
   * @return Polar angle @f$ \theta=2\tan^{-1}[\exp(-\eta)]@f$
   */  
  Double_t Eta2Theta(Double_t eta) const;
protected:
  ClassDef(AliFMDESDRevertexer,0) // Revertex and FMD ESD Object.
};

#endif
//
// Local Variables:
//  mode: C++
// End:
//



