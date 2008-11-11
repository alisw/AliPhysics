#ifndef ALIFMDESDREVERTEXER_H
#define ALIFMDESDREVERTEXER_H
# include <Rtypes.h>
class ::AliESDFMD;

//
// Class to recaculate quantities in an AliESDFMD object based on new
// a value for the z-coordinate of the primary interaction point. 
//
// This allows us, in case we get a better determination of the
// z-coordinate of the primary interaction, to recalibrate the signals
// in the FMD ESD object, without having to redo the reconstruction. 
//
class AliFMDESDRevertexer 
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
   * @return Polar angle @f$ \theta=2\atan[\exp(-\eta)]@f$
   */  
  Double_t Eta2Theta(Double_t eta) const;
protected:
};

# ifndef __CINT__
#  ifndef ALIFMDGEOMETRY_H
#   include <AliFMDGeometry.h>
#  endif
#  ifndef ALIESDFMD_H
#   include <AliESDFMD.h>
#  endif

//____________________________________________________________________
inline 
AliFMDESDRevertexer::AliFMDESDRevertexer()
{
  // Constructor 
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  geom->Init();
  geom->InitTransformations();
}

//____________________________________________________________________
inline Bool_t
AliFMDESDRevertexer::Revertex(AliESDFMD* fmdEsd, Double_t vz) const
{
  // Revertex the passed ESD.   The passed ESD object will be modified
  // directly. 
  // 
  // Parameters:
  //  	fmdEsd ESD object to revertex.
  // 	vz     New Z coordinate of primary vertex. 
  //
  // Return: 
  //    kTRUE on success, kFALSE failure.
  if (!fmdEsd) return kFALSE;
  
  Bool_t         ret  = kTRUE;
  const UShort_t sec0 = 0;
  for (UShort_t det = 1; det <= 3; det++) { 
    UShort_t nrng = (det == 1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nrng; ir++) {
      Char_t   rng  = (ir == 0 ? 'I' : 'O');
      UShort_t nsec = (ir == 0 ?  20 :  40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      for (UShort_t str = 0; str < nstr; str++) { 
	Double_t phi, r, theta;
	Double_t eta      = AliESDFMD::kInvalidEta;
	Double_t oldEta   = fmdEsd->Eta(det, rng, 0u, str);
	Double_t oldTheta = Eta2Theta(oldEta);
	Bool_t   ret1     = PhysicalCoordinates(det, rng, sec0, str, vz, 
						eta, phi, r, theta);
	fmdEsd->SetEta(det, rng, sec0, str, eta);

	if (!ret1) {
	  // If the was an error, then there's no reason to go on with
	  // this strip-ring.  Note, that the eta is correctly set to
	  // AliESDFMD::kInvalidMult. 
	  ret = kFALSE;
	  continue;
	}
	
	Double_t corr = 1; 
	if (fmdEsd->IsAngleCorrected()) 
	  corr = TMath::Abs(TMath::Cos(theta) / TMath::Cos(oldTheta));
	for (UShort_t sec = 0; sec < nsec; sec++) { 
	  Double_t mult = fmdEsd->Multiplicity(det, rng, sec, str);
	  fmdEsd->SetMultiplicity(det, rng, sec, str, corr * mult);
	}
      }
    }
  }

  return ret;
}

//____________________________________________________________________
inline Double_t
AliFMDESDRevertexer::Eta2Theta(Double_t eta) const
{
  // Calculate the polar angle @f$ \theta@f$ corresponding to the
  // psuedo-rapidity @f$ \eta@f$ 
  // 
  // Parameters:
  //	eta Psuedo rapidity @f$ \eta=-\log[\tan(\theta/2)]@f$ 
  // Return:
   //	Polar angle @f$ \theta=2\atan[\exp(-\eta)]@f$
  return 2 * TMath::ATan(TMath::Exp(-eta));
}


//____________________________________________________________________
inline Bool_t
AliFMDESDRevertexer::PhysicalCoordinates(UShort_t  det, 
					 Char_t    rng, 
					 UShort_t  sec, 
					 UShort_t  str,
					 Double_t  vz,
					 Double_t& eta, 
					 Double_t& phi,
					 Double_t& r,
					 Double_t& theta) const
{
  // Calculate the physical coordinates (@a eta, @a phi) corresponding
  // to the detector coordinates (@a det, @a rng, @a sec, @a str).
  // 
  // Parameters:
  //	det  The detector identifier 
  //	rng  The ring identifier 
  //	sec  The sector identifier 
  //	str  The strip identifier 
  //	vz   The z coordinate of the current primary interation vertex
  //	eta  On return, the psuedo-rapidity
  //	phi  On return, the azimuthal angle
  //	r    On return, the radius
  //	phi  On return, the polar angle
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  Double_t x=0, y=0, z=0;
  geom->Detector2XYZ(det, rng, sec, str, x, y, z);

  // Check that the conversion succeeded
  if (x == 0 && y == 0 && z == 0) return kFALSE;
  
  // Correct for vertex offset. 
  z     += vz;
  phi   =  TMath::ATan2(y, x);
  r     =  TMath::Sqrt(y * y + x * x);
  theta =  TMath::ATan2(r, z);
  eta   = -TMath::Log(TMath::Tan(theta / 2));

  return kTRUE;
}
# endif // __CINT__
#endif
//
// Local Variables:
//  mode: C++
// End:
//



