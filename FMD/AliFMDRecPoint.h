#ifndef ALIFMDRECPOINT_H
#define ALIFMDRECPOINT_H

// Reconstracted Particles Class: has number of reconstructed
// particles in sectors from NumOfMinSector to NumberOfMaxSector()
// rings from NumOfMinRing to NumOfMaxRing for each FMDvolume 
//
#ifndef ROOT_TObject
# include <TObject.h>
#endif

/** Reconstructed FMD points.  It contains the pseudo-inclusive
    multiplicity 
    @ingroup FMD_rec
 */
class AliFMDRecPoint: public TObject
{
public:
  /** CTOR */
  AliFMDRecPoint();
  /** Constrctor 
      @param detector Detector 
      @param ring     Ring
      @param sector   Sector
      @param strip    Strip 
      @param eta      Psuedo-rapidity @f$ \eta@f$ 
      @param phi      Azimuthal angle @f$ \varphi@f$ 
      @param edep     Energy deposited 
      @param particles Psuedo-inclusive multiplicity */
  AliFMDRecPoint(UShort_t detector,  Char_t   ring, 
		 UShort_t sector,    UShort_t strip, 
		 Float_t  eta,       Float_t  phi,
		 Float_t  edep,      Float_t  particles);
  /** DTOR */
  virtual ~AliFMDRecPoint() {};

  /** @return Detector # */
  UShort_t     Detector()	   const { return fDetector; }
  /** @return Ring ID */
  Char_t       Ring()	           const { return fRing;     }
  /** @return sector # */
  UShort_t     Sector()	           const { return fSector;   }
  /** @return strip # */
  UShort_t     Strip()	           const { return fStrip;    }
  /** @return Psuedo-rapidity @f$ \eta@f$ */
  Float_t      Eta() const             { return fEta; }
  /** @return phi      Azimuthal angle @f$ \varphi@f$ */
  Float_t      Phi() const             { return fPhi; }
  /** @return edep     Energy deposited */
  Float_t      Edep() const            { return fEdep; }
  /** @return particles Psuedo-inclusive multiplicity */
  Float_t      Particles() const       { return fParticles; }
  /** Print information 
      @param opt Not used */
  virtual void Print(Option_t* opt="D") const;
  /** @return Name */
  const char*  GetName()                const;
  /** @return Title */
  const char*  GetTitle()               const;
protected:
  UShort_t fDetector;        // Detector #
  Char_t   fRing;            // Ring ID
  UShort_t fSector;          // Sector #
  UShort_t fStrip;           // Strip # 
  Float_t  fEta;             // Eta value 
  Float_t  fPhi;             // Phi value
  Float_t  fEdep;            // Energy deposited 
  Float_t  fParticles;       // Quasi-number of particles 

  ClassDef(AliFMDRecPoint,1)     // Base class for multiplicity data
};
#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//
