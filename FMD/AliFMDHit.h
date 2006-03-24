#ifndef AliFMDhit_H
#define AliFMDhit_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
//___________________________________________________________________
//
#ifndef ALIHIT_H
# include <AliHit.h>
#endif

//___________________________________________________________________
/** AliFMDhit is the hit class for the FMD. Hits are the information
    that comes from a Monte Carlo at each step as a particle mass
    through sensitive detector elements as particles are transported
    through a detector. 
    @ingroup FMD_sim
*/
class AliFMDHit : public AliHit 
{
public:
  /** Default CTOR */
  AliFMDHit();
  /** Normal FMD hit ctor
      @param shunt     ???
      @param track     Track #
      @param detector  Detector # (1, 2, or 3)                      
      @param ring      Ring ID ('I' or 'O')
      @param sector    Sector # (For inner/outer rings: 0-19/0-39)
      @param strip     Strip # (For inner/outer rings: 0-511/0-255)
      @param x	       Track's X-coordinate at hit
      @param y	       Track's Y-coordinate at hit
      @param z	       Track's Z-coordinate at hit
      @param px	       X-component of track's momentum 
      @param py	       Y-component of track's momentum
      @param pz	       Z-component of track's momentum
      @param edep      Energy deposited by track
      @param pdg       Track's particle Id #
      @param t	       Time when the track hit 
      @param l         Track lenght through medium 
      @param stop      Whether track is stopped in medium */
  AliFMDHit(Int_t    shunt, 
	    Int_t    track, 
	    UShort_t detector, 
	    Char_t   ring, 
	    UShort_t sector, 
	    UShort_t strip, 
	    Float_t  x=0, 
	    Float_t  y=0, 
	    Float_t  z=0,
	    Float_t  px=0, 
	    Float_t  py=0, 
	    Float_t  pz=0,
	    Float_t  edep=0,
	    Int_t    pdg=0,
	    Float_t  t=0, 
	    Float_t  l=0, 
	    Bool_t   stop=kFALSE);
  /** DTOR  */
  virtual ~AliFMDHit() {}
  
  /** @return Detector # */
  UShort_t Detector()	const { return fDetector; }
  /** @return Ring ID */
  Char_t   Ring()	const { return fRing;     }
  /** @return Sector # */
  UShort_t Sector()	const { return fSector;   }
  /** @return Strip # */
  UShort_t Strip()	const { return fStrip;    }
  /** @return Energy deposited (MeV) */
  Float_t  Edep()       const { return fEdep;     }
  /** @return Track @f$ p_x@f$ - momentum in @f$ x@f$ (GeV) */
  Float_t  Px()         const { return fPx;       }
  /** @return Track @f$ p_y@f$ - momentum in @f$ y@f$ (GeV) */
  Float_t  Py()         const { return fPy;       }
  /** @return Track @f$ p_z@f$ - momentum in @f$ z@f$ (GeV) */
  Float_t  Pz()         const { return fPz;       } 
  /** @return Track @f$ |p|@f$ - momentum (GeV) */
  Float_t  P()          const;
  /** @return Track @f$ m@f$ - mass (GeV) */
  Float_t  M()          const;
  /** @return Track @f$ q@f$ - charge (1/3) */
  Float_t  Q()          const;
  /** @return Track PDG id number */
  Int_t    Pdg()        const { return fPdg;      }
  /** @return Time of hit in seconds */
  Float_t  Time()       const { return fTime;     }
  /** @return Path length through silicon */
  Float_t  Length()     const { return fLength;   }
  /** @return Whether track was stopped in silicon */
  Bool_t   IsStop()     const { return fStop;     }

  /** Print info 
      @param opt Not used */
  void        Print(Option_t* opt="") const;
  /** @return Get Name */
  const char* GetName()               const;
  /** @return Get title */
  const char* GetTitle()              const;

  /** Set enenrgy deposited
      @param edep Energy deposited */
  void     SetEdep(Float_t edep) { fEdep = edep; }
protected:
  UShort_t fDetector;  // (Sub) Detector # (1,2, or 3)
  Char_t   fRing;      // Ring ID ('I' or 'O')
  UShort_t fSector;    // Sector # (phi division)
  UShort_t fStrip;     // Strip # (radial division)
  Float_t  fPx;        // Particle's X momentum X
  Float_t  fPy;        // Particle's Y momentum Y
  Float_t  fPz;        // Particle's Z momentum Z
  Int_t    fPdg;       // Particles PDG code 
  Float_t  fEdep;      // Energy deposition
  Float_t  fTime;      // Particle's time of flight
  Float_t  fLength;    // Track length through material. 
  Bool_t   fStop;      // Whether track was stopped. 
  
  ClassDef(AliFMDHit,2)  //Hits for detector FMD
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
