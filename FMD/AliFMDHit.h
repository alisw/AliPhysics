#ifndef AliFMDhit_H
#define AliFMDhit_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
//___________________________________________________________________
//
// AliFMDhit is the hit class for the FMD. Hits are the information
// that comes from a Monte Carlo at each step as a particle mass
// through sensitive detector elements as particles are transported
// through a detector.
//
#ifndef ALIHIT_H
# include <AliHit.h>
#endif

//___________________________________________________________________
class AliFMDHit : public AliHit 
{
public:
  AliFMDHit();
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
  virtual ~AliFMDHit() {}

  UShort_t Detector()	const { return fDetector; }
  Char_t   Ring()	const { return fRing;     }
  UShort_t Sector()	const { return fSector;   }
  UShort_t Strip()	const { return fStrip;    }
  Float_t  Edep()       const { return fEdep;     }
  Float_t  Px()         const { return fPx;       }
  Float_t  Py()         const { return fPy;       }
  Float_t  Pz()         const { return fPz;       } 
  Float_t  P()          const;
  Float_t  M()          const;
  Float_t  Q()          const;
  Int_t    Pdg()        const { return fPdg;      }
  Float_t  Time()       const { return fTime;     }
  Float_t  Length()     const { return fLength;   }
  Bool_t   IsStop()     const { return fStop;     }
  void     Print(Option_t* opt="") const;

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
