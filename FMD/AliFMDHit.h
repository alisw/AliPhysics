// -*- mode: c++ -*-
#ifndef AliFMDhit_H
#define AliFMDhit_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
////////////////////////////////////////////////
//
//  Manager and hits classes for set:FMD     
//
////////////////////////////////////////////////
#ifndef ALIHIT_H
# include "AliHit.h"
#endif
 

///////////////////////////////////////////////////////////////////////
// AliFMDhit is the hit class for the FMD. Hits are the information
// that comes from a Monte Carlo at each step as a particle mass through
// sensitive detector elements as particles are transported through a
// detector.
//
// Data members:
//
// Int_t fTrack
//     See AliHit for a full description. The track number of the track
// that made this hit.
//
// Float_t fX
//     See AliHit for a full description. The global x position of the
// hit (in the standard units of the Monte Carlo).
//
// Float_t fY
//     See AliHit for a full description. The global y position of the
// hit (in the standard units of the Monte Carlo).
//
// Float_t fZ
//     See AliHit for a full description. The global z position of the
// hit (in the standard units of the Monte Carlo).
//
// Int_t fStatus
//     The track status flag. This flag indicates the track status
// at the time of creating this hit. It is made up of the following 8
// status bits from highest order to lowest order bits
// 0           :  IsTrackAlive():    IsTrackStop():IsTrackDisappeared():
// IsTrackOut():IsTrackExiting():IsTrackEntering():IsTrackInside()     .
// See AliMC for a description of these functions. If the function is
// true then the bit is set to one, otherwise it is zero.

// Int_t fVolume
//     The number of the FMD detector that contains this hit. 

// Float_t fEdep
//     The energy lost by the particle during the step ending in this
// hit. The units are those determined by the Monte Carlo.
//
// Float_t fPx
//     The x momentum, in global coordinates, of the particle that
// "created" the hit at the time and position of the hit. The units
// are those determined by the Monte Carlo.
//
// Float_t fPy
//     The y momentum, in global coordinates, of the particle that
// "created" the hit at the time and position of the hit. The units
// are those determined by the Monte Carlo.
//
// Float_t fPz
//     The z momentum, in global coordinates, of the particle that
// "created" the hit at the time and position of the hit. The units
// are those determined by the Monte Carlo.
//
///
// Float_t fTime
//     The time of flight associated with the particle ending in this
// hit. The time is typically measured from the point of creation of the
// original particle (if this particle is a daughter).  The units
// are those determined by the Monte Carlo. 

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
	    Float_t  t=0);
  virtual ~AliFMDHit() {}

  UShort_t Detector()	const { return fDetector; }
  Char_t   Ring()	const { return fRing;     }
  UShort_t Sector()	const { return fSector;   }
  UShort_t Strip()	const { return fStrip;    }
  Float_t  Edep()       const { return fEdep;     }
  Float_t  Px()         const { return fPx;       }
  Float_t  Py()         const { return fPy;       }
  Float_t  Pz()         const { return fPz;       } 
  Int_t    Pdg()        const { return fPdg;      }
  Float_t  Time()       const { return fTime;     }
  void     Print(Option_t* opt="") const;

  void     SetEdep(Float_t edep) { fEdep = edep; }
private:
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

  ClassDef(AliFMDHit,1)  //Hits for detector FMD
};
#endif
//____________________________________________________________________
//
// EOF
//
