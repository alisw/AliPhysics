#ifndef AliFMDhit_H
#define AliFMDhit_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
////////////////////////////////////////////////
//  Manager and hits classes for set:FMD     //
////////////////////////////////////////////////
#include "AliHit.h"
 
 
class AliFMDhit : public AliHit {
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

private:
  Int_t      fVolume;       //Volume copy identifier
  Int_t      fParticle;     //Particle identificator
  Float_t    fEdep;         //Energy deposition
  Float_t    fPx;            // Particle's momentum X
  Float_t    fPy;            // Particle's momentum Y
  Float_t    fPz;            // Particle's momentum Z
  Float_t    fTime;         // Particle's time of flight
  
public:
  AliFMDhit() {}
  AliFMDhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  virtual ~AliFMDhit() {}
  Int_t Volume();
  Float_t Particle();
  Float_t Edep();
  Float_t Px();
  Float_t Py();
  Float_t Pz();
  Float_t Time();

  ClassDef(AliFMDhit,1)  //Hits for detector FMD
};
inline Int_t AliFMDhit::Volume(){return fVolume;} 
inline Float_t AliFMDhit::Particle(){return fParticle;} 
inline Float_t AliFMDhit::Edep(){return fEdep;} 
inline Float_t AliFMDhit::Px(){return fPx;} 
inline Float_t AliFMDhit::Py(){return fPy;} 
inline Float_t AliFMDhit::Pz(){return fPz;} 
inline Float_t AliFMDhit::Time(){return fTime;} 
#endif
