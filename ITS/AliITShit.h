#ifndef ITSHIT_H
#define ITSHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliDetector.h"
#include "AliHit.h" 
#include "AliDigit.h"
#include "AliITSgeom.h"


class AliITShit : public AliHit {
////////////////////////////////////////////////////////////////////////
// Version: 0
// Written by Rene Brun, Federico Carminati, and Roberto Barbera
//
// Version: 1
// Modified and documented by Bjorn S. Nilsen
// July 11 1999
//
// AliITShit is the hit class for the ITS. Hits are the information
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
//
// Int_t fLayer
//     The layer number of the detector that contains this hit. See
// AliITSgeom and AliITSv? for a description of the geometry.
//
// Int_t fLadder
//     The ladder number of the detector that contains this hit. See
// AliITSgeom and AliITSv? for a description of the geometry.
//
// Int_t fDet
//     The detector number of the detector that contains this hit. See
// AliITSgeom and AliITSv? for a description of the geometry.
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
// Float_t fDestep
//     The energy lost by the particle during the step ending in this
// hit. The units are those determined by the Monte Carlo.
//
// Float_t fTof
//     The time of flight associated with the particle ending in this
// hit. The time is typically measured from the point of creation of the
// original particle (if this particle is a daughter).  The units
// are those determined by the Monte Carlo.
//
//
// Member functions:
//
// AliITShit()
//     The default creator of the AliITShit class.
//
// AliITShit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits)
//     The creator of the AliITShit class. The variables shunt and
// track are passed to the creator of the AliHit class. See the AliHit
// class for a full description. the integer array *vol contains, in order,
// fLayer = vol[0], fDet = vol[1], fLadder = vol[2], fStatus = vol[3].
// The array *hits contains, in order, fX = hits[0], fY = hits[1], 
// fZ = hits[2], fPx = hits[3], fPy = hits[4], fPz = hits[5],
// fDestep = hits[6], and fTof = hits[7].
//
// ~AliITShit()
//     The default destructor of the AliITShit class.
//
// int GetTrack()
//     See AliHit for a full description. Returns the track number fTrack
// for this hit.
//
// SetTrack(int track)
//     See AliHit for a full description. Sets the track number fTrack
// for this hit.
//
// Int_t GetTrackStatus()
//     Returns the value of the track status flag fStatus. This flag
// indicates the track status at the time of creating this hit. It is
// made up of the following 8 status bits from highest order to lowest
// order bits
// 0           :  IsTrackAlive():    IsTrackStop():IsTrackDisappeared():
// IsTrackOut():IsTrackExiting():IsTrackEntering():IsTrackInside()     .
// See AliMC for a description of these functions. If the function is
// true then the bit is set to one, otherwise it is zero.
//
// Int_t GetLayer()
//     Returns the layer number, fLayer, for this hit.
//
// Int_t GetLadder()
//     Returns the ladder number, fLadder, for this hit.
//
// Int_t GetDetector()
//     Returns the detector number, fDet, for this hit.
//
// GetDetectorID(Int_t &layer, Int_t &ladder, Int_t &detector)
//     Returns the layer, ladder, and detector numbers, fLayer fLadder fDet,
// in one call.
//
// Float_t GetIonization()
//     Returns the energy lost, fDestep, by the particle creating this hit,
// in the units defined by the Monte Carlo.
//
// GetPoositionG(Float_t &x, Float_t &y, Float_t &z)
//     Returns the global position, fX fY fZ, of this hit, in the units
// define by the Monte Carlo.
//
// Float_t GetTOF()
//     Returns the time of flight, fTof, of this hit, in the units defined
// by the Monte Carlo.
//
// GetPositionG(Float_t &x, Float_t &y, Float_t &z, Float_t &tof)
//     Returns the global position and time of flight, fX fY fZ fTof, of
// this hit, in the units define by the Monte Carlo.
//
// GetPositioonP(Float_t &px, Float_t &py, Float_t &pz)
//     Returns the global momentum, fPx fPy fPz, of the particle that made
// this hit, in the units define by the Monte Carlo.
////////////////////////////////////////////////////////////////////////
    // public;       // defined in AliHit
    // Int_t fTrack  // defined in AliHit
    // Float_t fX;   // defined in AliHit
    // Float_t fY;   // defined in AliHit
    // Float_t fZ;   // defined in AliHit

 public:
//private:
    Int_t     fStatus; // Track Status
    Int_t     fLayer;  // Layer number
    Int_t     fLadder; // Ladder number
    Int_t     fDet;    // Detector number  
    Float_t   fPx;     // PX of particle at the point of the hit
    Float_t   fPy;     // PY of particle at the point of the hit
    Float_t   fPz;     // PZ of particle at the point of the hit
    Float_t   fDestep; // Energy deposited in the current step
    Float_t   fTof;    // Time of flight at the point of the hit

 public:
    AliITShit() {}
    AliITShit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
    virtual ~AliITShit() {}
    // Get Hit information functions.
    // inline virtual int GetTrack(){return fTrack;} // define in AliHit
    // inline virtual void SetTrack(int track){fTrack=track;) // AliHit
    inline virtual Int_t GetTrackStatus() {return fStatus;}
    inline virtual Int_t GetLayer() {return fLayer;}
    inline virtual Int_t GetLadder() {return fLadder;}
    inline virtual Int_t GetDetector() {return fDet;}
    inline virtual void  GetDetectorID(Int_t &layer,Int_t &ladder,
	 			       Int_t &detector)
                     {layer=fLayer;ladder=fLadder;detector=fDet;return;};
    inline virtual Float_t GetIonization() {return fDestep;}
    //
    inline virtual void GetPositionG(Float_t &x,Float_t &y,Float_t &z)
                                    {x=fX;y=fY;z=fZ;return;};
    inline virtual Float_t GetTOF() {return fTof;}
    inline virtual void GetPositionG(Float_t &x,Float_t &y,Float_t &z,
				    Float_t &tof)
                                    {x=fX;y=fY;z=fZ,tof=fTof;return;};
    inline virtual Float_t GetXG(){return fX;}
    inline virtual Float_t GetYG(){return fY;}
    inline virtual Float_t GetZG(){return fZ;}
           virtual void GetPositionL(Float_t &x,Float_t &y,Float_t &z);
           virtual void GetPositionL(Float_t &x,Float_t &y,Float_t &z,
				     Float_t &tof);
           virtual Float_t GetXL();
           virtual Float_t GetYL();
           virtual Float_t GetZL();
    // Get Monti Carlo information about hit.
    inline virtual void GetMomentumG(Float_t &px,Float_t &py,Float_t &pz)
                                    {px=fPx;py=fPy;pz=fPz;return;};
           virtual void GetMomentumL(Float_t &px,Float_t &py,Float_t &pz);
    ClassDef(AliITShit,1)  //Hits object for set:ITS
};

#endif
