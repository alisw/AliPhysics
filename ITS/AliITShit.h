#ifndef ALIITSHIT_H
#define ALIITSHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliDetector.h"
#include "TParticle.h"
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
////////////////////////////////////////////////////////////////////////
    // public;       // defined in AliHit
    // Int_t fTrack  // defined in AliHit
    // Float_t fX;   // defined in AliHit
    // Float_t fY;   // defined in AliHit
    // Float_t fZ;   // defined in AliHit

 public:
//protected:
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
    // virtual int GetTrack() const {return fTrack;} // define in AliHit
    // virtual void SetTrack(int track) const {fTrack=track;) // AliHit
    virtual Int_t GetTrackStatus() const {return fStatus;}
    virtual Int_t GetLayer() const {return fLayer;}
    virtual Int_t GetLadder() const {return fLadder;}
    virtual Int_t GetDetector() const {return fDet;}
    virtual void  GetDetectorID(Int_t &layer,Int_t &ladder,
	 			       Int_t &detector)
                     const {layer=fLayer;ladder=fLadder;detector=fDet;return;};
    virtual Int_t GetModule();
    virtual Float_t GetIonization() const {return fDestep;}
    //
    virtual void GetPositionG(Float_t &x,Float_t &y,Float_t &z)
                                    const {x=fX;y=fY;z=fZ;return;};
    virtual void GetPositionG(Double_t &x,Double_t &y,Double_t &z)
                                    const {x=fX;y=fY;z=fZ;return;};
    virtual Float_t GetTOF() const {return fTof;}
    // Returns particle 3 position at this hit in global coordinates.
    virtual void GetPositionG(Float_t &x,Float_t &y,Float_t &z,
				    Float_t &tof)
                                    const {x=fX;y=fY;z=fZ,tof=fTof;return;};
    virtual void GetPositionG(Double_t &x,Double_t &y,Double_t &z,
				    Double_t &tof)
                                    const {x=fX;y=fY;z=fZ,tof=fTof;return;};
    // Returns particle 3 position and the time of flight at this hit
    // in global coordinates.
    virtual Float_t GetXG()const {return fX;}
    // Returns particle X position at this hit in global coordinates.
    virtual Float_t GetYG()const {return fY;}
    // Returns particle Y position at this hit in global coordinates.
    virtual Float_t GetZG()const {return fZ;}
    // Returns particle Z position at this hit in global coordinates.
    virtual void GetPositionL(Float_t &x,Float_t &y,Float_t &z);
    // Returns particle 3 position at this hit in local coordinates.
    virtual void GetPositionL(Float_t &x,Float_t &y,Float_t &z,
				     Float_t &tof);
    virtual void GetPositionL(Double_t &x,Double_t &y,Double_t &z){
	Float_t xf,yf,zf;GetPositionL(xf,yf,zf);x=xf,y=yf;z=zf;}
    // Returns particle 3 position at this hit in local coordinates.
    virtual void GetPositionL(Double_t &x,Double_t &y,Double_t &z,
				     Double_t &tof){
	Float_t xf,yf,zf,tf;GetPositionL(xf,yf,zf,tf);x=xf,y=yf;z=zf;tof=tf;}
    // Returns particle 3 position and the time of flight at this hit
    // in local coordinates.
    virtual Float_t GetXL();
    // Returns particle X position at this hit in local coordinates.
    virtual Float_t GetYL();
    // Returns particle Y position at this hit in local coordinates.
    virtual Float_t GetZL();
    // Returns particle Z position at this hit in local coordinates.
    // Get Monti Carlo information about hit.
    virtual void GetMomentumG(Float_t &px,Float_t &py,Float_t &pz)
                                    const {px=fPx;py=fPy;pz=fPz;return;};
    virtual void GetMomentumG(Double_t &px,Double_t &py,Double_t &pz)
                                    const {px=fPx;py=fPy;pz=fPz;return;};
    // Returns particle 3 momentum at this hit in global coordinates.
    virtual Float_t GetPXG()const {return fPx;}
    // Returns particle X momentum at this hit in global coordinates.
    virtual Float_t GetPYG()const {return fPy;}
    // Returns particle Y momentum at this hit in global coordinates.
    virtual Float_t GetPZG()const {return fPz;}
    // Returns particle Z momentum at this hit in global coordinates.
    virtual void GetMomentumL(Float_t &px,Float_t &py,Float_t &pz);
    virtual void GetMomentumL(Double_t &px,Double_t &py,Double_t &pz){
	Float_t x,y,z;GetMomentumL(x,y,z);px=x,py=y,pz=z;}
    // Returns particle 3 momentum at this hit in local coordinates.
    virtual Float_t GetPXL();
    // Returns particle X momentum at this hit in local coordinates.
    virtual Float_t GetPYL();
    // Returns particle Y momentum at this hit in local coordinates.
    virtual Float_t GetPZL();
    // Returns particle Z momentum at this hit in local coordinates.
    virtual TParticle * GetParticle(); // Returns pointer to this particle.
    Bool_t StatusInside() {if((fStatus&0x0001)==0) return kFALSE;
                                              else return kTRUE;}
    Bool_t StatusEntering() {if((fStatus&0x0002)==0) return kFALSE;
                                                else return kTRUE;}
    Bool_t StatusExiting() {if((fStatus&0x0004)==0) return kFALSE;
                                               else return kTRUE;}
    Bool_t StatusOut() {if((fStatus&0x0008)==0) return kFALSE;
                                           else return kTRUE;}
    Bool_t StatusDisappeared() {if((fStatus&0x00010)==0) return kFALSE;
                                                    else return kTRUE;}
    Bool_t StatusStop() {if((fStatus&0x00020)==0) return kFALSE;
                                             else return kTRUE;}
    Bool_t StatusAlive() {if((fStatus&0x00030)==0) return kFALSE;
                                              else return kTRUE;}

    ClassDef(AliITShit,1)  //Hits object for set:ITS
};

#endif
