#ifndef ALIITSHIT_H
#define ALIITSHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

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
#include <TLorentzVector.h>

#include "AliHit.h" 

class TParticle;

class AliITShit : public AliHit {

 public:
    AliITShit();// Default consrtructor
    // Old Standard Constructor
    AliITShit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
    // New Standard Constructor
    AliITShit(Int_t shunt,Int_t track,Int_t *vol,Float_t edep,Float_t tof,
              TLorentzVector &x,TLorentzVector &x0,TLorentzVector &p);
    // Copy Constructor
    AliITShit(const AliITShit &h);
    // Equals operator
    AliITShit& operator=(const AliITShit &h);
    // Default destructor
    virtual ~AliITShit() {};
    // Get Hit information functions.
    // virtual int GetTrack() const {return fTrack;} // define in AliHit
    // virtual void SetTrack(int track) const {fTrack=track;) // AliHit
    virtual void SetModule(Int_t mod){fModule=mod;};
    virtual void SetShunt(Int_t shunt);
    virtual void SetPosition(TLorentzVector &x){fX=x.X();fY=x.Y();fZ=x.Z();}
    virtual void SetStartPosition(TLorentzVector &x){fx0=x.X();fy0=x.Y();
                                                     fz0=x.Z();}
    virtual void SetTime(Float_t t){fTof = t;}
    virtual void SetStartTime(Float_t t){ft0 = t;}
    virtual void SetStatus(Int_t stat){fStatus = stat;}
    virtual void SetStartStatus(Int_t stat){fStatus0 = stat;}
    virtual void SetEdep(Float_t de){fDestep = de;}
    virtual void SetMomentum(TLorentzVector &p){fPx=p.Px();fPy=p.Py();
                                                fPz=p.Pz();}
    virtual Int_t GetTrackStatus() const {//returns the status code
        return fStatus;}
    virtual Int_t GetTrackStatus0() const {//returns the status code
        return fStatus0;}
    virtual Int_t GetLayer() const{Int_t a,b,c;// returns the layer number
    GetDetectorID(a,b,c); return a;}
    virtual Int_t GetLadder() const{Int_t a,b,c;// returns the ladder number
    GetDetectorID(a,b,c); return b;}
    virtual Int_t GetDetector() const{Int_t a,b,c;//returns the detector number
    GetDetectorID(a,b,c); return c;}
	// returns the detector ID
    virtual void  GetDetectorID(Int_t &layer,Int_t &ladder,Int_t &det)const;
    virtual Int_t GetModule(){return fModule;};
    virtual Float_t GetIonization() const {return fDestep;}//returns the Destep
    //
    virtual void GetPositionG(Float_t &x,Float_t &y,Float_t &z)const {
        // returns the position in the Global frame
        x=fX;y=fY;z=fZ;return;};
    virtual void GetPositionG(Double_t &x,Double_t &y,Double_t &z)const {
        // returns the position in the Global frame
        x=fX;y=fY;z=fZ;return;};
    virtual Float_t GetTOF() const {// returns the time of flight
        return fTof;}
    // Returns particle 3 position at this hit in global coordinates.
    virtual void GetPositionG(Float_t &x,Float_t &y,Float_t &z,Float_t &tof)
        const {// returns the position in the Global frame and the time of
        // flight
        x=fX;y=fY;z=fZ,tof=fTof;return;};
    virtual void GetPositionG(Double_t &x,Double_t &y,Double_t &z,Double_t &t)
        const {// Returns particle 3 position and the time of flight at this
        // hit in global coordinates.
        x=fX;y=fY;z=fZ,t=fTof;return;};
    virtual Float_t GetXG()const {// Returns particle X position at this hit
        // in global coordinates.
        return fX;}
    virtual Float_t GetYG()const {// Returns particle X position at this hit
        // in global coordinates.
        return fY;}
    virtual Float_t GetZG()const {// Returns particle Z position at this hit
        // in global coordinates.
        return fZ;}
    // Returns particle 3 position at this hit in global coordinates.
    virtual void GetPositionG0(Float_t &x,Float_t &y,Float_t &z,Float_t &tof)
        const {// returns the initial position in the Global frame and the
        // time of flight
        x=fx0;y=fy0;z=fz0,tof=fTof;return;};
    // Returns particle 3 position at this hit in global coordinates.
    virtual void GetPositionG0(Double_t &x,Double_t &y,Double_t &z,
                               Double_t &tof)const {
        // returns the initial position in the Global frame and the
        // time of flight
        x=fx0;y=fy0;z=fz0,tof=fTof;return;};
    virtual void GetPositionL(Float_t &x,Float_t &y,Float_t &z){
        // Returns particle 3 position at this hit in local coordinates.
        Float_t tf;GetPositionL(x,y,z,tf);}
    virtual void GetPositionL(Float_t &x,Float_t &y,Float_t &z,Float_t &tof);
    virtual void GetPositionL(Double_t &x,Double_t &y,Double_t &z){
        // Returns particle 3 position at this hit in local coordinates.
        Float_t xf,yf,zf,tf;GetPositionL(xf,yf,zf,tf);x=xf,y=yf;z=zf;}
    virtual void GetPositionL(Double_t &x,Double_t &y,Double_t &z,Double_t &t){
        // Returns particle 3 position and the time of flight at this hit
        // in local coordinates.
        Float_t xf,yf,zf,tf;GetPositionL(xf,yf,zf,tf);x=xf,y=yf;z=zf;t=tf;}
    // Returns particle 3 initial position and the time of flight at this hit
    // in local coordinates.
    virtual void GetPositionL0(Double_t &x,Double_t &y,Double_t &z,
                               Double_t &t);
    virtual Float_t GetXL(){
        // Returns particle X position at this hit in local coordinates.
        Float_t xf,yf,zf,tf;GetPositionL(xf,yf,zf,tf);return xf;}
    virtual Float_t GetYL(){
        // Returns particle Y position at this hit in local coordinates.
        Float_t xf,yf,zf,tf;GetPositionL(xf,yf,zf,tf);return yf;}
    virtual Float_t GetZL(){
        // Returns particle Z position at this hit in local coordinates.
        Float_t xf,yf,zf,tf;GetPositionL(xf,yf,zf,tf);return zf;}
  
    // Get Monti Carlo information about hit.
    virtual void GetMomentumG(Float_t &px,Float_t &py,Float_t &pz)const {
        // returns the particle momentum in the Global frame
        px=fPx;py=fPy;pz=fPz;return;};
    virtual void GetMomentumG(Double_t &px,Double_t &py,Double_t &pz)const {
        // returns the particle momentum in the Global frame
        px=fPx;py=fPy;pz=fPz;return;};
    virtual Float_t GetPXG()const {// Returns particle X momentum at this hit
        // in global coordinates.
        return fPx;}
    virtual Float_t GetPYG()const {// Returns particle Y momentum at
        // this hit in global coordinates.
        return fPy;}
    virtual Float_t GetPZG()const {// Returns particle Z momentum at
        // this hit in global coordinates.
        return fPz;}
    // Returns particle Z momentum at this hit in local coordinates.
    virtual TParticle * GetParticle() const; // Returns ptr to this particle.
    Bool_t StatusInside() const {// checks if the particle is "inside"
        if((fStatus&0x0001)==0) return kFALSE;else return kTRUE;}
    Bool_t StatusEntering() const {// checks if the particle is "entering"
        if((fStatus&0x0002)==0) return kFALSE;else return kTRUE;}
    Bool_t StatusExiting() const {// checks if the particle is "exiting"
        if((fStatus&0x0004)==0) return kFALSE;else return kTRUE;}
    Bool_t StatusOut() const {// checks if the particle is "out"
        if((fStatus&0x0008)==0) return kFALSE;else return kTRUE;}
    Bool_t StatusDisappeared() const {// checks if the part. is "disappeared"
        if((fStatus&0x00010)==0) return kFALSE;else return kTRUE;}
    Bool_t StatusStop() const {// checks if the particle is "stopped"
        if((fStatus&0x00020)==0) return kFALSE;else return kTRUE;}
    Bool_t StatusAlive() const {// checks if the particle is "alive"
        if((fStatus&0x00030)==0) return kFALSE;else return kTRUE;}
    // Prints out the content of this class in ASCII format.
    void Print(ostream *os) const; 
    // Reads in the content of this class in the format of Print
    void Read(istream *is);
    virtual void Print(Option_t *option="") const {TObject::Print(option);}
    virtual Int_t Read(const char *name) {return TObject::Read(name);}


 protected:
    Int_t     fStatus; // Track Status
    Int_t     fModule; // Module number 
    Float_t   fPx;     // PX of particle at the point of the hit
    Float_t   fPy;     // PY of particle at the point of the hit
    Float_t   fPz;     // PZ of particle at the point of the hit
    Float_t   fDestep; // Energy deposited in the current step
    Float_t   fTof;    // Time of flight at the point of the hit
    Int_t     fStatus0;// Track Status of Starting point
    Float_t   fx0;     // Starting point of this step
    Float_t   fy0;     // Starting point of this step
    Float_t   fz0;     // Starting point of this step
    Float_t   ft0;     // Starting point of this step

    ClassDef(AliITShit,3)  //Hits object for set:ITS
	 
}; 
// Input and output function for standard C++ input/output.
ostream& operator<<(ostream &os,AliITShit &source);
istream& operator>>(istream &os,AliITShit &source);

#endif
