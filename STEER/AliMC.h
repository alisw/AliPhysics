#ifndef ALIMC_H
#define ALIMC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
//   Abstract Monte Carlo interface                                          //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include <AliRndm.h>
#include "AliMCProcess.h"

class TLorentzVector;
class AliMC;
class AliDecayer;
class TArrayI;

R__EXTERN AliMC *gMC;

class AliMC : public TNamed, public AliRndm
{
  public:
    AliMC(const char *name, const char *title);
    AliMC() {fRandom=0;}
    virtual ~AliMC() {fgMC=gMC=0;fRandom=0;}
  
    // static access method
    static AliMC* GetMC() { return fgMC; }

    //
    // methods for building/management of geometry
    // ------------------------------------------------
    //

    // functions from GCONS 
    virtual void  Gfmate(Int_t imat, char *name, Float_t &a, Float_t &z,  
  		         Float_t &dens, Float_t &radl, Float_t &absl,
		         Float_t* ubuf, Int_t& nbuf) = 0; 

    // detector composition
    virtual void  Material(Int_t& kmat, const char* name, Float_t a, 
                     Float_t z, Float_t dens, Float_t radl, Float_t absl,
                     Float_t* buf, Int_t nwbuf) = 0;
    virtual void  Mixture(Int_t& kmat, const char *name, Float_t *a, 
                     Float_t *z, Float_t dens, Int_t nlmat, Float_t *wmat) = 0;
    virtual void  Medium(Int_t& kmed, const char *name, Int_t nmat, 
                     Int_t isvol, Int_t ifield, Float_t fieldm, Float_t tmaxfd, 
                     Float_t stemax, Float_t deemax, Float_t epsil, 
		     Float_t stmin, Float_t* ubuf, Int_t nbuf) = 0;
    virtual void  Matrix(Int_t& krot, Float_t thetaX, Float_t phiX, 
                     Float_t thetaY, Float_t phiY, Float_t thetaZ, 
		     Float_t phiZ) = 0;
    virtual void  Gstpar(Int_t itmed, const char *param, Float_t parval) = 0; 

    // functions from GGEOM 
    virtual Int_t Gsvolu(const char *name, const char *shape, Int_t nmed,  
                         Float_t *upar, Int_t np) = 0; 
    virtual void  Gsdvn(const char *name, const char *mother, Int_t ndiv, 
                        Int_t iaxis) = 0; 
    virtual void  Gsdvn2(const char *name, const char *mother, Int_t ndiv, 
                         Int_t iaxis, Float_t c0i, Int_t numed) = 0; 
    virtual void  Gsdvt(const char *name, const char *mother, Float_t step, 
                        Int_t iaxis, Int_t numed, Int_t ndvmx) = 0; 
    virtual void  Gsdvt2(const char *name, const char *mother, Float_t step, 
                         Int_t iaxis, Float_t c0, Int_t numed, Int_t ndvmx) = 0; 
    virtual void  Gsord(const char *name, Int_t iax) = 0; 
    virtual void  Gspos(const char *name, Int_t nr, const char *mother,  
                        Float_t x, Float_t y, Float_t z, Int_t irot, 
                        const char *konly="ONLY") = 0; 
    virtual void  Gsposp(const char *name, Int_t nr, const char *mother,  
                         Float_t x, Float_t y, Float_t z, Int_t irot,
                         const char *konly, Float_t *upar, Int_t np) = 0;

    virtual void  SetCerenkov(Int_t itmed, Int_t npckov, Float_t *ppckov,
                  Float_t *absco, Float_t *effic, Float_t *rindex) = 0;
  
    
    // functions for drawing
    virtual void  DrawOneSpec(const char* name) = 0;
    virtual void  Gsatt(const char* name, const char* att, Int_t val) = 0;
    virtual void  Gdraw(const char*,Float_t theta = 30, Float_t phi = 30,
		        Float_t psi = 0, Float_t u0 = 10, Float_t v0 = 10,
		        Float_t ul = 0.01, Float_t vl = 0.01) = 0;

    // Euclid
    virtual void WriteEuclid(const char*, const char*, Int_t, Int_t) = 0;
		               
    // get methods
    virtual Int_t VolId(const Text_t* volName) const = 0;
    virtual const char* VolName(Int_t id) const = 0;
    virtual Int_t NofVolumes() const = 0;
    virtual Int_t VolId2Mate(Int_t id) const = 0;

    //
    // methods for physics management
    // ------------------------------------------------
    //
 
    // set methods
    virtual void SetCut(const char* cutName, Float_t cutValue) = 0;
    virtual void SetProcess(const char* flagName, Int_t flagValue) = 0;
    virtual Float_t Xsec(char*, Float_t, Int_t, Int_t) = 0; 
 
        // particle table usage         
    virtual Int_t   IdFromPDG(Int_t id) const =0;  
    virtual Int_t   PDGFromId(Int_t pdg) const =0;  
    virtual void    DefineParticles() = 0;      
  
    //
    // methods for step management
    // ------------------------------------------------
    //

    // action methods
    virtual void StopTrack() =0;
    virtual void StopEvent() =0;   

    // set methods
    virtual void SetMaxStep(Float_t) = 0;
    virtual void SetMaxNStep(Int_t) = 0;
    virtual void SetUserDecay(Int_t) =0;  

    // get methods
         // tracking volume(s) 
    virtual Int_t CurrentVolID(Int_t& copyNo) const =0;
    virtual Int_t CurrentVolOffID(Int_t off, Int_t& copyNo) const =0;
    virtual const char* CurrentVolName() const =0;
    virtual const char* CurrentVolOffName(Int_t off) const =0;
    virtual Int_t CurrentMaterial(Float_t &a, Float_t &z, 
                    Float_t &dens, Float_t &radl, Float_t &absl) const =0;  
    virtual Int_t CurrentEvent() const =0; 
    virtual void  Gmtod(Float_t* xm, Float_t* xd, Int_t iflag) = 0;
    virtual void  Gdtom(Float_t* xd, Float_t* xm, Int_t iflag)= 0 ;
    virtual Float_t MaxStep() const =0;
    virtual Int_t GetMaxNStep() const = 0;
    virtual Int_t GetMedium() const =0;

        // tracking particle 
        // dynamic properties
    virtual void    TrackPosition(TLorentzVector& position) const =0;
    virtual void    TrackMomentum(TLorentzVector& momentum) const =0;
    virtual Float_t TrackStep() const =0;
    virtual Float_t TrackLength() const =0; 
    virtual Float_t TrackTime() const =0;
    virtual Float_t Edep() const =0;
        // static properties
    virtual Int_t   TrackPid() const =0;
    virtual Float_t TrackCharge() const =0;
    virtual Float_t TrackMass() const =0;
    virtual Float_t Etot() const =0;

        // track status
    virtual Bool_t  IsNewTrack() const =0;
    virtual Bool_t  IsTrackInside() const =0;
    virtual Bool_t  IsTrackEntering() const =0;
    virtual Bool_t  IsTrackExiting() const =0;
    virtual Bool_t  IsTrackOut() const =0;
    virtual Bool_t  IsTrackDisappeared() const =0;
    virtual Bool_t  IsTrackStop() const =0;
    virtual Bool_t  IsTrackAlive() const=0;

        // secondaries
    virtual Int_t NSecondaries() const=0;
    virtual void  GetSecondary(Int_t isec, Int_t& particleId, 
                    TLorentzVector& position, TLorentzVector& momentum) =0;
    virtual AliMCProcess ProdProcess(Int_t isec) const =0; 
    virtual Int_t StepProcesses(TArrayI &proc) const = 0;

    //
    // other (then geometry/step/run management) methods
    // ----------------------------------------------
    //
    
    //
    // Geant3 specific methods
    // !!! need to be transformed to common interface
    //
    virtual void Gdopt(const char*,const char*) = 0;
    virtual void SetClipBox(const char*,Float_t=-9999,Float_t=0, Float_t=-9999,Float_t=0,Float_t=-9999,Float_t=0) = 0;
    virtual void DefaultRange() = 0;
    virtual void Gdhead(Int_t, const char*, Float_t=0) = 0;   
    virtual void Gdman(Float_t, Float_t, const char*) = 0;
    virtual void SetColors() = 0;
    virtual void Gtreve() = 0;
    virtual void GtreveRoot() = 0;
    virtual void Gckmat(Int_t, char*) = 0;
    virtual void InitLego() = 0;
    virtual void Gfpart(Int_t, char*, Int_t&, Float_t&, Float_t&, Float_t&) = 0; 
    virtual void Gspart(Int_t, const char*, Int_t, Float_t, Float_t, Float_t) = 0; 

  // Control Methods

  virtual void Init() = 0;
  virtual void FinishGeometry() = 0;
  virtual void BuildPhysics() = 0;
  virtual void ProcessEvent() = 0;
  virtual void ProcessRun(Int_t nevent) = 0;
  // External Decayer
  virtual void SetExternalDecayer(AliDecayer* decayer) =0;
  virtual AliDecayer* Decayer() const =0;

  private:
  static AliMC*  fgMC;    // Pointer to the virtual MonteCarlo object
  AliMC(const AliMC &) {}
  AliMC & operator=(const AliMC &) {return (*this);}

  ClassDef(AliMC,1)  //Virtual MonteCarlo Interface
};

#endif 

