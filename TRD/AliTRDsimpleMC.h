#ifndef ALITRDSIMPLEMC_H
#define ALITRDSIMPLEMC_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/* $Id$ */
 
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Simple TRD Monte Carlo class                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TMCProcess.h"
#include "TVirtualMC.h"
 
#include "AliDecayer.h"

class AliTRDv1;
class AliTRDparameter;
 
class AliTRDsimpleMC : public TVirtualMC {
 
 public:     

  enum {
      kPdgElectron = 11
    , kPdgPion     = 211
  };

  AliTRDsimpleMC();
  AliTRDsimpleMC(const char *name, const char *title);
  AliTRDsimpleMC(const AliTRDsimpleMC &m); 
                                                                                
  virtual ~AliTRDsimpleMC();
  AliTRDsimpleMC &operator=(const AliTRDsimpleMC &m);    

  virtual void          Copy(TObject &m);

  //
  // Methods for building / management of geometry
  //

  // Functions from GCONS 
    virtual void  Gfmate(Int_t imat, char *name, Float_t &a, Float_t &z,  
  		         Float_t &dens, Float_t &radl, Float_t &absl,
		         Float_t* ubuf, Int_t& nbuf)  {}
    virtual void  Gfmate(Int_t imat, char *name, Double_t &a, Double_t &z,  
  		         Double_t &dens, Double_t &radl, Double_t &absl,
		         Double_t* ubuf, Int_t& nbuf) {}

  // Detector composition
    virtual void  Material(Int_t& kmat, const char* name, Double_t a, 
                     Double_t z, Double_t dens, Double_t radl, Double_t absl,
                     Float_t* buf, Int_t nwbuf) {}
    virtual void  Material(Int_t& kmat, const char* name, Double_t a, 
                     Double_t z, Double_t dens, Double_t radl, Double_t absl,
                     Double_t* buf, Int_t nwbuf) {}
    virtual void  Mixture(Int_t& kmat, const char *name, Float_t *a, 
                     Float_t *z, Double_t dens, Int_t nlmat, Float_t *wmat) {}
    virtual void  Mixture(Int_t& kmat, const char *name, Double_t *a, 
                     Double_t *z, Double_t dens, Int_t nlmat, Double_t *wmat) {}
    virtual void  Medium(Int_t& kmed, const char *name, Int_t nmat, 
                     Int_t isvol, Int_t ifield, Double_t fieldm, Double_t tmaxfd, 
                     Double_t stemax, Double_t deemax, Double_t epsil, 
		     Double_t stmin, Float_t* ubuf, Int_t nbuf) {}
    virtual void  Medium(Int_t& kmed, const char *name, Int_t nmat, 
                     Int_t isvol, Int_t ifield, Double_t fieldm, Double_t tmaxfd, 
                     Double_t stemax, Double_t deemax, Double_t epsil, 
		     Double_t stmin, Double_t* ubuf, Int_t nbuf) {}
    virtual void  Matrix(Int_t& krot, Double_t thetaX, Double_t phiX, 
                     Double_t thetaY, Double_t phiY, Double_t thetaZ, 
		     Double_t phiZ) {}
    virtual void  Gstpar(Int_t itmed, const char *param, Double_t parval) {} 

  // Functions from GGEOM 
    virtual Int_t  Gsvolu(const char *name, const char *shape, Int_t nmed,  
                          Float_t *upar, Int_t np)  { return 0; }
    virtual Int_t  Gsvolu(const char *name, const char *shape, Int_t nmed,  
                          Double_t *upar, Int_t np) { return 0; } 
    virtual void  Gsdvn(const char *name, const char *mother, Int_t ndiv, 
                         Int_t iaxis) {} 
    virtual void  Gsdvn2(const char *name, const char *mother, Int_t ndiv, 
                         Int_t iaxis, Double_t c0i, Int_t numed) {} 
    virtual void  Gsdvt(const char *name, const char *mother, Double_t step, 
                         Int_t iaxis, Int_t numed, Int_t ndvmx) {} 
    virtual void  Gsdvt2(const char *name, const char *mother, Double_t step, 
                         Int_t iaxis, Double_t c0, Int_t numed, Int_t ndvmx) {} 
    virtual void  Gsord(const char *name, Int_t iax) {} 
    virtual void  Gspos(const char *name, Int_t nr, const char *mother,  
                         Double_t x, Double_t y, Double_t z, Int_t irot, 
                         const char *konly="ONLY") {} 
    virtual void  Gsposp(const char *name, Int_t nr, const char *mother,  
                         Double_t x, Double_t y, Double_t z, Int_t irot,
                         const char *konly, Float_t *upar, Int_t np) {}
    virtual void  Gsposp(const char *name, Int_t nr, const char *mother,  
                         Double_t x, Double_t y, Double_t z, Int_t irot,
                         const char *konly, Double_t *upar, Int_t np) {}
    virtual void  Gsbool(const char* onlyVolName, const char* manyVolName) {}

    virtual void  SetCerenkov(Int_t itmed, Int_t npckov, Float_t *ppckov,
                               Float_t *absco, Float_t *effic, Float_t *rindex) {}
    virtual void  SetCerenkov(Int_t itmed, Int_t npckov, Double_t *ppckov,
                               Double_t *absco, Double_t *effic, Double_t *rindex) {}

  // Functions for drawing
    virtual void  DrawOneSpec(const char* name) {}
    virtual void  Gsatt(const char* name, const char* att, Int_t val) {}
    virtual void  Gdraw(const char*,Double_t theta = 30, Double_t phi = 30,
		        Double_t psi = 0, Double_t u0 = 10, Double_t v0 = 10,
		        Double_t ul = 0.01, Double_t vl = 0.01) {}

  // Euclid
  virtual void          WriteEuclid(const char *a, const char *b, Int_t c, Int_t d) {}
		               
  // Get methods
    virtual Int_t VolId(const Text_t* volName) const;
    virtual const char* VolName(Int_t id) const { return ""; }
    virtual Int_t NofVolumes() const { return 0; }
    virtual Int_t VolId2Mate(Int_t id) const { return 0; }

  //
  // Methods for physics management
  //
 
  // Set methods
    virtual void     SetCut(const char* cutName, Double_t cutValue) {}
    virtual void     SetProcess(const char* flagName, Int_t flagValue) {}
    virtual Double_t Xsec(char*, Double_t, Int_t, Int_t) { return 0.; } 
 
  // Particle table usage         
    virtual Int_t   IdFromPDG(Int_t id) const  { return 0; }  
    virtual Int_t   PDGFromId(Int_t pdg) const { return 0; }  
    virtual void    DefineParticles() {}      
  
  //
  // Methods for step management
  //

  // Action methods
  virtual void          StopTrack() { };
  virtual void          StopEvent() { };   

  // Set methods
  virtual void          SetMaxStep(Double_t step)                                         { fMaxStep = step; };
  virtual void          SetMaxNStep(Int_t n)                                              { };
  virtual void          SetUserDecay(Int_t d)                                             { };  

  virtual void          NewTrack(Int_t iTrack, Int_t pdg, Double_t px, Double_t py, Double_t pz);

  // Tracking volume(s) 
  virtual Int_t         CurrentVolID(Int_t& copyNo) const;
  virtual Int_t         CurrentVolOffID(Int_t off, Int_t& copyNo) const;
  virtual const char*   CurrentVolName() const                                            { return ""; };
  virtual const char*   CurrentVolOffName(Int_t off) const                                { return ""; };
  virtual Int_t         CurrentMaterial(Float_t &a, Float_t &z, 
                                        Float_t &dens, Float_t &radl, 
			  	        Float_t &absl) const                              { return 0;  };  
  virtual Int_t         CurrentEvent() const                                              { return 0;  }; 
  virtual void          Gmtod(Float_t* xm, Float_t* xd, Int_t iflag)   {}
  virtual void          Gmtod(Double_t* xm, Double_t* xd, Int_t iflag) {}
  virtual void          Gdtom(Float_t* xd, Float_t* xm, Int_t iflag)   {}
  virtual void          Gdtom(Double_t* xd, Double_t* xm, Int_t iflag) {}
  virtual Double_t      MaxStep() const                                                   { return fMaxStep; };
  virtual Int_t         GetNStep() const                                                  { return fNStep;   };
  virtual Int_t         GetMaxNStep() const                                               { return 0;  };
  virtual Int_t         GetMedium() const                                                 { return 0;  };

  // Dynamic properties
  virtual void          TrackPosition(TLorentzVector& position) const;
  virtual void          TrackPosition(Double_t &x, Double_t &y, Double_t &z) const;
  virtual void          TrackMomentum(TLorentzVector& momentum) const;
  virtual void          TrackMomentum(Double_t &px, Double_t &py, Double_t &pz, Double_t &etot) const;
  virtual Double_t      TrackStep() const                                                 { return fTrackStep; };
  virtual Double_t      TrackLength() const                                               { return 0.0; };
  virtual Double_t      TrackTime() const                                                 { return 0.0; };
  virtual Double_t      Edep() const                                                      { return 0.0; };
  
  // Static properties
  virtual Int_t         TrackPid() const                                                  { return fTrackPid;    };
  virtual Double_t      TrackCharge() const                                               { return fTrackCharge; };
  virtual Double_t      TrackMass() const                                                 { return fTrackMass;   };
  virtual Double_t      Etot() const                                                      { return fTrackEtot;   };

  // Track status
  virtual Bool_t        IsNewTrack() const                                                { return kFALSE; };
  virtual Bool_t        IsTrackInside() const                                             { return kFALSE; };
  virtual Bool_t        IsTrackEntering() const                                           { return fTrackEntering; };
  virtual Bool_t        IsTrackExiting() const                                            { return kFALSE; };
  virtual Bool_t        IsTrackOut() const                                                { return kFALSE; };
  virtual Bool_t        IsTrackDisappeared() const                                        { return kFALSE; };
  virtual Bool_t        IsTrackStop() const                                               { return kFALSE; };
  virtual Bool_t        IsTrackAlive() const                                              { return kFALSE; };

  // Secondaries
  virtual Int_t         NSecondaries() const                                              { return 0; };
  virtual void          GetSecondary(Int_t isec, Int_t& particleId, 
                                     TLorentzVector& position, 
                                     TLorentzVector& momentum)                            { };
  virtual TMCProcess  ProdProcess(Int_t isec) const                                     { return kPNoProcess; }; 
  virtual Int_t         StepProcesses(TArrayI &proc) const                                { return 0; };

  //
  // Other (then geometry/step/run management) methods
  //
    
  // Geant3 specific methods
    virtual void Gdopt(const char*,const char*) {}
    virtual void SetClipBox(const char*,Double_t=-9999,Double_t=0, Double_t=-9999,
                             Double_t=0,Double_t=-9999,Double_t=0) {}
    virtual void DefaultRange() {}
    virtual void Gdhead(Int_t, const char*, Double_t=0) {}   
    virtual void Gdman(Double_t, Double_t, const char*) {}
    virtual void SetColors() {}
    virtual void Gtreve() {}
    virtual void GtreveRoot() {}
    virtual void Gckmat(Int_t, char*) {}
    virtual void InitLego() {}
    virtual void Gfpart(Int_t, char*, Int_t&, Float_t&, Float_t&, Float_t&) {} 
    virtual void Gspart(Int_t, const char*, Int_t, Double_t, Double_t, Double_t) {} 

  // Control Methods
  virtual void          Init()                                                            { };
  virtual void          FinishGeometry()                                                  { };
  virtual void          BuildPhysics()                                                    { };
  virtual void          ProcessEvent();
  virtual void          ProcessRun(Int_t nevent)                                          { };
  //virtual TMCGeomType   GetMCGeomType() const                                             { return kGeant3; }

  // External Decayer
  virtual void          SetExternalDecayer(AliDecayer* decayer)                           { };
  virtual AliDecayer   *Decayer() const                                                   { return 0; };

 protected:

  enum {
      kVolDrRg
    , kVolAmRg
    , kVolDrCh
  };

  Float_t          fMaxStep;            //  Maximum step size
  Int_t            fNStep;              //  Number of steps
  Int_t            fTrack;              //  Track number
  Double_t         fTrackPx;            //  Track px
  Double_t         fTrackPy;            //  Track py
  Double_t         fTrackPz;            //  Track pz
  Double_t         fTrackPtot;          //  Track total momentum
  Double_t         fTrackEtot;          //  Track total energy
  Double_t         fTrackX;             //  Track x position
  Double_t         fTrackY;             //  Track y position
  Double_t         fTrackZ;             //  Track z position
  Double_t         fX0;                 //  X position of the beginning of the chamber
  Double_t         fTrackStep;          //  Track step size
  Int_t            fTrackPid;           //  Track PID
  Float_t          fTrackCharge;        //  Track charge
  Float_t          fTrackMass;          //  Track particle mass
  Bool_t           fTrackEntering;      //  Track entering chamber

  AliTRDv1        *fTRD;                //! TRD detector object
  AliTRDparameter *fPar;                //! TRD parameter object

  ClassDef(AliTRDsimpleMC,2)            //  Simple TRD Monte Carlo class
 
};
#endif                                                                          

