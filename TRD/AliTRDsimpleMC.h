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
 
#include "AliMC.h"
#include "AliMCProcess.h"

class AliTRDv1;
 
class AliTRDsimpleMC : public AliMC {
 
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
  virtual void          Gfmate(Int_t imat, char *name, Float_t &a, Float_t &z,  
  		               Float_t &dens, Float_t &radl, Float_t &absl,
		               Float_t* ubuf, Int_t& nbuf)                                { }; 

  // Detector composition
  virtual void          Material(Int_t& kmat, const char* name, Float_t a, 
                                 Float_t z, Float_t dens, Float_t radl, Float_t absl,
                                 Float_t* buf, Int_t nwbuf)                               { };
  virtual void          Mixture(Int_t& kmat, const char *name, Float_t *a, 
                                Float_t *z, Float_t dens, Int_t nlmat, Float_t *wmat)     { };
  virtual void          Medium(Int_t& kmed, const char *name, Int_t nmat, 
                               Int_t isvol, Int_t ifield, Float_t fieldm, Float_t tmaxfd, 
                               Float_t stemax, Float_t deemax, Float_t epsil, 
		               Float_t stmin, Float_t* ubuf, Int_t nbuf)                  { };
  virtual void          Matrix(Int_t& krot, Float_t thetaX, Float_t phiX, 
                               Float_t thetaY, Float_t phiY, Float_t thetaZ, 
	  	               Float_t phiZ)                                              { };
  virtual void          Gstpar(Int_t itmed, const char *param, Float_t parval)            { }; 

  // Functions from GGEOM 
  virtual Int_t         Gsvolu(const char *name, const char *shape, Int_t nmed,  
		               Float_t *upar, Int_t np)                                   { return 0; }; 
  virtual void          Gsdvn(const char *name, const char *mother, Int_t ndiv, 
                              Int_t iaxis)                                                { }; 
  virtual void          Gsdvn2(const char *name, const char *mother, Int_t ndiv, 
                               Int_t iaxis, Float_t c0i, Int_t numed)                     { }; 
  virtual void          Gsdvt(const char *name, const char *mother, Float_t step, 
                              Int_t iaxis, Int_t numed, Int_t ndvmx)                      { }; 
  virtual void          Gsdvt2(const char *name, const char *mother, Float_t step, 
                               Int_t iaxis, Float_t c0, Int_t numed, Int_t ndvmx)         { }; 
  virtual void          Gsord(const char *name, Int_t iax)                                { }; 
  virtual void          Gspos(const char *name, Int_t nr, const char *mother,  
                              Float_t x, Float_t y, Float_t z, Int_t irot, 
                              const char *konly="ONLY")                                   { }; 
  virtual void          Gsposp(const char *name, Int_t nr, const char *mother,  
                               Float_t x, Float_t y, Float_t z, Int_t irot,
                               const char *konly, Float_t *upar, Int_t np)                { };
  virtual void          Gsbool(const char* onlyVolName, const char* manyVolName) {}

  virtual void          SetCerenkov(Int_t itmed, Int_t npckov, Float_t *ppckov,
                                    Float_t *absco, Float_t *effic, Float_t *rindex)      { };
    
  // Functions for drawing
  virtual void          DrawOneSpec(const char* name)                                     { };
  virtual void          Gsatt(const char* name, const char* att, Int_t val)               { };
  virtual void          Gdraw(const char*,Float_t theta = 30, Float_t phi = 30,
	                      Float_t psi = 0, Float_t u0 = 10, Float_t v0 = 10,
	                      Float_t ul = 0.01, Float_t vl = 0.01)                       { };

  // Euclid
  virtual void          WriteEuclid(const char *a, const char *b, Int_t c, Int_t d)       { };
		               
  // Get methods
  virtual Int_t         VolId(const Text_t* volName) const;
  virtual const char   *VolName(Int_t id) const                                           { return " "; };
  virtual Int_t         NofVolumes() const                                                { return 0; };
  virtual Int_t         VolId2Mate(Int_t id) const                                        { return 0; };

  //
  // Methods for physics management
  //
 
  // Set methods
  virtual void          SetCut(const char* cutName, Float_t cutValue)                     { };
  virtual void          SetProcess(const char* flagName, Int_t flagValue)                 { };
  virtual Float_t       Xsec(char *a, Float_t b, Int_t c, Int_t d)                        { return 0.0; }; 
 
  // Particle table usage         
  virtual Int_t         IdFromPDG(Int_t id) const                                         { return 0;   };  
  virtual Int_t         PDGFromId(Int_t pdg) const                                        { return 0;   };  
  virtual void          DefineParticles()                                                 { };      
  
  //
  // Methods for step management
  //

  // Action methods
  virtual void          StopTrack() { };
  virtual void          StopEvent() { };   

  // Set methods
  virtual void          SetMaxStep(Float_t step)                                          { fMaxStep = step; };
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
  virtual void          Gmtod(Float_t* xm, Float_t* xd, Int_t iflag)                      { };
  virtual void          Gdtom(Float_t* xd, Float_t* xm, Int_t iflag)                      { } ;
  virtual Float_t       MaxStep() const                                                   { return fMaxStep; };
  virtual Int_t         GetNStep() const                                                  { return fNStep;   };
  virtual Int_t         GetMaxNStep() const                                               { return 0;  };
  virtual Int_t         GetMedium() const                                                 { return 0;  };

  // Dynamic properties
  virtual void          TrackPosition(TLorentzVector& position) const;
  virtual void          TrackMomentum(TLorentzVector& momentum) const;
  virtual Float_t       TrackStep() const                                                 { return fTrackStep; };
  virtual Float_t       TrackLength() const                                               { return 0.0; }; 
  virtual Float_t       TrackTime() const                                                 { return 0.0; };
  virtual Float_t       Edep() const                                                      { return 0.0; };
  
  // Static properties
  virtual Int_t         TrackPid() const                                                  { return fTrackPid;    };
  virtual Float_t       TrackCharge() const                                               { return fTrackCharge; };
  virtual Float_t       TrackMass() const                                                 { return fTrackMass;   };
  virtual Float_t       Etot() const                                                      { return fTrackEtot;   };

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
  virtual AliMCProcess  ProdProcess(Int_t isec) const                                     { return kPNoProcess; }; 
  virtual Int_t         StepProcesses(TArrayI &proc) const                                { return 0; };

  //
  // Other (then geometry/step/run management) methods
  //
    
  // Geant3 specific methods
  virtual void          Gdopt(const char *c1,const char*c2)                               { };
  virtual void          SetClipBox(const char* cc,Float_t a=-9999,Float_t b=0, 
                                   Float_t c=-9999,Float_t d=0,
                                   Float_t e=-9999,Float_t f=0)                           { };
  virtual void          DefaultRange()                                                    { };
  virtual void          Gdhead(Int_t, const char *c, Float_t a=0)                         { };   
  virtual void          Gdman(Float_t a, Float_t b, const char *c)                        { };
  virtual void          SetColors()                                                       { };
  virtual void          Gtreve()                                                          { };
  virtual void          GtreveRoot()                                                      { };
  virtual void          Gckmat(Int_t, char*)                                              { };
  virtual void          InitLego()                                                        { };
  virtual void          Gfpart(Int_t a, char *b, Int_t &c, Float_t &d, 
                               Float_t &e, Float_t &g)                                    { }; 
  virtual void          Gspart(Int_t a, const char *b, Int_t c, Float_t d, 
                               Float_t e, Float_t g)                                      { }; 

  // Control Methods
  virtual void          Init()                                                            { };
  virtual void          FinishGeometry()                                                  { };
  virtual void          BuildPhysics()                                                    { };
  virtual void          ProcessEvent();
  virtual void          ProcessRun(Int_t nevent)                                          { };
  virtual AliMCGeomType     GetMCGeomType() const { return kGeant3; }

  // External Decayer
  virtual void          SetExternalDecayer(AliDecayer* decayer)                           { };
  virtual AliDecayer   *Decayer() const                                                   { return 0; };

 protected:

  enum {
      kVolDrRg
    , kVolAmRg
    , kVolDrCh
  };

  Float_t         fMaxStep;            //  Maximum step size
  Int_t           fNStep;              //  Number of steps
  Int_t           fTrack;              //  Track number
  Double_t        fTrackPx;            //  Track px
  Double_t        fTrackPy;            //  Track py
  Double_t        fTrackPz;            //  Track pz
  Double_t        fTrackPtot;          //  Track total momentum
  Double_t        fTrackEtot;          //  Track total energy
  Double_t        fTrackX;             //  Track x position
  Double_t        fTrackY;             //  Track y position
  Double_t        fTrackZ;             //  Track z position
  Double_t        fX0;                 //  X position of the beginning of the chamber
  Double_t        fTrackStep;          //  Track step size
  Int_t           fTrackPid;           //  Track PID
  Float_t         fTrackCharge;        //  Track charge
  Float_t         fTrackMass;          //  Track particle mass
  Bool_t          fTrackEntering;      //  Track entering chamber

  AliTRDv1       *fTRD;                //! TRD detector object

  ClassDef(AliTRDsimpleMC,1)           //  Simple TRD Monte Carlo class
 
};
#endif                                                                          
