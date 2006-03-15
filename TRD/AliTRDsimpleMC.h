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

#include <RVersion.h>
#include <TMCProcess.h>
#include <TVirtualMC.h>
#include "AliDecayer.h"

class AliTRDv1;
class TArrayD;

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

  virtual void          Copy(TObject &m) const;

  //
  // Methods for building / management of geometry
  //

  // Functions from GCONS
    virtual void  Gfmate(Int_t , char* , Float_t& , Float_t& ,
  		         Float_t& , Float_t& , Float_t& ,
		         Float_t* , Int_t& )  {}
    virtual void  Gfmate(Int_t , char* , Double_t& , Double_t& ,
  		         Double_t& , Double_t& , Double_t& ,
		         Double_t* , Int_t& ) {}

  // Detector composition
    virtual void  Material(Int_t& , const char* , Double_t ,
                     Double_t , Double_t , Double_t , Double_t ,
                     Float_t* , Int_t ) {}
    virtual void  Material(Int_t& , const char* , Double_t ,
                     Double_t , Double_t , Double_t , Double_t ,
                     Double_t* , Int_t ) {}
    virtual void  Mixture(Int_t& , const char* , Float_t* ,
                     Float_t *, Double_t , Int_t , Float_t* ) {}
    virtual void  Mixture(Int_t& , const char* , Double_t *,
                     Double_t *, Double_t , Int_t , Double_t* ) {}
    virtual void  Medium(Int_t& , const char* , Int_t ,
                     Int_t , Int_t , Double_t , Double_t ,
                     Double_t , Double_t , Double_t ,
		     Double_t , Float_t* , Int_t ) {}
    virtual void  Medium(Int_t& , const char* , Int_t ,
                     Int_t , Int_t , Double_t , Double_t ,
                     Double_t , Double_t , Double_t ,
		     Double_t , Double_t* , Int_t ) {}
    virtual void  Matrix(Int_t& , Double_t , Double_t ,
                     Double_t , Double_t , Double_t ,
		     Double_t ) {}
    virtual void  Gstpar(Int_t , const char* , Double_t ) {}

  // Functions from GGEOM
    virtual Int_t  Gsvolu(const char* , const char* , Int_t ,
                          Float_t* , Int_t )  { return 0; }
    virtual Int_t  Gsvolu(const char* , const char* , Int_t ,
                          Double_t* , Int_t) { return 0; }
    virtual void  Gsdvn(const char* , const char* , Int_t ,
                         Int_t ) {}
    virtual void  Gsdvn2(const char* , const char* , Int_t ,
                         Int_t , Double_t , Int_t ) {}
    virtual void  Gsdvt(const char* , const char* , Double_t ,
                         Int_t , Int_t , Int_t ) {}
    virtual void  Gsdvt2(const char* , const char* , Double_t ,
                         Int_t , Double_t , Int_t , Int_t ) {}
    virtual void  Gsord(const char* , Int_t ) {}
    virtual void  Gspos(const char* , Int_t , const char* ,
                         Double_t , Double_t , Double_t , Int_t ,
                         const char* ) {}
    virtual void  Gsposp(const char* , Int_t , const char* ,
                         Double_t, Double_t, Double_t, Int_t ,
                         const char* , Float_t* , Int_t ) {}
    virtual void  Gsposp(const char* , Int_t , const char* ,
                         Double_t , Double_t , Double_t , Int_t ,
                         const char* , Double_t* , Int_t ) {}
    virtual void  Gsbool(const char* , const char* ) {}

    virtual void  SetCerenkov(Int_t , Int_t , Float_t* ,
                               Float_t* , Float_t* , Float_t* ) {}
    virtual void  SetCerenkov(Int_t , Int_t , Double_t* ,
                               Double_t* , Double_t* , Double_t* ) {}

  // Functions for drawing
    virtual void  DrawOneSpec(const char* ) {}
    virtual void  Gsatt(const char* , const char* , Int_t ) {}
    virtual void  Gdraw(const char* , Double_t , Double_t ,
		        Double_t , Double_t , Double_t ,
		        Double_t , Double_t ) {}

  // Euclid
  virtual void          WriteEuclid(const char* , const char* , Int_t , Int_t ) {}

  // Get methods
    virtual Int_t VolId(const Text_t* volName) const;
    virtual const char* VolName(Int_t ) const { return ""; }
    virtual Int_t NofVolumes() const { return 0; }
    virtual Int_t VolId2Mate(Int_t ) const { return 0; }

  //
  // Methods for physics management
  //

  // Set methods
#if ROOT_VERSION_CODE > 262150
    virtual Bool_t   SetCut(const char* , Double_t ) { return kTRUE; }
    virtual Bool_t   SetProcess(const char* , Int_t ) { return kTRUE; }
#else
    virtual void     SetCut(const char* , Double_t ) {}
    virtual void     SetProcess(const char* , Int_t ) {}
#endif
    virtual void     DefineParticles() {}
    virtual Double_t Xsec(char*, Double_t, Int_t, Int_t) { return 0.; }

  // Particle table usage
    virtual Int_t   IdFromPDG(Int_t ) const { return 0; }
    virtual Int_t   PDGFromId(Int_t ) const { return 0; }

  //
  // Methods for step management
  //

  // Action methods
  virtual void          StopTrack() { };
  virtual void          StopEvent() { };
#if ROOT_VERSION_CODE >= 262150
  virtual void          StopRun()   { }
#endif

  // Set methods
  virtual void          SetMaxStep(Double_t step)                                         { fMaxStep = step; };
  virtual void          SetMaxNStep(Int_t )                                              { };
  virtual void          SetUserDecay(Int_t )                                             { };

  virtual void          NewTrack(Int_t iTrack, Int_t pdg, Double_t px, Double_t py, Double_t pz);

  // Tracking volume(s)
  virtual Int_t         CurrentVolID(Int_t& copyNo) const;
  virtual Int_t         CurrentVolOffID(Int_t off, Int_t& copyNo) const;
  virtual const char*   CurrentVolName() const;
  virtual const char*   CurrentVolOffName(Int_t ) const                                { return ""; };
  virtual Int_t         CurrentMaterial(Float_t& , Float_t& ,
                                        Float_t& , Float_t& ,
			  	        Float_t& ) const                               { return 0;  };
  virtual Int_t         CurrentEvent() const                                              { return 0;  };
  virtual void          Gmtod(Float_t* , Float_t* , Int_t )   {}
  virtual void          Gmtod(Double_t* , Double_t* , Int_t ) {}
  virtual void          Gdtom(Float_t* , Float_t* , Int_t )   {}
  virtual void          Gdtom(Double_t* , Double_t* , Int_t ) {}
  virtual Double_t      MaxStep() const                                                   { return fMaxStep; };
  virtual Int_t         GetNStep() const                                                  { return fNStep;   };
  virtual Int_t         GetMaxNStep() const                                               { return 0;  };
  virtual Int_t         GetMedium() const                                                 { return 0;  };
  virtual Bool_t GetMedium(const TString& /*volumeName*/,
			   TString& /*name*/, Int_t& /*imed*/,
			   Int_t& /*nmat*/, Int_t& /*isvol*/, Int_t& /*ifield*/,
			   Double_t& /*fieldm*/, Double_t& /*tmaxfd*/, Double_t& /*stemax*/,
			   Double_t& /*deemax*/, Double_t& /*epsil*/, Double_t& /*stmin*/,
			   TArrayD& /*par*/) {
   return kFALSE;
}   

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
  virtual void          GetSecondary(Int_t , Int_t& ,
                                     TLorentzVector& ,
                                     TLorentzVector& )                                    { };
  virtual Bool_t        SecondariesAreOrdered() const {return kTRUE;}
  virtual TMCProcess    ProdProcess(Int_t ) const                                           { return kPNoProcess; };
  virtual Int_t         StepProcesses(TArrayI& ) const                                    { return 0; };

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
    // Dummy methods
#if ROOT_VERSION_CODE > 197895
#if ROOT_VERSION_CODE > 262150
    virtual Bool_t DefineParticle(int, const char*, TMCParticleType, double, double, double){ return kTRUE; }
    virtual Bool_t DefineIon(const char*, int, int, int, double, double) { return kTRUE; }
#else
    virtual void DefineParticle(int, const char*, TMCParticleType, double, double, double){;}
    virtual void DefineIon(const char*, int, int, int, double, double){;}
#endif
    virtual TString  ParticleName(int) const {return "";}
    virtual Double_t ParticleMass(int) const {return 0.;}
    virtual Double_t ParticleCharge(int) const {return 0.;}
    virtual Double_t ParticleLifeTime(int) const {return 0.;}
    virtual TMCParticleType ParticleMCType(int) const {return (TMCParticleType) 0;}
#endif
    //
  // Control Methods
  virtual void          Init()                                                            { };
  virtual void          FinishGeometry()                                                  { };
  virtual void          BuildPhysics()                                                    { };
  virtual void          ProcessEvent();
#if ROOT_VERSION_CODE >= 262150
  virtual Bool_t        ProcessRun(Int_t )                                                { return kTRUE; }
#else
  virtual void          ProcessRun(Int_t )                                                { };
#endif
  //virtual TMCGeomType   GetMCGeomType() const                                             { return kGeant3; }

  // External Decayer
  virtual void          SetExternalDecayer(AliDecayer* )                                  { };
  virtual AliDecayer   *Decayer() const                                                   { return 0; };

#if ROOT_VERSION_CODE>=262913
  virtual void SetRootGeometry() {}
  virtual Int_t NofVolDaughters(const char*) const {return 0;}
  virtual const char* VolDaughterName(const char*, Int_t) const {return 0x0;}
  virtual Int_t VolDaughterCopyNo(const char*, Int_t) const {return 0;}
  virtual void ForceDecayTime(Float_t) {}
  virtual const char* CurrentVolPath() {return 0x0;}
#endif

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

  ClassDef(AliTRDsimpleMC, 3)            //  Simple TRD Monte Carlo class

};
#endif

