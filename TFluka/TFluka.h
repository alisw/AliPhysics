#ifndef TFLUKA_H
#define TFLUKA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
// FLUKA implementation of the AliMC Interface                               //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliMC.h"

class TFluka : public AliMC
{
 public:
    TFluka(const char *title) 
	:AliMC("TFluka", title) {;}
    TFluka() {;}
    virtual ~TFluka() {;}
  
    //
    // methods for building/management of geometry
    // ------------------------------------------------
    //
    // Not yet implemented !!!
    //
    // functions from GCONS 
    virtual void  Gfmate(Int_t imat, char *name, Float_t &a, Float_t &z,  
  		         Float_t &dens, Float_t &radl, Float_t &absl,
		         Float_t* ubuf, Int_t& nbuf)
	{printf("Gfmate not yet implemented !\n");}
    

    // detector composition
    virtual void  Material(Int_t& kmat, const char* name, Float_t a, 
                     Float_t z, Float_t dens, Float_t radl, Float_t absl,
			   Float_t* buf, Int_t nwbuf)
	{printf("Material not yet implemented !\n");}
    
    virtual void  Mixture(Int_t& kmat, const char *name, Float_t *a, 
			  Float_t *z, Float_t dens, Int_t nlmat, Float_t *wmat)
	{printf("Mixture not yet implemented !\n");}
    
    virtual void  Medium(Int_t& kmed, const char *name, Int_t nmat, 
			 Int_t isvol, Int_t ifield, Float_t fieldm,
			 Float_t tmaxfd, 
			 Float_t stemax, Float_t deemax, Float_t epsil, 
			 Float_t stmin, Float_t* ubuf, Int_t nbuf)
	{printf("Medium not yet implemented !\n");}
    virtual void  Matrix(Int_t& krot, Float_t thetaX, Float_t phiX, 
			 Float_t thetaY, Float_t phiY, Float_t thetaZ, 
			 Float_t phiZ)
	{printf("Matrix not yet implemented !\n");}
    virtual void  Gstpar(Int_t itmed, const char *param, Float_t parval)
	{printf("Gstpar not yet implemented !\n");} 
    
    // functions from GGEOM 
    virtual Int_t Gsvolu(const char *name, const char *shape, Int_t nmed,  
                         Float_t *upar, Int_t np)
	{printf("Gsvolu not yet implemented !\n"); return -1;} 
    virtual void  Gsdvn(const char *name, const char *mother, Int_t ndiv, 
                        Int_t iaxis)
	{printf("Gsdvn not yet implemented !\n");} 
    virtual void  Gsdvn2(const char *name, const char *mother, Int_t ndiv, 
                         Int_t iaxis, Float_t c0i, Int_t numed)
	{printf("Gsdvn2 not yet implemented !\n");} 
    virtual void  Gsdvt(const char *name, const char *mother, Float_t step, 
                        Int_t iaxis, Int_t numed, Int_t ndvmx)
	{printf("Gsdvt not yet implemented !\n");} 
    virtual void  Gsdvt2(const char *name, const char *mother, Float_t step, 
                         Int_t iaxis, Float_t c0, Int_t numed, Int_t ndvmx)
	{printf("Gsdvt2 not yet implemented !\n");} 
    virtual void  Gsord(const char *name, Int_t iax)
	{printf("Gsord not yet implemented !\n");} 
    virtual void  Gspos(const char *name, Int_t nr, const char *mother,  
                        Float_t x, Float_t y, Float_t z, Int_t irot, 
                        const char *konly="ONLY")
	{printf("Gspos not yet implemented !\n");} 
    virtual void  Gsposp(const char *name, Int_t nr, const char *mother,  
                         Float_t x, Float_t y, Float_t z, Int_t irot,
                         const char *konly, Float_t *upar, Int_t np)
	{printf("Gsposp not yet implemented !\n");}
    virtual void  Gsbool(const char* onlyVolName, const char* manyVolName)
	{printf("Gsbool not yet implemented !\n");}

    virtual void  SetCerenkov(Int_t itmed, Int_t npckov, Float_t *ppckov,
			      Float_t *absco, Float_t *effic, Float_t *rindex)
	{printf("SetCerenkov not yet implemented !\n");}
    
    
    // functions for drawing
    virtual void  DrawOneSpec(const char* name)
	{printf("DrawOneSpec not yet implemented !\n");}
    virtual void  Gsatt(const char* name, const char* att, Int_t val)
	{printf("Gsatt not yet implemented !\n");}
    virtual void  Gdraw(const char*,Float_t theta = 30, Float_t phi = 30,
		        Float_t psi = 0, Float_t u0 = 10, Float_t v0 = 10,
		        Float_t ul = 0.01, Float_t vl = 0.01)
	{printf("Gdraw not yet implemented !\n");}

    // Euclid
    virtual void WriteEuclid(const char*, const char*, Int_t, Int_t)
	{printf("WriteEuclid not yet implemented !\n");}
		               
    // get methods
    virtual Int_t VolId(const Text_t* volName) const
	{printf("VolId not yet implemented !\n"); return -1;}
    virtual const char* VolName(Int_t id) const
	{printf("VolName not yet implemented !\n"); return "void";}
    virtual Int_t NofVolumes() const
	{printf("NofVolumes not yet implemented !\n"); return -1;}
    virtual Int_t VolId2Mate(Int_t id) const
	{printf("VolIdMate not yet implemented !\n"); return -1;}

    //
    // methods for physics management
    // ------------------------------------------------
    //
 
    // set methods
    virtual void SetCut(const char* cutName, Float_t cutValue)
	{printf("Set Cut not yet implemented !\n");}
    virtual void SetProcess(const char* flagName, Int_t flagValue)
	{printf("Set Process not yet implemented !\n");}
    virtual Float_t Xsec(char*, Float_t, Int_t, Int_t)
	{printf("Xsec not yet implemented !\n"); return -1.;} 
 
        // particle table usage         
    virtual Int_t   IdFromPDG(Int_t id)  const
	 {printf("IdFromPDG not yet implemented !\n"); return -1;}
    virtual Int_t   PDGFromId(Int_t pdg) const
	{printf("PDGFromId not yet implemented !\n"); return -1;}  
    virtual void    DefineParticles()
	{printf("DefineParticles not yet implemented !\n");}      
  
    //
    // methods for step management
    // ------------------------------------------------
    //

    // action methods
    virtual void StopTrack()
	{printf("StopTrack not yet implemented !\n");}
    
    virtual void StopEvent()
	{printf("StopEvent not yet implemented !\n");}   

    // set methods
    virtual void SetMaxStep(Float_t)
	{printf("SetMaxStep not yet implemented !\n");}
    virtual void SetMaxNStep(Int_t)
	{printf("SetMaxNStep not yet implemented !\n");}
    virtual void SetUserDecay(Int_t)
	{printf("SetUserDecay not yet implemented !\n");}  

    // get methods
         // tracking volume(s) 
    virtual Int_t CurrentVolID(Int_t& copyNo) const
	{printf("CurrentVolID not yet implemented !\n"); return -1;}
    virtual Int_t CurrentVolOffID(Int_t off, Int_t& copyNo) const
	{printf("CurrentVolOffID not yet implemented !\n"); return -1;}
    virtual const char* CurrentVolName() const
	{printf("CurrentVolName not yet implemented !\n"); return "void";}
    virtual const char* CurrentVolOffName(Int_t off) const
	{printf("CurrentVolOffName not yet implemented !\n"); return "void";}
    virtual Int_t CurrentMaterial(Float_t &a, Float_t &z, 
                    Float_t &dens, Float_t &radl, Float_t &absl) const
	{printf("CurrentMaterial not yet implemented !\n"); return -1;}  
    virtual Int_t CurrentEvent() const
	{printf("CurrentEvent not yet implemented !\n"); return -1;} 
    virtual void  Gmtod(Float_t* xm, Float_t* xd, Int_t iflag)	
	{printf("Gmtod not yet implemented !\n");}
    virtual void  Gdtom(Float_t* xd, Float_t* xm, Int_t iflag)
	{printf("Gdtom not yet implemented !\n");}
    virtual Float_t MaxStep() const
	{printf("MaxStep not yet implemented !\n"); return -1.;}
    virtual Int_t GetMaxNStep() const
	{printf("GetMaxNStep not yet implemented !\n"); return -1;}
    virtual Int_t GetMedium() const
	{printf("GetMedium not yet implemented !\n"); return -1;}

        // tracking particle 
        // dynamic properties
    virtual void    TrackPosition(TLorentzVector& position) const
	{printf("TrackPosition not yet implemented !\n");}
    virtual void    TrackMomentum(TLorentzVector& momentum) const
	{printf("TrackMomentum not yet implemented !\n");}
    virtual Float_t TrackStep() const
	{printf("TrackStep not yet implemented !\n"); return -1.;}
    virtual Float_t TrackLength() const
	{printf("TrackLength not yet implemented !\n"); return -1.;} 
    virtual Float_t TrackTime() const
	{printf("TrackTimenot yet implemented !\n"); return -1.;}
    virtual Float_t Edep() const
	{printf("Edep not yet implemented !\n"); return -1.;}
        // static properties
    virtual Int_t   TrackPid() const
	{printf("TrackPid not yet implemented !\n"); return -1;}
    virtual Float_t TrackCharge() const
	{printf("TrackCharge not yet implemented !\n"); return -1.;}
    virtual Float_t TrackMass() const
	{printf("TrackMasd not yet implemented !\n"); return -1.;}
    virtual Float_t Etot() const
	{printf("Etot not yet implemented !\n"); return -1.;}
    
        // track status
    virtual Bool_t  IsNewTrack() const
	{printf("IsNewTrack not yet implemented !\n"); return 0;}
    virtual Bool_t  IsTrackInside() const
	{printf("IsTrackInside not yet implemented !\n"); return 0;}
    virtual Bool_t  IsTrackEntering() const
	{printf("IsTrackEntering not yet implemented !\n"); return 0;}
    virtual Bool_t  IsTrackExiting() const
	{printf("IsTrackExiting not yet implemented !\n"); return 0;}
    virtual Bool_t  IsTrackOut() const
	{printf("IsTrackOut not yet implemented !\n"); return 0;}
    virtual Bool_t  IsTrackDisappeared() const
	{printf("IsTrackDisappeared not yet implemented !\n"); return 0;}
    virtual Bool_t  IsTrackStop() const
	{printf("IsTrackStop not yet implemented !\n"); return 0;}
    virtual Bool_t  IsTrackAlive() const
	{printf("IsTrackAlive not yet implemented !\n"); return 0;}

        // secondaries
    virtual Int_t NSecondaries() const
	{printf("NSecondaries not yet implemented !\n"); return -1;}
    virtual void  GetSecondary(Int_t isec, Int_t& particleId, 
                    TLorentzVector& position, TLorentzVector& momentum)
	{printf("GetSecondary not yet implemented !\n");}
    virtual AliMCProcess ProdProcess(Int_t isec) const
	{printf("ProdProcess not yet implemented !\n"); return kPNoProcess;} 
    virtual Int_t StepProcesses(TArrayI &proc) const
	{printf("StepProcess not yet implemented !\n"); return -1;}

    //
    // other (then geometry/step/run management) methods
    // ----------------------------------------------
    //
    virtual AliMCGeomType GetMCGeomType() const
	{printf("AliMCGeomType not yet implemented !\n"); return kFluka;}
    
    //
    // Geant3 specific methods
    // !!! need to be transformed to common interface
    //
    virtual void Gdopt(const char*,const char*)
	{printf("Gdopt not yet implemented !\n");}
    virtual void SetClipBox(const char*,Float_t=-9999,Float_t= 0, Float_t=-9999,Float_t=0,Float_t=-9999,Float_t=0)
	{printf("SetClipBox not yet implemented !\n");}
    virtual void DefaultRange()
	{printf("DefaultRange not yet implemented !\n");}
    virtual void Gdhead(Int_t, const char*, Float_t=0)
	{printf("Gdhead not yet implemented !\n");}   
    virtual void Gdman(Float_t, Float_t, const char*)
	{printf("Gdman not yet implemented !\n");}
    virtual void SetColors()
	{printf("SetColors not yet implemented !\n");}
    virtual void Gtreve()
	{printf("Gtreve not yet implemented !\n");}
    virtual void GtreveRoot()
	{printf("GtreveRoot not yet implemented !\n");}
    virtual void Gckmat(Int_t, char*)
	{printf("Gckmat not yet implemented !\n");}
    virtual void InitLego()
	{printf("InitLego not yet implemented !\n");}
    virtual void Gfpart(Int_t, char*, Int_t&, Float_t&, Float_t&, Float_t&)
	{printf("Gfpart not yet implemented !\n");} 
    virtual void Gspart(Int_t, const char*, Int_t, Float_t, Float_t, Float_t)
	{printf("Gspart not yet implemented !\n");} 

  // Control Methods

  
  virtual void FinishGeometry()
      {printf("FinishGeometry not yet implemented !\n");}
  virtual void BuildPhysics()
      {printf("BuildPhysics not yet implemented !\n");}
  // External Decayer
  virtual void SetExternalDecayer(AliDecayer* decayer)
      {printf("SetExternalDecayer not yet implemented !\n");}
  virtual AliDecayer* Decayer() const
      {printf("AliDecayer not yet implemented !\n"); return 0;}

  private:
  TFluka(const TFluka &) {;}
  TFluka & operator=(const TFluka &) {return (*this);}

  ClassDef(TFluka,1)  //Fluka implementation of the Monte Carlo Interface
};

#endif 

