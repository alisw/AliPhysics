#ifndef AliGEODB_H 
#define AliGEODB_H 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////// 
//  C++ interface to Geant3 basic routines    // 
//////////////////////////////////////////////// 
 
#include <AliMC.h> 
#include "AliGNode.h"
  
class AliGEODB : public AliMC { 

public: 
  AliGEODB() {}
  AliGEODB(const char *title, Int_t nwgeant=0); 
  virtual ~AliGEODB() {} 

///////////////////////////////////////////////////////////////////////
//                                                                   //
//                                                                   //
//     Here are the service routines from the geometry               //
//     which could be implemented also in other geometries           //
//                                                                   //
//                                                                   //
///////////////////////////////////////////////////////////////////////

  void  GeomIter();
  Int_t CurrentMaterial(Float_t &a, Float_t &z, Float_t &dens, Float_t &radl, Float_t &absl) const;
  Int_t NextVolUp(Text_t *name, Int_t &copy);
  Int_t CurrentVol(Text_t *name, Int_t &copy) const;
  Int_t CurrentVolOff(Int_t off, Text_t *name, Int_t &copy) const;
  Int_t VolId(Text_t *name) const;
  const char* VolName(Int_t id) const;
  void  TrackPosition(Float_t *xyz) const;
  void  TrackMomentum(Float_t *xyz) const;  
  Int_t NofVolumes() const;
  Float_t TrackTime() const;  
  Float_t TrackCharge() const;
  Float_t TrackMass() const;
  Float_t TrackStep() const;
  Float_t TrackLength() const;
  Int_t   TrackPid() const;
  Bool_t TrackInside() const;
  Bool_t TrackEntering() const;
  Bool_t TrackExiting() const;
  Bool_t TrackOut() const;
  Bool_t TrackDisappear() const;
  Bool_t TrackStop() const;
  Bool_t TrackAlive() const;
  Int_t   NSecondaries() const;
  Int_t   CurrentEvent() const;
  void    ProdProcess(char*) const;
  void    GetSecondary(Int_t, Int_t&, Float_t*, Float_t*);
  void   StopTrack();
  void   StopEvent();
  Float_t MaxStep() const;
  void   SetColors();
  void  SetMaxStep(Float_t maxstep);
  void  SetMaxNStep(Int_t maxnstp);
  Int_t GetMaxNStep() const;
  void GetParticle(const Int_t ipart, char *name, Float_t &mass) const;
  virtual Int_t GetMedium() const;
  virtual Float_t Edep() const;
  virtual Float_t Etot() const;
  virtual void    Rndm(Float_t* r, const Int_t n) const;
  virtual void    Material(Int_t&, const char*, Float_t, Float_t, Float_t, Float_t,
			    Float_t, Float_t* buf=0, Int_t nwbuf=0);
  virtual void    Mixture(Int_t&, const char*, Float_t*, Float_t*, Float_t, Int_t, Float_t*);
  virtual void    Medium(Int_t&, const char*, Int_t, Int_t, Int_t, Float_t, Float_t, 
		   Float_t, Float_t, Float_t, Float_t, Float_t* ubuf=0, Int_t nbuf=0);
  virtual void    Matrix(Int_t&, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t);

      // functions from GBASE 
   virtual  void  Gpcxyz(); 
   virtual  void  Ggclos(); 

   virtual  void  SetVisibility(Text_t* name, Int_t val);
   virtual  void  DrawTree( AliGNode* topnode, Int_t tabs );

   virtual  void  Gfile(const char *filename, const char *option="I"); 
   virtual  void  Glast(); 
   virtual  void  Gprint(const char *name); 
   virtual  void  Grun(); 
   virtual  void  Gtrig(); 
   virtual  void  Gtrigc(); 
   virtual  void  Gtrigi(); 
   virtual  void  Gwork(Int_t nwork); 
   virtual  void  Gzinit(); 
 
      // functions from GCONS 
   virtual  void  Gfmate(Int_t imat, char *name, Float_t &a, Float_t &z, Float_t &dens, 
                         Float_t &radl, Float_t &absl, Float_t* ubuf, Int_t& nbuf); 
   virtual  void  Gfpart(Int_t ipart, char *name, Int_t &itrtyp,  
                         Float_t &amass, Float_t &charge, Float_t &tlife); 
   virtual  void  Gftmed(Int_t numed, char *name, Int_t &nmat, Int_t &isvol,  
                         Int_t &ifield, Float_t &fieldm, Float_t &tmaxfd, 
                         Float_t &stemax, Float_t &deemax, Float_t &epsil, 
                         Float_t &stmin, Float_t *buf=0, Int_t *nbuf=0); 
   virtual  void  Gmate(); 
   virtual  void  Gpart(); 
   virtual  void  Gsckov(Int_t itmed, Int_t npckov, Float_t *ppckov,
			 Float_t *absco, Float_t *effic, Float_t *rindex); 
   virtual  void  Gsdk(Int_t ipart, Float_t *bratio, Int_t *mode); 
   virtual  void  Gsmate(Int_t imat, const char *name, Float_t a, Float_t z,  
                         Float_t dens, Float_t radl, Float_t absl); 
   virtual  void  Gsmixt(Int_t imat, const char *name, Float_t *a, Float_t *z,  
                         Float_t dens, Int_t nlmat, Float_t *wmat); 
   virtual  void  Gspart(Int_t ipart, const char *name, Int_t itrtyp,  
                         Float_t amass, Float_t charge, Float_t tlife); 
   virtual  void  Gstmed(Int_t numed, const char *name, Int_t nmat, Int_t isvol,  
                         Int_t ifield, Float_t fieldm, Float_t tmaxfd, 
                         Float_t stemax, Float_t deemax, Float_t epsil, 
                         Float_t stmin); 
   virtual  void  Gstpar(Int_t itmed, const char *param, Float_t parval); 
 
      // functions from GKINE 
   virtual  void  Gfkine(Int_t itra, Float_t *vert, Float_t *pvert, 
                         Int_t &ipart, Int_t &nvert); 
   virtual  void  Gfvert(Int_t nvtx, Float_t *v, Int_t &ntbeam, Int_t &nttarg, Float_t &tofg); 
   virtual  Int_t Gskine(Float_t *plab, Int_t ipart, Int_t nv, Float_t *ubuf=0, Int_t nwbuf=0); 
   virtual  Int_t Gsvert(Float_t *v, Int_t ntbeam, Int_t nttarg, Float_t *ubuf=0, Int_t nwbuf=0); 
 
      // functions from GPHYS 
   virtual  void  Gphysi(); 
 
      // functions from GTRAK 
   virtual  void  Gdebug(); 
   virtual  void  Gekbin(); 
   virtual  void  Gfinds(); 
   virtual  void  Gsking(Int_t igk); 
   virtual  void  Gskpho(Int_t igk); 
   virtual  void  Gsstak(Int_t iflag); 
   virtual  void  Gsxyz(); 
   virtual  void  Gtrack(); 
   virtual  void  Gtreve(); 
   virtual  void  Grndm(Float_t *rvec, const Int_t len) const; 
   virtual  void  Grndmq(Int_t &is1, Int_t &is2, const Int_t iseq, const Text_t *chopt); 
 
      // functions from GGEOM 
   virtual  void  Gdxyz(Int_t ); 
   virtual  void  Gdcxyz(); 

      // functions from GGEOM 
   virtual  void  Gdtom(Float_t *xd, Float_t *xm, Int_t iflag); 
   virtual  void  Glmoth(const char* iudet, Int_t iunum, Int_t &nlev, 
                         Int_t *lvols, Int_t *lindx); 
   virtual  void  Gmedia(Float_t *x, Int_t &numed); 
   virtual  void  Gmtod(Float_t *xm, Float_t *xd, Int_t iflag); 
   virtual  void  Gsdvn(const char *name, const char *mother, Int_t ndiv, Int_t iaxis); 
   virtual  void  Gsdvn2(const char *name, const char *mother, Int_t ndiv, Int_t iaxis, Float_t c0i, Int_t numed); 
   virtual  void  Gsdvs(const char *name, const char *mother, Float_t step, Int_t iaxis, Int_t numed); 
   virtual  void  Gsdvs2(const char *name, const char *mother, Float_t step, Int_t iaxis, Float_t c0, Int_t numed); 
   virtual  void  Gsdvt(const char *name, const char *mother, Float_t step, Int_t iaxis, Int_t numed, Int_t ndvmx); 
   virtual  void  Gsdvt2(const char *name, const char *mother, Float_t step, Int_t iaxis,
			 Float_t c0, Int_t numed, Int_t ndvmx); 
   virtual  void  Gsord(const char *name, Int_t iax); 
   virtual  void  Gspos(const char *name, Int_t nr, const char *mother,  
                         Float_t x, Float_t y, Float_t z, Int_t irot, const char *konly="ONLY"); 
   virtual  void  Gsposp(const char *nodename, Int_t nr, const char *mother,  
                         Float_t x, Float_t y, Float_t z, Int_t irot, const char *konly, Float_t *upar, Int_t np); 
   virtual  void  Gsrotm(Int_t nmat, Float_t theta1, Float_t phi1, Float_t theta2, Float_t phi2, 
                         Float_t theta3, Float_t phi3); 
   virtual  void  Gprotm(Int_t nmat=0); 
   virtual  Int_t Gsvolu(const char *name, const char *shape, Int_t nmed,  
                         Float_t *upar, Int_t np); 
   virtual  void  Gsatt(const char *name, const char *att, Int_t val);
   virtual  void  Gfpara(const char *name, Int_t number, Int_t intext, Int_t& npar,
			 Int_t& natt, Float_t* par, Float_t* att);
   virtual  void  Gckpar(Int_t, Int_t, Float_t*);
   virtual  void  Gckmat(Int_t, char*);
    
      // functions from GDRAW 
   virtual  void  DefaultRange();
   virtual  void  InitHIGZ();
   virtual  void  Gdopen(Int_t view);
   virtual  void  Gdclose();
   virtual  void  Gdelete(Int_t view);
   virtual  void  Gdshow(Int_t view);
   virtual  void  Gdopt(const char *name,const char *value);
   virtual  void  Gdraw(const char *name,Float_t theta=30, Float_t phi=30, Float_t psi=0,Float_t u0=10,Float_t v0=10,Float_t ul=0.01,Float_t vl=0.01);
   virtual  void  Gdrawc(const char *name,Int_t axis=1, Float_t cut=0,Float_t u0=10,Float_t v0=10,Float_t ul=0.01,Float_t vl=0.01);
   virtual  void  Gdrawx(const char *name,Float_t cutthe, Float_t cutphi, Float_t cutval,
                         Float_t theta=30, Float_t phi=30,Float_t u0=10,Float_t v0=10,Float_t ul=0.01,Float_t vl=0.01);
   virtual  void  Gdhead(Int_t isel, const char *name, Float_t chrsiz=0.6);   
   virtual  void  Gdman(Float_t u0, Float_t v0, const char *type="MAN");
   virtual  void  Gdspec(const char *name);
   virtual  void  DrawOneSpec(const char *name);
   virtual  void  Gdtree(const char *name,Int_t levmax=15,Int_t ispec=0);
   virtual  void  GdtreeParent(const char *name,Int_t levmax=15,Int_t ispec=0);

   virtual  void  WriteEuclid(const char*, const char*, Int_t, Int_t);

   virtual  void  SetABAN(Int_t par=1);
   virtual  void  SetANNI(Int_t par=1);
   virtual  void  SetAUTO(Int_t par=1);
   virtual  void  SetBOMB(Float_t bomb=1);
   virtual  void  SetBREM(Int_t par=1);
   virtual  void  SetCKOV(Int_t par=1);
   virtual  void  SetClipBox(const char *name,Float_t xmin=-9999,Float_t xmax=0, Float_t ymin=-9999,Float_t ymax=0,Float_t zmin=-9999,Float_t zmax=0);
   virtual  void  SetCOMP(Int_t par=1);
   virtual  void  SetCUTS(Float_t cutgam,Float_t cutele,Float_t cutneu,Float_t cuthad,
                      Float_t cutmuo ,Float_t bcute ,Float_t bcutm ,Float_t dcute ,
                      Float_t dcutm ,Float_t ppcutm, Float_t tofmax);
   virtual  void  SetDCAY(Int_t par=1);
   virtual  void  SetDEBU(Int_t emin=1, Int_t emax=999, Int_t emod=1);
   virtual  void  SetDRAY(Int_t par=1);
   virtual  void  SetHADR(Int_t par=1);
   virtual  void  SetKINE(Int_t kine, Float_t xk1=0, Float_t xk2=0, Float_t xk3=0, Float_t xk4=0,
                         Float_t xk5=0, Float_t xk6=0, Float_t xk7=0, Float_t xk8=0, Float_t xk9=0,
                         Float_t xk10=0);
   virtual  void  SetLOSS(Int_t par=2);
   virtual  void  SetMULS(Int_t par=1);
   virtual  void  SetMUNU(Int_t par=1);
   virtual  void  SetOPTI(Int_t par=2);
   virtual  void  SetPAIR(Int_t par=1);
   virtual  void  SetPFIS(Int_t par=1);
   virtual  void  SetPHOT(Int_t par=1);
   virtual  void  SetRAYL(Int_t par=1);
   virtual  void  SetSWIT(Int_t sw, Int_t val=1);
   virtual  void  SetTRIG(Int_t nevents=1);

   virtual  void  Vname(const char *name, char *vname);

   virtual  void  InitLego();
        
   ClassDef(AliGEODB,1)  //C++ interface to Geant basic routines 
}; 

#endif 
