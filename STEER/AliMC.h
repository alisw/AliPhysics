#ifndef ALIMC_H
#define ALIMC_H
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                                                                           //
//    Generic interface to MC for AliRoot                                    //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include <TLorentzVector.h>

class AliMC : public TNamed 
{

private:
  static AliMC* fgMC;

public:
  AliMC(const char *name, const char *title);
  AliMC();
  virtual ~AliMC() {}
  //Generic access functions
  static inline AliMC* GetMC() {return fgMC;}
  //
  virtual Int_t   CurrentMaterial(Float_t &a, Float_t &z, Float_t &dens, Float_t &radl, Float_t &absl) const =0;
  virtual Int_t   CurrentVolID(Int_t&) const =0;
  virtual Int_t   CurrentVolOffID(Int_t, Int_t& ) const =0;
  virtual const char* CurrentVolName() const =0;
  virtual const char* CurrentVolOffName(Int_t) const =0;
  virtual Int_t   NofVolumes() const =0;
  virtual Int_t   VolId(Text_t*) const =0;
  virtual Int_t   IdFromPDG(Int_t) const =0;
  virtual Int_t   PDGFromId(Int_t) const =0;
  virtual void    DefineParticles() = 0;
  virtual void    TrackPosition(TLorentzVector&) const =0;
  virtual void    TrackMomentum(TLorentzVector&) const =0;
  virtual Float_t TrackCharge() const =0;
  virtual Float_t TrackMass() const =0;
  virtual Float_t TrackStep() const =0;
  virtual Int_t   TrackPid() const =0;
  virtual Bool_t  IsTrackInside() const =0;
  virtual Bool_t  IsTrackEntering() const =0;
  virtual Bool_t  IsTrackExiting() const =0;
  virtual Bool_t  IsTrackOut() const =0;
  virtual Bool_t  IsTrackDisappeared() const =0;
  virtual Bool_t  IsTrackStop() const =0;
  virtual Float_t TrackLength() const =0;
  virtual Float_t TrackTime() const =0;
  virtual Bool_t  IsTrackAlive() const=0;
  virtual Int_t   NSecondaries() const=0;
  virtual Int_t   CurrentEvent() const=0;
  virtual void    ProdProcess(char*) const=0;
  virtual void    GetSecondary(Int_t, Int_t&, Float_t*, Float_t*)=0;
  virtual void    StopTrack() =0;
  virtual void    StopEvent() =0;
  virtual Float_t MaxStep() const =0;
  virtual void    SetMaxStep(Float_t ) =0;
  virtual void    SetMaxNStep(Int_t) =0;
  virtual Int_t   GetMaxNStep() const =0;
  virtual void    SetUserDecay(Int_t) =0;
  virtual void    GetParticle(const Int_t, char*, Float_t&) const =0;
  virtual Int_t   GetMedium() const =0;
  virtual void    DrawOneSpec(const char*)=0;
  virtual Float_t Edep() const =0;
  virtual Float_t Etot() const =0;
  virtual const char*   VolName(Int_t) const=0;
  virtual void    Gstpar(Int_t, const char *, Float_t)=0;
  virtual Int_t   Gsvolu(const char*, const char*, Int_t, Float_t*, Int_t)=0;
  virtual void    Gsdvn(const char*, const char*, Int_t, Int_t)=0;
  virtual void    Gsdvn2(const char*, const char*, Int_t, Int_t,
			 Float_t, Int_t)=0;
  virtual void    Gsdvt(const char*, const char*, Float_t, Int_t,
			Int_t, Int_t)=0; 
  virtual  void   Gsdvt2(const char*, const char*, Float_t, Int_t,
			 Float_t, Int_t, Int_t)=0; 
  virtual  void   Gspos(const char*, Int_t, const char*, Float_t, Float_t,
			Float_t, Int_t, const char *konly="ONLY")=0; 
  virtual  void   Gsposp(const char*, Int_t, const char*, Float_t, Float_t,
			 Float_t, Int_t, const char*, Float_t*, Int_t)=0;
  virtual  void   Gfmate(Int_t, char*, Float_t&, Float_t&, Float_t&, 
                         Float_t&, Float_t&, Float_t*, Int_t&)=0; 
  virtual  void   Gsatt(const char*, const char*, Int_t)=0;
  virtual  void   Gmtod(Float_t*, Float_t*, Int_t)=0;
  virtual  void   Gdtom(Float_t*, Float_t*, Int_t)=0;
  virtual  void   Gdraw(const char*,Float_t theta=30, Float_t phi=30,
			Float_t psi=0,Float_t u0=10,Float_t v0=10,
			Float_t ul=0.01,Float_t vl=0.01)=0;
  virtual  void   WriteEuclid(const char*, const char*, Int_t, Int_t)=0;
  virtual void    Rndm(Float_t*, const Int_t) const=0;
  virtual void    Material(Int_t&, const char*, Float_t, Float_t, Float_t, Float_t,
			    Float_t, Float_t* buf=0, Int_t nwbuf=0) =0;
  virtual void    Mixture(Int_t&, const char*, Float_t*, Float_t*, Float_t, Int_t, Float_t*) =0;
  virtual void    Medium(Int_t&, const char*, Int_t, Int_t, Int_t, Float_t, Float_t, 
		   Float_t, Float_t, Float_t, Float_t, Float_t *ubuf=0, Int_t nbuf=0) =0;
  virtual void    Matrix(Int_t&, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t) =0;

  //________________________________________________________________________
  //
  //                     Geant3 specific prototypes
  //
  //------------------------------------------------------------------------

  virtual void Gdopt(const char*,const char*)=0;
  virtual  void  SetClipBox(const char*,Float_t=-9999,Float_t=0, Float_t=-9999,Float_t=0,Float_t=-9999,Float_t=0)=0;
  virtual  void  DefaultRange()=0;
  virtual  void  Gdhead(Int_t, const char*, Float_t=0)=0;   
  virtual  void  Gdman(Float_t, Float_t, const char*)=0;
  virtual  void  Gsord(const char *name, Int_t iax)=0;
  virtual  void  Ggclos()=0;
  virtual  void  SetColors()=0;
  virtual  void  Gphysi()=0;
  virtual  void  Gtrigi()=0;
  virtual  void  Gtreve()=0;
  virtual  void  Gtreve_root()=0;
  virtual  void  Gtrigc()=0;
  virtual  void  Gtrig()=0;
  virtual  void  Gckmat(Int_t, char*)=0;
  virtual  void  InitLego()=0;
  virtual  void  Gfpart(Int_t, char*, Int_t&, Float_t&, Float_t&, Float_t&)=0; 
  virtual  void  Gspart(Int_t, const char*, Int_t, Float_t, Float_t, Float_t)=0; 

  ClassDef(AliMC,1) //Generic MonteCarlo Class

};

R__EXTERN AliMC *gMC;

#endif

