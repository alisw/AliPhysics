//////////////////////////////////////////////////////
//  C++ dummy interface to Geant3 basic routines    //
//////////////////////////////////////////////////////

#include "TGeant3.h"

//---------------------------------------------------------

ClassImp(TGeant3)

//___________________________________________
TGeant3::TGeant3(){}

//___________________________________________
TGeant3::TGeant3(const char *title, Int_t)
       :AliMC("TGeant3",title){}

//___________________________________________ 
void TGeant3::DefaultRange() {}
 
//___________________________________________
void TGeant3::LoadAddress() {}

//=======================functions from GBASE

//___________________________________________
void    TGeant3::Gfile(const char*, const char*) {}
void    TGeant3::GeomIter() {}
Int_t   TGeant3::CurrentMaterial(Float_t &, Float_t &, Float_t &, Float_t &, Float_t &) const {return 0;}
Int_t   TGeant3::NextVolUp(Text_t*, Int_t&) {return 0;}
Int_t   TGeant3::CurrentVolID(Int_t&) const {return 0;}
const char*  TGeant3::CurrentVolName() const {return 0;}
Int_t   TGeant3::NofVolumes() const {return 0;}
Int_t   TGeant3::CurrentVolOffID(Int_t, Int_t&) const {return 0;}
const char *TGeant3::CurrentVolOffName(Int_t) const {return 0;}
void    TGeant3::TrackPosition(TLorentzVector&) const {}
void    TGeant3::TrackMomentum(TLorentzVector&) const {}
Int_t   TGeant3::IdFromPDG(Int_t pdg) const {return -1;}
Int_t   TGeant3::PDGFromId(Int_t pdg) const {return -1;}
void    TGeant3::DefineParticles() {}
Int_t   TGeant3::VolId(Text_t*) const {return 0;}
const char*   TGeant3::VolName(Int_t ) const {return 0;}
Float_t TGeant3::TrackCharge() const {return 0;}
Float_t TGeant3::TrackMass() const {return 0;}
Bool_t  TGeant3::IsTrackInside() const {return 0;}
Bool_t  TGeant3::IsTrackEntering() const {return 0;}
Bool_t  TGeant3::IsTrackExiting() const {return 0;}
Bool_t  TGeant3::IsTrackOut() const {return 0;}
Bool_t  TGeant3::IsTrackDisappeared() const {return 0;}
Bool_t  TGeant3::IsTrackStop() const {return 0;}
Int_t   TGeant3::NSecondaries() const {return 0;}
void    TGeant3::ProdProcess(char*) const {}
void    TGeant3::GetSecondary(Int_t, Int_t&, Float_t*, Float_t*){}
Float_t TGeant3::MaxStep() const {return 0;}
void    TGeant3::SetMaxStep(Float_t ) {}
void    TGeant3::GetParticle(const Int_t, char*, Float_t&) const {}
Int_t   TGeant3::CurrentEvent() const {return 0;}
Int_t   TGeant3::GetMedium() const {return 0;}
Float_t TGeant3::Edep() const {return 0;}
Float_t TGeant3::Etot() const {return 0;}
void    TGeant3::Rndm(Float_t*, const Int_t) const {}
Float_t TGeant3::TrackStep() const {return 0;}
Float_t TGeant3::TrackLength() const {return 0;}
Float_t TGeant3::TrackTime() const {return 0;}
Int_t   TGeant3::TrackPid() const {return 0;}
Bool_t  TGeant3::IsTrackAlive() const {return 0;}
void    TGeant3::StopTrack() {}
void    TGeant3::StopEvent() {}
void    TGeant3::SetMaxNStep(Int_t) {}
void    TGeant3::SetColors() {}
Int_t   TGeant3::GetMaxNStep() const {return 0;}
void    TGeant3::Material(Int_t&, const char*, Float_t, Float_t, Float_t, Float_t,
			  Float_t, Float_t*, Int_t){}
void    TGeant3::Mixture(Int_t&, const char*, Float_t*, Float_t*, Float_t, Int_t, Float_t*){}
void    TGeant3::Medium(Int_t&, const char*, Int_t, Int_t, Int_t, Float_t, Float_t, 
			Float_t, Float_t, Float_t, Float_t, Float_t*, Int_t){}
void    TGeant3::Matrix(Int_t&, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t){}
void    TGeant3::WriteEuclid(const char*, const char*, Int_t, Int_t) {}


void  TGeant3::Gpcxyz() {}
void  TGeant3::Ggclos() {}
void  TGeant3::Glast() {}
void  TGeant3::Gprint(const char*) {}
void  TGeant3::Grun() {}
void  TGeant3::Gtrig() {}
void  TGeant3::Gtrigc() {}
void  TGeant3::Gtrigi() {}
void  TGeant3::Gwork(Int_t) {}
void  TGeant3::Gzinit() {}


//=======================functions from GCONS
//___________________________________________
void  TGeant3::Gfmate(Int_t, char*, Float_t&, Float_t&, Float_t&, Float_t&,
		      Float_t&, Float_t*, Int_t&){}

//___________________________________________
void  TGeant3::Gfpart(Int_t, char*, Int_t&, Float_t&, Float_t&, Float_t&){}

//___________________________________________
void  TGeant3::Gftmed(Int_t, char*, Int_t&, Int_t&, Int_t&, Float_t&, Float_t&,
                    Float_t&, Float_t&, Float_t&, Float_t&, Float_t*, Int_t*){}

//___________________________________________
void  TGeant3::Gmate() {}
void  TGeant3::Gpart() {}
void  TGeant3::Gsdk(Int_t, Float_t*, Int_t*) {}
void  TGeant3::Gsmate(Int_t, const char*, Float_t, Float_t, Float_t, Float_t, Float_t) {}
void  TGeant3::Gsmixt(Int_t, const char*, Float_t*, Float_t*, Float_t, Int_t, Float_t*) {}
void  TGeant3::Gspart(Int_t, const char*, Int_t,   Float_t, Float_t, Float_t) {}
void  TGeant3::Gstmed(Int_t, const char*, Int_t, Int_t, 
                   Int_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t) {}
void  TGeant3::Gstpar(Int_t, const char*, Float_t) {}
void  TGeant3::Gsckov(Int_t, Int_t, Float_t *, Float_t *, Float_t *, Float_t *) {}


//=======================functions from GKINE
//___________________________________________
void  TGeant3::Gfkine(Int_t, Float_t*, Float_t*, Int_t&, Int_t&) {}
void  TGeant3::Gfvert(Int_t, Float_t*, Int_t&, Int_t&, Float_t&) {}
Int_t TGeant3::Gskine(Float_t*, Int_t, Int_t, Float_t*, Int_t){return 0;}

//___________________________________________
Int_t TGeant3::Gsvert(Float_t*, Int_t, Int_t, Float_t*, Int_t){return 0;}

//=======================functions from GPHYS
//___________________________________________
void  TGeant3::Gphysi() {}


//=======================functions from GTRAK
//___________________________________________
void  TGeant3::Gdebug() {}
void  TGeant3::Gekbin() {}
void  TGeant3::Gfinds() {}
void  TGeant3::Gsking(Int_t) {}
void  TGeant3::Gskpho(Int_t) {}
void  TGeant3::Gsstak(Int_t) {}
void  TGeant3::Gsxyz() {}
void  TGeant3::Gtrack() {}
void  TGeant3::Gtreve() {}
void  TGeant3::Gtreve_root() {}
void  TGeant3::Grndm(Float_t*, const Int_t) const {}
void  TGeant3::Grndmq(Int_t&, Int_t&, const Int_t, const Text_t*) {}

//=======================functions from GDRAW
//___________________________________________
void  TGeant3::Gdxyz(Int_t ) {}
void  TGeant3::Gdcxyz() {}


//=======================functions from GGEOM
//___________________________________________
void  TGeant3::Gdtom(Float_t*, Float_t*, Int_t) {}
void  TGeant3::Glmoth(const char*, Int_t, Int_t&, Int_t*, Int_t*) {}
void  TGeant3::Gmedia(Float_t*, Int_t&) {}
void  TGeant3::Gmtod(Float_t*, Float_t*, Int_t) {}
void  TGeant3::Gsdvn(const char*, const char*, Int_t, Int_t) {}
void  TGeant3::Gsdvn2(const char*, const char*, Int_t, Int_t, Float_t, Int_t) {}
void  TGeant3::Gsdvs(const char*, const char*, Float_t, Int_t, Int_t) {}
void  TGeant3::Gsdvs2(const char*, const char*, Float_t, Int_t, Float_t, Int_t) {}
void  TGeant3::Gsdvt(const char*, const char*, Float_t, Int_t, Int_t, Int_t) {}
void  TGeant3::Gsdvt2(const char *, const char *, Float_t, Int_t, Float_t, Int_t, Int_t) {}
void  TGeant3::Gsord(const char*, Int_t) {}
void  TGeant3::Gspos(const char*, Int_t, const char*, Float_t, Float_t, Float_t, Int_t, const char*) {}
void  TGeant3::Gsposp(const char*, Int_t, const char*, Float_t, Float_t, Float_t, Int_t, const char*, Float_t*, Int_t) {}
void  TGeant3::Gsrotm(Int_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t)  {}
void  TGeant3::Gprotm(Int_t)  {}
Int_t TGeant3::Gsvolu(const char*, const char*, Int_t, Float_t*, Int_t){return 0;}
void  TGeant3::Gsatt(const char*, const char*, Int_t) {}
void  TGeant3::Gfpara(const char*, Int_t, Int_t, Int_t&, Int_t&, Float_t*, Float_t*) {}
void  TGeant3::Gckpar(Int_t, Int_t, Float_t*) {}
void  TGeant3::Gckmat(Int_t, char*) {}

//______________________________________________________________________________
void TGeant3::Streamer(TBuffer &R__b)
{
   // Stream an object of class TGeant3.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliMC::Streamer(R__b);
      R__b >> fNextVol;
   } else {
      R__b.WriteVersion(TGeant3::IsA());
      AliMC::Streamer(R__b);
      R__b << fNextVol;
   }
}

//=======================DRAW functions
//____________________________________
void  TGeant3::InitHIGZ() {}
void  TGeant3::Gdopen(Int_t) {}
void  TGeant3::Gdclose() {}
void  TGeant3::Gdelete(Int_t) {}
void  TGeant3::Gdshow(Int_t) {}
void  TGeant3::Gdopt(const char *,const char *) {}
void  TGeant3::Gdraw(const char *,Float_t, Float_t, Float_t,Float_t,Float_t,Float_t,Float_t) {}
void  TGeant3::Gdrawc(const char *,Int_t, Float_t,Float_t,Float_t,Float_t,Float_t) {}
void  TGeant3::Gdrawx(const char *,Float_t, Float_t, Float_t,
                         Float_t, Float_t,Float_t,Float_t,Float_t,Float_t) {}
void  TGeant3::Gdhead(Int_t, const char *, Float_t) {}
void  TGeant3::Gdman(Float_t, Float_t, const char *) {}
void  TGeant3::Gdspec(const char *) {}
void  TGeant3::DrawOneSpec(const char *) {}
void  TGeant3::Gdtree(const char *,Int_t,Int_t) {}
void  TGeant3::GdtreeParent(const char *,Int_t,Int_t) {}

//=======================Set functions
//____________________________________
void TGeant3::SetABAN(Int_t)  {}
void TGeant3::SetANNI(Int_t)  {}
void TGeant3::SetAUTO(Int_t)  {}
void TGeant3::SetBOMB(Float_t)  {}
void TGeant3::SetBREM(Int_t)  {}
void TGeant3::SetCKOV(Int_t)  {}
void TGeant3::SetClipBox(const char *,Float_t,Float_t, Float_t,Float_t,Float_t,Float_t) {}
void TGeant3::SetCOMP(Int_t)  {}
void TGeant3::SetCUTS(Float_t,Float_t,Float_t,Float_t,
                      Float_t,Float_t,Float_t,Float_t,
                      Float_t,Float_t, Float_t)  {}
void TGeant3::SetDCAY(Int_t)  {}
void TGeant3::SetDEBU(Int_t, Int_t, Int_t)  {}
void TGeant3::SetDRAY(Int_t)  {}
void TGeant3::SetHADR(Int_t)  {}
void TGeant3::SetKINE(Int_t, Float_t, Float_t, Float_t, Float_t,
                     Float_t, Float_t, Float_t, Float_t, Float_t,
                     Float_t)  {}
void TGeant3::SetLOSS(Int_t)  {}
void TGeant3::SetMULS(Int_t)  {}
void TGeant3::SetMUNU(Int_t)  {}
void TGeant3::SetOPTI(Int_t)  {}
void TGeant3::SetPAIR(Int_t)  {}
void TGeant3::SetPFIS(Int_t)  {}
void TGeant3::SetPHOT(Int_t)  {}
void TGeant3::SetRAYL(Int_t)  {}
void TGeant3::SetSWIT(Int_t , Int_t)  {}
void TGeant3::SetTRIG(Int_t)  {}
void TGeant3::SetUserDecay(Int_t) {}
void TGeant3::Vname(const char *, char *) {}
void TGeant3::InitLego() {}


#include "AliCallf77.h"
 

extern "C" type_of_call {

  // This is only for ITS

  void type_of_call gsvolu_(){}
  void type_of_call gspos_(){}
  void type_of_call gsatt_(){}
  void type_of_call gsdvn_(){}
  void type_of_call gsposp_(){}
  void type_of_call atg_(){}
  void type_of_call sxsrot_(){}

  // All this for ZDC

  void type_of_call grndm_(){}
  void type_of_call gphits_(){}
  void type_of_call gfhits_() {}
  void type_of_call gsvert_(){}
  void type_of_call gskine_(){}
  void type_of_call gfpart_(){}
  void type_of_call lorenf_(){}
  void type_of_call gpvert_(){}
  void type_of_call gpkine_(){}
  void type_of_call rnorml_(){}
  void type_of_call gfpath_() {}
  void type_of_call uctoh_() {}
  void type_of_call glvolu_() {}
  void type_of_call gmtod_() {}
  void type_of_call gfkine_() {}
  void type_of_call vmod_() {}
  void type_of_call gsahit_() {}
  void type_of_call gschit_() {}
  void type_of_call gdtom_() {}
  void type_of_call rnpssn_() {}
  void type_of_call ucopy_() {}
}
