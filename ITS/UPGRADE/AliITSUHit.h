#ifndef ALIITSUHIT_H
#define ALIITSUHIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//
// At the moment the same functionality/data-members as parent AliITShit 
// except the geometry transformation uses AliITSgeomTGeoUp 
//
////////////////////////////////////////////////////////////////////////

#include "AliITShit.h" 
#include "AliRun.h"

class AliITSUHit : public AliITShit {

 public:
  //
  AliITSUHit() {}
  AliITSUHit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  AliITSUHit(Int_t shunt,Int_t track,Int_t *vol,Float_t edep,Float_t tof,TLorentzVector &x,TLorentzVector &x0,TLorentzVector &p);
  AliITSUHit(const AliITSUHit &h);
  AliITSUHit& operator=(const AliITSUHit &h);
  virtual ~AliITSUHit() {}

  void SetChip(Int_t chip) {SetModule(chip);}
  Int_t GetChip()          {return GetModule();}

  virtual Int_t GetLayer() const;
  virtual Int_t GetStave() const;
  virtual Int_t GetHalfStave() const;
  virtual Int_t GetModule() const;  
  virtual Int_t GetChipInModule() const;
  virtual void  GetChipID(Int_t &layer,Int_t &stave,Int_t &sstave, Int_t &mod, Int_t &det) const;
  virtual void  GetPositionL(Float_t &x,Float_t &y,Float_t &z,Float_t &tof);
  virtual void  GetPositionL(Float_t &x,Float_t &y,Float_t &z) {Float_t tf;GetPositionL(x,y,z,tf);}
  virtual void  GetPositionL(Double_t &x,Double_t &y,Double_t &z,Double_t &t) {Float_t xf,yf,zf,tf;GetPositionL(xf,yf,zf,tf);x=xf,y=yf;z=zf;t=tf;}
  virtual void  GetPositionL(Double_t &x,Double_t &y,Double_t &z) {Float_t xf,yf,zf,tf;GetPositionL(xf,yf,zf,tf);x=xf,y=yf;z=zf;}
  virtual void  GetPositionL0(Double_t &x,Double_t &y,Double_t &z,Double_t &t);
  //
  virtual void Print(Option_t *option="") const;
  //
 protected:
  virtual void SetModule(Int_t mod){fModule=mod;};
  virtual Int_t GetModule(){return fModule;};

  ClassDef(AliITSUHit,1)  //Hits object
	 
}; 

#endif
