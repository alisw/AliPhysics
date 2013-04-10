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
  virtual Int_t GetLayer() const;
  virtual Int_t GetLadder() const;
  virtual Int_t GetDetector() const;
  virtual void  GetDetectorID(Int_t &layer,Int_t &ladder,Int_t &det) const;
  virtual void  GetPositionL(Float_t &x,Float_t &y,Float_t &z,Float_t &tof);
  virtual void  GetPositionL(Float_t &x,Float_t &y,Float_t &z) {Float_t tf;GetPositionL(x,y,z,tf);}
  virtual void  GetPositionL(Double_t &x,Double_t &y,Double_t &z,Double_t &t) {Float_t xf,yf,zf,tf;GetPositionL(xf,yf,zf,tf);x=xf,y=yf;z=zf;t=tf;}
  virtual void  GetPositionL(Double_t &x,Double_t &y,Double_t &z) {Float_t xf,yf,zf,tf;GetPositionL(xf,yf,zf,tf);x=xf,y=yf;z=zf;}
  virtual void  GetPositionL0(Double_t &x,Double_t &y,Double_t &z,Double_t &t);
  //
  virtual void Print(Option_t *option="") const;
  //
 protected:

  ClassDef(AliITSUHit,1)  //Hits object
	 
}; 

#endif
