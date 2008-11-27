// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2007

//=============================================================================
// parent class of all events; analyzers access data via this class
//=============================================================================

#ifndef ALIDEVENT_H
#define ALIDEVENT_H

#include <cmath>
#include <TObject.h>

class TTree;

//=============================================================================
class AliDEvent : public TObject {
 public:
  AliDEvent() : TObject(){printf("%s object named %s created\n",ClassName(),GetName());}
  virtual ~AliDEvent()   {printf("%s object named %s deleted\n",ClassName(),GetName());}

  // interface part

  virtual void        AttachTree(TTree *tr) = 0;
  virtual Double_t    Etamin() const = 0;     // experiment's acceptance
  virtual Double_t    Etamax() const = 0;
  virtual Bool_t      Good() const = 0;  
  virtual Double_t    Centrality() = 0;       // centrality (0,1); 0 is most central
  virtual void        RP(Double_t &qx, Double_t &qy) const = 0;
  virtual Double_t    RPphi() const = 0;
  virtual Double_t    Zver() const = 0;       // z-vertex (-1,1)
  virtual Int_t       NParticles() const = 0;

  virtual Bool_t      ParticleGood(Int_t i, Int_t pidi) const = 0;
  virtual Double_t    ParticleP(Int_t i) const = 0;
  virtual Double_t    ParticleTheta(Int_t i) const = 0;
  virtual Double_t    ParticlePhi(Int_t i) const = 0;
  virtual Double_t    ParticleDedx(Int_t i) const = 0;
  virtual Bool_t      PairGood(Double_t p0, Double_t the0, Double_t phi0, 
			       Double_t p1, Double_t the1, Double_t phi1) const = 0;

  // toolkit part

  void     RP(Double_t &qx, Double_t &qy, Int_t harmonic) const;
  Double_t ParticlePt(Int_t i) const {return ParticleP(i)*sin(ParticleTheta(i));}
  Double_t ParticlePz(Int_t i) const {return ParticleP(i)*cos(ParticleTheta(i));}
  Double_t ParticleEta(Int_t i) const;
  Double_t ParticleY(Int_t i, Double_t mass) const; 

  ClassDef(AliDEvent,0)
};
#endif 
//=============================================================================
