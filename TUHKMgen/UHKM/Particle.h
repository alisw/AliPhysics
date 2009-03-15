/*                                                                            
                                                                            
        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru 
                           November. 2, 2005                                

*/

#ifndef PARTICLE_INCLUDED
#define PARTICLE_INCLUDED

#include <list>

#include <TLorentzRotation.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#ifndef PARTICLE_PDG
#include "ParticlePDG.h"
#endif
#include <iostream>

//class ParticlePDG;

class Particle {
 protected:
  TLorentzVector   fPosition;
  TLorentzVector   fMomentum;
  TLorentzVector   fLastMotherDecayCoor;
  TLorentzVector   fLastMotherDecayMom;
  ParticlePDG     *fParticleProperties;
  Double_t         fLastInteractionTime;
  Int_t            fInteractionNumber;
  Int_t            fLastMotherPdg;
  Int_t            fType; //0-hydro, 1-jets
  Int_t            fIndex;                    // index (0 based) of particle in the final particle list which will contain both primaries and secondaries
  Int_t            fMotherIndex;              // index of the mother (-1 if its a primary particle)
  Int_t            fNDaughters;               // number of daughter particles (0 if the particle had not decayed)
  Int_t            fDaughterIndex[3];         // array of indexes for the daughter particles (the indexes are -1 for non-existing daughters)
  static Int_t     fLastIndex;                // the last index assigned

 public:
  Particle(const TLorentzVector &, const TLorentzVector &);
  Particle(ParticlePDG *pdg = 0);
  Particle(ParticlePDG *pdg, const TLorentzVector &pos, const TLorentzVector &mom,
	   Double_t lastInterTime = 0., Int_t lastInterNum = 0, Int_t type=0);
  Particle(ParticlePDG *pdg, const TLorentzVector &pos, const TLorentzVector &mom,
	   Double_t lastInterTime, Int_t lastInterNum, Int_t type, Int_t motherPdg, 
	   const TLorentzVector &motherPos, const TLorentzVector &motherMom);

  Double_t X()const{return fPosition.X();}
  Double_t X(Double_t val){fPosition.SetX(val); return val;}
  Double_t Y()const{return fPosition.Y();}
  Double_t Y(Double_t val){fPosition.SetY(val); return val;}
  Double_t Z()const{return fPosition.Z();}
  Double_t Z(Double_t val){fPosition.SetZ(val); return val;}
  Double_t T()const{return fPosition.T();}
  Double_t T(Double_t val){fPosition.SetT(val); return val;}
  Double_t Px()const{return fMomentum.Px();}
  Double_t Px(Double_t val){fMomentum.SetPx(val); return val;}
  Double_t Py()const{return fMomentum.Py();}
  Double_t Py(Double_t val){fMomentum.SetPy(val); return val;}
  Double_t Pz()const{return fMomentum.Pz();}
  Double_t Pz(Double_t val){fMomentum.SetPz(val); return val;}
  Double_t E()const{return fMomentum.E();}
  Double_t E(Double_t val){fMomentum.SetE(val); return val;}

  TLorentzVector &Pos(){return fPosition;}
  const TLorentzVector &Pos()const{return fPosition;}
  TLorentzVector &Pos(const TLorentzVector &val){return fPosition = val;}
  TLorentzVector &Mom(){return fMomentum;}
  const TLorentzVector &Mom()const{return fMomentum;}
  TLorentzVector &Mom(const TLorentzVector &val){return fMomentum = val;}
   		
  void Boost(const TVector3 &val){fMomentum.Boost(val);}
  void Boost(const TLorentzVector &val){fMomentum.Boost(val.BoostVector());}
  void TransformMomentum(const TRotation &rotator){fMomentum *= rotator;}
  void TransformPosition(const TRotation &rotator){fPosition *= rotator;}
  void Shift(const TVector3 &val){fPosition += TLorentzVector(val, 0.);}

  //Pseudorapidity
  Double_t Eta ()const;
  //Rapidity
  Double_t Rapidity()const;
  Double_t Phi()const;
  Double_t Theta()const;
  Double_t Pt()const;

  Int_t Encoding() const;
  Double_t TableMass() const;
  ParticlePDG *Def() const {return fParticleProperties;}
  ParticlePDG *Def(ParticlePDG *newProp) {return fParticleProperties = newProp;}
  //mother   
  void SetLastMotherPdg(Int_t value){fLastMotherPdg = value;}
  Int_t GetLastMotherPdg() const {return fLastMotherPdg;}

  // aic(2008/08/08): functions added in order to enable tracking of mother/daughter particles by a unique index
  // The index coincides with the position of the particle in the secondaries list.
  Int_t SetIndex() {fIndex = ++fLastIndex; return fIndex;}
  Int_t GetIndex() {return fIndex;}
  static Int_t GetLastIndex() {return fLastIndex;}
  void InitIndexing() {fLastIndex = -1;}
  void SetMother(Int_t value) {fMotherIndex = value;}
  Int_t GetMother() {return fMotherIndex;}
  void SetDaughter(Int_t value) {
    if(fNDaughters==3) {
      std::cout << "Warning in Particle::SetDaughter() Already 3 daughters are set!! Check it out!!" << std::endl;
      return;
    }
    fDaughterIndex[fNDaughters++] = value;
  }
  Int_t GetNDaughters() {return fNDaughters;}
  Int_t GetDaughter(Int_t value) {
    if(value<0 || value>fNDaughters-1) {
      std::cout << "Warning in Particle::GetDaughter(Int_t) This particle has " << fNDaughters 
		<< " daughters. The argument must range from 0 to " << fNDaughters-1 << std::endl;
      return -1;
    }
    return fDaughterIndex[value];
  }

  // void SetLastMotherDecayCoor(TLorentzVector fLastMotherDecayCoor);
  TLorentzVector &SetLastMotherDecayCoor(const TLorentzVector &val){return fLastMotherDecayCoor = val;}
  const TLorentzVector &GetLastMotherDecayCoor()const{return fLastMotherDecayCoor;}
  //  void SetLastMotherDecayMom(TLorentzVector fLastMotherDecayMom);
  TLorentzVector &SetLastMotherDecayMom(const TLorentzVector &val){return fLastMotherDecayMom = val;}
  const TLorentzVector &GetLastMotherDecayMom()const{return fLastMotherDecayMom;}

  void SetLastInterTime(Double_t value){fLastInteractionTime = value;}
  Double_t GetLastInterTime()const{return fLastInteractionTime;}
  void SetLastInterNumber(Int_t value){fInteractionNumber = value;}
  Int_t GetLastInterNumber()const{return fInteractionNumber;}
  void IncInter(){++fInteractionNumber;}

  void SetType(Int_t value){fType = value;}
  Int_t GetType()const{return fType;}


};

Double_t S(const TLorentzVector &, const TLorentzVector &);
Double_t T(const TLorentzVector &, const TLorentzVector &);

typedef std::list<Particle> List_t;
typedef std::list<Particle>::iterator LPIT_t;

class ParticleAllocator {
 public:
  void AddParticle(const Particle & particle, List_t & list);
  void FreeListNode(List_t & list, LPIT_t it);
  void FreeList(List_t & list);

 private:
  List_t fFreeNodes;
};

#endif
