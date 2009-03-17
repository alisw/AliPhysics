/*                                                                           
                                                                            
        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru 
                           November. 2, 2005                                

*/

#include <TMath.h>
#ifndef PARTICLE_PDG
#include "ParticlePDG.h"
#endif
#ifndef PARTICLE_INCLUDED
#include "Particle.h"
#endif

Particle::Particle(ParticlePDG *prop): 
  fPosition(),
  fMomentum(),
  fLastMotherDecayCoor(),
  fLastMotherDecayMom(),
  fParticleProperties(prop),
  fLastInteractionTime(0.),
  fInteractionNumber(0),
  fLastMotherPdg(0),
  fType(0),
  fIndex(-1),
  fMotherIndex(-1),
  fNDaughters(0)
{
  for(Int_t i=0; i<3; i++) fDaughterIndex[i] = -1;
}

Particle::Particle(ParticlePDG *prop, const TLorentzVector &pos, 
		   const TLorentzVector &mom, Double_t lit, Int_t lin, Int_t type):
  fPosition(pos),
  fMomentum(mom),
  fLastMotherDecayCoor(),
  fLastMotherDecayMom(),
  fParticleProperties(prop),
  fLastInteractionTime(lit),
  fInteractionNumber(lin),
  fLastMotherPdg(0),
  fType(type),
  fIndex(-1),
  fMotherIndex(-1),
  fNDaughters(0)
{
  for(Int_t i=0; i<3; i++) fDaughterIndex[i] = -1;
}

Particle::Particle(ParticlePDG *prop, const TLorentzVector &pos, const TLorentzVector &mom,
                   Double_t t, Int_t n, Int_t ty, Int_t motherPdg, const TLorentzVector &mPos, 
		   const TLorentzVector &mMom):
  fPosition(pos),
  fMomentum(mom),
  fLastMotherDecayCoor(mPos),
  fLastMotherDecayMom(mMom),
  fParticleProperties(prop),
  fLastInteractionTime(t),
  fInteractionNumber(n),
  fLastMotherPdg(motherPdg),
  fType(ty),
  fIndex(-1),
  fMotherIndex(-1),
  fNDaughters(0)
{
  for(Int_t i=0; i<3; i++) fDaughterIndex[i] = -1;
}

Int_t Particle::Encoding() const {
  return fParticleProperties->GetPDG();
}

Double_t Particle::TableMass() const {
  return fParticleProperties->GetMass();
}

Double_t Particle::Eta() const {
  if(fMomentum.P() != fMomentum.Pz())
    return 0.5 * TMath::Log((fMomentum.P() + fMomentum.Pz()) / (fMomentum.P()-fMomentum.Pz()));
  else return 1.e30;
}

Double_t Particle::Rapidity() const {
  if (fMomentum.E() != fMomentum.Pz())
    return 0.5 * TMath::Log((fMomentum.E() + fMomentum.Pz()) / (fMomentum.E() - fMomentum.Pz()));
  else return 1.e30;
}

Double_t Particle::Phi() const {
  return TMath::Pi()+TMath::ATan2(-fMomentum.Py(), -fMomentum.Px());
}

Double_t Particle::Theta() const {
  return !fMomentum.Pz() ? TMath::Pi() / 2 : TMath::ACos(fMomentum.Pz() / fMomentum.P());
}

Double_t Particle::Pt() const {
  return TMath::Sqrt(fMomentum.Px() * fMomentum.Px() + fMomentum.Py() * fMomentum.Py());
}

Double_t S(const TLorentzVector &v1, const TLorentzVector &v2) {
  return TMath::Power(v1.T() + v2.T(), 2) - TMath::Power(v1.X() + v2.X(), 2) -
    TMath::Power(v1.Y() + v2.Y(), 2) - TMath::Power(v1.Z() + v2.Z(), 2);
}

Double_t T(const TLorentzVector & v1, const TLorentzVector & v2) {
  return TMath::Power(v1.T() - v2.T(), 2) - TMath::Power(v1.X() - v2.X(), 2) - 
    TMath::Power(v1.Y() - v2.Y(), 2) - TMath::Power(v1.Z() - v2.Z(), 2);
}

void ParticleAllocator::AddParticle(const Particle & p, List_t &list) {
  if(fFreeNodes.empty())
    list.push_back(p);
  else {
    list.splice(list.end(), fFreeNodes, fFreeNodes.begin());
    list.back() = p;
  }
}

void ParticleAllocator::FreeListNode(List_t & list, LPIT_t it) {
  fFreeNodes.splice(fFreeNodes.end(), list, it);      
}

void ParticleAllocator::FreeList(List_t & list) {
  fFreeNodes.splice(fFreeNodes.end(), list);
}
