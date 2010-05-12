//                                                                            
//                                                                            
//        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
//      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru 
//                           November. 2, 2005                                
//
//   Implementation of class Particle
//   Contains particle PDG, 4-coordinates and 4-momentum vectors, history information

#include <TMath.h>
#include "ParticlePDG.h"
#include "Particle.h"

//___________________________________________________________________
Particle::Particle(ParticlePDG *prop): 
  fPosition(),
  fMomentum(),
  fLastMotherDecayCoor(),
  fLastMotherDecayMom(),
  fParticleProperties(prop),
  fLastInteractionTime(0.),
  fInteractionNumber(0),
  fPythiaStatusCode(-1),
  fLastMotherPdg(0),
  fType(0),
  fIndex(-1),
  fMotherIndex(-1),
  fNDaughters(0),
  fFirstDaughterIndex(-1),
  fLastDaughterIndex(-1),
  fDecayed(kFALSE)
{
//
// constructor
//
}

//___________________________________________________________________
Particle::Particle(ParticlePDG *prop, const TLorentzVector &pos, 
		   const TLorentzVector &mom, Double_t lit, Int_t lin, Int_t type):
  fPosition(pos),
  fMomentum(mom),
  fLastMotherDecayCoor(),
  fLastMotherDecayMom(),
  fParticleProperties(prop),
  fLastInteractionTime(lit),
  fInteractionNumber(lin),
  fPythiaStatusCode(-1),
  fLastMotherPdg(0),
  fType(type),
  fIndex(-1),
  fMotherIndex(-1),
  fNDaughters(0),
  fFirstDaughterIndex(-1),
  fLastDaughterIndex(-1),
  fDecayed(kFALSE)
{
  //
  // constructor
  //
}

//___________________________________________________________________
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
  fPythiaStatusCode(-1),
  fLastMotherPdg(motherPdg),
  fType(ty),
  fIndex(-1),
  fMotherIndex(-1),
  fNDaughters(0),
  fFirstDaughterIndex(-1),
  fLastDaughterIndex(-1),
  fDecayed(kFALSE)
{
  //
  // constructor
  //
}

//___________________________________________________________________
Particle::Particle(const Particle& copy) :
  fPosition(copy.Pos()),
  fMomentum(copy.Mom()),
  fLastMotherDecayCoor(copy.GetLastMotherDecayCoor()),
  fLastMotherDecayMom(copy.GetLastMotherDecayMom()),
  fParticleProperties(copy.Def()),
  fLastInteractionTime(copy.GetLastInterTime()),
  fInteractionNumber(copy.GetLastInterNumber()),
  fPythiaStatusCode(copy.GetPythiaStatusCode()),
  fLastMotherPdg(copy.GetLastMotherPdg()),
  fType(copy.GetType()),
  fIndex(copy.GetIndex()),
  fMotherIndex(copy.GetMother()),
  fNDaughters(copy.GetNDaughters()),
  fFirstDaughterIndex(copy.GetFirstDaughterIndex()),
  fLastDaughterIndex(copy.GetLastDaughterIndex()),
  fDecayed(copy.GetDecayed())
{
  //
  // copy constructor
  //
}

//___________________________________________________________________
Particle & Particle::operator=(const Particle& /*copy*/) {
  //
  // assignment operator
  //
  return *this;
}

//___________________________________________________________________
Int_t Particle::Encoding() const {
  //
  // particle code
  //
  return fParticleProperties->GetPDG();
}

//___________________________________________________________________
Double_t Particle::TableMass() const {
  //
  // particle mass
  //
  return fParticleProperties->GetMass();
}

//___________________________________________________________________
Double_t Particle::Eta() const {
  //
  // pseudo-rapidity
  //
  if(fMomentum.P() != fMomentum.Pz())
    return 0.5 * TMath::Log((fMomentum.P() + fMomentum.Pz()) / (fMomentum.P()-fMomentum.Pz()));
  else return 1.e30;
}

//___________________________________________________________________
Double_t Particle::Rapidity() const {
  //
  // rapidity
  //
  if (fMomentum.E() != fMomentum.Pz())
    return 0.5 * TMath::Log((fMomentum.E() + fMomentum.Pz()) / (fMomentum.E() - fMomentum.Pz()));
  else return 1.e30;
}

//___________________________________________________________________
Double_t Particle::Phi() const {
  //
  // azimuthal angle
  //
  return TMath::Pi()+TMath::ATan2(-fMomentum.Py(), -fMomentum.Px());
}

//___________________________________________________________________
Double_t Particle::Theta() const {
  //
  // polar angle
  //
  return !fMomentum.Pz() ? TMath::Pi() / 2 : TMath::ACos(fMomentum.Pz() / fMomentum.P());
}

//___________________________________________________________________
Double_t Particle::Pt() const {
  //
  // pt
  //
  return TMath::Sqrt(fMomentum.Px() * fMomentum.Px() + fMomentum.Py() * fMomentum.Py());
}

//___________________________________________________________________
Double_t S(const TLorentzVector &v1, const TLorentzVector &v2) {
  //
  // Mandelstam s
  //
  return TMath::Power(v1.T() + v2.T(), 2) - TMath::Power(v1.X() + v2.X(), 2) -
    TMath::Power(v1.Y() + v2.Y(), 2) - TMath::Power(v1.Z() + v2.Z(), 2);
}

//___________________________________________________________________
Double_t T(const TLorentzVector & v1, const TLorentzVector & v2) {
  //
  //  Mandelstam t
  //
  return TMath::Power(v1.T() - v2.T(), 2) - TMath::Power(v1.X() - v2.X(), 2) - 
    TMath::Power(v1.Y() - v2.Y(), 2) - TMath::Power(v1.Z() - v2.Z(), 2);
}

//___________________________________________________________________
void ParticleAllocator::AddParticle(const Particle & p, List_t &list) {
  //
  // Add a particle to the list
  //
  if(fFreeNodes.empty())
    list.push_back(p);
  else {
    list.splice(list.end(), fFreeNodes, fFreeNodes.begin());
    list.back() = p;
  }
}

//___________________________________________________________________
void ParticleAllocator::FreeListNode(List_t & list, LPIT_t it) {
  //
  // remove a particle from list
  //
  fFreeNodes.splice(fFreeNodes.end(), list, it);      
}

//___________________________________________________________________
void ParticleAllocator::FreeList(List_t & list) {
  //
  // free all list
  //
  fFreeNodes.splice(fFreeNodes.end(), list);
}
