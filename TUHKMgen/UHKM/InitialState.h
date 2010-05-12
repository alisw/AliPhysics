//
//  Ludmila Malinina  malinina@lav01.sinp.msu.ru,   SINP MSU/Moscow and JINR/Dubna
//  Ionut Arsene  i.c.arsene@fys.uio.no,            Oslo University and ISS-Bucharest
//  Date        : 2007/05/30
//  Updated     : 2008/08/11
//
// Virtual class for the initial state classes
// Include here common methods, but always declare them as virtual

#ifndef INITIALSTATE_H
#define INITIALSTATE_H

#include "Particle.h"
#include "DatabasePDG.h"

class InitialState {
 public:
  InitialState() : fDatabase(new DatabasePDG()) {};
  virtual ~InitialState() {
    if(fDatabase)
      delete fDatabase;
  };

  virtual void SetPDGParticleFilename(Char_t *filename) {fDatabase->SetParticleFilename(filename);}
  virtual void SetPDGDecayFilename(Char_t *filename) {fDatabase->SetDecayFilename(filename);}
  virtual void SetUseCharmParticles(Bool_t flag) {fDatabase->SetUseCharmParticles(flag);}
  virtual void SetMinimumWidth(Double_t value) {fDatabase->SetMinimumWidth(value);}
  virtual void SetMaximumWidth(Double_t value) {fDatabase->SetMaximumWidth(value);}
  virtual void SetMinimumMass(Double_t value) {fDatabase->SetMinimumMass(value);}
  virtual void SetMaximumMass(Double_t value) {fDatabase->SetMaximumMass(value);}
  virtual void SetWidthRange(Double_t min, Double_t max) {fDatabase->SetWidthRange(min, max);}
  virtual void SetMassRange(Double_t min, Double_t max) {fDatabase->SetMassRange(min, max);}
  virtual void LoadPDGInfo() {fDatabase->LoadData();}
  virtual void SetPDGParticleStable(Int_t pdg, Bool_t value) {fDatabase->SetStable(pdg, value);}
  virtual Bool_t GetPDGParticleStableStatus(Int_t pdg) {return fDatabase->GetStableStatus(pdg);}
  virtual Bool_t GetUseCharmParticles() {return fDatabase->GetUseCharmParticles();}
  virtual Double_t GetMinimumWidth() {return fDatabase->GetMinimumWidth();}
  virtual Double_t GetMaximumWidth() {return fDatabase->GetMaximumWidth();}
  virtual Double_t GetMinimumMass() {return fDatabase->GetMinimumMass();}
  virtual Double_t GetMaximumMass() {return fDatabase->GetMaximumMass();}

  virtual DatabasePDG* PDGInfo() const {return fDatabase;}
 
  virtual void Initialize(List_t &source, ParticleAllocator &allocator) = 0;
  virtual Bool_t ReadParams() = 0;
  virtual Bool_t MultIni() = 0;
  virtual Bool_t RunDecays() = 0;
  virtual Int_t GetNev() = 0;
  virtual Double_t GetWeakDecayLimit() = 0;
    
  //  virtual void Evolve(List_t &source, List_t &secondaries, ParticleAllocator &allocator, Double_t weakDecayLimit);
  virtual void Evolve(List_t &secondaries, ParticleAllocator &allocator, Double_t weakDecayLimit);
 protected:
   DatabasePDG *fDatabase;        // PDG database
 private:
   InitialState(const InitialState&);
   InitialState& operator=(const InitialState&);
};

#endif
