// DatabasePDG stores and handles PDG information
// The PDG particle definitions and decay channels are read
// in the begining from ASCII files

/*
  Copyright   : The FASTMC and SPHMC Collaboration
  Author      : Ionut Cristian Arsene 
  Affiliation : Oslo University, Norway & Institute for Space Sciences, Bucharest, Romania
  e-mail      : i.c.arsene@fys.uio.no
  Date        : 2007/05/30

  This class is using the particle and decay lists provided by the 
  THERMINATOR (Computer Physics Communications 174 669 (2006)) and
  SHARE (Computer Physics Communications 167 229 (2005)) collaborations.
*/

#ifndef DATABASEPDG_H
#define DATABASEPDG_H

#include <Rtypes.h>
#ifndef PARTICLE_PDG
#include "ParticlePDG.h"
#endif

const Int_t kMaxParticles = 1000;

class DatabasePDG {
 public:
  DatabasePDG();
  ~DatabasePDG();

  // Load the particle PDG information from the particle and decay files
  Bool_t LoadData();                        
  
  // Set particle and decay filenames
  void SetParticleFilename(Char_t *filename);
  void SetDecayFilename(Char_t *filename);
  // Set criteria for using particles. Those particle which do not match
  // these criteria will be flagged with FALSE in the fStatus array.
  void SetUseCharmParticles(Bool_t flag);
  void SetMinimumWidth(Double_t value);
  void SetMaximumWidth(Double_t value);
  void SetMinimumMass(Double_t value);
  void SetMaximumMass(Double_t value);
  void SetWidthRange(Double_t min, Double_t max);
  void SetMassRange(Double_t min, Double_t max);
  
  // Read a list of pdg codes from a specified file. The corresponding particles
  // will be flagged as good particles. If the exclusive flag is TRUE than
  // only this criteria will be used in selecting particles and, in consequence,
  // all the other particles will be flagged as NOT good. If the exclusive flag
  // is FALSE than we will take into account all the previous applied criterias
  // and we will flag as good only particles in this list which match also the mass, width and
  // charmness criteria.
  // Note: In order for the exclusive=FALSE to be effective, this function must be called after
  // calling all the width, mass and charmness criteria functions.
  void UseThisListOfParticles(Char_t *filename, Bool_t exclusive = kTRUE);
  
  const Char_t* GetParticleFilename() const {return fParticleFilename;}
  const Char_t* GetDecayFilename() const {return fDecayFilename;}
  Int_t GetNParticles(Bool_t all = kFALSE) const;      // true - no. of all particles; false - no. of good status particles
  ParticlePDG* GetPDGParticleByIndex(Int_t index) const;
  Bool_t GetPDGParticleStatusByIndex(Int_t index) const;
  ParticlePDG* GetPDGParticle(Int_t pdg) const;
  Bool_t GetPDGParticleStatus(Int_t pdg) const;
  ParticlePDG* GetPDGParticle(Char_t *name) const;
  Bool_t GetPDGParticleStatus(Char_t *name) const; 
  Bool_t GetUseCharmParticles() const {return fUseCharmParticles;};
  Double_t GetMinimumWidth() const {return fMinimumWidth;};
  Double_t GetMaximumWidth() const {return fMaximumWidth;};
  Double_t GetMinimumMass() const {return fMinimumMass;};
  Double_t GetMaximumMass() const {return fMaximumMass;};
  void DumpData(Bool_t dumpAll = kFALSE) const; // print the PDG information in the console
  Int_t CheckImpossibleDecays(Bool_t dump = kFALSE) const;   // print all impossible decays included in the database
  Bool_t IsChannelAllowed(DecayChannel *channel, Double_t motherMass) const;
  Int_t GetNAllowedChannels(ParticlePDG *particle, Double_t motherMass) const;
  void SetStable(Int_t pdg, Bool_t value) {GetPDGParticle(pdg)->SetStable(value);}
  Bool_t GetStableStatus(Int_t pdg) const {return GetPDGParticle(pdg)->GetStableStatus();}

 private:
  DatabasePDG(const DatabasePDG&);
  DatabasePDG& operator=(const DatabasePDG&);

  Int_t fNParticles;                        // no. of particles in database
  ParticlePDG *fParticles[kMaxParticles];   // array of particle pointers
  Bool_t fStatus[kMaxParticles];            // status of each particle
  Char_t fParticleFilename[256];            // particle list filename
  Char_t fDecayFilename[256];               // decay channels filename
  Bool_t fUseCharmParticles;                // flag for using (or not) charm particles
  Double_t fMinimumWidth;                   // minimum allowed width for resonances
  Double_t fMaximumWidth;                   // maximum allowed width for resonances
  Double_t fMinimumMass;                    // minimum allowed mass for resonances
  Double_t fMaximumMass;                    // maximum allowed mass for resonances

  Bool_t LoadParticles();
  Bool_t LoadDecays();
  void SortParticles();                     // put the good status particles at the beggining of the list
};

#endif
