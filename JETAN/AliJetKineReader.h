#ifndef ALIJETKINEREADER_H
#define ALIJETKINEREADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
// Jet Kine Reader 
// MC Kinematics reader for jet analysis
// Author: Andreas Morsch (andreas.morsch@cern.ch)

#include "AliJetReader.h"

class AliRunLoader;
class AliHeader;

class AliJetKineReader : public AliJetReader
{
 public: 
  AliJetKineReader();
  virtual ~AliJetKineReader();

  // Getters
  Float_t GetParticleMass() const {return fMass;}        // returns mass of the Track
  Int_t   GetParticlePdgCode() const {return fPdgC;}     // returns Pdg code
  // Setters
  Bool_t  FillMomentumArray(Int_t event);
  void    OpenInputFiles();
  // Fast Simulation
  Float_t SmearMomentum(Int_t ind, Float_t p);
  Bool_t  Efficiency(Float_t pt, Float_t eta, Float_t phi);
  // Others
  virtual Bool_t GetGenJets(AliJet* /*genJets*/);
  virtual AliHeader* GetAliHeader() {return fAliHeader;}
  
 protected:
  AliJetKineReader(const AliJetKineReader& rJetKine);
  AliJetKineReader& operator = (const AliJetKineReader& rkr);

  AliRunLoader *fRunLoader;       //! Pointer to the run loader
  AliHeader    *fAliHeader;       //! Header
  
  Float_t fMass;                  //!Particle mass
  Int_t   fPdgC;                  //!Pdg code
 
  ClassDef(AliJetKineReader,1)
};
 
#endif
