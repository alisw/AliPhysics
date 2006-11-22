#ifndef ALIJETESDMCREADER_H
#define ALIJETESDMCREADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
// Jet ESD Reader 
// ESD reader for jet analysis (it reads the esd and the MC trees)
// Author: Mercedes Lopez Noriega (mercedes.lopez.noriega@cern.ch)

#include "AliJetReader.h"
class AliJetESDReaderHeader;
class AliHeader;



class AliJetESDmcReader : public AliJetReader
{
 public: 
  AliJetESDmcReader();
  virtual ~AliJetESDmcReader();

  // Getters
  Float_t GetTrackMass() const {return fMass;}  // returns mass of the track
  Int_t   GetTrackSign() const {return fSign;}  // returns sign of the track

  // Setters
  Bool_t FillMomentumArray(Int_t event); 
  void   OpenInputFiles();
   
 protected:
  AliHeader* fAliHeader; //! Event header
  Float_t    fMass;      //! Particle mass
  Int_t      fSign;      //! Particle sign

  ClassDef(AliJetESDmcReader,1)
};
 
#endif
