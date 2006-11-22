#ifndef ALIJETESDREADER_H
#define ALIJETESDREADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// Jet ESD Reader 
// ESD reader for jet analysis
// Author: Mercedes Lopez Noriega (mercedes.lopez.noriega@cern.ch)
//---------------------------------------------------------------------

#include "AliJetReader.h"
class AliJetESDReaderHeader;


class AliJetESDReader : public AliJetReader
{
 public: 
  AliJetESDReader();
  virtual ~AliJetESDReader();

  // Getters
  Float_t GetTrackMass() const {return fMass;}  // returns mass of the track
  Int_t   GetTrackSign() const {return fSign;}  // returns sign of the track

  // Setters
  Bool_t FillMomentumArray(Int_t event); 
  void   OpenInputFiles();
  void   ConnectTree(TTree* tree);
  
 protected:
  Float_t fMass;    // Particle mass
  Int_t   fSign;    // Particle sign

  ClassDef(AliJetESDReader,1)
};
 
#endif
