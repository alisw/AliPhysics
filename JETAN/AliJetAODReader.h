#ifndef ALIJETAODREADER_H
#define ALIJETAODREADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// Jet AOD Reader
// AOD reader for jet analysis
// Author: Davide Perrino (davide.perrino@cern.ch)
//---------------------------------------------------------------------

#include "AliJetReader.h"
class AliJetAODReaderHeader;
class AliJetReaderHeader;
class AliAODEvent;
class TRefArray;

class AliJetAODReader : public AliJetReader
{
 public: 
  AliJetAODReader();
  virtual ~AliJetAODReader();

  TRefArray*   GetReferences() const {return fRef;}

  Bool_t FillMomentumArray(); 
  void   OpenInputFiles();
  void   ConnectTree(TTree* tree, TObject* data);
  void   SetInputEvent(TObject* /*esd*/, TObject* aod, TObject* /*mc*/) {fAOD = (AliAODEvent*) aod;}
 private:
  AliJetAODReader(const AliJetAODReader &det);
  AliJetAODReader &operator=(const AliJetAODReader &det);

 private:
  TChain                     *fChain;  //! chain for reconstructed tracks
  AliAODEvent                *fAOD;    //! pointer to aod
  TRefArray                  *fRef;    // pointer to array of references to tracks
  Int_t                       fDebug;  // Debug option
  Int_t                       fOpt;    // Detector to be used for jet reconstruction
  ClassDef(AliJetAODReader,1)
};
 
#endif
