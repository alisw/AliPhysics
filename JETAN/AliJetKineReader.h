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
class AliMCEvent;
class TClonesArray;


class AliJetKineReader : public AliJetReader
{
 public: 
  AliJetKineReader();
  virtual ~AliJetKineReader();
  // Setters
  Bool_t  FillMomentumArray();
  void    OpenInputFiles();
  void    SetInputEvent(TObject* esd, TObject* aod, TObject* mc);
  // Fast Simulation
  Float_t SmearMomentum(Int_t ind, Float_t p);
  Bool_t  Efficiency(Float_t pt, Float_t eta, Float_t phi);
  // Others
  TClonesArray*      GetGeneratedJets();
  virtual AliHeader* GetAliHeader() const {return fAliHeader;}
  
 protected:
  AliJetKineReader(const AliJetKineReader& rJetKine);
  AliJetKineReader& operator = (const AliJetKineReader& rkr);

  AliHeader          *fAliHeader;       //! Header
  AliMCEvent         *fMCEvent;  //! Monte Carlo Event Handler
  TClonesArray       *fGenJets;         //! List of generated jets
  ClassDef(AliJetKineReader,1)
};
 
#endif
