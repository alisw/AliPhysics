#ifndef ALIJETREADER_H
#define ALIJETREADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
// Jet reader base class
// manages the reading of input for jet algorithms
// Author: jgcn@mda.cinvestav.mx
  
#include <TObject.h>
#include <TChain.h>
#include <TArrayI.h>
class TTree;
class TClonesArray;
class AliJetReaderHeader;
class AliESD;
class AliJet;

class AliJetReader : public TObject 
{
 public: 
  AliJetReader();
  virtual ~AliJetReader();

  // Getters
  virtual TClonesArray *GetMomentumArray() {return fMomentumArray;}
  virtual Int_t GetChainEntries() {return fChain->GetEntries();} 
  virtual AliJetReaderHeader* GetReaderHeader() { return fReaderHeader;}
  virtual Int_t GetSignalFlag(Int_t i) const {return fSignalFlag[i];}
  virtual Int_t GetCutFlag(Int_t i) const {return fCutFlag[i];}
  
  // Setters
  virtual Bool_t FillMomentumArray(Int_t) {return kTRUE;}
  virtual void   SetReaderHeader(AliJetReaderHeader* header) 
    {fReaderHeader = header;}
	  
  // Others
  virtual void   OpenInputFiles() {}
  virtual void   ConnectTree(TTree* /*tree*/) {}
  virtual Bool_t GetGenJets(AliJet* /*genJets*/) {return kFALSE;}
  
  void ClearArray();
 
 protected:
  AliJetReader(const AliJetReader& rJetReader);
  AliJetReader& operator = (const AliJetReader& rhsr);

  TChain                  *fChain;         // chain for reconstructed tracks
  TChain                  *fChainMC;       // chain for mc information
  TClonesArray            *fMomentumArray; // array of particle momenta
  TClonesArray            *fArrayMC;       // array of mc particles
  AliESD                  *fESD;           // pointer to esd
  AliJetReaderHeader      *fReaderHeader;  // pointer to header
  TArrayI fSignalFlag;   // to flag if a particle comes from pythia or 
                         // from the underlying event
  TArrayI fCutFlag;      // to flag if a particle passed the pt cut or not

  ClassDef(AliJetReader,1)
};
 
#endif
