#ifndef ALIJETREADER_H
#define ALIJETREADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// Jet reader base class
// manages the reading of input for jet algorithms
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------
  
#include <TObject.h>
#include <TChain.h>
#include <TArrayI.h>

class TClonesArray;
class AliJetReaderHeader;
class AliESD;
class AliHeader;


class AliJetReader : public TObject 
{
 public: 
  AliJetReader();
  virtual ~AliJetReader();

  // Getters
  virtual TClonesArray *GetMomentumArray() {return fMomentumArray;}
  virtual Int_t GetChainEntries() {return fChain->GetEntries();} 
  virtual AliJetReaderHeader* GetReaderHeader() { return fReaderHeader;}
  virtual AliHeader* GetAliHeader() { return fAliHeader;}
  virtual Int_t GetSignalFlag(Int_t i) const {return fSignalFlag[i];}

  // Setters
  virtual void FillMomentumArray(Int_t) {}
  virtual void SetReaderHeader(AliJetReaderHeader* header) 
    {fReaderHeader = header;}
	  
  // others
  virtual void OpenInputFiles() {}
  void ClearArray();
 
 protected:
  TChain                  *fChain;         // chain for reconstructed tracks
  TChain                  *fChainMC;       // chain for mc information
  TClonesArray            *fMomentumArray; // array of particle momenta
  TClonesArray            *fArrayMC;       // array of mc particles
  AliESD                  *fESD;           // pointer to esd
  AliJetReaderHeader      *fReaderHeader;  // pointer to header
  AliHeader               *fAliHeader;     // pointer to event header
  TArrayI fSignalFlag;   // to flag if a particle comes from pythia or 
                        // from the underlying event

  ClassDef(AliJetReader,1)
};
 
#endif
