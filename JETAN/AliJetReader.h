#ifndef ALIJETREADER_H
#define ALIJETREADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
// Jet reader base class
// manages the reading of input for jet algorithms
// Authors: jgcn@mda.cinvestav.mx
//          Magali Estienne <magali.estienne@IReS.in2p3.fr>  

#include <TObject.h>
#include <TChain.h>
#include <TArrayI.h>
#ifndef ROOT_TTask
#include "TTask.h"
#endif

class TTree;
class TTask;
class TClonesArray;
class AliJetReaderHeader;
class AliJetUnitArray;
class AliJetHadronCorrectionv1;
class AliJet;

class AliJetReader : public TObject 
{
 public: 
  AliJetReader();
  virtual ~AliJetReader();

  // Getters
  virtual TClonesArray *GetMomentumArray() {return fMomentumArray;}

  virtual AliJetUnitArray     *GetUnitArray() const {return fUnitArray;}  
  virtual AliJetUnitArray     *GetUnitArrayNoCuts() const {return fUnitArrayNoCuts;}
  
  virtual AliJetReaderHeader* GetReaderHeader() { return fReaderHeader;}
  virtual Int_t GetSignalFlag(Int_t i) const {return fSignalFlag[i];}
  virtual Int_t GetCutFlag(Int_t i)    const {return fCutFlag[i];}
  virtual Int_t GetArrayInitialised() {return fArrayInitialised;}
  
  // Setters
  virtual Bool_t FillMomentumArray(Int_t) {return kTRUE;}
  virtual void   FillUnitArrayFromTPCTracks(Int_t) {}     // temporarily not used
  virtual void   FillUnitArrayFromEMCALHits() {}          // temporarily not used
  virtual void   FillUnitArrayFromEMCALDigits(Int_t) {}   // temporarily not used
  virtual void   FillUnitArrayFromEMCALClusters(Int_t) {} // temporarily not used
  virtual void   InitUnitArray() {}
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
  TClonesArray            *fMomentumArray;    // array of particle momenta
  TClonesArray            *fArrayMC;          // array of mc particles
  TTask                   *fFillUnitArray;    // task list for filling the UnitArray
  AliJetReaderHeader      *fReaderHeader;     // pointer to header
  TArrayI                  fSignalFlag;       // to flag if a particle comes from pythia or
                                              // from the underlying event
  TArrayI                  fCutFlag;          // to flag if a particle passed the pt cut or not
  AliJetUnitArray         *fUnitArray;        // array of digit position and energy 
  AliJetUnitArray         *fUnitArrayNoCuts;  // array of digit position and energy 
  Bool_t                   fArrayInitialised; // To check that array of units is initialised  
  ClassDef(AliJetReader,1)
};
 
#endif
