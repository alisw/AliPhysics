#ifndef ALIEVENTPOOLSPARSE_H
#define ALIEVENTPOOLSPARSE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Event pool based on THnSparseI
// This class is needed by the AnalysisManager to steer a mixing analysis.
// Author: Peter Hristov
// Peter.Hristov@cern.ch

#include <THnSparse.h>
#include <TEntryList.h>
#include "AliVEventPool.h"

class TChain;
class TTreeFormula;

//_____________________________________________________________________________
class AliEventPoolSparse : public AliVEventPool {
 public:

  AliEventPoolSparse();
  AliEventPoolSparse(const char* name, const char* title, TChain * tagchain, Int_t dim,
		     const char ** vars, const Int_t* nbins, const Double_t* xmin = 0,
		     const Double_t* xmax = 0, Int_t chunksize = 1024 * 16);

  virtual ~AliEventPoolSparse();

  // Interface from AiVEventPool, to be overloaded
  virtual TChain* GetNextChain();
  virtual void  GetCurrentBin(Float_t* xbin);
  virtual Int_t GetDimension(){return fN;}
  virtual void  Init();

  TEntryList * GetNextPool(Int_t i) {
    // Returns the array associated with bin "i"
    return fPool>0 ? fPool[i] : 0x0;
  }

  TEntryList * GetEvents(const Double_t * x) {
    // Returns the array associated with the bin
    // that corresponds to vector "x"
    Int_t bin = fHnSparseI.GetBin(x,kFALSE);
    return fPool>0 ? fPool[bin] : 0x0;
  }

  void SetTagChain(TChain * chain){
    // Input tag chain
    fTagChain = chain;
  }

  TChain * GetTagChain() const {
    // Return the input tag chain
    return fTagChain;
  }

  // Cuts
  void SetRunCut(const char * cut);
  void SetLHCCut(const char * cut);
  void SetDetCut(const char * cut);
  void SetEventCut(const char * cut);

  TTreeFormula ** GetPoolVars() const {return fVars;}
  TTreeFormula * GetRunCut() const {return fRunCut;}
  TTreeFormula * GetLHCCut() const {return fLHCCut;}
  TTreeFormula * GetDetCut() const {return fDetCut;}
  TTreeFormula * GetEventCut() const {return fEvCut;}
  Int_t BinNumber() const {return fBinNumber;}
	  
 protected:

  void Set(Int_t n);


 private:

  AliEventPoolSparse(const AliEventPoolSparse& source); // Not implemented
  AliEventPoolSparse& operator = (const AliEventPoolSparse& source); // Not implemented

  THnSparseI fHnSparseI; // Sparse histogram to 
  Int_t fChunkSize;      //! Cached chunk size since the getter is protected
  Int_t fN;              // Size of the array fPool
  TEntryList ** fPool;   // Arrays of pointers to the TEntryList containing the event IDs
  Int_t fCurrentBin;     //! Current bin
  TChain * fTagChain;    //! Input chain of tags

  TTreeFormula ** fVars; // Array of variables used to create the pools 
  TTreeFormula * fRunCut;// Run selection
  TTreeFormula * fLHCCut;// LNC-based selection
  TTreeFormula * fDetCut;// Detector-based selection
  TTreeFormula * fEvCut; // Event-based selection
  Int_t fBinNumber;      // Current bin
  
  ClassDef(AliEventPoolSparse,1)  // 
};

#endif


