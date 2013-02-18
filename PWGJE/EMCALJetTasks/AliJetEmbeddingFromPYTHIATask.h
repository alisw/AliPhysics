#ifndef ALIJETEMBEDDINGFROMPYTHIATASK_H
#define ALIJETEMBEDDINGFROMPYTHIATASK_H

// $Id: AliJetEmbeddingFromPYTHIATask.h  $

#include "AliJetEmbeddingFromAODTask.h"
#include <TArrayD.h>

class TString;

class AliJetEmbeddingFromPYTHIATask : public AliJetEmbeddingFromAODTask {
 public:
  AliJetEmbeddingFromPYTHIATask();
  AliJetEmbeddingFromPYTHIATask(const char *name); 
  virtual ~AliJetEmbeddingFromPYTHIATask();

  Bool_t         UserNotify();

  void           SetPYTHIAPath(const char* p)                      { fPYTHIAPath                                = p ; }
  void           SetPtHardBinScaling(Int_t n, Double_t *scaling)   { new (&fPtHardBinScaling) TArrayD(n, scaling)   ; }
  void           SetAnchorRun(Int_t r)                             { fAnchorRun                                 = r ; }
  void           SetLHC11hAnchorRuns(Bool_t a=kTRUE)               { fLHC11hAnchorRun                           = a ; }

 protected:
  Bool_t         ExecOnce()           ;// intialize task
  Bool_t         GetNextEntry()       ;// get next entry in current tree
  Int_t          GetRandomPtHardBin() ;// get a radnom pt hard bin according to fPtHardBinScaling
  TString        GetNextFileName()    ;// get next file name

  TString        fPYTHIAPath          ;// Path of the PYTHIA production
  TArrayD        fPtHardBinScaling    ;// Pt hard bin scaling
  Bool_t         fLHC11hAnchorRun     ;// LHC11h anchor runs
  Int_t          fAnchorRun           ;// Anchor run
  Int_t          fCurrentPtHardBin    ;//!Pt hard bin of the current open file


 private:
  AliJetEmbeddingFromPYTHIATask(const AliJetEmbeddingFromPYTHIATask&);            // not implemented
  AliJetEmbeddingFromPYTHIATask &operator=(const AliJetEmbeddingFromPYTHIATask&); // not implemented

  ClassDef(AliJetEmbeddingFromPYTHIATask, 1) // Jet embedding from PYTHIA task
};
#endif
