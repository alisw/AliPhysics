/// \class AliJetEmbeddingFromPYTHIATask
/// \brief Class for Pythia event embedding into data
///
/// The class inherits from AliJetEmbeddingFromAODTask and takes care of the implemetation of the track embedding into the original or a copy of the track array using the method AddTrack (see AliJetModelBaseTask)

#ifndef ALIJETEMBEDDINGFROMPYTHIATASK_H
#define ALIJETEMBEDDINGFROMPYTHIATASK_H

#include "AliJetEmbeddingFromAODTask.h"
#include <TArrayD.h>

template<class T> 
class TParameter;

class TString;
class TH1;
class THashTable;

class AliJetEmbeddingFromPYTHIATask : public AliJetEmbeddingFromAODTask {
 public:
  AliJetEmbeddingFromPYTHIATask();
  AliJetEmbeddingFromPYTHIATask(const char *name, Bool_t drawqa=kFALSE); 
  virtual ~AliJetEmbeddingFromPYTHIATask();

  Bool_t         UserNotify();
  void           UserCreateOutputObjects();

  void           SetPYTHIAPath(const char* p)                      { fPYTHIAPath                                = p ; }
  void           SetPtHardBinScaling(Int_t n, Double_t *scaling)   { new (&fPtHardBinScaling) TArrayD(n, scaling)   ; }
  void           SetAnchorRun(Int_t r)                             { fAnchorRun                                 = r ; }
  void           SetLHC11hAnchorRuns(Bool_t a=kTRUE)               { fLHC11hAnchorRun                           = a ; }
  void           SetFileTable(THashTable *t)                       { fFileTable                                 = t ; }
  void           SetUseAsVetoTable(Bool_t v)                       { fUseAsVetoTable                            = v ; }
  void           SetMinEntriesPerPtHardBin(Int_t r)                { fMinEntriesPerPtHardBin                    = r ; }
  void           SetDebugEmbedding(Bool_t d=kTRUE)                 { fDebugEmbedding                            = d ; }

 protected:
  Bool_t           ExecOnce()               ;// intialize task
  Bool_t           GetNextEntry()           ;// get next entry in current tree
  Int_t            GetRandomPtHardBin()     ;// get a radnom pt hard bin according to fPtHardBinScaling
  TFile           *GetNextFile()            ;// get next file

  TString          fPYTHIAPath              ;// Path of the PYTHIA production
  TArrayD          fPtHardBinScaling        ;// Pt hard bin scaling
  Bool_t           fLHC11hAnchorRun         ;// LHC11h anchor runs
  Int_t            fAnchorRun               ;// Anchor run
  THashTable      *fFileTable               ;// Table of allowed/vetoed files
  Bool_t           fUseAsVetoTable          ;// Use fFileTable as a veto table
  Int_t            fMinEntriesPerPtHardBin  ;// Minimum number of embedded events before changing pt hard bin, if < 0 change pt hard bin only when reach eof 
  Int_t            fCurrentPtHardBin        ;//!Pt hard bin of the current open file
  TParameter<int> *fPtHardBinParam          ;//!Pt hard bin param
  Int_t            fPtHardBinCount          ;//!Number of event embedded from the current pt hard bin
  Bool_t           fDebugEmbedding          ;// Debug embedding by embedding files in order. Do _not_ use on the grid! (It will embed the same files for each job)

  TH1             *fHistPtHardBins          ;//!Embeded pt hard bin distribution

 private:
  AliJetEmbeddingFromPYTHIATask(const AliJetEmbeddingFromPYTHIATask&);            // not implemented
  AliJetEmbeddingFromPYTHIATask &operator=(const AliJetEmbeddingFromPYTHIATask&); // not implemented

  ClassDef(AliJetEmbeddingFromPYTHIATask, 4) // Jet embedding from PYTHIA task
};
#endif
