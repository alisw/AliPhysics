#ifndef ALIJETMODELMERGEBRANCHES_H
#define ALIJETMODELMERGEBRANCHES_H

// $Id$

class TClonesArray;

#include "AliJetModelBaseTask.h"

class AliJetModelMergeBranches : public AliJetModelBaseTask {
 public:
  AliJetModelMergeBranches();
  AliJetModelMergeBranches(const char *name); 
  virtual ~AliJetModelMergeBranches();

  void                   SetTracksMergeName(const char *n)          { fTracksMergeName   = n;    }

 protected:
  Bool_t                 ExecOnce();
  void                   Run();

  void                   MergeTracks();

  TString                fTracksMergeName;        // name of track collection to be merged
  TClonesArray          *fTracksMerge;            //!track collection to be merged

 private:
  AliJetModelMergeBranches(const AliJetModelMergeBranches&);            // not implemented
  AliJetModelMergeBranches &operator=(const AliJetModelMergeBranches&); // not implemented

  ClassDef(AliJetModelMergeBranches, 1) // Jet model merge branches tas
};
#endif
