#ifndef ALIJETEMBEDDINGFROMGENTASK_H
#define ALIJETEMBEDDINGFROMGENTASK_H

// $Id$

class TClonesArray;
class TProfile;
class AliEMCALGeometry;

#include "AliJetModelBaseTask.h"
class AliGenerator;

class AliJetEmbeddingFromGenTask : public AliJetModelBaseTask {
 public:
  AliJetEmbeddingFromGenTask();
  AliJetEmbeddingFromGenTask(const char *name, Bool_t drawqa);
  virtual ~AliJetEmbeddingFromGenTask();

  void           UserCreateOutputObjects();
  void           FillPythiaHistograms();

  void           SetGen(AliGenerator *gen) { fGen = gen; }

 protected:
  Bool_t         ExecOnce();
  void           Run();

  AliGenerator  *fGen;                    //generator

  TH1F          *fHistTrials;             //!trials from generator
  TProfile      *fHistXsection;           //!x-section from generator
  TH1           *fHistPtHard;             //!pt hard distribution

 private:
  AliJetEmbeddingFromGenTask(const AliJetEmbeddingFromGenTask&);            // not implemented
  AliJetEmbeddingFromGenTask &operator=(const AliJetEmbeddingFromGenTask&); // not implemented

  ClassDef(AliJetEmbeddingFromGenTask, 2) // Jet embedding task
};
#endif
