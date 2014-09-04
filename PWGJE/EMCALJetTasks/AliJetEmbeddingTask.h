#ifndef ALIJETEMBEDDINGTASK_H
#define ALIJETEMBEDDINGTASK_H

// $Id$

#include "AliJetModelBaseTask.h"

class AliJetEmbeddingTask : public AliJetModelBaseTask {
 public:
  AliJetEmbeddingTask();
  AliJetEmbeddingTask(const char *name); 
  virtual ~AliJetEmbeddingTask();

  void           SetMasslessParticles(Bool_t b) { fMassless        = b ; }
  void           SetNeutralFraction(Double_t f) { fNeutralFraction = f ; }

 protected:
  void           Run();

 private:
  Bool_t         fMassless;               //make particles massless
  Double_t       fNeutralFraction;        //assign charge==0 to fraction of particles

  AliJetEmbeddingTask(const AliJetEmbeddingTask&);            // not implemented
  AliJetEmbeddingTask &operator=(const AliJetEmbeddingTask&); // not implemented

  ClassDef(AliJetEmbeddingTask, 3) // Jet embedding task
};
#endif
