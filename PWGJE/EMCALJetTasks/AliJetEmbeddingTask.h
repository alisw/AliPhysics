#ifndef ALIJETEMBEDDINGTASK_H
#define ALIJETEMBEDDINGTASK_H

// $Id$

#include "AliJetModelBaseTask.h"

class AliJetEmbeddingTask : public AliJetModelBaseTask {
 public:
  AliJetEmbeddingTask();
  AliJetEmbeddingTask(const char *name); 
  virtual ~AliJetEmbeddingTask();

  void           SetMasslessParticles(Bool_t b) { fMassless = b ; }

 protected:
  void           Run();

  Bool_t         fMassless;               //make particles massless

 private:
  AliJetEmbeddingTask(const AliJetEmbeddingTask&);            // not implemented
  AliJetEmbeddingTask &operator=(const AliJetEmbeddingTask&); // not implemented

  ClassDef(AliJetEmbeddingTask, 2) // Jet embedding task
};
#endif
