#ifndef ALIJETEMBEDDINGTASK_H
#define ALIJETEMBEDDINGTASK_H

// $Id$

#include "AliJetModelBaseTask.h"

class AliJetEmbeddingTask : public AliJetModelBaseTask {
 public:
  AliJetEmbeddingTask();
  AliJetEmbeddingTask(const char *name); 
  virtual ~AliJetEmbeddingTask();

 protected:
  void           Run();

 private:
  AliJetEmbeddingTask(const AliJetEmbeddingTask&);            // not implemented
  AliJetEmbeddingTask &operator=(const AliJetEmbeddingTask&); // not implemented

  ClassDef(AliJetEmbeddingTask, 2) // Jet embedding task
};
#endif
