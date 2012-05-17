#ifndef ALIJETEMBEDDINGTASK_H
#define ALIJETEMBEDDINGTASK_H

// $Id$

class TClonesArray;
class AliEMCALGeometry;

#include "AliJetModelBaseTask.h"

class AliJetEmbeddingTask : public AliJetModelBaseTask {
 public:
  AliJetEmbeddingTask();
  AliJetEmbeddingTask(const char *name); 
  virtual ~AliJetEmbeddingTask();

  virtual void           UserExec(Option_t* /*option*/);

 protected:

  virtual void           Run();                 // do embedding

 private:
  AliJetEmbeddingTask(const AliJetEmbeddingTask&);            // not implemented
  AliJetEmbeddingTask &operator=(const AliJetEmbeddingTask&); // not implemented

  ClassDef(AliJetEmbeddingTask, 2) // Jet embedding task
};
#endif
