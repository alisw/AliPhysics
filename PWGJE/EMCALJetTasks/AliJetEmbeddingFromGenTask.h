#ifndef ALIJETEMBEDDINGFROMGENTASK_H
#define ALIJETEMBEDDINGFROMGENTASK_H

// $Id: AliJetEmbeddingFromGenTask.h 57324 2012-06-21 04:33:52Z loizides $

class TClonesArray;
class AliEMCALGeometry;

#include "AliJetModelBaseTask.h"
class AliGenerator;

class AliJetEmbeddingFromGenTask : public AliJetModelBaseTask {
 public:
  AliJetEmbeddingFromGenTask();
  AliJetEmbeddingFromGenTask(const char *name); 
  virtual ~AliJetEmbeddingFromGenTask();

  void           SetGen(AliGenerator *gen) { fGen = gen; }

 protected:
  void           ExecOnce();
  void           Run();

  AliGenerator  *fGen;    //generator

 private:
  AliJetEmbeddingFromGenTask(const AliJetEmbeddingFromGenTask&);            // not implemented
  AliJetEmbeddingFromGenTask &operator=(const AliJetEmbeddingFromGenTask&); // not implemented

  ClassDef(AliJetEmbeddingFromGenTask, 2) // Jet embedding task
};
#endif
