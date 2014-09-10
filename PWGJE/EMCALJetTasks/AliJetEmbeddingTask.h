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
  void           SetNeutralMass(Double_t m)     { fNeutralMass     = m ; }

 protected:
  void           Run();

 private:
  Bool_t         fMassless;               //make particles massless
  Double_t       fNeutralFraction;        //assign charge==0 to fraction of particles
  Double_t       fNeutralMass;            //assign this mass to neutral particles

  AliJetEmbeddingTask(const AliJetEmbeddingTask&);            // not implemented
  AliJetEmbeddingTask &operator=(const AliJetEmbeddingTask&); // not implemented

  ClassDef(AliJetEmbeddingTask, 4) // Jet embedding task
};
#endif
