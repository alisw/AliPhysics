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

  void           SetGen(AliGenerator *gen)      { fGen = gen       ; }
  void           SetMasslessParticles(Bool_t b) { fMassless = b    ; }
  void           SetChargedOnly(Bool_t b)       { fChargedOnly = b ; }

 protected:
  Bool_t         ExecOnce();
  void           Run();

  AliGenerator  *fGen;                    //generator
  Bool_t         fMassless;               //make particles massless
  Bool_t         fChargedOnly;            //accept only charged particles

  TH1F          *fHistPt;                 //!pT spectrum of generated particles
  TH2F          *fHistEtaPhi;             //!eta-phi of generated particles
  TH1F          *fHistTrials;             //!trials from generator
  TProfile      *fHistXsection;           //!x-section from generator
  TH1           *fHistPtHard;             //!pt hard distribution

 private:
  AliJetEmbeddingFromGenTask(const AliJetEmbeddingFromGenTask&);            // not implemented
  AliJetEmbeddingFromGenTask &operator=(const AliJetEmbeddingFromGenTask&); // not implemented

  ClassDef(AliJetEmbeddingFromGenTask, 4) // Jet embedding task
};
#endif
