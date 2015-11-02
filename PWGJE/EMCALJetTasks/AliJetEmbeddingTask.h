/// \class AliJetEmbeddingTask
/// \brief Class for track embedding into an event
///
/// The class inherits from AliJetModelBaseTask and takes care of the implemetation of the track embedding into the original or a copy of the track array using the method AddTrack (see AliJetModelBaseTask)
/// Several choices on the track mass are possible: 
/// 1) pion mass
/// 2) massless
/// 3) a value set by the user
/// 4) a random value from a user-defined distribution
///
/// \author S.Aiola, 
/// \author C.Loizides
/// \date

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
  void           SetMass(Double_t mass)         { fMass = mass ; }
  void           SetNeutralFraction(Double_t f) { fNeutralFraction = f ; }
  void           SetNeutralMass(Double_t m)     { fNeutralMass     = m ; }
  void           SetMassDistribution(TH1F *hM);
  void           SetMassDistributionFromFile(TString filename, TString histoname);
  void           SetpTDistributionFromFile(TString filename, TString histoname);
  void           SetMassAndPtDistributionFromFile(TString filenameM, TString filenamepT, TString histonameM, TString histonamepT);
  
  
 protected:
  void           Run();

 private:
  Bool_t         fMassless;               ///< make particles massless
  Bool_t         fMassFromDistr;          ///< draw the particle mass from fHMassDistrib
  Double_t       fNeutralFraction;        ///< assign charge==0 to fraction of particles
  Double_t       fNeutralMass;            ///< assign this mass to neutral particles
  Double_t       fMass;                   ///< assign this mass to particle
  TH1F*          fHMassDistrib;           ///< shape of mass distribution of embedded tracks
  
  AliJetEmbeddingTask(const AliJetEmbeddingTask&);            // not implemented
  AliJetEmbeddingTask &operator=(const AliJetEmbeddingTask&); // not implemented
  
  /// \cond CLASSIMP
  ClassDef(AliJetEmbeddingTask, 5) /// Jet embedding task
  /// \endcond
};
#endif
