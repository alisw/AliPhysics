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

  void           UserCreateOutputObjects();
  void           SetMasslessParticles(Bool_t b) { fMassless        = b ; }
  void           SetMass(Double_t mass)         { fMass = mass ; }
  void           SetNeutralFraction(Double_t f) { fNeutralFraction = f ; }
  void           SetNeutralMass(Double_t m)     { fNeutralMass     = m ; }
  void           SetNamesForTree(TString path, TString treename, TString branchname) { fPathTreeinputFile = path; fTreeinputName = treename; fBranchJDetName = branchname; }
  void           SetPathAndNameInputpTDistribution(TString path, TString name) { fPathpTinputFile = path; fpTinputName = name;}
  void           SetPathAndNameInputMDistribution(TString path, TString name) { fPathMinputFile = path; fMinputName = name;}
  
 protected:
  void           Run();
  
  void           SetMassDistribution(TH1F *hM);
  void           SetMassDistributionFromFile(TString filename, TString histoname);
  void           SetpTDistributionFromFile(TString filename, TString histoname);
  void           SetMassAndPtDistributionFromFile(TString filenameM, TString filenamepT, TString histonameM, TString histonamepT);
  void           SetTree(TTree *tree);
  void           SetTreeFromFile(TString filenameM, TString treename);
 
 private:
  Bool_t         fMassless;               ///< make particles massless
  Bool_t         fMassFromDistr;          ///< draw the particle mass from fHMassDistrib
  Double_t       fNeutralFraction;        ///< assign charge==0 to fraction of particles
  Double_t       fNeutralMass;            ///< assign this mass to neutral particles
  Double_t       fMass;                   ///< assign this mass to particle
  TH1F*          fHMassDistrib;           //!<! shape of mass distribution of embedded tracks
  TString        fPathMinputFile;         ///< path to the file where the external input Mass distribution is (can be from alien)
  TString        fPathpTinputFile;        ///< path to the file where the external input pT distribution is (can be from alien)
  TString        fMinputName;             ///< name of the external input Mass histogram
  TString        fpTinputName;            ///< name of the external input pT histogram
  Bool_t         fFromTree;               ///< if true use TTree of jet 4-vectors as input single tracks
  TString        fPathTreeinputFile;      ///< path to the file where the external input Tree is (can be from alien)
  TString        fTreeinputName;          ///< name of the external input Tree
  TString        fBranchJDetName;         ///< name of the detector level jet branch in the TTree
  TTree*         fTreeJet4Vect;           //!<! tree containing the jet 4-vectors (input for embed.)
  Int_t          fCurrentEntry;           ///< Current TTree entry
  
  AliJetEmbeddingTask(const AliJetEmbeddingTask&);            // not implemented
  AliJetEmbeddingTask &operator=(const AliJetEmbeddingTask&); // not implemented
  
  /// \cond CLASSIMP
  ClassDef(AliJetEmbeddingTask, 6) /// Jet embedding task
  /// \endcond
};
#endif
