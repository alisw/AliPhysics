#ifndef ALIJETMODELCOPYTRACKS_H
#define ALIJETMODELCOPYTRACKS_H

// $Id$

class TClonesArray;
class TRandom3;
class AliVParticle;
class AliPicoTrack;

#include "AliAnalysisTaskEmcal.h"

class AliJetModelCopyTracks : public AliAnalysisTaskEmcal {
 public:
  enum ParticleMass {
    kMassive  = 0,
    kMassless = 1,
    kPionMass = 2
  };


  AliJetModelCopyTracks();
  AliJetModelCopyTracks(const char *name); 
  virtual ~AliJetModelCopyTracks();

  virtual void           UserCreateOutputObjects();

  void                   SetTracksOutName(const char *n)          { fTracksOutName   = n;    }
  void                   SetParticleMassType(ParticleMass m)      { fParticleMass    = m;    }

  void                   ExecOnce();
  Bool_t                 Run();

  void                   CopyTracks();

  TString                fTracksOutName;       // name of output track collection
  TClonesArray          *fTracksOut;           //!output track collection
  ParticleMass           fParticleMass;        // particle mass to use

  //Output objects
  TH1F     *fHistPtOut;                        //!pT spectrum of output particles
  
 private:
  AliJetModelCopyTracks(const AliJetModelCopyTracks&);            // not implemented
  AliJetModelCopyTracks &operator=(const AliJetModelCopyTracks&); // not implemented

  ClassDef(AliJetModelCopyTracks, 1) // copy tracks class
};
#endif
