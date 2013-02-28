#ifndef ALIJETCONSTITUENTTAGCOPIER_H
#define ALIJETCONSTITUENTTAGCOPIER_H

// $Id: AliJetConstituentTagCopier.h  $

#include "AliAnalysisTaskEmcal.h"

class TClonesArray;
class TString;

class AliJetConstituentTagCopier : public AliAnalysisTaskEmcal {
 public:
  AliJetConstituentTagCopier();
  AliJetConstituentTagCopier(const char *name);
  virtual ~AliJetConstituentTagCopier();

  void                        SetMCParticlesName(const char *n)      { fMCParticlesName       = n         ; }

 protected:
  void                        ExecOnce();
  Bool_t                      Run();
  void                        DoClusterLoop(TClonesArray *array);
  void                        DoTrackLoop(TClonesArray *array);
  void                        DoEmcalParticleLoop(TClonesArray *array);

  TString                     fMCParticlesName;                       // name of MC particle collection
  TClonesArray               *fMCParticles;                           //!MC particle collection
  TH1I                       *fMCParticlesMap;                        //!MC particle map

 private:
  AliJetConstituentTagCopier(const AliJetConstituentTagCopier&);            // not implemented
  AliJetConstituentTagCopier &operator=(const AliJetConstituentTagCopier&); // not implemented

  ClassDef(AliJetConstituentTagCopier, 1) // Copy tags from particle level constituent to detector level
};

#endif
