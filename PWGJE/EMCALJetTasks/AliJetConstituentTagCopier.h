#ifndef ALIJETCONSTITUENTTAGCOPIER_H
#define ALIJETCONSTITUENTTAGCOPIER_H

// $Id$

#include "AliAnalysisTaskEmcal.h"

class TString;
class AliNamedArrayI;
class AliParticleContainer;
class AliClusterContainer;

class AliJetConstituentTagCopier : public AliAnalysisTaskEmcal {
 public:
  AliJetConstituentTagCopier();
  AliJetConstituentTagCopier(const char *name);
  virtual ~AliJetConstituentTagCopier();

  void                        ConnectMCParticleContainerID(AliParticleContainer *cont)  { fMCParticleContainer   = cont      ; }
  void                        SetCleanBeforeCopy(Bool_t c)                              { fCleanBeforeCopy       = c         ; }
  void                        SetMCLabelShift(Int_t s)                                  { fMCLabelShift          = s         ; }

 protected:
  Bool_t                      Run();
  void                        DoClusterLoop(AliClusterContainer *cont);
  void                        DoParticleLoop(AliParticleContainer *cont);

  Bool_t                      fCleanBeforeCopy;                       // clean bit map before copying
  Int_t                       fMCLabelShift;                          // if MC label > fMCLabelShift, MC label -= fMCLabelShift
  AliParticleContainer       *fMCParticleContainer;                   // MC particle container

 private:
  AliJetConstituentTagCopier(const AliJetConstituentTagCopier&);            // not implemented
  AliJetConstituentTagCopier &operator=(const AliJetConstituentTagCopier&); // not implemented

  ClassDef(AliJetConstituentTagCopier, 4) // Copy tags from particle level constituent to detector level
};

#endif
