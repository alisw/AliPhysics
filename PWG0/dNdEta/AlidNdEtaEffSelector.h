#ifndef ALIDNDETAEFFSELECTOR_H
#define ALIDNDETAEFFSELECTOR_H

#include <TH2F.h>

#include "../AliSelector.h"
#include "../esdTrackCuts/AliESDtrackCuts.h"
#include "dNdEtaCorrection.h"

class TParticle;

class AlidNdEtaEffSelector : public AliSelector {
  public:
    AlidNdEtaEffSelector(TTree *tree=0);
    virtual ~AlidNdEtaEffSelector();

    virtual void    Begin(TTree *tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual Bool_t  Notify();
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
  Bool_t IsPrimary(const TParticle* aParticle, Int_t aTotalPrimaries);

  AliESDtrackCuts* fEsdTrackCuts;
  dNdEtaCorrection* fdNdEtaCorrection;

 private:

  ClassDef(AlidNdEtaEffSelector,0);
};

#endif
