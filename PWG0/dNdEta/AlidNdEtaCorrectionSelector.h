#ifndef ALIDNDETACORRECTIONSELECTOR_H
#define ALIDNDETACORRECTIONSELECTOR_H

#include "AliSelector.h"

class AliESDtrackCuts;
class dNdEtaCorrection;

class TParticle;

class AlidNdEtaCorrectionSelector : public AliSelector {
  public:
    AlidNdEtaCorrectionSelector(TTree *tree=0);
    virtual ~AlidNdEtaCorrectionSelector();

    virtual void    Begin(TTree *tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual Bool_t  Notify();
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
  AliESDtrackCuts*  fEsdTrackCuts;          // Object containing the parameters of the esd track cuts

  dNdEtaCorrection* fdNdEtaCorrection;      // contains the intermediate histograms (on each slave)
  dNdEtaCorrection* fdNdEtaCorrectionFinal; // contains the final histograms

 private:

  ClassDef(AlidNdEtaCorrectionSelector, 0);
};

#endif
