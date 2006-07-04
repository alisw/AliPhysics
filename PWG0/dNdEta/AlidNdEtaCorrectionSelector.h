/* $Id$ */

#ifndef ALIDNDETACORRECTIONSELECTOR_H
#define ALIDNDETACORRECTIONSELECTOR_H

#include "AliSelectorRL.h"

class AliESDtrackCuts;
class AlidNdEtaCorrection;

class AlidNdEtaCorrectionSelector : public AliSelectorRL {
  public:
    AlidNdEtaCorrectionSelector();
    virtual ~AlidNdEtaCorrectionSelector();

    virtual void    Begin(TTree *tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    Bool_t SignOK(Double_t charge);

    AliESDtrackCuts*  fEsdTrackCuts;          // Object containing the parameters of the esd track cuts

    AlidNdEtaCorrection* fdNdEtaCorrection;      // contains the intermediate histograms (on each slave)

    Int_t fSignMode;  // if 0 process all particles, if +-1 process only particles with that sign

 private:

  ClassDef(AlidNdEtaCorrectionSelector, 0);
};

#endif
