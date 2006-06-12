/* $Id$ */

#ifndef ALIDNDETACORRECTIONSELECTOR_H
#define ALIDNDETACORRECTIONSELECTOR_H

#include "AliSelectorRL.h"

class AliESDtrackCuts;
class dNdEtaCorrection;

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
    Bool_t  CheckVertex();
  
    AliESDtrackCuts*  fEsdTrackCuts;          // Object containing the parameters of the esd track cuts

    dNdEtaCorrection* fdNdEtaCorrection;      // contains the intermediate histograms (on each slave)
    dNdEtaCorrection* fdNdEtaCorrectionFinal; // contains the final histograms

 private:

  ClassDef(AlidNdEtaCorrectionSelector, 0);
};

#endif
