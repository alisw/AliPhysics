/* $Id$ */

#ifndef ALIMULTIPLICITYCORRECTIONSELECTOR_H
#define ALIMULTIPLICITYCORRECTIONSELECTOR_H

#include "AliSelectorRL.h"
#include "AliMultiplicityCorrection.h"

class AliESDtrackCuts;
class TH1F;
class TH2F;

class AliMultiplicityCorrectionSelector : public AliSelectorRL {
  public:
    AliMultiplicityCorrectionSelector();
    virtual ~AliMultiplicityCorrectionSelector();

    virtual void    Begin(TTree* tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    void ReadUserObjects(TTree* tree);

    AliMultiplicityCorrection* fMultiplicityCorrection;  // multiplicity histogram

    AliESDtrackCuts*           fEsdTrackCuts;  // Object containing the parameters of the esd track cuts

 private:
    AliMultiplicityCorrectionSelector(const AliMultiplicityCorrectionSelector&);
    AliMultiplicityCorrectionSelector& operator=(const AliMultiplicityCorrectionSelector&);

  ClassDef(AliMultiplicityCorrectionSelector, 0);
};

#endif
