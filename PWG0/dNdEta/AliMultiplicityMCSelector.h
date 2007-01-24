/* $Id$ */

#ifndef ALIMULTIPLICITYMCSELECTOR_H
#define ALIMULTIPLICITYMCSELECTOR_H

#include "AliSelectorRL.h"

class AliESDtrackCuts;
class TH2F;
class TH3F;
class AliMultiplicityCorrection;

class AliMultiplicityMCSelector : public AliSelectorRL {
  public:
    AliMultiplicityMCSelector();
    virtual ~AliMultiplicityMCSelector();

    virtual void    Begin(TTree* tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    void ReadUserObjects(TTree* tree);

    AliMultiplicityCorrection* fMultiplicity; // object containing the extracted data
    AliESDtrackCuts* fEsdTrackCuts;     // Object containing the parameters of the esd track cuts

 private:
    AliMultiplicityMCSelector(const AliMultiplicityMCSelector&);
    AliMultiplicityMCSelector& operator=(const AliMultiplicityMCSelector&);

  ClassDef(AliMultiplicityMCSelector, 0);
};

#endif
