/* $Id$ */

#ifndef ALIMULTIPLICITYMCSELECTOR_H
#define ALIMULTIPLICITYMCSELECTOR_H

#include "AliSelectorRL.h"

class AliESDtrackCuts;
class TH1F;
class TH2F;

class AliMultiplicityMCSelector : public AliSelectorRL {
  public:
    AliMultiplicityMCSelector();
    virtual ~AliMultiplicityMCSelector();

    virtual void    Begin(TTree* tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    void ReadUserObjects(TTree* tree);

    TH1F* fMultiplicityESD; // multiplicity histogram
    TH1F* fMultiplicityMC; // multiplicity histogram

    TH2F* fCorrelation; // (gene multiplicity) vs (meas multiplicity)

    AliESDtrackCuts*  fEsdTrackCuts;     // Object containing the parameters of the esd track cuts

 private:
    AliMultiplicityMCSelector(const AliMultiplicityMCSelector&);
    AliMultiplicityMCSelector& operator=(const AliMultiplicityMCSelector&);

  ClassDef(AliMultiplicityMCSelector, 0);
};

#endif
