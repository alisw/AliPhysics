/* $Id$ */

#ifndef ALIMULTIPLICITYESDSELECTOR_H
#define ALIMULTIPLICITYESDSELECTOR_H

#include "AliSelector.h"

class AliESDtrackCuts;
class AliMultiplicityCorrection;

class AliMultiplicityESDSelector : public AliSelector {
  public:
    AliMultiplicityESDSelector();
    virtual ~AliMultiplicityESDSelector();

    virtual void    Begin(TTree* tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    void ReadUserObjects(TTree* tree);

    AliMultiplicityCorrection* fMultiplicity; // object containing the extracted data
    AliESDtrackCuts* fEsdTrackCuts;           // Object containing the parameters of the esd track cuts

 private:
    AliMultiplicityESDSelector(const AliMultiplicityESDSelector&);
    AliMultiplicityESDSelector& operator=(const AliMultiplicityESDSelector&);

  ClassDef(AliMultiplicityESDSelector, 0);
};

#endif
