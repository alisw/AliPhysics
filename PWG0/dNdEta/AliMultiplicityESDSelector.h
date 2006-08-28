/* $Id$ */

#ifndef ALIMULTIPLICITYESDSELECTOR_H
#define ALIMULTIPLICITYESDSELECTOR_H

#include "AliSelector.h"

class AliESDtrackCuts;
class TH1F;

class TMonaLisaWriter;

class AliMultiplicityESDSelector : public AliSelector {
  public:
    AliMultiplicityESDSelector();
    virtual ~AliMultiplicityESDSelector();

    virtual void    Begin(TTree* tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    void ReadUserObjects(TTree* tree);

    TH1F* fMultiplicity; // multiplicity histogram

    AliESDtrackCuts*  fEsdTrackCuts;     // Object containing the parameters of the esd track cuts

 private:
    AliMultiplicityESDSelector(const AliMultiplicityESDSelector&);
    AliMultiplicityESDSelector& operator=(const AliMultiplicityESDSelector&);

    TMonaLisaWriter* fMonaLisaWriter; //! ML instance for monitoring

  ClassDef(AliMultiplicityESDSelector, 0);
};

#endif
