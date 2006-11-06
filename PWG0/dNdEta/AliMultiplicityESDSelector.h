/* $Id$ */

#ifndef ALIMULTIPLICITYESDSELECTOR_H
#define ALIMULTIPLICITYESDSELECTOR_H

#include "AliSelector.h"

// uncomment this to enable mona lisa monitoring
//#define ALISELECTOR_USEMONALISA

class AliESDtrackCuts;
class TH1F;

#ifdef ALISELECTOR_USEMONALISA
  class TMonaLisaWriter;
#endif

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

#ifdef ALISELECTOR_USEMONALISA
    TMonaLisaWriter* fMonaLisaWriter; //! ML instance for monitoring
#endif

  ClassDef(AliMultiplicityESDSelector, 0);
};

#endif
