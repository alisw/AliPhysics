/* $Id$ */

#ifndef ALISELECTORRL_H
#define ALISELECTORRL_H

// Use this selector, if you have AliROOT installed (at the moment only ESD, STEER + deps)

#include "AliSelector.h"

class AliRun;
class AliRunLoader;

class AliSelectorRL : public AliSelector {
  public:
    AliSelectorRL();
    virtual ~AliSelectorRL();

    virtual Bool_t  Notify();
    virtual void    SlaveTerminate();

 protected:
    AliRun* GetAliRun();

 private:
    void DeleteRunLoader();
    AliRunLoader* fRunLoader;    //! pointer to the RunLoader if galice.root was opened

    ClassDef(AliSelectorRL,0);
};

#endif
