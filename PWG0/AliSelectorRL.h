/* $Id$ */

#ifndef ALISELECTORRL_H
#define ALISELECTORRL_H

// Use this selector, if you have AliROOT installed (at the moment only ESD, STEER + deps)

#include "AliSelector.h"

class AliRunLoader;
class AliHeader;

class AliSelectorRL : public AliSelector {
  public:
    AliSelectorRL();
    virtual ~AliSelectorRL();

    virtual Bool_t  Notify();
    virtual void    SlaveTerminate();

 protected:
    AliRunLoader* GetAliRunLoader();
    AliHeader* GetHeader();

 private:
    void DeleteRunLoader();
    void DeleteHeaderFile();

    AliRunLoader* fRunLoader;    //! pointer to the RunLoader if galice.root was opened

    TFile*        fHeaderFile; //! pointer to galice.root, if the file was opened
    TTree*        fHeaderTree; //! holds TE tree of current galice.root
    AliHeader*    fHeader;     //! holds pointer to current header

    ClassDef(AliSelectorRL,0);
};

#endif
