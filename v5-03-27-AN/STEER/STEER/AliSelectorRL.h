/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#ifndef ALISELECTORRL_H
#define ALISELECTORRL_H

// Use this selector, if you have AliROOT installed (at the moment only ESD, STEER + deps)

#include "AliSelector.h"

class AliRunLoader;
class AliHeader;
class AliStack;

class AliSelectorRL : public AliSelector {
  public:
    AliSelectorRL();
    virtual ~AliSelectorRL();

    virtual Bool_t  Notify();
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();

 protected:
    AliRunLoader* GetRunLoader();
    AliHeader* GetHeader();
    AliStack* GetStack();

 private:
    void DeleteRunLoader();

    AliRunLoader* fRunLoader;    //! pointer to the RunLoader if galice.root was opened
    Bool_t fKinematicsLoaded;    // determines if the stack is properly loaded (AliRunLoader::LoadKinematics() succeeded), this is needed because the GetStack returnes a invalid stack object when the function failed
    Bool_t fHeaderLoaded;        // determines if the header is properly loaded

    AliSelectorRL(const AliSelectorRL&);
    AliSelectorRL& operator=(const AliSelectorRL&);

    ClassDef(AliSelectorRL,0);
};

#endif
