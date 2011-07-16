/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#ifndef ALISELECTOR_H
#define ALISELECTOR_H

// This selector is only dependent on the ESD library, if you need the whole of AliROOT use AliSelectorRL

#include <TSelector.h>

class TTree;
class AliESD;

class AliSelector : public TSelector {
  public:
    AliSelector();
    virtual ~AliSelector();

    virtual Int_t   Version() const {return 1;}
    virtual void    Begin(TTree*);
    virtual void    SlaveBegin(TTree* tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Notify();
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    void CheckOptions();

    TTree          *fTree;     //! pointer to the TTree containing the events
    AliESD*          fESD;     //! "ESD" branch in fChain
    Int_t fCountFiles;         // number of processed file

 private:
    AliSelector(const AliSelector&);
    AliSelector& operator=(const AliSelector&);

  ClassDef(AliSelector,0);
};

#endif
