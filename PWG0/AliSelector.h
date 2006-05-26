/* $Id$ */

#ifndef ALISELECTOR_H
#define ALISELECTOR_H

// This selector is only dependent on the ESD library, if you need the whole of AliROOT use AliSelectorRL

#include <TFile.h>
#include <TSelector.h>
#include <TChain.h>

#include <AliESD.h>
#include <AliHeader.h>

class TParticle;

class AliSelector : public TSelector {
  public:
    AliSelector();
    virtual ~AliSelector();

    virtual Int_t   Version() const {return 1;}
    virtual void    Begin(TTree *tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Notify();
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    TTree*  GetKinematics();
    AliHeader* GetHeader();
    Bool_t IsPrimaryCharged(TParticle* aParticle, Int_t aTotalPrimaries) const;

    TChain          *fChain;   //! pointer to the analyzed TTree or TChain

    AliESD*          fESD;     //! "ESD" branch in fChain

    Int_t fCountFiles;         // Nr. of current file

 private:
    void DeleteKinematicsFile();
    void DeleteHeaderFile();

    TFile*        fKineFile;   //! pointer to Kinematics.root if the file was opened

    TFile*        fHeaderFile; //! pointer to galice.root, if the file was opened
    TTree*        fHeaderTree; //!
    AliHeader*    fHeader;     //!

  ClassDef(AliSelector,0);
};

#endif
