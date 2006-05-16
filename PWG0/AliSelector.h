#ifndef ALISELECTOR_H
#define ALISELECTOR_H

#include <TFile.h>
#include <TSelector.h>
#include <TChain.h>

#include <AliESD.h>
#include <AliHeader.h>
#include <AliRun.h>
#include <AliRunLoader.h>

class AliSelector : public TSelector {
  public:
    AliSelector(TTree *tree=0);
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
    AliRun* GetAliRun();

    TChain          *fChain;   //! pointer to the analyzed TTree or TChain

    AliESD*          fESD;     //! "ESD" branch in fChain
    AliHeader*       fHeader;  //! "TE" branch in fChain, contains event header

 private:
  void DeleteKinematicsFile();
  void DeleteRunLoader();

  TFile*        fKineFile;  //! pointer to Kinematics.root if the file was opened
  AliRunLoader* fRunLoader;    //! pointer to the RunLoader if galice.root was opened

  ClassDef(AliSelector,0);
};

#endif
