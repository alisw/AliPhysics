//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec  27 19:14:54 2005 by ROOT version 5.12/01
// from TTree esdTree/Tree with ESD objects
// found on file: AliESDs.root
//////////////////////////////////////////////////////////

#ifndef EsdQa_h
#define EsdQa_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <AliESD.h>           //dtor deletes fEsd

class TH2F;
class TH1F;

class EsdQa : public TSelector {

 public :
           EsdQa(TTree *tree=0) : TSelector(), fChain(0), fEsd(0), fCkovMomH(0), fdist(0), fX(0), fY(0), fsigma2(0)  
                                                                                                         {for(Int_t i=0;i<5;i++) fProbH[i]=0;}
  virtual ~EsdQa(             )                                                                          {delete fEsd;}


  virtual Int_t   Version        () const {return 2;}
  virtual void    Begin          (TTree *tree);
  virtual void    SlaveBegin     (TTree *tree);
  virtual void    Init           (TTree *tree);
  virtual Bool_t  Notify         ();
  virtual Bool_t  Process        (Long64_t entry);
  virtual void    SetOption      (const char *option) { fOption = option; }
  virtual void    SetObject      (TObject *obj) { fObject = obj; }
  virtual void    SetInputList   (TList *input) {fInput = input;}
  virtual TList  *GetOutputList  () const { return fOutput; }
  virtual void    SlaveTerminate ();
  virtual void    Terminate      ();

 private: 
  TTree          *fChain ;   //!pointer to the analyzed TTree or TChain
  AliESD         *fEsd ;     //!

  TH2F           *fCkovMomH,*fdist; //!
  TH1F           *fX, *fY, *fsigma2;
  TH1F           *fProbH[5]; //!

  ClassDef(EsdQa,0);
};

#endif
