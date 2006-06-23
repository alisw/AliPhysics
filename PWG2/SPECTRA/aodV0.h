//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug 17 17:16:33 2005 by ROOT version 5.02/00
// from TTree esdTree/Tree with ESD objects
// found on file: AliESDs.root
// and modified by P.Hristov
//
// updated for using class ANALYSIS/AliAODv0 by B.Hippolyte
//////////////////////////////////////////////////////////

#ifndef aodV0_h
#define aodV0_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TColor.h>
#include <TPaveText.h>

class AliESD;


class aodV0 : public TSelector {
  public:

  aodV0(TTree *tree=0);
  virtual ~aodV0();

  virtual Int_t   Version() const {return 1;}
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) {fInput = input;}
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();

 private:
  TTree          *fChain;   //! pointer to the analyzed TTree or TChain

  // Declaration of leave types
  AliESD          *fESD;                     //! 

  // Histograms
  TH1F*           fHistV0PerEvent;           //!
  TH1F*           fHistMassK0;               //!
  TH1F*           fHistMassLambda;           //!
  TH1F*           fHistMassAntiLambda;       //!
  TH2F*           fHistMassLambdaVsProb;     //!
  TH1F*           fHistMassLambdaCut;        //!
  TH1F*           fHistMassAntiLambdaCut;    //!

  // AOD Histograms
  TH2F*           fHistPtVsRapK0Short;       //!
  TH2F*           fHistPtVsRapLambda;        //!
  TH2F*           fHistArmenterosPodolanski; //!

  ClassDef(aodV0,0);
};

#endif

