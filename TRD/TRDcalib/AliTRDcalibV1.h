//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 27 15:22:59 2007 by ROOT version 5.16/00
// from TTree esdTree/Tree with ESD objects
// found on file: AliESDs.root
//////////////////////////////////////////////////////////

#ifndef AliTRDcalibV1_h
#define AliTRDcalibV1_h

#include <TROOT.h>
#include <TH2I.h>
#include <TH2F.h>
#include <TProfile2D.h>
#include <TProfile.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TSelector.h>
#include <TObject.h>

class AliESD;
class AliESDEvent;
class AliESDfriend;
class AliESDtrack;
class AliESDfriendTrack;

class AliTRDCalibraFillHisto;
class AliTRDtrackV1;


class AliTRDcalibV1 : public TSelector {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain

  //variables
  AliESD      *fESD;
  AliESDEvent *fev;
  AliESDfriend *fevf;
  TObject *fo;
  AliTRDtrackV1 *ft;
  const AliESDtrack *fesdTrack;
  AliESDfriendTrack *ffriendTrack;

  //calibration class
  AliTRDCalibraFillHisto *fcalib;

  //Store infos
  TH2I *fCH2d;
  TProfile2D *fPH2d;
  TProfile2D *fPRF2d;
  TH2F     *fVdriftLinear[540];



  // For the case no proof
  Int_t    fFileNo;

  AliTRDcalibV1(TTree *tree=0);
  virtual ~AliTRDcalibV1() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();
  void            CleanESD();
  Int_t           ReadEvent(Long64_t entry);


   ClassDef(AliTRDcalibV1,0);
};

#endif
