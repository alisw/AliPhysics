//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug 17 17:16:33 2005 by ROOT version 5.02/00
// from TTree esdTree/Tree with ESD objects
// found on file: AliESDs.root
// and modified by M. Ivanov
//////////////////////////////////////////////////////////

#ifndef PointMaker_h
#define PointMaker_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1F.h>
#include <TH2F.h>

#include  "AliESDEvent.h"
#include  "AliESDfriend.h"
class AliTrackPointArray;

class PointMaker : public TSelector {
  public:

  PointMaker(char *outfil="AliTrackPoints.root");
  virtual ~PointMaker();

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
  Bool_t          IsIdenticalWithOneOf(AliTrackPoint *p, AliTrackPointArray *parray, int nmax=kMaxInt); 
  TTree          *fChain;           //! pointer to the analyzed TTree or TChain
  // Declaration of leave types
  AliESDEvent    *fESD;             //! 
  AliESDfriend   *fESDfriend;       //! 
  TFile          *fFile;            //! output file
  TTree          *fTree;            //! pointer to the output TTree 
  AliTrackPointArray *fArray;       // pointer to the track points
  Int_t          fNevents;          // number of events
  Int_t          fNtracks;          // number of tracks
  Int_t          fNAcceptedTracks;  // number of accepted tracks
  TString        fOutfil;           // output filename
 public:
  TH1D           *fCuttra;          // histogram for monitoring track cuts
  TH1D           *fCutpoi;          // histogram for monitoring point cuts
  TH2D           *fModpop;          // histogram for monitoring module population

  ClassDef(PointMaker,0);
};

#endif

