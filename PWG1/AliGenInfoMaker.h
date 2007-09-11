#ifndef ALIGENINFOMAKER_H
#define ALIGENINFOMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



//////////////////////////////////////////////////////////////////////////////
//                          Class AliGenInfo                               //
//   collect together MC info for comparison purposes - effieciency studies and so on//                                                                 //
//   marian.ivanov@cern.ch                                                  //
//////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////
//
// Start of implementation of the class AliTPCdigitRow
//
////////////////////////////////////////////////////////////////////////

#include <TParticle.h>
#include "AliTrackReference.h"

class TFile;
class AliRunLoader;
class AliStack;
class AliTPCParam;


////////////////////////////////////////////////////////////////////////
// 
// Start of implementation of the class AliGenInfoMaker
//
////////////////////////////////////////////////////////////////////////

class AliGenInfoMaker {

public:
  AliGenInfoMaker();
  AliGenInfoMaker(const char * fnGalice, const char* fnRes    ="genTracks.root",
		   Int_t nEvents=1, Int_t firstEvent=0);
  virtual ~AliGenInfoMaker();
  Int_t Exec();
  Int_t Exec(Int_t nEvents, Int_t firstEventNr);
  void CreateTreeGenTracks();
  void CloseOutputFile();
  Int_t TreeKLoop();
  Int_t TreeTRLoop();
  Int_t TreeTRLoopNew(); 
  Int_t TreeDLoop();
  Int_t BuildKinkInfo();  // build information about MC kinks
  Int_t BuildV0Info();  // build information about MC kinks
  Int_t BuildHitLines();  // build information about MC kinks
  void SetFirstEventNr(Int_t i) {fFirstEventNr = i;}
  void SetNEvents(Int_t i) {fNEvents = i;}
  void SetDebug(Int_t level) {fDebug = level;}
  Int_t SetIO(Int_t eventNr);
  Int_t CloseIOEvent();
  Int_t CloseIO();
  Int_t SetIO();
  Float_t TR2LocalX(AliTrackReference *trackRef,
		    AliTPCParam *paramTPC) const;
  AliMCInfo * GetInfo(UInt_t i) const {return (i<fNParticles)? fGenInfo[i]:0;}
  AliMCInfo * MakeInfo(UInt_t i);

private:
  AliTPCParam * GetTPCParam();
  Float_t TPCBetheBloch(Float_t bg);
  Int_t  fDebug;                   //! debug flag  
  Int_t  fEventNr;                 //! current event number
  Int_t  fLabel;                   //! track label
  Int_t  fNEvents;                 //! number of events to process
  Int_t  fFirstEventNr;            //! first event to process
  UInt_t fNParticles;              //! number of particles in TreeK
  TTree *fTreeGenTracks;          //! output tree with generated tracks
  TTree *fTreeKinks;             //!  output tree with Kinks
  TTree *fTreeV0;                //!  output tree with V0
  TTree *fTreeHitLines;          //! tree with hit lines
  char   fFnRes[1000];             //! output file name with stored tracks
  TFile *fFileGenTracks;          //! output file with stored fTreeGenTracks
  //
  AliRunLoader * fLoader;         //! pointer to the run loader
  TTree * fTreeD;                 //! current tree with digits
  TTree * fTreeTR;                //! current tree with TR
  AliStack *fStack;               //! current stack
  // 
  AliMCInfo **   fGenInfo;    //! array with pointers to gen info
  Int_t   fNInfos;                  //! number of tracks with infos
  //
  AliTPCParam* fParamTPC;         //! AliTPCParam
  Float_t  fVPrim[3];             //! primary vertex position
                                  // the fVDist[3] contains size of the 3-vector
  // cuts
  //
  Double_t fTPCPtCut; // do not store particles with generated pT less than this
  Double_t fITSPtCut; // do not store particles with generated pT less than this
  Double_t fTRDPtCut; // do not store particles with generated pT less than this
  Double_t fTOFPtCut; // do not store particles with generated pT less than this
 
  ClassDef(AliGenInfoMaker,0)    // class which creates and fills tree with TPCGenTrack objects
};






#endif
