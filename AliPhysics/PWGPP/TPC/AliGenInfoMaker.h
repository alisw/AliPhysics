#ifndef ALIGENINFOMAKER_H
#define ALIGENINFOMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



//////////////////////////////////////////////////////////////////////////////
//                          Class AliGenInfoMaker                           //
//   collect together MC info for comparison purposes - effieciency studies and so on//                                                                 //
//   marian.ivanov@cern.ch                                                  //
//////////////////////////////////////////////////////////////////////////////




#include <TParticle.h>
#include "AliAnalysisTask.h"
#include "AliTrackReference.h"

class TFile;
class AliRunLoader;
class AliStack;
class AliTPCParam;
class AliMCEventHandler;
class AliMCInfo;

////////////////////////////////////////////////////////////////////////
// 
// Start of implementation of the class AliGenInfoMaker
//
////////////////////////////////////////////////////////////////////////

class AliGenInfoMaker : public TObject {

public:
  AliGenInfoMaker();
  virtual ~AliGenInfoMaker();
  // 
  //
  AliGenInfoMaker(const char * fnGalice, const char* fnRes,
		   Int_t nEvents=1, Int_t firstEvent=0);
  //event by event function - used in the analysis task
  Int_t ProcessEvent(AliMCEventHandler* mcinfo);

  Int_t ProcessEvent();   // process event
  Int_t TreeKLoop();      // process kinamatics
  Int_t TreeTRLoop();     // process track refereces
  Int_t TreeDLoop();      // process digits tree
  Int_t BuildKinkInfo();  // build information about MC kinks
  Int_t BuildV0Info();    // build information about MC kinks
  //
  //
  Int_t Exec();  
  void CreateTreeGenTracks();
  void CloseOutputFile();
  void SetFirstEventNr(Int_t i) {fFirstEventNr = i;}
  void SetNEvents(Int_t i) {fNEvents = i;}
  void SetDebug(Int_t level) {fDebug = level;}
  Int_t SetIO(Int_t eventNr);
  Int_t CloseIOEvent();
  Int_t CloseIO();
  Int_t SetIO();

protected:
  AliGenInfoMaker(const AliGenInfoMaker& /*info*/);
  AliGenInfoMaker& operator=(const AliGenInfoMaker& /*info*/) { return *this;}

  AliMCInfo * MakeInfo(UInt_t i);
  AliMCInfo * GetInfo(UInt_t i) const {return (i<fNParticles)? fGenInfo[i]:0;}
  Float_t TR2LocalX(AliTrackReference *trackRef,
		    AliTPCParam *paramTPC) const;
  AliTPCParam * GetTPCParam();
  //
  TObjArray *fGenTracksArray;  //clones array with filtered particles
  TObjArray *fGenKinkArray;    //clones array with filtered Kinks
  TObjArray *fGenV0Array;      //clones array with filtered V0s
  //
  Int_t  fDebug;                   //! debug flag  
  Int_t  fEventNr;                 //! current event number
  Int_t  fLabel;                   //! track label
  Int_t  fNEvents;                 //! number of events to process
  Int_t  fFirstEventNr;            //! first event to process
  UInt_t fNParticles;              //! number of particles in TreeK
  TTree *fTreeGenTracks;          //! output tree with generated tracks
  TTree *fTreeKinks;             //!  output tree with Kinks
  TTree *fTreeV0;                //!  output tree with V0
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
  Float_t  fVPrim[3];             //! primary vertex position                                  // the fVDist[3] contains size of the 3-vector
  // cuts
  //
  Double_t fTPCPtCut; // do not store particles with generated pT less than this
  Double_t fITSPtCut; // do not store particles with generated pT less than this
  Double_t fTRDPtCut; // do not store particles with generated pT less than this
  Double_t fTOFPtCut; // do not store particles with generated pT less than this
 
  ClassDef(AliGenInfoMaker,0)    // class which creates and fills tree with TPCGenTrack objects
};






#endif
