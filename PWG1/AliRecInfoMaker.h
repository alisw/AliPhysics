#ifndef ALIRECINFOMAKER_H
#define ALIRECINFOMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



//////////////////////////////////////////////////////////////////////////////
//                          Class AliRecInfoMaker                           //
//   collect together MC info and Rec info for comparison purposes 
//                                           - effieciency studies and so on//                                                                 //
//   marian.ivanov@cern.ch                                                  //
//////////////////////////////////////////////////////////////////////////////


#include "TObject.h"
#include "AliESD.h"
#include "AliESDEvent.h"

#include "AliESDtrack.h"
#include "AliV0.h"
#include "AliESDfriendTrack.h"
#include "AliITStrackMI.h"
#include "AliTRDtrack.h"
class AliTPCseed;




////////////////////////////////////////////////////////////////////////
// 
// Start of implementation of the class AliRecInfoMaker
//
////////////////////////////////////////////////////////////////////////

class AliRecInfoMaker {

public:
  AliRecInfoMaker(const char* fnGenTracks,
	   const char* fnCmpRes      ="cmpTracks.root", 
	   const char* fnGalice      ="galice.root",
	   Int_t nEvents=1, Int_t firstEvent=0);
  static void MakeAliases(TTree *tree); 
  virtual ~AliRecInfoMaker();
  void Reset();
  Int_t Exec();
  Int_t Exec(Int_t nEvents, Int_t firstEventNr);
  Int_t SetIO();
  Int_t SetIO(Int_t eventNr );
  void CreateTreeCmp();
  void CloseOutputFile();
  Bool_t ConnectGenTree();
  Int_t TreeGenLoop(Int_t eventNr);
  Int_t TreeTLoop();
  Int_t BuildKinkInfo0(Int_t eventNr); // build kink info 0
  Int_t BuildV0Info(Int_t eventNr); // build kink info 0
  void SetFirstEventNr(Int_t i) {fFirstEventNr = i;}
  void SetNEvents(Int_t i) {fNEvents = i;}
  void SetDebug(Int_t level) {fDebug = level;}

// tmp method, should go to TrackReferenceESD
  static TVector3 TR2Local(AliTrackReference *trackRef,
		    AliTPCParam *paramTPC);
  static AliTPCParam * GetTPCParam();
private:

  Int_t fEventNr;                 //! current event number
  Int_t fNEvents;                 //! number of events to process
  Int_t fFirstEventNr;            //! first event to process
  //
  char  fFnCmp[1000];                   //! output file name with cmp tracks
  TFile *fFileCmp;                //! output file with cmp tracks
  TTree *fTreeCmp;                //! output tree with cmp tracks
  TTree *fTreeCmpKinks;                //! output tree with cmp Kinks
  TTree *fTreeCmpV0;                //! output tree with cmp V0
  //
  char  fFnGenTracks[1000];             //! input file name with gen tracks
  TFile *fFileGenTracks;                //! input files with generated tracks   
  TTree *fTreeGenTracks;           //! tree with generated tracks
  TTree *fTreeGenKinks;            // tree with gen kinks
  TTree *fTreeGenV0;            // tree with gen V0
  //
  //
  //
  AliRunLoader * fLoader;         //! pointer to the run loader
  //TTree *fTreeRecTracks;          //! tree with reconstructed tracks
  //
  Short_t *fIndexRecTracks;         //! index of particle label in the TreeT_ESD
  Short_t *fFakeRecTracks;          //! number of fake tracks
  Short_t *fMultiRecTracks;         //! number of multiple reconstructions
  //
  Short_t *fIndexRecKinks;         //! index of particle label in treeesd
  Short_t *fMultiRecKinks;         //! number of multiple reconstructions
  Short_t *fSignedKinks;           //! indicator that kink was not fake
  //
  Short_t *fIndexRecV0;         //! index of particle label in treeesd
  Short_t *fMultiRecV0;         //! number of multiple reconstructions
  Short_t *fSignedV0;                //! indicator that kink was not fake
  //
  TObjArray *fRecArray;           // container with rec infos
  AliESDEvent *fEvent;             //!event
  AliESDfriend *fESDfriend;              //!event friend
  //
  AliTPCParam* fParamTPC;         //! AliTPCParam
  Int_t fNParticles;              //! number of particles in the input tree genTracks
  Int_t fDebug;                   //! debug flag  
  Int_t fNextTreeGenEntryToRead;    //! last entry already read from genTracks tree
  Int_t fNextKinkToRead;            //! last entry already read from genKinks tree
  Int_t fNextV0ToRead;            //! last entry already read from genV0 tree
  //
  AliMCInfo*  fMCInfo;           //! MC information writen per particle
  AliGenKinkInfo* fGenKinkInfo;      //! MC information writen per Kink
  AliGenV0Info* fGenV0Info;      //! MC information writen per Kink
  AliESDRecInfo*  fRecInfo;          //! Rec. information writen per particle
  AliESDfriendTrack*  fFriend;          //! friend track
  AliESDRecKinkInfo* fRecKinkInfo;    //! reconstructed kink info
  AliESDRecV0Info* fRecV0Info;    //! reconstructed kink info
  //

  ClassDef(AliRecInfoMaker,1)    // class which creates and fills tree with ESDGenTrack objects
};


#endif
