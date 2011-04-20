#ifndef ALITRDGTUSIM_H
#define ALITRDGTUSIM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDgtuSim.h 27496 2008-07-22 08:35:45Z cblume $ */

// --------------------------------------------------------
// 
// GTU simulation
//
// --------------------------------------------------------

#include "TObject.h"

class AliRunLoader;
class AliLoader;
class AliESDEvent;

class AliTRDgtuTMU;
class AliTRDfeeParam;
class TTree;
class TList;

class AliTRDgtuSim : public TObject {
 public:
  AliTRDgtuSim(AliRunLoader *rl = 0x0);
  ~AliTRDgtuSim();

  Bool_t LoadTracklets(AliLoader * const loader);

  Bool_t RunGTU(AliLoader *loader, AliESDEvent *esd = 0x0);
  Bool_t RunGTUFromTrackletFile(TString filename, Int_t event, Int_t noev = 1);

  TTree* GetTreeOfTracks() { return fTrackTree; }
  Bool_t WriteTracksToTree(TList *listOfTracks, Int_t event = 0); 
  Bool_t WriteTracksToDataFile(TList *listOfTracks, Int_t event);
  Bool_t WriteTreesToFile() const;
  Bool_t WriteTracksToESD(const TList *const listOfTracks, AliESDEvent *esd);
  Bool_t WriteTracksToLoader(const TList *const listOfTracks);

 protected:
  AliRunLoader 	*fRunLoader;  	//!
  AliTRDfeeParam *fFeeParam;    //!
  AliTRDgtuTMU 	*fTMU; 		// pointer to TMU simulation class
  TClonesArray 	*fTrackletArray;	// array of tracklets
  TTree 	*fTrackTree; 	// tree to hold the tracks of one event, used for writing in WriteTracksToFile()
  TTree         *fTrackletTree; // tree to hold the gtu tracklets

 private:
  AliTRDgtuSim& operator=(const AliTRDgtuSim &rhs); // not implemented
  AliTRDgtuSim(const AliTRDgtuSim &rhs); // not implemented

  ClassDef(AliTRDgtuSim, 1);
};

#endif
