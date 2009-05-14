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
class TTree;
class TList;

class AliTRDgtuSim : public TObject {
 public:
  AliTRDgtuSim(AliRunLoader *rl = 0x0);
  ~AliTRDgtuSim();

  Bool_t LoadTracklets(AliLoader *loader);

  Bool_t RunGTU(AliLoader *loader, AliESDEvent *esd = 0x0);
  Bool_t RunGTUFromTrackletFile(TString filename, Int_t event, Int_t noev = 1);

  TTree* GetTreeOfTracks() { return fTrackTree; }
  Bool_t WriteTracksToTree(TList *ListOfTracks, Int_t event = 0); 
  Bool_t WriteTracksToDataFile(TList *ListOfTracks, Int_t event);
  Bool_t WriteTreesToFile();
  Bool_t WriteTracksToESD(TList *ListOfTracks, AliESDEvent *esd);
  Bool_t WriteTracksToLoader();

 protected:
  AliRunLoader 	*fRunLoader;  	//!
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
