#ifndef ALITRDGTUTMU_H
#define ALITRDGTUTMU_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDgtuTMU.h 27496 2008-07-22 08:35:45Z cblume $ */

//--------------------------------------------------------------------
//
// This class simulates the tracklet processing in a TMU
//
//--------------------------------------------------------------------

#include "TObject.h"
#include "TList.h"

#include "AliTRDtrackletGTU.h"
#include "AliTRDgtuParam.h"

class TTree;
class TBranch;
class AliTRDtrackGTU;
class AliESDEvent;

class AliTRDgtuTMU : public TObject {
 public:
  AliTRDgtuTMU(Int_t stack = -1, Int_t sector = -1);
  ~AliTRDgtuTMU();

  Bool_t SetSector(Int_t sector);
  Bool_t SetStack(Int_t stack);

  Bool_t AddTracklet(AliTRDtrackletGTU *tracklet, Int_t link);

  Bool_t RunTMU(TList *ListOfTracks = 0x0, AliESDEvent *esd = 0x0);
  Bool_t Reset();

  // ----- successive stages of the processing in the TMU -----
  Bool_t RunInputUnit(Int_t layer);
  Bool_t RunZChannelUnit(Int_t layer);
  Bool_t RunTrackFinder(Int_t zchannel, TList* ListOfTracks);
  Bool_t RunTrackMerging(TList* ListOfTracks);
  Bool_t RunTrackReconstruction(TList* ListOfTracks);

  Bool_t CalculateTrackParams(AliTRDtrackGTU *track);
  Bool_t Uniquifier(TList* inlist, TList *outlist);
  Bool_t CalculatePID(AliTRDtrackGTU *track);

protected:
  TObjArray **fTracklets; // holding all tracklets per link
  TObjArray **fTrackletsPostInput; // holding all tracklets of a layer
				   // after sorting/calculation in input units
  TList **fZChannelTracklets; // holding all tracklets for layer and z-channel
  TList **fTracks; // lists of tracks
  AliTRDgtuParam *fGtuParam; // pointer to the instance of the GtuParam class

  Int_t fStack;			// Stack of this TMU
  Int_t fSector;		// Sector of this TMU

 private:
  AliTRDgtuTMU(const AliTRDgtuTMU &rhs);              // not implemented
  AliTRDgtuTMU& operator=(const AliTRDgtuTMU &rhs);   // not implemented

  ClassDef(AliTRDgtuTMU, 1);
};

#endif
