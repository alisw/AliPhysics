/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliFlowEventSimple_H
#define AliFlowEventSimple_H


#include "AliFlowVector.h" //needed as include
class AliFlowTrackSimple;

// AliFlowEventSimple:
// A simple event for flow analysis
// authors: N. van der Kolk (kolk@nikhef.nl), A. Bilandzic (anteb@nikhef.nl)


class AliFlowEventSimple: public TObject {

 public:
  AliFlowEventSimple(Int_t aLenght);
  AliFlowEventSimple(const AliFlowEventSimple& anEvent);
  AliFlowEventSimple& operator=(const AliFlowEventSimple& anEvent);
  virtual  ~AliFlowEventSimple();
  
  Int_t NumberOfTracks() const              { return this->fNumberOfTracks; }
  Int_t GetEventNSelTracksIntFlow() const   { return this->fEventNSelTracksIntFlow; }
  void SetNumberOfTracks(Int_t nr)          { this->fNumberOfTracks = nr; }
  void SetEventNSelTracksIntFlow(Int_t nr)  { this->fEventNSelTracksIntFlow = nr; }
  AliFlowTrackSimple* GetTrack(Int_t i);
  TObjArray* TrackCollection() const        { return this->fTrackCollection; }
  AliFlowVector GetQ(Int_t n=2);//default harmonic n=2
  
 private:
  TObjArray*           fTrackCollection;         // collection of tracks
  Int_t                fNumberOfTracks;          // number of tracks
  Int_t                fEventNSelTracksIntFlow;  // number of tracks selected for integrated flow calculation
  
  ClassDef(AliFlowEventSimple,1)                 // macro for rootcint
};

#endif


