#ifndef ALIMUONTRACK_H
#define ALIMUONTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

#include <TROOT.h>
#include <TClonesArray.h>
#include "AliMUONTrackHit.h"
#include "AliMUONTrackParam.h"

class TClonesArray;
class AliMUONEventReconstructor;
class AliMUONHitForRec;
class AliMUONSegment;

class AliMUONTrack : public TObject {
 public:
  AliMUONTrack(){
    // Constructor
    ;} // Constructor
  virtual ~AliMUONTrack(){
    // Destructor
    ;} // Destructor
  AliMUONTrack (const AliMUONTrack& AliMUONTrack); // copy constructor
  AliMUONTrack& operator=(const AliMUONTrack& AliMUONTrack); // assignment operator

  AliMUONTrack(AliMUONSegment* BegSegment, AliMUONSegment* EndSegment, AliMUONEventReconstructor* EventReconstructor); // Constructor from two Segment's
  AliMUONTrack(AliMUONSegment* Segment, AliMUONHitForRec* HitForRec, AliMUONEventReconstructor* EventReconstructor); // Constructor from one Segment and one HitForRec

  AliMUONEventReconstructor* GetEventReconstructor(void) {return fEventReconstructor;};
  AliMUONTrackParam* GetTrackParamAtVertex(void);
  void SetTrackParamAtVertex(void); // Set track parameters at vertex from last stations 4 & 5
  void SetTrackParamAtVertex(AliMUONTrackParam* TrackParam) {fTrackParamAtVertex = *TrackParam;};

  TClonesArray* GetTrackHitsPtr(void);
  Int_t GetNTrackHits(void);
  Int_t GetFitMCS(void) {return fFitMCS;};
  void SetFitMCS(Int_t FitMCS);

  AliMUONTrackParam* GetTrackParamAtFirstHit(void);

  void RecursiveDump(void); // Recursive dump (with track hits)
  void Fit(AliMUONTrackParam *TrackParam, Int_t NParam); // Fit
  void AddSegment(AliMUONSegment* Segment); // Add Segment
  void AddHitForRec(AliMUONHitForRec* HitForRec); // Add HitForRec
  void SetTrackParamAtHit(Int_t indexHit, AliMUONTrackParam *TrackParam);

 protected:
 private:
  AliMUONEventReconstructor* fEventReconstructor; // Pointer to EventReconstructor
  AliMUONTrackParam fTrackParamAtVertex; // Track parameters at vertex
  TClonesArray *fTrackHitsPtr; // Pointer to array of TrackHit's
  Int_t fNTrackHits; // Number of TrackHit's
  Int_t fFitMCS; // 0(1) for fit without(with) multiple Coulomb scattering
  
  ClassDef(AliMUONTrack, 1) // Class definition in ROOT context
    };
	
#endif
