#ifndef ALIMUONTRACK_H
#define ALIMUONTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

#include "AliMUONTrackParam.h" // object belongs to the class

class TObjArray;
class TClonesArray;
class TVirtualFitter;
class AliMUONEventReconstructor;
class AliMUONHitForRec;
class AliMUONSegment;

class AliMUONTrack : public TObject {
 public:
  AliMUONTrack(){
    // Constructor
    ;} // Constructor
  virtual ~AliMUONTrack(); // Destructor
  AliMUONTrack (const AliMUONTrack& AliMUONTrack); // copy constructor
  AliMUONTrack& operator=(const AliMUONTrack& AliMUONTrack); // assignment operator

  AliMUONTrack(AliMUONSegment* BegSegment, AliMUONSegment* EndSegment, AliMUONEventReconstructor* EventReconstructor); // Constructor from two Segment's
  AliMUONTrack(AliMUONSegment* Segment, AliMUONHitForRec* HitForRec, AliMUONEventReconstructor* EventReconstructor); // Constructor from one Segment and one HitForRec
  void Remove(void);

  AliMUONEventReconstructor* GetEventReconstructor(void) {return fEventReconstructor;}
  AliMUONTrackParam* GetTrackParamAtVertex(void) { return &fTrackParamAtVertex;}
  void SetTrackParamAtVertex(void); // Set track parameters at vertex from last stations 4 & 5
  void SetTrackParamAtVertex(AliMUONTrackParam* TrackParam) {fTrackParamAtVertex = *TrackParam;}

  TObjArray* GetTrackHitsPtr(void) { return fTrackHitsPtr;}
  Int_t GetNTrackHits(void) { return fNTrackHits;}
  Int_t GetFitMCS(void) {return fFitMCS;}
  Int_t GetFitNParam(void) {return fFitNParam;}
  Int_t GetFitStart(void) {return fFitStart;}
  Double_t GetFitFMin(void) {return fFitFMin;}
  void SetFitMCS(Int_t FitMCS);
  void SetFitNParam(Int_t FitNParam);
  void SetFitStart(Int_t FitStart);

  AliMUONTrackParam* GetTrackParamAtFirstHit(void);

  void RecursiveDump(void); // Recursive dump (with track hits)
  void Fit(); // Fit
  void AddSegment(AliMUONSegment* Segment); // Add Segment
  void AddHitForRec(AliMUONHitForRec* HitForRec); // Add HitForRec
  void SetTrackParamAtHit(Int_t indexHit, AliMUONTrackParam *TrackParam);
  Int_t HitsInCommon(AliMUONTrack* Track);

  static TVirtualFitter* Fitter(void) {return fgFitter;}

 protected:
 private:
  static TVirtualFitter* fgFitter; // Pointer to track fitter
  AliMUONEventReconstructor* fEventReconstructor; // Pointer to EventReconstructor
  AliMUONTrackParam fTrackParamAtVertex; // Track parameters at vertex
  TObjArray *fTrackHitsPtr; // Pointer to array of pointers to TrackHit's
  Int_t fNTrackHits; // Number of TrackHit's
  Int_t fFitMCS; // 0(1) for fit without(with) multiple Coulomb scattering
  Int_t fFitNParam; // 3(5) for fit with 3(5) parameters
  Int_t fFitStart; // 0 or 1 for fit starting from parameters at vertex (0) or at first TrackHit(1)
  Double_t fFitFMin; // minimum value of the function minimized by the fit
  
  ClassDef(AliMUONTrack, 1) // Reconstructed track in ALICE dimuon spectrometer
    };
	
#endif
