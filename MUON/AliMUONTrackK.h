#ifndef ALIMUONTRACKK_H
#define ALIMUONTRACKK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TROOT.h>
class TArrayD;
class TMatrixD;
class AliMUONEventReconstructor;
class TClonesArray;
class AliMUONSegment;
class AliMUON;
class AliMUONHitForRec;

class AliMUONTrackK : public TObject {

 public:

  AliMUONTrackK(); // Default constructor
  virtual ~AliMUONTrackK(); // Destructor
  AliMUONTrackK (const AliMUONTrackK& source); // copy constructor
  AliMUONTrackK& operator=(const AliMUONTrackK& source); // assignment operator

  //AliMUONTrackK(const AliMUONEventReconstructor *EventReconstructor, const AliMUONHitForRec *hitForRec); // Constructor
  AliMUONTrackK(AliMUONEventReconstructor *EventReconstructor, TClonesArray *hitForRec); // Constructor
  AliMUONTrackK(AliMUONSegment *segment); // Constructor from a segment

  // Pointer to hits on track
  TObjArray* GetHitOnTrack(void) const {return fTrackHitsPtr;} // ptr. to hits on track
  Int_t GetNTrackHits(void) const {return fNTrackHits;} // hits on track
  Double_t GetTrackQuality(void) const {return fChi2;} // track quality
  TMatrixD* GetTrackParameters(void) const {return fTrackPar;} // track parameters
  Double_t GetZ(void) const {return fPosition;} // Z-coordinate of track 
  TMatrixD* GetCovariance(void) const {return fCovariance;} // covariance matrix
  Int_t GetTrackDir(void) const {return fTrackDir;} // get track propagation direction
  void SetTrackDir(Int_t iDir) {fTrackDir = iDir;} // set track propagation direction  
  Bool_t GetBPFlag(void) const {return fBPFlag;} // get backpropagation flag
  void SetBPFlag(Bool_t BPFlag) {fBPFlag = BPFlag;} // set backpropagation flag
  Int_t GetRecover(void) const {return fRecover;} // return recover flag 
  void SetRecover(Int_t iRecover) {fRecover = iRecover;} // set recover flag
  AliMUONSegment* GetStartSegment(void) const {return fStartSegment;} // return seed segment
  Bool_t KalmanFilter(Int_t ichamBeg, Int_t ichamEnd, Bool_t Back, Double_t zDipole1, Double_t zDipole2); // Kalman filter
  void StartBack(void); // start backpropagator
  void SetTrackQuality(Int_t iChi2); // compute track quality or Chi2
  Bool_t KeepTrack(AliMUONTrackK* track0); // keep or discard track 
  void Kill(void); // kill track candidate
  void Branson(void); // Branson correction
  void GoToZ(Double_t zEnd); // propagate track to given Z
  void GoToVertex(void); // propagate track to the vertex

  // What is necessary for sorting TClonesArray's
  Bool_t IsSortable() const { return kTRUE; }
  Int_t Compare(const TObject* trackK) const; // "Compare" function for sorting



 protected:

 private:
 
  static Int_t fgNOfPoints; // number of points in event
  static AliMUON *fgMUON; // pointer to MUON module  
  static AliMUONEventReconstructor *fgEventReconstructor; // pointer to event reconstructor
  static TClonesArray *fgHitForRec; // pointer to hits

  AliMUONSegment *fStartSegment; // seed segment  
  Double_t fPosition; // Z-coordinate of track
  Double_t fPositionNew; //! Z-coordinate of track
  Double_t fChi2; // Chi2 of track
  TObjArray *fTrackHitsPtr; // pointer to hits on track
  Int_t fNTrackHits; // number of points on track
  Int_t fTrackDir; // track direction (+(-) towards high (low) z)
  Bool_t fBPFlag; // backpropagation flag (TRUE if backpropagation)
  Int_t fRecover; // recover flag (!=0 if recovery procedure was applied)
  AliMUONHitForRec *fSkipHit; // hit to skip during recovery procedure

  TMatrixD *fTrackPar; // track parameters
  TMatrixD *fTrackParNew; //! track parameters
  TMatrixD *fCovariance; // covariance matrix
  TMatrixD *fWeight; //! weight matrix (inverse of covariance)

  // Functions

  void EvalCovariance(Double_t dZ);
  void ParPropagation(Double_t zEnd);
  void WeightPropagation(Double_t zEnd);
  void MSThin(Int_t sign);
  void MSLine(Double_t dZ, Double_t X0);
  Bool_t FindPoint(Int_t ichamb, Double_t zEnd, Int_t currIndx, Int_t iFB, AliMUONHitForRec *&hitAdd);
  void TryPoint(TMatrixD &point, const TMatrixD &pointWeight, TMatrixD &trackParTmp, Double_t &dChi2);
  void SetGeantParam(Double_t *VGeant3, Int_t iFB);
  void GetFromGeantParam(Double_t *VGeant3, Int_t iFB);
  void Recover(void);

  private:
   // Some constants
   static const Int_t fgkSize; // number of track parameters
   static const Int_t fgkNSigma; //4; // acceptance window width in sigmas
   static const Int_t fgkTriesMax; // max number of attempts to find exact position during tracking
   static const Double_t fgkEpsilon; // tracking precision (cm)

  ClassDef(AliMUONTrackK,0) // Kalman track in MUON arm of ALICE
    };
#endif
