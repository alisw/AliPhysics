#ifndef ALIMUONTRACKK_H
#define ALIMUONTRACKK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMUONTrackK
/// \brief Kalman track in MUON arm of ALICE
///
/// \author Alexander Zinchenko, JINR Dubna

class TArrayD;
class TClonesArray;
class TObjArray;
#include <TMatrixDfwd.h>
#include <TObject.h>

class AliMUONEventRecoCombi;
class AliMUONHitForRec;
class AliMUONSegment;
class AliMUONTrackReconstructorK;
#include "AliMUONTrack.h" 

class AliMUONTrackK : public AliMUONTrack {

 public:

  AliMUONTrackK(); // Default constructor
  virtual ~AliMUONTrackK(); // Destructor

  AliMUONTrackK(AliMUONTrackReconstructorK *TrackReconstructor, TClonesArray *hitForRec); // Constructor
  AliMUONTrackK(AliMUONSegment *segment); // Constructor from a segment

  // Pointer to hits on track
  TObjArray* GetTrackHits(void) const {return fTrackHits;} // ptr. to hits on track
  Int_t GetNTrackHits(void) const {return fNmbTrackHits;} // hits on track
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
  Bool_t KeepTrack(AliMUONTrackK* track0) const; // keep or discard track 
  void Kill(void); // kill track candidate
  void Branson(void); // Branson correction
  void GoToZ(Double_t zEnd); // propagate track to given Z
  void GoToVertex(Int_t iflag); // propagate track to the vertex
  Bool_t Smooth(void); // apply smoother
  Double_t GetChi2PerPoint(Int_t iPoint) const; // return Chi2 at point
  void Print(FILE *lun) const; // print track information
  void Print(const char* /*opt*/) const {return;} // print track information
  AliMUONHitForRec* GetHitLastOk(void); // get hit before the skipped one
  Int_t GetStation0(void); // return seed station number
  Int_t DebugLevel(void) const {return fgDebug;} // return debug level
  void SetDebugLevel(Int_t iDebug) {fgDebug = iDebug;} // set debug level
  void FillMUONTrack(void); // set track parameters as for AliMUONTrack
  void SetTrackParam(AliMUONTrackParam *trackParam, TMatrixD *par, Double_t z); // fill AliMUONTrackParam object

  // What is necessary for sorting TClonesArray's
  Bool_t IsSortable() const { return kTRUE; }
  Int_t Compare(const TObject* trackK) const; // "Compare" function for sorting



 protected:

  AliMUONTrackK& operator=(const AliMUONTrackK& source); // assignment operator

 private:
 
  static Int_t fgDebug; ///< debug level
  static Int_t fgNOfPoints; ///< number of points in event
  //static AliMUON *fgMUON; ///< pointer to MUON module  
  static AliMUONTrackReconstructorK *fgTrackReconstructor; ///< pointer to event reconstructor
  static TClonesArray *fgHitForRec; ///< pointer to hits
  static AliMUONEventRecoCombi *fgCombi; ///< pointer to combined cluster/track finder

  AliMUONSegment *fStartSegment; ///< seed segment  
  Double_t fPosition; ///< Z-coordinate of track
  Double_t fPositionNew; //!< Z-coordinate of track
  Double_t fChi2; ///< Chi2 of track
  TObjArray *fTrackHits; ///< pointer to hits on track
  Int_t fNmbTrackHits; ///< number of points on track
  Int_t fTrackDir; ///< track direction (+(-) towards high (low) z)
  Bool_t fBPFlag; ///< backpropagation flag (TRUE if backpropagation)
  Int_t fRecover; ///< recover flag (!=0 if recovery procedure was applied)
  AliMUONHitForRec *fSkipHit; ///< hit to skip during recovery procedure

  TMatrixD *fTrackPar; ///< track parameters
  TMatrixD *fTrackParNew; //!< track parameters
  TMatrixD *fCovariance; ///< covariance matrix
  TMatrixD *fWeight; //!< weight matrix (inverse of covariance)

  // For smoother
  TObjArray *fParExtrap; //!< extrapolated track parameters
  TObjArray *fParFilter; //!< filtered track parameters
  TObjArray *fParSmooth; //!< smoothed track parameters

  TObjArray *fCovExtrap; //!< extrapolated covariance matrices
  TObjArray *fCovFilter; //!< filtered covariance matrices

  TObjArray *fJacob; //!< Jacobian matrices
  Int_t fNSteps;     //!< number of the track propagation points
  TArrayD *fSteps;   //!< Z-coordinates of the track propagation points
  TArrayD *fChi2Array; //!< measurements' contributions to the track Chi2
  TArrayD *fChi2Smooth;//!< measurements' contributions to the smoothed track Chi2

  // Functions

  AliMUONTrackK (const AliMUONTrackK& source); // copy constructor
  void EvalCovariance(Double_t dZ);
  void ParPropagation(Double_t zEnd);
  void WeightPropagation(Double_t zEnd, Bool_t smooth);
  void MSThin(Int_t sign);
  void MSLine(Double_t dZ, Double_t X0);
  Bool_t FindPoint(Int_t ichamb, Double_t zEnd, Int_t currIndx, Int_t iFB, AliMUONHitForRec *&hitAdd, Int_t iz);
  void TryPoint(TMatrixD &point, const TMatrixD &pointWeight, TMatrixD &trackParTmp, Double_t &dChi2);
  void SetGeantParam(Double_t *VGeant3, Int_t iFB);
  void GetFromGeantParam(Double_t *VGeant3, Int_t iFB);
  Bool_t Recover(void);
  void AddMatrices(AliMUONTrackK *trackK, Double_t dChi2, AliMUONHitForRec *hitAdd);
  void CreateMatrix(TObjArray *objArray) const;
  void RemoveMatrices(Double_t zEnd);
  void RemoveMatrices(AliMUONTrackK* trackK);
  void Outlier();
  void SortHits(Int_t iflag, TObjArray *array);
  void DropBranches(Int_t imax, TObjArray *hits);
  void DropBranches(AliMUONSegment *segment);
  Bool_t ExistDouble(AliMUONHitForRec *hit);
  Bool_t ExistDouble(void);

  private:
   // Some constants
   static const Int_t fgkSize; ///< number of track parameters
   static const Int_t fgkNSigma; ///< 4; acceptance window width in sigmas
   static const Double_t fgkChi2max; ///< 25; chi2 cut in smoother for outlier detection
   static const Int_t fgkTriesMax; ///< max number of attempts to find exact position during tracking
   static const Double_t fgkEpsilon; ///< tracking precision (cm)

  ClassDef(AliMUONTrackK,0) // Kalman track in MUON arm of ALICE
    };
#endif
