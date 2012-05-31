#ifndef ALIHLTMUONFULLTRACKER_H
#define ALIHLTMUONFULLTRACKER_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

///
///  @file   AliHLTMUONFullTracker.h
///  @author Indranil Das <indra.das@saha.ac.in> | <indra.ehep@gmail.com>
///  @date   09 Feb 2010
///  @brief  For full tracking in the dimuon HLT.
///


/**********************************************************************
 Created on : 08/12/2009
 Purpose    : First version implementation of the Full tracker for dHLT.
 Author     : Indranil Das, HEP Division, SINP
 Email      : indra.das@saha.ac.in | indra.ehep@gmail.com
**********************************************************************/

#include <iostream>
#include <map>

#include "TMatrixD.h"
#include "TMath.h"
#include "AliHLTMUONDataBlockReader.h"
#include "AliHLTMUONDataBlockWriter.h"


class	AliHLTMUONUtils;
class	AliHLTMUONDataTypes;
class	AliHLTMUONConstants;
class	AliRunInfo;
class	AliLog;
class	AliCDBEntry;
class	AliMpDEIterator;
class	AliMpCDB;
class	AliMpSegmentation;
class	AliMpDDLStore;
class	AliMUONTrackParam;
class	AliMUONGeometryTransformer;
class	AliMUONConstants;
class	AliMUONTrackExtrap;
class	AliMUONTrackParam;
class	AliMUONTrackExtrap;
class	AliGRPObject;
class	AliGeomManager;
class	AliCDBStorage;
class	AliCDBManager;
class	AliMagF;
class	TGeoGlobalMagField;
class	TVector3;
class	TString;
class	TMap;

typedef std::map<Int_t, Int_t> DetElemList;

class AliHLTMUONFullTracker : public AliHLTLogging
{

 public :
  ///Constructor
  AliHLTMUONFullTracker() ;
  ///Destructor
  virtual ~AliHLTMUONFullTracker();

  ///Print message
  void Print();
  ///Set the input of trigrec blocks
  Bool_t SetInput(AliHLTInt32_t ddl, const AliHLTMUONTriggerRecordStruct  *data, AliHLTInt32_t size);
  ///Set the input of rechit blocks
  Bool_t SetInput(AliHLTInt32_t ddl, const AliHLTMUONRecHitStruct  *data, AliHLTInt32_t size);
  ///Main method to run and compute the tracks
  Bool_t Run(AliHLTEventID_t iEvent,AliHLTMUONTrackStruct *data, AliHLTUInt32_t& size);
  ///To be called once from DoInit method of component
  Bool_t Init();
  ///Max number of points per chamber
  int MaxNofPointsPerCh(){return fgkMaxNofPointsPerCh;}
  ///Set for fast tracking bypass  default Kalman tracking
  void FastTracking(Bool_t isFast){fFastTracking = isFast;}
  ///Getter for fast tracking
  Bool_t FastTracking(){return fFastTracking;}

 protected:

  /// copy constructor
  AliHLTMUONFullTracker(const AliHLTMUONFullTracker& rhs); 
  /// assignment operator
  AliHLTMUONFullTracker& operator=(const AliHLTMUONFullTracker& rhs); 

 private :

  /// intger pair needed for QuadTrackSeg method
  struct IntPair{
    Int_t fFirst,fSecond; /// First and second
  };
  ///Structure for internal track segments
  struct TrackSeg{
    Int_t fIndex[4]; /// index array for cluster address
    AliHLTInt32_t fTrigRec; /// trigrec
  };

  ///Sructure for clusters
  struct Cluster{
    Float_t fX,fY,fZ;  /// position
    Float_t fErrX2,fErrY2;  /// error in position
  };
 
  ///Sructure for clusters
  struct HalfTrack{
    Float_t fPx,fPy,fPz; /// momentum
    Int_t fCharge;       /// charge
  }; 

  static const Float_t fgkTrackDetCoordinate[3]; /// set the constant value for third station position and size
  
  static const Double_t fgkAbsoedge[4] ;     /// edge of the absorber
  static const Double_t fgkRadLen[3] ;       /// radiation length of the main three matirials of the front absorber
  static const Double_t fgkRho[3] ;          /// density of the main three matirials of the front absorber
  static const Double_t fgkAtomicZ[3] ;      /// atomic number the main three matirials of the front absorber
  static const Double_t fgkAtomicA[3] ;      /// atomic mass of the main three matirials of the front absorber

  static const Int_t fgkMaxNofCellsPerCh ;      /// maximum number of cell are allowed to create
  static const Int_t fgkMaxNofPointsPerCh ;     /// maximim number of points per chamber
  static const Int_t fgkMaxNofCh ;              /// maximum number of chambrs
  static const Int_t fgkMaxNofTracks;           /// maximum number of allowed tracks
  static const Int_t fgkMaxNofConnectedTracks;  /// maximum number of back to front connected tracks
  static const Int_t fgkMaxNofTriggers;         /// maximum number of triggers (condition comes from simulation prediction)
  
  AliMUONGeometryTransformer *fChamberGeometryTransformer; /// Pointer to AliMUONGeometryTransformer
  AliHLTMUONRecHitStruct ***fChPoint; /// array of pointer to rechit data
  AliHLTMUONTriggerRecordStruct **fChPoint11; ///array of pointer to trigrec data
  TrackSeg *fBackTrackSeg; /// track segments at the rear part of the spectrometer
  TrackSeg *fFrontTrackSeg; /// track segments close the part of interaction point  of ALICE
  Float_t *fExtrapSt3X ; /// Extrapolated x position in third station
  Float_t *fExtrapSt3Y ; /// Extrapolated y position in third station
  Float_t *fInclinationBack; /// values of inclination angle of back track segments

  Int_t *fNofConnectedfrontTrackSeg ; /// nof connected tracks in front direction for each back track segments
  Int_t **fBackToFront; /// Pointer to back to front segment mapping
  Int_t *fNofPoints ; /// Number of points for each stations
  AliMUONTrackParam *fTrackParam ; /// track parameters;
  HalfTrack *fHalfTrack; /// momentum parameters for the tracks which doesnot have tracksegment in quadrants

  Int_t fTotNofPoints; /// Total number of points received from all rechit source
  Int_t fTotTrackSeg; /// Total number of track segments
  Int_t fNofCells[2]; // nof cell count per station
  Int_t fNofbackTrackSeg; /// number of back track segments
  Int_t fNoffrontTrackSeg; /// number of front track segments
  Int_t fNofConnected ; /// number of connected track segments
  AliHLTUInt32_t fNofTracks; /// number of connected track segments
  DetElemList fDetElemList; ///Map for valid detelem
  Bool_t fFastTracking ; ///flag for fast tracking avoiding kalman
  Int_t   fNofInputs; /// Nof inputs
  Int_t   fNofTriggerInputs; /// Nof inputs
  Int_t   fNofTrackerInputs; /// Nof inputs
  Bool_t  fIsMagfield ; /// checks the status of magfield

  ///  Cross Check the inputs
  Bool_t CheckInput(AliHLTEventID_t iEvent);
  /// Slat Track segments 
  Bool_t SlatTrackSeg();
  /// Calculate preliminary momentum
  Bool_t PrelimMomCalc();
  /// Quad Track segments 
  Bool_t QuadTrackSeg();
  /// Kalman Chi2 test
  Bool_t KalmanChi2Test();
  /// track extrapolation through  dipole magnet to connect front and back track seg
  Bool_t SelectFront();
  /// Propagate tracks
  void PropagateTracks(Double_t charge, Float_t& px, Float_t& py, Float_t& pz, 
		       Float_t& xr, Float_t& yr, Float_t& zr, Float_t zprop);
  /// extrapolate to origin
  Bool_t ExtrapolateToOrigin();
  /// Clean after each run
  Bool_t Clear();


  /// Angle calculate
  inline Double_t Angle(const AliHLTMUONRecHitStruct *v1, const AliHLTMUONRecHitStruct *v2);
  /// Subtracktion of two point
  inline void Sub(const AliHLTMUONRecHitStruct *v1, const AliHLTMUONRecHitStruct *v2, AliHLTMUONRecHitStruct *v3) const;
  /// Kalman Filter
  inline Double_t KalmanFilter(AliMUONTrackParam &trackParamAtCluster, Cluster *cluster);
  /// Try onecluster
  inline Double_t TryOneCluster(const AliMUONTrackParam &trackParam, Cluster* cluster,
				AliMUONTrackParam &trackParamAtCluster, Bool_t updatePropagator);
  inline Bool_t TryOneClusterFast(const AliMUONTrackParam &trackParam, const Cluster* cluster);

  /// MCS effect correction
  inline void CorrectMCSEffectInAbsorber(AliMUONTrackParam* param,
					 Double_t xVtx, Double_t yVtx, Double_t zVtx,
					 Double_t absZBeg, 
					 Double_t f1, Double_t f2);
  /// Covariant handling function
  inline void Cov2CovP(const TMatrixD &param, TMatrixD &cov);
  /// Covariant handling function
  inline void CovP2Cov(const TMatrixD &param, TMatrixD &covP);
  /// Energy loss coreection in front absorber
  inline void CorrectELossEffectInAbsorber(AliMUONTrackParam* param, Double_t eLoss);
  /// Linear Extrapolation to Z position
  inline void LinearExtrapToZ(AliMUONTrackParam* trackParam, Double_t zEnd);
  /// Energy loss
  inline Double_t EnergyLossFluctuation2(Double_t pTotal, Double_t pathLength, Double_t rho, Double_t atomicA, Double_t atomicZ);
  /// Bethe Bloch formula of enrgy loss
  inline Double_t BetheBloch(Double_t pTotal, Double_t pathLength, Double_t rho, Double_t atomicA, Double_t atomicZ);
  
  /// Runge Kutta method of track extrapolation through mag field
  inline void OneStepRungekutta(Double_t charge, Double_t step, const Double_t* vect, Double_t* vout);
  /// Helix3 method of track extrapolation through mag field
  inline void OneStepHelix3(Double_t field, Double_t step, const Double_t *vect, Double_t *vout) const;				  
  /// Fill the tracks to output pointer
  Bool_t FillOutData(AliHLTMUONTrackStruct *data, AliHLTUInt32_t& size);
  
};
#endif // ALIHLTMUONMANSOTRACKERFSM_H
